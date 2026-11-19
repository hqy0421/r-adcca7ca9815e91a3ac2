#' Specify response and explanatory variables
#'
#' @description
#'
#' `specify()` is used to specify which columns in the supplied data frame are
#' the relevant response (and, if applicable, explanatory) variables. Note that
#' character variables are converted to `factor`s.
#'
#' Learn more in `vignette("infer")`.
#'
#' @param x A data frame that can be coerced into a [tibble][tibble::tibble].
#' @param formula A formula with the response variable on the left and the
#'   explanatory on the right. Alternatively, a `response` and `explanatory`
#'   argument can be supplied.
#' @param response The variable name in `x` that will serve as the response.
#'   This is an alternative to using the `formula` argument.
#' @param explanatory The variable name in `x` that will serve as the
#'   explanatory variable. This is an alternative to using the formula argument.
#' @param success The level of `response` that will be considered a success, as
#'   a string. Needed for inference on one proportion, a difference in
#'   proportions, and corresponding z stats.
#' @param clusters For cluster bootstrap inference. An unquoted column name in
#'   `x` giving the cluster (e.g. customer) that each row belongs to. Rows
#'   sharing a cluster value will never be split across resamples. When
#'   supplied, the column is retained through the pipeline so that it can be
#'   used by [generate()]. Defaults to `NULL`, in which case rows are treated
#'   as independent.
#'
#' @return A tibble containing the response (and explanatory, if specified)
#'   variable data.
#'
#' @examples
#' # specifying for a point estimate on one variable
#' gss |>
#'    specify(response = age)
#'
#' # specify a relationship between variables as a formula...
#' gss |>
#'   specify(age ~ partyid)
#'
#' # ...or with named arguments!
#' gss |>
#'   specify(response = age, explanatory = partyid)
#'
#' # more in-depth explanation of how to use the infer package
#' \dontrun{
#' vignette("infer")
#' }
#'
#' @importFrom rlang f_lhs f_rhs get_expr caller_env
#' @importFrom dplyr select any_of across
#' @importFrom methods hasArg
#' @family core functions
#' @export
specify <- function(
  x,
  formula,
  response = NULL,
  explanatory = NULL,
  success = NULL,
  clusters = NULL
) {
  check_type(x, is.data.frame)

  # Standardize variable types
  x <- standardize_variable_types(x)

  # Parse response and explanatory variables
  response <- enquo(response)
  explanatory <- enquo(explanatory)

  x <- parse_variables(x, formula, response, explanatory)

  # Parse the (optional) cluster variable. Parsed before column selection so
  # that the cluster column can be carried through the pipeline.
  cluster_name <- parse_clusters(x, enquo(clusters))
  attr(x, "clusters") <- cluster_name

  # Add attributes
  attr(x, "success") <- success
  attr(x, "generated") <- FALSE
  attr(x, "hypothesized") <- FALSE
  attr(x, "fitted") <- FALSE

  # Check the success argument
  check_success_arg(x, success)

  # Select variables, retaining the cluster column if one was supplied
  x <- x |>
    select(any_of(c(
      response_name(x),
      explanatory_name(x),
      cluster_name
    )))

  is_complete <- stats::complete.cases(x)
  if (!all(is_complete)) {
    x <- dplyr::filter(x, is_complete)
    cli_warn("Removed {sum(!is_complete)} rows containing missing values.")
  }

  # Add "infer" class
  append_infer_class(x)
}

# Parse the `clusters` argument of specify() into a single column name
# (a character string), or NULL when no cluster variable was supplied.
parse_clusters <- function(x, clusters_quo, call = caller_env()) {
  clusters_expr <- rlang::get_expr(clusters_quo)

  if (is.null(clusters_expr)) {
    return(NULL)
  }

  if (!rlang::is_symbolic(clusters_expr)) {
    cli_abort(
      "The `clusters` argument should be an unquoted column name (not a \\
       string in quotation marks).",
      call = call
    )
  }

  cluster_names <- all.vars(clusters_expr)

  if (length(cluster_names) != 1) {
    cli_abort(
      "The `clusters` argument should reference exactly one column.",
      call = call
    )
  }

  if (!cluster_names %in% colnames(x)) {
    cli_abort(
      'The cluster variable `{cluster_names}` cannot be found in this \\
       dataframe.',
      call = call
    )
  }

  cluster_names
}

parse_variables <- function(
  x,
  formula,
  response,
  explanatory,
  call = caller_env()
) {
  if (methods::hasArg(formula)) {
    tryCatch(
      rlang::is_formula(formula),
      error = function(e) {
        cli_abort(
          c(
            "The argument you passed in for the formula does not exist.",
            i = "Were you trying to pass in an unquoted column name?",
            i = "Did you forget to name one or more arguments?"
          ),
          call = call
        )
      }
    )
    if (!rlang::is_formula(formula)) {
      cli_abort(
        c(
          "The first unnamed argument must be a formula.",
          i = "You passed in '{get_type(formula)}'.",
          x = "Did you forget to name one or more arguments?"
        ),
        call = call
      )
    }
  }

  attr(x, "response") <- get_expr(response)
  attr(x, "explanatory") <- get_expr(explanatory)
  attr(x, "formula") <- NULL

  if (methods::hasArg(formula)) {
    attr(x, "response") <- f_lhs(formula)
    attr(x, "explanatory") <- f_rhs(formula)
    attr(x, "formula") <- formula
  }

  # Check response and explanatory variables to be appropriate for later use
  if (!has_response(x)) {
    cli_abort(
      "Please supply a response variable that is not `NULL`.",
      call = call
    )
  }

  check_var_correct(x, "response", call = call)
  check_var_correct(x, "explanatory", call = call)

  # If there's an explanatory var
  check_vars_different(x, call = call)

  if (!has_attr(x, "response")) {
    attr(x, "response_type") <- NULL
  } else {
    attr(x, "response_type") <- class(response_variable(x))
  }

  if (!has_attr(x, "explanatory")) {
    attr(x, "explanatory_type") <- NULL
  } else {
    attr(x, "explanatory_type") <-
      purrr::map_chr(as.data.frame(explanatory_variable(x)), class)
  }

  attr(x, "type_desc_response") <- determine_variable_type(x, "response")
  attr(x, "type_desc_explanatory") <- determine_variable_type(x, "explanatory")

  # Determine params for theoretical fit
  x <- set_params(x)

  x
}

check_success_arg <- function(x, success, call = caller_env()) {
  response_col <- response_variable(x)

  if (!is.null(success)) {
    if (!is.character(success)) {
      cli_abort("`success` must be a string.", call = call)
    }
    if (!is.factor(response_col)) {
      cli_abort(
        "`success` should only be specified if the response is a categorical \\
         variable.",
        call = call
      )
    }
    if (!(success %in% levels(response_col))) {
      cli_abort(
        '{success} is not a valid level of {response_name(x)}.',
        call = call
      )
    }
    if (sum(table(response_col) > 0) > 2) {
      cli_abort(
        "`success` can only be used if the response has two levels. \\
         `filter()` can reduce a variable to two levels.",
        call = call
      )
    }
  }

  if (
    (attr(x, "response_type") == "factor" &&
      is.null(success) &&
      length(levels(response_variable(x))) == 2) &&
      ((!has_attr(x, "explanatory_type") ||
        length(levels(explanatory_variable(x))) == 2))
  ) {
    cli_abort(
      'A level of the response variable `{response_name(x)}` needs to be \\
        specified for the `success` argument in `specify()`.',
      call = call
    )
  }
}

check_var_correct <- function(x, var_name, call = caller_env()) {
  var <- attr(x, var_name)

  # Variable (if present) should be a symbolic column name
  if (!is.null(var)) {
    if (!rlang::is_symbolic(var)) {
      cli_abort(
        "The {var_name} should be a bare variable name (not a string in \\
         quotation marks).",
        call = call
      )
    }

    if (any(!(all.vars(var) %in% names(x)))) {
      cli_abort(
        'The {var_name} variable `{var}` cannot be found in this dataframe.',
        call = call
      )
    }
  }

  TRUE
}

check_vars_different <- function(x, call = caller_env()) {
  if (has_response(x) && has_explanatory(x)) {
    if (identical(response_name(x), explanatory_name(x))) {
      cli_abort(
        "The response and explanatory variables must be different from one \\
         another.",
        call = call
      )
    }
  }

  TRUE
}
