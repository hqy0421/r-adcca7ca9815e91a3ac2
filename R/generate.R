#' Generate resamples, permutations, or simulations
#'
#' @description
#'
#' Generation creates a simulated distribution from `specify()`.
#' In the context of confidence intervals, this is a bootstrap distribution
#' based on the result of `specify()`. In the context of hypothesis testing,
#' this is a null distribution based on the result of `specify()` and
#' `hypothesize().`
#'
#' Learn more in `vignette("infer")`.
#'
#' @param x A data frame that can be coerced into a [tibble][tibble::tibble].
#' @param reps The number of resamples to generate.
#' @param type The method used to generate resamples of the observed
#'   data reflecting the null hypothesis. Currently one of
#'   `"bootstrap"`, `"permute"`, or `"draw"` (see below).
#' @param variables If `type = "permute"`, a set of unquoted column names in the
#'   data to permute (independently of each other). Defaults to only the
#'   response variable. Note that any derived effects that depend on these
#'   columns (e.g., interaction effects) will also be affected.
#' @param clusters For cluster bootstrap inference. An unquoted column name in
#'   the data giving the cluster (e.g. customer) that each row belongs to.
#'   Supplying this argument switches `type = "bootstrap"` from row-wise to
#'   cluster resampling: each bootstrap replicate draws the original number of
#'   clusters with replacement and includes every row belonging to each
#'   selected cluster, duplicating rows as many times as their cluster is
#'   selected. Rows sharing a cluster are never split. The same column must
#'   also have been supplied to [specify()] via its `clusters` argument so
#'   that it is retained through the pipeline. Defaults to `NULL`, in which
#'   case individual rows are resampled as before. Ignored (with a warning)
#'   for generation types other than `"bootstrap"`.
#' @param ... Currently ignored.
#'
#' @return A tibble containing `reps` generated datasets, indicated by the
#'   `replicate` column.
#'
#' @section Generation Types:
#'
#' The `type` argument determines the method used to create the null
#' distribution.
#'
#' \itemize{
#'   \item `bootstrap`: A bootstrap sample will be drawn for each replicate,
#'   where a sample of size equal to the input sample size is drawn (with
#'   replacement) from the input sample data.
#'   \item `permute`: For each replicate, each input value will be randomly
#'   reassigned (without replacement) to a new output value in the sample.
#'   \item `draw`: A value will be sampled from a theoretical distribution
#'   with parameter `p` specified in [hypothesize()] for each replicate. This
#'   option is currently only applicable for testing on one proportion. This
#'   generation type was previously called `"simulate"`, which has been
#'   superseded.
#' }
#'
#' @includeRmd man-roxygen/seeds.Rmd
#'
#' @examples
#' # generate a null distribution by taking 200 bootstrap samples
#' gss |>
#'  specify(response = hours) |>
#'  hypothesize(null = "point", mu = 40) |>
#'  generate(reps = 200, type = "bootstrap")
#'
#' # generate a null distribution for the independence of
#' # two variables by permuting their values 200 times
#' gss |>
#'  specify(partyid ~ age) |>
#'  hypothesize(null = "independence") |>
#'  generate(reps = 200, type = "permute")
#'
#' # generate a null distribution via sampling from a
#' # binomial distribution 200 times
#' gss |>
#' specify(response = sex, success = "female") |>
#'   hypothesize(null = "point", p = .5) |>
#'   generate(reps = 200, type = "draw") |>
#'   calculate(stat = "z")
#'
#' # more in-depth explanation of how to use the infer package
#' \dontrun{
#' vignette("infer")
#' }
#'
#' @importFrom dplyr group_by
#' @family core functions
#' @export
generate <- function(
  x,
  reps = 1,
  type = NULL,
  variables = !!response_expr(x),
  clusters = NULL,
  ...
) {
  # Check type argument, warning if necessary
  type <- sanitize_generation_type(type)
  auto_type <- sanitize_generation_type(attr(x, "type"))
  type <- if (!is.null(type)) {
    compare_type_vs_auto_type(type, auto_type, x)
  } else {
    use_auto_type(auto_type)
  }
  # For MLR workflows no auto generation type is stored (see set_params()),
  # so `use_auto_type()` may return NULL or an empty string there. Resolve
  # either to bootstrap, the only generation type available without a
  # hypothesis.
  if (is.null(type) || !nzchar(type)) {
    type <- "bootstrap"
  }
  attr(x, "type") <- type

  check_cols(x, rlang::enquo(variables), type, missing(variables))

  # Resolve the cluster variable, possibly inherited from specify()
  clusters <- resolve_clusters_arg(x, rlang::enquo(clusters), type)

  attr(x, "generated") <- TRUE

  switch(
    type,
    bootstrap = bootstrap(x, reps, clusters = clusters, ...),
    permute = {
      check_permutation_attributes(x)
      permute(x, reps, rlang::enquo(variables), ...)
    },
    draw = draw(x, reps, ...),
    simulate = draw(x, reps, ...)
  )
}

# Resolve the `clusters` argument of generate(). Cluster resampling is only
# activated when the user explicitly supplies `clusters` to generate(); a
# cluster column stored by specify() is otherwise merely retained in the data
# and ordinary row-wise resampling keeps its previous behavior. Returns a
# single column name as a character string, or NULL for row-wise resampling.
resolve_clusters_arg <- function(x, clusters_quo, type, call = caller_env()) {
  arg_expr <- rlang::get_expr(clusters_quo)
  stored <- attr(x, "clusters")

  if (is.null(arg_expr)) {
    return(NULL)
  }

  if (!rlang::is_symbolic(arg_expr)) {
    cli_abort(
      "The `clusters` argument should be an unquoted column name (not a \\
       string in quotation marks).",
      call = call
    )
  }

  arg_names <- all.vars(arg_expr)

  if (length(arg_names) != 1) {
    cli_abort(
      "The `clusters` argument should reference exactly one column.",
      call = call
    )
  }

  if (!arg_names %in% colnames(x)) {
    cli_abort(
      c(
        'The cluster variable `{arg_names}` cannot be found in this \\
         dataframe.',
        i = "Cluster columns must be supplied to `specify()` via its \\
             `clusters` argument so that they are retained in the pipeline."
      ),
      call = call
    )
  }

  if (!is.null(stored) && !identical(arg_names, stored)) {
    cli_abort(
      "The `clusters` argument to `generate()` (`{arg_names}`) does not \\
       match the cluster variable supplied to `specify()` (`{stored}`).",
      call = call
    )
  }

  cluster_name <- arg_names

  if (type != "bootstrap") {
    cli_warn(
      "The {.arg clusters} argument is only used for the \\
       {.val bootstrap} generation type and will be ignored for \\
       {.val {type}} resampling."
    )
    cluster_name <- NULL
  }

  cluster_name
}

# Check that type argument is an implemented type
sanitize_generation_type <- function(x, call = caller_env()) {
  if (is.null(x)) return(x)

  check_type(x, is.character, call = call)

  if (!x %in% c("bootstrap", "permute", "simulate", "draw")) {
    cli_abort(
      'The `type` argument should be one of "bootstrap", "permute", \\
       or "draw". See {.help [{.fun generate}](infer::generate)} for more details.',
      call = call
    )
  }

  if (x == "simulate") {
    cli_inform(
      'The `"simulate"` generation type has been renamed to `"draw"`. \\
       Use `type = "draw"` instead to quiet this message.'
    )
  }

  x
}

# Ensure that the supplied type matches what would be assumed from input
compare_type_vs_auto_type <- function(type, auto_type, x) {
  if (is.null(auto_type)) {
    return(type)
  }

  if (
    (type == "bootstrap" && has_p_param(x)) ||
      (type != "bootstrap" &&
        auto_type != type &&
        # make sure auto_type vs type difference isn't just an alias
        (any(!c(auto_type, type) %in% c("draw", "simulate"))))
  ) {
    cli_warn(
      "You have given `type = \"{type}\"`, but `type` is expected \\
         to be `\"{auto_type}\"`. This workflow is untested and \\
         the results may not mean what you think they mean."
    )
  }

  type
}

has_p_param <- function(x) {
  if (!has_attr(x, "params")) {
    return(FALSE)
  }

  if (all(grepl("^p\\.", names(attr(x, "params"))))) {
    return(TRUE)
  }

  FALSE
}

use_auto_type <- function(auto_type) {
  if (!suppress_infer_messages()) {
    cli_inform('Setting `type = "{auto_type}"` in `generate()`.')
  }
  auto_type
}

check_permutation_attributes <- function(x, call = caller_env()) {
  if (
    any(!has_attr(x, "response"), !has_attr(x, "explanatory")) &&
      !identical(attr(x, "null"), "paired independence")
  ) {
    cli_abort(
      "Please {.fun specify} an explanatory and a response variable \\
        when permuting.",
      call = call
    )
  }
}

check_cols <- function(
  x,
  variables,
  type,
  missing,
  arg_name = "variables",
  call = caller_env()
) {
  if (!rlang::is_symbolic(rlang::get_expr(variables)) && type == "permute") {
    cli_abort(
      "The {.arg {arg_name}} argument should be one or more unquoted variable names \\
        (not strings in quotation marks).",
      call = call
    )
  }

  if (!missing && type != "permute") {
    cli_warn(
      'The {.arg {arg_name}} argument is only relevant for the "permute" \\
       generation type and will be ignored.'
    )

    should_prompt <- FALSE
  } else {
    should_prompt <- TRUE
  }

  col_names <- process_variables(variables, should_prompt)

  if (any(!col_names %in% colnames(x))) {
    bad_cols <- col_names[!col_names %in% colnames(x)]

    cli_abort(
      '{qty(bad_cols)}The column{?s} {.field {bad_cols}} provided to \\
       the {.arg {arg_name}} argument{qty(bad_cols)} {?is/are} not in the supplied data.',
      call = call
    )
  }
}

bootstrap <- function(x, reps = 1, clusters = NULL, ...) {
  # Check if hypothesis test chosen
  if (is_hypothesized(x)) {
    # If so, shift the variable chosen to have a mean corresponding
    # to that specified in `hypothesize`
    if (!is.null(attr(attr(x, "params"), "names"))) {
      if (identical(attr(attr(x, "params"), "names"), "mu")) {
        col <- response_name(x)
        x[[col]] <- x[[col]] - mean(x[[col]], na.rm = TRUE) + attr(x, "params")
      } # Similarly for median
      else if (identical(attr(attr(x, "params"), "names"), "med")) {
        col <- response_name(x)
        x[[col]] <- x[[col]] -
          stats::median(x[[col]], na.rm = TRUE) +
          attr(x, "params")
      }
    }
  }

  if (!is.null(clusters)) {
    return(cluster_bootstrap(x, reps, clusters))
  }

  # Set variables for use in calculate()
  result <- rep_sample_n(x, size = nrow(x), replace = TRUE, reps = reps)
  result <- copy_attrs(to = result, from = x)
  # Row-wise resampling never carries cluster bootstrap bookkeeping; the
  # cluster column itself (if specify() retained one) stays in the data, but
  # the absence of `cluster_sampling` marks this as a row-wise resample.
  attr(result, "cluster_sampling") <- NULL
  attr(result, "cluster_levels") <- NULL
  attr(result, "cluster_data") <- NULL

  append_infer_class(result)
}

#' @importFrom dplyr bind_rows group_by
permute <- function(x, reps = 1, variables, ..., call = caller_env()) {
  nrow_x <- nrow(x)
  df_out <- replicate(
    reps,
    permute_once(x, variables, call = call),
    simplify = FALSE
  ) |>
    dplyr::bind_rows() |>
    dplyr::mutate(replicate = rep(1:reps, each = !!nrow_x)) |>
    group_by_replicate(reps, nrow_x)

  df_out <- copy_attrs(to = df_out, from = x)

  append_infer_class(df_out)
}

permute_once <- function(x, variables, ..., call = caller_env()) {
  dots <- list(...)
  null <- attr(x, "null")

  if (
    !is_hypothesized(x) ||
      !null %in% c("independence", "paired independence")
  ) {
    cli_abort(
      "Permuting should be done only when doing an independence \\
        hypothesis test. See {.help [{.fun hypothesize}](infer::hypothesize)}.",
      call = call
    )
  }

  variables <- process_variables(variables, FALSE)
  if (null == "independence") {
    # for each column, determine whether it should be permuted
    needs_permuting <- colnames(x) %in% variables

    # pass each to permute_column with its associated logical
    out <- purrr::map2(x, needs_permuting, permute_column)
    out <- tibble::new_tibble(out)
  } else {
    out <- x
    signs <- sample(c(-1, 1), nrow(x), replace = TRUE, prob = c(.5, .5))
    out[[variables]] <- x[[variables]] * signs
  }

  copy_attrs(out, x)

  return(out)
}

process_variables <- function(variables, should_prompt) {
  # extract the expression and convert each element to string
  out <- rlang::get_expr(variables)

  if (length(out) == 1) {
    out <- as.character(out)
  } else {
    out <- as.list(out)
    out <- purrr::map(out, as.character)
  }

  # drop c()
  out[out == "c"] <- NULL

  # drop interactions and message
  interactions <- purrr::map_lgl(out, `%in%`, x = "*")

  if (any(interactions) && should_prompt) {
    cli_inform(
      "Message: Please supply only data columns to the {.arg variables} argument. \\
       Note that any derived effects that depend on these columns will also \\
       be affected."
    )
  }

  out <- out[!interactions]

  out
}

permute_column <- function(col, permute) {
  if (permute) {
    sample(col, size = length(col), replace = FALSE)
  } else {
    col
  }
}

#' @importFrom dplyr pull
#' @importFrom tibble tibble
#' @importFrom rlang :=
draw <- function(x, reps = 1, ...) {
  fct_levels <- as.character(unique(response_variable(x)))

  probs <- format_params(x)
  col_simmed <- unlist(replicate(
    reps,
    sample(fct_levels, size = nrow(x), replace = TRUE, prob = probs),
    simplify = FALSE
  ))

  x_nrow <- nrow(x)
  rep_tbl <- tibble::tibble(
    !!response_expr(x) := as.factor(col_simmed),
    replicate = as.factor(rep(1:reps, rep(x_nrow, reps)))
  )

  rep_tbl <- copy_attrs(to = rep_tbl, from = x)

  rep_tbl <- group_by_replicate(rep_tbl, reps, nrow(x))

  append_infer_class(rep_tbl)
}
