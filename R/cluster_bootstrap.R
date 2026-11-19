# Cluster bootstrap and cluster-aware BCa confidence intervals -------------
#
# This file implements resampling, fitting, and BCa interval machinery in
# which whole clusters (e.g. customers) are the resampling unit. Every fit in
# the pipeline -- the observed fit, each bootstrap replicate fit, and each
# leave-one-cluster-out fit -- goes through the same canonical fitting helper
# so that factor encodings and coefficient identities never drift between
# fits.

# Generation ---------------------------------------------------------------

#' Resample whole clusters with replacement
#'
#' For each replicate, draw as many clusters with replacement as there are
#' clusters in the original data. Every row belonging to a selected cluster is
#' included, duplicated as many times as the cluster was selected. Rows that
#' share a cluster are never split across resamples.
#'
#' @return A grouped `infer` tibble with a `replicate` column. The original
#'   observed data, the cluster levels, and a table recording which clusters
#'   were drawn in each replicate are stored as attributes for downstream
#'   fitting and confidence interval functions.
#'
#' @noRd
cluster_bootstrap <- function(x, reps, cluster_name, call = caller_env()) {
  check_type(
    reps,
    is_single_number,
    "single number not less than 1",
    min_val = 1,
    call = call
  )

  cluster_values <- as.vector(x[[cluster_name]])
  cluster_levels <- unique(cluster_values)
  n_clusters <- length(cluster_levels)

  if (n_clusters < 3) {
    cli_abort(
      c(
        "Cluster bootstrap requires at least three clusters.",
        i = "The data contain {n_clusters} cluster{?s} in `{cluster_name}`; \\
             leave-one-cluster-out acceleration is not defined for fewer \\
             than three clusters."
      ),
      call = call
    )
  }

  # Row indices belonging to each cluster, in original row order. Building a
  # factor with the observed levels keeps empty clusters (none here, but the
  # same structure is reused for leave-one-cluster-out data) represented.
  rows_by_cluster <- split(
    seq_len(nrow(x)),
    factor(cluster_values, levels = cluster_levels)
  )

  # One column per cluster draw; row r gives the cluster indices of replicate r.
  sampled <- matrix(
    sample.int(n_clusters, size = n_clusters * reps, replace = TRUE),
    nrow = reps,
    ncol = n_clusters,
    byrow = TRUE
  )

  rep_tbls <- vector("list", reps)
  for (r in seq_len(reps)) {
    row_idx <- unlist(
      rows_by_cluster[sampled[r, ]],
      use.names = FALSE,
      recursive = TRUE
    )
    rep_tbls[[r]] <- vctrs::vec_slice(x, row_idx)
    rep_tbls[[r]]$replicate <- r
  }

  result <- dplyr::bind_rows(rep_tbls)
  result <- dplyr::group_by(result, replicate)

  # Record the resampled population: every cluster draw for every replicate.
  cluster_sampling <- tibble::tibble(
    replicate = rep(seq_len(reps), each = n_clusters),
    cluster = cluster_levels[as.vector(t(sampled))]
  )

  result <- copy_attrs(to = result, from = x)
  attr(result, "clusters") <- cluster_name
  attr(result, "cluster_levels") <- cluster_levels
  attr(result, "cluster_sampling") <- cluster_sampling
  # The observed data are needed for the leave-one-cluster-out fits, which
  # cannot be reconstructed from the resampled replicates.
  attr(result, "cluster_data") <- x

  append_infer_class(result)
}

# Canonical fitting --------------------------------------------------------

# Capture factor levels from the observed data so that the same reference
# encoding can be imposed on every resampled and leave-one-cluster-out data set.
capture_factor_levels <- function(data) {
  fct_cols <- names(data)[vapply(data, is.factor, logical(1))]
  lapply(data[fct_cols], levels)
}

restore_factor_levels <- function(data, level_ref) {
  for (nm in names(level_ref)) {
    data[[nm]] <- factor(as.character(data[[nm]]), levels = level_ref[[nm]])
  }

  data
}

# Coefficient names in glm's raw encoding, e.g. "(Intercept)", "collegedegree".
canonical_raw_terms <- function(data, formula) {
  mf <- stats::model.frame(formula, data)
  mm <- stats::model.matrix(stats::terms(mf), mf)

  colnames(mm)
}

map_raw_term <- function(x) {
  dplyr::case_when(x == "(Intercept)" ~ "intercept", TRUE ~ x)
}

# Fit one model and return estimates aligned, by name, to `raw_terms`.
# Non-estimable coefficients (singular/aliased columns or a failed fit) are
# returned as NA_real_ rather than dropped, preserving coefficient identity.
# A list with the estimate vector and a status description is returned.
fit_canonical_once <- function(
  data,
  formula,
  raw_terms,
  level_ref,
  dots
) {
  data <- restore_factor_levels(data, level_ref)

  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      do.call(stats::glm, c(list(formula = formula, data = data), dots)),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  status <- list(
    converged = FALSE,
    error = NA_character_,
    warning = if (length(warnings) == 0) NA_character_ else warnings[[1]]
  )

  estimates <- setNames(rep(NA_real_, length(raw_terms)), raw_terms)

  if (inherits(fit, "error")) {
    status$error <- conditionMessage(fit)

    return(list(estimate = estimates, status = status))
  }

  status$converged <- is.null(fit$converged) || isTRUE(fit$converged)
  if (!status$converged && (is.na(status$error) || is.null(status$error))) {
    status$error <- "model did not converge"
  }

  coefs <- stats::coef(fit)
  matched <- intersect(names(coefs), raw_terms)
  estimates[matched] <- unname(coefs[matched])

  # Singular/aliased columns (e.g. a factor level absent from a resample)
  # yield NA coefficients even when the fit otherwise converges. Record these
  # explicitly as non-estimable rather than letting them pass silently.
  missing_terms <- raw_terms[is.na(estimates)]
  if (length(missing_terms) > 0) {
    status$error <- if (!is.na(status$error) && !is.null(status$error)) {
      paste0(
        status$error,
        "; coefficient(s) not estimable: ",
        paste(map_raw_term(missing_terms), collapse = ", ")
      )
    } else {
      paste0(
        "coefficient(s) not estimable: ",
        paste(map_raw_term(missing_terms), collapse = ", ")
      )
    }
  }

  list(estimate = estimates, status = status)
}

# Fit the observed model, every bootstrap replicate, and a model with each
# complete cluster removed in turn, all through fit_canonical_once() so that
# encoding and coefficient identity are shared across all three families.
fit_cluster_bootstrap <- function(object, formula, dots) {
  # The observed data stored at generation time carry the same infer
  # attributes as the pipeline object, so apply the same response releveling
  # the ordinary fit path uses (see fit.infer). Factor levels captured below
  # then impose that exact reference encoding on every resampled and
  # leave-one-cluster-out fit.
  obs_data <- relevel_response(attr(object, "cluster_data"))
  cluster_name <- attr(object, "clusters")
  cluster_levels <- attr(object, "cluster_levels")

  level_ref <- capture_factor_levels(obs_data)
  raw_terms <- canonical_raw_terms(obs_data, formula)
  out_terms <- map_raw_term(raw_terms)
  n_terms <- length(raw_terms)

  # Observed fit -----------------------------------------------------------
  observed <- fit_canonical_once(obs_data, formula, raw_terms, level_ref, dots)

  if (!is.null(observed$status$error) && !is.na(observed$status$error)) {
    cli_abort(
      c(
        "The observed model could not be fit to the full data.",
        x = "{observed$status$error}"
      )
    )
  }

  # Replicate fits ---------------------------------------------------------
  nested <- object |>
    dplyr::ungroup() |>
    tidyr::nest(data = -replicate)

  reps <- nrow(nested)
  rep_estimates <- matrix(
    NA_real_,
    nrow = reps,
    ncol = n_terms,
    dimnames = list(NULL, raw_terms)
  )
  fit_status <- tibble::tibble(
    replicate = integer(reps),
    converged = logical(reps),
    error = character(reps),
    warning = character(reps)
  )

  for (i in seq_len(reps)) {
    one <- fit_canonical_once(
      nested$data[[i]],
      formula,
      raw_terms,
      level_ref,
      dots
    )
    rep_estimates[i, ] <- one$estimate
    fit_status$replicate[i] <- nested$replicate[[i]]
    fit_status$converged[i] <- one$status$converged
    fit_status$error[i] <- if (is.null(one$status$error)) {
      NA_character_
    } else {
      one$status$error
    }
    fit_status$warning[i] <- one$status$warning
  }

  failed_reps <- fit_status$replicate[
    vapply(
      seq_len(reps),
      function(i) {
        !fit_status$converged[i] ||
          (!is.na(fit_status$error[i])) ||
          any(!is.finite(rep_estimates[i, ]))
      },
      logical(1)
    )
  ]

  if (length(failed_reps) > 0) {
    cli_warn(
      c(
        "{length(failed_reps)} bootstrap replicate{?s} did not yield finite \\
         estimates for at least one coefficient and {?is/are} recorded as \\
         `NA` in the output.",
        i = "See {.code attr(<fit>, \"fit_status\")} for details. BCa \\
             confidence intervals will not be reported until every replicate \\
             is estimable."
      )
    )
  }

  replicate_tbl <- tibble::tibble(
    replicate = rep(nested$replicate, each = n_terms),
    term = rep(out_terms, times = reps),
    estimate = as.vector(t(rep_estimates))
  )

  # Leave-one-cluster-out fits --------------------------------------------
  cluster_values <- as.vector(obs_data[[cluster_name]])
  loo_estimates <- matrix(
    NA_real_,
    nrow = length(cluster_levels),
    ncol = n_terms,
    dimnames = list(NULL, raw_terms)
  )
  loo_status <- character(length(cluster_levels))

  for (j in seq_along(cluster_levels)) {
    keep <- cluster_values != cluster_levels[[j]]
    one <- fit_canonical_once(
      obs_data[keep, , drop = FALSE],
      formula,
      raw_terms,
      level_ref,
      dots
    )
    loo_estimates[j, ] <- one$estimate
    loo_status[[j]] <- if (is.null(one$status$error) ||
      is.na(one$status$error)) {
      NA_character_
    } else {
      one$status$error
    }
  }

  failed_clusters <- cluster_levels[
    vapply(
      seq_along(cluster_levels),
      function(j) any(!is.finite(loo_estimates[j, ])),
      logical(1)
    )
  ]

  if (length(failed_clusters) > 0) {
    cli_warn(
      c(
        "Leave-one-cluster-out fits were not estimable for {length(failed_clusters)} \\
         cluster{?s}: {.val {failed_clusters}}.",
        i = "BCa acceleration (and thus BCa intervals) for coefficients \\
             depending on these fits is not available; see \\
             {.code attr(<fit>, \"jackknife\")}."
      )
    )
  }

  jackknife_tbl <- tibble::tibble(
    cluster = rep(cluster_levels, each = n_terms),
    term = rep(out_terms, times = length(cluster_levels)),
    estimate = as.vector(t(loo_estimates)),
    error = rep(loo_status, each = n_terms)
  )

  observed_tbl <- tibble::tibble(
    term = out_terms,
    estimate = unname(observed$estimate)
  )

  replicate_tbl <- copy_attrs(to = replicate_tbl, from = object)
  attr(replicate_tbl, "fitted") <- TRUE
  attr(replicate_tbl, "clusters") <- cluster_name
  attr(replicate_tbl, "cluster_levels") <- cluster_levels
  attr(replicate_tbl, "cluster_sampling") <- attr(object, "cluster_sampling")
  attr(replicate_tbl, "cluster_data") <- obs_data
  attr(replicate_tbl, "jackknife") <- jackknife_tbl
  attr(replicate_tbl, "fit_status") <- fit_status
  attr(replicate_tbl, "observed_fit_canonical") <- observed_tbl

  replicate_tbl
}

# BCa confidence intervals -------------------------------------------------

# Translate the `contrast` argument into a matrix of contrast weights with one
# row per requested estimand and canonical term names as columns.
process_contrast <- function(
  contrast,
  terms,
  call = caller_env()
) {
  p <- length(terms)

  if (is.null(contrast)) {
    weights <- diag(1, nrow = p, ncol = p)
    dimnames(weights) <- list(terms, terms)

    return(weights)
  }

  if (is.character(contrast)) {
    unknown <- setdiff(contrast, terms)
    if (length(unknown) > 0) {
      cli_abort(
        c(
          "All term names supplied to `contrast` must be coefficients of \\
           the fitted model.",
          x = "Unknown term{?s}: {.val {unknown}}.",
          i = "Available terms: {.val {terms}}."
        ),
        call = call
      )
    }

    weights <- matrix(0, nrow = length(contrast), ncol = p)
    dimnames(weights) <- list(contrast, terms)
    for (j in seq_along(contrast)) {
      weights[j, contrast[[j]]] <- 1
    }

    return(weights)
  }

  if (is.numeric(contrast)) {
    if (!is.null(names(contrast))) {
      unknown <- setdiff(names(contrast), terms)
      if (length(unknown) > 0) {
        cli_abort(
          c(
            "All names of a numeric `contrast` must be coefficients of the \\
             fitted model.",
            x = "Unknown term{?s}: {.val {unknown}}.",
            i = "Available terms: {.val {terms}}."
          ),
          call = call
        )
      }

      weights <- matrix(0, nrow = 1, ncol = p)
      dimnames(weights) <- list("contrast", terms)
      weights[1, names(contrast)] <- contrast

      return(weights)
    }

    if (length(contrast) != p) {
      cli_abort(
        "An unnamed numeric `contrast` must have length equal to the number \\
         of coefficients ({p}), not {length(contrast)}.",
        call = call
      )
    }

    weights <- matrix(contrast, nrow = 1, ncol = p)
    dimnames(weights) <- list("contrast", terms)

    return(weights)
  }

  cli_abort(
    "`contrast` must be `NULL`, a vector of coefficient term names, or a \\
     numeric vector (optionally named) defining a linear combination.",
    call = call
  )
}

# BCa acceleration from leave-one-cluster-out estimates.
bca_acceleration <- function(loo) {
  lbar <- mean(loo)
  diffs <- lbar - loo
  denom <- 6 * sum(diffs^2)^(3 / 2)

  if (denom == 0) {
    return(list(a = 0, note = "acceleration undefined: constant \\
                    leave-one-cluster-out estimates; acceleration set to 0"))
  }

  list(a = sum(diffs^3) / denom, note = NA_character_)
}

# Compute a BCa interval for one estimand. Every failure is reported rather
# than silently dropping resamples: the function errors for non-estimable
# cases instead of returning a partial interval.
bca_one <- function(theta_obs, boot, loo, level, label, call = caller_env()) {
  bad_boot <- which(!is.finite(boot))
  if (length(bad_boot) > 0) {
    cli_abort(
      c(
        "BCa interval for {.val {label}} is not estimable: {length(bad_boot)} \\
         bootstrap replicate estimate{?s} {?is/are} non-finite.",
        i = "No replicates were removed; inspect \\
             {.code attr(<fit>, \"fit_status\")} and the bootstrap estimates."
      ),
      call = call
    )
  }

  bad_loo <- which(!is.finite(loo))
  if (length(bad_loo) > 0) {
    cli_abort(
      "BCa interval for {.val {label}} is not estimable: leave-one-cluster-out \\
       estimates were non-finite for {length(bad_loo)} cluster{?s}.",
      call = call
    )
  }

  prop_below <- mean(boot < theta_obs)
  z0 <- stats::qnorm(prop_below)

  if (!is.finite(z0)) {
    cli_abort(
      c(
        "BCa interval for {.val {label}} is not estimable: the bias \\
         correction is infinite.",
        i = "{round(100 * prop_below, 2)}% of bootstrap estimates fall below \\
             the observed estimate; the BCa adjustment lies outside the \\
             range of the bootstrap distribution."
      ),
      call = call
    )
  }

  accel <- bca_acceleration(loo)
  a <- accel$a

  z_alpha <- stats::qnorm((1 + c(-level, level)) / 2)
  adj <- z0 + z_alpha
  new_probs <- stats::pnorm(z0 + adj / (1 - a * adj))

  if (any(new_probs <= 0 | new_probs >= 1)) {
    cli_abort(
      c(
        "BCa interval for {.val {label}} is not estimable: the adjusted \\
         quantile probabilities fall outside (0, 1).",
        i = "Adjusted probabilities: {.val {round(new_probs, 4)}}."
      ),
      call = call
    )
  }

  ci <- stats::quantile(boot, probs = new_probs, type = 7, names = FALSE)

  list(
    lower_ci = ci[[1]],
    upper_ci = ci[[2]],
    z0 = z0,
    acceleration = a,
    note = accel$note
  )
}

# Entry point used by get_confidence_interval() for `type = "bca"`.
ci_bca_cluster <- function(
  x,
  level,
  point_estimate,
  contrast,
  call = caller_env()
) {
  if (!is_fitted(x) || is.null(attr(x, "cluster_sampling"))) {
    cli_abort(
      c(
        '`type = "bca"` requires a fitted cluster bootstrap distribution.',
        i = "Generate one with {.code generate(..., clusters = <cluster col>)} \\
             followed by {.fun fit}.",
        i = "For ordinary row-wise bootstrap objects, use \\
             {.val percentile}, {.val se}, or {.val bias-corrected}; a \\
             row-wise jackknife acceleration is intentionally not provided."
      ),
      call = call
    )
  }

  if (is.null(point_estimate) || !is_fitted(point_estimate)) {
    cli_abort(
      'A numeric value needs to be given for `point_estimate` for \\
       `type` "bca". Pass the observed `fit()` output.',
      call = call
    )
  }

  check_mlr_x_and_obs_stat(
    x,
    point_estimate,
    "get_confidence_interval",
    "point_estimate",
    call = call
  )

  terms <- unique(x$term)
  p <- length(terms)

  # Replicate estimates as a reps x terms matrix, columns in canonical term
  # order (term order as first seen in the fitted output).
  rep_wide <- tidyr::pivot_wider(
    dplyr::ungroup(x),
    names_from = term,
    values_from = estimate,
    id_cols = replicate
  )
  rep_mat <- as.matrix(rep_wide[, terms, drop = FALSE])
  storage.mode(rep_mat) <- "double"

  jack <- attr(x, "jackknife")
  jack_wide <- tidyr::pivot_wider(
    jack[, c("cluster", "term", "estimate")],
    names_from = term,
    values_from = estimate,
    id_cols = cluster
  )
  jack_mat <- as.matrix(jack_wide[, terms, drop = FALSE])
  storage.mode(jack_mat) <- "double"

  obs <- point_estimate$estimate[match(terms, point_estimate$term)]
  if (length(obs) != p || any(is.na(obs))) {
    cli_abort(
      "The observed point estimate could not be aligned with the bootstrap \\
       coefficient terms.",
      call = call
    )
  }

  weights <- process_contrast(contrast, terms, call = call)
  labels <- rownames(weights)

  boot_c <- rep_mat %*% t(weights)
  loo_c <- jack_mat %*% t(weights)
  obs_c <- as.vector(weights %*% obs)

  out <- tibble::tibble(
    term = labels,
    lower_ci = NA_real_,
    upper_ci = NA_real_
  )
  diagnostics <- tibble::tibble(
    term = labels,
    point_estimate = obs_c,
    z0 = NA_real_,
    acceleration = NA_real_,
    note = NA_character_
  )

  notes <- character()
  for (k in seq_along(labels)) {
    res <- bca_one(
      theta_obs = obs_c[[k]],
      boot = boot_c[, k],
      loo = loo_c[, k],
      level = level,
      label = labels[[k]],
      call = call
    )
    out$lower_ci[[k]] <- res$lower_ci
    out$upper_ci[[k]] <- res$upper_ci
    diagnostics$z0[[k]] <- res$z0
    diagnostics$acceleration[[k]] <- res$acceleration
    diagnostics$note[[k]] <- res$note
    if (!is.na(res$note)) notes <- c(notes, paste0(labels[[k]], ": ", res$note))
  }

  if (length(notes) > 0) cli_warn(notes)

  boot_estimates <- tibble::as_tibble(as.data.frame(boot_c))
  names(boot_estimates) <- labels
  boot_estimates <- tibble::tibble(
    replicate = rep_wide$replicate,
    boot_estimates
  )

  jack_estimates <- tibble::as_tibble(as.data.frame(loo_c))
  names(jack_estimates) <- labels
  jack_estimates <- tibble::tibble(
    cluster = jack_wide$cluster,
    jack_estimates
  )

  attr(out, "bca") <- TRUE
  attr(out, "level") <- level
  attr(out, "point_estimate") <- diagnostics[, c("term", "point_estimate")]
  attr(out, "bootstrap_estimates") <- boot_estimates
  attr(out, "jackknife_estimates") <- jack_estimates
  attr(out, "cluster_sampling") <- attr(x, "cluster_sampling")
  attr(out, "bca_diagnostics") <- diagnostics
  attr(out, "clusters") <- attr(x, "clusters")

  out
}
