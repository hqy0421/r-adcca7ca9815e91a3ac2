# Tests for cluster bootstrap and cluster-aware BCa confidence intervals

# Build a small clustered data set: several rows per customer (cluster)
set.seed(2026)
n_clust <- 30
clust_dat <- do.call(rbind, lapply(seq_len(n_clust), function(id) {
  ni <- sample(2:6, 1)
  age <- sample(20:60, 1)
  b0 <- rnorm(1, 40, 4)
  data.frame(
    customer = paste0("c", id),
    hours = b0 + 0.05 * age + rnorm(ni, 0, 2),
    age = age,
    college = sample(c("degree", "no degree"), 1)
  )
}))
rownames(clust_dat) <- NULL

clust_spec <- clust_dat |>
  specify(hours ~ age + college, clusters = customer)

# specify() ----------------------------------------------------------------

test_that("specify() retains and records the cluster column", {
  expect_true("customer" %in% names(clust_spec))
  expect_equal(attr(clust_spec, "clusters"), "customer")

  # without clusters, behavior is unchanged
  plain <- clust_dat |> specify(hours ~ age + college)
  expect_false("customer" %in% names(plain))
  expect_null(attr(plain, "clusters"))
})

test_that("specify() validates the clusters argument", {
  expect_snapshot(
    error = TRUE,
    specify(clust_dat, hours ~ age, clusters = "customer")
  )
  expect_snapshot(
    error = TRUE,
    specify(clust_dat, hours ~ age, clusters = not_a_col)
  )
})

# generate() ---------------------------------------------------------------

test_that("cluster bootstrap resamples whole clusters, never rows", {
  set.seed(1)
  reps <- 20
  gen <- clust_spec |>
    generate(reps = reps, type = "bootstrap", clusters = customer)

  sampling <- attr(gen, "cluster_sampling")
  cluster_levels <- attr(gen, "cluster_levels")

  # resampled population: same number of clusters as the original each rep
  expect_equal(length(cluster_levels), n_clust)
  expect_equal(
    as.integer(table(sampling$replicate)),
    rep(n_clust, reps)
  )

  # every sampled cluster value is one of the original customers
  expect_true(all(sampling$cluster %in% cluster_levels))

  # replicate 1: rows are every row of each drawn cluster, in draw order
  drawn <- sampling$cluster[sampling$replicate == 1]
  manual <- do.call(
    rbind,
    lapply(drawn, function(cl) {
      clust_spec[clust_spec$customer == cl, ]
    })
  )
  got <- as.data.frame(gen[gen$replicate == 1, names(clust_spec)])
  rownames(manual) <- rownames(got) <- NULL
  expect_equal(got, manual, ignore_attr = TRUE)

  # clusters are never split: row counts per (replicate, cluster) are integer
  # multiples of that cluster's original size
  tab <- as.data.frame(table(gen$replicate, gen$customer))
  names(tab) <- c("rep", "cluster", "n")
  orig <- as.data.frame(table(clust_spec$customer))
  names(orig) <- c("cluster", "orig_n")
  merged <- merge(tab, orig, by = "cluster")
  expect_true(all(merged$n %% as.integer(merged$orig_n) == 0))
})

test_that("generate() requires explicit clusters for cluster resampling", {
  # no clusters argument -> ordinary row-wise bootstrap, fully backward
  # compatible even though the data still carry a cluster column
  row_gen <- clust_spec |> generate(reps = 5, type = "bootstrap")

  expect_null(attr(row_gen, "cluster_sampling"))
  expect_equal(nrow(row_gen), nrow(clust_spec) * 5)

  row_fit <- row_gen |> fit()
  expect_null(attr(row_fit, "cluster_sampling"))
})

test_that("generate() validates cluster arguments", {
  expect_snapshot(
    error = TRUE,
    clust_spec |>
      generate(reps = 5, type = "bootstrap", clusters = not_a_col)
  )

  few <- clust_dat[clust_dat$customer %in% c("c1", "c2"), ] |>
    specify(hours ~ age, clusters = customer)
  expect_snapshot(
    error = TRUE,
    few |> generate(reps = 5, clusters = customer)
  )

  # clusters are ignored with a warning under non-bootstrap generation types
  expect_warning(
    clust_spec |>
      hypothesize(null = "independence") |>
      generate(reps = 3, type = "permute", clusters = customer),
    "ignored"
  )
})

# fit() --------------------------------------------------------------------

test_that("fit() shares encoding and coefficient identity across all fits", {
  set.seed(2)
  obs <- clust_spec |> fit()
  boot_fit <- clust_spec |>
    generate(reps = 25, clusters = customer) |>
    fit()

  terms <- unique(obs$term)

  # identical term identities in observed, replicate, and jackknife fits
  expect_equal(unique(boot_fit$term), terms)
  expect_equal(unique(attr(boot_fit, "jackknife")$term), terms)

  # the canonical observed fit equals the ordinary fit() output
  expect_equal(
    attr(boot_fit, "observed_fit_canonical")$estimate,
    obs$estimate
  )

  # one row per (replicate, term) and one jackknife row per (cluster, term)
  expect_equal(nrow(boot_fit), 25 * length(terms))
  expect_equal(
    nrow(attr(boot_fit, "jackknife")),
    n_clust * length(terms)
  )

  # factor encoding matches the row-wise bootstrap path as well
  row_terms <- clust_spec |>
    generate(reps = 3, type = "bootstrap") |>
    fit() |>
    unique()
  expect_equal(unique(row_terms$term), terms)
})

test_that("fit() aligns binomial reference levels with the ordinary path", {
  set.seed(8)
  bin_dat <- do.call(rbind, lapply(1:25, function(id) {
    ni <- sample(2:5, 1)
    x <- rnorm(ni)
    data.frame(
      customer = paste0("c", id),
      y = factor(
        rbinom(ni, 1, 1 / (1 + exp(-(0.3 * x + rnorm(1, 0, 0.3))))),
        levels = 0:1
      ),
      x = x
    )
  }))

  sp <- bin_dat |>
    specify(y ~ x, clusters = customer, success = "1")
  obs <- sp |> fit()
  boot_fit <- sp |>
    generate(reps = 30, clusters = customer) |>
    fit()

  # the canonical observed fit must equal the ordinary fit() output even
  # though the latter relevels the binary response
  expect_equal(
    attr(boot_fit, "observed_fit_canonical")$estimate,
    obs$estimate
  )
  expect_equal(unique(boot_fit$term), obs$term)
})

test_that("fit() records cluster-deletion estimates and fit status", {
  set.seed(3)
  boot_fit <- clust_spec |>
    generate(reps = 15, clusters = customer) |>
    fit()

  jack <- attr(boot_fit, "jackknife")
  expect_setequal(unique(jack$cluster), paste0("c", seq_len(n_clust)))
  expect_true(all(is.finite(jack$estimate)))

  status <- attr(boot_fit, "fit_status")
  expect_equal(nrow(status), 15)
  expect_true(all(status$converged))
})

# BCa confidence intervals -------------------------------------------------

test_that("BCa intervals match an independent manual implementation", {
  set.seed(42)
  reps <- 400
  boot_fit <- clust_spec |>
    generate(reps = reps, clusters = customer) |>
    fit()
  obs <- clust_spec |> fit()

  ci <- get_confidence_interval(
    boot_fit,
    type = "bca",
    point_estimate = obs,
    contrast = "age"
  )

  # independent recomputation from the recorded resampled population
  cluster_col <- attr(boot_fit, "clusters")
  cdata <- attr(boot_fit, "cluster_data")
  ids <- attr(boot_fit, "cluster_levels")
  sampling <- attr(boot_fit, "cluster_sampling")

  fit_age <- function(d) {
    unname(stats::coef(stats::glm(
      get_formula(clust_spec),
      data = d,
      family = stats::gaussian
    ))["age"])
  }

  boot_est <- vapply(seq_len(reps), function(r) {
    drawn <- sampling$cluster[sampling$replicate == r]
    d2 <- do.call(rbind, lapply(drawn, function(k) {
      cdata[cdata[[cluster_col]] == k, ]
    }))
    fit_age(d2)
  }, numeric(1))

  loo_est <- vapply(ids, function(k) {
    fit_age(cdata[cdata[[cluster_col]] != k, ])
  }, numeric(1))

  theta <- fit_age(cdata)
  z0 <- stats::qnorm(mean(boot_est < theta))
  lbar <- mean(loo_est)
  a <- sum((lbar - loo_est)^3) /
    (6 * sum((lbar - loo_est)^2)^(3 / 2))
  z_a <- stats::qnorm((1 + c(-0.95, 0.95)) / 2)
  probs <- stats::pnorm(z0 + (z0 + z_a) / (1 - a * (z0 + z_a)))
  manual <- stats::quantile(boot_est, probs = probs, type = 7, names = FALSE)

  expect_equal(c(ci$lower_ci, ci$upper_ci), unname(manual))
  expect_equal(attr(ci, "bca_diagnostics")$z0, z0)
  expect_equal(attr(ci, "bca_diagnostics")$acceleration, a)

  # returned resampled population and estimate collections are consistent
  expect_equal(
    sort(attr(ci, "bootstrap_estimates")$age),
    sort(unname(boot_est))
  )
  expect_equal(
    sort(attr(ci, "jackknife_estimates")$age),
    sort(unname(loo_est))
  )
  expect_equal(attr(ci, "cluster_sampling"), sampling)
})

test_that("BCa supports term names and named or unnamed linear contrasts", {
  set.seed(4)
  boot_fit <- clust_spec |>
    generate(reps = 30, clusters = customer) |>
    fit()
  obs <- clust_spec |> fit()
  terms <- unique(boot_fit$term)

  # default: every coefficient
  ci_all <- get_confidence_interval(
    boot_fit,
    type = "bca",
    point_estimate = obs
  )
  expect_equal(ci_all$term, terms)
  expect_true(all(is.finite(ci_all$lower_ci), is.finite(ci_all$upper_ci)))

  # single term by name
  ci_age <- get_confidence_interval(
    boot_fit,
    type = "bca",
    point_estimate = obs,
    contrast = "age"
  )
  expect_equal(ci_age$term, "age")
  expect_equal(
    ci_age$lower_ci,
    ci_all$lower_ci[ci_all$term == "age"]
  )

  # named numeric contrast
  w <- setNames(rep(0, length(terms)), terms)
  w[["age"]] <- 2
  college_term <- setdiff(terms, c("intercept", "age"))
  w[[college_term]] <- -1

  ci_c <- get_confidence_interval(
    boot_fit,
    type = "bca",
    point_estimate = obs,
    contrast = c(age = 2, rlang::set_names(-1, college_term))
  )
  expect_equal(ci_c$term, "contrast")

  wide <- tidyr::pivot_wider(
    boot_fit,
    names_from = term,
    values_from = estimate,
    id_cols = replicate
  )
  manual_boot <- as.matrix(wide[, terms]) %*% w
  expect_equal(
    attr(ci_c, "bootstrap_estimates")$contrast,
    as.numeric(manual_boot)
  )

  # unnamed numeric contrast of the right length gives the same interval
  ci_u <- get_confidence_interval(
    boot_fit,
    type = "bca",
    point_estimate = obs,
    contrast = unname(w)
  )
  expect_equal(ci_u$lower_ci, ci_c$lower_ci)
  expect_equal(ci_u$upper_ci, ci_c$upper_ci)
})

test_that("BCa validates its arguments", {
  set.seed(5)
  obs <- clust_spec |> fit()
  boot_fit <- clust_spec |>
    generate(reps = 15, clusters = customer) |>
    fit()
  row_fit <- clust_spec |>
    generate(reps = 15, type = "bootstrap") |>
    fit()

  # BCa is not defined for ordinary row-wise bootstrap objects
  expect_snapshot(
    error = TRUE,
    get_confidence_interval(row_fit, type = "bca", point_estimate = obs)
  )

  # point estimate is required
  expect_snapshot(
    error = TRUE,
    get_confidence_interval(boot_fit, type = "bca")
  )

  # unknown term name
  expect_snapshot(
    error = TRUE,
    get_confidence_interval(
      boot_fit,
      type = "bca",
      point_estimate = obs,
      contrast = "nope"
    )
  )

  # wrong-length unnamed contrast
  expect_snapshot(
    error = TRUE,
    get_confidence_interval(
      boot_fit,
      type = "bca",
      point_estimate = obs,
      contrast = c(1, 2)
    )
  )

  # contrast is ignored (with a warning) for non-BCa interval types
  expect_warning(
    get_confidence_interval(
      row_fit,
      type = "percentile",
      point_estimate = obs,
      contrast = "age"
    ),
    "contrast"
  )
})

test_that("non-estimable replicates and jackknife fits are not silently dropped", {
  set.seed(6)
  # one cluster is the sole source of the "rare" factor level
  rare <- do.call(rbind, lapply(1:15, function(id) {
    ni <- sample(2:4, 1)
    g <- if (id == 1) {
      "rare"
    } else {
      sample(c("a", "b"), ni, replace = TRUE)
    }
    data.frame(
      customer = paste0("c", id),
      y = rnorm(ni) + ifelse(g == "rare", 5, 0),
      g = factor(g, levels = c("a", "b", "rare"))
    )
  }))

  sp <- rare |> specify(y ~ g, clusters = customer)

  expect_warning(
    boot_fit <- sp |>
      generate(reps = 30, clusters = customer) |>
      fit(),
    "did not yield finite"
  )

  # failed replicates are retained as NA, not removed
  expect_true(any(!is.finite(boot_fit$estimate)))
  expect_true(any(!is.na(attr(boot_fit, "fit_status")$error)))

  # jackknife failure is recorded per cluster
  jack <- attr(boot_fit, "jackknife")
  expect_true(any(!is.finite(jack$estimate)))
  expect_true(any(!is.na(jack$error)))

  # BCa must error rather than silently subsetting replicates
  obs <- sp |> fit()
  expect_snapshot(
    error = TRUE,
    get_confidence_interval(boot_fit, type = "bca", point_estimate = obs)
  )
})

test_that("BCa errors when the bias correction is at the distribution edge", {
  bca_one <- infer:::bca_one

  expect_snapshot(
    error = TRUE,
    bca_one(
      theta_obs = 0,
      boot = runif(50, 1, 2),
      loo = rnorm(50),
      level = 0.95,
      label = "g"
    )
  )

  # constant jackknife values -> acceleration of zero with a note
  res <- bca_one(
    theta_obs = 1,
    boot = rnorm(50, 1, 0.1),
    loo = rep(2, 50),
    level = 0.95,
    label = "k"
  )
  expect_equal(res$acceleration, 0)
  expect_false(is.na(res$note))
})

# Existing modes remain available -----------------------------------------

test_that("existing sampling and interval modes are unchanged", {
  set.seed(7)

  # row-wise bootstrap with calculate() still works alongside cluster columns
  calc_dist <- clust_dat |>
    specify(response = hours, clusters = customer) |>
    generate(reps = 20, type = "bootstrap") |>
    calculate(stat = "mean")
  expect_equal(nrow(calc_dist), 20)
  expect_equal(
    get_confidence_interval(calc_dist, type = "percentile") |> nrow(),
    1
  )

  # percentile, se, and bias-corrected intervals on fitted cluster objects
  obs <- clust_spec |> fit()
  boot_fit <- clust_spec |>
    generate(reps = 20, clusters = customer) |>
    fit()

  for (ty in c("percentile", "se", "bias-corrected")) {
    ci <- get_confidence_interval(
      boot_fit,
      type = ty,
      point_estimate = obs
    )
    expect_equal(names(ci), c("term", "lower_ci", "upper_ci"))
    expect_equal(nrow(ci), length(unique(boot_fit$term)))
  }
})
