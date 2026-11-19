# specify() validates the clusters argument

    Code
      specify(clust_dat, hours ~ age, clusters = "customer")
    Condition
      Error in `specify()`:
      ! The `clusters` argument should be an unquoted column name (not a string in quotation marks).

---

    Code
      specify(clust_dat, hours ~ age, clusters = not_a_col)
    Condition
      Error in `specify()`:
      ! The cluster variable `not_a_col` cannot be found in this dataframe.

# generate() validates cluster arguments

    Code
      generate(clust_spec, reps = 5, type = "bootstrap", clusters = not_a_col)
    Condition
      Error in `generate()`:
      ! The cluster variable `not_a_col` cannot be found in this dataframe.
      i Cluster columns must be supplied to `specify()` via its `clusters` argument so that they are retained in the pipeline.

---

    Code
      generate(few, reps = 5, clusters = customer)
    Condition
      Error in `bootstrap()`:
      ! Cluster bootstrap requires at least three clusters.
      i The data contain 2 clusters in `customer`; leave-one-cluster-out acceleration is not defined for fewer than three clusters.

# BCa validates its arguments

    Code
      get_confidence_interval(row_fit, type = "bca", point_estimate = obs)
    Condition
      Error in `get_confidence_interval()`:
      ! `type = "bca"` requires a fitted cluster bootstrap distribution.
      i Generate one with `generate(..., clusters = <cluster col>)` followed by `fit()`.
      i For ordinary row-wise bootstrap objects, use "percentile", "se", or "bias-corrected"; a row-wise jackknife acceleration is intentionally not provided.

---

    Code
      get_confidence_interval(boot_fit, type = "bca")
    Condition
      Error in `get_confidence_interval()`:
      ! A numeric value needs to be given for `point_estimate` for `type` "bca". Pass the observed `fit()` output.

---

    Code
      get_confidence_interval(boot_fit, type = "bca", point_estimate = obs, contrast = "nope")
    Condition
      Error in `get_confidence_interval()`:
      ! All term names supplied to `contrast` must be coefficients of the fitted model.
      x Unknown term: "nope".
      i Available terms: "intercept", "age", and "collegeno degree".

---

    Code
      get_confidence_interval(boot_fit, type = "bca", point_estimate = obs, contrast = c(
        1, 2))
    Condition
      Error in `get_confidence_interval()`:
      ! An unnamed numeric `contrast` must have length equal to the number of coefficients (3), not 2.

# non-estimable replicates and jackknife fits are not silently dropped

    Code
      get_confidence_interval(boot_fit, type = "bca", point_estimate = obs)
    Condition
      Error in `get_confidence_interval()`:
      ! BCa interval for "intercept" is not estimable: 16 bootstrap replicate estimates are non-finite.
      i No replicates were removed; inspect `attr(<fit>, "fit_status")` and the bootstrap estimates.

# BCa errors when the bias correction is at the distribution edge

    Code
      bca_one(theta_obs = 0, boot = runif(50, 1, 2), loo = rnorm(50), level = 0.95,
      label = "g")
    Condition
      Error:
      ! BCa interval for "g" is not estimable: the bias correction is infinite.
      i 0% of bootstrap estimates fall below the observed estimate; the BCa adjustment lies outside the range of the bootstrap distribution.

