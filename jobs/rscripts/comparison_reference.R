# Reference values for the comparisons across records, from R's `survival` and `metafor`.
#
# The suite reads the two CSVs this script writes and never invokes R, so pytest stays
# runnable with no R installed. Re-run when a case is added or a package changes version:
#
#     Rscript jobs/rscripts/comparison_reference.R
#
# Inputs are written beside the values, so Python reads the exact numbers R saw. Every case
# is a set of records, each a list of right-censored ages; ages are rounded to 0.1 so ties
# are exact floats on both sides.
#
# Conventions pinned here, each one a choice that changes a number:
#   - Two ages tie only when they are the same float, as in `event_table.risk_table`: survfit
#     runs with timefix = FALSE. survdiff 3.8.6 forwards `timefix` into model.frame and
#     fails on it, so it runs with the default, and the script asserts that `aeqSurv` (what
#     timefix = TRUE applies) leaves every age of the case unchanged.
#   - The comparison tau is the smallest, over the records of a case, of each record's largest
#     age; the record that sets it is the first one attaining it.
#   - survdiff(rho = 0): the log-rank test; its p is the chi-square(1) upper tail.
#   - summary(fit, rmean = tau): each record's restricted mean and its standard error. The
#     difference of two records is a - b with SE sqrt(se_a^2 + se_b^2).
#   - rma(method = "DL"): DerSimonian-Laird on the per-record restricted means, with the
#     squared standard errors as known variances; normal CI on the pooled mean.
#   - predict(predtype = "Riley"): the prediction interval on t with k - 2 df.
#   - confint(): the Q-profile interval for tau^2, with the search widened and tightened
#     (the default stops at max(100, 10 tau^2) with tolerance 1.2e-4).
#
# Each pooling case asserts the branch it exists for, so a re-run that drifts off it fails.

suppressPackageStartupMessages({
  library(survival)
  library(metafor)
})

out_dir <- file.path("jobs", "reference")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

inputs <- list()
values <- list()

add_value <- function(case, quantity, i, value) {
  values[[length(values) + 1]] <<- data.frame(
    case = case, quantity = quantity, i = i, value = as.numeric(value),
    stringsAsFactors = FALSE
  )
}

add_records <- function(case, records) {
  for (name in names(records)) {
    r <- records[[name]]
    inputs[[length(inputs) + 1]] <<- data.frame(
      case = case, record = name, i = seq_along(r$age), age = r$age,
      event = as.integer(r$event), stringsAsFactors = FALSE
    )
  }
}

gen <- function(n, mean, cens) {
  list(age = round(rexp(n, rate = 1 / mean), 1), event = as.integer(runif(n) > cens))
}

comparison_tau <- function(records) {
  largest <- sapply(records, function(r) max(r$age))
  list(tau = min(largest), set_by = which.min(largest))
}

restricted_mean <- function(r, tau) {
  fit <- survfit(Surv(r$age, r$event) ~ 1, timefix = FALSE)
  tab <- summary(fit, rmean = tau)$table
  c(tab[["rmean"]], tab[["se(rmean)"]])
}

pair_case <- function(case, a, b) {
  records <- list(a = a, b = b)
  add_records(case, records)
  ct <- comparison_tau(records)
  add_value(case, "tau", 1, ct$tau)
  add_value(case, "tau_set_by_index", 1, ct$set_by)
  d <- data.frame(
    record = rep(c("a", "b"), c(length(a$age), length(b$age))),
    age = c(a$age, b$age), event = c(a$event, b$event)
  )
  stopifnot(identical(as.numeric(aeqSurv(Surv(d$age, d$event))[, 1]), as.numeric(d$age)))
  lr <- survdiff(Surv(age, event) ~ record, data = d, rho = 0)
  add_value(case, "logrank_chisq", 1, lr$chisq)
  add_value(case, "logrank_p", 1, pchisq(lr$chisq, df = 1, lower.tail = FALSE))
  add_value(case, "logrank_o_minus_e_a", 1, lr$obs[1] - lr$exp[1])
  ra <- restricted_mean(a, ct$tau)
  rb <- restricted_mean(b, ct$tau)
  add_value(case, "rmean_a", 1, ra[1])
  add_value(case, "rmean_se_a", 1, ra[2])
  add_value(case, "rmean_b", 1, rb[1])
  add_value(case, "rmean_se_b", 1, rb[2])
  add_value(case, "rmean_diff", 1, ra[1] - rb[1])
  add_value(case, "rmean_diff_se", 1, sqrt(ra[2]^2 + rb[2]^2))
}

pool_case <- function(case, records, branch) {
  names(records) <- paste0("r", seq_along(records))
  add_records(case, records)
  ct <- comparison_tau(records)
  add_value(case, "tau", 1, ct$tau)
  add_value(case, "tau_set_by_index", 1, ct$set_by)
  m <- t(sapply(records, restricted_mean, tau = ct$tau))
  yi <- m[, 1]
  vi <- m[, 2]^2
  stopifnot(all(vi > 0))
  k <- length(yi)
  add_value(case, "yi", seq_len(k), yi)
  add_value(case, "vi", seq_len(k), vi)
  res <- rma(yi, vi, method = "DL")
  add_value(case, "dl_mu", 1, res$beta[1])
  add_value(case, "dl_se", 1, res$se)
  add_value(case, "dl_ci_lb", 1, res$ci.lb)
  add_value(case, "dl_ci_ub", 1, res$ci.ub)
  add_value(case, "dl_tau2", 1, res$tau2)
  add_value(case, "dl_q", 1, res$QE)
  add_value(case, "dl_i2_percent", 1, res$I2)
  pr <- predict(res, predtype = "Riley")
  add_value(case, "pi_lb", 1, pr$pi.lb)
  add_value(case, "pi_ub", 1, pr$pi.ub)
  ci <- confint(res, control = list(tol = 1e-12, tau2.max = 1e6))$random["tau^2", ]
  add_value(case, "qp_tau2_lb", 1, ci[["ci.lb"]])
  add_value(case, "qp_tau2_ub", 1, ci[["ci.ub"]])
  ok <- switch(branch,
    both_bounds_zero = res$tau2 == 0 && ci[["ci.lb"]] == 0 && ci[["ci.ub"]] == 0,
    truncated_lower_zero = res$tau2 == 0 && ci[["ci.lb"]] == 0 && ci[["ci.ub"]] > 0,
    positive_lower_zero = res$tau2 > 0 && ci[["ci.lb"]] == 0,
    positive_lower_positive = res$tau2 > 0 && ci[["ci.lb"]] > 0
  )
  stopifnot(isTRUE(ok))
}

r_version <- c(as.integer(R.version$major), as.integer(strsplit(R.version$minor, ".",
  fixed = TRUE)[[1]]))
add_value("meta", "r_major", 1, r_version[1])
add_value("meta", "r_minor", 1, r_version[2])
add_value("meta", "r_patch", 1, r_version[3])
for (pkg in c("survival", "metafor")) {
  v <- unlist(packageVersion(pkg))
  add_value("meta", paste0(pkg, "_major"), 1, v[1])
  add_value("meta", paste0(pkg, "_minor"), 1, v[2])
  add_value("meta", paste0(pkg, "_patch"), 1, v[3])
}

# ---------------------------------------------------------------- pairs
# Deaths and censorings tied at the same ages within and across records; a sets tau with a
# censored largest age.
pair_case("pair_ties_across_records",
  list(age = c(1, 2, 2, 3, 4, 4, 5, 6), event = c(1, 1, 0, 1, 1, 0, 1, 0)),
  list(age = c(2, 2, 3, 4, 5, 5, 7, 8), event = c(1, 0, 1, 0, 1, 1, 1, 0)))
# Both records end at the same age: the first one sets tau.
pair_case("pair_equal_largest_ages",
  list(age = c(1, 3, 3, 5), event = c(1, 1, 0, 1)),
  list(age = c(2, 2, 4, 5), event = c(1, 1, 1, 0)))
# The second record sets tau with a death at its largest age.
pair_case("pair_tau_set_by_a_death",
  list(age = c(1, 3, 5, 6, 9), event = c(1, 1, 0, 1, 0)),
  list(age = c(1, 2, 4), event = c(1, 0, 1)))
# One record has no death: its restricted mean is tau with standard error 0.
pair_case("pair_one_record_without_deaths",
  list(age = c(2, 3, 5, 7), event = c(0, 0, 0, 0)),
  list(age = c(1, 2, 4, 6, 8), event = c(1, 1, 1, 0, 1)))
# A larger mixed pair, away from the edges.
set.seed(20261030)
pair_case("pair_mixed", gen(40, 12, 0.2), gen(35, 8, 0.25))

# ---------------------------------------------------------------- pooling
# One record copied three times: Q = 0, tau^2 = 0, and both Q-profile bounds are 0.
set.seed(20261031)
one <- gen(20, 10, 0.2)
pool_case("pool_identical_k3", list(one, one, one), "both_bounds_zero")
# One law, five records: Q below its df, tau^2 truncated at 0, upper bound above 0.
set.seed(2)
pool_case("pool_near_homogeneous_k5", lapply(rep(10, 5), gen, n = 25, cens = 0.2),
  "truncated_lower_zero")
# Three records: tau^2 above 0, but Q below the upper chi-square quantile.
set.seed(2)
pool_case("pool_moderate_k3", lapply(c(6, 10, 16), gen, n = 20, cens = 0.2),
  "positive_lower_zero")
# Six records from laws far apart: the Q-profile lower bound is above 0.
set.seed(1)
pool_case("pool_heterogeneous_k6", lapply(c(5, 8, 12, 18, 25, 40), gen, n = 30, cens = 0.2),
  "positive_lower_positive")

write.csv(do.call(rbind, inputs), file.path(out_dir, "r_comparison_inputs.csv"),
  row.names = FALSE)
write.csv(do.call(rbind, values), file.path(out_dir, "r_comparison_values.csv"),
  row.names = FALSE)
cat("wrote", file.path(out_dir, "r_comparison_inputs.csv"), "and r_comparison_values.csv\n")
