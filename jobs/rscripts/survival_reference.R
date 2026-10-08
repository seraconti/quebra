# Reference values for the survival estimators, from R's `survival` package.
#
# The suite reads the two CSVs this script writes and never invokes R, so pytest stays
# runnable with no R installed. Re-run when a case is added or `survival` changes version:
#
#     Rscript jobs/rscripts/survival_reference.R
#
# Inputs are written beside the values, so Python reads the exact numbers R saw.
#
# Conventions pinned here, each one a choice that changes a number:
#   - survfit(conf.type = "log-log"). R's default is "log"; the estimators use log-log.
#   - ctype = 1: the Nelson-Aalen cumulative hazard, with std.chaz = sqrt(sum d / n^2).
#   - quantile(): R's rule, the midpoint of a stretch where S sits exactly at 1 - p; the
#     interval is the same rule applied to the band's two curves.
#   - summary(fit, rmean = tau): restricted mean and its standard error.

suppressPackageStartupMessages(library(survival))

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

add_case <- function(case, age, event, tau) {
  inputs[[length(inputs) + 1]] <<- data.frame(
    case = case, i = seq_along(age), age = age, event = as.integer(event), tau = tau,
    stringsAsFactors = FALSE
  )
  fit <- survfit(Surv(age, event) ~ 1, conf.type = "log-log", ctype = 1)
  s <- summary(fit)
  n_steps <- length(s$time)
  if (n_steps > 0) {
    idx <- seq_len(n_steps)
    add_value(case, "km_time", idx, s$time)
    add_value(case, "km_surv", idx, s$surv)
    add_value(case, "km_lower", idx, s$lower)
    add_value(case, "km_upper", idx, s$upper)
    add_value(case, "na_cumhaz", idx, s$cumhaz)
    add_value(case, "na_std", idx, s$std.chaz)
  }
  # The probabilities are written out, so Python asks for the ones R was asked for.
  probs <- c(0.25, 0.5, 0.75)
  q <- quantile(fit, probs)
  add_value(case, "quantile_p", seq_along(probs), probs)
  add_value(case, "quantile", seq_along(probs), q$quantile)
  add_value(case, "quantile_lower", seq_along(probs), q$lower)
  add_value(case, "quantile_upper", seq_along(probs), q$upper)
  tab <- summary(fit, rmean = tau)$table
  add_value(case, "rmean", 1, tab[["rmean"]])
  add_value(case, "rmean_se", 1, tab[["se(rmean)"]])
}

# The full version, so a docstring citing it can be checked against the fixture.
survival_version <- unlist(packageVersion("survival"))
add_value("meta", "survival_major", 1, survival_version[1])
add_value("meta", "survival_minor", 1, survival_version[2])
add_value("meta", "survival_patch", 1, survival_version[3])

# Deaths and censorings tied at the same ages: the risk-set convention matters.
add_case("ties_and_censoring",
         c(1, 2, 2, 3, 4, 4, 5, 6, 7, 8), c(1, 1, 0, 1, 1, 0, 1, 0, 1, 0), 8)
# No censoring, S lands exactly on 0.75, 0.5 and 0.25: the quantiles are the sample ones.
add_case("flat_at_every_quartile", c(10, 20, 30, 40), c(1, 1, 1, 1), 40)
# S sits at 0.5 until the end of follow-up, which ends censored.
add_case("flat_to_censored_end", c(10, 20, 30, 40), c(1, 1, 0, 0), 40)
# A censoring falls inside the stretch where S sits at 0.5.
add_case("censoring_inside_flat", c(1, 2, 3, 4, 5, 6), c(1, 1, 1, 0, 1, 1), 6)
# The last two windows die together: n == d at the last age, which is tau.
add_case("risk_set_consumed_at_tau", c(1, 2, 3, 3), c(1, 1, 1, 1), 3)
# The same, with tau past the last age: the n == d term has no area left after it.
add_case("risk_set_consumed_before_tau", c(1, 2, 3, 3), c(1, 1, 1, 1), 5)
# tau far past the last age, which is censored: R extends the last value of S flat to tau.
add_case("tau_past_follow_up", c(1, 2), c(1, 0), 100)
# Duplicate timestamps give age-0 windows, one of them an observed death.
add_case("age_zero", c(0, 0, 1, 2, 3, 5), c(1, 0, 1, 1, 0, 1), 4)
# A larger mixed case, so the band and the restricted-mean SE are exercised away from edges.
set.seed(20261006)
age <- round(rexp(40, rate = 1 / 12), 1)
event <- as.integer(runif(40) > 0.2)
add_case("mixed_forty", age, event, 20)

# ---------------------------------------------------------------- Turnbull, interval2
# Surv(L, R, type = "interval2"): L == R is an exact length, R = NA is right-censored.
# survfit reports a time at the midpoint of each innermost interval (at the point itself
# for an exact length), and can report other times that carry zero mass (2 and 6 in
# tb_read_grid). Its EM has no exposed tolerance and can stop short of the maximum on a
# flat likelihood, so the tests pin the POSITIONS and compare log-likelihoods; they
# compare masses only where the estimate has a closed form.
tb_inputs <- list()
add_turnbull_case <- function(case, L, R) {
  tb_inputs[[length(tb_inputs) + 1]] <<- data.frame(
    case = case, i = seq_along(L), lo = L, hi = ifelse(is.na(R), Inf, R),
    stringsAsFactors = FALSE
  )
  fit <- survfit(Surv(L, R, type = "interval2") ~ 1)
  mass <- -diff(c(1, fit$surv))
  idx <- seq_along(fit$time)
  add_value(case, "tb_time", idx, fit$time)
  add_value(case, "tb_mass", idx, mass)
}
# Overlapping intervals, an exact length, a right-censored one.
add_turnbull_case("tb_mixed", c(0, 1, 2, 3, 2, 5), c(2, 3, 4, NA, 2, 7))
# Intervals of a regular read grid: [a - 1, a + 1] per death, [a, Inf) per censoring.
add_turnbull_case("tb_read_grid",
                  c(0, 1, 1, 2, 2, 3, 4, 5, 2, 6), c(2, 3, 3, 4, 4, 5, 6, 7, NA, NA))
# Disjoint intervals: the estimate is the count in each over n.
add_turnbull_case("tb_disjoint", c(0, 0, 3, 3, 3, 7), c(1, 1, 4, 4, 4, 9))

write.csv(do.call(rbind, tb_inputs), file.path(out_dir, "r_turnbull_inputs.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, inputs), file.path(out_dir, "r_survival_inputs.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, values), file.path(out_dir, "r_survival_values.csv"),
          row.names = FALSE)
cat("wrote", file.path(out_dir, "r_survival_inputs.csv"), "and r_survival_values.csv\n")
