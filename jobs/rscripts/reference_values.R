# Tier 4: reference values from the R implementations our Python routines were transcribed
# from, or transcribed to agree with.
#
# WHY THIS EXISTS. Tier 2 pins arithmetic against a published table; Tier 3 measures whether
# a p-value holds its level. Neither catches a routine that is self-consistently wrong - one
# that computes SOME well-behaved statistic that is not the statistic its source defines.
# Only an independent implementation catches that, and for these estimators the independent
# implementations are in R, several of them written by the authors of the papers.
#
# HOW IT IS USED. This script writes two committed CSV fixtures. The test suite reads the
# FIXTURES and never invokes R, so `pytest` stays runnable on a machine with no R at all.
# Re-run this script to regenerate them when a package version changes:
#
#     Rscript jobs/rscripts/reference_values.R
#
# The script writes the INPUT DATA into the fixture alongside the reference values. That is
# deliberate: matching R's RNG stream from Python is fragile and would make a spurious
# mismatch look like a statistical disagreement. Python reads the exact numbers R saw.
#
# CSV rather than JSON because `jsonlite` is not installed here and a long-format CSV diffs
# cleanly under review, which a nested JSON blob does not.

suppressPackageStartupMessages({
  library(XICOR)
  library(energy)
  library(randtests)
  library(copula)
})

# TRACKED, unlike tests/ which this repo gitignores. The fixture must survive a
# fresh clone or the suite would silently require R - the one thing the plan forbids.
out_dir <- file.path("jobs", "reference")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

inputs <- list()
values <- list()

add_value <- function(case, quantity, value) {
  values[[length(values) + 1]] <<- data.frame(
    case = case, quantity = quantity, value = as.numeric(value),
    stringsAsFactors = FALSE
  )
}

add_input <- function(case, x, y) {
  inputs[[length(inputs) + 1]] <<- data.frame(
    case = case, i = seq_along(x), x = as.numeric(x), y = as.numeric(y),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------- provenance of the run
add_value("meta", "r_major", as.numeric(R.version$major))
add_value("meta", "r_minor", as.numeric(strsplit(R.version$minor, "[.]")[[1]][1]))
for (pkg in c("XICOR", "energy", "randtests", "copula")) {
  v <- as.character(packageVersion(pkg))
  parts <- as.numeric(strsplit(v, "[.-]")[[1]])
  add_value("meta", paste0(pkg, "_major"), parts[1])
  add_value("meta", paste0(pkg, "_minor"), parts[2])
}

# ------------------------------------------------------------- Chatterjee's xi, TIE-FREE
# The reduction `1 - 3*sum|dr|/(n^2-1)` is only valid here. This case pins that our
# tie-corrected xi (Chatterjee 2021) still agrees with the reference implementation where
# the two forms coincide.
set.seed(101)
n <- 60
x <- rnorm(n)
y <- x^2 + rnorm(n, sd = 0.3)   # non-monotone on purpose: xi sees it, Spearman does not
add_input("xi_tie_free", x, y)
add_value("xi_tie_free", "xicor", XICOR::xicor(x, y, ties = TRUE))
add_value("xi_tie_free", "xicor_ties_false", XICOR::xicor(x, y, ties = FALSE))
add_value("xi_tie_free", "spearman", cor(x, y, method = "spearman"))
add_value("xi_tie_free", "dcor", energy::dcor(x, y))

# ------------------------------------------------------------------ Chatterjee's xi, TIED
# THE CASE THAT MATTERS. Our estimator uses the tie-corrected form, not the tie-free
# reduction, and no shipped window exercises it - all 309 are tie-free. This is the only
# external evidence that the tie-corrected path is right. `ties = TRUE` is XICOR's own
# tie-corrected path.
set.seed(202)
n <- 60
x <- round(rnorm(n), 1)              # ties on x
y <- round(x + rnorm(n, sd = 0.5))   # heavy ties on y
add_input("xi_tied", x, y)
# `xicor` is RANDOM when x has ties: Chatterjee breaks them uniformly at random, so repeated
# calls on identical data give different answers - measured, 7 distinct values in 8 calls,
# spanning 0.425 to 0.563. A single draw is therefore NOT a reference value, and pinning it
# would pin R's RNG state rather than the estimator. What IS well defined is the
# DISTRIBUTION over tie-breaks, so the fixture carries its summary and our deterministic
# stable-sort value is checked for membership in it.
set.seed(2024)
xi_draws <- replicate(2000, XICOR::xicor(x, y, ties = TRUE))
add_value("xi_tied", "xicor_draws_n", length(xi_draws))
add_value("xi_tied", "xicor_mean", mean(xi_draws))
add_value("xi_tied", "xicor_sd", sd(xi_draws))
add_value("xi_tied", "xicor_q001", unname(quantile(xi_draws, 0.001)))
add_value("xi_tied", "xicor_q01", unname(quantile(xi_draws, 0.01)))
add_value("xi_tied", "xicor_q99", unname(quantile(xi_draws, 0.99)))
add_value("xi_tied", "xicor_q999", unname(quantile(xi_draws, 0.999)))
add_value("xi_tied", "xicor_min", min(xi_draws))
add_value("xi_tied", "xicor_max", max(xi_draws))
add_value("xi_tied", "xicor_ties_on_x", sum(duplicated(x)))
add_value("xi_tied", "xicor_ties_on_y", sum(duplicated(y)))
add_value("xi_tied", "spearman", cor(x, y, method = "spearman"))
add_value("xi_tied", "dcor", energy::dcor(x, y))

# A second tied case with ties on y ONLY, which is the case Chatterjee's Thm 2.1 excludes
# and therefore the case our asymptotic null is not entitled to.
set.seed(303)
n <- 45
x <- rnorm(n)
y <- as.numeric(cut(x + rnorm(n, sd = 0.8), breaks = 5))
add_input("xi_tied_y_only", x, y)
# x is continuous here, so there are no x-ties to break and `xicor` IS deterministic - a
# single value is a genuine reference. Confirmed rather than assumed: the draws below must
# all be identical, and the fixture records that they are.
set.seed(3030)
y_only_draws <- replicate(50, XICOR::xicor(x, y, ties = TRUE))
add_value("xi_tied_y_only", "xicor", y_only_draws[1])
add_value("xi_tied_y_only", "xicor_distinct_draws", length(unique(y_only_draws)))
add_value("xi_tied_y_only", "xicor_ties_on_x", sum(duplicated(x)))
add_value("xi_tied_y_only", "spearman", cor(x, y, method = "spearman"))

# --------------------------------------------------------------- distance correlation
# Szekely-Rizzo-Bakirov. `energy` is the authors' own package, so this is as close to a
# definitional reference as exists.
set.seed(404)
n <- 50
x <- runif(n, -2, 2)
y <- sin(3 * x) + rnorm(n, sd = 0.2)   # dependent but rank-uncorrelated by construction
add_input("dcor_nonmonotone", x, y)
add_value("dcor_nonmonotone", "dcor", energy::dcor(x, y))
add_value("dcor_nonmonotone", "spearman", cor(x, y, method = "spearman"))
add_value("dcor_nonmonotone", "xicor", XICOR::xicor(x, y, ties = TRUE))

# ------------------------------------------------------- a duration series for the checks
# One series, several instruments, so a disagreement can be localised to an instrument
# rather than to the data.
set.seed(505)
n <- 80
durations <- rexp(n, rate = 1 / 30)
add_input("durations_iid", durations, rep(NA_real_, n))
add_value("durations_iid", "lag1_rank_autocorr",
          cor(rank(durations)[-n], rank(durations)[-1], method = "pearson"))
bt <- randtests::bartels.rank.test(durations)
add_value("durations_iid", "bartels_statistic", unname(bt$statistic))
add_value("durations_iid", "bartels_p_value", bt$p.value)

# A trended series, so the same instruments are pinned somewhere they should react.
set.seed(606)
n <- 80
durations_trend <- rexp(n, rate = 1 / 30) * seq(0.4, 2.0, length.out = n)
add_input("durations_trend", durations_trend, rep(NA_real_, n))
add_value("durations_trend", "lag1_rank_autocorr",
          cor(rank(durations_trend)[-n], rank(durations_trend)[-1], method = "pearson"))
bt2 <- randtests::bartels.rank.test(durations_trend)
add_value("durations_trend", "bartels_statistic", unname(bt2$statistic))
add_value("durations_trend", "bartels_p_value", bt2$p.value)

# ------------------------------------------------------------------- C3, copula::serialIndepTest
# Pinned with an explicit simulation seed. The statistic is seed-invariant; the p-value is
# NOT, because the null is simulated - the same distinction the C3 module docstring records.
set.seed(707)
sim <- copula::serialIndepTestSim(n = 80, lag.max = 5, N = 1000, verbose = FALSE)
sit <- copula::serialIndepTest(durations, sim)
add_value("durations_iid", "serial_indep_global_statistic", sit$global.statistic)
add_value("durations_iid", "serial_indep_global_p_value", sit$global.statistic.pvalue)
add_value("durations_iid", "serial_indep_sim_seed", 707)
add_value("durations_iid", "serial_indep_N", 1000)
add_value("durations_iid", "serial_indep_lag_max", 5)

# ------------------------------------------------------------------------------ write out
inputs_df <- do.call(rbind, inputs)
values_df <- do.call(rbind, values)

write.csv(inputs_df, file.path(out_dir, "r_reference_inputs.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(values_df, file.path(out_dir, "r_reference_values.csv"),
          row.names = FALSE, quote = FALSE)

cat(sprintf("wrote %d input rows and %d reference values\n",
            nrow(inputs_df), nrow(values_df)))
