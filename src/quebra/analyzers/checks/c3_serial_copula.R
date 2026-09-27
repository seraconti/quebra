#!/usr/bin/env Rscript
# C3 - copula::serialIndepTest, called from checks/c3_serial_copula.py over a file bridge.
#
# Usage: Rscript c3_serial_copula.R <in.csv> <out.csv> <lag.max> <n.sim> <seed>
#
# Reads a one-column CSV of durations, runs the empirical-copula serial independence test,
# and writes a one-row CSV with `statistic` and `p_value`.
#
# It exits non-zero on any surprise rather than writing a plausible number, because the
# Python side turns a non-zero exit into a traceback and a written file into a reported
# p-value.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("expected 5 arguments: in.csv out.csv lag.max n.sim seed", call. = FALSE)
}

in_path  <- args[[1L]]
out_path <- args[[2L]]
lag_max  <- as.integer(args[[3L]])
n_sim    <- as.integer(args[[4L]])
seed     <- as.integer(args[[5L]])

if (!requireNamespace("copula", quietly = TRUE)) {
  stop("the 'copula' package is required for C3", call. = FALSE)
}

frame <- utils::read.csv(in_path)
if (!"duration_s" %in% names(frame)) {
  stop("input CSV must have a 'duration_s' column", call. = FALSE)
}
x <- as.numeric(frame$duration_s)
if (anyNA(x) || any(!is.finite(x))) {
  stop("input durations must all be finite", call. = FALSE)
}
if (length(x) <= lag_max + 1L) {
  stop("series is too short for the requested lag.max", call. = FALSE)
}

# Pinned so a rerun with the same seed reproduces the same simulated null. The null
# depends only on (n, lag.max, seed).
set.seed(seed)

sim <- copula::serialIndepTestSim(
  n       = length(x),
  lag.max = lag_max,
  N       = n_sim,
  verbose = FALSE
)
result <- copula::serialIndepTest(x, d = sim)

statistic <- as.numeric(result$global.statistic)
p_value   <- as.numeric(result$global.statistic.pvalue)
if (length(statistic) != 1L || length(p_value) != 1L) {
  stop("serialIndepTest returned an unexpected shape", call. = FALSE)
}

utils::write.csv(
  data.frame(statistic = statistic, p_value = p_value),
  file = out_path,
  row.names = FALSE
)
