## Study 6: exact recovery (Corollary 3.12). n=200, break at b=75 (d=75), trimming 0.1.
## Oracle-exponent design: N=10, gamma = Delta (2000 replications). Direction-exponent design: N=8000, Delta chosen so that
## Delta^2 d/N = gamma, between 3.6 and 15 times below Delta (1000 replications).
## Reported: frequency of a missed date, P(bhat != b), with 95% Wilson intervals; in the direction design also the known-direction oracle.
source("common.R"); n <- 200; b <- 75; d <- 75
args <- commandArgs(TRUE); dir.create("parts", showWarnings = FALSE)
wilson <- function(k, R) { p <- k / R; z <- 1.96; c0 <- p + z^2 / (2 * R); h <- z * sqrt(p * (1 - p) / R + z^2 / (4 * R^2)); c((c0 - h), (c0 + h)) / (1 + z^2 / R) }
## run one cell per call: Rscript s6_exact.R <design> <gamma>; then Rscript s6_exact.R combine
if (args[1] == "combine") { f <- Sys.glob("parts/s6_*.csv"); x <- do.call(rbind, lapply(f, read.csv)); x <- x[order(x$design != "oracle-exponent", x$gamma), ]; write_csv(x, "s6_exact.csv"); quit() }
for (design in args[1]) for (g in as.numeric(args[2])) {
  set.seed(round(1000 * g) + if (design == "oracle-exponent") 1 else 2); out <- NULL
  N <- if (design == "oracle-exponent") 10 else 8000; R <- if (design == "oracle-exponent") 2000 else 1000
  Delta <- if (design == "oracle-exponent") g else sqrt(g * N / d)
  miss <- miss_or <- 0
  for (r in 1:R) {
    delta <- sample(c(-1, 1), N, TRUE) * sqrt(Delta / N)
    Y <- matrix(rnorm(n * N), n); Y[(b + 1):n, ] <- sweep(Y[(b + 1):n, , drop = FALSE], 2, delta, "+")
    miss <- miss + (pooled(Y, 0.1)$k != b)
    if (design == "direction-exponent") miss_or <- miss_or + (pooled(Y %*% (delta / sqrt(Delta)), 0.1)$k != b)
  }
  ci <- wilson(miss, R); ci_or <- wilson(miss_or, R)
  out <- rbind(out, data.frame(design = design, N = N, Delta = Delta, gamma = g, dir_exponent = Delta^2 * d / N, R = R,
    miss = miss / R, lo = ci[1], hi = ci[2], miss_oracle = if (design == "direction-exponent") miss_or / R else NA,
    lo_or = if (design == "direction-exponent") ci_or[1] else NA, hi_or = if (design == "direction-exponent") ci_or[2] else NA))
  cat(sprintf("%-19s gamma=%.1f Delta=%.2f  miss=%.4f [%.4f, %.4f]%s\n", design, g, Delta, miss / R, ci[1], ci[2],
      if (design == "direction-exponent") sprintf("  oracle=%.4f", miss_or / R) else ""))
  write_csv(out, sprintf("parts/s6_%s_%g.csv", design, g))
}
