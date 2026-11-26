## Study 5: local versus global spacing for seeded NOT (Theorem 4.2).
## T=1000, N=20, sigma=1. A strong pair of breaks at 200 and 350 (Delta=8 each) fixes d_min=150 in both layouts.
## A weak break with jump Delta_w sits either isolated at 700 (local spacing 300) or crowded at 500 (local spacing 150).
## The global condition (Delta_j d_min) treats the two layouts alike; the local condition (Delta_j d_j) favours the isolated one.
source("common.R"); args <- commandArgs(TRUE); mode <- if (length(args)) args[1] else "run"
set.seed(501); T <- 1000; N <- 20; fam <- seeded_family(T); tol <- 15
jump <- function(Delta) sqrt(Delta / N) * sample(c(-1, 1), N, replace = TRUE)
make_mean <- function(dates, Deltas) {
  M <- matrix(0, T, N); level <- rep(0, N)
  ext <- c(0, dates, T)
  for (j in seq_along(ext)[-1]) {
    if (j > 2) level <- level + jump(Deltas[j - 2])
    rows <- (ext[j - 1] + 1):ext[j]; M[rows, ] <- matrix(level, length(rows), N, byrow = TRUE)
  }
  M
}
if (mode == "time") {
  t0 <- proc.time(); for (r in 1:5) { Y <- matrix(rnorm(T * N), T); sc <- seeded_scan(cumfull(Y), fam, N) }
  cat("seconds per replication:", round((proc.time() - t0)[3] / 5, 3), " intervals:", length(fam$a), "\n"); quit()
}
## threshold: 95% null quantile of the maximal seeded excess, from 2000 null panels
R0 <- 2000; nm <- numeric(R0)
for (r in 1:R0) nm[r] <- max(seeded_scan(cumfull(matrix(rnorm(T * N), T)), fam, N)$M)
tau <- unname(quantile(nm, 0.95)); cat("threshold:", round(tau, 2), "\n")
layouts <- list(isolated = c(200, 350, 700), crowded = c(200, 350, 500))
grid <- c(0.4, 0.6, 0.8, 1.0, 1.4, 2.0); R <- 200; out <- NULL
for (lay in names(layouts)) for (Dw in grid) {
  dates <- layouts[[lay]]; Mn <- make_mean(dates, c(8, 8, Dw))
  weak <- all3 <- nd <- numeric(R); errw <- numeric(0)
  for (r in 1:R) {
    Y <- Mn + matrix(rnorm(T * N), T); f <- seeded_not(seeded_scan(cumfull(Y), fam, N), fam, tau, "intersect")
    hit <- sapply(dates, function(b) any(abs(f - b) <= tol))
    weak[r] <- hit[3]; all3[r] <- length(f) == 3 && all(hit); nd[r] <- length(f)
    if (hit[3]) errw <- c(errw, min(abs(f - dates[3])))
  }
  dw <- if (lay == "isolated") 300 else 150
  out <- rbind(out, data.frame(layout = lay, Delta_w = Dw, d_w = dw, local = Dw * dw, global = Dw * 150, R = R,
    weak = mean(weak), all3 = mean(all3), dates = mean(nd), mae_weak = if (length(errw)) mean(errw) else NA))
  cat(sprintf("%-8s Dw=%.1f  local=%5.0f global=%5.0f  weak=%.3f all3=%.3f dates=%.2f\n", lay, Dw, Dw * dw, Dw * 150, mean(weak), mean(all3), mean(nd)))
}
## null false-positive rate at this threshold
fp <- mean(replicate(1000, length(seeded_not(seeded_scan(cumfull(matrix(rnorm(T * N), T)), fam, N), fam, tau, "intersect")) > 0))
cat("null false-positive rate:", fp, "\n")
write_csv(cbind(out, tau = tau, null_fp = fp), "s5_local.csv")
