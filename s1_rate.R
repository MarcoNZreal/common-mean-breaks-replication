## Study 1a/1b: date error of the pooled squared CUSUM and the known-direction oracle. One break at b=250 in (0,500], lambda=0.1, sigma=1.
source("common.R"); set.seed(101)
n <- 500; b <- 250; d <- 250; lam <- 0.1
one <- function(N, Delta) {
  X <- matrix(rnorm(n * N), n); X[(b + 1):n, ] <- X[(b + 1):n, ] + sqrt(Delta / N)
  c(abs(pooled(X, lam)$k - b), abs(pooled(rowSums(X) / sqrt(N), lam)$k - b))
}
design <- rbind(data.frame(panel = "a", N = c(1, 10, 100, 400, 1000, 2000), Delta = 0.5),
                data.frame(panel = "b", N = 400, Delta = c(0.5, 1, 2, 4, 8)))
R <- 500; out <- NULL
for (i in seq_len(nrow(design))) {
  N <- design$N[i]; D <- design$Delta[i]; r <- t(replicate(R, one(N, D)))
  out <- rbind(out, data.frame(design[i, ], R = R, q = N / (D * d), r_or = 1 / D, r_dir = N / (D^2 * d),
    oracle_mean = mean(r[, 2]), oracle_se = sd(r[, 2]) / sqrt(R), pooled_mean = mean(r[, 1]), pooled_se = sd(r[, 1]) / sqrt(R),
    oracle_med = median(r[, 2]), pooled_med = median(r[, 1])))
  cat(sprintf("%s N=%d Delta=%.2f pooled=%.2f oracle=%.2f\n", design$panel[i], N, D, mean(r[, 1]), mean(r[, 2])))
}
write_csv(out, "s1_rate.csv")
