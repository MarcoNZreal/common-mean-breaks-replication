## Study 2: power of the pooled maximum-CUSUM test against rho = Delta d / (sigma^2 sqrt N). Break at 100 in (0,200], lambda=0.1.
source("common.R"); set.seed(103)
n <- 200; b <- 100; d <- 100; out <- NULL
stat <- function(X) max(pooled(X, 0.1)$v)
for (N in c(25, 100, 400)) {
  null <- replicate(2000, stat(matrix(rnorm(n * N), n))); crit <- unname(quantile(null, 0.95))
  for (rho in c(0.5, 1, 2, 4, 6, 8, 12)) {
    Delta <- rho * sqrt(N) / d
    pw <- mean(replicate(500, { X <- matrix(rnorm(n * N), n); X[(b + 1):n, ] <- X[(b + 1):n, ] + sqrt(Delta / N); stat(X) > crit }))
    out <- rbind(out, data.frame(N = N, rho = rho, Delta = Delta, power = pw, R = 500, crit = crit, null_mean = mean(null)))
  }
  cat("N =", N, "done\n")
}
write_csv(out, "s2_detect.csv")
