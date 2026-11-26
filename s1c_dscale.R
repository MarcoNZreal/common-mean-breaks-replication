## Study 1c: out-of-sample test of the 1/d scaling. N=400, Delta=1, b=d=T/2, lambda=0.1. Predictions use constants fitted to Study 1a/1b.
source("common.R"); set.seed(102)
fit <- read.csv("s1_rate.csv"); use <- fit$pooled_mean < 40
cf <- coef(lm(pooled_mean ~ 0 + r_or + r_dir, data = fit[use, ]))
N <- 400; Delta <- 1; out <- NULL
for (dd in list(c(200, 400), c(400, 300), c(800, 200), c(1600, 150))) {
  d <- dd[1]; R <- dd[2]; n <- 2 * d; b <- d
  r <- t(replicate(R, { X <- matrix(rnorm(n * N), n); X[(b + 1):n, ] <- X[(b + 1):n, ] + sqrt(Delta / N)
                        c(abs(pooled(X, 0.1)$k - b), abs(pooled(rowSums(X) / sqrt(N), 0.1)$k - b)) }))
  out <- rbind(out, data.frame(N = N, Delta = Delta, d = d, R = R, r_or = 1, r_dir = N / d, q = N / d,
    pooled_mean = mean(r[, 1]), pooled_se = sd(r[, 1]) / sqrt(R), oracle_mean = mean(r[, 2]), oracle_se = sd(r[, 2]) / sqrt(R),
    predicted = cf[1] * 1 + cf[2] * N / d, c1 = cf[1], c2 = cf[2]))
  cat(sprintf("d=%d pooled=%.2f (%.2f) oracle=%.2f predicted=%.2f\n", d, mean(r[, 1]), sd(r[, 1]) / sqrt(R), mean(r[, 2]), cf[1] + cf[2] * N / d))
}
write_csv(out, "s1c_dscale.csv")
