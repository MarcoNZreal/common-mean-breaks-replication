## Numerical checks behind the results marked [P] (Supplement S1-S2). Writes verify_*.csv.
source("common.R"); set.seed(201)
## V1 exact interval criterion and margin (Lemma 2.1)
v1 <- NULL
for (n in c(101, 400)) {
  delta <- rnorm(5); D <- sum(delta^2); maxrel <- 0
  for (b in 1:(n - 1)) {
    mu <- matrix(0, n, 5); mu[(b + 1):n, ] <- matrix(delta, n - b, 5, byrow = TRUE)
    S <- colcum(mu); ks <- 1:(n - 1); C <- S[ks, , drop = FALSE] - outer(ks / n, S[n, ]); vd <- n / (ks * (n - ks)) * rowSums(C^2)
    f <- ifelse(ks <= b, D * ks * (n - b)^2 / (n * (n - ks)), D * b^2 * (n - ks) / (n * ks))
    maxrel <- max(maxrel, max(abs(vd - f) / pmax(abs(f), 1e-12)))
  }
  mg <- sapply(c(0.05, 0.1, 0.25), function(lam) {
    K <- Kset(n, lam); w <- Inf
    for (b in K) { k <- K[K != b]; Vk <- ifelse(k <= b, D * k * (n - b)^2 / (n * (n - k)), D * b^2 * (n - k) / (n * k))
      w <- min(w, min((D * b * (n - b) / n - Vk) / (lam * D * abs(k - b)))) }
    w })
  v1 <- rbind(v1, data.frame(n = n, maxrel = maxrel, m05 = mg[1], m10 = mg[2], m25 = mg[3]))
}
write_csv(v1, "verify_v1.csv"); print(v1)
## V2 OU correlation identity on the grid (Lemma 3.2)
T <- 200; ks <- 1:(T - 1); u <- ks / T; r <- log(u / (1 - u))
lo <- outer(ks, ks, pmin); hi <- outer(ks, ks, pmax); cv <- lo * (T - hi) / T; sdv <- sqrt(ks * (T - ks) / T)
v2 <- max(abs(cv / outer(sdv, sdv) - exp(-0.5 * abs(outer(r, r, "-"))))); write_csv(data.frame(T = T, maxerr = v2), "verify_v2.csv"); cat("V2", v2, "\n")
## V3 chi-square identity of the dense detection lower bound
v3 <- NULL; M <- 400000
for (pr in list(c(3, 0.4, 0.5), c(5, 0.3, 0.4), c(8, 0.25, 0.3))) {
  N <- pr[1]; a <- pr[2]; s2 <- pr[3]; X <- matrix(rnorm(M * N, sd = sqrt(s2)), M)
  L <- apply(cosh(a * X / s2), 1, prod) * exp(-N * a^2 / (2 * s2))
  v3 <- rbind(v3, data.frame(N = N, a = a, s2 = s2, mc = mean(L^2), se = sd(L^2) / sqrt(M), exact = cosh(a^2 / s2)^N, bound = exp(N * a^4 / (2 * s2^2))))
}
xs <- seq(0, 6, length.out = 601); cat("cosh bound holds:", all(cosh(xs) <= exp(xs^2 / 2) + 1e-12), "\n"); write_csv(v3, "verify_v3.csv"); print(v3)
## V4 seeded covering property, exhaustive over break locations (Lemma 4.1)
v4 <- NULL
for (T in c(500, 2000, 8192)) {
  fam <- seeded_family(T); Ln <- fam$c - fam$a
  for (d in c(16, 23, 50, 100, 257)) { if (2 * d > T) next
    okl <- Ln >= d / 4 & Ln <= d / 2; A <- fam$a[okl]; C <- fam$c[okl]; L <- Ln[okl]; fails <- 0; worst <- Inf
    for (b in d:(T - d)) { m <- A < b & b < C; dist <- pmin(b - A[m], C[m] - b) / L[m]; g <- dist >= 0.25
      if (!any(g)) fails <- fails + 1 else worst <- min(worst, max(dist[g])) }
    v4 <- rbind(v4, data.frame(T = T, intervals = length(fam$a), total_length = sum(Ln), d = d, locations = T - 2 * d + 1, failures = fails, worst = worst)) }
}
write_csv(v4, "verify_v4.csv"); print(v4)
## V5 robust variance estimator (Lemma 5.1)
v5 <- NULL
for (pr in list(c(1, 200, 0), c(1, 200, 3), c(10, 200, 3), c(10, 2000, 0), c(10, 2000, 10), c(50, 2000, 10))) {
  N <- pr[1]; T <- pr[2]; m <- pr[3]
  est <- replicate(400, { Y <- matrix(rnorm(T * N, sd = 1.7), T)
    if (m > 0) for (b in floor(seq(T / (m + 1), T - T / (m + 1), length.out = m))) Y[(b + 1):T, ] <- sweep(Y[(b + 1):T, , drop = FALSE], 2, rnorm(N, sd = 5), "+")
    sig_hat(Y) / 1.7 - 1 })
  v5 <- rbind(v5, data.frame(N = N, T = T, m = m, bias = mean(est), sd = sd(est), rate = (N * T)^(-1 / 2) + m / T))
}
write_csv(v5, "verify_v5.csv"); print(v5)
