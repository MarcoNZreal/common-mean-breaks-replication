## Numerical checks behind the results marked [P] (Supplement S1-S2). Writes verify_*.csv.
source("common.R"); set.seed(202)
## Gauss-Hermite nodes (probabilists') by Golub-Welsch
gh <- function(k) { J <- matrix(0, k, k); off <- sqrt(1:(k - 1)); J[cbind(1:(k - 1), 2:k)] <- off; J[cbind(2:k, 1:(k - 1))] <- off
  e <- eigen(J, symmetric = TRUE); list(x = e$values, w = e$vectors[1, ]^2) }
q <- gh(200); logcosh <- function(t) abs(t) + log1p(exp(-2 * abs(t))) - log(2)
Fq <- function(m, s) sum(q$w * logcosh(m + s * q$x))
gc <- function(T, b) { v <- numeric(T); v[(b + 1):T] <- 1; v - mean(v) }
Dexact <- function(a, g, h) { s <- a * sqrt(sum(g^2)); Fq(a^2 * sum(g^2), s) - Fq(a^2 * abs(sum(g * h)), s) }
## direct Monte Carlo for three symmetric pairs (T=40), compared with the quadrature reduction
mcmp <- NULL
for (bb in c(19, 17, 12)) for (a in c(0.15, 0.4, 1.0)) {
  g <- gc(40, bb); h <- gc(40, 40 - bb); M <- 200000; xi <- sample(c(-1, 1), M, TRUE); Z <- matrix(rnorm(M * 40), M)
  Y <- xi * a * matrix(g, M, 40, byrow = TRUE) + Z; vals <- logcosh(a * (Y %*% g)) - logcosh(a * (Y %*% h))
  mcmp <- rbind(mcmp, data.frame(b = bb, a = a, mc = mean(vals), mc_se = sd(vals) / sqrt(M), quad = Dexact(a, g, h)))
}
write_csv(mcmp, "verify_mc_vs_quad.csv"); print(mcmp)
worst <- 0; ncase <- 0
for (T in c(40, 200, 1000)) for (b in ceiling(0.1 * T):(T / 2 - 1)) { g <- gc(T, b); h <- gc(T, T - b)
  fro <- sum(g^2)^2 + sum(h^2)^2 - 2 * sum(g * h)^2
  for (a in exp(seq(log(1e-3), log(3), length.out = 25))) { B <- a^4 / 4 * fro; if (B > 0) { worst <- max(worst, Dexact(a, g, h) / B); ncase <- ncase + 1 } } }
cat("quadrature cases:", ncase, " max ratio:", worst, "\n")
## Frobenius bound (25)/(22): ||g_b g_b' - g_{b+h} g_{b+h}'||_F^2 / (T h)
wf <- 0
for (T in c(100, 400, 1600)) for (b in ceiling(0.1 * T):floor(0.9 * T)) for (hh in c(1, 2, 5, 10, 30)) { if (b + hh > 0.9 * T) next
  g <- gc(T, b); h <- gc(T, b + hh); wf <- max(wf, (sum(g^2)^2 + sum(h^2)^2 - 2 * sum(g * h)^2) / (T * hh)) }
cat("Frobenius max ratio:", wf, "\n")
## repartition intervals (Theorem 4.3)
bad <- 0; ws <- Inf; wr <- Inf; nint <- 0
for (it in 1:200000) {
  m <- sample.int(5, 1); gaps <- sample(16:399, m + 1, TRUE); bb <- cumsum(c(0, gaps)); T <- tail(bb, 1); br <- bb[2:(m + 1)]
  d <- sapply(1:m, function(j) min(br[j] - bb[j], bb[j + 2] - br[j]))
  bt <- br + sapply(d, function(x) { k <- ceiling(x / 8) - 1; sample(-k:k, 1) }); ext <- c(0, bt, T)
  for (j in 1:m) { dt <- min(ext[j + 1] - ext[j], ext[j + 2] - ext[j + 1]); h <- dt %/% 2; a <- bt[j] - h; c <- bt[j] + h; nint <- nint + 1
    oth <- br[-j]
    if (!(a < br[j] && br[j] < c) || any(a < oth & oth < c) || a < 0 || c > T) { bad <- bad + 1; next }
    side <- min(br[j] - a, c - br[j]); ws <- min(ws, side / d[j]); wr <- min(wr, side / (c - a)) }
}
cat("repartition intervals:", nint, "failures:", bad, "min side/d:", ws, "min side/|I|:", wr, "\n")
write_csv(data.frame(quad_cases = ncase, quad_max = worst, frob_max = wf, rep_intervals = nint, rep_failures = bad, rep_min_side_d = ws, rep_min_side_len = wr), "verify_misc.csv")
## V6 unequal-norm mixture divergence (Supplement S2)
v6 <- NULL; Tm <- 12
for (c0 in c(0.05, 0.2)) { a <- sqrt(c0 / Tm)
  for (pr in list(c(5, 6), c(4, 7), c(3, 9))) { g <- gc(Tm, pr[1]); h <- gc(Tm, pr[2]); M <- 400000
    sg <- sample(c(-1, 1), M, TRUE); Y <- sg * a * matrix(g, M, Tm, byrow = TRUE) + matrix(rnorm(M * Tm), M)
    lm1 <- function(v) { l1 <- -0.5 * rowSums((Y - a * matrix(v, M, Tm, byrow = TRUE))^2); l2 <- -0.5 * rowSums((Y + a * matrix(v, M, Tm, byrow = TRUE))^2)
      mx <- pmax(l1, l2); mx + log(0.5 * exp(l1 - mx) + 0.5 * exp(l2 - mx)) }
    vals <- lm1(g) - lm1(h); rhs <- a^4 * (sum(g^2)^2 + sum(h^2)^2 - 2 * sum(g * h)^2)
    v6 <- rbind(v6, data.frame(c0 = c0, bg = pr[1], bh = pr[2], kl = mean(vals), se2 = 2 * sd(vals) / sqrt(M), rhs = rhs)) } }
write_csv(v6, "verify_v6.csv"); print(v6)
