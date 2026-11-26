## Shared helpers for all simulation studies. Data matrices are n x p (time in rows, series in columns).
colcum <- function(X) {                      # column-wise cumulative sums via one vectorized cumsum
  n <- nrow(X); p <- ncol(X); S <- matrix(cumsum(X), n)
  if (p > 1) S <- S - rep(c(0, S[n, -p]), each = n)
  S
}
Kset <- function(n, lam) seq(ceiling(lam * n), floor(n - lam * n))
## pooled standardized squared CUSUM on (0,n]: returns list(k = argmax, v = values over ks, ks)
pooled <- function(X, lam) {
  X <- as.matrix(X); n <- nrow(X); ks <- Kset(n, lam); S <- colcum(X)
  C <- S[ks, , drop = FALSE] - outer(ks / n, S[n, ])
  v <- n / (ks * (n - ks)) * rowSums(C^2)
  list(k = ks[which.max(v)], v = v, ks = ks)
}
## per-coordinate max over k of the standardized CUSUM squared
coordmax <- function(X, lam) {
  n <- nrow(X); ks <- Kset(n, lam); S <- colcum(X)
  C <- S[ks, , drop = FALSE] - outer(ks / n, S[n, ])
  apply(C^2 * (n / (ks * (n - ks))), 2, max)
}
## interval statistic from cumulative sums Sfull ((T+1) x p, first row zero): Vhat_(a,c](k) for integer ks
vint <- function(Sfull, a, c, ks) {
  n <- c - a; Sa <- Sfull[a + 1, ]; Sc <- Sfull[c + 1, ]
  C <- Sfull[ks + 1, , drop = FALSE] - matrix(Sa, length(ks), length(Sa), byrow = TRUE) - outer((ks - a) / n, Sc - Sa)
  n / ((ks - a) * (c - ks)) * rowSums(C^2)
}
cumfull <- function(Y) rbind(0, colcum(Y))
central <- function(a, c, lam) {
  n <- c - a; lo <- a + max(1, ceiling(lam * n)); hi <- c - max(1, ceiling(lam * n))
  if (hi >= lo) lo:hi else a + n %/% 2
}
## seeded family (26): dyadic lengths, shift L/8, plus a right-aligned interval per scale
seeded_family <- function(T) {
  A <- integer(0); C <- integer(0); L <- 4
  while (L <= T) {
    q <- max(1, L %/% 8); st <- seq(0, T - L, by = q); if (tail(st, 1) + L < T) st <- c(st, T - L)
    A <- c(A, st); C <- c(C, st + L); L <- 2 * L
  }
  list(a = A, c = C)
}
## max over the central half of each seeded interval of Vhat - N*s2
seeded_scan <- function(Sfull, fam, N, s2 = 1) {
  m <- length(fam$a); M <- numeric(m); K <- integer(m)
  for (i in seq_len(m)) {
    ks <- central(fam$a[i], fam$c[i], 0.25); v <- vint(Sfull, fam$a[i], fam$c[i], ks) - N * s2
    j <- which.max(v); M[i] <- v[j]; K[i] <- ks[j]
  }
  list(M = M, K = K)
}
## narrowest-over-threshold selection. rule = "basic" (remove intervals containing the date) or "intersect" (analysed rule)
seeded_not <- function(scan, fam, tau, rule = "intersect") {
  act <- rep(TRUE, length(fam$a)); Ln <- fam$c - fam$a; found <- integer(0)
  repeat {
    cand <- which(act & scan$M > tau); if (!length(cand)) break
    i <- cand[order(Ln[cand], fam$a[cand])[1]]; k <- scan$K[i]; found <- c(found, k)
    if (rule == "basic") act <- act & !(fam$a < k & k < fam$c) else act <- act & !(fam$a < fam$c[i] & fam$c > fam$a[i])
  }
  sort(found)
}
## pooled binary segmentation: splits at least max(2, floor(lam*n)) from either end; root threshold at every node
bs_stat <- function(Sfull, a, c, lam, N) {
  n <- c - a; g <- max(2, floor(lam * n)); if (c - a < 4 || a + g > c - g) return(list(m = -Inf, k = NA))
  ks <- (a + g):(c - g); v <- vint(Sfull, a, c, ks) - N; j <- which.max(v); list(m = v[j], k = ks[j])
}
bs <- function(Sfull, T, tau, lam, N) {
  found <- integer(0); stack <- list(c(0, T))
  while (length(stack)) {
    iv <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL
    r <- bs_stat(Sfull, iv[1], iv[2], lam, N)
    if (r$m > tau) { found <- c(found, r$k); stack[[length(stack) + 1]] <- c(iv[1], r$k); stack[[length(stack) + 1]] <- c(r$k, iv[2]) }
  }
  sort(found)
}
## repartition on symmetric intervals (Theorem 4.3), splits in the central half
repartition <- function(Sfull, dates, T) {
  ext <- c(0, dates, T); out <- integer(length(dates))
  for (j in seq_along(dates)) {
    dj <- min(ext[j + 1] - ext[j], ext[j + 2] - ext[j + 1]); h <- dj %/% 2
    a <- max(0, ext[j + 1] - h); c <- min(T, ext[j + 1] + h); ks <- central(a, c, 0.25)
    out[j] <- ks[which.max(vint(Sfull, a, c, ks))]
  }
  out
}
## robust variance estimate (32) from non-overlapping first differences
sig_hat <- function(Y) {
  T <- nrow(Y); R <- T %/% 2; D <- Y[2 * (1:R), , drop = FALSE] - Y[2 * (1:R) - 1, , drop = FALSE]
  median(abs(D)) / (sqrt(2) * qnorm(0.75))
}
write_csv <- function(df, f) write.csv(df, f, row.names = FALSE)
