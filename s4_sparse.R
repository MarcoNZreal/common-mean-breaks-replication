## Study 4: sparse changes with an unknown support. Off-centre breaks of both parities (alternating over replications).
## Estimators: oracle projection; pooled statistic on the known support; endpoint-block screening (Theorem 3.14, alpha=0.05);
## parity screen and split as previously analysed (A=1, midpoint fallback); dense pooling. Trimming lambda = 0.25 throughout,
## so that the end blocks have l = n/8 observations.
## run one cell per call: Rscript s4_sparse.R <A|B> <Delta>; then Rscript s4_sparse.R combine
source("common.R"); args <- commandArgs(TRUE); part <- args[1]; dir.create("parts", showWarnings = FALSE)
if (part == "combine") { for (p in c("A", "B")) { f <- Sys.glob(sprintf("parts/s4_%s_*.csv", p)); x <- do.call(rbind, lapply(f, read.csv)); write_csv(x[order(x$Delta), ], sprintf("s4_sparse_%s.csv", p)) }; quit() }
Del <- as.numeric(args[2]); set.seed(round(100 * Del) + if (part == "A") 205 else 206)
lam <- 0.25
endpoint_block <- function(Y, alpha) {
  n <- nrow(Y); N <- ncol(Y); l <- floor(lam * n / 2)
  D <- (colMeans(Y[(n - l + 1):n, , drop = FALSE]) - colMeans(Y[1:l, , drop = FALSE])) / sqrt(2 / l)
  sel <- which(abs(D) > sqrt(2 * log(2 * N / alpha)))
  if (!length(sel)) return(list(b = n %/% 2, sel = sel))
  Sf <- cumfull(Y[, sel, drop = FALSE]); ks <- Kset(n, lam)        # original admissible dates
  v <- vint(Sf, l, n - l, ks); list(b = ks[which.max(v)], sel = sel)
}
parity_split <- function(Y, A = 1) {
  n <- nrow(Y); N <- ncol(Y); Ye <- Y[seq(2, n, 2), , drop = FALSE]; Yo <- Y[seq(1, n, 2), , drop = FALSE]; np <- nrow(Ye)
  sel <- which(coordmax(Ye, lam / 2) > 2 * (1 + A) * log(N * np))
  k <- if (length(sel)) pooled(Yo[, sel, drop = FALSE], lam / 2)$k else np %/% 2
  list(b = 2 * k, sel = sel)
}
study <- function(tag, N, s, n, bs, grid, R) {
  out <- NULL; meth <- c("oracle", "known", "endpoint", "parity", "dense")
  for (Delta in grid) {
    e <- matrix(NA, R, length(meth), dimnames = list(NULL, meth)); nsel <- matrix(NA, R, 4)
    for (r in 1:R) {
      b <- bs[1 + (r %% 2)]
      S <- sample.int(N, s); delta <- numeric(N); delta[S] <- sample(c(-1, 1), s, TRUE) * sqrt(Delta / s)
      Y <- matrix(rnorm(n * N), n); Y[(b + 1):n, ] <- sweep(Y[(b + 1):n, , drop = FALSE], 2, delta, "+")
      e[r, "oracle"] <- pooled(Y %*% (delta / sqrt(Delta)), lam)$k - b
      e[r, "known"] <- pooled(Y[, S, drop = FALSE], lam)$k - b
      e[r, "dense"] <- pooled(Y, lam)$k - b
      eb <- endpoint_block(Y, 0.05); ps <- parity_split(Y)
      e[r, "endpoint"] <- eb$b - b; e[r, "parity"] <- ps$b - b
      nsel[r, ] <- c(length(intersect(eb$sel, S)), length(eb$sel), length(intersect(ps$sel, S)), length(ps$sel))
    }
    row <- data.frame(design = tag, N = N, s = s, n = n, Delta = Delta, R = R, per_coord_snr = Delta * min(bs[1], n - bs[1]) / s,
      ep_true = mean(nsel[, 1]), ep_sel = mean(nsel[, 2]), pa_true = mean(nsel[, 3]), pa_sel = mean(nsel[, 4]))
    for (m in meth) { row[[paste0(m, "_mae")]] <- mean(abs(e[, m])); row[[paste0(m, "_se")]] <- sd(abs(e[, m])) / sqrt(R)
                      row[[paste0(m, "_exact")]] <- mean(e[, m] == 0) }
    out <- rbind(out, row)
    cat(sprintf("%s Delta=%g  MAE: %s | exact: %s\n", tag, Delta, paste(meth, round(colMeans(abs(e)), 2), collapse = " "), paste(round(colMeans(e == 0), 2), collapse = " ")))
  }
  out
}
if (part == "A") write_csv(study("A", 1000, 10, 500, c(175, 176), Del, 300), sprintf("parts/s4_A_%g.csv", Del))   # Delta in 2, 4, 8, 16, 32
if (part == "B") write_csv(study("B", 10000, 5, 200, c(70, 71), Del, 150), sprintf("parts/s4_B_%g.csv", Del))   # Delta in 10, 20, 40, 80
