## Study 3: masking. T=400, N=10, sigma=1, pulse on (150,190]. All methods are applied to the same datasets.
## Methods: pooled BS with lambda=0.05 and 1/4; seeded NOT with the basic and the analysed (intersection) removal rule;
## seeded NOT (analysed rule) with the plug-in variance; repartition after seeded NOT (analysed rule).
source("common.R"); part <- as.integer(commandArgs(TRUE)[1]); set.seed(104 + 10 * part)
T <- 400; N <- 10; truth <- c(150, 190); fam <- seeded_family(T)
ok <- function(f) length(f) == 2 && all(abs(f - truth) <= 10)
if (part == 1) {
  nb1 <- nb4 <- ns <- numeric(1000)
  for (r in 1:1000) { Sf <- cumfull(matrix(rnorm(T * N), T))
    nb1[r] <- bs_stat(Sf, 0, T, 0.05, N)$m; nb4[r] <- bs_stat(Sf, 0, T, 0.25, N)$m; ns[r] <- max(seeded_scan(Sf, fam, N)$M) }
  tau <- c(bs1 = unname(quantile(nb1, 0.95)), bs4 = unname(quantile(nb4, 0.95)), nt = unname(quantile(ns, 0.95)))
  write_csv(data.frame(t(tau)), "s3_thresholds.csv")
} else { tt <- read.csv("s3_thresholds.csv"); tau <- c(bs1 = tt$bs1, bs4 = tt$bs4, nt = tt$nt) }
print(round(tau, 2))
run <- function(Y) {
  Sf <- cumfull(Y); sc <- seeded_scan(Sf, fam, N); s2 <- sig_hat(Y)^2
  scp <- list(M = sc$M + N * (1 - s2), K = sc$K)               # same maximizers, centring N * s2
  list(bs1 = bs(Sf, T, tau["bs1"], 0.05, N), bs4 = bs(Sf, T, tau["bs4"], 0.25, N),
       basic = seeded_not(sc, fam, tau["nt"], "basic"), inter = seeded_not(sc, fam, tau["nt"], "intersect"),
       plug = seeded_not(scp, fam, tau["nt"] * s2, "intersect"), Sf = Sf)
}
out <- NULL; R <- 500
grid <- list(c(2, 4), c(6, 10), numeric(0))[[part]]
for (Delta in grid) {
  succ <- matrix(0, R, 5, dimnames = list(NULL, c("bs1", "bs4", "basic", "inter", "plug"))); ndates <- succ; before <- after <- numeric(0)
  for (r in 1:R) {
    Y <- matrix(rnorm(T * N), T); Y[151:190, ] <- Y[151:190, ] + sqrt(Delta / N); f <- run(Y)
    for (m in colnames(succ)) { succ[r, m] <- ok(f[[m]]); ndates[r, m] <- length(f[[m]]) }
    if (ok(f$inter)) { rp <- repartition(f$Sf, f$inter, T); before <- c(before, abs(f$inter - truth)); after <- c(after, abs(rp - truth)) }
  }
  out <- rbind(out, data.frame(design = "pulse", Delta = Delta, R = R, t(colMeans(succ)), setNames(as.data.frame(t(colMeans(ndates))), paste0("dates_", colnames(succ))),
    mae_before = mean(before), mae_before_se = sd(before) / sqrt(length(before)), mae_after = mean(after), mae_after_se = sd(after) / sqrt(length(after)), n_scored = length(before)))
  cat(sprintf("Delta=%g: %s\n", Delta, paste(names(colMeans(succ)), round(colMeans(succ), 3), collapse = " ")))
}
if (part == 3) {
fp <- matrix(0, 1000, 5, dimnames = list(NULL, c("bs1", "bs4", "basic", "inter", "plug")))
for (r in 1:1000) { f <- run(matrix(rnorm(T * N), T)); for (m in colnames(fp)) fp[r, m] <- length(f[[m]]) > 0 }
nullrow <- data.frame(design = "null", Delta = 0, R = 1000, t(colMeans(fp)))
cat("null FP:", round(colMeans(fp), 3), "\n")
write_csv(nullrow, "s3_null.csv") }
if (part < 3) write_csv(out, sprintf("s3_mask_part%d.csv", part))
