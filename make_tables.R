dir.create("../verify", showWarnings = FALSE)   # supplement tables are written to ../verify/
## Builds every LaTeX table used by the paper (../tab_*.tex) and the supplement (../verify/tab_*.tex) from the CSV output of the R studies.
fmt <- function(m, s) if (m >= 10) sprintf("%.1f (%.1f)", m, s) else sprintf("%.2f (%.2f)", m, s)
f1 <- function(x) if (x >= 10) sprintf("%.1f", x) else sprintf("%.2f", x)
w <- function(lines, f) writeLines(lines, f)
## Table 1: two scales
r <- read.csv("s1_rate.csv"); use <- r$pooled_mean < 40; cf <- coef(lm(pooled_mean ~ 0 + r_or + r_dir, data = r[use, ]))
r2 <- summary(lm(pooled_mean ~ 0 + r_or + r_dir, data = r[use, ]))$r.squared
oc <- sum(r$r_or * r$oracle_mean) / sum(r$r_or^2)
L <- c("\\begin{tabular}{rrrrrrrrrr}\\toprule", "Panel & $N$ & $\\Delta$ & $R$ & $q$ & $r_{\\rm or}$ & $r_{\\rm dir}$ & oracle & pooled & fit \\\\ \\midrule")
for (i in seq_len(nrow(r))) { x <- r[i, ]; fit <- min(cf[1] * x$r_or + cf[2] * x$r_dir, 100)
  L <- c(L, sprintf("%s & %d & %g & %d & %.3g & %.3g & %.3g & %s & %s%s & %.1f \\\\", x$panel, x$N, x$Delta, x$R, x$q, x$r_or, x$r_dir,
    fmt(x$oracle_mean, x$oracle_se), fmt(x$pooled_mean, x$pooled_se), if (x$pooled_mean >= 40) "$^\\ast$" else "", fit)) }
w(c(L, "\\bottomrule\\end{tabular}"), "../tab_rate.tex")
w(sprintf("c1=%.2f c2=%.2f R2=%.4f oracle_c=%.2f nfit=%d", cf[1], cf[2], r2, oc, sum(use)), "rate_fit.txt")
## Table 2: d scaling
d <- read.csv("s1c_dscale.csv")
L <- c("\\begin{tabular}{rrrrrrr}\\toprule", "$d$ & $R$ & $r_{\\rm dir}$ & $q$ & oracle & pooled & prediction \\\\ \\midrule")
for (i in seq_len(nrow(d))) { x <- d[i, ]; L <- c(L, sprintf("%d & %d & %g & %g & %s & %s & %.1f \\\\", x$d, x$R, x$r_dir, x$q, fmt(x$oracle_mean, x$oracle_se), fmt(x$pooled_mean, x$pooled_se), x$predicted)) }
w(c(L, "\\bottomrule\\end{tabular}"), "../tab_dscale.tex")
## Table 3: detection
p <- read.csv("s2_detect.csv"); Ns <- sort(unique(p$N)); rh <- sort(unique(p$rho))
L <- c("\\begin{tabular}{rrrr}\\toprule", paste0("$\\rho=\\Delta d/(\\sigma^2\\sqrt N)$ & ", paste(sprintf("$N=%d$", Ns), collapse = " & "), " \\\\ \\midrule"))
for (x in rh) L <- c(L, paste0(sprintf("%g & ", x), paste(sprintf("%.2f", sapply(Ns, function(n) p$power[p$N == n & p$rho == x])), collapse = " & "), " \\\\"))
w(c(L, "\\bottomrule\\end{tabular}"), "../tab_detect.tex")
## Table 4: masking; Table S: repartition
m <- read.csv("s3_mask.csv"); nl <- read.csv("s3_null.csv")
L <- c("\\begin{tabular}{rccccc}\\toprule", " & \\multicolumn{2}{c}{binary segmentation} & \\multicolumn{3}{c}{seeded NOT} \\\\ \\cmidrule(lr){2-3}\\cmidrule(lr){4-6}",
       "$\\Delta$ per jump & $\\lambda=0.05$ & $\\lambda=1/4$ & basic rule & analysed rule & analysed, $\\widehat\\sigma$ \\\\ \\midrule")
for (i in seq_len(nrow(m))) { x <- m[i, ]; L <- c(L, sprintf("%g & %.2f & %.2f & %.2f & %.2f & %.2f \\\\", x$Delta, x$bs1, x$bs4, x$basic, x$inter, x$plug)) }
L <- c(L, "\\midrule", sprintf("null (no breaks) & %.2f & %.2f & %.2f & %.2f & %.2f \\\\", nl$bs1, nl$bs4, nl$basic, nl$inter, nl$plug))
w(c(L, "\\bottomrule\\end{tabular}"), "../tab_mask.tex")
L <- c("\\begin{tabular}{rrcc}\\toprule", "$\\Delta$ per jump & dates scored & before repartition & after repartition \\\\ \\midrule")
for (i in seq_len(nrow(m))) { x <- m[i, ]; L <- c(L, sprintf("%g & %d & %.2f (%.2f) & %.2f (%.2f) \\\\", x$Delta, x$n_scored, x$mae_before, x$mae_before_se, x$mae_after, x$mae_after_se)) }
w(c(L, "\\bottomrule\\end{tabular}"), "../tab_repart.tex")
## Table S: sparse (endpoint-block versus parity split), trimming 0.25, off-centre breaks of both parities
sp <- rbind(read.csv("s4_sparse_A.csv"), read.csv("s4_sparse_B.csv"))
g2 <- function(x) sprintf("%.2f", x)
L <- c("{\\setlength{\\tabcolsep}{3.5pt}\\begin{tabular}{crrccccccccc}\\toprule",
       " & & & & \\multicolumn{5}{c}{mean absolute error} & \\multicolumn{3}{c}{exact recovery} \\\\ \\cmidrule(lr){5-9}\\cmidrule(lr){10-12}",
       "Design & $\\Delta$ & $\\Delta d/(s\\sigma^2)$ & selected & oracle & known $S$ & endpoint & parity & dense & endpoint & parity & dense \\\\ \\midrule")
for (i in seq_len(nrow(sp))) { x <- sp[i, ]
  if (i > 1 && sp$design[i] != sp$design[i - 1]) L <- c(L, "\\midrule")
  L <- c(L, sprintf("%s & %g & %g & %.2f/%.2f & %s & %s & %s & %s & %s & %s & %s & %s \\\\", x$design, x$Delta, x$per_coord_snr, x$ep_true, x$ep_sel,
    g2(x$oracle_mae), g2(x$known_mae), g2(x$endpoint_mae), g2(x$parity_mae), g2(x$dense_mae), g2(x$endpoint_exact), g2(x$parity_exact), g2(x$dense_exact))) }
w(c(L, "\\bottomrule\\end{tabular}}"), "../tab_sparse.tex")
## Table S: exact recovery (Study 6)
ex <- read.csv("s6_exact.csv"); o <- ex[ex$design == "oracle-exponent", ]; dd <- ex[ex$design == "direction-exponent", ]
ci <- function(p, lo, hi) sprintf("%.3f [%.3f, %.3f]", p, lo, hi)
L <- c("\\begin{tabular}{rcccc}\\toprule",
       " & oracle exponent ($N=10$) & \\multicolumn{3}{c}{direction exponent ($N=8000$)} \\\\ \\cmidrule(lr){2-2}\\cmidrule(lr){3-5}",
       "$\\gamma$ & pooled & $\\Delta/\\sigma^2$ & pooled & known direction \\\\ \\midrule")
for (i in seq_len(nrow(o))) L <- c(L, sprintf("%g & %s & %.1f & %s & %s \\\\", o$gamma[i], ci(o$miss[i], o$lo[i], o$hi[i]), dd$Delta[i],
    ci(dd$miss[i], dd$lo[i], dd$hi[i]), ci(dd$miss_oracle[i], dd$lo_or[i], dd$hi_or[i])))
w(c(L, "\\bottomrule\\end{tabular}"), "../tab_exact.tex")
## supplement tables
v1 <- read.csv("verify_v1.csv"); sci <- function(x) { e <- floor(log10(x)); sprintf("$%.1f\\times10^{%d}$", x / 10^e, e) }
L <- c("\\begin{tabular}{rrccc}\\toprule", "$n$ & max.\\ relative error of formula & \\multicolumn{3}{c}{min.\\ of $\\{V_I(b)-V_I(k)\\}/(\\lambda\\Delta|k-b|)$} \\\\ \\cmidrule(l){3-5}",
       " & (all $b,k$) & $\\lambda=0.05$ & $\\lambda=0.1$ & $\\lambda=0.25$ \\\\ \\midrule")
for (i in seq_len(nrow(v1))) L <- c(L, sprintf("%d & %s & %.3f & %.3f & %.3f \\\\", v1$n[i], sci(v1$maxrel[i]), v1$m05[i], v1$m10[i], v1$m25[i]))
w(c(L, "\\bottomrule\\end{tabular}"), "../verify/tab_v1.tex")
v3 <- read.csv("verify_v3.csv")
L <- c("\\begin{tabular}{rrrccc}\\toprule", "$N$ & $a$ & $s^2$ & Monte Carlo $1+\\chi^2$ (s.e.) & $\\cosh(a^2/s^2)^N$ & $\\exp\\{Na^4/(2s^4)\\}$ \\\\ \\midrule")
for (i in seq_len(nrow(v3))) L <- c(L, sprintf("%d & %g & %g & %.4f (%.4f) & %.4f & %.4f \\\\", v3$N[i], v3$a[i], v3$s2[i], v3$mc[i], v3$se[i], v3$exact[i], v3$bound[i]))
w(c(L, "\\bottomrule\\end{tabular}"), "../verify/tab_v3.tex")
v4 <- read.csv("verify_v4.csv")
L <- c("\\begin{tabular}{rrrrrrc}\\toprule", "$T$ & intervals & total length & $d$ & break locations checked & failures & worst centring \\\\ \\midrule")
for (i in seq_len(nrow(v4))) L <- c(L, sprintf("%d & %d & %d & %d & %d & %d & %.4f \\\\", v4$T[i], v4$intervals[i], v4$total_length[i], v4$d[i], v4$locations[i], v4$failures[i], v4$worst[i]))
w(c(L, "\\bottomrule\\end{tabular}"), "../verify/tab_v4.tex")
v5 <- read.csv("verify_v5.csv")
L <- c("\\begin{tabular}{rrrccc}\\toprule", "$N$ & $T$ & $m$ & mean of $\\hat\\sigma/\\sigma-1$ & s.d. & $(NT)^{-1/2}+m/T$ \\\\ \\midrule")
for (i in seq_len(nrow(v5))) L <- c(L, sprintf("%d & %d & %d & %+.4f & %.4f & %.4f \\\\", v5$N[i], v5$T[i], v5$m[i], v5$bias[i], v5$sd[i], v5$rate[i]))
w(c(L, "\\bottomrule\\end{tabular}"), "../verify/tab_v5.tex")
v6 <- read.csv("verify_v6.csv"); se <- function(x) { e <- floor(log10(abs(x))); sprintf("%.2f\\mathrm{e}{%d}", x / 10^e, e) }
L <- c("\\begin{tabular}{rrccc}\\toprule", "$a^2T/\\sigma^2$ & $(b_g,b_h)$ & Monte Carlo $D(p_g\\Vert p_h)$ ($\\pm2$ s.e.) & $a^4\\Vert gg'-hh'\\Vert_F^2/\\sigma^4$ & ratio \\\\ \\midrule")
for (i in seq_len(nrow(v6))) L <- c(L, sprintf("%g & (%d,%d) & $%s \\pm %s$ & $%s$ & %.2f \\\\", v6$c0[i], v6$bg[i], v6$bh[i], se(v6$kl[i]), se(v6$se2[i]), se(v6$rhs[i]), v6$kl[i] / v6$rhs[i]))
w(c(L, "\\bottomrule\\end{tabular}"), "../verify/tab_v6.tex")
L <- c("\\begin{tabular}{rrrrrrr}\\toprule", "Panel & $N$ & $\\Delta$ & $q$ & oracle median & pooled median & pooled/oracle mean \\\\ \\midrule")
for (i in seq_len(nrow(r))) { x <- r[i, ]; L <- c(L, sprintf("%s & %d & %g & %.3g & %g & %g & %.2f \\\\", x$panel, x$N, x$Delta, x$q, x$oracle_med, x$pooled_med, x$pooled_mean / x$oracle_mean)) }
w(c(L, "\\bottomrule\\end{tabular}"), "../verify/tab_s1extra.tex")
L <- c("\\begin{tabular}{rrrr}\\toprule", "$N$ & $N\\sigma^2$ & null mean of $\\max_k\\widehat V_I(k)$ & 95\\% critical value \\\\ \\midrule")
for (n in Ns) { x <- p[p$N == n, ][1, ]; L <- c(L, sprintf("%d & %d & %.1f & %.1f \\\\", n, n, x$null_mean, x$crit)) }
w(c(L, "\\bottomrule\\end{tabular}"), "../verify/tab_s2crit.tex")
cat(readLines("rate_fit.txt"), "\n")
## Table S: local versus global spacing (Study 5)
lc <- read.csv("s5_local.csv"); iso <- lc[lc$layout == "isolated", ]; cro <- lc[lc$layout == "crowded", ]
L <- c("\\begin{tabular}{rcccc}\\toprule",
       " & \\multicolumn{2}{c}{isolated, $d_w=300$} & \\multicolumn{2}{c}{crowded, $d_w=150$} \\\\ \\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
       "$\\Delta_w$ & weak break & all three & weak break & all three \\\\ \\midrule")
for (i in seq_len(nrow(iso))) L <- c(L, sprintf("%.1f & %.2f & %.2f & %.2f & %.2f \\\\", iso$Delta_w[i], iso$weak[i], iso$all3[i], cro$weak[i], cro$all3[i]))
L <- c(L, "\\midrule", sprintf("null (no breaks) & \\multicolumn{4}{c}{%.3f} \\\\", lc$null_fp[1]))
w(c(L, "\\bottomrule\\end{tabular}"), "../tab_local.tex")
