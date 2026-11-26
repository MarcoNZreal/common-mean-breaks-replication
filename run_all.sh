#!/bin/sh
# Reproduces every simulation and numerical check in the paper and supplement, then rebuilds all tables. Run from this directory.
set -e
Rscript s1_rate.R            # Table 1 (two localization scales)
Rscript s1c_dscale.R         # Table 2 (out-of-sample 1/d test; uses s1_rate.csv)
Rscript s2_detect.R          # Table 3 (detection boundary)
Rscript s3_mask.R 1          # masking: null calibration and Delta = 2, 4
Rscript s3_mask.R 2          #          Delta = 6, 10
Rscript s3_mask.R 3          #          null false-positive rates
Rscript -e 'write.csv(rbind(read.csv("s3_mask_part1.csv"), read.csv("s3_mask_part2.csv")), "s3_mask.csv", row.names = FALSE)'
for D in 2 4 8 16 32; do Rscript s4_sparse.R A $D; done      # sparse design A, one cell per call
for D in 10 20 40 80; do Rscript s4_sparse.R B $D; done      # sparse design B
Rscript s4_sparse.R combine
for g in 0.5 1 2 4 6 8; do Rscript s6_exact.R oracle-exponent $g; Rscript s6_exact.R direction-exponent $g; done   # exact recovery
Rscript s6_exact.R combine
Rscript s5_local.R           # local versus global spacing for seeded NOT (Supplement S1)
Rscript verify_part1.R       # Supplement S1 checks (V1-V5)
Rscript verify_part2.R       # Supplement S1-S2 checks (quadrature, Frobenius, repartition, unequal norms)
Rscript make_tables.R        # writes ../tab_*.tex and ../verify/tab_*.tex
