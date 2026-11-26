Simulation code for "Common Mean Breaks in Panels: Non-asymptotic Localization, Exact Recovery and the Cost of an
Unknown Sparse Support". Base R only; no packages are needed.

run_all.sh reproduces every simulation and rebuilds the tables; make_tables.R writes ../tab_*.tex, so run it from a
folder placed next to the LaTeX sources. The .csv files hold the results reported in the paper and supplement.

  s1_rate.R, s1c_dscale.R   two localization scales and the out-of-sample 1/d test
  s2_detect.R               detection boundary
  s3_mask.R                 masking, seeded NOT and repartition
  s4_sparse.R               sparse changes: endpoint-block screening, parity splitting, dense pooling
                            (one cell per call: Rscript s4_sparse.R A 8; then Rscript s4_sparse.R combine)
  s5_local.R                local versus global spacing for seeded NOT
  s6_exact.R                exact recovery (one cell per call: Rscript s6_exact.R direction-exponent 4;
                            then Rscript s6_exact.R combine)
  verify_part1.R, verify_part2.R   additional numerical checks of the population margin and related identities
