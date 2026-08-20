test_that("multipcf() reproduces the vignette's lymphoma workflow", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)

  multi.seg <- multipcf(data = lymph.wins, verbose = FALSE)

  # multipcf finds common breakpoints, so segments are wide: one row per
  # segment with a value column per sample, rather than pcf's long format.
  sample_cols <- c("X01.B1", "X01.B2", "X01.B3")
  expect_named(
    multi.seg,
    c("chrom", "arm", "start.pos", "end.pos", "n.probes", sample_cols)
  )
  # a real invariant: segments partition every probe in the input exactly
  # once (breakpoints are shared across samples, so there is one n.probes
  # column, not one per sample).
  expect_equal(sum(multi.seg$n.probes), nrow(lymph.wins))
})
