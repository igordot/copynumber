test_that("multipcf()'s documented example produces stable output", {
  # Values pinned from an R 4.6.1 run, not derived independently.
  data(lymphoma)
  sub.lymphoma <- subsetData(lymphoma, sample = 1:3)
  wins.lymph <- winsorize(sub.lymphoma, verbose = FALSE)

  multi.segments <- multipcf(
    data = wins.lymph,
    gamma = 12,
    Y = sub.lymphoma,
    verbose = FALSE
  )

  expect_equal(nrow(multi.segments), 54)
  expect_equal(sum(multi.segments$n.probes), 3091)
  expect_equal(
    as.numeric(multi.segments[1, c("X01.B1", "X01.B2", "X01.B3")]),
    c(-0.0439, -0.0324, -0.0643),
    tolerance = 0.01
  )
})

test_that("multipcf() rejects an invalid assembly and names every valid build", {
  data(lymphoma)
  err <- expect_error(multipcf(lymphoma, assembly = "bogus", verbose = FALSE))
  for (build in validAssemblies()) {
    expect_true(grepl(build, conditionMessage(err), fixed = TRUE))
  }
})

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
