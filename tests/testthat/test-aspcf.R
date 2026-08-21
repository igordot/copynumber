test_that("aspcf()'s documented example produces stable output", {
  # Values pinned from an R 4.6.1 run, not derived independently.
  data(logR)
  data(BAF)
  wins.logR <- winsorize(logR, verbose = FALSE)

  aspcf.segments <- aspcf(wins.logR, BAF, verbose = FALSE)

  expect_equal(nrow(aspcf.segments), 218)
  expect_equal(sum(aspcf.segments$n.probes), 20000)
  expect_equal(aspcf.segments$logR.mean[1], -0.2323, tolerance = 0.01)
  expect_equal(aspcf.segments$BAF.mean[1], 0.6916, tolerance = 0.01)
})

test_that("aspcf() rejects an invalid assembly and names every valid build", {
  data(logR)
  data(BAF)
  err <- expect_error(aspcf(logR, BAF, assembly = "bogus", verbose = FALSE))
  for (build in validAssemblies()) {
    expect_true(grepl(build, conditionMessage(err), fixed = TRUE))
  }
})

test_that("aspcf() reproduces the vignette's SNP-array workflow", {
  data(logR)
  data(BAF)
  logR.wins <- winsorize(logR, verbose = FALSE)

  allele.seg <- aspcf(logR.wins, BAF, verbose = FALSE)

  # aspcf output matches pcf's, plus an 8th column of segment BAF-values
  expect_named(
    allele.seg,
    c(
      "sampleID",
      "chrom",
      "arm",
      "start.pos",
      "end.pos",
      "n.probes",
      "logR.mean",
      "BAF.mean"
    )
  )
  # a real invariant, known without running aspcf(): segments partition
  # every probe in the input exactly once, per sample.
  for (id in unique(allele.seg$sampleID)) {
    n <- sum(allele.seg$n.probes[allele.seg$sampleID == id])
    expect_equal(n, nrow(logR.wins))
  }
})

test_that("aspcf() with return.est=TRUE returns logR estimates alongside segments", {
  data(logR)
  data(BAF)
  logR.wins <- winsorize(logR, verbose = FALSE)

  res <- aspcf(logR.wins, BAF, verbose = FALSE, return.est = TRUE)

  expect_named(res, c("logR_estimates", "segments"))
  expect_named(res$logR_estimates, c("chrom", "pos", colnames(logR.wins)[-c(1, 2)]))
  expect_equal(nrow(res$logR_estimates), nrow(logR.wins))
})
