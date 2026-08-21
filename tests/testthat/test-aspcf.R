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
