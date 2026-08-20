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
