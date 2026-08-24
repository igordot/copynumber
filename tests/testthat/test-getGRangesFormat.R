test_that("getGRangesFormat() converts pcf() segments to a GRanges", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  gr <- getGRangesFormat(singleSeg)

  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), nrow(singleSeg))
  seqnames <- as.character(GenomicRanges::seqnames(gr))
  expect_equal(seqnames, as.character(singleSeg$chrom))
  expect_equal(names(gr), singleSeg$sampleID)
  expect_equal(colnames(S4Vectors::mcols(gr)), c("arm", "n.probes", "mean"))
})
