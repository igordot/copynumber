test_that("getGRangesFormat() converts pcf() segments to a GRanges", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  gr <- getGRangesFormat(single.seg)

  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), nrow(single.seg))
  seqnames <- as.character(GenomicRanges::seqnames(gr))
  expect_equal(seqnames, as.character(single.seg$chrom))
  expect_equal(names(gr), single.seg$sampleID)
  expect_equal(colnames(S4Vectors::mcols(gr)), c("arm", "n.probes", "mean"))
})
