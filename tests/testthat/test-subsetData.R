test_that("subsetData() filters by chrom and by sample", {
  data(lymphoma)

  by_chrom <- subsetData(data = lymphoma, chrom = c(1, 2))
  expect_equal(unique(by_chrom[, 1]), c(1, 2))
  expect_equal(ncol(by_chrom), ncol(lymphoma))

  by_sample <- subsetData(data = lymphoma, sample = 1:3)
  expect_equal(ncol(by_sample), 2 + 3)
  expect_equal(nrow(by_sample), nrow(lymphoma))
})
