test_that("subsetSegments() filters segments by chromosome", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  ss <- subsetSegments(segments = single.seg, chrom = 1)

  expect_s3_class(ss, "data.frame")
  expect_equal(unique(ss$chrom), 1)
  expect_true(nrow(ss) < nrow(single.seg))
  expect_true(nrow(ss) > 0)
})
