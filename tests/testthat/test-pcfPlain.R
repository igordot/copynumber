test_that("pcfPlain() segments a single track given only position + values", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)

  posData <- lymphWins[lymphWins[, 1] == 1, c(2, 3)]
  res <- pcfPlain(pos.data = posData, gamma = 12, verbose = FALSE)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("sampleID", "start.pos", "end.pos", "n.probes", "mean"))
  expect_equal(unique(res$sampleID), "X01.B1")
  expect_equal(sum(res$n.probes), nrow(posData))
})
