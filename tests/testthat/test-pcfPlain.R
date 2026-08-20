test_that("pcfPlain() segments a single track given only position + values", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)

  pos.data <- lymph.wins[lymph.wins[, 1] == 1, c(2, 3)]
  res <- pcfPlain(pos.data = pos.data, gamma = 12, verbose = FALSE)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("sampleID", "start.pos", "end.pos", "n.probes", "mean"))
  expect_equal(unique(res$sampleID), "X01.B1")
  expect_equal(sum(res$n.probes), nrow(pos.data))
})
