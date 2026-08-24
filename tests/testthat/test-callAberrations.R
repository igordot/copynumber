test_that("callAberrations() classifies pcf (single-sample) segments", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  called <- callAberrations(segments = singleSeg, thres.gain = 0.2)

  expect_equal(colnames(called), c(colnames(singleSeg)[-7], "call"))
  expect_true(all(called$call %in% c("normal", "gain", "loss")))
  expect_true(all(called$call[called$mean > 0.2] == "gain"))
  expect_true(all(called$call[called$mean < -0.2] == "loss"))
})

test_that("callAberrations() classifies multipcf (wide) segments per sample", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  multiSeg <- multipcf(data = lymphWins, verbose = FALSE)

  called <- callAberrations(segments = multiSeg, thres.gain = 0.2)

  sample_cols <- c("X01.B1", "X01.B2", "X01.B3")
  expect_true(all(unlist(called[sample_cols]) %in% c("normal", "gain", "loss")))
  expect_equal(nrow(called), nrow(multiSeg))
})
