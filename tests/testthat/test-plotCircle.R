test_that("plotCircle() runs without error or warning on valid input", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plotCircle(segments = singleSeg, thres.gain = 0.15))
  expect_no_warning(plotCircle(segments = singleSeg, thres.gain = 0.15))
})
