test_that("plotFreq() runs without error on valid input", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(
    plotFreq(segments = singleSeg, thres.gain = 0.2, thres.loss = -0.1)
  )
})

test_that("plotFreq() rejects a bad pos.unit", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_error(
    plotFreq(segments = singleSeg, thres.gain = 0.2, pos.unit = "bogus")
  )
})
