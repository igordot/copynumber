test_that("plotSample() runs without error on valid input", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(
    plotSample(data = subLymphoma, segments = singleSeg, sample = 1)
  )
})

test_that("plotSample() rejects missing data/segments", {
  pdf(NULL)
  on.exit(dev.off())
  expect_error(plotSample())
})
