test_that("plotHeatmap() runs without error on valid input", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plotHeatmap(segments = singleSeg, upper.lim = 0.3))
})

test_that("plotHeatmap() rejects a bad pos.unit and a non-positive upper.lim", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_error(
    plotHeatmap(segments = singleSeg, upper.lim = 0.3, pos.unit = "bogus")
  )
  expect_error(plotHeatmap(segments = singleSeg, upper.lim = -1))
})
