test_that("plotGenome() runs without error on valid input", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(
    plotGenome(data = subLymphoma, segments = singleSeg, sample = 1)
  )
})

test_that("plotGenome() rejects missing data/segments and a bad pos.unit", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_error(plotGenome())
  expect_error(plotGenome(segments = singleSeg, pos.unit = "bogus"))
})
