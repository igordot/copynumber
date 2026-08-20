test_that("plotGenome() runs without error on valid input", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(
    plotGenome(data = sub.lymphoma, segments = single.seg, sample = 1)
  )
})

test_that("plotGenome() rejects missing data/segments and a bad pos.unit", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_error(plotGenome())
  expect_error(plotGenome(segments = single.seg, pos.unit = "bogus"))
})
