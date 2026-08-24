test_that("plotChrom() runs without error on valid input", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  multiSeg <- multipcf(data = lymphWins, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(
    plotChrom(data = lymphWins, segments = multiSeg, chrom = 1)
  )
})

test_that("plotChrom() rejects missing data/segments", {
  pdf(NULL)
  on.exit(dev.off())
  expect_error(plotChrom())
})

test_that("getSeglim() rejects equalRange = FALSE with no k or sampleID", {
  # Internal helper; every current caller always supplies k or sampleID, so tested directly.
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  expect_error(getSeglim(singleSeg, equalRange = FALSE))
})
