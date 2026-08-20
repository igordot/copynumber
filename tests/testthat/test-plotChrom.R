test_that("plotChrom() runs without error on valid input", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  multi.seg <- multipcf(data = lymph.wins, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(
    plotChrom(data = lymph.wins, segments = multi.seg, chrom = 1)
  )
})

test_that("plotChrom() rejects missing data/segments", {
  pdf(NULL)
  on.exit(dev.off())
  expect_error(plotChrom())
})
