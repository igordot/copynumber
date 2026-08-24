test_that("plotAllele() runs without error on valid input", {
  data(logR)
  data(BAF)
  logRWins <- winsorize(logR, verbose = FALSE)
  alleleSeg <- aspcf(logRWins, BAF, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(
    plotAllele(logR, BAF, alleleSeg, sample = 1, chrom = 1)
  )
})

test_that("plotAllele() rejects missing logR/BAF/segments and unpaired input", {
  pdf(NULL)
  on.exit(dev.off())
  expect_error(plotAllele())
  expect_error(plotAllele(logR = data.frame(a = 1, b = 2, c = 3)))
})
