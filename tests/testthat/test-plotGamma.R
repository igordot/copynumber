test_that("plotGamma() runs without error on valid input", {
  data(micma)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plotGamma(micma, chrom = 17))
})

test_that("plotGamma(cv = TRUE) returns a gamma grid that matches gammaRange", {
  data(micma)

  pdf(NULL)
  on.exit(dev.off())
  res <- plotGamma(
    micma,
    chrom = 17,
    gammaRange = c(10, 100),
    cv = TRUE,
    K = 2
  )

  # These are known ahead of running plotGamma(): 10 evenly spaced gamma
  # values across the given range, one pred.error per gamma, and the
  # reported optGamma must be one of the evaluated gamma values.
  expect_equal(res$gamma, seq(10, 100, length.out = 10))
  expect_equal(length(res$pred.error), length(res$gamma))
  expect_true(res$optGamma %in% res$gamma)
})
