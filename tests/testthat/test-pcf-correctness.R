test_that("pcf() finds the one true break in a noiseless step function", {
  # This input has a flat block of 0s, then a flat block of 5s.
  # There is one correct segmentation: two 10-probe segments, means 0
  # and 5, split between probe 10 and probe 11. We know this before we
  # run pcf(). We do not read it from pcf()'s own output.
  d <- data.frame(chrom = 1, pos = 1:20, S1 = c(rep(0, 10), rep(5, 10)))

  res <- pcf(d, gamma = 10, normalize = FALSE, verbose = FALSE)

  expect_equal(nrow(res), 2)
  expect_equal(res$n.probes, c(10, 10))
  expect_equal(res$mean, c(0, 5))
  expect_equal(res$end.pos[1], 10)
  expect_equal(res$start.pos[2], 11)
})

test_that("pcf() default normalize=TRUE needs noise, not just a jump", {
  # normalize=TRUE scales the penalty by a MAD noise estimate. A
  # noiseless step has a MAD of exactly 0. So the penalty blocks every
  # split, at any gamma. Real measurements always have some noise. So
  # this case is real, but not a problem in practice. We test it here
  # so a future change to the normalization logic does not hide it.
  d <- data.frame(chrom = 1, pos = 1:20, S1 = c(rep(0, 10), rep(5, 10)))

  expect_equal(getMad(d$S1), 0)
  res <- pcf(d, gamma = 0.01, normalize = TRUE, verbose = FALSE)
  expect_equal(nrow(res), 1)
})
