test_that("plotCircle() runs without error or warning on valid input", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plotCircle(segments = single.seg, thres.gain = 0.15))
  expect_no_warning(plotCircle(segments = single.seg, thres.gain = 0.15))
})
