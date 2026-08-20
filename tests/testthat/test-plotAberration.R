test_that("plotAberration() runs without error on valid input", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plotAberration(segments = single.seg, thres.gain = 0.2))
})

test_that("plotAberration() rejects a bad pos.unit", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_error(
    plotAberration(segments = single.seg, thres.gain = 0.2, pos.unit = "bogus")
  )
})
