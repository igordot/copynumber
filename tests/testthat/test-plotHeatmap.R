test_that("plotHeatmap() runs without error on valid input", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plotHeatmap(segments = single.seg, upper.lim = 0.3))
})

test_that("plotHeatmap() rejects a bad pos.unit and a non-positive upper.lim", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_error(
    plotHeatmap(segments = single.seg, upper.lim = 0.3, pos.unit = "bogus")
  )
  expect_error(plotHeatmap(segments = single.seg, upper.lim = -1))
})
