test_that("callAberrations() classifies pcf (single-sample) segments", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  called <- callAberrations(segments = single.seg, thres.gain = 0.2)

  expect_equal(colnames(called), c(colnames(single.seg)[-7], "call"))
  expect_true(all(called$call %in% c("normal", "gain", "loss")))
  expect_true(all(called$call[called$mean > 0.2] == "gain"))
  expect_true(all(called$call[called$mean < -0.2] == "loss"))
})

test_that("callAberrations() classifies multipcf (wide) segments per sample", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  multi.seg <- multipcf(data = lymph.wins, verbose = FALSE)

  called <- callAberrations(segments = multi.seg, thres.gain = 0.2)

  sample_cols <- c("X01.B1", "X01.B2", "X01.B3")
  expect_true(all(unlist(called[sample_cols]) %in% c("normal", "gain", "loss")))
  expect_equal(nrow(called), nrow(multi.seg))
})
