test_that("callAberrations() classifies pcf (single-sample) segments", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  called <- callAberrations(segments = singleSeg, thres.gain = 0.2)

  expect_equal(colnames(called), c(colnames(singleSeg)[-7], "call"))
  expect_true(all(called$call %in% c("normal", "gain", "loss")))
  expect_true(all(called$call[called$mean > 0.2] == "gain"))
  expect_true(all(called$call[called$mean < -0.2] == "loss"))
})

test_that("callAberrations() classifies multipcf (wide) segments per sample", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  multiSeg <- multipcf(data = lymphWins, verbose = FALSE)

  called <- callAberrations(segments = multiSeg, thres.gain = 0.2)

  sample_cols <- c("X01.B1", "X01.B2", "X01.B3")
  expect_true(all(unlist(called[sample_cols]) %in% c("normal", "gain", "loss")))
  expect_equal(nrow(called), nrow(multiSeg))
})

test_that("callAberrations() keeps missing pcf segment values unknown", {
  segments <- data.frame(
    sampleID = "S1",
    chrom = 1,
    arm = "p",
    start.pos = 1:7,
    end.pos = 1:7,
    n.probes = 1,
    mean = c(NA_real_, NaN, -0.3, -0.2, 0, 0.2, 0.3)
  )

  called <- callAberrations(segments, thres.gain = 0.2)

  expect_identical(
    called$call,
    c(
      NA_character_,
      NA_character_,
      "loss",
      "normal",
      "normal",
      "normal",
      "gain"
    )
  )
  expect_equal(called[1:6], segments[1:6])
})

test_that("callAberrations() keeps missing multipcf values unknown per sample", {
  segments <- data.frame(
    chrom = 1,
    arm = "p",
    start.pos = 1:7,
    end.pos = 1:7,
    n.probes = 1,
    S1 = c(NA_real_, NaN, -0.3, -0.2, 0, 0.2, 0.3),
    S2 = c(0.3, 0.2, 0, -0.2, -0.3, NaN, NA_real_)
  )

  called <- callAberrations(segments, thres.gain = 0.2)

  expect_identical(
    called$S1,
    c(
      NA_character_,
      NA_character_,
      "loss",
      "normal",
      "normal",
      "normal",
      "gain"
    )
  )
  expect_identical(
    called$S2,
    c(
      "gain",
      "normal",
      "normal",
      "normal",
      "loss",
      NA_character_,
      NA_character_
    )
  )
  expect_equal(called[1:5], segments[1:5])
})
