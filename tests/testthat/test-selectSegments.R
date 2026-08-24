test_that("selectSegments() picks the highest-variance multipcf segments", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:5)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  multiSeg <- multipcf(data = lymphWins, verbose = FALSE)

  sel <- selectSegments(segments = multiSeg, what = "variance", nseg = 5)

  # Independent top-5-by-variance check, matched by chrom+start.pos (row order isn't guaranteed).
  sample_cols <- c("X01.B1", "X01.B2", "X01.B3", "X03.B1", "X03.B2")
  true_var <- apply(multiSeg[, sample_cols], 1, var)
  want <- order(true_var, decreasing = TRUE)[1:5]
  want_key <- paste(multiSeg$chrom[want], multiSeg$start.pos[want])
  got_key <- paste(sel$sel.seg$chrom, sel$sel.seg$start.pos)

  expect_equal(nrow(sel$sel.seg), 5)
  expect_setequal(got_key, want_key)
})

test_that("selectSegments() picks low-variance segments below a threshold", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:5)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  multiSeg <- multipcf(data = lymphWins, verbose = FALSE)

  # thres picked by inspecting real output, not derived independently.
  sel <- selectSegments(
    segments = multiSeg,
    what = "variance",
    thres = 0.001,
    large = FALSE
  )

  # Independently compute which segments truly have variance below thres.
  sample_cols <- c("X01.B1", "X01.B2", "X01.B3", "X03.B1", "X03.B2")
  true_var <- apply(multiSeg[, sample_cols], 1, var)
  want <- which(true_var < 0.001)
  want_key <- paste(multiSeg$chrom[want], multiSeg$start.pos[want])
  got_key <- paste(sel$sel.seg$chrom, sel$sel.seg$start.pos)

  expect_setequal(got_key, want_key)
})

test_that("selectSegments() picks the longest multipcf segments", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:5)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  multiSeg <- multipcf(data = lymphWins, verbose = FALSE)

  sel <- selectSegments(segments = multiSeg, what = "length", nseg = 5)

  # Independently compute the true top-5 longest segments.
  true_length <- multiSeg$end.pos - multiSeg$start.pos + 1
  want <- order(true_length, decreasing = TRUE)[1:5]
  want_key <- paste(multiSeg$chrom[want], multiSeg$start.pos[want])
  got_key <- paste(sel$sel.seg$chrom, sel$sel.seg$start.pos)

  expect_equal(nrow(sel$sel.seg), 5)
  expect_setequal(got_key, want_key)
})

test_that("selectSegments() picks segments aberrant in a minimum proportion of samples", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:5)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  multiSeg <- multipcf(data = lymphWins, verbose = FALSE)

  # thres/p picked by inspecting real output, not derived independently.
  sel <- selectSegments(
    segments = multiSeg,
    what = "aberration",
    thres = 0.2,
    p = 0.2
  )

  # Independently compute which segments truly clear the proportion cutoff.
  sample_cols <- c("X01.B1", "X01.B2", "X01.B3", "X03.B1", "X03.B2")
  true_prop <- rowMeans(multiSeg[, sample_cols] > 0.2)
  want <- which(true_prop >= 0.2)
  want_key <- paste(multiSeg$chrom[want], multiSeg$start.pos[want])
  got_key <- paste(sel$sel.seg$chrom, sel$sel.seg$start.pos)

  expect_setequal(got_key, want_key)
})

test_that("selectSegments() picks segments with the highest aberration quantile", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:5)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  multiSeg <- multipcf(data = lymphWins, verbose = FALSE)

  sel <- selectSegments(
    segments = multiSeg,
    what = "aberration",
    nseg = 20,
    p = 0.5
  )

  # Independently compute the true top-20 by (1-p)-quantile across samples.
  sample_cols <- c("X01.B1", "X01.B2", "X01.B3", "X03.B1", "X03.B2")
  true_q <- apply(multiSeg[, sample_cols], 1, quantile, probs = 0.5, type = 1)
  want <- order(true_q, decreasing = TRUE)[1:20]
  want_key <- paste(multiSeg$chrom[want], multiSeg$start.pos[want])
  got_key <- paste(sel$sel.seg$chrom, sel$sel.seg$start.pos)

  expect_equal(nrow(sel$sel.seg), 20)
  expect_setequal(got_key, want_key)
})

test_that("selectSegments() returns a list when nseg exceeds available segments", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:5)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  multiSeg <- multipcf(data = lymphWins, verbose = FALSE)

  expect_warning(
    sel <- selectSegments(segments = multiSeg, what = "variance", nseg = 10^6)
  )

  expect_named(sel, c("sel.seg", "seg.var"))
  expect_equal(nrow(sel$sel.seg), nrow(multiSeg))
})

test_that("selectSegments() returns 0 rows, not a fabricated row, for 0-row input", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:5)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  multiSeg <- multipcf(data = lymphWins, verbose = FALSE)
  emptySeg <- multiSeg[0, ]

  expect_warning(
    sel <- selectSegments(segments = emptySeg, what = "variance")
  )

  expect_equal(nrow(sel$sel.seg), 0)
})

test_that("selectSegments() rejects segments not from multipcf", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:5)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  expect_error(selectSegments(segments = singleSeg))
})
