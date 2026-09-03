test_that("multipcf()'s documented example produces stable output", {
  # Values pinned from an R 4.6.1 run, not derived independently.
  data(lymphoma)
  subLymphoma <- subsetData(lymphoma, sample = 1:3)
  winsLymph <- winsorize(subLymphoma, verbose = FALSE)

  multiSegments <- multipcf(
    data = winsLymph,
    gamma = 12,
    Y = subLymphoma,
    verbose = FALSE
  )

  expect_equal(nrow(multiSegments), 54)
  expect_equal(sum(multiSegments$n.probes), 3091)
  expect_equal(
    as.numeric(multiSegments[1, c("X01.B1", "X01.B2", "X01.B3")]),
    c(-0.0439, -0.0324, -0.0643),
    tolerance = 0.01
  )
})

test_that("multipcf() rejects an invalid assembly and names every valid build", {
  data(lymphoma)
  err <- expect_error(multipcf(lymphoma, assembly = "bogus", verbose = FALSE))
  for (build in validAssemblies()) {
    expect_true(grepl(build, conditionMessage(err), fixed = TRUE))
  }
})

test_that("multipcf() reproduces the vignette's lymphoma workflow", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)

  multiSeg <- multipcf(data = lymphWins, verbose = FALSE)

  # multipcf finds common breakpoints, so segments are wide: one row per
  # segment with a value column per sample, rather than pcf's long format.
  sample_cols <- c("X01.B1", "X01.B2", "X01.B3")
  expect_named(
    multiSeg,
    c("chrom", "arm", "start.pos", "end.pos", "n.probes", sample_cols)
  )
  # a real invariant: segments partition every probe in the input exactly
  # once (breakpoints are shared across samples, so there is one n.probes
  # column, not one per sample).
  expect_equal(sum(multiSeg$n.probes), nrow(lymphWins))
})

test_that("multipcf() handles a chromosome arm with a single probe", {
  # Regression test: a 1-probe arm used to crash inside doMultiPCF().
  data <- data.frame(
    chrom = c(1, 1, 2),
    pos = c(1000, 2000, 1000),
    sample1 = c(0.1, 0.15, 0.5),
    sample2 = c(0.2, 0.22, 0.6)
  )
  arms <- c("p", "p", "p")

  segments <- multipcf(data = data, arms = arms, verbose = FALSE)

  expect_equal(sum(segments$n.probes), nrow(data))
  expect_true(1 %in% segments$n.probes)
})

test_that("multipcf() handles a single-probe arm when a sample has zero variance", {
  # Regression test: a constant sample column forces the sd==0 short-circuit
  # branch, which used to corrupt output shape for a 1-probe arm.
  data <- data.frame(
    chrom = c(1, 1, 2),
    pos = c(1000, 2000, 1000),
    sample1 = c(0.5, 0.5, 0.5),
    sample2 = c(0.2, 0.22, 0.6)
  )
  arms <- c("p", "p", "p")

  result <- multipcf(
    data = data,
    arms = arms,
    return.est = TRUE,
    verbose = FALSE
  )

  expect_equal(sum(result$segments$n.probes), nrow(data))
  expect_equal(nrow(result$estimates), nrow(data))
})
