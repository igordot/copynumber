test_that("pcf()'s documented example produces stable output", {
  # Values pinned from an R 4.6.1 run, not derived independently.
  data(lymphoma)
  sub.lymphoma <- subsetData(lymphoma, sample = 1:3)
  wins.lymph <- winsorize(sub.lymphoma, verbose = FALSE)

  pcf.segments <- pcf(
    data = wins.lymph,
    gamma = 12,
    Y = sub.lymphoma,
    verbose = FALSE
  )

  expect_equal(nrow(pcf.segments), 175)
  expect_equal(sum(pcf.segments$n.probes), 9273)
  expect_equal(mean(pcf.segments$mean), 0.004773714, tolerance = 1e-6)
  expect_equal(pcf.segments$mean[1], -0.0439)
})

test_that("pcf() segments a single sample into the expected columns", {
  data(lymphoma)
  cols <- c("Chrom", "Median.bp", "X01.B1")
  sub <- lymphoma[lymphoma$Chrom == 1, cols][1:200, ]
  colnames(sub) <- c("chrom", "pos", "S1")

  res <- pcf(sub, gamma = 40, verbose = FALSE)

  expect_s3_class(res, "data.frame")
  expect_named(
    res,
    c("sampleID", "chrom", "arm", "start.pos", "end.pos", "n.probes", "mean")
  )
  expect_true(nrow(res) > 0)
  expect_equal(sum(res$n.probes), nrow(sub))
})

test_that("pcf() accepts matrix input, not just data frames", {
  data(lymphoma)
  cols <- c("Chrom", "Median.bp", "X01.B1")
  sub <- lymphoma[lymphoma$Chrom == 1, cols][1:200, ]
  colnames(sub) <- c("chrom", "pos", "S1")

  res_df <- pcf(sub, gamma = 40, verbose = FALSE)
  res_mat <- pcf(as.matrix(sub), gamma = 40, verbose = FALSE)

  # row.names differ (matrix dimnames carry through); values must not.
  expect_equal(res_mat, res_df, ignore_attr = "row.names")
})

test_that("pcf() rejects an invalid assembly and names every valid build", {
  data(lymphoma)
  err <- expect_error(pcf(lymphoma, assembly = "bogus", verbose = FALSE))
  for (build in validAssemblies()) {
    expect_true(grepl(build, conditionMessage(err), fixed = TRUE))
  }
})

test_that("pcf() reproduces the vignette's lymphoma workflow", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)

  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  ids <- c("X01.B1", "X01.B2", "X01.B3")
  expect_equal(sort(unique(single.seg$sampleID)), ids)
  # a real invariant, known without running pcf(): segments partition every
  # probe in the input exactly once, per sample.
  for (id in ids) {
    n <- sum(single.seg$n.probes[single.seg$sampleID == id])
    expect_equal(n, nrow(lymph.wins))
  }
})
