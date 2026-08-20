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
