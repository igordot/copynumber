test_that("interpolate.pcf() reads off segment means at arbitrary positions", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  x <- subLymphoma[subLymphoma[, 1] == 1, 1:2][1:5, ]
  interp <- interpolate.pcf(segments = singleSeg, x = x)

  expect_named(interp, c("chr", "pos", unique(singleSeg$sampleID)))
  expect_equal(nrow(interp), nrow(x))

  # every interpolated value should equal the mean of some chrom-1 segment
  is_chrom1 <- singleSeg$chrom == 1 & singleSeg$sampleID == "X01.B1"
  chrom1 <- singleSeg[is_chrom1, ]
  expect_true(all(interp$X01.B1 %in% chrom1$mean))
})

test_that("interpolate.pcf() returns NA, not an error, for a chromosome absent from segments", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  x <- cbind(chr = c(99, 98), pos = c(1000, 2000))
  interp <- expect_no_error(interpolate.pcf(segments = singleSeg, x = x))

  expect_true(all(is.na(interp[, -c(1, 2)])))
})

test_that("interpolate.pcf() accepts a single-row matrix x", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  x <- cbind(chr = 1, pos = 2000000)
  interp <- expect_no_error(interpolate.pcf(segments = singleSeg, x = x))

  expect_equal(nrow(interp), 1)
  expect_named(interp, c("chr", "pos", unique(singleSeg$sampleID)))
})
