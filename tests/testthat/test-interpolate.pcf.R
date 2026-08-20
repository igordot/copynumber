test_that("interpolate.pcf() reads off segment means at arbitrary positions", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  x <- sub.lymphoma[sub.lymphoma[, 1] == 1, 1:2][1:5, ]
  interp <- interpolate.pcf(segments = single.seg, x = x)

  expect_named(interp, c("chr", "pos", unique(single.seg$sampleID)))
  expect_equal(nrow(interp), nrow(x))

  # every interpolated value should equal the mean of some chrom-1 segment
  is_chrom1 <- single.seg$chrom == 1 & single.seg$sampleID == "X01.B1"
  chrom1 <- single.seg[is_chrom1, ]
  expect_true(all(interp$X01.B1 %in% chrom1$mean))
})
