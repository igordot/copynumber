test_that("subsetSegments() filters segments by chromosome", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  ss <- subsetSegments(segments = single.seg, chrom = 1)

  expect_s3_class(ss, "data.frame")
  expect_equal(unique(ss$chrom), 1)
  expect_true(nrow(ss) < nrow(single.seg))
  expect_true(nrow(ss) > 0)
})

test_that("subsetSegments() accepts a character matrix, not just a data frame", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  # as.matrix() on pcf() output coerces it to a character matrix.
  seg_mat <- as.matrix(single.seg)
  expect_true(is.character(seg_mat))

  ss <- subsetSegments(segments = seg_mat, chrom = 1)
  ss_from_df <- subsetSegments(segments = single.seg, chrom = 1)

  expect_s3_class(ss, "data.frame")
  # reparse: as.matrix() pads numbers for column alignment.
  tidy <- function(df) {
    as.data.frame(lapply(df, \(col) type.convert(trimws(as.character(col)), as.is = TRUE)))
  }
  expect_equal(tidy(ss), tidy(ss_from_df))
})
