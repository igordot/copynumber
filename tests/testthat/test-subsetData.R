test_that("subsetData() filters by chrom and by sample", {
  data(lymphoma)

  by_chrom <- subsetData(data = lymphoma, chrom = c(1, 2))
  expect_equal(unique(by_chrom[, 1]), c(1, 2))
  expect_equal(ncol(by_chrom), ncol(lymphoma))

  by_sample <- subsetData(data = lymphoma, sample = 1:3)
  expect_equal(ncol(by_sample), 2 + 3)
  expect_equal(nrow(by_sample), nrow(lymphoma))
})

test_that("subsetData() keeps file values paired with selected sample names", {
  d <- data.frame(
    chrom = c(1, 2, 1),
    pos = c(10, 20, 30),
    S1 = 101:103,
    S2 = 201:203,
    S3 = 301:303
  )
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path))
  write.table(d, path, sep = "\t", row.names = FALSE, quote = FALSE)

  selections <- list(NULL, c(1, 3), c(3, 1), 3, c(3, 1, 3), integer())
  for (sample in selections) {
    selected <- if (is.null(sample)) 1:3 else sample
    expected <- d[c(1, 3), c(1:2, selected + 2), drop = FALSE]
    names(expected) <- c("chrom", "pos", names(d)[selected + 2])
    rownames(expected) <- NULL

    expect_equal(subsetData(path, chrom = 1, sample = sample), expected)
  }
})
