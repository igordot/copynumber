test_that("imputeMissing() replaces NAs with a constant", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  d <- subLymphoma[1:20, ]
  d[3, 3] <- NA
  d[7, 4] <- NA

  imp <- imputeMissing(data = d, method = "constant", c = 0)

  expect_equal(sum(is.na(imp)), 0)
  expect_equal(imp[3, 3], 0)
  expect_equal(imp[7, 4], 0)
  # untouched cells are unaffected
  expect_equal(imp[-c(3, 7), ], d[-c(3, 7), ])
})

test_that("imputeMissing() rejects an unknown method", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)

  expect_error(imputeMissing(data = subLymphoma, method = "bogus"))
})

test_that("imputeMissing() leaves complete data unchanged with the pcf method", {
  d <- data.frame(chrom = 1, pos = 1:10, S1 = rep(c(0.1, 0.2), 5), S2 = 0.3)

  expect_identical(imputeMissing(d, method = "pcf"), d)
})

test_that("imputeMissing() rejects missing coordinates with the pcf method", {
  d <- data.frame(chrom = 1, pos = 1:10, S1 = 0.2)
  for (column in 1:2) {
    invalid <- d
    invalid[1, column] <- NA
    expect_error(
      imputeMissing(invalid, method = "pcf"),
      "missing values are not allowed in data columns 1 and 2",
      fixed = TRUE
    )
  }
})

test_that("imputeMissing() fills missing samples with pcf estimates", {
  d <- data.frame(chrom = 1, pos = 1:10, S1 = 0.2, S2 = rep(c(0.1, 0.3), 5))
  expected <- d
  d[5, 3] <- NA

  expect_equal(
    imputeMissing(d, method = "pcf", normalize = FALSE, verbose = FALSE),
    expected
  )
})
