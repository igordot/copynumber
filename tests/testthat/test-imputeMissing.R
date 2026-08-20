test_that("imputeMissing() replaces NAs with a constant", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  d <- sub.lymphoma[1:20, ]
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
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)

  expect_error(imputeMissing(data = sub.lymphoma, method = "bogus"))
})
