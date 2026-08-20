test_that("winsorize() flags outliers and preserves data shape", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)

  wins.res <- winsorize(
    data = sub.lymphoma,
    return.outliers = TRUE,
    verbose = FALSE
  )

  expect_named(wins.res, c("wins.data", "wins.outliers"))
  expect_equal(dim(wins.res$wins.data), dim(sub.lymphoma))
  expect_equal(dim(wins.res$wins.outliers), dim(sub.lymphoma))
  expect_true(all(unlist(wins.res$wins.outliers[, -c(1, 2)]) %in% c(-1, 0, 1)))
})
