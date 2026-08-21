test_that("winsorize()'s documented example produces stable output", {
  # Values pinned from an R 4.6.1 run, not derived independently.
  data(lymphoma)
  sub.lymphoma <- subsetData(lymphoma, sample = 1:3)

  wins.data <- winsorize(data = sub.lymphoma, verbose = FALSE)

  expect_equal(dim(wins.data), c(3091, 5))
  expect_equal(
    as.numeric(wins.data[1, c("X01.B1", "X01.B2", "X01.B3")]),
    c(0.0284, -0.0144, 0.0255),
    tolerance = 0.01
  )
})

test_that("winsorize() rejects an invalid assembly and names every valid build", {
  data(lymphoma)
  err <- expect_error(winsorize(lymphoma, assembly = "bogus", verbose = FALSE))
  for (build in validAssemblies()) {
    expect_true(grepl(build, conditionMessage(err), fixed = TRUE))
  }
})

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
