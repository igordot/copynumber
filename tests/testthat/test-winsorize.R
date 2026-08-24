test_that("winsorize()'s documented example produces stable output", {
  # Values pinned from an R 4.6.1 run, not derived independently.
  data(lymphoma)
  subLymphoma <- subsetData(lymphoma, sample = 1:3)

  winsData <- winsorize(data = subLymphoma, verbose = FALSE)

  expect_equal(dim(winsData), c(3091, 5))
  expect_equal(
    as.numeric(winsData[1, c("X01.B1", "X01.B2", "X01.B3")]),
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
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)

  winsRes <- winsorize(
    data = subLymphoma,
    return.outliers = TRUE,
    verbose = FALSE
  )

  expect_named(winsRes, c("wins.data", "wins.outliers"))
  expect_equal(dim(winsRes$wins.data), dim(subLymphoma))
  expect_equal(dim(winsRes$wins.outliers), dim(subLymphoma))
  expect_true(all(unlist(winsRes$wins.outliers[, -c(1, 2)]) %in% c(-1, 0, 1)))
})
