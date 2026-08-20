test_that("selectSegments() picks the highest-variance multipcf segments", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  multi.seg <- multipcf(data = lymph.wins, verbose = FALSE)

  sel <- selectSegments(segments = multi.seg, what = "variance", nseg = 5)

  # Compute the "true" top-5 by variance independently of selectSegments(),
  # from the same 3 sample columns, and check the same segments come back
  # (identified by chrom + start.pos, since row order isn't guaranteed).
  sample_cols <- c("X01.B1", "X01.B2", "X01.B3")
  true_var <- apply(multi.seg[, sample_cols], 1, var)
  want <- order(true_var, decreasing = TRUE)[1:5]
  want_key <- paste(multi.seg$chrom[want], multi.seg$start.pos[want])
  got_key <- paste(sel$sel.seg$chrom, sel$sel.seg$start.pos)

  expect_equal(nrow(sel$sel.seg), 5)
  expect_setequal(got_key, want_key)
})

test_that("selectSegments() rejects segments not from multipcf", {
  data(lymphoma)
  sub.lymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymph.wins <- winsorize(data = sub.lymphoma, verbose = FALSE)
  single.seg <- pcf(data = lymph.wins, gamma = 12, verbose = FALSE)

  expect_error(selectSegments(segments = single.seg))
})
