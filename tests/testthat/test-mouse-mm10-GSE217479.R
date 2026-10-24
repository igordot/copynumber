test_that("the mouse fixture has finite values at ordered mm10 positions", {
  mouse <- read.delim(test_path("fixtures", "mouse-mm10-GSE217479-sparse.tsv"))

  expect_named(mouse, c("chrom", "pos", "WT4", "WT6"))
  expect_equal(dim(mouse), c(2032L, 4L))
  expect_identical(unique(mouse$chrom), c(as.character(1:19), "X", "Y"))
  expect_equal(anyDuplicated(mouse[c("chrom", "pos")]), 0L)
  expect_equal(sum(!is.finite(as.matrix(mouse[-1]))), 0L)
  for (chrom in unique(mouse$chrom)) {
    positions <- mouse$pos[mouse$chrom == chrom]
    bands <- mm10[as.character(mm10[[1]]) == paste0("chr", chrom), ]
    expect_equal(length(positions), if (chrom == "Y") 32L else 100L)
    expect_gt(min(diff(positions)), 0)
    expect_gt(min(positions), 0)
    expect_lte(max(positions), max(as.numeric(as.character(bands[[3]]))))
  }
})

test_that("the dense mouse fixture keeps all usable bins on chromosomes 18 and 19", {
  mouse <- read.delim(test_path("fixtures", "mouse-mm10-GSE217479-dense.tsv"))

  expect_named(mouse, c("chrom", "pos", "WT4", "WT6"))
  expect_equal(dim(mouse), c(13909L, 4L))
  expect_equal(as.integer(table(mouse$chrom)), c(6126L, 7783L))
  expect_identical(unique(mouse$chrom), c(18L, 19L))
  expect_equal(anyDuplicated(mouse[c("chrom", "pos")]), 0L)
  expect_equal(sum(!is.finite(as.matrix(mouse))), 0L)
  for (chrom in unique(mouse$chrom)) {
    expect_gt(min(diff(mouse$pos[mouse$chrom == chrom])), 0)
  }
})

for (fixture in c(
  "mouse-mm10-GSE217479-sparse.tsv",
  "mouse-mm10-GSE217479-dense.tsv"
)) {
  test_that(
    paste(fixture, "supports the mm10 analysis and plotting workflow"),
    {
      mouse <- read.delim(test_path("fixtures", fixture))
      wins <- winsorize(mouse, assembly = "mm10", verbose = FALSE)
      expect_identical(wins[1:2], mouse[1:2])
      expect_equal(dim(wins), dim(mouse))
      expect_equal(sum(!is.finite(as.matrix(wins[-c(1, 2)]))), 0L)

      single <- pcf(wins, assembly = "mm10", verbose = FALSE)
      joint <- multipcf(wins, assembly = "mm10", verbose = FALSE)
      # Values from the first mm10 commit, ef99c55 (2020-07-15).
      if (fixture == "mouse-mm10-GSE217479-sparse.tsv") {
        expect_equal(dim(single), c(50L, 7L))
        expect_equal(dim(joint), c(30L, 7L))
        expect_equal(wins$WT4[c(1, 1016, 2032)], c(-0.4985, -0.0982, 1.1289))
        expect_equal(wins$WT6[c(1, 1016, 2032)], c(-0.1664, -0.1434, 1.5674))
        expect_equal(single$mean[c(1, 25, 50)], c(-0.1221, -0.2116, 1.4691))
        expect_equal(single$end.pos[c(1, 25)], c(195176582, 120018546))
        expect_equal(joint$WT4[c(1, 30)], c(-0.1221, 1.1289))
        expect_equal(joint$WT6[c(1, 30)], c(-0.0497, 1.5674))
      } else {
        expect_equal(dim(single), c(51L, 7L))
        expect_equal(dim(joint), c(26L, 7L))
        expect_equal(wins$WT4[c(1, 6955, 13909)], c(0.6093, 0.3501, 0.9679))
        expect_equal(wins$WT6[c(1, 6955, 13909)], c(0.6999, 0.9729, 1.29))
        expect_equal(single$mean[c(1, 26, 51)], c(0.6584, 0.5712, 0.6968))
        expect_equal(single$end.pos[c(1, 26)], c(3267660, 12322708))
        expect_equal(joint$WT4[c(1, 26)], c(0.901, 0.2827))
        expect_equal(joint$WT6[c(1, 26)], c(1.032, 0.762))
      }
      expect_setequal(single$sampleID, c("WT4", "WT6"))
      expect_named(
        joint,
        c("chrom", "arm", "start.pos", "end.pos", "n.probes", "WT4", "WT6")
      )
      for (sample in c("WT4", "WT6")) {
        for (segments in list(single[single$sampleID == sample, ], joint)) {
          coverage <- integer(nrow(mouse))
          for (i in seq_len(nrow(segments))) {
            segment <- segments[i, ]
            inside <- mouse$chrom == segment$chrom &
              mouse$pos >= segment$start.pos &
              mouse$pos <= segment$end.pos
            coverage[inside] <- coverage[inside] + 1L
            expect_equal(sum(inside), segment$n.probes)
            value <- if ("mean" %in% names(segment)) {
              segment$mean
            } else {
              segment[[sample]]
            }
            expect_lte(abs(value - mean(wins[[sample]][inside])), 0.00005001)
          }
          expect_identical(coverage, rep(1L, nrow(mouse)))
        }
      }

      pdf(NULL)
      on.exit(dev.off())
      expect_no_error(plotGenome(
        mouse,
        segments = single,
        assembly = "mm10",
        sample = 1
      ))
    }
  )
}
