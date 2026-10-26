test_that("plotCircle() runs without error or warning on valid input", {
  data(lymphoma)
  subLymphoma <- subsetData(data = lymphoma, sample = 1:3)
  lymphWins <- winsorize(data = subLymphoma, verbose = FALSE)
  singleSeg <- pcf(data = lymphWins, gamma = 12, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plotCircle(segments = singleSeg, thres.gain = 0.15))
  expect_no_warning(plotCircle(segments = singleSeg, thres.gain = 0.15))
})

test_that("mouse global positions include X and Y with one gap per chromosome", {
  cyto <- data.frame(
    chrom = c("chrY", "chr1", "chrX"),
    start = 0,
    end = c(30, 10, 20),
    band = "q1",
    stain = "gneg"
  )

  expect_equal(
    getGlobPos(c(1, 23, 24, 23), c(1, 1, 1, 25), "bp", cyto, delta = 3),
    c(1, 14, 37, 33)
  )
})

test_that("plotCircle() supports mouse assemblies and sex chromosome arcs", {
  segments <- data.frame(
    sampleID = "S1",
    chrom = c("1", "19", "X", "Y"),
    arm = "q",
    start.pos = 1,
    end.pos = 100,
    n.probes = 10,
    mean = c(0.3, -0.3, 0.4, -0.4)
  )
  arcs <- matrix(c(19, 50, 23, 50, 1, 23, 50, 24, 50, 1), ncol = 5, byrow = TRUE)

  pdf(NULL)
  on.exit(dev.off())
  for (assembly in c("mm7", "mm8", "mm9", "mm10", "mm39")) {
    expect_no_error(expect_no_warning(plotCircle(
      segments,
      thres.gain = 0.2,
      arcs = arcs,
      assembly = assembly
    )))
  }
})
