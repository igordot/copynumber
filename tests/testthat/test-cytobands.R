test_that("bundled cytoband data is structurally sound for every genome build", {
  for (build in validAssemblies()) {
    cyto <- get(build)
    label <- paste("build:", build)

    expect_equal(ncol(cyto), 5, info = label)
    expect_false(anyNA(cyto), info = label)

    chrom <- as.character(cyto[[1]])
    start <- as.numeric(as.character(cyto[[2]]))
    end <- as.numeric(as.character(cyto[[3]]))
    band <- as.character(cyto[[4]])

    expect_true(all(grepl("^chr([0-9]+|X|Y)$", chrom)), info = label)
    expect_true(all(end > start), info = label)
    expect_true(all(grepl("^[pq]", band)), info = label)

    for (thisChrom in unique(chrom)) {
      idx <- chrom == thisChrom
      chromLabel <- paste(label, "chrom:", thisChrom)

      # coverage must start at the telomere (position 0)
      expect_equal(min(start[idx]), 0, info = chromLabel)

      # bands must tile the chromosome with no gaps or overlaps
      chromStart <- start[idx]
      chromEnd <- end[idx]
      if (length(chromStart) > 1) {
        expect_equal(
          chromStart[-1],
          chromEnd[-length(chromEnd)],
          info = chromLabel
        )
      }
    }
  }
})
