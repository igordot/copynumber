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

test_that("getArms() never assigns a p-arm boundary from a different chromosome", {
  for (build in validAssemblies()) {
    cyto <- get(build)
    label <- paste("build:", build)

    chrom <- as.character(cyto[[1]])
    end <- as.numeric(as.character(cyto[[3]]))
    armChar <- substring(as.character(cyto[[4]]), 1, 1)

    for (thisChrom in unique(chrom)) {
      idx <- chrom == thisChrom
      chromNum <- numericChrom(sub("^chr", "", thisChrom))
      chromLabel <- paste(label, "chrom:", thisChrom)
      hasP <- any(armChar[idx] == "p")
      chromEnd <- max(end[idx])

      # a position at the telomere must be "q"
      startArm <- getArms(chromNum, 1, pos.unit = "bp", cyto.data = cyto)
      expect_equal(startArm, if (hasP) "p" else "q", info = chromLabel)

      # a position at the far end of the chromosome must always be "q"
      endArm <- getArms(
        chromNum,
        chromEnd - 1,
        pos.unit = "bp",
        cyto.data = cyto
      )
      expect_equal(endArm, "q", info = chromLabel)
    }
  }
})

test_that("getArmandChromStop() reproduces hg19's known p-arm/chromosome stops", {
  l <- getArmandChromStop(hg19, "bp")

  expect_equal(
    l$pstop,
    c(
      125000000,
      93300000,
      91000000,
      50400000,
      48400000,
      61000000,
      59900000,
      45600000,
      49000000,
      40200000,
      53700000,
      35800000,
      17900000,
      17600000,
      19000000,
      36600000,
      24000000,
      17200000,
      26500000,
      27500000,
      13200000,
      14700000,
      60600000,
      12500000
    ),
    tolerance = 1e-6
  )
  expect_equal(
    l$chromstop,
    c(
      249250621,
      243199373,
      198022430,
      191154276,
      180915260,
      171115067,
      159138663,
      146364022,
      141213431,
      135534747,
      135006516,
      133851895,
      115169878,
      107349540,
      102531392,
      90354753,
      81195210,
      78077248,
      59128983,
      63025520,
      48129895,
      51304566,
      155270560,
      59373566
    )
  )
})
