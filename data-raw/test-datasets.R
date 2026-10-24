# Run from the package root to make the mouse test fixtures.
# Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217479
# CNVkit CNR files contain bin log2 ratios before segmentation.

samples <- c(WT4 = "GSM6720269", WT6 = "GSM6720270")
chromosomes <- c(as.character(1:19), "X", "Y")
tracks <- lapply(names(samples), function(sample) {
  accession <- samples[[sample]]
  filename <- paste0(accession, "_", sample, ".mapped.sorted.cnr.txt.gz")
  path <- file.path(tempdir(), filename)
  if (!file.exists(path)) {
    url <- paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6720nnn/",
      accession,
      "/suppl/",
      filename
    )
    download.file(url, path, mode = "wb")
  }
  bins <- read.delim(gzfile(path))
  stopifnot(all(
    c(
      "chromosome",
      "start",
      "end",
      "log2",
      "depth",
      "weight"
    ) %in%
      names(bins)
  ))
  bins <- bins[
    bins$chromosome %in%
      paste0("chr", chromosomes) &
      is.finite(bins$log2) &
      bins$depth > 0 &
      bins$weight > 0,
    c("chromosome", "start", "end", "log2")
  ]
  stopifnot(!anyNA(bins), all(bins$end > bins$start))
  names(bins)[4] <- sample
  bins
})

bins <- merge(tracks[[1]], tracks[[2]], by = c("chromosome", "start", "end"))
bins <- bins[
  order(match(bins$chromosome, paste0("chr", chromosomes)), bins$start),
]
selected <- unlist(
  lapply(paste0("chr", chromosomes), function(chromosome) {
    indices <- which(bins$chromosome == chromosome)
    stopifnot(length(indices) > 0)
    indices[round(seq(
      1,
      length(indices),
      length.out = min(100, length(indices))
    ))]
  }),
  use.names = FALSE
)

# CNR intervals use zero-based starts and exclusive ends.
mouse <- data.frame(
  chrom = sub("^chr", "", bins$chromosome),
  pos = floor((bins$start + bins$end) / 2) + 1,
  WT4 = bins$WT4,
  WT6 = bins$WT6
)
stopifnot(length(selected) == 2032, !anyDuplicated(mouse[c("chrom", "pos")]))
dir.create("tests/testthat/fixtures", showWarnings = FALSE)
write.table(
  mouse[selected, ],
  "tests/testthat/fixtures/mouse-mm10-GSE217479-sparse.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Chromosomes 18 and 19 are the two smallest mm10 autosomes.
dense <- mouse[mouse$chrom %in% c("18", "19"), ]
stopifnot(nrow(dense) == 13909)
write.table(
  dense,
  "tests/testthat/fixtures/mouse-mm10-GSE217479-dense.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
