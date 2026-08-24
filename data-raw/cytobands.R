# cytoband data for various genome builds

# download and clean UCSC cytoBand data for one genome build
prepare_cytoband <- function(build, chroms) {
  url <- paste0(
    "http://hgdownload.cse.ucsc.edu/goldenpath/",
    build,
    "/database/cytoBand.txt.gz"
  )
  cb <- data.table::fread(url, header = FALSE, data.table = FALSE)
  cb <- cb[cb$V1 %in% chroms, ]
  rownames(cb) <- seq_len(nrow(cb))
  cb$V1 <- factor(cb$V1)
  cb$V4 <- factor(cb$V4)
  cb$V5 <- factor(cb$V5)
  cb
}

human_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
mouse_chrs <- c(paste0("chr", 1:19), "chrX", "chrY")

# load existing sysdata.rda
original_sysdata <- load("R/sysdata.rda")

# fix mm7, mm8, and mm9
# originally consumed first cytoband row as a column header
# there's no p band for mouse builds, so row 1 is not needed
for (build in c("mm7", "mm8", "mm9")) {
  assign(build, prepare_cytoband(build, mouse_chrs))
}

# add hg38 and mm10 genome builds
hg38 <- prepare_cytoband("hg38", human_chrs)
mm10 <- prepare_cytoband("mm10", mouse_chrs)

# save updated sysdata.rda with new genome builds
save(list = c(original_sysdata, "hg38", "mm10"), file = "R/sysdata.rda")
