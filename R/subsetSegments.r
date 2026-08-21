####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# Get a subset of segments for given chromosomes and/or sampleIDs

## Input:
### segments : either a data frame or the name of a file containing the segments found by either pcf, multipcf or aspcf
### chrom : a numeric vector with chromosome number(s) for which data or segments should be selected. If unspecified, all chromosomes in data or segments will be selected
### sample : a numeric vector indicating which sample(s) is to selected. The number(s) should correspond to the sample's place (in order of appearance) in the data
### sep : the separator of the input files if \code{data} is a file. Default is "\t"
### ... : optional parameters to be passed to \code{read.table} in the case where data are to be read from files

## Output:
### sel.segments : the data frame with segments selected on chromosomes and sampleIds

## Required by:
### none

## Requires:
### is.multiseg
### numericChrom
### pullOutContent



#' Retrieve a subset of segments
#' 
#' This function returns a subset of segments according to the input and the
#' specified chromosomes and/or samples.
#' 
#' 
#' @param segments either a data frame or the name of a tab-separated file from
#' which segmentation results can be read. Segmentation results may come from
#' \code{\link{pcf}}, \code{\link{multipcf}} or \code{\link{aspcf}}.
#' @param chrom a numeric or character vector with chromosome(s) for which
#' segments should be selected. If unspecified, all chromosomes in
#' \code{segments} will be selected.
#' @param sample a numeric vector indicating for which sample(s) segments
#' should be selected. The number(s) should correspond to the sample's place
#' (in order of appearance) in \code{segments}.
#' @param sep the separator of the input files if \code{segments} is a file.
#' Default is tab.
#' @param \dots optional parameters to be passed to \code{read.table} in the
#' case where \code{segments} is to be read from file.
#' @return A data frame containing the desired subset of segments.
#' @author Gro Nilsen
#' @examples
#' 
#' #Load lymphoma data
#' data(lymphoma)
#' 
#' #Select segments only for samples 1 and 6 and chromosomes 1:9:
#' segments <- pcf(lymphoma,gamma=12)
#' sub.segments <- subsetSegments(segments=segments,chrom=c(1:9),sample=c(1,6))
#' 
#' 
#' @export
subsetSegments <- function(
  segments,
  chrom = NULL,
  sample = NULL,
  sep = "\t",
  ...
) {
  # Check if segments is a file:
  isfile <- is.character(segments) && is.vector(segments)

  # get header and chrom from data
  if (isfile) {
    # read segment-file
    segments <- read.table(segments, header = TRUE, sep = sep, as.is = TRUE)
  } else {
    # Make sure segments is a data frame
    segments <- pullOutContent(res = segments, what = "segments")
  }

  # Check whether we have multiPCF segments:
  multi <- is.multiseg(segments)

  if (multi) {
    seg.sampleid <- colnames(segments)[-c(1:5)]
    seg.chrom <- segments[, 1]
  } else {
    seg.sampleid <- segments[, 1]
    seg.chrom <- segments[, 2]
  }

  # make sure chrom in segments are numeric
  seg.chrom <- numericChrom(seg.chrom)

  # Pick out relevant chromosomes:
  keepchrom <- 1:length(seg.chrom)
  if (!is.null(chrom)) {
    # make sure selected chrom are numeric:
    chrom <- numericChrom(chrom)
    # Check that these chrom are found in data.chrom
    use <- chrom %in% unique(seg.chrom)
    use.chrom <- chrom[use]
    if (length(use.chrom) == 0) {
      msg <- "None of the specified chromosomes are found in segments"
      stop(msg, call. = FALSE)
    } else if (length(use.chrom) != length(chrom)) {
      not.use <- paste(chrom[!use], sep = "", collapse = ",")
      msg <- paste(
        "The following chromosome(s) are not found in segments:",
        not.use,
        sep = " "
      )
      warning(msg, call. = FALSE, immediate. = TRUE)
    }
    keepchrom <- which(seg.chrom %in% use.chrom)
  }

  # Select segments for the desired chromosome number:
  sel.segments <- segments[keepchrom, , drop = FALSE]

  # Check that specified sampleIDs are found in segments:
  if (!is.null(sample)) {
    # Get sampleID for selected samples:
    id <- as.character(unique(seg.sampleid))
    sampleID <- id[sample]
    if (any(is.na(sampleID))) {
      stop(
        "Input in 'sample' is outside the number of samples represented in segments",
        call. = FALSE
      )
    }

    # selct segments for this sampleID(s)
    if (multi) {
      keepsample <- which(colnames(segments) %in% sampleID)
      sel.segments <- sel.segments[, c(1:5, keepsample)]
      colnames(sel.segments) <- c(
        colnames(segments)[1:5],
        colnames(segments)[keepsample]
      )
    } else {
      keepsample <- which(sel.segments[, 1] %in% sampleID)
      sel.segments <- sel.segments[keepsample, , drop = FALSE]
    }
  }

  return(as.data.frame(sel.segments))
}
