#' Get segments on the GRanges format
#' 
#' The segments data frame obtained e.g. by \code{pcf}, \code{multipcf} or
#' \code{aspcf} is converted to the GRanges format.
#' 
#' GRanges, in the GenomicRanges package, is the standard BioConductor
#' containers for range data. For some applications it may therefore be useful
#' to convert segmentation results to this format.
#' 
#' @param segments a data frame containing segmentation results found by e.g.
#' \code{\link{pcf}}, \code{\link{multipcf}} or \code{\link{aspcf}}.
#' @return The segments converted to the GRanges container class.
#' @author Gro Nilsen
#' @examples
#' 
#' #load lymphoma data
#' data(lymphoma)
#' #Run pcf
#' seg <- pcf(data=lymphoma,gamma=12)
#' #Obtain the GRanges format
#' gr <- getGRangesFormat(seg)
#' 
#' @import S4Vectors
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges mcols<-
#' @export
getGRangesFormat <- function(segments) {
  nc <- ncol(segments)
  if (!is.multiseg(segments)) {
    gr <- GRanges(
      seqnames = segments$chrom,
      ranges = IRanges(
        start = segments$start.pos,
        end = segments$end.pos,
        names = segments$sampleID
      )
    )
    mcols(gr) <- DataFrame(segments[, c(3, 6:nc)])
  } else {
    gr <- GRanges(
      seqnames = segments$chrom,
      ranges = IRanges(start = segments$start.pos, end = segments$end.pos)
    )
    mcols(gr) <- DataFrame(segments[, c(2, 5:nc)])
  }

  return(gr)
}
