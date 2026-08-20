####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

### Requires:
### pullOutContent

# Function with interpolates pcf-segments


#' Interpolation of pcf-estimates.
#' 
#' Given a segmentation by \code{pcf}, interpolate pcf-estimates for specific
#' positions.
#' 
#' Pcf-estimates are interpolated for the chromosomes and postions specified in
#' \code{x}.
#' 
#' @param segments a data frame containing the segmentation results from
#' \code{\link{pcf}}.
#' @param x matrix or data.frame where the first column gives chrosomomes and
#' the second gives positions.
#' @return A data frame where the first two columns give the chromsomes and
#' positions specified in the input \code{x} and the remaining columns give the
#' interpolated pcf-estimate for each sample represented in \code{segments}.
#' @note The positions in \code{segments} and \code{x} must be of the same unit
#' (bp, kbp, or mbp).
#' @author Gro Nilsen, Ole Christian Lingjaerde.
#' @seealso \code{\link{pcf}}
#' @examples
#' 
#' #Load the lymphoma data set:
#' data(lymphoma)
#' 
#' #Take out a smaller subset of 3 samples (using subsetData):
#' sub.lymphoma <- subsetData(lymphoma,sample=1:3)
#' 
#' #Run pcf:
#' seg <- pcf(data=sub.lymphoma,gamma=12)
#' 
#' #Make a matrix with two positions and chromosomes for which we want to 
#' #interpolate the pcf-estimate:
#' pos <-  c(2000000,50000000)
#' chr <- c(1,2)
#' x <- cbind(chr,pos)
#' 
#' #Interpolate
#' int.pcf <- interpolate.pcf(seg,x)
#' 
#' 
#' @rawNamespace export(interpolate.pcf)
interpolate.pcf <- function(segments, x) {
  # Make sure segments is a data frame
  segments <- pullOutContent(res = segments, what = "segments")

  usamp <- unique(segments$sampleID)
  nsamp <- length(usamp)
  chrom <- unique(x[, 1])
  z <- data.frame(cbind(x[, c(1:2)], matrix(0, nrow(x), nsamp)))
  # z = data.frame(x[,c(1,2)], matrix(0, nrow(x), nsamp))
  names(z) <- c("chr", "pos", usamp)
  for (i in 1:nsamp) {
    for (j in 1:length(chrom)) {
      fitij <- segments[
        segments$sampleID == usamp[i] & segments$chrom == chrom[j],
      ]
      v <- (c(fitij$start.pos[-1], 10^9) + fitij$end.pos) / 2
      xj <- x[x[, 1] == chrom[j], 2]
      kj <- rep(0, length(xj))
      for (k in rev(1:length(v))) {
        kj[xj <= v[k]] <- k
      }
      z[z$chr == chrom[j], 2 + i] <- fitij$mean[kj]
    }
  }
  return(z)
}
