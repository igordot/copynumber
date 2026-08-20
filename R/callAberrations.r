####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# Function that calls segments as gain, normal or loss

## Input:
### segments: segmentation results from pcf or multipcf
### thres.gain,thres.loss: threshold(s) to be applied for abberration calling

## Required by: none

## Requires:
### is.multiseg
### pullOutContent



#' Call aberrations in segmented data
#' 
#' Segments, obtained by \code{pcf} or \code{multipcf}, are classified as
#' "gain", "normal" or "loss" given the specified thresholds.
#' 
#' Each region found in \code{segments} is classified as "gain", "normal" or
#' "loss". Regions with gain or loss will be those segments where the segment
#' value is above or below the value given in \code{thres.gain} or
#' \code{thres.loss}, respectively.
#' 
#' @param segments a data frame containing the segmentation results found by
#' either \code{\link{pcf}} or \code{\link{multipcf}}.
#' @param thres.gain a numeric value giving the threshold to be applied for
#' calling gains.
#' @param thres.loss a numeric value giving the threshold to be applied for
#' calling losses. Default is to use the negative value of \code{thres.gain}.
#' @return A new segment data frame where the segment values have been replaced
#' by the classification "gain", "normal" or "loss".
#' @author Gro Nilsen
#' @examples
#' 
#' #load lymphoma data
#' data(lymphoma)
#' #Run pcf
#' seg <- pcf(data=lymphoma,gamma=12)
#' 
#' #Call gains as segments whose value is > 0.2, and losses as segments whose
#' # value < -0.1
#' ab.seg <- callAberrations(seg,thres.gain=0.2,thres.loss=-0.1)
#' 
#' 
#' @export
callAberrations <- function(segments, thres.gain, thres.loss = -thres.gain) {
  # Make sure segments is a data frame
  segments <- pullOutContent(res = segments, what = "segments")

  if (is.multiseg(segments)) {
    call.seg <- matrix(
      "normal",
      nrow = nrow(segments),
      ncol = ncol(segments) - 5
    )
    colnames(call.seg) <- colnames(segments)[-c(1:5)]
    gain <- segments[, -c(1:5), drop = FALSE] > thres.gain
    call.seg[gain] <- "gain"
    loss <- segments[, -c(1:5), drop = FALSE] < thres.loss
    call.seg[loss] <- "loss"

    call.seg <- as.data.frame(cbind(segments[, c(1:5)], call.seg))
  } else {
    call.seg <- rep("normal", nrow(segments))
    gain <- segments[, 7] > thres.gain
    call.seg[gain] <- "gain"
    loss <- segments[, 7] < thres.loss
    call.seg[loss] <- "loss"

    call.seg <- as.data.frame(cbind(segments[, c(1:6)], call.seg))
    colnames(call.seg)[7] <- "call"
  }

  return(call.seg)
}
