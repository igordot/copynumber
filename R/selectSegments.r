####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# FUNCTION THAT RETURNS A SELECTION OF MULTIPCF SEGMENTS BASED ON A DESIRED CHARACHTERISTIC

## Input:
### segments: data frame with segmentation results from multipcf
### what: selection criteria; one of "variance","length" and "aberration"
### thres: optional threshold for selcting, default is NULL
### nseg: desired number of segments, only used if thres=NULL
### p: minimum proportion when what="aberration"
### large: logical - select segments where variance, length or mean value is large (TRUE) or small (FALSE) for what equal to variance, length and aberration respectively.

## Required by:
###  none

## Requires:
### is.multiseg
### pullOutContent



#' Select multipcf segments
#' 
#' Selects multipcf segments based on a desired characteristic.
#' 
#' The input in \code{what} determines how the segments are selected. Three
#' options are available:
#' 
#' If \code{what="variance"} the variance of the segment values across all
#' samples is calculated for each segment. If \code{thres} is specified, the
#' subset of segments for which the variance is above (if \code{large=TRUE}) or
#' below (if \code{large=FALSE}) the threshold is returned. If \code{thres} is
#' not given by the user, a given number of segments determined by the input in
#' \code{nseg} is selected; if \code{large=TRUE} this will be the \code{nseg}
#' segments with the highest variance, whereas if \code{large=FALSE} the subset
#' will consist of the \code{nseg} segments with the lowest variance.
#' 
#' If \code{what="length"} selection is based on the genomic length of the
#' segments (end position minus start position). If \code{thres} is specified,
#' the subset of segments for which the length is above (if \code{large=TRUE})
#' or below (if \code{large=FALSE}) this threshold is returned. If \code{thres}
#' is left unspecified, a given number of segments determined by the input in
#' \code{nseg} is selected; if \code{large=TRUE} this will be the \code{nseg}
#' longest segments, whereas if \code{large=FALSE} it will be the \code{nseg}
#' shortest segments.
#' 
#' If \code{what="aberration"} the aberration frequency is used to select the
#' subset of segments. If \code{thres} is specified, the proportion of samples
#' for which the segment value is above (if \code{large=TRUE}) or below (if
#' \code{large=FALSE}) the threshold is calculated for each segment. The subset
#' of segments where this frequency is above or equal to the proportion set by
#' the parameter \code{p} is returned. If \code{thres} is not specified, the
#' \code{nseg} segments with the highest (1-p)-quantile (if \code{large=TRUE})
#' or the lowest p-quantile (if \code{large=FALSE}) is returned.
#' 
#' @param segments a data frame containing segments found by
#' \code{\link{multipcf}}.
#' @param what the desired characteristic to base selection on. Must be one of
#' "variance" (default),"length" and "aberration". See details below.
#' @param thres an optional numeric threshold to be applied in the selection.
#' @param nseg the desired number of segments to be selected, default is 10.
#' Only used if \code{thres=NULL}.
#' @param large logical value indicating whether segments with large (TRUE) or
#' small (FALSE) variance, length or mean value should be selected when
#' \code{what} is "variance", "length" or "aberration", respectively.
#' @param p a number between 0 and 1 giving the minimum proportion of samples
#' for which an aberration must be detected, default is 0.1. Only applicable if
#' \code{what="aberration"}.
#' @return A list containing: \item{sel.seg}{data frame containing the selected
#' segments.} In addition, depending on the value of \code{what}:
#' \item{seg.var}{a vector giving the variance for each segment. Only returned
#' if \code{what = "variance"}.} \item{seg.length}{a vector giving the length
#' of each segment. Only returned if \code{what = "length"}.}
#' \item{seg.ab.prop}{a vector giving the aberration proportion for each
#' segment. Only returned if \code{what = "aberration"} and \code{thres} is
#' specified.} \item{seg.quantile}{a vector giving the (1-p)- or p-quantile for
#' each segment. Only returned if \code{what = "aberration"} and
#' \code{thres=NULL}.}
#' @author Gro Nilsen
#' @seealso \code{\link{multipcf}}
#' @examples
#' 
#' #Lymphoma data
#' data(lymphoma)
#' 
#' #Run multipcf
#' segments <- multipcf(lymphoma,gamma=12)
#' 
#' #Select the 10 segments with the highest variance:
#' sel.seg1 <- selectSegments(segments)
#' 
#' #Select the segments where the variance is below 0.001
#' sel.seg2 <- selectSegments(segments,thres=0.001,large=FALSE)
#' 
#' #Select the 5 longest segments:
#' sel.seg3 <- selectSegments(segments,what="length",nseg=5)
#' 
#' #Select the segments where 20 % of the samples have segment value of 0.2 or more:
#' sel.seg4 <- selectSegments(segments,what="aberration",thres=0.2,p=0.2)
#' 
#' #Select the 20 segments with the largest median:
#' sel.seg5 <- selectSegments(segments,what="aberration",nseg=20,p=0.5)
#' 
#' @export
selectSegments <- function(
  segments,
  what = "variance",
  thres = NULL,
  nseg = 10,
  large = TRUE,
  p = 0.1
) {
  # Check input
  if (!what %in% c("variance", "length", "aberration")) {
    stop("'what' must be one of variance, length and aberration")
  }
  # Make sure segments is a data frame
  segments <- pullOutContent(res = segments, what = "segments")

  if (!is.multiseg(segments)) {
    stop("'segments' must be on the format resulting from running multipcf!")
  }

  if (is.null(thres) && nseg > nrow(segments)) {
    nseg <- nrow(segments)
    warning(
      "'nseg' is larger than number of rows in 'segments'. Returning all segments.",
      call. = FALSE
    )
  }

  sel.res <- switch(
    what,
    variance = subset.var(segments, nseg, thres, large),
    length = subset.length(segments, nseg, thres, large),
    aberration = subset.abe(segments, nseg, thres, p, large)
  )

  # Sort sel.seg according to chromosome numbers:
  sel.res$sel.seg <- sel.res$sel.seg[order(sel.res$sel.seg[, 1]), ]

  return(sel.res)
}

subset.var <- function(segments, nseg, thres, large) {
  # calculate variance across samples for each segment:
  seg.var <- apply(segments[, -c(1:5)], 1, var)

  if (!is.null(thres)) {
    # Find the segments with variance above thres
    if (large) {
      sel.seg <- segments[seg.var > thres, ]
      if (nrow(sel.seg) == 0) {
        warning(
          paste(
            "none of the segments have variance above ",
            thres,
            ". Returning empty data frame.",
            sep = ""
          ),
          call. = FALSE
        )
      }
    } else {
      sel.seg <- segments[seg.var < thres, ]
      if (nrow(sel.seg) == 0) {
        warning(
          paste(
            "none of the segments have variance below ",
            thres,
            ". Returning empty data frame.",
            sep = ""
          ),
          call. = FALSE
        )
      }
    }
  } else {
    # Find the nseg segments with the highest variance
    if (large) {
      sel.seg <- segments[order(seg.var, decreasing = TRUE)[seq_len(nseg)], ]
    } else {
      sel.seg <- segments[order(seg.var, decreasing = FALSE)[seq_len(nseg)], ]
    }
  }

  return(list(sel.seg = sel.seg, seg.var = seg.var))
}

subset.length <- function(segments, nseg, thres, large) {
  # Find length of each segment:
  L <- segments[, 4] - segments[, 3] + 1
  if (!is.null(thres)) {
    if (large) {
      # Pick out long segments:
      sel.seg <- segments[L > thres, ]
      if (nrow(sel.seg) == 0) {
        warning(
          paste(
            "none of the segments are longer than ",
            thres,
            ". Returning empty data frame.",
            sep = ""
          ),
          call. = FALSE
        )
      }
    } else {
      # Pick out short segments:
      sel.seg <- segments[L < thres, ]
      if (nrow(sel.seg) == 0) {
        warning(
          paste(
            "none of the segments are shorter than ",
            thres,
            ". Returning empty data frame.",
            sep = ""
          ),
          call. = FALSE
        )
      }
    }
  } else {
    if (large) {
      sel.seg <- segments[order(L, decreasing = TRUE)[seq_len(nseg)], ]
    } else {
      sel.seg <- segments[order(L, decreasing = FALSE)[seq_len(nseg)], ]
    }
  }

  return(list(sel.seg = sel.seg, seg.length = L))
}

subset.abe <- function(segments, nseg, thres, p, large) {
  if (!is.null(thres)) {
    if (large) {
      prop.ab <- rowMeans(segments[, -c(1:5)] > thres)
    } else {
      prop.ab <- rowMeans(segments[, -c(1:5)] < thres)
    }

    sel.seg <- segments[prop.ab >= p, ]
    if (nrow(sel.seg) == 0) {
      if (large) {
        warning(
          paste(
            "none of the segments have mean value above ",
            thres,
            "for minimum ",
            p * 100,
            "% of the samples. Returning empty data frame.",
            sep = ""
          ),
          call. = FALSE
        )
      } else {
        warning(
          paste(
            "none of the segments have mean value below ",
            thres,
            "for minimum ",
            p * 100,
            "% of the samples. Returning empty data frame.",
            sep = ""
          ),
          call. = FALSE
        )
      }
    }
    return(list(sel.seg = sel.seg, seg.ab.prop = prop.ab))
  } else {
    if (large) {
      q <- apply(segments[, -c(1:5)], 1, quantile, probs = 1 - p, type = 1)
      q.ord <- order(q, decreasing = TRUE)
    } else {
      q <- apply(segments[, -c(1:5)], 1, quantile, probs = p, type = 1)
      q.ord <- order(q, decreasing = FALSE)
    }
    sel.seg <- segments[q.ord[seq_len(nseg)], ]

    return(list(sel.seg = sel.seg, seg.quantile = q))
  }
}
