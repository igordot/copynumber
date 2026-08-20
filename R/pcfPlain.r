####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

## Required by:

## Requires:
### findNN
### getMad
### pcf(not the main function, but the rest of the help functions found in the same document)
### handleMissing
### pullOutContent

## Main function for pcf-analysis to be called by the user



#' Plain single-sample copy number segmentation.
#' 
#' A basic single-sample pcf segmentation which does not take chromosome
#' borders into account
#' 
#' A piecewise constant segmentation curve is fitted to the copy number
#' observations as described in the PCF algorithm in Nilsen and Liestoel et al.
#' (2012). Unlike the regular \code{\link{pcf}} function, \code{pcfPlain} does
#' not make independent segmentations for each chromosome arm (i.e. breakpoints
#' are not automatically inserted at the beginning and end of chromosome arms).
#' The segmentation can thus be performed independently of assembly.
#' 
#' @param pos.data a data frame where the rows represent the probes, column 1
#' holds probe positions, and subsequent column(s) give the numeric copy number
#' measurements for one or more samples. The header of copy number columns
#' should give sample IDs.
#' @param kmin minimum number of probes in each segment, default is 5.
#' @param gamma penalty for each discontinuity in the curve, default is 40.
#' @param normalize logical value indicating whether the copy number
#' measurements should be scaled by the sample residual standard error. Default
#' is TRUE.
#' @param fast a logical value indicating whether a fast (not guaranteed to be
#' exact) version should be run if the number of probes are > 400.
#' @param digits the number of decimals to be applied when reporting results.
#' Default is 4.
#' @param return.est logical value indicating whether a data frame holding copy
#' number estimates (pcf values) should be returned along with the segments.
#' Default is FALSE, which means that only segments are returned.
#' @param verbose logical value indicating whether or not to print a progress
#' message each time pcf analysis is finished for a sample.
#' @return If \code{return.est = TRUE} a list with the following components:
#' \item{estimates}{a data frame where the first column gives the probe
#' positions, while subsequent column(s) give the copy number estimates for
#' each sample. The estimate for a given probe equals the mean of the segment
#' where the probe is located.} \item{segments}{a data frame describing each
#' segment found in the data. Each row represents a segment, while columns give
#' the sampleID, start position, end position, number of probes in the segment
#' and mean value, respectively.}
#' 
#' If \code{return.est = FALSE}, only the data frame containing the segments is
#' returned.
#' @note If probe positions are not available, the first column in \code{data}
#' may, e.g., contain the values \code{1:nrow(data)}.
#' @author Gro Nilsen, Knut Liestoel, Ole Christian Lingjaerde.
#' @seealso \code{\link{pcf}}
#' @references Nilsen and Liestoel et al., "Copynumber: Efficient algorithms
#' for single- and multi-track copy number segmentation", BMC Genomics 13:591
#' (2012), doi:10.1186/1471-2164-13-59
#' @examples
#' 
#' #Load the lymphoma data set:
#' data(lymphoma)
#' 
#' #Take out a smaller subset of 3 samples (using subsetData):
#' sub.lymphoma <- subsetData(lymphoma,sample=1:3)
#' 
#' #Run pcfPlain (remove first column of chromosome numbers):
#' plain.segments <- pcfPlain(pos.data=sub.lymphoma[,-1],gamma=12)
#' 
#' 
#' 
#' @export
pcfPlain <- function(
  pos.data,
  kmin = 5,
  gamma = 40,
  normalize = TRUE,
  fast = TRUE,
  digits = 4,
  return.est = FALSE,
  verbose = TRUE
) {
  # Input could come from winsorize and thus be a list; check and possibly retrieve data frame wins.data
  pos.data <- pullOutContent(pos.data, what = "wins.data")

  # Make sure all data columns are numeric:
  if (any(!sapply(pos.data, is.numeric))) {
    stop("All input in pos.data must be numeric", call. = FALSE)
  }

  # Check data input:
  stopifnot(ncol(pos.data) >= 2) # something is missing in input data

  # Extract information from pos.data:
  position <- pos.data[, 1]
  cn.data <- as.matrix(pos.data[, -1])
  nSample <- ncol(pos.data) - 1
  sampleid <- colnames(pos.data)[-1]
  nProbe <- length(position)

  # save user's gamma
  gamma0 <- gamma

  sd <- rep(1, nSample) # sd is only used if normalize=TRUE, and then these values are replaced by MAD-sd
  # If number of probes in entire data set is less than 100K, the MAD sd-estimate is calculated using all obs for every sample
  # Only required if normalize=T
  if (nProbe < 100000 && normalize) {
    # calculate MAD-sd for each sample:
    for (j in 1:nSample) {
      sample.data <- pos.data[, j + 1]
      sd[j] <- getMad(sample.data[!is.na(sample.data)], k = 25) # Take out missing values before calculating mad
    }
  } # endif

  # Initialize
  pcf.names <- c("pos", sampleid)
  seg.names <- c("sampleID", "start.pos", "end.pos", "n.probes", "mean")

  segments <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(segments) <- seg.names
  if (return.est) {
    pcf.est <- matrix(nrow = 0, ncol = nSample)
  }

  # Run PCF separately for each sample:
  for (i in 1:nSample) {
    if (return.est) {
      # Initialize:
      yhat <- rep(NA, length(nProbe))
    }

    sample.data <- cn.data[, i]

    # Remove probes with missing obs; Only run pcf on non-missing values
    obs <- !is.na(sample.data)
    obs.data <- sample.data[obs]

    if (length(obs.data) > 0) {
      ## Make sure there are observations for this sample! If not, estimates are left NA as well

      # If number of probes in entire data set is >= 100K, the MAD sd-estimate is calculated using obs in this arm for this sample.
      # Only required if normalize=T
      if (nProbe >= 100000 && normalize) {
        sd[i] <- getMad(obs.data, k = 25)
      }

      # Scale gamma by variance if normalize is TRUE
      use.gamma <- gamma0
      if (normalize) {
        use.gamma <- gamma0 * (sd[i])^2
      }

      # Must check that sd!=0 and sd!!=NA -> use.gamma=0/NA. If not, simply calculate mean of observations
      if (use.gamma == 0 || is.na(use.gamma)) {
        if (return.est) {
          res <- list(
            Lengde = length(obs.data),
            sta = 1,
            mean = mean(obs.data),
            nIntervals = 1,
            yhat = rep(mean(obs.data))
          )
        } else {
          res <- list(
            Lengde = length(obs.data),
            sta = 1,
            mean = mean(obs.data),
            nIntervals = 1
          )
        }
      } else {
        # Compute piecewise constant fit
        # run fast approximate PCF if fast=TRUE and number of probes>400, or exact PCF otherwise
        if (!fast || length(obs.data) < 400) {
          # Exact PCF:
          res <- exactPcf(
            y = obs.data,
            kmin = kmin,
            gamma = use.gamma,
            yest = return.est
          )
        } else {
          # Run fast PCF:
          res <- selectFastPcf(
            x = obs.data,
            kmin = kmin,
            gamma = use.gamma,
            yest = return.est
          )
        } # endif
      } # endif

      # Retrieve segment info from results:
      seg.start <- res$sta
      seg.stop <- c(seg.start[-1] - 1, length(obs.data))
      seg.npos <- res$Lengde
      seg.mean <- res$mean
      nSeg <- res$nIntervals
      if (return.est) {
        yhat[obs] <- res$yhat
      }
      # Find genomic positions for start and stop of each segment:
      pos.start <- position[seg.start]
      pos.stop <- position[seg.stop]

      # Handle missing values:
      if (any(!obs)) {
        # first find nearest non-missing neighbour for missing probes:
        nn <- findNN(pos = position, obs = obs)

        # Include probes with missing values in segments where their nearest neighbour probes are located
        new.res <- handleMissing(
          nn = nn,
          pos = position,
          obs = obs,
          pos.start = pos.start,
          pos.stop = pos.stop,
          seg.npos = seg.npos
        )
        pos.start <- new.res$pos.start
        pos.stop <- new.res$pos.stop
        seg.npos <- new.res$seg.npos

        if (return.est) {
          yhat[!obs] <- yhat[nn]
        }
      }
    } else {
      warning(
        paste(
          "pcf is not run for sample ",
          i,
          " because all observations are missing. NA is returned.",
          sep = ""
        ),
        immediate. = TRUE,
        call. = FALSE
      )
      seg.start <- 1
      seg.stop <- nProbe
      pos.start <- position[seg.start]
      pos.stop <- position[seg.stop]
      nSeg <- 1
      seg.mean <- NA
      seg.npos <- nProbe
    }

    # Round:
    seg.mean <- round(seg.mean, digits = digits)

    # Add results for this sample to results for other samples in data frame:
    seg <- data.frame(
      rep(sampleid[i], nSeg),
      pos.start,
      pos.stop,
      seg.npos,
      seg.mean,
      stringsAsFactors = FALSE
    )
    colnames(seg) <- seg.names
    segments <- rbind(segments, seg)

    if (return.est) {
      # Rounding:
      yhat <- round(yhat, digits = digits)
      pcf.est <- cbind(pcf.est, yhat)
    }
  } # endfor

  if (verbose) {
    cat(paste("pcf finished for sample ", i, sep = ""), "\n")
  }

  if (return.est) {
    pcf.est <- data.frame(position, pcf.est, stringsAsFactors = FALSE)
    colnames(pcf.est) <- pcf.names
    return(list(estimates = pcf.est, segments = segments))
  } else {
    return(segments)
  }
} # endfunction
