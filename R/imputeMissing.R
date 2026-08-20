####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# FUNCTION THAT IMPUTES MISSING VALUES IN COPY NUMBER DATA; EITHER A CONSTANT
# VALUE OR THE PCF-VALUE OF NEAREST NON-MISSING NEIGHBOUR

## Input:
### data : a numeric matrix with chromosome numbers in the first column, local probe positions in the second, and copy number data for one or more samples in subsequent columns
### method : the imputation method to be used. Must be one of "constant" and "pcf"
### c : a numerical value to be imputed if method is "constant". Default is 0
### pcf.est : a numeric matrix of the same size as data, with chromosome numbers and positions in the first two columns, and copy number estimates from \code{pcf} in the subsequent columns. Only applicable if \code{method="pcf"}. If unspecified, \code{pcf} is run internally to find estimates
### ... : other optional parameters to be passed to \code{pcf}

## Output:
### data matrix of same size as data with missing values imputed

## Required by:
### none

## Requires:
###  pcf
### pullOutContent



#' Impute missing copy number values
#' 
#' Missing copy number values are imputed by a constant value or pcf-estimates.
#' 
#' The available imputation methods are: \describe{
#' \item{\code{constant}:}{all missing values in \code{data} are replaced by
#' the specified value \code{c}.} \item{\code{pcf}:}{the estimates from
#' pcf-segmentation (see \code{\link{pcf}}) are used to impute missing values.
#' If \code{pcf} has already been run, these estimates may be specified in
#' \code{pcf.est}. If \code{pcf.est} is unspecified, \code{pcf} is run on the
#' input data. In \code{pcf} the analysis is done on the observed values, and
#' estimates for missing observations are set to be the estimate of the
#' nearest observed probe.} }
#' 
#' @param data a data frame with numeric or character chromosome numbers in the
#' first column, numeric local probe positions in the second, and numeric copy
#' number data for one or more samples in subsequent columns.
#' @param method the imputation method to be used. Must be one of "constant"
#' and "pcf".
#' @param c a numerical value to be imputed if method is "constant". Default is
#' 0.
#' @param pcf.est a data frame of same size as \code{data}, with chromosome
#' numbers and positions in the first two columns, and copy number estimates
#' obtained from \code{pcf} in the subsequent columns. Only applicable if
#' \code{method="pcf"}. If unspecified and \code{method="pcf"}, \code{pcf} is
#' run internally to find estimates.
#' @param \dots other relevant parameters to be passed on to \code{pcf}
#' @return A data frame of the same size and format as \code{data} with all
#' missing values imputed.
#' @author Gro Nilsen
#' @seealso \code{\link{pcf}}
#' @examples
#' 
#' #Load lymphoma data
#' data(lymphoma)
#' chrom <- lymphoma[,1]
#' pos <- lymphoma[,2]
#' #pick out data for the first six samples:
#' cn.data <- lymphoma[,3:8]
#' 
#' #Create missing values in cn.data at random positions:
#' n <- nrow(cn.data)*ncol(cn.data)
#' r <- matrix(rbinom(n=n,size=1,prob=0.95),nrow=nrow(cn.data),ncol=ncol(cn.data))
#' cn.data[r==0] <- NA    #matrix with approximately 5% missing values
#' mis.data <- data.frame(chrom,pos,cn.data)
#' 
#' #Impute missing values by constant, c=0:
#' imp.data <- imputeMissing(data=mis.data,method="constant")
#' 
#' #Impute missing values by obtained pcf-values:
#' pcf.est <- pcf(data=mis.data,return.est=TRUE)
#' imp.data <- imputeMissing(data=mis.data,method="pcf",pcf.est=pcf.est)
#' 
#' #Or run pcf within imputeMissing:
#' imp.data <- imputeMissing(data=mis.data,method="pcf")
#' 
#' 
#' @export
imputeMissing <- function(data, method, c = 0, pcf.est = NULL, ...) {
  # Check input
  # First check if data comes from winsorize, and if so make sure it is a data frame
  data <- pullOutContent(data, what = "wins.data")
  stopifnot(ncol(data) >= 3) # something is missing in input data

  # Check method input:
  if (!method %in% c("constant", "pcf")) {
    stop("method must be one of 'constant' and 'pcf'", call. = FALSE)
  }

  imp.data <- switch(
    method,
    constant = imp.constant(data, c),
    pcf = imp.pcf(data, yhat = pcf.est, ...)
  )

  return(imp.data)
} # endimputeMissing


# Two possible imputation methods:
imp.constant <- function(data, c) {
  cn.data <- data[, -c(1:2), drop = FALSE]
  na <- is.na(cn.data)
  imp.data <- cn.data
  imp.data[na] <- c

  imp.data <- cbind(data[, c(1:2), drop = FALSE], imp.data)
  return(imp.data)
}

imp.pcf <- function(data, yhat, ...) {
  na <- is.na(data)
  imp.data <- data

  if (is.null(yhat)) {
    # Run pcf:
    s <- apply(na, 2, any) # find samples with missing values
    if (any(s[1:2] == TRUE)) {
      stop(
        "missing values are not allowed in data columns 1 and 2",
        call. = FALSE
      )
    }
    pcf.res <- pcf(data = data[, c(1:2, which(s))], return.est = TRUE, ...) # Only run PCF on samples with missing values
    # Make sure yhat is the same size as cn.data to maintain location of na
    yhat <- data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
    yhat[, c(1:2, which(s))] <- pcf.res$estimates
  } else {
    # Make sure pcf.est is a data frame
    yhat <- pullOutContent(res = yhat, what = "estimates")
    if (ncol(yhat) != ncol(data) || nrow(yhat) != nrow(data)) {
      stop("pcf.est must be the same size as data", call. = FALSE)
    }
  }
  imp.data[na] <- yhat[na]

  return(imp.data)
}
