#' 3K aCGH data
#'
#' A subset of the aCGH data set taken from the reference below.
#'
#' @format Data frame containing 3091 probes with log2-ratio copy numbers for
#'   21 samples. The first column contains the chromosome numbers, the second
#'   gives the local probe positions (in base pairs), while the subsequent
#'   columns contain the copy number measurements for the individual samples.
#' @source Eide et al., "Genomic alterations reveal potential for higher
#'   grade transformation in follicular lymphoma and confirm parallel
#'   evolution of tumor cell clones", Blood 116:1489-1497, 2010
#' @usage data(lymphoma)
#' @examples
#' #Get data
#' data(lymphoma)
"lymphoma"

#' Subset of 244K aCGH data
#'
#' A subset of the 244K MicMa data set containing copy number measurements
#' for six samples on chromosome 17.
#'
#' @format Data frame containing 7658 probes with log2-ratio copy numbers
#'   for 6 samples on chromosome 17. The first column contains the
#'   chromosome numbers, the second gives the local probe positions (in
#'   base pairs), while the subsequent columns contain the copy number
#'   measurements for the individual samples.
#' @source Mathiesen et al., "High resolution analysis of copy number
#'   changes in disseminated tumor cells of patients with breast cancer",
#'   Int J Cancer 131(4):E405:E415, 2011
#' @usage data(micma)
#' @examples
#' #Get data
#' data(micma)
"micma"

#' Artificial SNP array data
#'
#' Artificial SNP array data containing a logR track and a BAF track
#'
#' @format Two corresponding data sets containing 10000 probes with logR
#'   and BAF measurements, respectively, for 2 samples. The two first
#'   columns in both data sets contain chromosome numbers and local probe
#'   positions (in base pairs), while the subsequent columns contain
#'   logR-values and BAF-values in the two data sets, respectively.
#' @usage
#' data(logR)
#' data(BAF)
#' @examples
#' #Get data
#' data(logR)
#' data(BAF)
#' @rdname SNPdata
"logR"

#' @format NULL
#' @usage NULL
#' @rdname SNPdata
"BAF"
