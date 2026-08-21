####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# Function that plots the frequency of deletions and amplifications given a threshold on the genome represented by a circle. Can also add arc to represent connections between different
# part of the genome

## Input:
### segments: segmentations results from pcf or multipcf
### thres.gain,thres.loss: a vector with threshold(s) to be applied for abberration calling
### pos.unit: unit used for positions (bp,kbp,mbp)
### freq.colors: colors used to plot gain and loss frequencies
### alpha: a scalar between 0 and 1 that sets the amount of scaling of frequencies in order to fit in circle
### arcs: a data frame with 5 columns. The first four columns gives chromosome numbers and local positions for start point and end point of arc, respectively, while the last column should contain a vector of numbers 1,2,... indication that the arcs belong to different classes. Each class of arcs will then be plotted in a different color
### arc.colors: colors used to plot the lines representing the different classes in arcs. Must be longer or equal to the number of classes found in arcs.
### d: a scalar > 0 representing the distance from the genome circle to the starting points of the arcs. Set d=0 to make arcs start at the genome circle
### assembly: which assembly to use, must be one of hg19, hg18, hg17 or hg16

## Required by: none

### Requires:
## getGlobPos
## numericChrom
## circ
## c.lines
## pullOutContent
## getFreqData



#' Plot a circular genome with aberration frequencies and connections between
#' genomic loci added.
#' 
#' A circular genome is plotted and the percentage of samples that have a gain
#' or a loss at a genomic position is added in the middle of the circle.
#' Gains/losses correspond to copy number values that are above/below a
#' pre-defined threshold. In addition arcs representing some connection between
#' genomic loci may be added.
#' 
#' To zoom in on the observed aberration frequencies one may increase
#' \code{alpha}. However, the user should be aware that this implies that the
#' distance between the genome circle and the frequency zero-line does not
#' reflect an aberration frequency of 100 \%. Since the distance between the
#' two circles is always 1/7, the maximum plotted percentage will be
#' 100/(alpha*7) and any percentages that are higher than this will be
#' truncated to this value.
#' 
#' @param segments a data frame containing the segmentation results found by
#' either \code{\link{pcf}} or \code{\link{multipcf}}.
#' @param thres.gain a scalar giving the threshold value to be applied for
#' calling gains.
#' @param thres.loss a scalar giving the threshold value to be applied for
#' calling losses. Default is to use the negative value of \code{thres.gain}.
#' @param pos.unit the unit used to represent the probe positions. Allowed
#' options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base
#' pairs). By default assumed to be "bp".
#' @param freq.colors a vector giving the colors to be used for the
#' amplification and deletion frequencies, respectively. Default is
#' c("red","limegreen").
#' @param alpha a scalar in the range 0 to 1 determining the amount of scaling
#' of the aberration frequencies. For the default value of 1/7 the distance
#' between the genome circle and the zero-line of the frequency-circle
#' corresponds to an aberration percentage of 100 \%. See details.
#' @param arcs an optional matrix or data frame with 5 columns specifying
#' connections between genomic loci. The first two columns must give the
#' chromosome numbers and local positions for the start points of the arcs,
#' while the two next columns give the chromosome numbers and local positions
#' for the end point of arcs. The last column should contain a vector of
#' numbers 1,2,... indicating that the arcs belong to different classes. Each
#' class of arcs will then be plotted in a different color.
#' @param arc.colors a vector giving the colors to be used for plotting the
#' different classes of arcs. Cannot be shorter than the number of unique
#' classes in \code{arcs}. The first color will represent the first class in
#' \code{arcs}, the second color the second class and so on.
#' @param d a scalar > 0 representing the distance from the genome circle to
#' the starting points of the arcs. Set d=0 to make arcs start at the genome
#' circle.
#' @param assembly a string specifying which genome assembly version should be
#' applied to determine chromosome ideograms. Allowed options are "hg19",
#' "hg18", "hg17" and "hg16" (corresponding to the four latest human genome
#' annotations in the UCSC genome browser).
#' @author Gro Nilsen
#' @examples
#' 
#' #load lymphoma data
#' data(lymphoma)
#' #Run pcf
#' pcf.res <- pcf(data=lymphoma,gamma=12)
#' 
#' plotCircle(segments=pcf.res,thres.gain=0.1)
#' 
#' #Use alpha to view the frequencies in more detail:
#' plotCircle(segments=pcf.res,thres.gain=0.1,alpha=1/5)
#' 
#' #An example of how to specify arcs
#' #Using multipcf, we compute the correlation between all segments and then
#' #retrieve those that have absolute inter-chromosomal correlation > 0.7
#' multiseg <- multipcf(lymphoma)
#' nseg = nrow(multiseg)
#' cormat = cor(t(multiseg[,-c(1:5)]))
#' chr.from <- c()
#' pos.from <- c()
#' chr.to <- c()
#' pos.to <- c()
#' cl <- c()
#' 
#' thresh = 0.7
#' for (i in 1:(nseg-1)) {
#'   for (j in (i+1):nseg) {
#'     #Check if segment-correlation is larger than threshold and that the two 
#'     #segments are located on different chromosomes
#'     if (abs(cormat[i,j]) > thresh && multiseg$chrom[i] != multiseg$chrom[j]) {
#'       chr.from = c(chr.from,multiseg$chrom[i])
#'       chr.to = c(chr.to,multiseg$chrom[j])
#'       pos.from = c(pos.from,(multiseg$start.pos[i] + multiseg$end.pos[i])/2)
#'       pos.to = c(pos.to,(multiseg$start.pos[j] + multiseg$end.pos[j])/2)
#'       if(cormat[i,j] > thresh){
#'         cl <- c(cl,1)           #class 1 for those with positive correlation
#'       }else{
#'         cl <- c(cl,2)           #class 2 for those with negative correlation
#'       }    
#'     }
#'   }
#' }
#'   
#' arcs <- cbind(chr.from,pos.from,chr.to,pos.to,cl)  
#' 
#' #Plot arcs between segment with high correlations; positive correlation in
#' #orange, negative correlation in blue:
#' plotCircle(segments=pcf.res,thres.gain=0.15,arcs=arcs,d=0)
#' 
#' 
#'   
#' 
#' 
#' @export
plotCircle <- function(
  segments,
  thres.gain,
  thres.loss = -thres.gain,
  pos.unit = "bp",
  freq.colors = c("red", "limegreen"),
  alpha = 1 / 7,
  arcs = NULL,
  arc.colors = c("goldenrod1", "dodgerblue"),
  d = 0.3,
  assembly = "hg19"
) {
  delta <- 20000000 # defines the amount of space to be plotted between each chromosome

  # Empty plot
  par(mar = c(0, 0, 0, 0))
  plot(
    0,
    0,
    xlim = c(-1.2, 1.2),
    ylim = c(-1.2, 1.2),
    axes = F,
    xlab = "",
    ylab = "",
    asp = 1,
    type = "n"
  )

  outer.circ <- 0 # location of outer circle
  inner.circ <- 0.1 # location of inner circle

  # Plot cytoband

  cytoband <- get(assembly)
  cyto.chr <- substring(cytoband[, 1], 4)
  cyto.chr[grep("X", cyto.chr)] <- 23
  cyto.chr[grep("Y", cyto.chr)] <- 24
  cyto.chr <- as.numeric(cyto.chr)
  o <- order(cyto.chr)
  cyto.chr <- cyto.chr[o]
  cytoband <- cytoband[o, ]

  chr.end <- diff(cyto.chr) > 0
  cyto.start <- getGlobPos(
    cyto.chr,
    cytoband[, 2],
    pos.unit,
    cytoband,
    delta = delta
  )
  cyto.end <- getGlobPos(
    cyto.chr,
    cytoband[, 3],
    pos.unit,
    cytoband,
    delta = delta
  )
  cyto.end[chr.end[1]:length(cyto.end)] <- cyto.end[
    chr.end[1]:length(cyto.end)
  ] +
    delta # need to add delta between first and second chromosome

  cyto.stain <- cytoband[, 5]
  xmin <- cyto.start[1]
  xmax <- cyto.end[length(cyto.end)]
  ngrid <- 5000
  gridpts <- seq(xmin, xmax, length = ngrid)
  x0 <- gridpts
  x1 <- x0
  y0 <- rep(0, ngrid)
  y1 <- rep(inner.circ, ngrid)
  tmp0 <- circ(x0, y0, xmax)
  tmp1 <- circ(x1, y1, xmax)

  # Get colors for the cytoband
  c0 <- rep(0, ngrid)
  chr.stop.pos <- c(cyto.end[chr.end], cyto.end[length(cyto.end)])
  delta.gridpts <- c()
  for (i in 1:ngrid) {
    grid.chr <- which(chr.stop.pos >= gridpts[i])[1]
    # Check if we are within chromosome or in the delta-region added at the end of each chrom
    if (gridpts[i] <= (chr.stop.pos[grid.chr] - delta)) {
      band <- which(cyto.end >= gridpts[i])[1]
      c0[i] <- switch(
        as.character(cyto.stain[band]),
        "gneg" = "white",
        "gpos100" = "black",
        "gpos75" = "gray25",
        "gpos50" = "gray50",
        "gpos25" = "gray75",
        "stalk" = "gray90",
        "gvar" = "grey",
        "acen" = "yellow"
      )
    } else {
      delta.gridpts <- c(delta.gridpts, i)
      # Between chromosomes the color should be white
      c0[i] <- "white"
    }
  }
  # Plot cytobands:
  segments(
    tmp0$x[-delta.gridpts],
    tmp0$y[-delta.gridpts],
    tmp1$x[-delta.gridpts],
    tmp1$y[-delta.gridpts],
    col = c0[-delta.gridpts],
    lwd = 2
  )

  # Plot outer circle
  x <- seq(0, 1, len = 2000)
  y <- rep(outer.circ, 2000)
  c.lines(x, y, xmax = 1)

  # Plot inner circle
  x <- seq(0, 1, len = 2000)
  y <- rep(inner.circ, 2000)
  c.lines(x, y, xmax = 1)

  # Plot the white area between each chromosome:
  segments(
    tmp0$x[delta.gridpts],
    tmp0$y[delta.gridpts],
    tmp1$x[delta.gridpts],
    tmp1$y[delta.gridpts],
    col = c0[delta.gridpts],
    lwd = 2
  )

  # Plot chromosome borders
  borders <- c(chr.stop.pos, chr.stop.pos - delta)
  x0 <- borders
  x1 <- x0
  y0 <- rep(outer.circ - 0.04, length(borders)) # How much lines separating chromosomes will point out of the outer circle
  y1 <- rep(inner.circ, length(borders))
  tmp0 <- circ(x0, y0, xmax)
  tmp1 <- circ(x1, y1, xmax)
  segments(tmp0$x, tmp0$y, tmp1$x, tmp1$y, lwd = 3)

  # Plot chromosome numbers
  chrmiddle <- (chr.stop.pos[1:23] + (chr.stop.pos - delta)[2:24]) / 2
  chrmiddle <- c((chr.stop.pos[1] - delta) / 2, chrmiddle) # add middle of chrom 1
  x0 <- chrmiddle
  y0 <- rep(outer.circ - 0.1, length(chrmiddle))
  tmp <- circ(x0, y0, xmax)
  text(tmp$x, tmp$y, c(1:22, "X", "Y"))

  # Plot frequency of aberration
  # Check data input:
  # Make sure data input is a data frame (could be a list if it comes from segmentation algorithm or winsorize). Could also be a segments data frame or a data frame with original data.
  # data <- pullOutContent(res=data,what="estimates")
  data <- segments
  if (is.data.frame(data)) {
    # Check if segments or data:
    data <- getFreqData(data)
  } else {
    data <- pullOutContent(res = data, what = "estimates")
  }
  stopifnot(ncol(data) >= 3) # something is missing in input data

  # Make sure that chromosomes in data are numeric:
  chr <- numericChrom(data[, 1])
  pos <- data[, 2]
  yhat <- data[, -c(1:2)]
  freq.circ <- 1 / 7 + inner.circ # placement of frequencies-circle
  m <- freq.circ + (freq.circ - inner.circ) # max space for amp frequencies
  AmpFreq <- apply(yhat > thres.gain, 1, mean)
  glob.pos <- getGlobPos(chr, pos, pos.unit, cytoband, delta = delta)
  x0 <- glob.pos
  y0 <- rep(freq.circ, length(x0))
  x1 <- x0
  y1 <- freq.circ + 0.01 + AmpFreq * alpha # Add 0.01 to avoid overlapping colors
  y1[y1 > m] <- m # make sure that frequencies don't exceed the maximum space; truncate large values to this max
  tmp0 <- circ(x0, y0, xmax)
  tmp1 <- circ(x1, y1, xmax)
  segments(tmp0$x, tmp0$y, tmp1$x, tmp1$y, col = freq.colors[1], lwd = 2)

  m <- freq.circ - (freq.circ - inner.circ) # max space for del frequencies
  DelFreq <- apply(yhat < thres.loss, 1, mean)
  x0 <- glob.pos
  y0 <- rep(freq.circ, length(x0))
  x1 <- x0
  y1 <- freq.circ - 0.01 - DelFreq * alpha # Subtract 0.01 to avoid overlapping colors
  y1[y1 < m] <- m # make sure that frequencies don't exceed the maximum space; truncate large values to this max
  tmp0 <- circ(x0, y0, xmax)
  tmp1 <- circ(x1, y1, xmax)
  segments(tmp0$x, tmp0$y, tmp1$x, tmp1$y, col = freq.colors[2], lwd = 2)

  # Plot x-axis (x-circle)
  x <- seq(0, 1, len = 2000)
  y <- rep(freq.circ, 2000)
  c.lines(x, y, xmax = 1)

  # Plot arcs if specified
  if (!is.null(arcs)) {
    cl <- arcs[, 5]
    u.cl <- sort(unique(cl))
    if (length(arc.colors) < length(u.cl)) {
      arc.colors <- rep(arc.colors, length(u.cl))
      warning(
        "Number of colors in 'arc.colors' is fewer than number of unique classes in 'arcs'. Colors are reused.",
        call. = FALSE
      )
    }
    chr0 <- arcs[, 1]
    chr1 <- arcs[, 3]
    pos0 <- arcs[, 2]
    pos1 <- arcs[, 4]
    x0 <- getGlobPos(chr0, pos0, pos.unit, cytoband, delta = delta)
    x1 <- getGlobPos(chr1, pos1, pos.unit, cytoband, delta = delta)
    m <- d + inner.circ # determine where arcs should start
    for (i in 1:length(x0)) {
      tmp <- circ(c(x0[i], 0, x1[i]), c(0.1, 1.0, 0.1), xmax)
      tmp <- circ(c(x0[i], 0, x1[i]), c(m, 1.0, m), xmax)
      a0 <- xspline(tmp$x, tmp$y, shape = c(0, 1, 0), draw = F)
      lines(a0$x, a0$y, col = arc.colors[which(u.cl == cl[i])])
    }
  } # endif
}

# Circos transformation

circ <- function(x, y, xmax) {
  x <- x * 2 * pi / xmax
  xnew <- (1 - y) * cos(-x + pi / 2)
  ynew <- (1 - y) * sin(-x + pi / 2)
  list(x = xnew, y = ynew)
}


c.lines <- function(x, y, xmax, ...) {
  tmp <- circ(x, y, xmax)
  lines(tmp$x, tmp$y, ...)
}
