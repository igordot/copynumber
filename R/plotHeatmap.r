####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# Function that plots heatmap given limits - by genome og chromosomes

## Input:
### segments: segmentation results from pcf or multipcf
### upper.lim, lower.lim: a vector with limit(s) to be applied for abberration calling
### pos.unit: unit used for positions (bp,kbp,mbp)
### chrom: a vector with chromosomes to be plotted. If specified the frequencies are plotted with one panel for each chromosome
### layout: number of columns and rows in plot
### ... : other optional plot parameters

## Required by: none

## Requires:
### numericChrom
### checkChrom
### genomeHeat
### chromosomeHeat
### is.multiseg
### checkSegments
### pullOutContent

# Main function for heatmap plotting:


#' Plot copy number heatmap
#' 
#' Heatmap reflecting the magnitude of estimated copy numbers relative to some
#' pre-defined limits. Estimates may be obtained using \code{pcf} or
#' \code{multipcf}, and results may be visualized over the entire genome or by
#' chromosomes.
#' 
#' For each sample, the segments are represented by a rectangle plotted in a
#' color corresponding to the difference between the segment copy number value
#' and the limits. If the value is below \code{lower.lim}, the color of the
#' rectangle will equal the input in \code{colors[1]} (default dodgerblue). If
#' the value is above \code{lower.lim}, but below zero, the color of the
#' rectangle will be a nuance between the input in \code{colors[1]} and
#' \code{colors[2]} (default black). The closer the value is to zero, the
#' closer the nuance will be to \code{colors[2]}. Similary, if the value is
#' above \code{upper.lim}, the color of the rectangle will equal the input in
#' \code{colors[3]} (default red), whereas if the value is below
#' \code{upper.lim}, but above zero, the color will be a nuance between the
#' input in \code{colors[2]} and \code{colors[3]}. Again, the closer the value
#' is to zero, the closer the nuance will be to \code{colors[2]}.
#' 
#' Each row in the heatmap represents a sample, while probe positions are
#' reflected along the x-axis.
#' 
#' @param segments a data frame containing the segmentation results found by
#' either \code{\link{pcf}} or \code{\link{multipcf}}.
#' @param upper.lim a positive numeric vector giving the upper limits(s) to be
#' applied. The colors in the heatmap will reflect the magnitude of the
#' estimated copy numbers relative to this limit, see details.
#' @param lower.lim a negative numeric vector of same length as
#' \code{upper.lim} giving the lower limits(s) to be applied. Default is to use
#' the negative value of \code{upper.lim}.
#' @param pos.unit the unit used to represent the probe positions. Allowed
#' options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base
#' pairs). By default assumed to be "bp".
#' @param chrom a numeric or character vector with chromosome number(s) to
#' indicate which chromosome(s) is (are) to be plotted. If unspecified the
#' whole genome is plotted.
#' @param layout the vector of length two giving the number of rows and columns
#' in the plot window. Default is \code{c(1,1)}.
#' @param \dots other optional graphical parameters. These include the plot arguments
#' \code{xlab}, \code{ylab}, \code{main}, \code{cex.main}, \code{mgp},
#' \code{cex.lab}, \code{cex.axis}, \code{mar} and \code{title} (see
#' \code{\link{par}} on these), as well as \code{plot.size}, \code{plot.unit},
#' \code{plot.ideo}, \code{ideo.frac}, \code{cyto.text}, \code{assembly} and
#' \code{cex.cytotext} (see \code{\link{plotSample}} on these). In addition, a
#' range of graphical arguments specific for this plot function may be
#' specified: \describe{ \item{\code{colors}}{a character vector of length
#' three giving the colors to interpolate in the heatmap, default is
#' c("dodgerblue","black","red").} \item{\code{n.col}}{an integer giving the
#' number of color shades to be applied in the interpolation, default is 50.}
#' \item{\code{sample.labels}}{a logical value indicating whether sample
#' labels are to be plotted along the y-axis. Default is TRUE.}
#' \item{\code{sep.samples}}{a number in the range 0 to 0.4 used to create
#' some space between samples. Default is 0, which implies that there is no
#' space.} \item{\code{sample.line}}{a numeric scalar giving the margin line
#' where the sample labels should be written, starting at 0 counting outwards.
#' Default is 0.2.} \item{\code{sample.cex}}{the size of the sample labels.} }
#' @note This function applies \code{par(fig)}, and is therefore not compatible
#' with other setups for arranging multiple plots in one device such as
#' \code{par(mfrow,mfcol)}.
#' @author Gro Nilsen
#' @examples
#' 
#' #Load lymphoma data
#' data(lymphoma)
#' 
#' #Run pcf to obtain estimated copy number values
#' seg <- pcf(data=lymphoma,gamma=12)
#' 
#' #Heatmap for entire genome, two limit values:
#' plotHeatmap(segments=seg,upper.lim=c(0.1,0.5),layout=c(2,1))
#' 
#' #Heatmap for the first 4 chromosomes:
#' plotHeatmap(segments=seg,upper.lim=0.1,chrom=c(1:4),layout=c(2,2))
#' 
#' 
#' @export
plotHeatmap <- function(
  segments,
  upper.lim,
  lower.lim = -upper.lim,
  pos.unit = "bp",
  chrom = NULL,
  layout = c(1, 1),
  ...
) {
  # If chrom is unspecified, the whole genome is plotted. Otherwise, selected chromosomes are plotted
  type <- ifelse(is.null(chrom), "genome", "bychrom")

  # Check input in segments:
  segments <- pullOutContent(res = segments, what = "segments")

  # Check segments (convert multi segments to unisegments etc)
  segments <- checkSegments(segments, type)

  # Check and if necessary modify chrom to be plotted
  chrom <- checkChrom(data = NULL, segments = segments, chrom)

  # Convert segments back to data frame (is returned as list in checkSegments)
  segments <- as.data.frame(segments)

  # Get sampleIds
  sampleID <- unique(segments[, 1])

  # Check pos.unit input:
  if (!pos.unit %in% c("bp", "kbp", "mbp")) {
    stop("pos.unit must be one of bp, kbp and mbp", call. = FALSE)
  }

  # Make sure upper.lim is positive and lower.lim is negative:
  if (!all(upper.lim > 0) || !all(lower.lim < 0)) {
    stop(
      "upper.lim must be positive and lower.lim must be negative",
      call. = FALSE
    )
  }

  # Making sure number of upperlimits and lowerlimits are the same:
  nT <- min(length(upper.lim), length(lower.lim))
  upper.lim <- upper.lim[1:nT]
  lower.lim <- lower.lim[1:nT]

  # plot heatmap, either over genome or by chromosomes:
  switch(
    type,
    genome = genomeHeat(
      segments,
      upper.lim,
      lower.lim,
      pos.unit,
      sampleID,
      layout,
      ...
    ),
    bychrom = chromosomeHeat(
      segments,
      upper.lim,
      lower.lim,
      pos.unit,
      sampleID,
      chrom,
      layout,
      ...
    )
  )
}


# Function that plots heatmap by genome

## Input:
### segments: segmentation data frame
### upper.lim, lower.lim: a vector with limit(s) to be applied for abberration calling
### pos.unit: unit used for positions (bp,kbp,mbp)
### sampleID: ids for the samples to be plotted
### layout: number of columns and rows in plot
### ... : other optional plot parameters

## Required by: plotHeatmap

## Requires:
### getHeatParameters
### getGlobalXlim
### adjustSegPos
### framedim
### colorSetup
### getCol
### addChromlines

genomeHeat <- function(
  segments,
  upper.lim,
  lower.lim,
  pos.unit,
  sampleID,
  layout,
  ...
) {
  nT <- length(upper.lim)
  nr <- layout[1]
  nc <- layout[2]
  rc <- nr * nc
  nSample <- length(sampleID)

  op <- getHeatParameters(
    type = "genome",
    nc = nc,
    nr = nr,
    nSample = nSample,
    upper.lim = upper.lim,
    lower.lim = lower.lim,
    ...
  )

  # Set global xlimits if not specified by user:
  if (is.null(op$xlim)) {
    op$xlim <- getGlobalXlim(
      op = op,
      pos.unit = pos.unit,
      chrom = unique(segments[, 2])
    )
  }

  # Check if there should be more than one file/window with plot(s), and get file.name accordingly
  if (dev.cur() <= 1) {
    # to make Sweave work
    dev.new(width = op$plot.size[1], height = op$plot.size[2], record = TRUE)
  }

  # Initialize:
  row <- 1
  clm <- 1
  new <- FALSE

  # Division of plotting window:
  frames <- framedim(nr, nc)

  for (t in 1:nT) {
    # Frame dimensions for plot t:
    fig.t <- c(
      frames$left[clm],
      frames$right[clm],
      frames$bot[row],
      frames$top[row]
    )
    par(fig = fig.t, new = new, oma = c(0, 0, 0.5, 0), mar = op$mar)

    # Empty plot with correct dimensions:
    plot(
      1,
      1,
      type = "n",
      ylim = c(0, nSample),
      ylab = op$ylab,
      xlab = op$xlab,
      xlim = op$xlim,
      xaxs = "i",
      yaxt = "n",
      xaxt = "n",
      yaxs = "i",
      cex.lab = 0.9,
      mgp = op$mgp,
      main = ""
    )
    # Let background be black to avoid white parts in arms without probes:
    rect(
      par("usr")[1],
      par("usr")[3],
      par("usr")[2],
      par("usr")[4],
      col = op$colors[2]
    )

    # main title for this plot
    title(main = op$main[t], line = op$main.line, cex.main = op$cex.main)

    # Get colorsetup for these limits:
    cs <- colorSetup(
      upper.lim = upper.lim[t],
      lower.lim = lower.lim[t],
      op = op
    )

    # Plot heatmap for each sample:
    for (i in 1:nSample) {
      sample.segments <- segments[segments[, 1] == sampleID[i], ]

      # Adjust positions to be plotted along xaxis; i.e. get global positions, scale according to plot.unit, and get left and right pos for rectangles
      # to be plotted (either continuous or 1 probe long):
      x <- adjustSegPos(
        chrom = sample.segments[, 2],
        char.arms = sample.segments[, 3],
        start = sample.segments[, 4],
        stop = sample.segments[, 5],
        type = "genome",
        unit = pos.unit,
        op = op
      )
      xleft <- x$use.start
      xright <- x$use.stop

      # Find appropriate colour for each probe:
      heat.col <- sapply(
        sample.segments[, 7],
        getCol,
        colors = cs$colors,
        intervals = cs$intervals
      )

      # Plot rectangles with appropriate color for each probe:
      ytop <- i - op$sep.samples
      ybottom <- i - (1 - op$sep.samples)

      rect(xleft, ybottom, xright, ytop, col = heat.col, border = NA)

      # Add sampleid on yaxis
      if (op$sample.labels) {
        axis(
          side = 2,
          at = (ytop - (ytop - ybottom) / 2),
          labels = sampleID[i],
          line = op$sample.line,
          tcl = 0,
          cex.axis = op$sample.cex,
          las = 1,
          mgp = op$mgp,
          tick = FALSE
        )
      }
    } # endfor

    # Separate chromosomes by vertical lines:
    addChromlines(
      chromosomes = segments[, 2],
      xaxis = "pos",
      unit = pos.unit,
      cex = op$cex.chrom,
      op = op
    )

    # Box:
    abline(v = op$xlim)
    abline(h = c(0, nSample))

    # Get new page, or update column/row:
    if (t %% (nr * nc) == 0) {
      # Start new plot page (prompted by user)
      devAskNewPage(ask = TRUE)

      # Reset columns and row in layout:
      clm <- 1
      row <- 1
      new <- FALSE
    } else {
      # Update column and row index:
      if (clm < nc) {
        clm <- clm + 1
      } else {
        clm <- 1
        row <- row + 1
      } # endif
      new <- TRUE
    } # endif
  } # endfor
}

# Function that plots heatmap by chromosomes (each chrom in separate panel)

## Input:
### segments: segmentation data frame
### upper.lim, lower.lim: a vector with limits(s) to be applied for abberration calling
### pos.unit: unit used for positions (bp,kbp,mbp)
### sampleID: the IDs for samples to be plotted
### chrom: the chromosomes to be plotted
### layout: number of columns and rows in plot
### ... : other optional plot parameters

## Required by: plotHeatmap

## Requires:
### getHeatParameters
### adjustSegPos
### framedim
### plotIdeogram
### chromMax
### colorSetup
### getCol
### getXticks

chromosomeHeat <- function(
  segments,
  upper.lim,
  lower.lim,
  pos.unit,
  sampleID,
  chrom,
  layout,
  ...
) {
  nT <- length(upper.lim)
  nr <- layout[1]
  nc <- layout[2]
  nChrom <- length(chrom)
  nSample <- length(sampleID)

  # Default plot options:
  op <- getHeatParameters(
    type = "bychrom",
    nc = nc,
    nr = nr,
    nSample = nSample,
    upper.lim = upper.lim,
    lower.lim = lower.lim,
    chrom = chrom,
    ...
  )

  # Margins for entire plot in window:
  if (all(op$title == "")) {
    oma <- c(0, 0, 0, 0)
  } else {
    oma <- c(0, 0, 1, 0)
  }
  mar <- c(0.2, 0.2, 0.3, 0.2)

  # Divide the plotting window by the function "framedim":
  frames <- framedim(nr, nc)

  # make separate plots for each value of limits
  for (t in 1:nT) {
    # Start new window/file:
    if (dev.cur() <= 1) {
      # to make Sweave work
      dev.new(width = op$plot.size[1], height = op$plot.size[2], record = TRUE)
    }

    # Initialize row and column index:
    row <- 1
    clm <- 1
    new <- FALSE

    # Get colorsetup for these limits
    cs <- colorSetup(
      upper.lim = upper.lim[t],
      lower.lim = lower.lim[t],
      op = op
    )

    # Separate plots for each chromosome:
    for (c in 1:nChrom) {
      # Frame dimensions for plot c:
      fig.c <- c(
        frames$left[clm],
        frames$right[clm],
        frames$bot[row],
        frames$top[row]
      )
      par(fig = fig.c, new = new, oma = oma, mar = mar)

      # Make list with frame dimensions:
      frame.c <- list(
        left = frames$left[clm],
        right = frames$right[clm],
        bot = frames$bot[row],
        top = frames$top[row]
      )

      # Select relevant chromosome number
      k <- chrom[c]

      # Pick out indeces for this chromosome
      ind.c <- which(segments[, 2] == k)
      chrom.segments <- segments[ind.c, ]

      # Find maximum position for this chromosome (and scale according to plot unit)
      scale.fac <- convert.unit(unit1 = op$plot.unit, unit2 = pos.unit)
      xlim <- c(0, max(chrom.segments[, 5])) * scale.fac

      # Plot ideogram below heatmaps:
      if (op$plot.ideo) {
        # Ideogram frame:
        ideo.frame <- frame.c
        ideo.frame$top <- frame.c$bot +
          (frame.c$top - frame.c$bot) * op$ideo.frac

        par(fig = unlist(ideo.frame), new = new, mar = op$mar.i)
        # Plot ideogram and get maximum probe position in ideogram:
        plotIdeogram(
          chrom = k,
          cyto.text = op$cyto.text,
          cyto.data = op$assembly,
          cex = op$cex.cytotext,
          unit = op$plot.unit
        )

        # Get maximum position for this chromosome:
        xmaxI <- chromMax(
          chrom = k,
          cyto.data = op$assembly,
          pos.unit = op$plot.unit
        )
        xlim <- c(0, xmaxI)

        new <- TRUE
      }

      # Plot-dimensions:
      frame.c$bot <- frame.c$bot + (frame.c$top - frame.c$bot) * op$ideo.frac
      par(fig = unlist(frame.c), new = new, mar = op$mar)

      # Limits:
      if (!is.null(op$xlim)) {
        xlim <- op$xlim
      }

      # empty plot set up:
      plot(
        1,
        1,
        type = "n",
        ylim = c(0, nSample),
        xlim = xlim,
        ylab = op$ylab,
        xlab = op$xlab,
        xaxs = "i",
        yaxt = "n",
        xaxt = "n",
        yaxs = "i",
        mgp = op$mgp,
        main = ""
      )
      # Let background be black to avoid white parts between arms:
      rect(
        par("usr")[1],
        par("usr")[3],
        par("usr")[2],
        par("usr")[4],
        col = op$colors[2]
      )

      # main title for this plot
      title(main = op$main[c], line = op$main.line, cex.main = op$cex.main)

      # Plot heatmap for each sample:
      for (i in 1:nSample) {
        sample.segments <- chrom.segments[chrom.segments[, 1] == sampleID[i], ]
        # Adjust positions to be plotted along xaxis; i.e. scale according to plot.unit, and get left and right pos for rectangles to be plotted (either continuous
        # or 1 probe long):
        x <- adjustSegPos(
          chrom = sample.segments[, 2],
          char.arms = sample.segments[, 3],
          start = sample.segments[, 4],
          stop = sample.segments[, 5],
          type = "chromosome",
          unit = pos.unit,
          op = op
        )
        xleft <- x$use.start
        xright <- x$use.stop

        heat.col <- sapply(
          sample.segments[, 7],
          getCol,
          colors = cs$colors,
          intervals = cs$intervals
        )
        # Plot rectangles with appropriate color for each probe:
        ytop <- i - op$sep.samples
        ybottom <- i - (1 - op$sep.samples)
        rect(xleft, ybottom, xright, ytop, col = heat.col, border = NA)

        # Add sampleid on yaxis
        if (op$sample.labels) {
          axis(
            side = 2,
            at = (ytop - (ytop - ybottom) / 2),
            labels = sampleID[i],
            line = op$sample.line,
            tcl = 0,
            cex.axis = op$sample.cex,
            las = 1,
            mgp = op$mgp,
            tick = FALSE
          )
        }
      } # endfor

      if (!op$plot.ideo) {
        # Add xaxis:
        at.x <- getXticks(xlim[1], xlim[2], unit = op$plot.unit, ideal.n = 6)
        axis(
          side = 1,
          tcl = -0.2,
          at = at.x,
          cex.axis = op$cex.axis,
          mgp = op$mgp
        )
        title(xlab = op$xlab, cex.lab = op$cex.lab, line = op$mgp[1])
      }

      # Add box around plot:
      abline(v = xlim)
      abline(h = c(0, nSample))

      # If page is full; start plotting on new page
      if (c %% (nr * nc) == 0 && c != nChrom) {
        # Add main title to page:
        title(op$title[t], outer = TRUE)

        devAskNewPage(ask = TRUE)

        # Reset columns and row in layout:
        clm <- 1
        row <- 1
        new <- FALSE
      } else {
        # Update column and row index:
        if (clm < nc) {
          clm <- clm + 1
        } else {
          clm <- 1
          row <- row + 1
        } # endif
        new <- TRUE
      } # endif
    } # endfor

    # Add main title to page:
    title(op$title[t], outer = TRUE)
    if (t != nT) {
      devAskNewPage(ask = TRUE)
    }
  } # endfor
} # endfunction
