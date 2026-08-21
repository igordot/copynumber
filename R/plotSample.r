####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# Function to plot observed copy number data and/or segmentation result for one given sample with a separate figure for each chromosome

# Input:
### data: dataframe or matrix with chromosomes in first column, positions in second column and copy number data for one or more samples in subsequent columns
### segments: a data frame or a list with data frames containing segmentation results
### pos.unit: the unit used to represent the positions in data
### sample: a numeric vector indicating which samples are to be plotted (numerated such that sample=1 indicates the first sample found in data)
### chrom: numeric/character vector giving the chromosomes to be plotted
### xaxis: what is to be plotted along xaxis; either positions ("pos") or indeces ("index")
### layout: the grid layout for the plot (number of columns and rows)
### plot.ideo: should ideogram be plotted below plots
### ... : other optional plot parameters

## Required by:
### none

## Requires:
### checkAndRetrievePlotInput
### chromMax
### framedim
### getSeglim
### getFilename
### getPlotParameters
### plotIdeogram
### plotObs
### plotSegments



#' Plot copy number data and/or segmentation results by sample
#' 
#' Plot copy number data and/or segmentation results for each sample separately
#' with chromosomes in different panels.
#' 
#' Several plots may be produced on the same page with the \code{layout}
#' option. If the number of plots exceeds the desired page layout, the user is
#' prompted before advancing to the next page of output.
#' 
#' @param data a data frame with numeric or character chromosome numbers in the
#' first column, numeric local probe positions in the second, and numeric copy
#' number data for one or more samples in subsequent columns. The header of the
#' copy number columns should be the sample IDs.
#' @param segments a data frame or a list of data frames containing the
#' segmentation results found by either \code{\link{pcf}} or
#' \code{\link{multipcf}}.
#' @param pos.unit the unit used to represent the probe positions. Allowed
#' options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base
#' pairs). By default assumed to be "bp".
#' @param sample a numeric vector indicating which sample(s) is (are) to be
#' plotted. The number(s) should correspond to the sample's place (in order of
#' appearance) in \code{data}, or in \code{segments} in case \code{data} is
#' unspecified.
#' @param chrom a numeric or character vector with chromosome number(s) to
#' indicate which chromosome(s) is (are) to be plotted.
#' @param assembly a string specifying which genome assembly version should be
#' applied to define the chromosome ideogram. Allowed options are "hg19",
#' "hg18", "hg17" and "hg16" (corresponding to the four latest human genome
#' annotations in the UCSC genome browser).
#' @param winsoutliers an optional data frame of the same size as \code{data}
#' identifying observations classified as outliers by \code{\link{winsorize}}.
#' If specified, outliers will be marked by a different color and symbol than
#' the other observations (see \code{wins.col} and \code{wins.pch}).
#' @param xaxis either "pos" or "index". The former implies that the xaxis will
#' represent the genomic positions, whereas the latter implies that the xaxis
#' will represent the probe index. Default is "pos".
#' @param layout an integer vector of length two giving the number of rows and
#' columns in the plot. Default is \code{c(1,1)}.
#' @param plot.ideo a logical value indicating whether the chromosome ideogram
#' should be plotted. Only applicable when \code{xaxis="pos"}.
#' @param \dots other graphical parameters. These include the common plot arguments
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}, \code{col}
#' (default is "grey"), \code{pch} (default is 46, equivalent to "."),
#' \code{cex}, \code{cex.lab}, \code{cex.main}, \code{cex.axis}, \code{las},
#' \code{tcl}, \code{mar} and \code{mgp} (see \code{\link{par}} on these). In
#' addition, a range of graphical arguments specific for \code{plotSample} (as
#' well as the similar functions \code{plotChrom}, \code{plotGenome} and
#' \code{plotAllele}) may be specified:
#' @section Additional graphical parameters for \code{plotSample}:
#' \describe{
#' \item{\code{dir.print}:}{an optional directory where the plot(s) is (are)
#' to be saved as pdf file(s). Defaults to NULL which implies that the plot(s)
#' is (are) printed to screen instead.} \item{\code{file.name}}{an optional
#' character vector containing file name(s) for the pdf file(s) to be saved.}
#' \item{\code{onefile}:}{logical value indicating whether all plots should be
#' plotted in one device / saved in one file. Default is TRUE. If FALSE, a new
#' window is opened or a new file is saved for each sample (each chromosome
#' for \code{plotChrom}).} \item{\code{plot.size}:}{a numeric vector of length
#' 2 giving the width and height of the plotting window. Default is
#' \code{c(11.6,8.2)}.} \item{\code{title}:}{an overall title for all plots on
#' one page.} \item{\code{plot.unit}:}{the desired unit to be applied for
#' probe position tick marks along the x-axis. Only "mbp" (default) and "kbp"
#' is allowed.} \item{\code{equalRange}:}{logical value indicating whether the
#' range of the y-axis should be the same across all plots. Defaults to TRUE.}
#' \item{\code{q}:}{a numerical value in the range 0 to 1 indicating that
#' \code{ylim} will be set to only include observations between the (1-q/2)-
#' and the (q/2)-quantile. Observations that fall outside these quantiles are
#' truncated to the limits of the plot, and are by default marked by a special
#' symbol (see \code{q.pch}). Default is \code{q=0.01} when \code{data} is
#' specified, and \code{q=0} otherwise.} \item{\code{q.col},
#' \code{wins.col}:}{colors used to plot truncated observations and outliers.
#' Default is "grey" and "magenta", respectively.} \item{\code{q.pch},
#' \code{wins.pch}:}{symbols used to plot truncated observations and outliers.
#' Default is 42 (equivalent to "*") for both. Note that input must be of the
#' same class as \code{pch} (numeric or character).} \item{\code{q.cex},
#' \code{wins.cex}:}{magnification used for truncated observations and
#' outliers relative to \code{cex}. Default is 0.4 for both.}
#' \item{\code{h}:}{a numerical value indicating that a horizontal reference
#' line should be plotted at \code{y=h}. Default is \code{h=0}. \code{h=NULL}
#' suppresses the plotting of a reference line.} \item{\code{at.x}:}{the
#' points at which tick-marks on x-axis are to be drawn.}
#' \item{\code{at.y}:}{the points at which tick-marks on y-axis are to be
#' drawn.} \item{\code{main.line}:}{the margin line for the main title.}
#' \item{\code{legend}:}{either a logical value indicating whether legends
#' should be added to the plot if there is more than one segmentation result
#' present in \code{segments}, or a character vector giving the legend texts
#' to be used for the segmentation results. Default is TRUE, in which case the
#' legend will be plotted in the topright corner of each plot.}
#' \item{\code{seg.col}:}{color(s) used to plot the segmentation result(s).
#' The default colors are found using the function \code{rainbow(n)}, where
#' \code{n} is the number of segmentation results found in \code{segments}
#' (see \code{\link{rainbow}} for details).} \item{\code{seg.lty}:}{the line
#' type(s) used to plot the segmentation result(s). Default is 1.}
#' \item{\code{seg.lwd}:}{the line width(s) used to plot the segmentation
#' result(s).} \item{\code{connect}:}{logical value indicating whether
#' segments should be connected by vertical lines, default is TRUE.}
#' \item{\code{ideo.frac}:}{a numerical value in the range 0 to 1 indicating
#' the fraction of the plot to be occupied by the chromosome ideogram.}
#' \item{\code{cyto.text}:}{a logical value indicating whether cytoband-names
#' should be plotted along with the ideogram. Not recommended when many plots
#' are plotted in the same grid, default is FALSE.}
#' \item{\code{cex.cytotext}:}{the magnification used for the plotting of the
#' cytoband-names.} \item{\code{cex.chrom}:}{the text size used to plot
#' chromosome numbers in \code{plotGenome}.}
#' }
#' @note These functions apply \code{par(fig)}, and are therefore not
#' compatible with other setups for arranging multiple plots in one device such
#' as \code{par(mfrow,mfcol)}.
#' @author Gro Nilsen
#' @seealso \code{\link{plotChrom}}, \code{\link{plotGenome}}
#' @examples
#' 
#' #Lymphoma data
#' data(lymphoma)
#' #Take out a smaller subset of 6 samples (using subsetData):
#' sub.lymphoma <- subsetData(lymphoma,sample=1:6)
#' 
#' #Winsorize data:
#' wins.data <- winsorize(data=sub.lymphoma)
#' 
#' #Use pcf to find segments:        
#' uni.segments <- pcf(data=wins.data,gamma=12)
#' 
#' #Use multipcf to find segments as well:
#' multi.segments <- multipcf(data=wins.data,gamma=12)
#' 
#' # plotSample() opens a new graphics device on its first call (needed for
#' # its multi-page/multi-window output modes), which is not safe to run
#' # under automated example checks (e.g. pkgdown's example capture).
#' \dontrun{
#' #Plot data and pcf-segments for one sample separately for each chromosome:
#' plotSample(data=sub.lymphoma,segments=uni.segments,sample=1,layout=c(5,5))
#' #Add cytoband text to ideogram (one page per chromosome to ensure sufficient
#' #space)
#' plotSample(data=sub.lymphoma,segments=uni.segments,sample=1,layout=c(1,1),
#'     cyto.text=TRUE)
#' #Add multipcf-segmentation results, drop legend
#' plotSample(data=sub.lymphoma,segments=list(uni.segments,multi.segments),sample=1,
#'     layout=c(5,5),seg.col=c("red","blue"),seg.lwd=c(3,2),legend=FALSE)
#' #Plot by chromosome for two samples, but only chromosome 1-9. One window per
#' #sample:
#' plotSample(data=sub.lymphoma,segments=list(uni.segments,multi.segments),sample=
#'     c(2,3),chrom=c(1:9),layout=c(3,3),seg.col=c("red","blue"),
#'     seg.lwd=c(3,2),onefile=FALSE)
#' #Zoom in on a particular region by setting xlim:
#' plotSample(data=sub.lymphoma,segments=uni.segments,sample=1,chrom=1,plot.ideo=
#'     FALSE,xlim=c(140,170))
#' }
#'
#' @export
plotSample <- function(
  data = NULL,
  segments = NULL,
  pos.unit = "bp",
  sample = NULL,
  chrom = NULL,
  assembly = "hg19",
  winsoutliers = NULL,
  xaxis = "pos",
  layout = c(1, 1),
  plot.ideo = TRUE,
  ...
) {
  # Check, modify and retrieve plot input:
  input <- checkAndRetrievePlotInput(
    data = data,
    segments = segments,
    winsoutliers = winsoutliers,
    type = "sample",
    xaxis = xaxis,
    pos.unit = pos.unit,
    sample = sample,
    chrom = chrom
  )
  data <- input$data
  segments <- input$segments
  sampleID <- input$sampleID
  chrom <- input$chrom
  winsoutliers <- input$winsoutliers

  nSample <- length(sampleID)
  nChrom <- length(chrom)
  nSeg <- length(segments) # will be 0 if segments=NULL
  sample.names <- colnames(data)[-c(1:2)] # will be NULL if data=NULL

  # Plot layout (number of columns and rows in plot grid) specified by user:
  nr <- layout[1]
  nc <- layout[2]

  # If xaxis is index; cannot plot ideogram:
  if (xaxis == "index") {
    plot.ideo <- FALSE
  }
  # Set default plot parameters and change these if user has specified other choices via ... :
  arg <- getPlotParameters(
    type = "sample",
    nSeg = nSeg,
    cr = nc * nr,
    sampleID = sampleID,
    plot.ideo = plot.ideo,
    xaxis = xaxis,
    assembly = assembly,
    ...
  )

  # Margins used for the plot window:
  if (arg$title[1] == "") {
    oma <- c(0, 0, 0, 0)
  } else {
    oma <- c(0, 0, 1, 0)
  }
  mar <- c(0.2, 0.2, 0.3, 0.2)

  # Check if there should be more than one file/window with plot(s), and get file.name accordingly
  nPage <- ifelse(arg$onefile, 1, nSample)
  file.name <- getFilename(nPage, arg$file.name, ID = sampleID, type = "sample")

  # Divide the plotting window by the function "framedim":
  frames <- framedim(nr, nc)

  # Separate plots for each sample:

  for (i in 1:nSample) {
    # Start new window/file:
    if (!arg$onefile || i == 1) {
      # Either print to file, or plot on screen
      if (!is.null(arg$dir.print)) {
        pdf(
          file = paste(arg$dir.print, "/", file.name[i], ".pdf", sep = ""),
          width = arg$plot.size[1],
          height = arg$plot.size[2],
          onefile = TRUE,
          paper = "a4"
        ) # a4-paper
      } else {
        if (dev.cur() <= i) {
          # to make Sweave work
          dev.new(
            width = arg$plot.size[1],
            height = arg$plot.size[2],
            record = TRUE
          )
        }
      }
    } else {
      # Start new page when prompted by user:
      if (is.null(arg$dir.print)) {
        devAskNewPage(ask = TRUE)
      }
    }

    # Initialize row and column index:
    row <- 1
    clm <- 1
    new <- FALSE

    # Which sample is this:
    id <- sampleID[i]
    ind.sample <- which(sample.names == id) # will be integer(0) if data=NULL

    # Get data limits for this sample if range should include all chromosomes (equalRange=TRUE)
    if (!is.null(data) && arg$equalRange) {
      all.chrom <- which(data[, 1] %in% chrom)
      data.lim <- quantile(
        data[all.chrom, ind.sample + 2],
        probs = c(arg$q / 2, (1 - arg$q / 2)),
        names = FALSE,
        type = 4,
        na.rm = TRUE
      )
    }
    # Picking out all segments where sampleid (in first column) is id, returns new list
    if (!is.null(segments)) {
      sample.segments <- lapply(
        segments,
        function(seg, id) {
          seg[seg[, 1] == id, ]
        },
        id = id
      )
    }

    # Make separate plots for each chromosome:

    for (c in 1:nChrom) {
      # Frame dimensions for plot c:
      fig.c <- c(
        frames$left[clm],
        frames$right[clm],
        frames$bot[row],
        frames$top[row]
      )
      par(fig = fig.c, new = new, oma = oma, mar = mar)
      frame.c <- list(
        left = frames$left[clm],
        right = frames$right[clm],
        bot = frames$bot[row],
        top = frames$top[row]
      )

      # Divide frame for this chromosome into a frame for the actual plot and a frame for the ideogram:
      plot.frame <- frame.c
      plot.frame$bot <- frame.c$bot +
        (frame.c$top - frame.c$bot) * arg$ideo.frac
      ideo.frame <- frame.c
      ideo.frame$top <- plot.frame$bot

      # Select relevant chromosome number
      k <- chrom[c]

      if (!is.null(segments)) {
        # Get min and max values in segments to make sure all are shown in plot
        seg.lim <- sapply(
          sample.segments,
          getSeglim,
          equalRange = arg$equalRange,
          k = k
        ) # matrix with limits for each segmentation for this sample, min in row 1, max in row2
        seg.lim <- c(
          min(seg.lim[1, ], na.rm = TRUE),
          max(seg.lim[2, ], na.rm = TRUE)
        ) # Get overall min and max over all segments
      } else {
        seg.lim <- NULL
      }

      # Get maximum position on chromosome from cytoband info
      if (plot.ideo) {
        xmax <- chromMax(
          chrom = k,
          cyto.data = arg$assembly,
          pos.unit = arg$plot.unit
        )
        # PLOT IDEOGRAM
        par(fig = unlist(ideo.frame), new = new, mar = arg$mar.i)
        plotIdeogram(
          chrom = k,
          arg$cyto.text,
          cyto.data = arg$assembly,
          cex = arg$cex.cytotext,
          unit = arg$plot.unit
        )
        new <- TRUE
      } else {
        xmax <- NULL
      }

      # PLOT DATA POINTS:
      add <- FALSE
      if (!is.null(data)) {
        ind.chrom <- which(data[, 1] == k)
        if (!arg$equalRange) {
          # Get data limits for this sample using just this chromosome (equalRange=FALSE)
          data.lim <- quantile(
            data[ind.chrom, ind.sample + 2],
            probs = c(arg$q / 2, (1 - arg$q / 2)),
            names = FALSE,
            type = 4,
            na.rm = TRUE
          )
        }

        plotObs(
          y = data[ind.chrom, ind.sample + 2],
          pos = data[ind.chrom, 2],
          unit = pos.unit,
          winsoutliers = winsoutliers[ind.chrom, ind.sample + 2],
          type = "sample",
          xaxis = xaxis,
          plot.ideo = plot.ideo,
          k = k,
          frame = plot.frame,
          new = new,
          op = arg,
          data.lim = data.lim,
          seg.lim = seg.lim,
          xmax = xmax
        )

        add <- TRUE ## segment plot will be added on top of data plot
      } else {
        data.lim <- NULL
      }

      # Plot segments:
      if (!is.null(segments)) {
        # Plot all segments in list:
        for (s in 1:nSeg) {
          use.segments <- sample.segments[[s]]
          plotSegments(
            use.segments[use.segments[, 2] == k, , drop = FALSE],
            type = "sample",
            k = k,
            xaxis = xaxis,
            add = add,
            plot.ideo = plot.ideo,
            col = arg$seg.col[s],
            lty = arg$seg.lty[s],
            lwd = arg$seg.lwd[s],
            frame = plot.frame,
            new = new,
            unit = pos.unit,
            seg.lim = seg.lim,
            data.lim = data.lim,
            op = arg
          )
          if (k %in% use.segments[, 2]) {
            add <- TRUE
          }
        } # endfor

        # Add segmentation legends:
        if (!is.null(arg$legend)) {
          legend(
            "topright",
            legend = arg$legend,
            col = arg$seg.col,
            lty = arg$seg.lty,
            cex = arg$cex.axis
          )
        }
      }

      # If page is full; plot on new page
      if (c %% (nr * nc) == 0) {
        # Add main title to page:
        title(arg$title[i], outer = TRUE)

        # Start new page when prompted by user:
        if (is.null(arg$dir.print)) {
          devAskNewPage(ask = TRUE)
        }

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

    # Plot sampleid as title

    title(arg$title[i], outer = TRUE)

    # Close graphcis:
    if (!is.null(arg$dir.print)) {
      if (!arg$onefile) {
        cat(
          "Plot was saved in ",
          paste(arg$dir.print, "/", file.name[i], ".pdf", sep = ""),
          "\n"
        )
        graphics.off()
      } else {
        if (i == nSample) {
          cat(
            "Plot was saved in ",
            paste(arg$dir.print, "/", file.name, ".pdf", sep = ""),
            "\n"
          )
          graphics.off()
        }
      }
    } # endif
  } # endfor
} # endfunction
