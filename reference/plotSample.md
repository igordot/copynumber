# Plot copy number data and/or segmentation results by sample

Plot copy number data and/or segmentation results for each sample
separately with chromosomes in different panels.

## Usage

``` r
plotSample(
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
)
```

## Arguments

- data:

  a data frame with numeric or character chromosome numbers in the first
  column, numeric local probe positions in the second, and numeric copy
  number data for one or more samples in subsequent columns. The header
  of the copy number columns should be the sample IDs.

- segments:

  a data frame or a list of data frames containing the segmentation
  results found by either
  [`pcf`](https://igordot.github.io/copynumber/reference/pcf.md) or
  [`multipcf`](https://igordot.github.io/copynumber/reference/multipcf.md).

- pos.unit:

  the unit used to represent the probe positions. Allowed options are
  "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs).
  By default assumed to be "bp".

- sample:

  a numeric vector indicating which sample(s) is (are) to be plotted.
  The number(s) should correspond to the sample's place (in order of
  appearance) in `data`, or in `segments` in case `data` is unspecified.

- chrom:

  a numeric or character vector with chromosome number(s) to indicate
  which chromosome(s) is (are) to be plotted.

- assembly:

  a string specifying which genome assembly version should be applied to
  define the chromosome ideogram. Allowed options are "hg19", "hg18",
  "hg17" and "hg16" (corresponding to the four latest human genome
  annotations in the UCSC genome browser).

- winsoutliers:

  an optional data frame of the same size as `data` identifying
  observations classified as outliers by
  [`winsorize`](https://igordot.github.io/copynumber/reference/winsorize.md).
  If specified, outliers will be marked by a different color and symbol
  than the other observations (see `wins.col` and `wins.pch`).

- xaxis:

  either "pos" or "index". The former implies that the xaxis will
  represent the genomic positions, whereas the latter implies that the
  xaxis will represent the probe index. Default is "pos".

- layout:

  an integer vector of length two giving the number of rows and columns
  in the plot. Default is `c(1,1)`.

- plot.ideo:

  a logical value indicating whether the chromosome ideogram should be
  plotted. Only applicable when `xaxis="pos"`.

- ...:

  other graphical parameters. These include the common plot arguments
  `xlab`, `ylab`, `main`, `xlim`, `ylim`, `col` (default is "grey"),
  `pch` (default is 46, equivalent to "."), `cex`, `cex.lab`,
  `cex.main`, `cex.axis`, `las`, `tcl`, `mar` and `mgp` (see
  [`par`](https://rdrr.io/r/graphics/par.html) on these). In addition, a
  range of graphical arguments specific for `plotSample` (as well as the
  similar functions `plotChrom`, `plotGenome` and `plotAllele`) may be
  specified:

## Details

Several plots may be produced on the same page with the `layout` option.
If the number of plots exceeds the desired page layout, the user is
prompted before advancing to the next page of output.

## Note

These functions apply `par(fig)`, and are therefore not compatible with
other setups for arranging multiple plots in one device such as
`par(mfrow,mfcol)`.

## Additional graphical parameters for `plotSample`

- `dir.print`::

  an optional directory where the plot(s) is (are) to be saved as pdf
  file(s). Defaults to NULL which implies that the plot(s) is (are)
  printed to screen instead.

- `file.name`:

  an optional character vector containing file name(s) for the pdf
  file(s) to be saved.

- `onefile`::

  logical value indicating whether all plots should be plotted in one
  device / saved in one file. Default is TRUE. If FALSE, a new window is
  opened or a new file is saved for each sample (each chromosome for
  `plotChrom`).

- `plot.size`::

  a numeric vector of length 2 giving the width and height of the
  plotting window. Default is `c(11.6,8.2)`.

- `title`::

  an overall title for all plots on one page.

- `plot.unit`::

  the desired unit to be applied for probe position tick marks along the
  x-axis. Only "mbp" (default) and "kbp" is allowed.

- `equalRange`::

  logical value indicating whether the range of the y-axis should be the
  same across all plots. Defaults to TRUE.

- `q`::

  a numerical value in the range 0 to 1 indicating that `ylim` will be
  set to only include observations between the (1-q/2)- and the
  (q/2)-quantile. Observations that fall outside these quantiles are
  truncated to the limits of the plot, and are by default marked by a
  special symbol (see `q.pch`). Default is `q=0.01` when `data` is
  specified, and `q=0` otherwise.

- `q.col`, `wins.col`::

  colors used to plot truncated observations and outliers. Default is
  "grey" and "magenta", respectively.

- `q.pch`, `wins.pch`::

  symbols used to plot truncated observations and outliers. Default is
  42 (equivalent to "\*") for both. Note that input must be of the same
  class as `pch` (numeric or character).

- `q.cex`, `wins.cex`::

  magnification used for truncated observations and outliers relative to
  `cex`. Default is 0.4 for both.

- `h`::

  a numerical value indicating that a horizontal reference line should
  be plotted at `y=h`. Default is `h=0`. `h=NULL` suppresses the
  plotting of a reference line.

- `at.x`::

  the points at which tick-marks on x-axis are to be drawn.

- `at.y`::

  the points at which tick-marks on y-axis are to be drawn.

- `main.line`::

  the margin line for the main title.

- `legend`::

  either a logical value indicating whether legends should be added to
  the plot if there is more than one segmentation result present in
  `segments`, or a character vector giving the legend texts to be used
  for the segmentation results. Default is TRUE, in which case the
  legend will be plotted in the topright corner of each plot.

- `seg.col`::

  color(s) used to plot the segmentation result(s). The default colors
  are found using the function `rainbow(n)`, where `n` is the number of
  segmentation results found in `segments` (see
  [`rainbow`](https://rdrr.io/r/grDevices/palettes.html) for details).

- `seg.lty`::

  the line type(s) used to plot the segmentation result(s). Default is
  1.

- `seg.lwd`::

  the line width(s) used to plot the segmentation result(s).

- `connect`::

  logical value indicating whether segments should be connected by
  vertical lines, default is TRUE.

- `ideo.frac`::

  a numerical value in the range 0 to 1 indicating the fraction of the
  plot to be occupied by the chromosome ideogram.

- `cyto.text`::

  a logical value indicating whether cytoband-names should be plotted
  along with the ideogram. Not recommended when many plots are plotted
  in the same grid, default is FALSE.

- `cex.cytotext`::

  the magnification used for the plotting of the cytoband-names.

- `cex.chrom`::

  the text size used to plot chromosome numbers in `plotGenome`.

## See also

[`plotChrom`](https://igordot.github.io/copynumber/reference/plotChrom.md),
[`plotGenome`](https://igordot.github.io/copynumber/reference/plotGenome.md)

## Author

Gro Nilsen

## Examples

``` r

#Lymphoma data
data(lymphoma)
#Take out a smaller subset of 6 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:6)

#Winsorize data:
wins.data <- winsorize(data=sub.lymphoma)
#> winsorize finished for chromosome arm 1p 
#> winsorize finished for chromosome arm 1q 
#> winsorize finished for chromosome arm 2p 
#> winsorize finished for chromosome arm 2q 
#> winsorize finished for chromosome arm 3p 
#> winsorize finished for chromosome arm 3q 
#> winsorize finished for chromosome arm 4p 
#> winsorize finished for chromosome arm 4q 
#> winsorize finished for chromosome arm 5p 
#> winsorize finished for chromosome arm 5q 
#> winsorize finished for chromosome arm 6p 
#> winsorize finished for chromosome arm 6q 
#> winsorize finished for chromosome arm 7p 
#> winsorize finished for chromosome arm 7q 
#> winsorize finished for chromosome arm 8p 
#> winsorize finished for chromosome arm 8q 
#> winsorize finished for chromosome arm 9p 
#> winsorize finished for chromosome arm 9q 
#> winsorize finished for chromosome arm 10p 
#> winsorize finished for chromosome arm 10q 
#> winsorize finished for chromosome arm 11p 
#> winsorize finished for chromosome arm 11q 
#> winsorize finished for chromosome arm 12p 
#> winsorize finished for chromosome arm 12q 
#> winsorize finished for chromosome arm 13q 
#> winsorize finished for chromosome arm 14q 
#> winsorize finished for chromosome arm 15q 
#> winsorize finished for chromosome arm 16p 
#> winsorize finished for chromosome arm 16q 
#> winsorize finished for chromosome arm 17p 
#> winsorize finished for chromosome arm 17q 
#> winsorize finished for chromosome arm 18p 
#> winsorize finished for chromosome arm 18q 
#> winsorize finished for chromosome arm 19p 
#> winsorize finished for chromosome arm 19q 
#> winsorize finished for chromosome arm 20p 
#> winsorize finished for chromosome arm 20q 
#> winsorize finished for chromosome arm 21q 
#> winsorize finished for chromosome arm 22q 
#> winsorize finished for chromosome arm 23p 
#> winsorize finished for chromosome arm 23q 

#Use pcf to find segments:        
uni.segments <- pcf(data=wins.data,gamma=12)
#> pcf finished for chromosome arm 1p 
#> pcf finished for chromosome arm 1q 
#> pcf finished for chromosome arm 2p 
#> pcf finished for chromosome arm 2q 
#> pcf finished for chromosome arm 3p 
#> pcf finished for chromosome arm 3q 
#> pcf finished for chromosome arm 4p 
#> pcf finished for chromosome arm 4q 
#> pcf finished for chromosome arm 5p 
#> pcf finished for chromosome arm 5q 
#> pcf finished for chromosome arm 6p 
#> pcf finished for chromosome arm 6q 
#> pcf finished for chromosome arm 7p 
#> pcf finished for chromosome arm 7q 
#> pcf finished for chromosome arm 8p 
#> pcf finished for chromosome arm 8q 
#> pcf finished for chromosome arm 9p 
#> pcf finished for chromosome arm 9q 
#> pcf finished for chromosome arm 10p 
#> pcf finished for chromosome arm 10q 
#> pcf finished for chromosome arm 11p 
#> pcf finished for chromosome arm 11q 
#> pcf finished for chromosome arm 12p 
#> pcf finished for chromosome arm 12q 
#> pcf finished for chromosome arm 13q 
#> pcf finished for chromosome arm 14q 
#> pcf finished for chromosome arm 15q 
#> pcf finished for chromosome arm 16p 
#> pcf finished for chromosome arm 16q 
#> pcf finished for chromosome arm 17p 
#> pcf finished for chromosome arm 17q 
#> pcf finished for chromosome arm 18p 
#> pcf finished for chromosome arm 18q 
#> pcf finished for chromosome arm 19p 
#> pcf finished for chromosome arm 19q 
#> pcf finished for chromosome arm 20p 
#> pcf finished for chromosome arm 20q 
#> pcf finished for chromosome arm 21q 
#> pcf finished for chromosome arm 22q 
#> pcf finished for chromosome arm 23p 
#> pcf finished for chromosome arm 23q 

#Use multipcf to find segments as well:
multi.segments <- multipcf(data=wins.data,gamma=12)
#> multipcf finished for chromosome arm 1p 
#> multipcf finished for chromosome arm 1q 
#> multipcf finished for chromosome arm 2p 
#> multipcf finished for chromosome arm 2q 
#> multipcf finished for chromosome arm 3p 
#> multipcf finished for chromosome arm 3q 
#> multipcf finished for chromosome arm 4p 
#> multipcf finished for chromosome arm 4q 
#> multipcf finished for chromosome arm 5p 
#> multipcf finished for chromosome arm 5q 
#> multipcf finished for chromosome arm 6p 
#> multipcf finished for chromosome arm 6q 
#> multipcf finished for chromosome arm 7p 
#> multipcf finished for chromosome arm 7q 
#> multipcf finished for chromosome arm 8p 
#> multipcf finished for chromosome arm 8q 
#> multipcf finished for chromosome arm 9p 
#> multipcf finished for chromosome arm 9q 
#> multipcf finished for chromosome arm 10p 
#> multipcf finished for chromosome arm 10q 
#> multipcf finished for chromosome arm 11p 
#> multipcf finished for chromosome arm 11q 
#> multipcf finished for chromosome arm 12p 
#> multipcf finished for chromosome arm 12q 
#> multipcf finished for chromosome arm 13q 
#> multipcf finished for chromosome arm 14q 
#> multipcf finished for chromosome arm 15q 
#> multipcf finished for chromosome arm 16p 
#> multipcf finished for chromosome arm 16q 
#> multipcf finished for chromosome arm 17p 
#> multipcf finished for chromosome arm 17q 
#> multipcf finished for chromosome arm 18p 
#> multipcf finished for chromosome arm 18q 
#> multipcf finished for chromosome arm 19p 
#> multipcf finished for chromosome arm 19q 
#> multipcf finished for chromosome arm 20p 
#> multipcf finished for chromosome arm 20q 
#> multipcf finished for chromosome arm 21q 
#> multipcf finished for chromosome arm 22q 
#> multipcf finished for chromosome arm 23p 
#> multipcf finished for chromosome arm 23q 

# plotSample() opens a new graphics device on its first call (needed for
# its multi-page/multi-window output modes), which is not safe to run
# under automated example checks (e.g. pkgdown's example capture).
if (FALSE) { # \dontrun{
#Plot data and pcf-segments for one sample separately for each chromosome:
plotSample(data=sub.lymphoma,segments=uni.segments,sample=1,layout=c(5,5))
#Add cytoband text to ideogram (one page per chromosome to ensure sufficient
#space)
plotSample(data=sub.lymphoma,segments=uni.segments,sample=1,layout=c(1,1),
    cyto.text=TRUE)
#Add multipcf-segmentation results, drop legend
plotSample(data=sub.lymphoma,segments=list(uni.segments,multi.segments),sample=1,
    layout=c(5,5),seg.col=c("red","blue"),seg.lwd=c(3,2),legend=FALSE)
#Plot by chromosome for two samples, but only chromosome 1-9. One window per
#sample:
plotSample(data=sub.lymphoma,segments=list(uni.segments,multi.segments),sample=
    c(2,3),chrom=c(1:9),layout=c(3,3),seg.col=c("red","blue"),
    seg.lwd=c(3,2),onefile=FALSE)
#Zoom in on a particular region by setting xlim:
plotSample(data=sub.lymphoma,segments=uni.segments,sample=1,chrom=1,plot.ideo=
    FALSE,xlim=c(140,170))
} # }
```
