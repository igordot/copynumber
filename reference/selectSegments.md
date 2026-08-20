# Select multipcf segments

Selects multipcf segments based on a desired characteristic.

## Usage

``` r
selectSegments(
  segments,
  what = "variance",
  thres = NULL,
  nseg = 10,
  large = TRUE,
  p = 0.1
)
```

## Arguments

- segments:

  a data frame containing segments found by
  [`multipcf`](https://igordot.github.io/copynumber/reference/multipcf.md).

- what:

  the desired characteristic to base selection on. Must be one of
  "variance" (default),"length" and "aberration". See details below.

- thres:

  an optional numeric threshold to be applied in the selection.

- nseg:

  the desired number of segments to be selected, default is 10. Only
  used if `thres=NULL`.

- large:

  logical value indicating whether segments with large (TRUE) or small
  (FALSE) variance, length or mean value should be selected when `what`
  is "variance", "length" or "aberration", respectively.

- p:

  a number between 0 and 1 giving the minimum proportion of samples for
  which an aberration must be detected, default is 0.1. Only applicable
  if `what="aberration"`.

## Value

A list containing:

- sel.seg:

  data frame containing the selected segments.

In addition, depending on the value of `what`:

- seg.var:

  a vector giving the variance for each segment. Only returned if
  `what = "variance"`.

- seg.length:

  a vector giving the length of each segment. Only returned if
  `what = "length"`.

- seg.ab.prop:

  a vector giving the aberration proportion for each segment. Only
  returned if `what = "aberration"` and `thres` is specified.

- seg.quantile:

  a vector giving the (1-p)- or p-quantile for each segment. Only
  returned if `what = "aberration"` and `thres=NULL`.

## Details

The input in `what` determines how the segments are selected. Three
options are available:

If `what="variance"` the variance of the segment values across all
samples is calculated for each segment. If `thres` is specified, the
subset of segments for which the variance is above (if `large=TRUE`) or
below (if `large=FALSE`) the threshold is returned. If `thres` is not
given by the user, a given number of segments determined by the input in
`nseg` is selected; if `large=TRUE` this will be the `nseg` segments
with the highest variance, whereas if `large=FALSE` the subset will
consist of the `nseg` segments with the lowest variance.

If `what="length"` selection is based on the genomic length of the
segments (end position minus start position). If `thres` is specified,
the subset of segments for which the length is above (if `large=TRUE`)
or below (if `large=FALSE`) this threshold is returned. If `thres` is
left unspecified, a given number of segments determined by the input in
`nseg` is selected; if `large=TRUE` this will be the `nseg` longest
segments, whereas if `large=FALSE` it will be the `nseg` shortest
segments.

If `what="aberration"` the aberration frequency is used to select the
subset of segments. If `thres` is specified, the proportion of samples
for which the segment value is above (if `large=TRUE`) or below (if
`large=FALSE`) the threshold is calculated for each segment. The subset
of segments where this frequency is above or equal to the proportion set
by the parameter `p` is returned. If `thres` is not specified, the
`nseg` segments with the highest (1-p)-quantile (if `large=TRUE`) or the
lowest p-quantile (if `large=FALSE`) is returned.

## See also

[`multipcf`](https://igordot.github.io/copynumber/reference/multipcf.md)

## Author

Gro Nilsen

## Examples

``` r

#Lymphoma data
data(lymphoma)

#Run multipcf
segments <- multipcf(lymphoma,gamma=12)
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

#Select the 10 segments with the highest variance:
sel.seg1 <- selectSegments(segments)

#Select the segments where the variance is below 0.001
sel.seg2 <- selectSegments(segments,thres=0.001,large=FALSE)

#Select the 5 longest segments:
sel.seg3 <- selectSegments(segments,what="length",nseg=5)

#Select the segments where 20 % of the samples have segment value of 0.2 or more:
sel.seg4 <- selectSegments(segments,what="aberration",thres=0.2,p=0.2)

#Select the 20 segments with the largest median:
sel.seg5 <- selectSegments(segments,what="aberration",nseg=20,p=0.5)
```
