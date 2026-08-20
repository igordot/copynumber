# Call aberrations in segmented data

Segments, obtained by `pcf` or `multipcf`, are classified as "gain",
"normal" or "loss" given the specified thresholds.

## Usage

``` r
callAberrations(segments, thres.gain, thres.loss = -thres.gain)
```

## Arguments

- segments:

  a data frame containing the segmentation results found by either
  [`pcf`](https://igordot.github.io/copynumber/reference/pcf.md) or
  [`multipcf`](https://igordot.github.io/copynumber/reference/multipcf.md).

- thres.gain:

  a numeric value giving the threshold to be applied for calling gains.

- thres.loss:

  a numeric value giving the threshold to be applied for calling losses.
  Default is to use the negative value of `thres.gain`.

## Value

A new segment data frame where the segment values have been replaced by
the classification "gain", "normal" or "loss".

## Details

Each region found in `segments` is classified as "gain", "normal" or
"loss". Regions with gain or loss will be those segments where the
segment value is above or below the value given in `thres.gain` or
`thres.loss`, respectively.

## Author

Gro Nilsen

## Examples

``` r

#load lymphoma data
data(lymphoma)
#Run pcf
seg <- pcf(data=lymphoma,gamma=12)
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

#Call gains as segments whose value is > 0.2, and losses as segments whose
# value < -0.1
ab.seg <- callAberrations(seg,thres.gain=0.2,thres.loss=-0.1)

```
