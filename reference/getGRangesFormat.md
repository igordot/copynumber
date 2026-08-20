# Get segments on the GRanges format

The segments data frame obtained e.g. by `pcf`, `multipcf` or `aspcf` is
converted to the GRanges format.

## Usage

``` r
getGRangesFormat(segments)
```

## Arguments

- segments:

  a data frame containing segmentation results found by e.g.
  [`pcf`](https://igordot.github.io/copynumber/reference/pcf.md),
  [`multipcf`](https://igordot.github.io/copynumber/reference/multipcf.md)
  or [`aspcf`](https://igordot.github.io/copynumber/reference/aspcf.md).

## Value

The segments converted to the GRanges container class.

## Details

GRanges, in the GenomicRanges package, is the standard BioConductor
containers for range data. For some applications it may therefore be
useful to convert segmentation results to this format.

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
#Obtain the GRanges format
gr <- getGRangesFormat(seg)
```
