# Retrieve a subset of segments

This function returns a subset of segments according to the input and
the specified chromosomes and/or samples.

## Usage

``` r
subsetSegments(segments, chrom = NULL, sample = NULL, sep = "\t", ...)
```

## Arguments

- segments:

  either a data frame or the name of a tab-separated file from which
  segmentation results can be read. Segmentation results may come from
  [`pcf`](https://igordot.github.io/copynumber/reference/pcf.md),
  [`multipcf`](https://igordot.github.io/copynumber/reference/multipcf.md)
  or [`aspcf`](https://igordot.github.io/copynumber/reference/aspcf.md).

- chrom:

  a numeric or character vector with chromosome(s) for which segments
  should be selected. If unspecified, all chromosomes in `segments` will
  be selected.

- sample:

  a numeric vector indicating for which sample(s) segments should be
  selected. The number(s) should correspond to the sample's place (in
  order of appearance) in `segments`.

- sep:

  the separator of the input files if `segments` is a file. Default is
  tab.

- ...:

  optional parameters to be passed to `read.table` in the case where
  `segments` is to be read from file.

## Value

A data frame containing the desired subset of segments.

## Author

Gro Nilsen

## Examples

``` r

#Load lymphoma data
data(lymphoma)

#Select segments only for samples 1 and 6 and chromosomes 1:9:
segments <- pcf(lymphoma,gamma=12)
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
sub.segments <- subsetSegments(segments=segments,chrom=c(1:9),sample=c(1,6))

```
