# Interpolation of pcf-estimates.

Given a segmentation by `pcf`, interpolate pcf-estimates for specific
positions.

## Usage

``` r
# S3 method for class 'pcf'
interpolate(segments, x)
```

## Arguments

- segments:

  a data frame containing the segmentation results from
  [`pcf`](https://igordot.github.io/copynumber/reference/pcf.md).

- x:

  matrix or data.frame where the first column gives chrosomomes and the
  second gives positions.

## Value

A data frame where the first two columns give the chromsomes and
positions specified in the input `x` and the remaining columns give the
interpolated pcf-estimate for each sample represented in `segments`.

## Details

Pcf-estimates are interpolated for the chromosomes and postions
specified in `x`.

## Note

The positions in `segments` and `x` must be of the same unit (bp, kbp,
or mbp).

## See also

[`pcf`](https://igordot.github.io/copynumber/reference/pcf.md)

## Author

Gro Nilsen, Ole Christian Lingjaerde.

## Examples

``` r

#Load the lymphoma data set:
data(lymphoma)

#Take out a smaller subset of 3 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:3)

#Run pcf:
seg <- pcf(data=sub.lymphoma,gamma=12)
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

#Make a matrix with two positions and chromosomes for which we want to 
#interpolate the pcf-estimate:
pos <-  c(2000000,50000000)
chr <- c(1,2)
x <- cbind(chr,pos)

#Interpolate
int.pcf <- interpolate.pcf(seg,x)

```
