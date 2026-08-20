# Impute missing copy number values

Missing copy number values are imputed by a constant value or
pcf-estimates.

## Usage

``` r
imputeMissing(data, method, c = 0, pcf.est = NULL, ...)
```

## Arguments

- data:

  a data frame with numeric or character chromosome numbers in the first
  column, numeric local probe positions in the second, and numeric copy
  number data for one or more samples in subsequent columns.

- method:

  the imputation method to be used. Must be one of "constant" and "pcf".

- c:

  a numerical value to be imputed if method is "constant". Default is 0.

- pcf.est:

  a data frame of same size as `data`, with chromosome numbers and
  positions in the first two columns, and copy number estimates obtained
  from `pcf` in the subsequent columns. Only applicable if
  `method="pcf"`. If unspecified and `method="pcf"`, `pcf` is run
  internally to find estimates.

- ...:

  other relevant parameters to be passed on to `pcf`

## Value

A data frame of the same size and format as `data` with all missing
values imputed.

## Details

The available imputation methods are:

- `constant`::

  all missing values in `data` are replaced by the specified value `c`.

- `pcf`::

  the estimates from pcf-segmentation (see
  [`pcf`](https://igordot.github.io/copynumber/reference/pcf.md)) are
  used to impute missing values. If `pcf` has already been run, these
  estimates may be specified in `pcf.est`. If `pcf.est` is unspecified,
  `pcf` is run on the input data. In `pcf` the analysis is done on the
  observed values, and estimates for missing observations are set to be
  the estimate of the nearest observed probe.

## See also

[`pcf`](https://igordot.github.io/copynumber/reference/pcf.md)

## Author

Gro Nilsen

## Examples

``` r

#Load lymphoma data
data(lymphoma)
chrom <- lymphoma[,1]
pos <- lymphoma[,2]
#pick out data for the first six samples:
cn.data <- lymphoma[,3:8]

#Create missing values in cn.data at random positions:
n <- nrow(cn.data)*ncol(cn.data)
r <- matrix(rbinom(n=n,size=1,prob=0.95),nrow=nrow(cn.data),ncol=ncol(cn.data))
cn.data[r==0] <- NA    #matrix with approximately 5% missing values
mis.data <- data.frame(chrom,pos,cn.data)

#Impute missing values by constant, c=0:
imp.data <- imputeMissing(data=mis.data,method="constant")

#Impute missing values by obtained pcf-values:
pcf.est <- pcf(data=mis.data,return.est=TRUE)
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
imp.data <- imputeMissing(data=mis.data,method="pcf",pcf.est=pcf.est)

#Or run pcf within imputeMissing:
imp.data <- imputeMissing(data=mis.data,method="pcf")
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

```
