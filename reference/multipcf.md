# Multi-sample copy number segmentation.

Joint segmentation resulting in piecewise constant curves with common
break points for all samples.

## Usage

``` r
multipcf(
  data,
  pos.unit = "bp",
  arms = NULL,
  Y = NULL,
  gamma = 40,
  normalize = TRUE,
  w = 1,
  fast = TRUE,
  assembly = "hg19",
  digits = 4,
  return.est = FALSE,
  save.res = FALSE,
  file.names = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  either a data frame or the name of a tab-separated file from which
  copy number data can be read. The rows of the data frame or file
  should represent the probes. Column 1 must hold numeric or character
  chromosome numbers, column 2 the numeric local probe positions, and
  subsequent columns the numeric copy number measurements for two or
  more samples. The header of copy number columns should give sample
  IDs.

- pos.unit:

  the unit used to represent the probe positions. Allowed options are
  "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs).
  By default assumed to be "bp".

- arms:

  optional character vector containing chromosome arms (denoted 'p' and
  'q') corresponding to the chromosomes and positions found in `data`.
  If not specified chromosome arms are found using the built-in genome
  assembly version determined by `assembly`.

- Y:

  either a data frame or the name of a tab-separated file containing
  original copy number data in the case where `data` contains Winsorized
  values. If provided, these values are used to calculate the mean of
  each segment, otherwise the copy number values in `data` are used. `Y`
  must be on the same form as `data`.

- gamma:

  penalty for each discontinuity in the curve, default is 40.

- normalize:

  a logical value indicating whether each sample's copy number
  measurements should be scaled by the sample specific residual standard
  error. Default is TRUE.

- w:

  a numeric vector giving an individual weight to be used for each
  sample. May be of length 1 if the same weight should be applied for
  each sample, default is 1 (no weighting).

- fast:

  a logical value indicating whether a fast (not guaranteed to be exact)
  version should be run on chromosome arms with \> 400 probes. Default
  is TRUE.

- assembly:

  a string specifying which genome assembly version should be applied to
  determine chromosome arms. Allowed options are "hg19", "hg18", "hg17"
  and "hg16" (corresponding to the four latest human genome annotations
  in the UCSC genome browser).

- digits:

  the number of decimals to be applied when reporting results. Default
  is 4.

- return.est:

  logical value indicating whether a data frame with copy number
  estimates (multipcf estimates)should be returned along with the
  segments. Default is FALSE, which means that only segments are
  returned.

- save.res:

  logical value indicating whether results should be saved in text
  files, default is FALSE.

- file.names:

  optional character vector of length two giving the name of the files
  where the multipcf estimates and segments, respectively, should be
  saved in case `save.res = TRUE`.

- verbose:

  logical value indicating whether or not to print a progress message
  each time multipcf analysis is finished for a new chromosome arm.

## Value

If `return.est = TRUE` a list with the following components:

- estimates:

  a data frame where the first two columns give the chromosome numbers
  and probe positions, respectively, while subsequent columns give the
  copy number estimates for each sample. The estimate for a given probe
  and sample equals the sample mean of the segment where the probe is
  located.

- segments:

  a data frame describing the segments found in the data. Each row
  represents a segment, and the first five columns give the chromosome
  numbers, arms, local start positions, local end positions, and the
  number of probes in the segments, respectively. Subsequent columns
  give the mean segment value for each sample, with sample IDs as column
  headers.

If `return.est = FALSE` only the data frame containing the segments is
returned.

If `save.res = TRUE` the results are also saved in text files with names
as specified in `file.names`. If `file.names=NULL`, a folder named
"multipcf_results" is created in the working directory, and the segments
and copy number estimates are saved in this folder as tab-separated
files named segments.txt and estimates.txt, respectively.

## Details

Piecewise constant curves are simultaneously fitted to the copy number
data for several samples as described in the multiPCF algorithm in
Nilsen and Liestoel et al. (2012). This implies that break points will
be the same for all segmentation curves, but the mean segment values
will differ among samples. Segmentation is done separately on each
chromosome arm.

## Note

It is usually advisable to Winsorize data before running pcf, see
[`winsorize`](https://igordot.github.io/copynumber/reference/winsorize.md)
on this.

The input data must be complete, see
[`imputeMissing`](https://igordot.github.io/copynumber/reference/imputeMissing.md)
for imputation of missing copy number values.

## References

Nilsen and Liestoel et al., "Copynumber: Efficient algorithms for
single- and multi-track copy number segmentation", BMC Genomics 13:591
(2012), doi:10.1186/1471-2164-13-59

## See also

[`imputeMissing`](https://igordot.github.io/copynumber/reference/imputeMissing.md),
[`pcf`](https://igordot.github.io/copynumber/reference/pcf.md)

## Author

Gro Nilsen, Knut Liestoel

## Examples

``` r

#Load lymphoma data:
data(lymphoma)

#Take out a subset of 3 biopsies from the first patient (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:3)

#Check for missing values in data:
any(is.na(sub.lymphoma))
#> [1] FALSE
#FALSE

#First winsorize data to handle outliers:
wins.lymph <- winsorize(sub.lymphoma)
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

#Run multipcf on subset lymphoma data (using a low gamma because of low-density data)
multi.segments <- multipcf(data=wins.lymph,gamma=12,Y=sub.lymphoma)
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

```
