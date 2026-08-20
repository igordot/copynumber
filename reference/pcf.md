# Single-sample copy number segmentation.

Fit a individual piecewise constant segmentation curve to each sample's
copy number data.

## Usage

``` r
pcf(
  data,
  pos.unit = "bp",
  arms = NULL,
  Y = NULL,
  kmin = 5,
  gamma = 40,
  normalize = TRUE,
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
  subsequent column(s) the numeric copy number measurements for one or
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

- kmin:

  minimum number of probes in each segment, default is 5.

- gamma:

  penalty for each discontinuity in the curve, default is 40.

- normalize:

  logical value indicating whether the copy number measurements should
  be scaled by the sample residual standard error. Default is TRUE.

- fast:

  a logical value indicating whether a fast (not guaranteed to be exact)
  version should be run on chromosome arms with \> 400 probes.

- assembly:

  a string specifying which genome assembly version should be applied to
  determine chromosome arms. Allowed options are "hg19", "hg18", "hg17"
  and "hg16" (corresponding to the four latest human genome annotations
  in the UCSC genome browser).

- digits:

  the number of decimals to be applied when reporting results. Default
  is 4.

- return.est:

  logical value indicating whether a data frame holding copy number
  estimates (pcf values) should be returned along with the segments.
  Default is FALSE, which means that only segments are returned.

- save.res:

  logical value indicating whether results should be saved in text
  files.

- file.names:

  optional character vector of length two giving the name of the files
  where the pcf estimates and segments, respectively, should be saved in
  case `save.res=TRUE`.

- verbose:

  logical value indicating whether or not to print a progress message
  each time pcf analysis is finished for a new chromosome arm.

## Value

If `return.est = TRUE` a list with the following components:

- estimates:

  a data frame where the first two columns give the chromosome numbers
  and probe positions respectively, while subsequent column(s) give the
  copy number estimates for each sample. The estimate for a given probe
  equals the mean of the segment where the probe is located.

- segments:

  a data frame describing each segment found in the data. Each row
  represents a segment, while columns give the sampleID, chromosome
  number, arm, local start position, local end position, number of
  probes in the segment and mean value, respectively.

If `return.est = FALSE`, only the data frame containing the segments is
returned.

If `save.res = TRUE` the results are also saved in text files with names
as specified in `file.names`. If `file.names=NULL`, a folder named
"pcf_results" is created in the working directory, and the pcf estimates
and segments are saved in this directory in tab-separated files named
estimates.txt and segments.txt, respectively.

## Details

A piecewise constant segmentation curve is fitted to the copy number
observations as described in the PCF algorithm in Nilsen and Liestoel et
al. (2012). Segmentation is done separately on each chromosome arm in
each sample.

## Note

It is usually advisable to Winsorize data before running pcf, see
[`winsorize`](https://igordot.github.io/copynumber/reference/winsorize.md)
on this.

Missing copy number values are allowed. These are kept out of the pcf
analysis, and copy number estimates for missing observations are later
set to be the same as the estimate of the nearest observed probe.

## References

Nilsen and Liestoel et al., "Copynumber: Efficient algorithms for
single- and multi-track copy number segmentation", BMC Genomics 13:591
(2012), doi:10.1186/1471-2164-13-59

## See also

[`multipcf`](https://igordot.github.io/copynumber/reference/multipcf.md)

## Author

Gro Nilsen, Knut Liestoel, Ole Christian Lingjaerde.

## Examples

``` r

#Load the lymphoma data set:
data(lymphoma)

#Take out a smaller subset of 3 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:3)

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

#Run pcf (using small gamma because of low-density data):
pcf.segments <- pcf(data=wins.lymph,gamma=12,Y=sub.lymphoma)
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
