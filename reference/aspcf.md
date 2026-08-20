# Allele-specific copy number segmentation.

Joint segmentation of SNP array data resulting in piecewise constant
curves with common break points for copy number data and B-allelle
frequency data.

## Usage

``` r
aspcf(
  logR,
  BAF,
  pos.unit = "bp",
  arms = NULL,
  kmin = 5,
  gamma = 40,
  baf.thres = c(0.1, 0.9),
  skew = 3,
  assembly = "hg19",
  digits = 4,
  return.est = FALSE,
  save.res = FALSE,
  file.names = NULL,
  verbose = TRUE
)
```

## Arguments

- logR:

  either a data frame or the name of a tab-separated file from which
  copy number data can be read. The rows of the data frame or file
  should represent the probes. Column 1 must hold numeric or character
  chromosome numbers, column 2 the numeric local probe positions, and
  subsequent columns the numeric copy number measurements for one or
  more samples. The header of copy number column(s) should give sample
  ID(s).

- BAF:

  either a data frame or the name of a tab-separated file from which
  B-allelle frequency data can be read. Must be on the same format and
  size as `logR`, with chromosomes and local probe positions in the two
  first columns, and numeric BAF-measurements for one or more samples in
  subsequent columns.

- pos.unit:

  the unit used to represent the probe positions. Allowed options are
  "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs).
  By default assumed to be "bp".

- arms:

  optional character vector containing chromosome arms (denoted 'p' and
  'q') corresponding to the chromosomes and positions found in `logR`
  and `BAF`. If not specified chromosome arms are found using the
  built-in genome assembly version determined by `assembly`.

- kmin:

  minimum number of probes in each segment, default is 5.

- gamma:

  penalty for each discontinuity in the curve, default is 40.

- baf.thres:

  a numeric vector of length two giving the thresholds below and above
  which BAF probes are considered germline homozygous. Must be in the
  range 0 to 1, default is 0.1 and 0.9 for the lower and upper limit,
  respectively.

- skew:

  a numeric value used to determine whether there is allelic skewness
  (one or two bands) in BAF. Default is 3. The larger the value the
  further the BAF measurements must be from 0.5 to imply two bands.

- assembly:

  a string specifying which genome assembly version should be applied to
  determine chromosome arms. Allowed options are "hg19", "hg18", "hg17"
  and "hg16" (corresponding to the four latest human genome annotations
  in the UCSC genome browser).

- digits:

  the number of decimals to be applied when reporting results. Default
  is 4.

- return.est:

  logical value indicating whether a data frame holding LogR estimates
  should be returned along with the segments. Default is FALSE, which
  means that only segments are returned.

- save.res:

  logical value indicating whether results should be saved in text
  files, default is FALSE.

- file.names:

  optional character vector of length two giving the name of the files
  where the logR estimates and segments, respectively, should be saved
  in case `save.res=TRUE`.

- verbose:

  logical value indicating whether or not to print a progress message
  each time aspcf analysis is finished for a new chromosome arm.

## Value

If `return.est = TRUE` a list with the following components:

- logR_estimates:

  a data frame where the first two columns give the chromosome numbers
  and probe positions, respectively, while subsequent column(s) give the
  LogR estimates for each sample. The estimate for a given probe equals
  the mean of the segment where the probe is located.

- segments:

  a data frame describing each segment found. Each row represents a
  segment, and columns give the sample IDs, chromosome numbers, arms,
  local start positions, local end positions, number of probes in the
  segments, mean LogR values and mean BAF values, respectively.

If `return.est = FALSE`, only the data frame containing the segments is
returned.

If `save.res = TRUE` the results are also saved in text files with names
as specified in `file.names`. If `file.names=NULL`, a folder named
"aspcf_results" is created in the working directory, and the LogR
estimates and the segmentation results are saved in this folder as
tab-separated files named logR_estimates.txt and segments.txt,
respectively.

## Details

Piecewise constant curves are simultaneously fitted to the LogR and BAF
data as described in Nilsen and Liestoel et al.(2012). This implies that
break points will be the same for the LogR and BAF segmentation curves,
while segment values differ. Segmentation is done separately on each
chromosome arm in each sample.

## Note

It will usually be advisable to Winsorize the logR data before running
`aspcf`, see
[`winsorize`](https://igordot.github.io/copynumber/reference/winsorize.md)
on this. Missing values are not allowed in `logR`, see `imputeMissing`
for imputation of missing copy number values.

## References

Nilsen and Liestoel et al., "Copynumber: Efficient algorithms for
single- and multi-track copy number segmentation", BMC Genomics 13:591
(2012), doi:10.1186/1471-2164-13-59

## See also

[`plotAllele`](https://igordot.github.io/copynumber/reference/plotAllele.md),
[`winsorize`](https://igordot.github.io/copynumber/reference/winsorize.md)

## Author

Gro Nilsen, Knut Liestoel, Ole Christian Lingjaerde

## Examples

``` r

#Load LogR and BAF data:
data(logR)
data(BAF)

#First winsorize logR to handle outliers:
wins.logR <- winsorize(logR)
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
#> winsorize finished for chromosome arm Xp 
#> winsorize finished for chromosome arm Xq 

#Run aspcf:
aspcf.segments <- aspcf(wins.logR,BAF)
#> aspcf finished for chromosome arm 1p 
#> aspcf finished for chromosome arm 1q 
#> aspcf finished for chromosome arm 2p 
#> aspcf finished for chromosome arm 2q 
#> aspcf finished for chromosome arm 3p 
#> aspcf finished for chromosome arm 3q 
#> aspcf finished for chromosome arm 4p 
#> aspcf finished for chromosome arm 4q 
#> aspcf finished for chromosome arm 5p 
#> aspcf finished for chromosome arm 5q 
#> aspcf finished for chromosome arm 6p 
#> aspcf finished for chromosome arm 6q 
#> aspcf finished for chromosome arm 7p 
#> aspcf finished for chromosome arm 7q 
#> aspcf finished for chromosome arm 8p 
#> aspcf finished for chromosome arm 8q 
#> aspcf finished for chromosome arm 9p 
#> aspcf finished for chromosome arm 9q 
#> aspcf finished for chromosome arm 10p 
#> aspcf finished for chromosome arm 10q 
#> aspcf finished for chromosome arm 11p 
#> aspcf finished for chromosome arm 11q 
#> aspcf finished for chromosome arm 12p 
#> aspcf finished for chromosome arm 12q 
#> aspcf finished for chromosome arm 13q 
#> aspcf finished for chromosome arm 14q 
#> aspcf finished for chromosome arm 15q 
#> aspcf finished for chromosome arm 16p 
#> aspcf finished for chromosome arm 16q 
#> aspcf finished for chromosome arm 17p 
#> aspcf finished for chromosome arm 17q 
#> aspcf finished for chromosome arm 18p 
#> aspcf finished for chromosome arm 18q 
#> aspcf finished for chromosome arm 19p 
#> aspcf finished for chromosome arm 19q 
#> aspcf finished for chromosome arm 20p 
#> aspcf finished for chromosome arm 20q 
#> aspcf finished for chromosome arm 21q 
#> aspcf finished for chromosome arm 22q 
#> aspcf finished for chromosome arm Xp 
#> aspcf finished for chromosome arm Xq 

```
