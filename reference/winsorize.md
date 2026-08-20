# Winsorization of copy number data

Outliers in copy number data are detected and modified using MAD or PCF
Winsorization.

## Usage

``` r
winsorize(
  data,
  pos.unit = "bp",
  arms = NULL,
  method = "mad",
  tau = 2.5,
  k = 25,
  gamma = 40,
  iter = 1,
  assembly = "hg19",
  digits = 4,
  return.outliers = FALSE,
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

- method:

  the Winsorization method to be applied, must be one of "mad" (default)
  or "pcf".

- tau:

  Winsorization threshold, default is 2.5.

- k:

  the half window size to be applied in median filtering, default is 25.

- gamma:

  penalty for each discontinuity in the pcf curve, default is 40. Only
  applicable when `method="pcf"`.

- iter:

  number of iterations in PCF Winsorization, default is 1.

- assembly:

  a string specifying which genome assembly version should be applied to
  determine chromosome arms. Allowed options are "hg19", "hg18", "hg17"
  and "hg16" (corresponding to the four latest human genome annotations
  in the UCSC genome browser).

- digits:

  the number of decimals to be applied when reporting results. Default
  is 4.

- return.outliers:

  logical value indicating whether a data frame identifying outliers
  should be returned, default is FALSE.

- save.res:

  logical value indicating whether results should be saved in text
  files, default is FALSE.

- file.names:

  optional character vector of length two giving the name of the files
  where the Winsorized data and outlier statuses, respectively, should
  be saved if `save.res=TRUE`.

- verbose:

  logical value indicating whether or not to print a progress message
  each time Winsorization is finished for a new chromosome arm.

## Value

If `return.outliers = TRUE` a list with the following components:

- wins.data:

  a data frame with chromosome numbers in the first column, probe
  positions in the second and the Winsorized copy number values for the
  sample(s) in subsequent column(s).

- wins.outliers:

  a data frame with chromosome numbers in the first column, probe
  positions in the second and outlier statuses for each sample in the
  subsequent column(s). The values +/- 1 indicate that the observation
  is an outlier, whereas the value 0 indicates that it is not.

If `return.outliers = FALSE` only the data frame containing the
winsorized data is returned.

If `save.res=TRUE` the results are saved in text files with names as
specified in `file.names`. If `file.names=NULL`, a folder named
"Wins_res" is created in the working directory and Winsorized data and
outlier statuses are saved in this directory in tab-separated files
named wins.data.txt and wins.outliers.txt, respectively.

## Details

The copy number data are either MAD Winsorized or PCF Winsorized as
described in Nilsen and Liestoel et al. (2012). Winsorization is done
separately on each chromosome arm in each sample.

## Note

Any missing values in `data` imply that the Winsorized value and outlier
status for this probe will be missing as well. Also, if the number of
probes within a chromosome arm is less than 2\*k, Winsorization cannot
be done and the data values are thus left unchanged.

## References

Nilsen and Liestoel et al., "Copynumber: Efficient algorithms for
single- and multi-track copy number segmentation", BMC Genomics 13:591
(2012), doi:10.1186/1471-2164-13-59

## Author

Gro Nilsen, Knut Liestoel, Ole Christian Lingjaerde

## Examples

``` r

#Lymphoma data
data(lymphoma)
#Take out a smaller subset of 3 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:3)

#Do MAD Winsorization:
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
         
```
