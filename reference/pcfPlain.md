# Plain single-sample copy number segmentation.

A basic single-sample pcf segmentation which does not take chromosome
borders into account

## Usage

``` r
pcfPlain(
  pos.data,
  kmin = 5,
  gamma = 40,
  normalize = TRUE,
  fast = TRUE,
  digits = 4,
  return.est = FALSE,
  verbose = TRUE
)
```

## Arguments

- pos.data:

  a data frame where the rows represent the probes, column 1 holds probe
  positions, and subsequent column(s) give the numeric copy number
  measurements for one or more samples. The header of copy number
  columns should give sample IDs.

- kmin:

  minimum number of probes in each segment, default is 5.

- gamma:

  penalty for each discontinuity in the curve, default is 40.

- normalize:

  logical value indicating whether the copy number measurements should
  be scaled by the sample residual standard error. Default is TRUE.

- fast:

  a logical value indicating whether a fast (not guaranteed to be exact)
  version should be run if the number of probes are \> 400.

- digits:

  the number of decimals to be applied when reporting results. Default
  is 4.

- return.est:

  logical value indicating whether a data frame holding copy number
  estimates (pcf values) should be returned along with the segments.
  Default is FALSE, which means that only segments are returned.

- verbose:

  logical value indicating whether or not to print a progress message
  each time pcf analysis is finished for a sample.

## Value

If `return.est = TRUE` a list with the following components:

- estimates:

  a data frame where the first column gives the probe positions, while
  subsequent column(s) give the copy number estimates for each sample.
  The estimate for a given probe equals the mean of the segment where
  the probe is located.

- segments:

  a data frame describing each segment found in the data. Each row
  represents a segment, while columns give the sampleID, start position,
  end position, number of probes in the segment and mean value,
  respectively.

If `return.est = FALSE`, only the data frame containing the segments is
returned.

## Details

A piecewise constant segmentation curve is fitted to the copy number
observations as described in the PCF algorithm in Nilsen and Liestoel et
al. (2012). Unlike the regular
[`pcf`](https://igordot.github.io/copynumber/reference/pcf.md) function,
`pcfPlain` does not make independent segmentations for each chromosome
arm (i.e. breakpoints are not automatically inserted at the beginning
and end of chromosome arms). The segmentation can thus be performed
independently of assembly.

## Note

If probe positions are not available, the first column in `data` may,
e.g., contain the values `1:nrow(data)`.

## References

Nilsen and Liestoel et al., "Copynumber: Efficient algorithms for
single- and multi-track copy number segmentation", BMC Genomics 13:591
(2012), doi:10.1186/1471-2164-13-59

## See also

[`pcf`](https://igordot.github.io/copynumber/reference/pcf.md)

## Author

Gro Nilsen, Knut Liestoel, Ole Christian Lingjaerde.

## Examples

``` r

#Load the lymphoma data set:
data(lymphoma)

#Take out a smaller subset of 3 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:3)

#Run pcfPlain (remove first column of chromosome numbers):
plain.segments <- pcfPlain(pos.data=sub.lymphoma[,-1],gamma=12)
#> pcf finished for sample 3 


```
