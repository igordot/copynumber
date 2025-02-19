# copynumber with hg38 and mm10

<!-- badges: start -->
[![R-CMD-check](https://github.com/igordot/copynumber/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/igordot/copynumber/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This is an updated version of the [copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html) R package.
It has been modified to support the human hg38 and mouse mm10 genome builds.
Specifically, this is used for the `assembly` parameter in `aspcf`, `multipcf`, `pcf`, and `winsorize` functions.

You can install this version of the package from GitHub using:

```r
BiocManager::install("igordot/copynumber")
```

This modification builds upon the previous work of [aroneklund](https://github.com/aroneklund/copynumber) and [ShixiangWang](https://github.com/ShixiangWang/copynumber).
There is an additional fork by [Irrationone](https://github.com/Irrationone/copynumber) that provides a species-agnostic approach, but the enhanced functionality adds a new parameter that makes it not backward compatible.
Any packages that depend on `copynumber` (such as `sequenza`, `PureCN`, or `scarHRD`) would need to be modified as well to take advantage of the new feature.

This fork is based on copynumber 1.29.0 (part of Bioconductor 3.12).
The package was removed from Bioconductor in version 3.18.
The specific version is likely not relevant for most applications as there have been no changes to the code since 2013.
