# copynumber with hg38 and mm10

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/igordot/copynumber.svg?branch=master)](https://travis-ci.org/igordot/copynumber)
![R-CMD-check](https://github.com/igordot/copynumber/workflows/R-CMD-check/badge.svg)
![R-CMD-check-bioc](https://github.com/igordot/copynumber/workflows/R-CMD-check-bioc/badge.svg)
<!-- badges: end -->

This is an updated version of the [copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html) R package. It has been modified to support the human hg38 and mouse mm10 genome builds (`assembly` parameter in `aspcf`, `multipcf`, `pcf`, and `winsorize` functions).

You can install the package from GitHub with:

```r
BiocManager::install("igordot/copynumber")
```

This modification builds upon the previous work of [aroneklund](https://github.com/aroneklund/copynumber) and [ShixiangWang](https://github.com/ShixiangWang/copynumber). There is also a fork by [Irrationone](https://github.com/Irrationone/copynumber) that provides a species-agnostic approach, but the enhanced functionality adds a new parameter that makes it not backward compatible. If you are using packages that depend on `copynumber`, such as `sequenza` and `PureCN`, they would need to be modified as well.

This fork is based on copynumber 1.29.0 (part of Bioconductor 3.12). However, versioning may not be particularly relevant as there have been no changes to the code since 2013.
