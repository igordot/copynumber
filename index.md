# copynumber with hg38 and mm10

This is an updated version of the
[copynumber](https://bioconductor.org/packages/3.17/bioc/html/copynumber.html)
R package. It has been modified to support the human hg38 and mouse mm10
genome builds.

## Installation

You can install the package from
[R-universe](https://igordot.r-universe.dev/copynumber) (pre-compiled
binary):

``` r

BiocManager::install("copynumber", site_repository = "https://igordot.r-universe.dev")
```

Alternatively, you can install from source (slower):

``` r

BiocManager::install("igordot/copynumber")
```

## About

This fork of copynumber adds the following improvements to the original
Bioconductor package:

- Support for more recent genome builds, used by the `assembly`
  parameter in
  [`aspcf()`](https://igordot.github.io/copynumber/reference/aspcf.md),
  [`multipcf()`](https://igordot.github.io/copynumber/reference/multipcf.md),
  [`pcf()`](https://igordot.github.io/copynumber/reference/pcf.md), and
  [`winsorize()`](https://igordot.github.io/copynumber/reference/winsorize.md)
- A [documentation website](https://igordot.github.io/copynumber/) (via
  roxygen2 and pkgdown)
- A clean, browsable [HTML
  vignette](https://igordot.github.io/copynumber/articles/copynumber.html)
  (via knitr and R Markdown)
- Automated tests to help ensure that the package continues to function

This modification builds upon the previous work of
[aroneklund](https://github.com/aroneklund/copynumber) and
[ShixiangWang](https://github.com/ShixiangWang/copynumber). There is an
additional fork by
[Irrationone](https://github.com/Irrationone/copynumber) that provides a
species-agnostic approach, but the enhanced functionality adds a new
parameter that makes it not backward compatible. Any packages that
depend on `copynumber` (such as `sequenza`, `PureCN`, or `scarHRD`)
would need to be modified as well to take advantage of the new feature.

This fork is based on copynumber 1.29.0. The package was removed in
Bioconductor 3.18 (released October 2023). There have been no changes to
the code besides Bioconductor-mandated version bumps since 2013.

Please use [GitHub issues](https://github.com/igordot/copynumber/issues)
to report any problems or request new features.
