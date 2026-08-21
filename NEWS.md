# copynumber 2.0.0

## Bug fixes

* `addChromlines()`, `get.seglim()`, and `plotCircle()` had minor bugs fixed; several error messages had incomplete or misspelled arguments.
* `aspcf(..., return.est = TRUE)` now returns named logR estimates.
* `aspcf()`, `multipcf()`, `pcf()`, `plotAberration()`, `plotHeatmap()`, `subsetData()`, `subsetSegments()`, and `winsorize()` no longer error on matrix input.
* `interpolate.pcf()` no longer errors on a missing chromosome or single-row input.
* `selectSegments()` now returns a list, not a bare data frame, for large `nseg`.

## Development

* Converted function documentation to roxygen2.
* Converted the vignette from Sweave to R Markdown.
* Added a pkgdown website.
* Added a testthat test suite (previously no automated tests).

# copynumber 1.29.0.9000

* Added `hg38` and `mm10` genome build support.
* Set up continuous integration via GitHub Actions.

# copynumber 1.2.0 to 1.41.0

Version-only bumps per Bioconductor release, except where noted.

* 1.41.0 — devel bump after Bioconductor 3.17 (2023-04-25); no further releases followed
* 1.40.0 — Bioconductor 3.17 (2023-04-25); final Bioconductor release
* 1.39.1 — marked `Deprecated` (unresponsive maintainer, invalid contact email)
* 1.38.0 — Bioconductor 3.16 (2022-11-01)
* 1.36.0 — Bioconductor 3.15 (2022-04-26)
* 1.34.0 — Bioconductor 3.14 (2021-10-26)
* 1.32.0 — Bioconductor 3.13 (2021-05-19)
* 1.30.0 — Bioconductor 3.12 (2020-10-27)
* 1.29.0 — devel bump after Bioconductor 3.11; parent of this fork's 1.29.0.9000 (2020-04-27)
* 1.28.0 — Bioconductor 3.11 (2020-04-27)
* 1.26.0 — Bioconductor 3.10 (2019-10-29)
* 1.24.0 — Bioconductor 3.9 (2019-05-02)
* 1.22.0 — Bioconductor 3.8 (2018-10-30)
* 1.20.0 — Bioconductor 3.7 (2018-04-30)
* 1.18.0 — Bioconductor 3.6 (2017-10-30)
* 1.16.0 — Bioconductor 3.5 (2017-04-24)
* 1.14.0 — Bioconductor 3.4 (2016-10-17)
* 1.12.0 — Bioconductor 3.3 (2016-05-03)
* 1.10.0 — Bioconductor 3.2 (2015-10-13)
* 1.8.0 — Bioconductor 3.1 (2015-04-16)
* 1.6.0 — Bioconductor 3.0 (2014-10-13)
* 1.4.0 — Bioconductor 2.14 (2014-04-11)
* 1.2.0 — Bioconductor 2.13 (2013-10-14)

# copynumber 1.1.1

* `plotHeatmap()` now uses different colors above and below zero, instead of assuming `upper.lim`/`lower.lim` are symmetrical.

# copynumber 1.0.0

* Initial version bump for Bioconductor 2.12 (2013-04-03).

# copynumber 0.99.1

* Package first submitted to Bioconductor (2013-03-29).
