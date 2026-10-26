# Mouse mm10 data

Bin log2 copy ratios from mouse tumor cell lines WT4 (GSM6720269) and WT6
(GSM6720270), from [GSE217479](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217479).
Reads were aligned directly to mm10. CNVkit used a flat exome reference.

| File | Bins | Selection |
|---|---:|---|
| `mouse-mm10-GSE217479-sparse.tsv` | 2,032 | 100 bins per chromosome (1–19 and X), plus 32 Y bins. |
| `mouse-mm10-GSE217479-dense.tsv` | 13,909 | All usable shared bins on chromosomes 18 and 19. |

Columns: `chrom`, `pos` (one-based bin midpoint), `WT4`, and `WT6`.

Run `Rscript data-raw/test-datasets.R` from the package root to make these fixtures.
The script downloads the CNR files, joins samples by bin coordinates, and removes
bins with non-finite ratios, zero depth, or zero weight in either sample.
It keeps the log2 ratios unchanged and sets `pos = floor((start + end) / 2) + 1`.
The sparse fixture selects bins at equal intervals in bin order.

Tests cover `winsorize()`, `pcf()`, `multipcf()`, `plotGenome()`, and `plotCircle()` with mm10.
The fixtures have exome gaps, no BAF track, no CNVkit weights, and no known true
breakpoints. Do not use them to measure segmentation accuracy.

Tests check output dimensions and selected values from the first mm10 commit,
`ef99c55fb69a94cdf4a672b93bb0ca079a8f1b5b` (2020-07-15). The values were
checked with historical code and cytobands under the current R environment.
They are a regression baseline, not known biological truth.
