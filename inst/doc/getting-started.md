# Getting started with InterMeta

This guide gives a short, offline installation check using the synthetic data
bundled with InterMeta. The synthetic cohort is for learning and software
verification only; it is not biological or clinical evidence.

## TCGA-PRAD and the synthetic check are different examples

TCGA-PRAD is the principal real-cohort example. For a complete recomputation
use `example09_tcga_prad_complete_workflow.R` and its documented configuration.
It selects gene pairs using `PV < 0.05` and pathway interactions for scoring
using `Adjust.P.Value < 0.05`. Permutation FDR is calculated and saved separately.

The synthetic example below uses 120 generated samples and 30 genes. In 1.1.4
its minimum permutation FDR is 0.09. Its `Permutation.FDR < 0.05` filter thus
returns an empty score/network and no selected pathway chord links. Expected
warnings explain this result; it does not mean installation failed. Do not
relax an inferential threshold solely to force a nonempty figure.

For a quick look at real-data plotting objects, run
`example05_tcga_prad_results.R`. It loads precomputed PRAD results rather than
refitting models. Its interaction display uses `Permutation.FDR < 0.2`, while
the saved score was built with `Adjust.P.Value < 0.05`; these are different
selection rules. Use the full workflow for one consistent analysis contract.

## 1. Inspect the bundled inputs

```r
library(InterMeta)

data("InterMetaExampleInfo")
data("InterMetaExampleRNA")
data("InterMetaExamplePathwayset")

dim(InterMetaExampleInfo)
dim(InterMetaExampleRNA)
table(InterMetaExampleInfo$event)
table(InterMetaExamplePathwayset$Pathway)
```

The same inputs are available as plain-text files:

```r
system.file("extdata", "example", package = "InterMeta")
```

Run the consistency check that compares the package objects with those text
files:

```r
source(system.file(
  "examples", "example00_bundled_data.R", package = "InterMeta"
))
```

## 2. Run the small end-to-end workflow

The complete synthetic example preprocesses expression values, constructs
gene pairs, fits Cox models, performs pathway-pair enrichment and permutation
FDR estimation, calculates pathway scores, and constructs plotting objects.
It runs without network access or credentials.

```r
source(system.file(
  "examples",
  "example01_metabolic_pathway_prognosis.R",
  package = "InterMeta"
))
```

The principal objects created by the example are `cox_res`, `enrich`, `score`,
and `net`. Consult their function help before adapting thresholds or profiles:

```r
?InterMetaCox
?InterMetaEnrich
?InterMetaPathwayScore
```

## 3. Choose the appropriate next guide

- Use `prad-reproduction.md` for the bundled PCaDB TCGA-PRAD inputs and the
  compatibility analysis.
- Use `ai-workflow.md` only after the statistical result has been inspected;
  the AI layer is optional and does not replace statistical or biological
  validation.

From an installed package, locate all guides and scripts with:

```r
list.files(system.file("doc", package = "InterMeta"), full.names = TRUE)
list.files(system.file("examples", package = "InterMeta"), full.names = TRUE)
```
