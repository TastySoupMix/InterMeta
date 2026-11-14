# Reproducing the bundled TCGA-PRAD analysis

## Data source and scope

The prognosis and normalized bulk transcriptomic input tables used for the
bundled TCGA-PRAD objects were downloaded from PCaDB:

- PCaDB historical website endpoint: recorded in `InterMetaPRADProvenance`
- PCaDB source repository: <https://github.com/rli012/PCaDB>
- PCaDB preprint: <https://doi.org/10.1101/2021.06.29.449134>

The package retains the 463-sample aligned prognosis cohort and the 1,260-gene
metabolic expression subset used by its examples. The manuscript and any
derived analysis should identify PCaDB as the download source and cite it.
Machine-readable source, transformation, analysis-lineage, and fingerprint
fields are stored in `InterMetaPRADProvenance`. The original download/access
date, database release version, and immutable snapshot were not present in the
available files, so they are explicitly recorded as unknown.

```r
library(InterMeta)
data("InterMetaPRADProvenance")
str(InterMetaPRADProvenance, max.level = 2)
```

The package-data generation and provenance-update scripts are shipped as
reproducibility material under `system.file("scripts", package = "InterMeta")`.
They require the source inputs described in their headers and are not run when
the package is installed.

## Inspect the bundled reference result

This path is fast and does not refit Cox models or rerun permutations:

```r
source(system.file(
  "examples", "example05_tcga_prad_results.R", package = "InterMeta"
))
```

The bundled reference preserves the historical analysis contract. Its
`Permutation.FDR` column was republished from a fresh InterMeta 1.1.3 run that
refitted Cox models and reran all 100 permutations with checkpoint reuse
disabled. The historical enrichment object is now used only for numerical
parity checks. Loading the bundled object remains a data-loading operation and
must not itself be described as running permutations.

The accepted run passed 12 full-sample parity checks and all three
prespecified 80% internal resamples. Its identifier and package, workflow, and
output-manifest SHA-256 values are stored under
`InterMetaPRADProvenance$permutation_fdr$fresh_reproduction`. Those three
resamples describe internal consistency; they are not an external validation
cohort.

## Recompute from the bundled PRAD inputs

The following script starts with `InterMetaPRADInfo` and `InterMetaPRADRNA`,
uses the explicit `TimiGP-legacy` compatibility profile, refits the gene-pair
Cox models, reruns pathway enrichment, and by default reruns 100 annotation
permutations:

```r
script <- system.file(
  "examples", "example06_tcga_prad_recompute.R", package = "InterMeta"
)
source(script)
```

The output directory is `InterMeta_PRAD_recomputed` below the working
directory. The full run is computationally intensive. The script contains an
explicit `run_permutation` switch for a faster structural check, but a run with
that switch disabled must not be reported as having newly estimated
`Permutation.FDR` values.

Before a manuscript run, record the installed package version and execution
environment alongside the result:

```r
packageVersion("InterMeta")
sessionInfo()
```

Also archive the exact script, parameter values, random seeds, input
fingerprints, and output fingerprints used for the paper. The complete
publication workflow is available separately:

```r
system.file(
  "examples", "example09_tcga_prad_complete_workflow.R", package = "InterMeta"
)
```

That script adds internal resampling and the optional evidence workflow. Read
its configuration block before running it; PubMed and model-service actions
have different network, privacy, and cost implications from the statistical
analysis.

For release-provenance publication, the source tree also contains the
fail-closed `data-raw/update_prad_provenance.R` wrapper and its public
`inst/scripts/update_prad_provenance_impl.R` implementation. The updater
accepts only a complete 1.1.3 output whose manifest hashes, 12 parity checks,
three successful resamples, current-input fingerprint, candidate archive
workflow-script hash, and recomputed FDR all verify. Because the original run
did not record an archive or core-code hash, the candidate archive is retained
as supporting evidence and is not claimed as a cryptographic identity of the
executing installation.

## Interpretation boundary

TimiGP-legacy is the sole analysis mode from 1.1.4 and preserves the historical analysis.
It is not automatically a confirmatory analysis. In particular, report the
gene-pair selection field (`PV` or `QV`), cutoff, pathway-level adjustment,
permutation count, score condition, and analysis profile explicitly. Keep
external-cohort validation distinct from internal resampling.
