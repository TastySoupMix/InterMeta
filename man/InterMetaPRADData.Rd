\name{InterMetaPRADData}
\alias{InterMetaPRADData}
\alias{InterMetaPRADInfo}
\alias{InterMetaPRADRNA}
\alias{InterMetaPRADEnrich}
\alias{InterMetaPRADScore}
\alias{InterMetaPRADProvenance}
\docType{data}
\title{PCaDB TCGA-PRAD Data for the InterMeta Workflow}
\description{
A PCaDB TCGA-PRAD biochemical-recurrence cohort, a metabolic-gene expression
subset, and converted historical PRAD reference results for demonstrating the
InterMeta workflow.
}
\usage{
data("InterMetaPRADInfo")
data("InterMetaPRADRNA")
data("InterMetaPRADEnrich")
data("InterMetaPRADScore")
data("InterMetaPRADProvenance")
}
\format{
Five objects are provided:
\describe{
  \item{\code{InterMetaPRADInfo}}{A data frame with 463 TCGA-PRAD samples as
  row names. The columns are event status and positive follow-up time in days.
  There are 99 observed events and 364 censored observations.}
  \item{\code{InterMetaPRADRNA}}{A numeric matrix of normalized bulk
  expression for 1,260 metabolic genes (rows) across the same 463 samples
  (columns).}
  \item{\code{InterMetaPRADEnrich}}{A data frame containing all 2,652 directed
  interactions among 52 pathways, BH-adjusted enrichment P-values, and FDR
  values recomputed with 100 \code{TimiGP-legacy} annotation permutations.}
  \item{\code{InterMetaPRADScore}}{A 50-row TimiGP-legacy pathway score
  table calculated using \code{Adjust.P.Value < 0.05}, with signed ranks and
  neutral pathways removed. Re-scored in 1.1.4 from the existing enrichment.}
  \item{\code{InterMetaPRADProvenance}}{A list containing PCaDB source links and
  DOI, the TCGA-PRAD cohort and day time unit, cohort dimensions, analysis and
  FDR-lineage metadata, fresh no-checkpoint computational-reproduction
  evidence, and SHA-256 fingerprints of the source inputs and historical
  result files.}
}
}
\details{
The original expression table contained 55,141 genes and 495 samples. The
package object retains the 463 samples with prognosis data and the 1,260 genes
that matched the 1,292 unique genes in the bundled \code{geneset-all}
annotation. To reproduce the bundled analysis, use the sole
\code{"TimiGP-legacy"} mode with \code{log = TRUE}. It retains genes
positive in at least half of samples, safely removes rows containing values
less than or equal to -1, applies \code{log10(x + 1)}, and median-centers each
gene. It retains 861 genes for the full cohort and yields 370,230 gene pairs.

The bundled enrichment result migrates the historical TimiGP PRAD analysis.
Its gene-pair query used unadjusted Cox \code{PV < 0.05}, BH adjustment for
pathway enrichment, and 100 annotation permutations. The public package-data
generator recomputes enrichment counts, P-values, and permutation FDR; the
supplied historical enrichment object is used only to check numerical
agreement. A fresh no-checkpoint full workflow, identified in provenance as
\code{prad-complete-full-1.1.3/run_20260917_101835_fresh/
InterMeta_TCGA_PRAD_complete_1_1_3}, then recomputed the Cox models,
enrichment, and 100-permutation FDR from the bundled PCaDB-derived expression
and outcome inputs with seed 1. It produced 370,230 Cox pairs, selected 101,323
pairs at \code{PV < 0.05}, and produced 2,652 enrichment rows and 50
non-neutral legacy score rows. All 12 parity checks and all three prespecified
80-percent internal resamples passed. The freshly recomputed FDR values agree
with the previous bundled values to a maximum absolute difference of about
\code{2.22e-16}; the bundled column was republished from the fresh run rather
than copied from the historical reference. The three resamples are a
descriptive internal-consistency check, not external validation. The legacy profile also
uses the original rounded gene-pair prevalence threshold, factor-based Cox
predictor for logical pairs, signed \code{Difference} rank, and removal of
neutral pathways. Its 50-pathway score table agrees with the original pathway
identity, order, ranks, groups, and numeric scores to floating-point error of
about \code{5.33e-15}. This preserves the historical analysis but is
exploratory. New confirmatory analyses should use \code{QV < 0.05}, or clearly
report both thresholds as a sensitivity analysis.

The 0.4.0 server result is not an equivalent legacy reproduction: it retained
1,260 genes, selected \code{QV}, sampled events and censored observations with
seeds 1--3, and used absolute score ranks while retaining neutral pathways.
Those choices affect the Cox universe and sampled patients, so exact historical
equivalence cannot be recovered by post-processing that output; recompute from
these inputs with \code{analysis_profile = "TimiGP-legacy"}.

The compact InterMeta implementation reproduced every legacy enrichment
count and P-value. Large gene-pair strings were omitted from the migrated
object to keep the installed package small. Consequently,
\code{InterMetaPRADEnrich} can be used for the pathway chord but not directly
for \code{InterMetaGeneChord()}; rerun enrichment with
\code{keep_gene_pairs = "shared"} to retain the gene pairs needed by that
plot.
}
\source{
The prognosis table and normalized bulk transcriptomic table used to construct
these TCGA-PRAD objects were downloaded from PCaDB. Its historical website
endpoint is recorded in \code{InterMetaPRADProvenance}; the stable source
repository is \url{https://github.com/rli012/PCaDB}. Cite the PCaDB preprint using
\doi{10.1101/2021.06.29.449134}. The package retains the aligned prognosis
cohort and metabolic-gene expression subset used by the examples, and records
source-input fingerprints in \code{InterMetaPRADProvenance}. The original
PCaDB access date, release version, and immutable snapshot were not present in
the available source files and are recorded as unknown rather than inferred.
The accompanying
manuscript should identify PCaDB as the download source and cite it explicitly.
}
\keyword{datasets}
