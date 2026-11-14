\name{InterMeta}
\alias{InterMeta}
\alias{InterMetaAvailableGenesets}
\alias{InterMetaLoadGeneset}
\alias{InterMetaCheckEvent}
\alias{InterMetaPreprocess}
\alias{InterMetaGenePair}
\alias{InterMetaCox}
\alias{InterMetaBG}
\alias{InterMetaPathwayPair}
\alias{InterMetaEnrich}
\alias{InterMetaPermFDR}
\alias{InterMetaPathwayNetwork}
\alias{InterMetaPathwayScore}
\alias{InterMetaDotplot}
\alias{InterMetaPathwayScoreBar}
\alias{InterMetaStability}
\alias{InterMetaStabilityScoreBar}
\alias{InterMetaSummarizeStability}
\title{InterMeta Metabolic Pathway Interaction Workflow}
\description{
InterMeta adapts the TimiGP directional gene-pair framework to metabolic
pathway interaction analysis using bulk transcriptomic and survival data.
}
\usage{
InterMetaAvailableGenesets()

InterMetaLoadGeneset(dataset = "all", file = NULL)

InterMetaCheckEvent(
  info = NULL,
  event_col = 1,
  time_col = 2,
  min_events = 1,
  min_censored = 1
)

InterMetaPreprocess(
  pathwayset = NULL,
  marker = NULL,
  rna = NULL,
  cohort = NULL,
  log = TRUE,
  GMNorm = TRUE,
  low_expression_filter = TRUE,
  analysis_profile = "TimiGP-legacy"
)

InterMetaGenePair(
  rna = NULL,
  both = FALSE,
  cont = FALSE,
  chunk_size = 10000,
  analysis_profile = "TimiGP-legacy"
)

InterMetaCox(
  mps = NULL,
  info = NULL,
  p.adj = "BH",
  parallel = FALSE,
  core = 1,
  keep_mps = TRUE,
  analysis_profile = "TimiGP-legacy"
)

InterMetaBG(marker.pair = NULL)

InterMetaPathwayPair(pathwayset = NULL, dataset = NULL, core = 1)

InterMetaEnrich(
  gene = NULL,
  background = NULL,
  pathway.pair = NULL,
  p.adj = "BH",
  core = 1,
  keep_gene_pairs = c("shared", "none", "all")
)

InterMetaPermFDR(
  resdata = NULL,
  pathwayset = NULL,
  gene = NULL,
  background = NULL,
  niter = 100,
  core = 1,
  seed = 1,
  analysis_profile = "TimiGP-legacy"
)

InterMetaPathwayNetwork(
  resdata = NULL,
  select = NULL,
  pathwayset = NULL,
  condition = "Adjust.P.Value",
  cutoff = 0.05,
  export = TRUE,
  path = NULL
)

InterMetaPathwayScore(
  resdata = NULL,
  condition = "Adjust.P.Value",
  cutoff = 0.05,
  analysis_profile = "TimiGP-legacy"
)

InterMetaDotplot(
  resdata = NULL,
  select = NULL,
  condition = "Adjust.P.Value",
  cutoff = 0.05,
  top = 10L
)

InterMetaPathwayScoreBar(score = NULL, select = NULL, top = 5L)

InterMetaStability(
  rna = NULL,
  info = NULL,
  pathwayset = NULL,
  pathway.pair = NULL,
  n_resample = 3,
  sample_fraction = 0.8,
  seeds = NULL,
  p.adj = "BH",
  pair_condition = "PV",
  pair_cutoff = 0.05,
  enrich_p.adj = "BH",
  score_condition = "Adjust.P.Value",
  score_cutoff = 0.05,
  run_permutation = FALSE,
  permutation_niter = 100,
  run_network = FALSE,
  log = TRUE,
  GMNorm = TRUE,
  parallel = FALSE,
  core = 1,
  outdir = NULL,
  reuse = TRUE,
  stratify_event = FALSE,
  min_events = 5,
  min_censored = 5,
  continue_on_error = TRUE,
  keep_details = NULL,
  analysis_profile = "TimiGP-legacy"
)

InterMetaSummarizeStability(
  stability = NULL,
  scores = NULL,
  pathways = NULL,
  top = c(5, 10, 20),
  high.frequency = 0.8,
  moderate.frequency = 0.5,
  high.direction = 0.8,
  moderate.direction = 0.6
)

InterMetaStabilityScoreBar(
  stability = NULL,
  select = NULL,
  top = 5L,
  min_selection_frequency = 1,
  min_direction_consistency = 1,
  error = c("sd", "se", "none")
)
}
\arguments{
  \item{dataset}{Bundled geneset key or a Dataset value used for filtering.}
  \item{file}{Optional path to a custom tab-delimited geneset.}
  \item{info}{Clinical data frame with sample IDs as row names.}
  \item{event_col}{Event column name or index; values must be 0 or 1.}
  \item{time_col}{Follow-up time column name or index.}
  \item{min_events}{Minimum number of observed events.}
  \item{min_censored}{Minimum number of censored observations.}
  \item{pathwayset}{Data frame containing Pathway (or CellType), Gene, and
  optionally Dataset columns.}
  \item{marker}{Optional character vector of genes used instead of pathwayset.}
  \item{rna}{Numeric gene-by-sample expression matrix.}
  \item{cohort}{Optional sample IDs to retain.}
  \item{log}{Whether to apply \code{log10(x + 1)}; defaults to TRUE.
  Genes containing values less than or equal to -1 are removed before logging.
  Set FALSE explicitly when the documented input scale must not be logged again.}
  \item{GMNorm}{Whether to subtract each gene's median across samples.}
  \item{low_expression_filter}{Whether to retain genes positive in at least
  half of samples. The rule also applies to matrices containing negative values.}
  \item{analysis_profile}{Only \code{"TimiGP-legacy"} is supported, and it is
  the default. The argument remains for explicit provenance. Other values,
  including the former \code{"InterMeta"} mode, raise an error.}
  \item{both}{Whether to retain both directions of each gene pair.}
  \item{cont}{Whether to use continuous expression differences instead of
  logical comparisons.}
  \item{chunk_size}{Number of gene pairs processed per memory-bounded chunk.}
  \item{mps}{Gene-pair score matrix.}
  \item{p.adj}{Multiple-testing adjustment method.}
  \item{parallel}{Whether Cox models should use parallel workers.}
  \item{core}{Positive number of worker processes.}
  \item{keep_mps}{Whether \code{InterMetaCox()} returns the direction-adjusted
  gene-pair score matrix. Set this to \code{FALSE} in memory-constrained
  workflows that only require Cox results.}
  \item{marker.pair}{Character vector of directional gene-pair identifiers.}
  \item{pathway.pair}{Compact object returned by
  \code{InterMetaPathwayPair()}. Expanded legacy data frames are also accepted
  by \code{InterMetaEnrich()}.}
  \item{gene}{Selected prognostic gene-pair identifiers. Use
  \code{character(0)} when no pairs were selected.}
  \item{background}{Background directional gene-pair identifiers.}
  \item{keep_gene_pairs}{Whether enrichment output stores shared pairs, no
  pair strings, or both shared and all annotated pairs.}
  \item{resdata}{Pathway-interaction enrichment result.}
  \item{niter}{Number of annotation permutations.}
  \item{seed}{Retained for API compatibility. Permutations always use
  the fixed iteration seeds \code{1:niter}; this argument does not change them.}
  \item{select}{Explicit enrichment indices, pathway-score indices, or
  pathway names to select, depending on the function. Explicit selections
  override automatic top-pathway selection.}
  \item{condition}{Significance column used to select interactions.}
  \item{cutoff}{Strict significance cutoff.}
  \item{export}{Whether network files should be written.}
  \item{path}{Output directory for network files.}
  \item{score}{Pathway score data frame.}
  \item{n_resample}{Number of repeated internal subsamples. The default of
  three implements the package's prespecified internal consistency check.
  Report the number, sampling fraction, seeds, and stratification rule; do not
  interpret this check as external validation or as a population-level
  uncertainty estimate.}
  \item{sample_fraction}{Fraction of samples retained in each subsample.}
  \item{seeds}{Deterministic seed vector with at least \code{n_resample}
  entries. When omitted, seeds start at 20260701.}
  \item{pair_condition}{Cox selection column. \code{InterMetaStability()}
  defaults to \code{"PV"} for exploratory unadjusted P-values. Set
  \code{"QV"} explicitly to use adjusted P-values for a confirmatory or
  sensitivity analysis.}
  \item{pair_cutoff}{Gene-pair significance cutoff.}
  \item{enrich_p.adj}{Adjustment method for pathway-interaction enrichment.}
  \item{score_condition}{Enrichment significance column used for scoring.}
  \item{score_cutoff}{Pathway-interaction cutoff used for scoring.}
  \item{run_permutation}{Whether to calculate permutation FDR in every run.}
  \item{permutation_niter}{Permutation count used in every run.}
  \item{run_network}{Whether to create a network for every run.}
  \item{outdir}{Optional checkpoint and result directory. Use a different
  directory for each concurrently running call; one output directory supports
  only one writer at a time.}
  \item{reuse}{Whether matching per-run checkpoints may be reused. Cached
  successful results are accepted only when their input, actual pathway-pair
  annotation, parameter, seed, and package-version fingerprints match. Failed
  checkpoints are retained for diagnosis but retried on the next call.}
  \item{stratify_event}{Whether subsampling preserves event/censoring balance.
  Defaults to FALSE; an explicit TRUE remains available as a sampling parameter.}
  \item{continue_on_error}{Whether a failed run is recorded while remaining
  runs continue.}
  \item{keep_details}{Whether large Cox and enrichment tables remain in the
  returned object. By default they remain in memory only when no output
  directory is supplied; complete per-run RDS files are always checkpointed.}
  \item{stability}{Object returned by \code{InterMetaStability()}.}
  \item{scores}{Optional list of pathway-score tables.}
  \item{pathways}{Complete candidate pathway vector, including pathways never
  selected in any run.}
  \item{top}{For \code{InterMetaDotplot()}, the total number of qualifying
  interactions to display in enrichment \code{Index} order. For the two score
  bars, the number of pathways per prognostic direction; use \code{NULL} to
  disable automatic limiting. For \code{InterMetaSummarizeStability()},
  magnitude-rank cutoffs used for top-set recurrence and Jaccard indices.}
  \item{high.frequency}{Selection and top-rank frequency required for high
  stability.}
  \item{moderate.frequency}{Selection and top-rank frequency required for
  moderate stability.}
  \item{high.direction}{Conditional direction consistency required for high
  stability.}
  \item{moderate.direction}{Conditional direction consistency required for
  moderate stability.}
  \item{min_selection_frequency}{Minimum fraction of successful resamples in
  which a pathway must receive a score before stability plotting.}
  \item{min_direction_consistency}{Minimum conditional agreement of the
  favorable or unfavorable direction before stability plotting.}
  \item{error}{Uncertainty shown by \code{InterMetaStabilityScoreBar()}:
  sample standard deviation, standard error, or no error bar. SD and SE require
  at least two successful score observations for every displayed pathway; use
  \code{"none"} for a single-run diagnostic summary.}
}
\value{
Functions return validated genesets, preprocessed matrices, Cox results,
compact pathway annotations, enrichment and score tables, plots, networks, or
a stability object. The stability object contains pathway summaries, pairwise
score comparisons, a long table of every resample score,
pathway-interaction recurrence, per-run status, the active-run manifest,
parameters, and provenance.
}
\details{
The only supported mode is \code{analysis_profile = "TimiGP-legacy"}.
It defaults to \code{log = TRUE}; choose the explicit log setting according to
the documented input scale. It retains genes positive in at least half of the selected
samples, removes any retained gene containing a value less than or equal to
-1 before applying \code{log10(x + 1)}, and median-centers each gene. It uses
\code{round(0.1 * n)} for the minimum prevalence of both gene-pair states and,
for logical matrices, fits the Cox predictor as a factor. These operations
reproduce the original result while avoiding the original propagation of
\code{NaN} values.

The compact pathway-pair representation avoids materializing millions of
repeated annotation rows. Enrichment counts are calculated with gene adjacency
and pathway membership matrices. Stability sampling defaults to unstratified 80-percent
sampling, seeds beginning at 20260701, \code{PV < 0.05}, and sorted sampled
IDs. Cached results are reused only after fingerprint validation.
Increasing the requested number of resamples can reuse existing compatible
per-seed checkpoints. Failed checkpoints are never treated as completed runs.

Stability Top-K sets are profile-independent: eligible pathways have finite
\code{Difference}, are ordered by decreasing \code{abs(Difference)}, and ties
are broken by ascending pathway name. Overall, favorable-only, and
unfavorable-only frequencies and Jaccard indices are reported separately.
Empty-versus-empty overlap is undefined (\code{NA}), not perfect agreement.
\code{Stability.Class} uses the largest requested K (Top 20 by default), and
is \code{"Not assessed"} when fewer than two runs succeeded. Use the explicit
\code{Top5.*.Frequency} fields for a Top-5 stability statement.

Pathway scoring removes pathways
whose \code{Difference} is zero and uses signed \code{Difference} for both
rank columns. Its permutation FDR also restores iteration seeds
\code{1:niter}, the original enrichment \code{Rank} denominator, the
shuffle-before-de-duplication order of the original marker slots, and the
unmodified (non-monotone) FDR values. The bundled all-pathway geneset exposes
1,999 unique pathway-gene rows for ordinary analysis while retaining 2,022 raw
marker slots internally for permutation.

The installed complete PRAD example and ordinary function calls use the same
\code{"TimiGP-legacy"} mode. With
the bundled TCGA-PRAD inputs this mode restores 861 full-cohort genes,
370,230 gene pairs, and 50 scored pathways; pathway identity, order, ranks,
groups, and scores match the original TimiGP result to about \code{5.33e-15}.
The earlier 0.4.0 server run instead retained 1,260 genes, used \code{QV},
event-stratified seeds 1--3, and absolute ranks with neutral pathways. Because
these differences affect both Cox fitting and sample membership, exact legacy
equivalence requires recomputation and cannot be recovered by post-processing
the 0.4.0 output alone.

\code{InterMetaDotplot()} defaults to the first ten qualifying interactions
in enrichment \code{Index} order. \code{InterMetaPathwayScoreBar()} defaults to
the five strongest favorable and five strongest unfavorable pathways.
\code{InterMetaStabilityScoreBar()} shows the mean Difference plus or minus
its standard deviation across resamples, after requiring complete selection
and direction consistency by default.
}
\examples{
\dontrun{
library(InterMeta)
pathwayset <- InterMetaLoadGeneset("all")
pathway_pair <- InterMetaPathwayPair(pathwayset)
}
}
\keyword{package}
