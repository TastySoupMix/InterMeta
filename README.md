# InterMeta 1.1.4

`InterMeta` infers directional prognostic interactions between metabolic
pathways from bulk transcriptomic profiles and survival data. It adapts the
TimiGP gene-pair framework using metabolic pathway gene sets.

## Start here

The main worked example is **TCGA-PRAD**. After installing the package, use
the [bundled example](#tcga-prad-example) to inspect its results or run the
complete analysis. For your own cohort, follow the
[full-sample workflow](#full-sample-workflow).

Detailed guides cover [getting started](inst/doc/getting-started.md),
[PRAD reproduction](inst/doc/prad-reproduction.md), and the
[optional PubMed/AI workflow](inst/doc/ai-workflow.md).

## Analysis overview

The workflow selects prognostic gene pairs, tests their enrichment in
directional pathway interactions, and calculates signed pathway scores.
Outputs include interaction networks, chord diagrams, enrichment dot plots,
and pathway score plots. Internal resampling summarizes the repeatability of
pathway scores and interactions. An optional PubMed/AI workflow retrieves and
summarizes literature for selected pathways.

## Installation

Install InterMeta 1.1.4 from GitHub in R or RStudio:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "TastySoupMix/InterMeta",
  ref = "v1.1.4",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE
)

library(InterMeta)
packageVersion("InterMeta")
```

A downloaded source folder can also be installed with
`remotes::install_local("/path/to/InterMeta")`.

## Bundled Metabolic Gene Sets

InterMeta includes six gene sets:

| Key | Scope | Provenance |
| --- | --- | --- |
| `glucose` | Glucose metabolism | KEGG-derived subset curated by Zeheng Tan |
| `lipid` | Lipid metabolism | KEGG-derived subset curated by Zeheng Tan |
| `amino` | Amino-acid metabolism | KEGG-derived subset curated by Zeheng Tan |
| `all` | Integrated energy metabolism | KEGG-derived selection, organization, and merging by Zeheng Tan |
| `kegg` | All KEGG metabolism-related pathways | Converted database resource |
| `reactome` | All Reactome metabolism-related pathways | Underlying Reactome data are CC0 |

The four energy-metabolism collections were selected, organized, and merged
by Zeheng Tan from KEGG-derived annotations. The full KEGG and Reactome files
are content-equivalent to the corresponding scMetabolism 0.2.1 GMT resources
after the recorded label conversions. Sources, curation, available release
information, and file fingerprints are documented in the
[gene-set provenance record](inst/extdata/metabolic-genesets/PROVENANCE.md).

Load a bundled set:

```r
pathwayset <- InterMetaLoadGeneset("all")
```

Or supply a custom tab-delimited file:

```r
pathwayset <- InterMetaLoadGeneset(file = "my_geneset.txt")
```

## TCGA-PRAD example

The package includes TCGA-PRAD prognosis and normalized expression data
downloaded from [PCaDB](https://github.com/rli012/PCaDB): 463 samples,
99 biochemical-recurrence events, and expression for 1,260 metabolic genes.
Follow-up time is measured in days. Precomputed results include 2,652 directed
pathway interactions and a 50-row pathway score table.
`InterMetaPRADProvenance` records the input fingerprints and result provenance.

```r
data("InterMetaPRADInfo")
data("InterMetaPRADRNA")
data("InterMetaPRADEnrich")
data("InterMetaPRADScore")
data("InterMetaPRADProvenance")
```

Three scripts provide different entry points for the same cohort:

| Script | Purpose |
| --- | --- |
| [example05_tcga_prad_results.R](inst/examples/example05_tcga_prad_results.R) | Inspect precomputed results and generate plots. |
| [example06_tcga_prad_recompute.R](inst/examples/example06_tcga_prad_recompute.R) | Recompute Cox models, enrichment, permutation FDR, and pathway scores. |
| [example09_tcga_prad_complete_workflow.R](inst/examples/example09_tcga_prad_complete_workflow.R) | Run the full analysis, internal resampling, and optional PubMed/AI workflow. |

For a quick look at the precomputed results:

```r
source(system.file(
  "examples",
  "example05_tcga_prad_results.R",
  package = "InterMeta"
))
```

The display script uses `Permutation.FDR < 0.2` for selected interaction
plots; its saved pathway score uses `Adjust.P.Value < 0.05`.

For the complete workflow, review the configuration block and run:

```r
source(system.file(
  "examples",
  "example09_tcga_prad_complete_workflow.R",
  package = "InterMeta"
))
```

This reproduces the bundled reference with `log = TRUE`, `GMNorm = TRUE`,
gene-pair selection at `PV < 0.05`, pathway scoring at `Adjust.P.Value < 0.05`,
100 annotation permutations per analysis, and three 80% internal resamples.
The full-cohort reference retains 861 genes and tests 370,230 gene pairs.
Model calls are disabled in the script until explicitly enabled.
For data-source and reproduction details, see the
[PRAD guide](inst/doc/prad-reproduction.md).

## Input requirements

Provide a numeric gene-by-sample expression matrix with unique gene row names
and sample column names. Clinical data use sample IDs as row names, event
status (`0` or `1`) in the first column, and follow-up time in the second.

Preprocessing defaults to `log = TRUE` and `GMNorm = TRUE`. Set `log = FALSE`
for expression values already on a logarithmic scale. The examples below use
`log = FALSE`; choose this setting from your input data's documented scale.

## Full-sample workflow

Starting with `bulk_expression` and `clinical_info`:

```r
library(InterMeta)

pathwayset <- InterMetaLoadGeneset("all")
info <- InterMetaCheckEvent(clinical_info)

rna <- InterMetaPreprocess(
  pathwayset = pathwayset,
  rna = bulk_expression,
  cohort = rownames(info),
  log = FALSE,
  GMNorm = TRUE
)

mps <- InterMetaGenePair(
  rna,
  chunk_size = 10000
)
cox <- InterMetaCox(
  mps = mps,
  info = info,
  p.adj = "BH",
  parallel = TRUE,
  core = 8,
  keep_mps = FALSE
)

cox_res <- cox$cox_res
query_pairs <- rownames(cox_res)[
  is.finite(cox_res$PV) & cox_res$PV < 0.05
]
background <- InterMetaBG(rownames(cox_res))

pathway_pair <- InterMetaPathwayPair(pathwayset)
enrich <- InterMetaEnrich(
  gene = query_pairs,
  background = background,
  pathway.pair = pathway_pair,
  p.adj = "BH",
  keep_gene_pairs = "shared"
)

score <- InterMetaPathwayScore(
  enrich,
  condition = "Adjust.P.Value",
  cutoff = 0.05
)
```

This example selects gene pairs using unadjusted `PV < 0.05` and scores
pathways using BH-adjusted enrichment P-values below 0.05. Use `QV` in the
gene-pair selection step for an adjusted sensitivity analysis.

### Permutation FDR

```r
enrich <- InterMetaPermFDR(
  resdata = enrich,
  pathwayset = pathwayset,
  gene = query_pairs,
  background = background,
  niter = 100,
  core = 8,
  seed = 1
)
```

`InterMetaPermFDR()` parallelizes annotation permutations. Iteration seeds are
fixed at `1:niter`, so the `seed` argument does not change that sequence.

## Visualization

Plot the leading pathway interactions and pathway scores:

```r
dotplot <- InterMetaDotplot(
  enrich,
  condition = "Adjust.P.Value",
  cutoff = 0.05,
  top = 10
)
score_bar <- InterMetaPathwayScoreBar(score, top = 5)
```

Draw a directional pathway chord diagram:

```r
pdf("pathway_chord.pdf", width = 14, height = 14, useDingbats = FALSE)
chord <- InterMetaPathwayChord(
  enrich,
  condition = "Adjust.P.Value",
  cutoff = 0.05,
  max_edges = NULL,
  link_weight = "No.Shared.MGP",
  label_cex = 0.7,
  start_degree = 90,
  adaptive_labels = FALSE,
  canvas_limit = 1.45,
  sector_floor_quantile = NULL
)
dev.off()

print(chord)
```

`max_edges = NULL` draws all interactions meeting the cutoff. Use `select`
with enrichment `Index` values for a specified set, or set `max_edges` for a
smaller diagram. The `group`, `color`, and `order` arguments control pathway
grouping, colors, and arrangement.

To inspect the gene pairs behind one pathway interaction:

```r
gene_chord_candidates <- which(
  !is.na(enrich$Shared.MGP) & nzchar(enrich$Shared.MGP) &
    enrich$No.Shared.MGP > 0 & enrich$No.Shared.MGP <= 50
)
gene_chord_candidates <- gene_chord_candidates[
  order(enrich$Adjust.P.Value[gene_chord_candidates])
]
gene_chord_index <- enrich$Index[gene_chord_candidates[1]]
pdf("pathway_gene_pairs.pdf", width = 9, height = 9)
gene_chord <- InterMetaGeneChord(
  enrich,
  select = gene_chord_index,
  condition = "Adjust.P.Value"
)
dev.off()
```

Gene-level chords require enrichment with `keep_gene_pairs = "shared"` or
`"all"`, and at least one interaction meeting the selection criteria. For
precomputed objects without `Shared.MGP`, recompute enrichment first.

## Resampling stability

`InterMetaStability()` repeats the analysis on subsets of the cohort and
summarizes pathway selection, score direction, Top-K overlap, and interaction
recurrence. Defaults are three unstratified 80% subsets with seeds
20260701--20260703. These results describe internal repeatability rather than
external validation or population-level uncertainty.

```r
stability <- InterMetaStability(
  rna = bulk_expression,
  info = clinical_info,
  pathwayset = pathwayset,
  n_resample = 3,
  sample_fraction = 0.8,
  seeds = 20260701:20260703,
  pair_condition = "PV",
  pair_cutoff = 0.05,
  score_condition = "Adjust.P.Value",
  score_cutoff = 0.05,
  run_permutation = FALSE,
  log = FALSE,
  GMNorm = TRUE,
  parallel = TRUE,
  core = 8,
  outdir = "InterMeta_stability",
  reuse = TRUE,
  stratify_event = FALSE,
  min_events = 5,
  min_censored = 5
)
```

Use `run_permutation = TRUE, permutation_niter = 100` to calculate permutation
FDR in each resample. The complete PRAD script enables this step.
Completed runs are checkpointed and reused when their data and settings match;
failed runs can be retried. Use a separate output directory for each analysis.

The result includes `summary` for pathway stability, `pairwise` for score
correlations and Top-K overlap, `edge_summary` for interaction recurrence,
and `run_status` for completion and sample counts. See `?InterMetaStability`
for the full result structure and Top-K definitions.

Plot mean pathway scores with sample standard deviations across resamples:

```r
stability_plot <- InterMetaStabilityScoreBar(
  stability,
  top = 5,
  min_selection_frequency = 1,
  min_direction_consistency = 1,
  error = "sd"
)
ggplot2::ggsave(
  "stability_pathway_score_bar.pdf",
  stability_plot,
  width = 8.5,
  height = 5.5
)
```

The plot above includes pathways selected in every resample with a consistent
score direction. Standard-deviation bars require at least two successful runs.

## Optional PubMed and AI workflow

The AI workflow builds a pathway-level context, retrieves PubMed records,
screens abstracts, verifies cited excerpts, and synthesizes a report and
figure. `InterMetaAIWorkflow()` uses `biology_exploration` to summarize
biological findings for the selected disease and pathways. For endpoint-specific
literature comparison, use `InterMetaAIInterpret()` separately with
`purpose = "prognostic_validation"`; see `?InterMetaAIInterpret`.

The following example uses the bundled PRAD pathway score. Set
`DEEPSEEK_API_KEY` in the R process environment before running it. Credentials
stay outside scripts and reports; model calls require
`allow_model_calls = TRUE` and may incur provider charges.

```r
library(InterMeta)

data("InterMetaPRADScore")

context <- InterMetaAIContext(
  score = InterMetaPRADScore,
  disease = "prostate adenocarcinoma",
  endpoint = "biochemical recurrence-free survival",
  time_unit = "days",
  cohort = "TCGA-PRAD",
  top_n = 5L,
  direction = "unfavorable"
)

if (!nzchar(Sys.getenv("DEEPSEEK_API_KEY"))) {
  stop("Set DEEPSEEK_API_KEY outside the analysis script, then restart R.")
}

run_dir <- "InterMeta_AI_PRAD_1_1_4"

ai_result <- InterMetaAIWorkflow(
  context = context,
  run_dir = run_dir,
  allow_model_calls = TRUE,
  language = "zh",
  purpose = "biology_exploration",
  model = "deepseek-v4-pro",
  screen_batch_size = 20L,
  max_findings_per_pathway = 5L,
  render = TRUE,
  figure_text_mode = "compact_english"
)

InterMetaAIStatus(run_dir)
ai_result$report$exploration_summary
```

Use a new `run_dir` for the first call. Inspect progress with
`InterMetaAIStatus(run_dir)` and continue an interrupted run with:

```r
ai_result <- InterMetaAIResume(
  run_dir,
  allow_model_calls = TRUE,
  retry_failed = FALSE
)
```

Review failed or delivery-uncertain calls before setting `retry_failed = TRUE`,
which permits another model request. Replay a completed run offline with:

```r
replay <- InterMetaAIReplay(
  run_dir,
  outdir = "InterMeta_AI_PRAD_replay",
  render = TRUE
)
```

Reports contain `exploration_summary.tsv`, `exploration_findings.tsv`,
retrieved evidence, validation results, and provenance. The model receives
pathway-level summaries and literature rather than patient-level expression
or survival records. PubMed searches are bounded; results that reach a query
limit are marked as truncated.

AI summaries require scientific review of the cited studies and pathway
scope. Excerpt matching verifies the text source, not the biological claim;
the output is exploratory and does not establish causality, clinical utility,
or independent validation. The
[AI guide](inst/doc/ai-workflow.md) covers interpretation, retrieval limits,
credentials, recovery, and replay. See `?InterMetaAIWorkflow` and
`?InterMetaAIScreening` for checkpoint and validation details.

### Model providers

DeepSeek is the default provider. Adapters are also available for OpenAI,
Anthropic, Google, Alibaba, and Moonshot. To inspect presets and construct
a provider:

```r
InterMetaAIModels()
provider <- InterMetaLLMProvider(provider = "deepseek")
# Or select another service and a model available to your account:
provider <- InterMetaLLMProvider(provider = "openai", model = "gpt-4.1")
```

Pass `provider` to `InterMetaAIWorkflow()` and again when resuming a model
stage. Configure the service's environment variable as documented in
`?InterMetaLLMProvider`. Presets are configuration examples; service and model
availability depend on the provider account.

The package also includes
[example08_ai_workflow.R](inst/examples/example08_ai_workflow.R), which runs
the AI workflow from precomputed PRAD pathway scores.

## Citation and license

Use `citation("InterMeta")` for the software citation. Cite the TimiGP method
and the original databases used in your analysis. The bundled PRAD data were
downloaded from PCaDB; cite the
[PCaDB preprint](https://doi.org/10.1101/2021.06.29.449134) when using them.

InterMeta is distributed under GPL-3; see [COPYING](COPYING). Zeheng Tan is the
package author and maintainer. Portions of `R/Chord.R`, `R/Cox.R`, and
`R/InterMeta.R` were adapted from [TimiGP](https://github.com/CSkylarL/TimiGP).
[inst/NOTICE](inst/NOTICE) records the upstream attribution and modifications.
Database-derived gene sets retain their source terms, documented in the
[gene-set provenance record](inst/extdata/metabolic-genesets/PROVENANCE.md).
