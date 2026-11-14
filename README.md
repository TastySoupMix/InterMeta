# InterMeta 1.1.4

`InterMeta` infers directional prognostic interactions between metabolic
pathways from bulk transcriptomic profiles and survival data. It adapts the
TimiGP gene-pair framework while replacing immune-cell annotations with
metabolic pathway gene sets.

Short, task-oriented guides are available for
[getting started](inst/doc/getting-started.md),
[reproducing the bundled PRAD analysis](inst/doc/prad-reproduction.md), and
[running the optional PubMed/AI workflow](inst/doc/ai-workflow.md).

## Start here

Install from GitHub with `remotes::install_github("TastySoupMix/InterMeta")`,
then load `library(InterMeta)`. See the installation section below for
dependencies and installation of a fixed version.

The principal worked research example is **TCGA-PRAD**. The small synthetic
cohort is a separate offline software check and does not represent TCGA data.

- Full TCGA-PRAD recomputation, internal resampling and optional AI:
  `inst/examples/example09_tcga_prad_complete_workflow.R`.
  Read its configuration block before sourcing it; model calls require an
  explicit opt-in and a credential supplied outside the script.
- Fast inspection of previously computed PRAD results:
  `inst/examples/example05_tcga_prad_results.R`. Its selected interaction
  plots use `Permutation.FDR < 0.2`; its saved score uses
  `Adjust.P.Value < 0.05`. It is not a new full analysis.
- Small offline installation check:
  `inst/examples/example01_metabolic_pathway_prognosis.R`.
  Under the fixed 1.1.4 example settings its minimum permutation FDR is 0.09;
  the 0.05 filter therefore produces an empty pathway score and no selected
  pathway chord links. This is an expected nonsignificant result.

Use `citation("InterMeta")` for the software citation, and cite the TimiGP
method and the data sources relevant to the analysis. See
[the example guide](inst/doc/getting-started.md) and
[AI interpretation limits](inst/doc/ai-workflow.md).

## Analysis overview

InterMeta 1.1.4 supports one analysis mode, `TimiGP-legacy`. The workflow
includes gene-pair Cox analysis, pathway-interaction enrichment, permutation
FDR, signed pathway scores, and network and chord plots. Internal resampling
defaults to three 80% subsets. The bundled TCGA-PRAD and synthetic examples
use this analysis mode.

An optional PubMed/AI workflow supports literature screening, synthesis,
checkpoint recovery, offline replay, and local report and figure generation.

## Optional PubMed/AI workflow

For a new workflow, the default PubMed cache is saved in fixed configuration
as the relative path `pubmed_cache`, anchored to the workflow run's parent
directory. Synthesis checkpoint references, replay-source metadata, and the
complete example's human-readable AI-run summary are likewise stored relative
to their documented anchors. This allows a workflow tree to be moved between
machines without retaining obsolete server paths. If a relocated run
stopped before PubMed retrieval completed,
`InterMetaAIResume(pubmed_cache_dir = ...)` can provide an operational cache
location without rewriting the fixed configuration or changing its hash. The
direct synthesis aggregate resolves `calls/...` from its `outdir`; the workflow
aggregate resolves `checkpoints/calls/...` from the directory containing
`05_synthesis/synthesis.rds`. Offline replay verifies the frozen configuration
SHA before publishing a report.

In the optional AI workflow, R builds and runs reproducible PubMed searches;
the first model stage
screens bounded, pathway-specific article batches; R then checks candidate IDs,
structured labels, disease-label consistency, and contiguous Abstract excerpts.
Only locally verified candidates enter the second model stage, which summarizes
one pathway at a time and selects existing candidate IDs. R maps citations and
quotes back from the frozen screening table, validates the merged report, and
renders the fixed publication figure locally.

The run directory records all eight stages: context, queries, PubMed retrieval,
screening, verified evidence, synthesis, report, and figure. Completed artifacts
are hash-bound and reusable. Model calls are disabled by default, failed or
delivery-uncertain calls are never retried automatically, and an intact
completed run can be revalidated and redrawn without PubMed or model access.
Invalid screening rows are quarantined individually as `unresolved` instead of
being silently repaired. When no whole batch failed, screening with such rows
is marked `complete_with_quarantine`: it remains usable, but quarantined rows
never enter synthesis. Failed or incomplete whole batches remain unusable. Raw
returned model responses are retained unchanged as pre-validation artifacts,
so a compatible new
validation revision can re-check them as a new immutable attempt with zero new
model calls. Report and figure publications use numbered attempt directories,
preventing a later revalidation from reusing stale output. These controls
improve traceability; they do not turn an AI summary into proof of causality,
pathway activity, clinical utility, or independent validation of the InterMeta
score direction.

A PubMed stage that reaches one or more configured per-query limits is
marked `complete_with_truncation`. It remains usable as an explicitly bounded
evidence snapshot, but must not be described as an exhaustive search. In
`biology_exploration`, the authoritative result tables are
`exploration_summary.tsv` and `exploration_findings.tsv`; prognostic
compatibility tables are labelled `not_applicable` and do not signal model
failure or insufficient evidence. Reports use schema 4.3.0.

The following example starts from the package's public, precomputed TCGA-PRAD
pathway score. Replace only the context labels and score object when adapting it
to another cohort. Set `DEEPSEEK_API_KEY` in the server environment before
starting R; never paste a credential into this script or the package.

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

`run_dir` must be new on the first call. After interruption or a normal pause,
reuse completed work with:

```r
ai_result <- InterMetaAIResume(
  run_dir,
  allow_model_calls = TRUE,
  retry_failed = FALSE
)
```

Only after moving an unfinished run whose cache is not the conventional
sibling directory, add
`pubmed_cache_dir = "/new/path/to/pubmed_cache"` to that call. It is an
operational override and is not saved into fixed workflow configuration. The
override is also useful when the saved cache path is ambiguous after relocation.

If status shows a failed or delivery-uncertain model attempt, inspect the saved
stage and provider status first. A new potentially billable attempt requires the
deliberate combination `allow_model_calls = TRUE, retry_failed = TRUE`. To
revalidate a completed run and draw a new figure without any external call:

```r
replay <- InterMetaAIReplay(
  run_dir,
  outdir = "InterMeta_AI_PRAD_replay",
  render = TRUE
)
```

For lower-level control, call `InterMetaAIScreenEvidence()`,
`InterMetaAIScreenedEvidence()`, and `InterMetaAISynthesizeScreened()` directly.
See `?InterMetaAIWorkflow` and `?InterMetaAIScreening` for checkpoint and retry
semantics.

When a package update changes only a compatible local validation revision,
`InterMetaAIResume()` first revalidates the saved screening or synthesis
response without a provider request. Newly published workflow files are kept
under `06_report/outputs/attempt-NNN` and
`07_figure/outputs/attempt-NNN`; earlier attempts remain auditable.

### Evidence interpretation and validation

Reports include finding-level quote quarantine, pathway-level
`Quote.Grounding`, AI contract and validation revisions in provenance, and a
16-file report bundle. Saved reports retain their original generation version
and assessment provenance; changing the installed package does not regenerate
a report.

The default post-analysis purpose is `biology_exploration`.
Disease/pathway retrieval is not gated by the statistical endpoint. Findings
describe functional experiments, omics/bioinformatics and patient studies as
parallel, non-exclusive sources, with biological processes, model context,
component/pathway scope, null results and testable hypotheses. The report
includes `exploration_summary.tsv` and `exploration_findings.tsv`.
New model calls request at most five non-redundant findings per pathway while
retaining the larger selected literature set as provider-visible evidence.
Each non-empty exploration quote is tested as one contiguous verbatim excerpt
and must contain at least 20 Unicode letter-or-number characters to qualify.
The provider-visible prompt requests eight characters; the independent local
validator applies the stricter 20-character rule.
Validation first collapses whitespace, then permits only documented Unicode
typography equivalents. A non-contiguous or otherwise source-mismatched quote is
not repaired or fuzzy-matched: its model text remains in the audit row with
`Quote.Verified = FALSE`, `Quote.Match.Method = source_mismatch`, and
`Eligible.For.Summary = FALSE`. A shorter non-empty quote is retained as
`insufficient_information`; a quote longer than 1,000 characters is retained as
`overlength`. None contributes to finding-derived evidence channels, PMID or
process aggregates, or source-linked observations. The summary's
`Quote.Grounding` is `complete`, `partial`, `none`, or `no_findings` and describes
strict quote matching only. `Evidence.Eligibility` uses the same four states,
derived separately from the eligible/all finding count after strict quote
matching and consistency of the model-supplied `Study.Disease` and
`Disease.Match` fields with the target label. The finding table exposes this
second gate as `Disease.Match.Verified`. It is structural and lexical; it does
not independently infer the paper's disease from the Abstract or MeSH. When only
some findings are eligible, the public `Biological.Summary` is built solely from
eligible `Source.Linked.Observations`, public `Context.Dependence` and
`Research.Hypothesis` are suppressed, and the original `Model.*` narratives
remain audit-only. No case
folding, transliteration, scientific-symbol substitution, word editing, or
fuzzy matching is used.
Unique-PMID counts are not quality scores or independent-cohort counts.
Explicit `purpose = "prognostic_validation"` retains the endpoint-specific mode.

Pathway-level biological syntheses are distinct from strict endpoint-specific
prognostic consistency. Mechanism findings retain their
supporting abstracts, disease applicability, pathway/component scope, and
contrary or null results; they cannot establish pathway activation or validate
the sign of a gene-pair network score. Reports without these assessments
remain not assessed. Model syntheses require human scientific review.

The optional PubMed evidence layer supports DeepSeek, OpenAI, Anthropic
(Claude), Google (Gemini), Alibaba (Qwen), and Moonshot (Kimi). Each user
provides their own credentials; no API key or shared account is bundled.
It deterministically selects the Top 5 unfavorable and/or favorable pathways,
keeps patient-level inputs out of model payloads, validates every model-cited PMID against
retrieved records, and publishes a traceable report bundle. Endpoint-aware
evidence grading and separately labelled component retrieval are available.

## Installation

Install from GitHub in R or RStudio:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "TastySoupMix/InterMeta",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE
)

library(InterMeta)
packageVersion("InterMeta")
```

To install the fixed 1.1.4 snapshot:

```r
remotes::install_github(
  "TastySoupMix/InterMeta",
  ref = "v1.1.4",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE
)
```

A downloaded source folder can also be installed with
`remotes::install_local("/path/to/InterMeta")`.

The expected package version is `1.1.4`.

## License, authorship, and upstream code

The package code is distributed under GPL-3; the complete license text is in
[COPYING](COPYING). Zeheng Tan is the package author and maintainer. The final
manuscript author order and affiliations will be added when confirmed and
should not be inferred from the software metadata.

InterMeta is a metabolic-pathway extension of TimiGP. Portions of
`R/Chord.R`, `R/Cox.R`, and `R/InterMeta.R` were adapted from the
[TimiGP project](https://github.com/CSkylarL/TimiGP); file headers and
[inst/NOTICE](inst/NOTICE) record the mapping, modifications, upstream author,
and article DOI. InterMeta adds metabolic pathway adaptation, stricter
validation and compatibility controls, internal resampling/stability analysis,
and the optional AI-assisted evidence workflow. The NOTICE preserves
attribution without adding a use restriction to GPL-3.

The supplied upstream tree was verified file-for-file against official TimiGP
tag `v1.3.0`, commit `9a527dd7dfe18c3155647d57c4392dad977ff13f`.

The code license does not automatically relicense database-derived files.
Their sources, curation, and file fingerprints are recorded
separately in
[the gene-set provenance record](inst/extdata/metabolic-genesets/PROVENANCE.md).

## Choosing an AI model (optional)

The orchestrated workflow uses DeepSeek by default. A custom
provider created below can instead be passed as `provider` to
`InterMetaAIWorkflow()` and must be passed again when resuming a model stage.
The following lower-level sequence remains useful for inspecting a deterministic
query plan, frozen PubMed evidence, or the evidence-only report skeleton; it is
not a substitute for the screening-then-synthesis workflow.

After constructing an aggregate `context`:

```r
queries <- InterMetaBuildPubMedQueries(context, purpose = "biology_exploration")
# Explicit public-literature request; no model is called here:
evidence <- InterMetaPubMedSearch(context, queries, retmax = 100L,
  page_size = 100L, sort = "relevance", cache_dir = "pubmed_cache")
report <- InterMetaAIReport(context, evidence = evidence) # Evidence-only, offline
# Supply provider explicitly to request a potentially billable interpretation.
```

The search and model-selection caps are separate. The default model selection
is at most 30 articles per pathway, with a shared 200,000-character abstract
budget and recorded prefix truncation. Selection/type hints improve diversity
but are not verified study classifications and cannot guarantee comprehensive
coverage. Source-matched quotations validate provenance, not semantic truth.
Exploration outputs expose `Quote.Match.Method` as `whitespace_exact`,
`typography_equivalent`, `empty`, `insufficient_information`, or
`source_mismatch`, with `overlength` for a quote exceeding 1,000 characters; the
second category changes presentation only, not words, case or scientific
symbols. Only rows with `Quote.Verified = TRUE` and
`Disease.Match.Verified = TRUE` are eligible for derived evidence summaries.
The latter checks lexical consistency of model-supplied disease fields rather
than independently inferring study disease from Abstract or MeSH text. `Quote.Grounding` shows
whether all, some, none, or no finding rows have verified quotes, while
`Evidence.Eligibility` exposes a separate four-state result derived from
eligible/all after both gates. The main figure
shows quote quarantine separately from target-disease exclusion. Read the
original papers before publishing an interpretation.

```r
InterMetaAIModels() # Offline presets; not a live service availability check.
provider <- InterMetaLLMProvider(provider = "deepseek")
# Or choose another service, and optionally an explicit supported model ID:
provider <- InterMetaLLMProvider(provider = "openai", model = "gpt-4.1")
# Constructing a provider does not read credentials, contact a service, or bill.
# A request happens only when explicitly passed to InterMetaAIReport/AIInterpret.
```

Configure the selected service's own environment variable in the R process;
see `?InterMetaLLMProvider` for the mapping. Do not embed credentials in the
package, analysis script, saved provider object, or report. The catalog is a
dated set of examples, not an exhaustive model list or recommendation ranking.
Set `model` explicitly to an ID supported by your account when needed. The
existing `InterMetaDeepSeekProvider()` and custom `InterMetaAIProvider()` remain
available. Each report records the requested and returned model when supplied
by the service. Mocked protocol tests do not establish live account access.

## Bundled Metabolic Gene Sets

InterMeta 1.1.4 includes six gene sets:

| Key | Scope | Provenance |
| --- | --- | --- |
| `glucose` | Glucose metabolism | KEGG-derived subset curated by Zeheng Tan |
| `lipid` | Lipid metabolism | KEGG-derived subset curated by Zeheng Tan |
| `amino` | Amino-acid metabolism | KEGG-derived subset curated by Zeheng Tan |
| `all` | Integrated energy metabolism | KEGG-derived selection, organization, and merging by Zeheng Tan |
| `kegg` | All KEGG metabolism-related pathways | Converted database resource; exact release pending |
| `reactome` | All Reactome metabolism-related pathways | Underlying Reactome data are CC0; immediate compiled source/version pending |

The four energy-metabolism collections were selected, organized and merged
by Zeheng Tan from KEGG-derived pathway annotations. The full KEGG and
Reactome files are content-equivalent to the corresponding scMetabolism
0.2.1 GMT resources after the recorded label conversions. Exact original
database releases and download dates were not retained. Source links,
curation notes and file fingerprints are provided in the
[gene-set source record](inst/extdata/metabolic-genesets/PROVENANCE.md).
Please cite the original databases as well as InterMeta when using these sets.

Load the integrated set without uploading a separate text file:

```r
pathwayset <- InterMetaLoadGeneset("all")
```

A custom tab-delimited file can still be used:

```r
pathwayset <- InterMetaLoadGeneset(file = "my_geneset.txt")
```

## Bundled Examples

### Synthetic Workflow Check

InterMeta includes a generated cohort with 120 samples, 30 genes, and six
metabolic pathways. It is deliberately small enough to run on a personal
computer and contains no patient records. The signal is intentionally strong,
so the example is suitable for checking installation and workflow behavior,
not for biological inference.

Load the input objects directly:

```r
library(InterMeta)
data("InterMetaExampleInfo")
data("InterMetaExampleRNA")
data("InterMetaExamplePathwayset")
```

Run the complete Cox, enrichment, permutation, network, score, and plotting
workflow:

```r
source(system.file(
  "examples",
  "example01_metabolic_pathway_prognosis.R",
  package = "InterMeta"
))
```

Precomputed results are also available for immediate inspection:

```r
data("InterMetaExampleCox")
data("InterMetaExampleEnrich")
data("InterMetaExampleScore")
```

The same inputs are included in the text format used by the server workflow:

```r
example_dir <- system.file("extdata", "example", package = "InterMeta")
list.files(example_dir, full.names = TRUE)
```

Installed example scripts cover data loading, a complete analysis, internal
resampling stability, a real-data server run, and precomputed results. Locate
them with `list.files(system.file("examples", package = "InterMeta"))`.

### Real TCGA-PRAD Example

The package also contains the TCGA-PRAD cohort used in the original metabolic
analysis. Its prognosis and normalized bulk transcriptomic input tables were
downloaded from PCaDB (historical website endpoint recorded in the provenance
object; [source repository](https://github.com/rli012/PCaDB)); analyses and the
accompanying manuscript should identify PCaDB as the download source and cite
the [PCaDB preprint](https://doi.org/10.1101/2021.06.29.449134). The bundled
objects contain the aligned prognosis cohort and metabolic-gene expression
subset used by the examples:

- 463 samples with prognosis data: 99 events and 364 censored observations;
- normalized expression for 1,260 genes matched to `geneset-all`; the legacy
  preprocessing profile retains 861 genes for the full cohort;
- all 2,652 directed pathway interactions with BH-adjusted P-values and
  100-iteration permutation FDR;
- a compact migrated enrichment result, its 50-row TimiGP-legacy score table,
  and SHA-256 provenance fingerprints.

Load the real input data and precomputed results directly:

```r
data("InterMetaPRADInfo")
data("InterMetaPRADRNA")
data("InterMetaPRADEnrich")
data("InterMetaPRADScore")
data("InterMetaPRADProvenance")
```

Inspect the data, plots, network, and precomputed results without rerunning
Cox or permutation:

```r
source(system.file(
  "examples",
  "example05_tcga_prad_results.R",
  package = "InterMeta"
))
```

Recompute the analysis from the bundled PRAD inputs on a server:

```r
source(system.file(
  "examples",
  "example06_tcga_prad_recompute.R",
  package = "InterMeta"
))
```

Run the publication-package acceptance workflow from the same bundled inputs.
This single script chains the legacy full-sample analysis, 100-permutation
enrichment, three 80% internal resamples, the Top 5 unfavorable pathway
context, PubMed retrieval, screened AI synthesis, report, and fixed figure:

```r
source(system.file(
  "examples",
  "example09_tcga_prad_complete_workflow.R",
  package = "InterMeta"
))
```

The three prespecified 80% resamples are the internal consistency analysis for
this PRAD workflow. They summarize repeatability across those three subsets and
must not be interpreted as external validation or as a population-level
uncertainty estimate. The script uses fingerprinted full-sample
checkpoints and the package-native stability/AI resume mechanisms. Paid model
calls and failed-call retries both default to `FALSE`; set their separate
switches to `TRUE` only after inspecting the saved state and deliberately
authorizing the relevant request. API credentials are read only from named
environment variables and are never stored in the script or checkpoints.

The bundled reference enrichment used unadjusted Cox `PV < 0.05` and is
retained for exact continuity with the original analysis. The installed
complete PRAD example selects `analysis_profile = "TimiGP-legacy"` and runs the
same 100-permutation step for the full cohort and all three internal
resamples. These settings are deliberately fixed in this reproduction script;
use the exported InterMeta functions directly when designing a different
analysis profile or permutation scheme.

## One analysis mode

In InterMeta 1.1.4, `TimiGP-legacy` is the only analysis mode. All statistical
functions use it when `analysis_profile` is omitted. The argument remains
available for explicit provenance, but accepts only the exact value
`"TimiGP-legacy"`.

Defaults are `log = TRUE`, gene-wise median centering, the rounded gene-pair
prevalence threshold, logical-factor Cox fitting, legacy annotation permutation,
signed pathway ranking with neutral rows removed, and three unstratified 80%
subsamples with seeds 20260701--20260703. Explicit preprocessing and sampling
parameters remain configurable. In particular, select `log` according to the
input scale; choosing this mode does not establish that another log transform
is appropriate for every cohort. Internal resampling and AI remain available.

For the bundled TCGA-PRAD inputs, this mode yields 861 genes
and 370,230 gene pairs in the full cohort, selects 101,323 pairs at
`PV < 0.05`, returns all 2,652 enrichment rows, and produces 50 scored
pathways. Pathways, order, ranks, groups, and scores agree with the original
TimiGP result to floating-point error of about `5.33e-15`.

The complete PRAD example enforces these dimensions and compares enrichment and
score identities plus numeric columns with the bundled migrated reference at
`tolerance = 1e-12`. It writes
`01_full_sample/tables/parity_checks.tsv` and stops before resampling if the
check fails. `Permutation.FDR` is included in the numeric comparison.

The bundled `Permutation.FDR` column comes from a full rerun with checkpoint
reuse disabled: all 370,230 Cox models were refitted, enrichment and the
100-permutation FDR were regenerated, all 12 full-sample parity checks passed,
and all three prespecified resamples completed. The original enrichment is
retained as a numerical parity reference. The generating package version,
run identifier, package/script/output hashes, and original PCaDB input
fingerprints remain recorded in `InterMetaPRADProvenance`.

The bundled `geneset-all` file contains 2,022 original marker slots and 1,999
unique pathway-gene rows. Ordinary preprocessing and observed enrichment use
the unique table. Historical TimiGP permutation shuffled all 2,022 slots first
and only then formed unique pathway-pair annotations. InterMeta preserves those
raw slots internally for its sole `"TimiGP-legacy"` mode. This ordering is required to reproduce the
original PRAD `Permutation.FDR` values and their 104/50 interactions below
0.2/0.05.

## Input Requirements

The expression object must be a numeric gene-by-sample matrix with unique gene
row names and sample column names. Clinical data must have sample IDs as row
names, event in the first column, and follow-up time in the second column.
Event values must be `0` and `1`.

The only mode defaults to `log = TRUE` and `GMNorm = TRUE`. It applies the
half-sample positivity filter and, when logging, removes rows with any value
`<= -1` before `log10(x + 1)`. Explicitly set `log = FALSE` for inputs whose
documented scale should not be logged again. The bundled synthetic example
uses that explicit setting. Record preprocessing for each cohort.

## Full-Sample Workflow

```r
library(InterMeta)

pathwayset <- InterMetaLoadGeneset("all")
info <- InterMetaCheckEvent(clinical_info)

rna <- InterMetaPreprocess(
  pathwayset = pathwayset,
  rna = bulk_expression,
  cohort = rownames(info),
  log = FALSE,
  GMNorm = TRUE,
  analysis_profile = "TimiGP-legacy"
)

mps <- InterMetaGenePair(
  rna,
  chunk_size = 10000,
  analysis_profile = "TimiGP-legacy"
)
cox <- InterMetaCox(
  mps = mps,
  info = info,
  p.adj = "BH",
  parallel = TRUE,
  core = 8,
  keep_mps = FALSE,
  analysis_profile = "TimiGP-legacy"
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
  cutoff = 0.05,
  analysis_profile = "TimiGP-legacy"
)
```

All analysis calls use `analysis_profile = "TimiGP-legacy"`, including
preprocessing, gene-pair generation, Cox, permutation FDR, and pathway scoring,
and select Cox pairs by `PV < 0.05`. The installed complete PRAD example wires
these choices together to avoid mixing profiles accidentally.

## Optional endpoint-specific validation

### Clinical assessment fields

This entire section describes the optional `purpose = "prognostic_validation"`
workflow. The current default is `biology_exploration`, documented at the top
of this README. The explicit purpose arguments below retain the clinical
workflow; reports use schema 4.3.0 and contain 16 bundle files.

`Literature.Consistency` assesses the cited subset's endpoint-comparable
prognostic direction, not a generic biological support vote. Only direct human
prognostic evidence with a matching endpoint, a verified abstract excerpt, and
an explicitly same/opposite direction contributes. Incidence, preclinical/tissue
studies, reviews and indirect pathway enrichment provide context but cannot
alone validate prognosis. Null associations are never opposing directions.
Unknown endpoints and responses without assessments remain insufficient.
Original model labels/narratives are retained in unendorsed `Model.*` audit
columns; biological context is presented separately. Human review of studies,
exposures, endpoint comparability and cohort overlap remains essential.

Providers return an `assessments` array per pathway: `pmid`, `evidence_type`,
`study_endpoint`, `endpoint_match`, `direction`, `directness`, `abstract_quote`.
The runtime schema supplies allowed enums. Matching a quote verifies grounding,
not semantic correctness. These fields are retained in `assessments.tsv`.

PubMed queries include OXPHOS/lysine-catabolism synonyms and separately labelled
component layers (e.g. taurine/hypotaurine). Components are retrieval expansions,
not strict synonyms of the whole pathway. BCR queries expand biochemical
recurrence/relapse. `include_components = FALSE` disables component layers
without changing ranks or scores.

The evidence layer runs only after the statistical analysis. It does not alter
scores, P-values, pathway ranks, or the existing result objects. No network or
model call occurs unless it is requested explicitly.

Create the restricted reporting context from aggregate enrichment and score
results. The default `direction = "both"` returns at most five pathways in each
direction; use `"unfavorable"` to report only the Top 5 unfavorable pathways:

```r
context <- InterMetaAIContext(
  enrichment = enrich,
  score = score,
  stability = stability,       # optional
  disease = "prostate adenocarcinoma",
  endpoint = "biochemical recurrence-free survival",
  time_unit = "days",
  cohort = "TCGA-PRAD",
  top_n = 5L,
  direction = "unfavorable"
)

context$pathways
context$warnings
```

`InterMetaAIContext()` uses a whitelist. It does not retain expression
matrices, clinical rows, survival times, sample identifiers, or resampling run
payloads. The context contains selected pathway scores, related aggregate
interactions, optional driver gene pairs, scientific warnings, and SHA-256
input fingerprints.

Build reproducible, field-tagged PubMed queries without accessing the network:

```r
queries <- InterMetaBuildPubMedQueries(context, purpose = "prognostic_validation")
queries[, c("pathway", "layer", "query")]
```

Run retrieval explicitly. `retmax = 20L` is a per-query safety bound, and every
query records its status and whether the result was truncated. A failed API
request remains `NA` plus an error status; it is never reported as a genuine
zero-literature result.

```r
evidence <- InterMetaPubMedSearch(
  context,
  queries = queries,
  purpose = "prognostic_validation",
  retmax = 20L,
  email = Sys.getenv("NCBI_EMAIL"),
  cache_dir = "pubmed_cache"
)
```

The optional packages `httr2`, `jsonlite`, and `xml2` are needed for live
retrieval. Set `NCBI_EMAIL` to the email address associated with your NCBI
requests. An NCBI API key can be placed in the `NCBI_API_KEY` environment
variable; the key is read at call time and is not written to results or cache.

A deterministic evidence skeleton does not require a model:

```r
report <- InterMetaAIReport(
  context,
  evidence = evidence,
  purpose = "prognostic_validation",
  language = "zh",
  outdir = "InterMeta_AI_report"
)
```

For model-written interpretation, configure `DEEPSEEK_API_KEY` securely in the
R process environment, then use the built-in provider. Do not put credentials
in scripts, notebooks, reports, or source control. The constructor does not
read the credential or contact the service; the key is read at invocation time.

```r
provider <- InterMetaDeepSeekProvider(model = "deepseek-v4-pro")

report <- InterMetaAIReport(
  context,
  evidence = evidence,
  provider = provider,
  purpose = "prognostic_validation",
  language = "zh",
  max_records_per_pathway = 20L,
  outdir = "InterMeta_AI_report_with_model"
)
```

The provider uses the official `https://api.deepseek.com/chat/completions`
endpoint, not the Anthropic compatibility endpoint. It requests JSON mode,
then validates the parsed response locally; JSON mode alone is not a schema
or biological-correctness guarantee. Defaults are thinking enabled, high
reasoning effort, `max_tokens = 16384L`, and a 180-second timeout. Use
`model = "deepseek-v4-flash"` to select Flash explicitly. A model request may
incur API charges. There is no automatic retry or fallback model: empty,
malformed, incomplete, filtered, or credential-containing responses fail
closed. For non-stop DeepSeek completions, the error exposes only an allowlisted
`finish_reason` (`length`, `content_filter`, `tool_calls`, or
`insufficient_system_resource`; otherwise `unknown`) so a cached-evidence retry
can be chosen deliberately. Review failures before explicitly trying again.

Safe response metadata (requested/returned model, token usage, response ID,
settings, and request/response hashes) is retained in
`report$provenance$provider_call`. Credentials, raw HTTP bodies, and reasoning
text are not retained. The separate public server example can preserve only the
parsed JSON object when a completed provider response later fails local
validation or report publication. That explicitly unvalidated, hash-bound
candidate is not a report and can only be replayed offline against the exact
saved input contract. The generic `InterMetaAIProvider(generate, name, model)`
interface remains available for custom clients returning a parsed R list.

The model may classify literature consistency only as `supported`, `mixed`,
`potentially_contradicted`, or `insufficient_evidence`. These labels describe
the retrieved literature subset, not biological truth. Invented pathways,
unsupported labels, uncited claims, and PMIDs absent from the corresponding
pathway's evidence sent to the model cause validation to fail. The three
non-insufficient labels require a cited claim. Source-membership validation
cannot establish that a paper supports every model assertion: human review is
still required. Article text is treated as untrusted evidence content. Model
output must use structured PMID fields for citations rather than inline identifiers
or URLs in narrative text; uncited insufficient-evidence summaries use a fixed
safe statement. Retrieved records
are selected round-robin across retrieval layers so prognosis results cannot
consume the entire default quota. Abstracts and aggregate tables are bounded,
and truncation is recorded.

Each published directory contains Markdown, RDS and TSV outputs, the exact
queries and evidence, per-query search status, response hashes, warnings,
provenance, a SHA-256 manifest, and a completion marker. Existing output
directories are never overwritten. Every report also
states that results are exploratory, non-causal, and unsuitable for clinical
decisions.

## Directional Pathway Chord Diagram

`InterMetaPathwayChord()` visualizes favorable-to-unfavorable pathway links
directly from an enrichment table. By default it draws every interaction
meeting `condition < cutoff`. The default classic layout keeps sectors
proportional to their link weights, uses fixed labels (`label_cex = 0.7`), and
sets the circular canvas extent to 1.45. Pathway colors are assigned
automatically. The server workflow and example below keep a 14 by 14 inch PDF
page to give long labels room.

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

Dense-layout adjustments remain available by explicitly setting
`adaptive_labels = TRUE`, `canvas_limit = 1.65`, and
`sector_floor_quantile = 0.75`. These options adjust the layout while preserving
the selected interactions, weights, palette, arrows, and analysis results.

Use `select = c(...)` with enrichment `Index` values for an exact edge set.
Use `max_edges = 20` only when an explicitly reduced diagram is desired.
Optional named `group`, `color`, and `order` mappings may include unused
pathways; the function safely crops them to the pathways actually drawn.
Default colors are stable per pathway across different selections; a named
`color` mapping can still be reused when a specific palette is required.

`InterMetaGeneChord()` displays the retained metabolic gene pairs inside one
selected pathway interaction:

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

This requires enrichment made with `keep_gene_pairs = "shared"` (the default)
or `"all"`. Compact migrated objects that deliberately omit `Shared.MGP` cannot
be used for this gene-level chord without recomputing enrichment.
The server script exposes `requested_gene_chord_index`; leave it `NULL` to
auto-select the most significant qualifying interaction with at most
`gene_chord_max_pairs = 50`, or set it to a specific Index to override the cap.

The reference-aligned summary plots use ten pathway interactions and five
pathways in each prognostic direction by default:

```r
dotplot <- InterMetaDotplot(
  enrich,
  condition = "Adjust.P.Value",
  cutoff = 0.05,
  top = 10
)
score_bar <- InterMetaPathwayScoreBar(score, top = 5)
```

`InterMetaPathwayPair()` returns a compact object, and `InterMetaEnrich()` uses
matrix counting for pathway-pair enrichment.

The package workflow defaults to `PV`. Because `PV` is
unadjusted at the gene-pair stage, report this choice explicitly and consider
an explicit `QV` sensitivity analysis for new inference.

## Permutation FDR

```r
enrich <- InterMetaPermFDR(
  resdata = enrich,
  pathwayset = pathwayset,
  gene = query_pairs,
  background = background,
  niter = 100,
  core = 8,
  seed = 1,
  analysis_profile = "TimiGP-legacy"
)
```

Permutation iterations are parallelized at the outer iteration level. The
compact annotation avoids rebuilding millions of text rows in every
permutation. The `TimiGP-legacy` mode uses the original iteration seeds `1:niter`,
rank denominator, raw marker-slot shuffle-before-de-duplication order, and
non-monotone FDR output; the `seed` argument therefore does not
alter this sequence. `InterMetaLoadGeneset()` prints and
returns the unique pathway-gene table; its raw multiplicity is retained as
internal provenance only for compatible legacy permutation and caching.

## Resampling Stability

The user-adjustable defaults are three runs with 80 percent of the cohort.
Three runs implement the package's prespecified internal consistency check.
Report the number of runs, sampling fraction, seeds, and stratification rule,
and do not present this check as external validation or as a population-level
uncertainty estimate. At least two runs are needed to draw an across-run sample
standard deviation; `error = "none"` is the explicit single-run diagnostic mode
of `InterMetaStabilityScoreBar()`. Gene-pair selection defaults
to exploratory `PV < 0.05`; set `pair_condition = "QV"` explicitly and use a
separate output directory for an adjusted sensitivity analysis.

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
  min_censored = 5,
  analysis_profile = "TimiGP-legacy"
)
```

The same workflow can omit `analysis_profile`, `seeds`, and `stratify_event`
because these settings are the public defaults. Set
`run_permutation = TRUE, permutation_niter = 100` when permutation FDR is needed.

This uses `PV < 0.05`, seeds 20260701--20260703, unstratified sampling,
legacy preprocessing, and signed pathway scoring. The complete PRAD example
also enables 100 annotation permutations per analysis. The permutation count
and the three internal resamples are separate design parameters: the former
estimates the pathway-interaction permutation FDR, whereas the latter records
internal repeatability across three prespecified 80% subsets.

The returned object contains:

- `summary`: pathway selection, signed score, legacy rank, magnitude rank,
  conditional direction stability, and overall/direction-specific Top 5, 10,
  and 20 frequencies;
- `pairwise`: pairwise score correlations plus overall, favorable, and
  unfavorable magnitude-Top-K Jaccard indices;
- `edge_summary`: recurrence of directional pathway interactions;
- `score_long`: every pathway score from every successful resample;
- `run_status`: success, failure, event balance, runtime, and result counts;
- `active_run_manifest`: the exact run IDs, seeds, status, and checkpoint paths
  that belong to the current request, so older higher-numbered directories are
  not accidentally included after reducing `n_resample`;
- `top_set_definition`: a machine-readable record that Top-K uses finite
  `abs(Difference)`, an ascending pathway-name tie-break, and exact set size.
  `Stability.Class` uses the largest requested K (Top 20 by default); report
  `Top5.*.Frequency` directly when the claim concerns Top 5. Empty-versus-empty
  Jaccard values are undefined (`NA`), and one successful run is `Not assessed`;
- `provenance`: package version, parameters, session information, and input
  fingerprint.

Each run is checkpointed as `resample_NNN/run_result.rds`. Only a successful
checkpoint whose data, geneset, actual pathway-pair annotation, parameters,
seed, and package-version fingerprints all match is reused. Failed checkpoints
remain available for diagnosis and are retried on the next call. Large per-run
details stay in these files by default, while the final
`stability_result.rds` remains compact.

When `outdir` is supplied, `active_run_manifest.txt` is written before the
first resample with `pending` status and refreshed after every reused or newly
completed round. This makes an interrupted directory self-identifying instead
of leaving a manifest from an earlier request. Manifest refreshes use a
complete same-directory staging file and atomic replacement on Unix-like
servers; the Windows fallback preserves `.previous` and restores it on the
next write if replacement was interrupted. Use one `InterMetaStability()`
writer per `outdir`; concurrent calls must use different output directories.

Plot the score differences across resamples: pathways must appear in every run
with a fully consistent favorable or unfavorable direction; bars show mean
Difference and error bars show its sample standard deviation.

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

The installed server example and the complete PRAD example execute the full
cohort first, write its Cox, enrichment, score, network, and PDF outputs, and
only then start internal subsampling with the same profile, pair-selection, and
pathway-scoring rules. The complete example defaults to `"TimiGP-legacy"`, 80%
unstratified sampling repeated three times from seed 20260701, and 100
permutations per analysis. By default, it requires all requested rounds to
succeed before it draws the stability summary; after a failure, rerunning
reuses completed rounds and retries the failed checkpoint. It invalidates old
aggregate plots at the start and writes `internal_resampling_complete.rds`
only after the current run manifest, tables, plot, and SHA-256 output manifest
are complete.
