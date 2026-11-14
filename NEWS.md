# InterMeta 1.1.4

- TimiGP-legacy is the only analysis mode. Omitted profile arguments use it;
  other values fail with migration guidance. Alternative statistical branches
  have been removed rather than silently mapped to legacy semantics.
- Legacy log, ranking, permutation and three-resample defaults are explicit.
  Preprocessing parameters remain configurable for the documented input scale.
- Re-scores the bundled PRAD enrichment to its 50-row legacy view; synthetic
  reference results are regenerated with the same sole mode. The historical
  1.1.3 full-run/FDR provenance remains intact with an added score migration record.
- Updates examples and help. AI algorithms and historical AI archive metadata
  are unchanged; this change does not address evidence-to-pathway interpretation.

# InterMeta 1.1.3

- Version 1.1.3 adds named maintainer metadata,
  GPL-3 COPYING and TimiGP attribution notices, gene-set provenance records,
  and PCaDB source/FDR-lineage metadata. It removes nine unused upstream image
  copies, adopts PRAD-specific reference filenames with warnings for historical
  names, and corrects the bundled PRAD time unit in examples to days. These
  provenance changes do not alter the statistical implementation or any
  substantive PRAD result; the republished FDR differs from the previous
  bundle only at machine-roundoff scale.
- Defines the default three 80% internal subsamples as a prespecified,
  descriptive internal-consistency check and removes the contradictory
  recommendation and runtime warning that treated the package default as a
  workflow-only setting. Documentation now distinguishes these three samples
  from the 100 annotation permutations used to estimate permutation FDR.
- The complete PRAD example supports environment-controlled output, stage, and
  checkpoint-reuse settings for clean acceptance runs. It records PCaDB data
  provenance, RNG kind, script hash, run timing, and a SHA-256 output manifest.
- A no-checkpoint 1.1.3 PRAD acceptance run refitted all Cox models, recomputed
  the 100-permutation FDR, passed all 12 full-sample parity checks, and
  completed the three prespecified internal resamples. The bundled FDR column
  is republished from that run; the historical enrichment object is used only
  for parity verification. A fail-closed provenance updater verifies every
  output-manifest hash before publishing the data and run fingerprints.
- Adds a Chinese manuscript methods/results template tied to the fixed PRAD
  example parameters. These changes affect reporting and reproducibility, not
  the fitted models or scoring rules.

# InterMeta 1.1.2

- Adds a portability-only maintenance update to the recoverable AI workflow.
  A new workflow records the default PubMed cache as the relative
  `pubmed_cache` path anchored to the run-directory parent, rather than
  serializing the machine-specific absolute path.
- `InterMetaAIResume()` accepts an optional operational `pubmed_cache_dir`
  override when an unfinished PubMed stage must be continued from a relocated
  run. The override is not written into the fixed workflow configuration and
  does not change its hash. Compatible legacy runs retain their original
  configuration; their conventional sibling `pubmed_cache` is mapped after a
  move, or an explicit override can identify a nonstandard cache.
- Synthesis checkpoint references are stored relative to their declared
  output base: direct synthesis results use their `outdir`, while the workflow
  stage artifact uses its own `05_synthesis` directory. Replay metadata records
  its source run relatively, and
  the complete-workflow human summary records the AI run relative to the
  workflow root. On synthesis-checkpoint reuse,
  `Model.Call.Attempted = FALSE` now describes the current invocation, while
  `Paid.May.Have.Occurred` preserves the checkpoint's historical billing-risk
  record.
- Offline replay now refuses a modified workflow configuration before reading
  it for report publication, closing the configuration-SHA verification loop.
- There is no change to statistical analysis, PubMed query definitions,
  retrieval limits, screening or synthesis prompts, model selection, evidence
  validation, report schema, or figure semantics. Completed InterMeta 1.1.1
  results remain valid and can be inspected, resumed when already complete,
  or replayed locally with InterMeta 1.1.2 without repeating PubMed or model
  calls.

# InterMeta 1.1.1

- Fixes a direction bias in resampling summaries under the historical
  `TimiGP-legacy` profile. Stability Top-K membership and pairwise Top-K
  Jaccard values are now defined by finite absolute signed difference with a
  deterministic pathway-name tie-break, so favorable and unfavorable
  pathways receive the same opportunity to enter a Top-K set. Per-run legacy
  scores and ranks remain unchanged.
- Adds Top 5 to the default stability summary, favorable- and
  unfavorable-specific Top-K audit frequencies/Jaccards, magnitude ranks, a
  machine-readable Top-K definition, and typed zero-/one-run results. Empty
  versus empty Top-K overlap is reported as undefined rather than perfect,
  and a single successful resample is explicitly `Not assessed`.
- PubMed query dictionaries use version 3.1 (biology exploration) or 2.1
  (endpoint-oriented retrieval). Lysine-degradation component searches now use
  pathway-specific intermediates instead of the nonspecific bare term
  `lysine`, based on the low-precision pattern observed in the 1.1.0 server
  acceptance run.
- A successful but capped PubMed stage is recorded as
  `complete_with_truncation`, with an explicit workflow-state message. Figure
  checkpoint metadata uses portable run-relative output paths.
- Biology-exploration reports explicitly identify
  `exploration_summary.tsv` and `exploration_findings.tsv` as authoritative.
  Legacy prognostic compatibility tables are marked `not_applicable` instead
  of looking like a failed AI assessment. This additive report metadata and
  clarified table semantics are recorded as report schema 4.3.0. Report schema
  participates in the local report-stage fingerprint, so resuming an intact
  older run can republish the report and figure without another model call.

# InterMeta 1.1.0

- Adds a package-native, resumable two-stage literature workflow. R builds and
  executes deterministic PubMed searches; a model screens articles in bounded
  pathway batches; R verifies identifiers, enums, disease labels and exact
  abstract excerpts; a second per-pathway model pass synthesizes only the
  verified structured evidence.
- The synthesis pass selects frozen `Candidate.ID` values and never recopies an
  abstract or invents citation identifiers. R maps the selected identifiers
  back to the verified evidence and applies the existing quote validator again.
  This removes the repeated quote-copy boundary that caused several 1.0.0
  server-trial failures.
- Adds explicit screening and synthesis checkpoints, compatible-success reuse,
  row-level quarantine, isolated pathway failures, and opt-in retry of a failed
  paid call. Individually invalid rows remain `unresolved`; structurally
  complete batches aggregate to `complete_with_quarantine`, which remains usable
  while quarantined rows are excluded from synthesis. Failed or incomplete
  whole batches are not treated as complete. Credentials remain environment-only
  and are never serialized.
- Preserves raw returned provider responses unchanged as pre-validation
  artifacts, separately from local validation results. A compatible
  validation-revision change can create a new immutable screening or synthesis
  attempt by local revalidation with zero new model calls. Report and evidence
  publication is hash-bound to its source records,
  screening audit and synthesis, and report/figure files use numbered attempt
  directories so stale outputs are not silently reused.
- Adds a package-native fixed exploration figure and an end-to-end workflow,
  status, resume and offline replay interface. The bundled TCGA-PRAD example
  remains public; no external or confidential cohort is included.
- Adds `example09_tcga_prad_complete_workflow.R`, a single
  fingerprinted acceptance script spanning the legacy full-sample analysis,
  100-permutation enrichment, internal resampling, Top 5 unfavorable pathway
  context, PubMed screening, AI synthesis, report and fixed figure. Model calls
  and failed-call retries remain separately opt-in.
- Adds `Biological.Plausibility` as an evidence-bounded exploratory label. It is
  not a pathway truth score, causal conclusion, clinical grade or validation of
  the InterMeta score direction.

# InterMeta 1.0.1

- Maintenance release of the schema 4.2.0 tumor-biology exploration workflow.
  Statistical analysis, PubMed query version 3, model presets, provider-visible
  AI contract `20260905-quote-replay-v1`, and local validation revision
  `20260906-quote-quarantine-v1` are unchanged.
- The versioned RStudio Server entry accepts hash-bound candidates produced by
  either 1.0.0 or 1.0.1 when their contract, validation revision, context,
  evidence and response hashes agree. This permits the preserved trial05
  candidate to be revalidated without another provider request.
- The separately distributed main-figure helper now defaults to an ASCII
  `compact_english` article mode. It retains structured pathways, signed scores,
  biological processes, evidence channels, scope, key PMIDs and grounding
  counts while keeping the full source-linked narrative in the report. Explicit
  `text_mode = "full"` preserves the prior full-text behavior and its strict
  non-ASCII font preflight.

# InterMeta 1.0.0

- Designates the existing 0.10.0 implementation as InterMeta 1.0.0.
- Package version, release filenames and server installation checks are updated.
  Statistical, retrieval and model-selection logic are unchanged; the current
  figure helper additionally displays the new quote and eligibility audit counts.
- A release-candidate server trial exposed an opaque DeepSeek non-stop
  completion. The adapter now reports an allowlisted `finish_reason` without
  retaining raw output. Biology-exploration prompts request at most five
  non-redundant, source-linked findings per pathway to prevent unbounded JSON.
- Later completed responses exposed non-contiguous or clipped Abstract excerpts.
  Report schema 4.2.0 quarantines these finding rows instead of aborting the
  otherwise auditable report. The original model text remains visible, while
  `Quote.Verified = FALSE`, `Quote.Match.Method = source_mismatch`, and
  `Eligible.For.Summary = FALSE` keep it out of finding-derived evidence
  channels, PMID/process aggregates, and source-linked observations. A non-empty
  quote with fewer than 20 Unicode letters or numbers is retained as
  `insufficient_information`; one longer than 1,000 characters is retained as
  `overlength`. None is repaired, paraphrase-matched, or counted as evidence.
  The unchanged provider-visible prompt retains its historical eight-character
  request for replay compatibility; the independent local validator applies
  the stricter 20-character rule.
- Valid exploration quotes still use two auditable matching levels:
  whitespace-exact, followed only by an explicit set of Unicode typography
  equivalents. Case, words, scientific symbols and ordering are never
  fuzzy-matched. Pathway summaries add `Quote.Grounding` with `complete`,
  `partial`, `none`, or `no_findings` so mixed grounding is visible rather than
  hidden by a report-level pass/fail outcome. This field measures strict quote
  matching only. `Evidence.Eligibility` separately reports `complete`,
  `partial`, `none`, or `no_findings`, derived from eligible/all findings after
  strict quote matching and consistency of the model-supplied disease fields
  with the target label (`Disease.Match.Verified`). This is a structural/lexical
  gate, not independent Abstract or MeSH disease verification. With partial eligibility, the
  public biological summary is derived only from eligible
  `Source.Linked.Observations`, public context/hypothesis text is suppressed,
  and the original `Model.*` narratives remain audit-only. The main figure
  separates quote quarantine from target-disease exclusion.
- The provider-visible AI contract revision remains
  `20260905-quote-replay-v1`; a compatible parsed candidate from the earlier
  trial05 run can therefore be replayed without another provider call. The
  independent validation revision is `20260906-quote-quarantine-v1` and is
  checked separately to reject a stale installed or in-memory 1.0.0 build.
  Provenance records both revisions. The complete report bundle remains 16
  files. Previous 0.10.0 reports and figures retain schema 4.0.0 and their
  original version and model/assistant assessment provenance.

# InterMeta 0.10.0

- Default AI purpose is now exploratory tumour biology. Explicit
  `prognostic_validation` retains the endpoint-specific mode and legacy query plan.
- Query version 3 uses broad disease/pathway and biological-process searches;
  components and reviews remain separately tagged. Optional paged retrieval
  records sorting, page size, response hashes, truncation and cache provenance.
- Schema 4 supports multiple findings and evidence types per PMID: functional
  experiments, omics/bioinformatics and patient studies are parallel sources,
  not a universal evidence ladder. Null and contextual findings are retained.
- Reports add exploration_summary.tsv and exploration_findings.tsv (16 files).
  Public cohort labels help identify reuse; unique PMID counts do not imply
  independent replication or quantify biological truth.
- Statistical code/data are unchanged. No credential is bundled; model calls
  remain explicit. Literature interpretation cannot establish pathway activity,
  causal importance or validate the original signed network score.

# InterMeta 0.9.0

- Added pathway-level biological synthesis separately from endpoint-specific
  prognostic consistency. Same-disease experimental pathway mechanisms,
  component or associative context, conflicts, and missing assessments are
  distinguished; a citation alone never upgrades biological interpretability.
- Report schema 3 adds grounded biological findings and cross-paper mechanism
  narratives to the same explicit model request. Biological evidence cannot
  bypass the existing strict clinical endpoint/direction gate.
- Reports add pathway_synthesis.tsv and biology_findings.tsv (14 files total).
  Legacy reports lacking synthesis remain explicitly not assessed. Source
  membership, quotation checks, model provenance, and human-review limits apply.
- Original statistical code and data remain unchanged. A signed pathway score
  is not an assay of pathway activity or metabolic flux.
- Revalidation binds interpretation to the original context and visible evidence
  payload, including selection settings and abstract truncation. Mixed-disease
  labels cannot match a target merely because they mention its organ. Clinical
  null associations do not suppress independent biological mechanism findings.

# InterMeta 0.8.0

- Added InterMetaLLMProvider() for explicit provider and model selection, and
  InterMetaAIModels() for an offline catalog of documented example model IDs.
  DeepSeek remains available through its existing backwards-compatible client.
- Added native OpenAI Responses, Anthropic Messages, and Gemini generateContent
  adapters, plus compatible Qwen and Moonshot chat-completions adapters. Users
  supply their own service credentials at invocation time; no key is bundled.
- All clients retain local schema, citation, endpoint, and direction checks.
  A different model does not relax evidence requirements. Requests have no
  automatic retries, redirects, or cross-provider/model fallback.
- Statistical analysis functions, bundled data, and report schema 2 are unchanged.

# InterMeta 0.7.0

- Query version 2 separates conservative pathway synonyms from component-level
  retrieval and expands BCR to biochemical recurrence/relapse. Unknown endpoint
  labels are no longer literal search terms. Query scope is retained throughout.
- Report schema 2 separates biological context from endpoint-comparable
  directional consistency. Structured cited-study assessments record study type,
  endpoint, directness, direction, and a verified abstract excerpt. Incidence,
  indirect mechanisms, mismatched endpoints, and null associations cannot alone
  support or contradict a prognostic direction. Legacy unassessed outputs are
  conservatively downgraded; original model labels remain auditable.
- Published reports add assessments.tsv (12 files total). These checks are
  evidence-grounding safeguards, not independent biological or clinical validation.
- Statistical analysis functions and bundled data are unchanged.

# InterMeta 0.6.1

- Added a ready-to-use DeepSeek Chat Completions provider with JSON output,
  explicit model and thinking controls, and invocation-time credentials from
  an environment variable. Provider construction never performs a request.
- Reject incomplete, empty, malformed, or credential-containing model output;
  do not automatically retry potentially billable requests or switch models.
- Preserve non-secret model call metadata, token usage, and request/response
  hashes in report provenance without storing raw responses or reasoning text.

# InterMeta 0.6.0

- Added an optional post-analysis evidence-reporting layer. Deterministic R
  code selects the five strongest favorable and/or unfavorable pathways,
  records the analysis profile and statistical limitations, and excludes
  patient-level expression, survival values, and sample identifiers from the
  model context.
- Added reproducible, tiered PubMed query construction and retrieval for
  pathway-disease prognosis, broader pathway-disease biology, reviews, and
  optional driver-gene bridge evidence. Searches retain their exact query,
  PMID set, retrieval status, raw-response hashes, and cache provenance.
- Added a provider-neutral AI interface that accepts structured model output,
  validates every cited PMID against the retrieved evidence, and limits
  literature conclusions to supported, mixed, potentially contradicted, or
  insufficient evidence.
- Added a deterministic Markdown evidence report with mandatory exploratory,
  non-causal, and non-clinical limitations. The AI and PubMed layers are
  optional and cannot change Cox, enrichment, score, stability, or existing
  analysis checkpoints.

# InterMeta 0.5.1

- Changed the `InterMetaStability()` gene-pair selection default from adjusted
  `QV < 0.05` to exploratory unadjusted `PV < 0.05`, while retaining explicit
  `pair_condition = "QV"` for adjusted sensitivity or confirmatory analyses.
  The existing unadjusted-P-value warning remains in place.

# InterMeta 0.5.0

- Added an explicit `analysis_profile` switch to preprocessing, gene-pair
  construction, Cox analysis, permutation FDR, pathway scoring, and internal
  resampling. The package API continues to default to the current
  `"InterMeta"` behavior; `"TimiGP-legacy"` is an opt-in compatibility mode.
- Reproduced the original TimiGP preprocessing without propagating invalid
  logarithms: retain genes positive in at least half of samples, safely remove
  genes containing values less than or equal to -1 before `log10(x + 1)`, and
  median-center each retained gene. Legacy gene-pair prevalence uses the
  original rounded 10-percent threshold, and logical gene pairs use the
  original factor-based Cox model.
- Restored the complete legacy selection and scoring rules: Cox `PV < 0.05`,
  signed `Difference` ranks, removal of neutral pathways, and the original
  permutation-FDR seed and rank-denominator behavior.
- Corrected the legacy permutation annotation order exposed by the first
  server run. The original geneset contains 2,022 marker slots but 1,999
  unique pathway-gene rows; TimiGP shuffled all slots before de-duplication.
  `InterMetaLoadGeneset()` still exposes the 1,999-row unique table for normal
  analysis, while legacy permutation now preserves and shuffles the hidden raw
  2,022 slots. Modern permutation remains unchanged. The release script uses a
  versioned permutation-annotation cache key, so an affected enrichment/full
  checkpoint is recomputed while the compatible Cox checkpoint is retained.
- Added legacy internal-resampling defaults: unstratified 80-percent sampling,
  seeds beginning at 20260701, and deterministic sample ordering. The new
  `InterMeta_0.5.0_server_trial.R` selects this profile by default, repeats the
  sample three times, and runs 100 permutations for the full cohort and every
  resample. Users can still change the profile, fraction, repeat count, seeds,
  and permutation setting at the top of the script.
- With the same TCGA-PRAD expression, prognosis, and geneset inputs, the legacy
  profile restores 861 full-cohort genes, 370,230 gene pairs, and 50 scored
  pathways. The verified pathway table, order, ranks, groups, and scores match
  the original analysis up to floating-point error of about `5.33e-15`.
- Restored classic pathway-chord defaults: proportional sectors, fixed
  `label_cex = 0.7`, a canvas extent of 1.45, and no sector-width floor. Dense
  layout controls remain opt-in with `adaptive_labels = TRUE`,
  `canvas_limit = 1.65`, and `sector_floor_quantile = 0.75`. All qualifying
  interactions are still drawn by default (`max_edges = NULL`); weights,
  pathway colors, arrows, and analysis results are unchanged.
- Kept the server pathway-chord PDF at 14 by 14 inches to leave room for long
  labels while restoring the classic layout. The separate gene-pair chord
  design is unchanged: it accepts an explicit enrichment Index, with automatic
  selection limited to interactions having at most 50 retained gene pairs.
- Hardened tabular exports so embedded tabs and line breaks in plot labels or
  error messages cannot split one record across multiple physical TSV lines;
  server exports are read back and structurally validated before completion.
- Documented why 0.4.0 did not match the historical run: it analyzed 1,260
  genes instead of the 861-gene legacy universe, selected `QV` instead of
  `PV`, used event-stratified seeds 1--3 instead of the legacy unstratified
  seeds, and ranked absolute differences while retaining neutral pathways.
  Because these choices change the Cox universe and sampled patients, an old
  0.4.0 result cannot be made fully equivalent by thresholding or plotting it
  after the fact; it must be recomputed with the compatibility profile.

# InterMeta 0.4.0

- Changed `InterMetaPathwayChord()` to draw all pathway interactions meeting
  the requested significance cutoff by default; edge truncation is now
  opt-in through `max_edges`.
- Added `InterMetaGeneChord()` to expand one selected pathway interaction into
  its retained favorable-to-unfavorable metabolic gene pairs, with
  deterministic side-specific colors and auditable returned data.
- Restyled the pathway enrichment dot plot to match the restored reference and
  changed the score bar default to the five strongest pathways in each
  prognostic direction.
- Added `InterMetaStabilityScoreBar()` for the mean and standard deviation of
  recurrent, direction-consistent pathway scores across internal resamples.
- Added a long-form resample-score output and made compatible early checkpoints
  reusable when the requested number of resamples is increased. The cache
  fingerprint now includes the actual pathway-pair object, and failed run
  checkpoints are diagnosed but retried instead of being reused. An explicit
   active-run manifest prevents stale higher-numbered directories from being
   mistaken for part of a later, smaller request; it is written with pending
   states before computation and refreshed after every round through a complete
   same-directory staging file. Each output directory is a single-writer
   checkpoint location. Resampling
  controls now use strict scalar/integer validation, and the versioned per-run
  permutation-seed scheme is part of the cache key without integer-overflow
  multiplication.
- Reordered the server workflow so full-sample analysis and visual outputs are
  completed before user-adjustable internal resampling, whose defaults remain
  80 percent of samples repeated three times.
- Added a 0.4.0 TCGA-PRAD server script with full-sample Cox and enrichment
  checkpoints, matched full/resample scoring rules, and reference-aligned PDF
  outputs. By default, internal resampling is blocked until a matching
  full-sample checkpoint and all key tables/PDFs are present. Old aggregate
   resampling plots are invalidated at the next request only after an
   in-progress state has been written; failed stale-file deletion stops the
   workflow. A completion marker is written only after the current tables,
   plot, and SHA-256 manifest exist.

# InterMeta 0.3.1

- Added `InterMetaPathwayChord()` for directional metabolic pathway chord
  diagrams using enrichment results.
- Added automatic deterministic pathway colors and optional pathway grouping,
  color, and order mappings without the legacy `dataset` or `geneset`
  requirements.
- Added strict cropping and factor-level cleanup for selected pathways, finite
  link weights (including zero P-values), readable edge limits, invisible
  reproducibility metadata, and graphics-state cleanup.
- Added a TCGA-PRAD chord example, server-trial output, and regression tests for
  the failure modes observed in the legacy `TimiCellChord()` workflow.

# InterMeta 0.3.0

- Added a user-curated TCGA-PRAD cohort containing prognosis data for 463
  samples and normalized expression for 1,260 genes matched to the bundled
  `geneset-all` annotation.
- Added migrated TCGA-PRAD pathway enrichment, 100-iteration permutation FDR,
  pathway scores, and input provenance with SHA-256 fingerprints.
- Verified that compact InterMeta enrichment reproduces all 2,652 legacy
  pathway-interaction counts and P-values.
- Added quick precomputed and full server recomputation examples for the real
  TCGA-PRAD data.

# InterMeta 0.2.1

- Added a reproducible synthetic prognosis cohort with 120 samples and 30
  metabolic genes.
- Added a six-pathway example annotation and precomputed Cox, enrichment,
  permutation-FDR, and pathway-score results.
- Added package-readable `conditions.txt`, `exp_tumor.txt`, and geneset files
  under `inst/extdata/example`.
- Updated the full workflow and stability examples to run directly with the
  bundled data.

# InterMeta 0.2.0

- Added strict expression and survival-data validation.
- Changed the default expression setting to `log = FALSE`.
- Added memory-bounded gene-pair generation.
- Replaced expanded pathway-pair tables with compact pathway membership data.
- Reimplemented enrichment with matrix-based counting.
- Parallelized permutation iterations and made their seeds configurable.
- Added event-stratified repeated subsampling.
- Added input- and parameter-validated checkpoint reuse.
- Added per-run error isolation and status reporting.
- Corrected direction consistency to be conditional on pathway selection.
- Added zero-frequency rows for pathways never selected.
- Added directional pathway-interaction recurrence summaries.
- Preserved neutral pathway scores instead of dropping them.
- Ranked pathway importance by absolute score difference while retaining a
  separate signed rank.
- Moved runnable examples under `inst/examples`.
- Removed legacy immune datasets and code from the installed package.
- Added automated tests and standard package metadata.
