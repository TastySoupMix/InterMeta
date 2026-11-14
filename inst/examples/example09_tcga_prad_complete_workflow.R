# End-to-end InterMeta 1.1.3 workflow using only the bundled public
# TCGA-PRAD expression, BCR, and reference-result objects.
#
# This example runs both analysis blocks from the original InterMeta/TimiGP
# workflow and then hands the newly computed aggregate results to the
# recoverable PubMed/AI workflow:
#
#   all 463 samples -> Cox -> enrichment -> 100 permutations -> pathway score
#   repeated 80% subsamples -> stability summary (100 permutations per run)
#   current enrichment + score + stability -> Top 5 unfavorable pathways
#   PubMed -> bounded AI screening -> AI synthesis -> report and fixed figure
#
# Cox, permutation, and resampling stages are computationally expensive. The
# three prespecified 80% resamples are the internal consistency analysis used
# by this PRAD example. They describe repeatability across those three subsets;
# they are not external validation or a population-level uncertainty estimate.
# Completed, fingerprinted checkpoints are reused when this script is sourced
# again.
#
# No API key belongs in this file. Configure DEEPSEEK_API_KEY outside R and
# change allow_model_calls to TRUE only when you deliberately authorize paid
# requests. With its default FALSE value, the workflow cannot call a model.

library(InterMeta)

## 0. User settings -----------------------------------------------------------

output_dir <- file.path(getwd(), "InterMeta_TCGA_PRAD_complete_1_1_3")

# These switches allow an interrupted workflow to be continued by stage.
run_full_sample <- TRUE
run_stability <- TRUE
run_ai_workflow <- TRUE
reuse_checkpoints <- TRUE

# A validation runner can override only these operational settings without
# editing the scientific parameters below. Accepted flag values are
# true/false, yes/no, and 1/0 (case-insensitive).
read_flag_override <- function(name, default) {
  value <- trimws(tolower(Sys.getenv(name, unset = "")))
  if (!nzchar(value)) {
    return(default)
  }
  if (value %in% c("true", "yes", "1")) {
    return(TRUE)
  }
  if (value %in% c("false", "no", "0")) {
    return(FALSE)
  }
  stop(name, " must be one of true/false, yes/no, or 1/0.")
}
output_dir_override <- Sys.getenv("INTERMETA_OUTPUT_DIR", unset = "")
if (nzchar(output_dir_override)) {
  output_dir <- output_dir_override
}
run_full_sample <- read_flag_override(
  "INTERMETA_RUN_FULL_SAMPLE", run_full_sample
)
run_stability <- read_flag_override(
  "INTERMETA_RUN_STABILITY", run_stability
)
run_ai_workflow <- read_flag_override(
  "INTERMETA_RUN_AI_WORKFLOW", run_ai_workflow
)
reuse_checkpoints <- read_flag_override(
  "INTERMETA_REUSE_CHECKPOINTS", reuse_checkpoints
)

# Three runs are the prespecified internal consistency design for this example.
n_resample <- 3L
sample_fraction <- 0.8

# This example deliberately reproduces the legacy TimiGP analysis contract.
analysis_profile <- "TimiGP-legacy"
pair_condition <- "PV"
pair_cutoff <- 0.05
score_condition <- "Adjust.P.Value"
score_cutoff <- 0.05
permutation_niter <- 100L

# Paid model requests require this explicit safety gate. Leave FALSE to run
# only stages that require no model request. Setting TRUE does not authorize a
# retry after a failed or delivery-uncertain request unless retry_failed is
# also deliberately changed to TRUE.
allow_model_calls <- FALSE
retry_failed_model_calls <- FALSE
ai_model <- "deepseek-v4-pro"
ai_key_env <- "DEEPSEEK_API_KEY"

if (!identical(analysis_profile, "TimiGP-legacy")) {
  stop("This public reproduction example requires TimiGP-legacy.")
}
if (!identical(pair_condition, "PV") || !identical(pair_cutoff, 0.05)) {
  stop("This reproduction example requires pair_condition = 'PV' and 0.05.")
}
if (!identical(score_condition, "Adjust.P.Value") ||
    !identical(score_cutoff, 0.05)) {
  stop(
    "This reproduction example requires Adjust.P.Value scoring at 0.05."
  )
}
if (!identical(permutation_niter, 100L)) {
  stop("This complete reproduction example requires 100 permutations.")
}
if (!is.numeric(n_resample) || length(n_resample) != 1L ||
    !is.finite(n_resample) || n_resample < 2L ||
    n_resample != floor(n_resample)) {
  stop("n_resample must be one integer of at least 2.")
}
n_resample <- as.integer(n_resample)
resample_seeds <- 20260701L + seq_len(n_resample) - 1L
if (!is.numeric(sample_fraction) || length(sample_fraction) != 1L ||
    !is.finite(sample_fraction) || sample_fraction <= 0 ||
    sample_fraction > 1) {
  stop("sample_fraction must be one finite number greater than 0 and <= 1.")
}
logical_settings <- list(
  run_full_sample = run_full_sample,
  run_stability = run_stability,
  run_ai_workflow = run_ai_workflow,
  reuse_checkpoints = reuse_checkpoints,
  allow_model_calls = allow_model_calls,
  retry_failed_model_calls = retry_failed_model_calls
)
valid_logical_setting <- vapply(
  logical_settings,
  function(value) is.logical(value) && length(value) == 1L && !is.na(value),
  logical(1)
)
if (!all(valid_logical_setting)) {
  stop(
    "These settings must each be exactly TRUE or FALSE: ",
    paste(names(logical_settings)[!valid_logical_setting], collapse = ", ")
  )
}
if (allow_model_calls &&
    !nzchar(Sys.getenv(ai_key_env, unset = ""))) {
  stop(
    ai_key_env, " is not configured. Set it outside this script, restart R, ",
    "and then explicitly set allow_model_calls <- TRUE."
  )
}

available_cores <- parallel::detectCores(logical = FALSE)
if (is.na(available_cores)) {
  available_cores <- 2L
}
core <- max(1L, min(8L, available_cores - 1L))
use_parallel_cox <- core > 1L
RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)
rng_kind <- RNGkind()
run_started_at <- Sys.time()
pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

metadata_dir <- file.path(output_dir, "00_metadata")
full_sample_dir <- file.path(output_dir, "01_full_sample")
full_table_dir <- file.path(full_sample_dir, "tables")
full_figure_dir <- file.path(full_sample_dir, "figures")
stability_dir <- file.path(output_dir, "02_internal_resampling")
stability_figure_dir <- file.path(stability_dir, "figures")
ai_root_dir <- file.path(output_dir, "03_ai_workflow")

for (path in c(
  metadata_dir, full_table_dir, full_figure_dir,
  stability_dir, stability_figure_dir, ai_root_dir
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The digest dependency is required for checkpoint fingerprints.")
}

# Checkpoints are written to a same-directory temporary file before rename.
# An incompatible or unreadable old checkpoint is preserved with a suffix
# instead of being silently accepted as a completed stage.
atomic_save_checkpoint <- function(value, file, stage, fingerprint) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  payload <- list(
    schema = "InterMeta-example09-checkpoint-v1",
    stage = stage,
    fingerprint = fingerprint,
    value_sha256 = digest::digest(
      value, algo = "sha256", serialize = TRUE
    ),
    saved_at = Sys.time(),
    value = value
  )
  temporary <- tempfile(
    pattern = paste0(".", basename(file), "-"),
    tmpdir = dirname(file)
  )
  on.exit({
    if (file.exists(temporary)) {
      unlink(temporary)
    }
  }, add = TRUE)
  saveRDS(payload, temporary)

  if (file.exists(file)) {
    suffix <- format(Sys.time(), "%Y%m%d-%H%M%S")
    preserved <- paste0(file, ".superseded-", suffix)
    if (!file.rename(file, preserved)) {
      stop("Could not preserve the previous checkpoint: ", file)
    }
  }
  if (!file.rename(temporary, file)) {
    stop("Could not publish the completed checkpoint: ", file)
  }
  invisible(file)
}

load_checkpoint <- function(file, stage, fingerprint) {
  if (!reuse_checkpoints || !file.exists(file)) {
    return(NULL)
  }
  payload <- tryCatch(readRDS(file), error = function(error) error)
  valid <- !inherits(payload, "error") && is.list(payload) &&
    identical(payload$schema, "InterMeta-example09-checkpoint-v1") &&
    identical(payload$stage, stage) &&
    identical(payload$fingerprint, fingerprint) &&
    "value" %in% names(payload) &&
    is.character(payload$value_sha256) &&
    length(payload$value_sha256) == 1L &&
    identical(
      payload$value_sha256,
      digest::digest(payload$value, algo = "sha256", serialize = TRUE)
    )
  if (!valid) {
    warning(
      "Ignoring an unreadable or incompatible checkpoint: ", file,
      immediate. = TRUE
    )
    return(NULL)
  }
  message("Reusing validated ", stage, " checkpoint: ", file)
  payload$value
}

write_tsv <- function(x, file) {
  utils::write.table(
    x,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

## 1. Load and fingerprint public inputs --------------------------------------

data("InterMetaPRADInfo")
data("InterMetaPRADRNA")
data("InterMetaPRADEnrich")
data("InterMetaPRADScore")
data("InterMetaPRADProvenance")

if (utils::packageVersion("InterMeta") < "1.1.3") {
  stop("InterMeta 1.1.3 or newer is required.")
}
if (!identical(colnames(InterMetaPRADRNA), row.names(InterMetaPRADInfo))) {
  stop("The bundled TCGA-PRAD expression and clinical samples are misaligned.")
}

pathwayset <- InterMetaLoadGeneset("all")
pathway_pair <- InterMetaPathwayPair(pathwayset)

analysis_parameters <- list(
  analysis_profile = analysis_profile,
  pair_condition = pair_condition,
  pair_cutoff = pair_cutoff,
  score_condition = score_condition,
  score_cutoff = score_cutoff,
  permutation_niter = permutation_niter,
  log = TRUE,
  GMNorm = TRUE
)
full_fingerprint <- digest::digest(
  list(
    schema = "InterMeta-example09-full-v1",
    package_version = as.character(utils::packageVersion("InterMeta")),
    clinical = InterMetaPRADInfo,
    expression = InterMetaPRADRNA,
    pathwayset = pathwayset,
    parameters = analysis_parameters
  ),
  algo = "sha256",
  serialize = TRUE
)
analysis_id <- substr(full_fingerprint, 1L, 16L)
full_checkpoint_dir <- file.path(
  full_sample_dir, "checkpoints", paste0("analysis-", analysis_id)
)
dir.create(full_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(
  capture.output(utils::sessionInfo()),
  file.path(metadata_dir, "sessionInfo.txt")
)
workflow_script <- system.file(
  "examples", "example09_tcga_prad_complete_workflow.R",
  package = "InterMeta"
)
workflow_script_sha256 <- if (nzchar(workflow_script) &&
                              file.exists(workflow_script)) {
  digest::digest(workflow_script, algo = "sha256", file = TRUE)
} else {
  NA_character_
}
writeLines(
  capture.output(dput(InterMetaPRADProvenance)),
  file.path(metadata_dir, "prad_data_provenance.txt")
)
writeLines(
  c(
    paste0(
      "Run started (UTC): ",
      format(run_started_at, tz = "UTC", usetz = TRUE)
    ),
    paste0("InterMeta version: ", utils::packageVersion("InterMeta")),
    paste0("Workflow script SHA-256: ", workflow_script_sha256),
    paste0("PRAD source database: ", InterMetaPRADProvenance$source$database),
    paste0("PRAD source URL: ", InterMetaPRADProvenance$source$url),
    paste0("Analysis fingerprint: ", full_fingerprint),
    paste0("Analysis profile: ", analysis_profile),
    paste0("Cox selection: ", pair_condition, " < ", pair_cutoff),
    paste0("Pathway scoring: ", score_condition, " < ", score_cutoff),
    paste0("Permutations per analysis: ", permutation_niter),
    paste0("Requested resamples: ", n_resample),
    paste0("Sample fraction: ", sample_fraction),
    paste0("Resample seeds: ", paste(resample_seeds, collapse = ",")),
    paste0("RNG kind: ", paste(rng_kind, collapse = " | ")),
    "Event-stratified resampling: FALSE",
    paste0("Worker cores: ", core),
    paste0("Reuse checkpoints: ", reuse_checkpoints),
    paste0("Run full sample: ", run_full_sample),
    paste0("Run stability: ", run_stability),
    paste0("Run AI workflow: ", run_ai_workflow),
    paste0("Model calls authorized: ", allow_model_calls)
  ),
  file.path(metadata_dir, "run_settings.txt")
)

## 2. Full-sample TimiGP-legacy analysis -------------------------------------

full_result_file <- file.path(full_checkpoint_dir, "full_result.rds")
cox_file <- file.path(full_checkpoint_dir, "cox_result.rds")
enrich_file <- file.path(full_checkpoint_dir, "enrichment_result.rds")
full_result <- load_checkpoint(
  full_result_file, "full_sample", full_fingerprint
)
if (!is.null(full_result)) {
  required_full_fields <- c(
    "enrich", "score", "parameters", "analysis_fingerprint",
    "preprocessed_dimension", "cox_gene_pairs"
  )
  valid_full_result <- is.list(full_result) &&
    all(required_full_fields %in% names(full_result)) &&
    is.data.frame(full_result$enrich) &&
    is.data.frame(full_result$score) &&
    identical(full_result$parameters, analysis_parameters) &&
    identical(full_result$analysis_fingerprint, full_fingerprint) &&
    is.numeric(full_result$preprocessed_dimension) &&
    length(full_result$preprocessed_dimension) == 2L &&
    identical(
      as.integer(full_result$preprocessed_dimension), c(861L, 463L)
    ) &&
    identical(as.integer(full_result$cox_gene_pairs), 370230L)
  if (!valid_full_result) {
    warning(
      "Ignoring a structurally incomplete full-sample checkpoint.",
      immediate. = TRUE
    )
    full_result <- NULL
  }
}

cox_res <- load_checkpoint(cox_file, "cox", full_fingerprint)
if (!is.null(cox_res)) {
  valid_cox_result <- is.data.frame(cox_res) &&
    nrow(cox_res) == 370230L &&
    all(c("PV", "QV") %in% names(cox_res)) &&
    !is.null(row.names(cox_res)) &&
    !anyDuplicated(row.names(cox_res))
  if (!valid_cox_result) {
    warning(
      "Ignoring a structurally incomplete Cox checkpoint.",
      immediate. = TRUE
    )
    cox_res <- NULL
  }
}

need_full_computation <- is.null(full_result) || is.null(cox_res)
if (need_full_computation && !run_full_sample) {
  stop(
    "run_full_sample is FALSE and the full-sample/Cox checkpoint set is ",
    "incomplete. Set run_full_sample <- TRUE to reconstruct only the ",
    "missing or invalid stages."
  )
}

if (need_full_computation) {
  message("Full sample: validating and preprocessing public TCGA-PRAD")
  info <- InterMetaCheckEvent(
    InterMetaPRADInfo,
    min_events = 5,
    min_censored = 5
  )
  rna <- InterMetaPreprocess(
    pathwayset = pathwayset,
    rna = InterMetaPRADRNA,
    cohort = row.names(info),
    log = TRUE,
    GMNorm = TRUE,
    analysis_profile = analysis_profile
  )
  if (!identical(dim(rna), c(861L, 463L))) {
    stop("Legacy preprocessing did not produce the expected 861 x 463 matrix.")
  }

  if (is.null(cox_res)) {
    message("Full sample: generating 370,230 pairs and fitting Cox models")
    mps <- InterMetaGenePair(
      rna,
      chunk_size = 10000,
      analysis_profile = analysis_profile
    )
    if (!identical(dim(mps), c(370230L, 463L))) {
      stop("Legacy gene-pair generation did not produce 370,230 x 463.")
    }
    cox <- InterMetaCox(
      mps = mps,
      info = info,
      p.adj = "BH",
      parallel = use_parallel_cox,
      core = core,
      keep_mps = FALSE,
      analysis_profile = analysis_profile
    )
    cox_res <- cox$cox_res
    rm(mps, cox)
    invisible(gc())
    atomic_save_checkpoint(cox_res, cox_file, "cox", full_fingerprint)
  }
  if (nrow(cox_res) != 370230L) {
    stop("The Cox checkpoint does not contain 370,230 gene pairs.")
  }
}

if (is.null(full_result)) {
  enrich <- load_checkpoint(enrich_file, "enrichment", full_fingerprint)
  if (!is.null(enrich)) {
    valid_enrichment <- is.data.frame(enrich) && nrow(enrich) == 2652L &&
      all(c(
        "Pathway.Interaction", "Favorable.Pathway", "Unfavorable.Pathway",
        "No.Total.MGP", "No.Shared.MGP", "P.Value", "Adjust.P.Value",
        "Permutation.FDR", "Shared.MGP"
      ) %in% names(enrich))
    if (!valid_enrichment) {
      warning(
        "Ignoring a structurally incomplete enrichment checkpoint.",
        immediate. = TRUE
      )
      enrich <- NULL
    }
  }
  if (is.null(enrich)) {
    message("Full sample: enrichment and 100 annotation permutations")
    query_pairs <- row.names(cox_res)[
      is.finite(cox_res[[pair_condition]]) &
        cox_res[[pair_condition]] < pair_cutoff
    ]
    background <- InterMetaBG(row.names(cox_res))
    enrich <- InterMetaEnrich(
      gene = query_pairs,
      background = background,
      pathway.pair = pathway_pair,
      p.adj = "BH",
      keep_gene_pairs = "shared"
    )
    enrich <- InterMetaPermFDR(
      resdata = enrich,
      pathwayset = pathwayset,
      gene = query_pairs,
      background = background,
      niter = permutation_niter,
      core = core,
      seed = 1L,
      analysis_profile = analysis_profile
    )
    atomic_save_checkpoint(
      enrich, enrich_file, "enrichment", full_fingerprint
    )
  }

  message("Full sample: calculating signed pathway scores")
  score <- InterMetaPathwayScore(
    enrich,
    condition = score_condition,
    cutoff = score_cutoff,
    analysis_profile = analysis_profile
  )
  full_result <- list(
    enrich = enrich,
    score = score,
    parameters = analysis_parameters,
    analysis_fingerprint = full_fingerprint,
    preprocessed_dimension = dim(rna),
    cox_gene_pairs = nrow(cox_res)
  )
  atomic_save_checkpoint(
    full_result, full_result_file, "full_sample", full_fingerprint
  )
  message("Full sample checkpoint completed")
}

enrich <- full_result$enrich
score <- full_result$score
if (is.null(cox_res) || nrow(cox_res) != 370230L) {
  stop(
    "The validated full-sample result lacks its compatible Cox checkpoint. ",
    "Set run_full_sample <- TRUE to reconstruct it."
  )
}

## 3. Structural and reference-parity acceptance checks ----------------------

required_enrich_columns <- c(
  "Index", "Pathway.Interaction", "Favorable.Pathway",
  "Unfavorable.Pathway", "No.Total.MGP", "No.Shared.MGP",
  "P.Value", "Adjust.P.Value", "Permutation.FDR"
)
required_score_columns <- c(
  "Index", "Rank", "Signed.Rank", "Pathway", "Favorable.Score",
  "Unfavorable.Score", "Difference", "Group"
)

reference_match <- match(
  enrich$Pathway.Interaction,
  InterMetaPRADEnrich$Pathway.Interaction
)
same_reference_keys <- nrow(enrich) == nrow(InterMetaPRADEnrich) &&
  !anyNA(reference_match) && !anyDuplicated(enrich$Pathway.Interaction)

reference_counts_match <- same_reference_keys && identical(
  as.integer(enrich$No.Total.MGP),
  as.integer(InterMetaPRADEnrich$No.Total.MGP[reference_match])
) && identical(
  as.integer(enrich$No.Shared.MGP),
  as.integer(InterMetaPRADEnrich$No.Shared.MGP[reference_match])
)
reference_p_match <- same_reference_keys && isTRUE(all.equal(
  enrich$P.Value,
  InterMetaPRADEnrich$P.Value[reference_match],
  tolerance = 1e-12
)) && isTRUE(all.equal(
  enrich$Adjust.P.Value,
  InterMetaPRADEnrich$Adjust.P.Value[reference_match],
  tolerance = 1e-12
))
reference_permutation_match <- same_reference_keys && isTRUE(all.equal(
  enrich$Permutation.FDR,
  InterMetaPRADEnrich$Permutation.FDR[reference_match],
  tolerance = 1e-12
))

reference_legacy_score <- InterMetaPathwayScore(
  InterMetaPRADEnrich,
  condition = score_condition,
  cutoff = score_cutoff,
  analysis_profile = analysis_profile
)
score_match <- match(score$Pathway, reference_legacy_score$Pathway)
same_score_keys <- nrow(score) == nrow(reference_legacy_score) &&
  !anyNA(score_match) && !anyDuplicated(score$Pathway)
score_parity <- same_score_keys && identical(
  as.character(score$Group),
  as.character(reference_legacy_score$Group[score_match])
) && isTRUE(all.equal(
  score$Difference,
  reference_legacy_score$Difference[score_match],
  tolerance = 1e-12
))

# The bundled score uses the same sole analysis mode from 1.1.4.
bundled_match <- match(score$Pathway, InterMetaPRADScore$Pathway)
bundled_legacy_score_parity <- nrow(score) == nrow(InterMetaPRADScore) &&
  !anyNA(bundled_match) && !anyDuplicated(score$Pathway) &&
  identical(as.character(score$Group),
            as.character(InterMetaPRADScore$Group[bundled_match])) &&
  isTRUE(all.equal(score$Difference,
                  InterMetaPRADScore$Difference[bundled_match], tolerance = 1e-12))

parity_checks <- data.frame(
  Check = c(
    "legacy_preprocessed_dimension_861x463",
    "legacy_cox_pair_count_370230",
    "enrichment_required_columns",
    "enrichment_row_count_2652",
    "enrichment_reference_keys",
    "enrichment_reference_counts",
    "enrichment_reference_p_values",
    "enrichment_reference_permutation_fdr",
    "score_required_columns",
    "legacy_score_row_count_50",
    "legacy_score_reference_parity",
    "bundled_legacy_score_reference_parity"
  ),
  Observed = c(
    paste(full_result$preprocessed_dimension, collapse = "x"),
    as.character(full_result$cox_gene_pairs),
    paste(intersect(required_enrich_columns, names(enrich)), collapse = ","),
    as.character(nrow(enrich)),
    as.character(same_reference_keys),
    as.character(reference_counts_match),
    as.character(reference_p_match),
    as.character(reference_permutation_match),
    paste(intersect(required_score_columns, names(score)), collapse = ","),
    as.character(nrow(score)),
    as.character(score_parity),
    as.character(bundled_legacy_score_parity)
  ),
  Expected = c(
    "861x463", "370230",
    paste(required_enrich_columns, collapse = ","),
    "2652", "TRUE", "TRUE", "TRUE", "TRUE",
    paste(required_score_columns, collapse = ","),
    "50", "TRUE", "TRUE"
  ),
  Pass = c(
    identical(as.integer(full_result$preprocessed_dimension), c(861L, 463L)),
    identical(as.integer(full_result$cox_gene_pairs), 370230L),
    all(required_enrich_columns %in% names(enrich)),
    nrow(enrich) == 2652L,
    same_reference_keys,
    reference_counts_match,
    reference_p_match,
    reference_permutation_match,
    all(required_score_columns %in% names(score)),
    nrow(score) == 50L,
    score_parity,
    bundled_legacy_score_parity
  ),
  stringsAsFactors = FALSE
)
write_tsv(parity_checks, file.path(full_table_dir, "parity_checks.tsv"))
if (!all(parity_checks$Pass)) {
  stop(
    "Full-sample acceptance checks failed. Review ",
    file.path(full_table_dir, "parity_checks.tsv"),
    " before running resampling or AI stages."
  )
}

write_tsv(enrich, file.path(full_table_dir, "enrichment.tsv"))
write_tsv(score, file.path(full_table_dir, "pathway_score.tsv"))
saveRDS(full_result, file.path(full_table_dir, "full_sample_result.rds"))
saveRDS(cox_res, file.path(full_table_dir, "cox_result.rds"))
write_tsv(
  cbind(Metabolic.Gene.Pair = row.names(cox_res), cox_res),
  file.path(full_table_dir, "cox_result.tsv")
)
query_pairs <- row.names(cox_res)[
  is.finite(cox_res[[pair_condition]]) &
    cox_res[[pair_condition]] < pair_cutoff
]
write_tsv(
  data.frame(Metabolic.Gene.Pair = query_pairs),
  file.path(full_table_dir, "selected_gene_pairs.tsv")
)
saveRDS(pathway_pair, file.path(full_table_dir, "pathway_pair.rds"))

InterMetaPathwayNetwork(
  resdata = enrich,
  pathwayset = pathwayset,
  condition = score_condition,
  cutoff = score_cutoff,
  export = TRUE,
  path = full_table_dir
)
dotplot <- InterMetaDotplot(
  enrich,
  condition = score_condition,
  cutoff = score_cutoff,
  top = 10L
)
score_bar <- InterMetaPathwayScoreBar(score, top = 5L)
ggplot2::ggsave(
  file.path(full_figure_dir, "pathway_enrichment_dotplot.pdf"),
  dotplot,
  width = 8,
  height = 7,
  device = pdf_device
)
ggplot2::ggsave(
  file.path(full_figure_dir, "pathway_score_bar.pdf"),
  score_bar,
  width = 8,
  height = 7,
  device = pdf_device
)

chord_file <- file.path(full_figure_dir, "pathway_interaction_chord.pdf")
grDevices::pdf(chord_file, width = 14, height = 14, useDingbats = FALSE)
tryCatch(
  InterMetaPathwayChord(
    enrich,
    condition = score_condition,
    cutoff = score_cutoff,
    max_edges = NULL,
    label_cex = 0.7,
    adaptive_labels = FALSE,
    canvas_limit = 1.45,
    sector_floor_quantile = NULL
  ),
  finally = grDevices::dev.off()
)

gene_chord_candidate <- which(
  is.finite(enrich[[score_condition]]) &
    enrich[[score_condition]] < score_cutoff &
    !is.na(enrich$Shared.MGP) & nzchar(enrich$Shared.MGP) &
    is.finite(enrich$No.Shared.MGP) & enrich$No.Shared.MGP > 0 &
    enrich$No.Shared.MGP <= 50
)
if (length(gene_chord_candidate) > 0L) {
  gene_chord_candidate <- gene_chord_candidate[order(
    enrich[[score_condition]][gene_chord_candidate],
    -enrich$No.Shared.MGP[gene_chord_candidate],
    enrich$Index[gene_chord_candidate]
  )]
  gene_chord_file <- file.path(
    full_figure_dir, "selected_pathway_gene_pair_chord.pdf"
  )
  grDevices::pdf(gene_chord_file, width = 9, height = 9)
  tryCatch(
    InterMetaGeneChord(
      enrich,
      select = enrich$Index[gene_chord_candidate[1L]],
      condition = score_condition
    ),
    finally = grDevices::dev.off()
  )
}

## 4. Repeated internal subsampling ------------------------------------------

if (run_stability) {
  if (reuse_checkpoints) {
    message("Starting internal resampling; compatible per-run caches are reused")
  } else {
    message("Starting a fresh internal resampling run; cache reuse is disabled")
  }
  stability <- InterMetaStability(
    rna = InterMetaPRADRNA,
    info = InterMetaPRADInfo,
    pathwayset = pathwayset,
    pathway.pair = pathway_pair,
    n_resample = n_resample,
    sample_fraction = sample_fraction,
    seeds = resample_seeds,
    p.adj = "BH",
    pair_condition = pair_condition,
    pair_cutoff = pair_cutoff,
    enrich_p.adj = "BH",
    score_condition = score_condition,
    score_cutoff = score_cutoff,
    run_permutation = TRUE,
    permutation_niter = permutation_niter,
    run_network = FALSE,
    log = TRUE,
    GMNorm = TRUE,
    parallel = use_parallel_cox,
    core = core,
    outdir = stability_dir,
    reuse = reuse_checkpoints,
    stratify_event = FALSE,
    min_events = 5,
    min_censored = 5,
    continue_on_error = TRUE,
    keep_details = FALSE,
    analysis_profile = analysis_profile
  )
} else {
  stability_file <- file.path(stability_dir, "stability_result.rds")
  if (!file.exists(stability_file)) {
    stop(
      "run_stability is FALSE and stability_result.rds is absent. ",
      "Set run_stability <- TRUE."
    )
  }
  stability <- readRDS(stability_file)
}

if (!inherits(stability, "InterMetaStability")) {
  stop("The stability result is not a valid InterMetaStability object.")
}
expected_stability_parameters <- list(
  analysis_profile = analysis_profile,
  n_resample = as.integer(n_resample),
  sample_fraction = sample_fraction,
  seeds = as.integer(resample_seeds),
  pair_condition = pair_condition,
  pair_cutoff = pair_cutoff,
  score_condition = score_condition,
  score_cutoff = score_cutoff,
  run_permutation = TRUE,
  permutation_niter = permutation_niter,
  log = TRUE,
  GMNorm = TRUE,
  stratify_event = FALSE
)
for (parameter_name in names(expected_stability_parameters)) {
  if (!identical(
    stability$parameters[[parameter_name]],
    expected_stability_parameters[[parameter_name]]
  )) {
    stop("The stability checkpoint has incompatible ", parameter_name, ".")
  }
}

successful_runs <- sum(stability$run_status$Status == "success")
if (successful_runs != n_resample) {
  stop(
    successful_runs, " of ", n_resample,
    " resamples succeeded. Review 02_internal_resampling/run_status.txt, ",
    "then rerun; successful checkpoints will be reused."
  )
}
stability_plot <- InterMetaStabilityScoreBar(
  stability,
  top = 5L,
  min_selection_frequency = 1,
  min_direction_consistency = 1,
  error = "sd"
)
ggplot2::ggsave(
  file.path(stability_figure_dir, "stability_pathway_score_bar.pdf"),
  stability_plot,
  width = 8.5,
  height = 5.5,
  device = pdf_device
)
stability_plot_data <- attr(stability_plot, "plot_data")
if (!"Display.Pathway" %in% names(stability_plot_data)) {
  stability_plot_data$Display.Pathway <- as.character(
    stability_plot_data$Pathway
  )
}
stability_plot_data_export <- stability_plot_data
text_columns <- vapply(
  stability_plot_data_export,
  function(column) is.character(column) || is.factor(column),
  logical(1)
)
stability_plot_data_export[text_columns] <- lapply(
  stability_plot_data_export[text_columns],
  function(column) gsub("[\t\r\n]+", " ", as.character(column))
)
stability_plot_data_file <- file.path(
  stability_dir, "stability_score_plot_data.tsv"
)
write_tsv(
  stability_plot_data_export,
  stability_plot_data_file
)
stability_plot_data_check <- utils::read.delim(
  stability_plot_data_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(dim(stability_plot_data_check),
               dim(stability_plot_data_export)) ||
    !identical(names(stability_plot_data_check),
               names(stability_plot_data_export)) ||
    !identical(as.character(stability_plot_data_check$Pathway),
               as.character(stability_plot_data_export$Pathway)) ||
    !identical(as.character(stability_plot_data_check$Display.Pathway),
               as.character(stability_plot_data_export$Display.Pathway))) {
  stop("The stability plot-data TSV failed its round-trip validation.")
}

## 5. Build the Top 5 aggregate AI context ----------------------------------

# Passing all three current aggregate results through analysis ensures the AI
# context records the exact thresholds/profile while excluding patient rows,
# expression matrices, sample IDs, and resampling payloads.
analysis_for_ai <- list(
  enrich = enrich,
  score = score,
  stability = stability,
  parameters = list(
    analysis_profile = analysis_profile,
    pair_condition = pair_condition,
    pair_cutoff = pair_cutoff,
    score_condition = score_condition,
    score_cutoff = score_cutoff
  )
)
context <- InterMetaAIContext(
  analysis = analysis_for_ai,
  disease = "prostate adenocarcinoma",
  endpoint = "biochemical recurrence (BCR)",
  time_unit = "days",
  cohort = "TCGA-PRAD",
  top_n = 5L,
  direction = "unfavorable"
)
if (nrow(context$pathways) != 5L ||
    !all(context$pathways$Direction == "unfavorable")) {
  stop("The AI context did not contain exactly five unfavorable pathways.")
}

context_id <- substr(context$provenance$input_hash, 1L, 16L)
prepared_dir <- file.path(ai_root_dir, paste0("prepared-", context_id))
dir.create(prepared_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(context, file.path(prepared_dir, "ai_context.rds"))
write_tsv(context$pathways, file.path(prepared_dir, "top5_pathways.tsv"))
queries <- InterMetaBuildPubMedQueries(context)
write_tsv(queries, file.path(prepared_dir, "pubmed_query_plan.tsv"))

## 6. Recoverable PubMed, AI, report, and figure workflow --------------------

# A content-affecting AI configuration change receives a separate run
# directory. Authorization flags and credential values are intentionally not
# hashed: they control whether a request may occur, not what the run means.
ai_configuration <- list(
  schema = "InterMeta-example09-ai-config-v1",
  context_input_hash = context$provenance$input_hash,
  query_plan_sha256 = digest::digest(
    queries, algo = "sha256", serialize = TRUE
  ),
  language = "zh",
  purpose = "biology_exploration",
  retmax = 100L,
  page_size = 100L,
  sort = "relevance",
  screen_batch_size = 20L,
  max_findings_per_pathway = 5L,
  model = ai_model,
  model_max_tokens = 32768L,
  model_timeout = 900,
  thinking = "enabled",
  reasoning_effort = "high",
  continue_on_error = TRUE,
  render = TRUE,
  figure_text_mode = "compact_english"
)
ai_configuration_hash <- digest::digest(
  ai_configuration, algo = "sha256", serialize = TRUE
)
writeLines(
  capture.output(dput(ai_configuration)),
  file.path(prepared_dir, "ai_workflow_configuration.txt")
)

ai_result <- NULL
ai_workflow_status <- "disabled"
ai_run_dir <- file.path(
  ai_root_dir,
  paste0("run-", context_id, "-", substr(ai_configuration_hash, 1L, 16L))
)
if (run_ai_workflow) {
  if (!allow_model_calls) {
    message(
      "allow_model_calls = FALSE: no model request can be made. The new ",
      "workflow may retrieve PubMed records and then pause at the model gate."
    )
  }
  if (!file.exists(ai_run_dir) && !dir.exists(ai_run_dir)) {
    ai_result <- InterMetaAIWorkflow(
      context = context,
      run_dir = ai_run_dir,
      queries = queries,
      allow_model_calls = allow_model_calls,
      language = ai_configuration$language,
      purpose = ai_configuration$purpose,
      retmax = ai_configuration$retmax,
      page_size = ai_configuration$page_size,
      sort = ai_configuration$sort,
      ncbi_email_env = "NCBI_EMAIL",
      ncbi_api_key_env = "NCBI_API_KEY",
      screen_batch_size = ai_configuration$screen_batch_size,
      max_findings_per_pathway = ai_configuration$max_findings_per_pathway,
      model = ai_configuration$model,
      api_key_env = ai_key_env,
      model_max_tokens = ai_configuration$model_max_tokens,
      model_timeout = ai_configuration$model_timeout,
      thinking = ai_configuration$thinking,
      reasoning_effort = ai_configuration$reasoning_effort,
      continue_on_error = ai_configuration$continue_on_error,
      render = ai_configuration$render,
      figure_text_mode = ai_configuration$figure_text_mode
    )
  } else if (dir.exists(ai_run_dir)) {
    ai_result <- InterMetaAIResume(
      run_dir = ai_run_dir,
      allow_model_calls = allow_model_calls,
      retry_failed = retry_failed_model_calls
    )
  } else {
    stop("The AI run path exists but is not a directory: ", ai_run_dir)
  }

  workflow_status <- InterMetaAIStatus(ai_run_dir)
  ai_workflow_status <- workflow_status$status
  print(workflow_status$status)
  print(workflow_status$stages)
  if (identical(workflow_status$status, "complete")) {
    print(ai_result$report$exploration_summary)
    print(lapply(ai_result$figure, function(x) x$output))
  } else {
    message(
      "AI workflow status is '", workflow_status$status,
      "'. Inspect InterMetaAIStatus(ai_run_dir), then source this script ",
      "again to continue from validated checkpoints."
    )
  }
} else {
  message(
    "AI execution is disabled; the aggregate context and deterministic ",
    "PubMed query plan were saved under: ", prepared_dir
  )
}

writeLines(
  c(
    "InterMeta TCGA-PRAD workflow status summary",
    paste0("Full-sample checks passed: ", all(parity_checks$Pass)),
    paste0("Successful resamples: ", successful_runs, "/", n_resample),
    paste0("AI context pathways: ", paste(context$pathways$Pathway, collapse = " | ")),
    paste0("AI workflow enabled: ", run_ai_workflow),
    paste0("AI workflow status: ", ai_workflow_status),
    paste0("Model calls authorized: ", allow_model_calls),
    paste0("AI configuration hash: ", ai_configuration_hash),
    paste0(
      "AI run directory (relative to workflow root): ",
      file.path("03_ai_workflow", basename(ai_run_dir))
    )
  ),
  file.path(metadata_dir, "workflow_summary.txt")
)

run_completed_at <- Sys.time()
writeLines(
  c(
    paste0(
      "Run started (UTC): ",
      format(run_started_at, tz = "UTC", usetz = TRUE)
    ),
    paste0(
      "Run completed (UTC): ",
      format(run_completed_at, tz = "UTC", usetz = TRUE)
    ),
    paste0(
      "Elapsed seconds: ",
      format(as.numeric(difftime(
        run_completed_at, run_started_at, units = "secs"
      )), digits = 12)
    )
  ),
  file.path(metadata_dir, "run_timing.txt")
)

# Hash every completed output except this self-referential manifest.
manifest_file <- file.path(metadata_dir, "output_manifest_sha256.tsv")
manifest_paths <- list.files(
  output_dir, recursive = TRUE, full.names = TRUE,
  all.files = TRUE, no.. = TRUE
)
manifest_paths <- manifest_paths[
  file.info(manifest_paths)$isdir %in% FALSE &
    normalizePath(manifest_paths, winslash = "/", mustWork = FALSE) !=
      normalizePath(manifest_file, winslash = "/", mustWork = FALSE)
]
output_root <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
manifest_normalized <- normalizePath(
  manifest_paths, winslash = "/", mustWork = TRUE
)
output_manifest <- data.frame(
  Relative.Path = substring(manifest_normalized, nchar(output_root) + 2L),
  Bytes = as.numeric(file.info(manifest_paths)$size),
  SHA256 = vapply(
    manifest_paths,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE
  ),
  stringsAsFactors = FALSE
)
output_manifest <- output_manifest[
  order(output_manifest$Relative.Path), , drop = FALSE
]
write_tsv(output_manifest, manifest_file)

message("Workflow outputs are available at: ", normalizePath(output_dir))

# After a completed AI run, a strictly offline revalidation/figure replay can
# be written to a NEW directory without PubMed or model calls:
#
# replay_result <- InterMetaAIReplay(
#   ai_run_dir,
#   outdir = file.path(ai_root_dir, paste0("replay-", context_id)),
#   render = TRUE,
#   text_mode = "compact_english"
# )
