# Implementation for update_prad_provenance.R. Run through the public wrapper.

update_prad_provenance_main <- function() {

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Provide the example09 output directory and acceptance candidate tarball.")
}
if (!file.exists("DESCRIPTION") || !dir.exists("data")) {
  stop("Run this script from an InterMeta source-tree root.")
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The digest package is required to validate release provenance.")
}

description <- read.dcf("DESCRIPTION")
if (!identical(unname(description[1L, "Package"]), "InterMeta") ||
    !identical(unname(description[1L, "Version"]), "1.1.3")) {
  stop("This updater is locked to the InterMeta 1.1.3 release contract.")
}

output_dir <- normalizePath(args[[1L]], mustWork = TRUE)
expected_output_name <- "InterMeta_TCGA_PRAD_complete_1_1_3"
if (!identical(basename(output_dir), expected_output_name)) {
  nested <- file.path(output_dir, expected_output_name)
  if (!dir.exists(nested)) {
    stop("The supplied directory is not an InterMeta 1.1.3 example09 output.")
  }
  output_dir <- normalizePath(nested, mustWork = TRUE)
}
analysis_archive <- normalizePath(args[[2L]], mustWork = TRUE)
if (!identical(basename(analysis_archive), "InterMeta_1.1.3.tar.gz")) {
  stop("The candidate archive must be named InterMeta_1.1.3.tar.gz.")
}

metadata_dir <- file.path(output_dir, "00_metadata")
full_dir <- file.path(output_dir, "01_full_sample", "tables")
stability_dir <- file.path(output_dir, "02_internal_resampling")
settings_file <- file.path(metadata_dir, "run_settings.txt")
summary_file <- file.path(metadata_dir, "workflow_summary.txt")
timing_file <- file.path(metadata_dir, "run_timing.txt")
manifest_file <- file.path(metadata_dir, "output_manifest_sha256.tsv")
run_provenance_file <- file.path(metadata_dir, "prad_data_provenance.txt")
parity_file <- file.path(full_dir, "parity_checks.tsv")
full_result_file <- file.path(full_dir, "full_sample_result.rds")
selected_file <- file.path(full_dir, "selected_gene_pairs.tsv")
stability_status_file <- file.path(stability_dir, "run_status.txt")
stability_result_file <- file.path(stability_dir, "stability_result.rds")
bundled_enrichment_file <- file.path("data", "InterMetaPRADEnrich.rda")
bundled_score_file <- file.path("data", "InterMetaPRADScore.rda")
provenance_file <- file.path("data", "InterMetaPRADProvenance.rda")
bundled_info_file <- file.path("data", "InterMetaPRADInfo.rda")
bundled_rna_file <- file.path("data", "InterMetaPRADRNA.rda")
geneset_file <- file.path(
  "inst", "extdata", "metabolic-genesets", "geneset-all.txt"
)
required_files <- c(
  settings_file, summary_file, timing_file, manifest_file,
  run_provenance_file, parity_file,
  full_result_file, selected_file, stability_status_file,
  stability_result_file, bundled_enrichment_file, bundled_score_file,
  provenance_file, bundled_info_file, bundled_rna_file, geneset_file,
  analysis_archive
)
if (any(!file.exists(required_files))) {
  stop(
    "Missing validation or package-data file(s): ",
    paste(required_files[!file.exists(required_files)], collapse = ", ")
  )
}

read_colon_file <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  divider <- regexpr(": ", lines, fixed = TRUE)
  if (!length(lines) || any(divider < 1L)) {
    stop(basename(path), " does not have the expected key/value schema.")
  }
  keys <- substr(lines, 1L, divider - 1L)
  values <- substring(lines, divider + 2L)
  if (any(!nzchar(keys)) || anyDuplicated(keys)) {
    stop(basename(path), " contains blank or duplicate keys.")
  }
  setNames(values, keys)
}

settings <- read_colon_file(settings_file)
expected_settings <- c(
  "InterMeta version" = "1.1.3",
  "PRAD source database" = "PCaDB",
  "PRAD source URL" = "http://bioinfo.jialab-ucr.org/PCaDB/",
  "Analysis profile" = "TimiGP-legacy",
  "Cox selection" = "PV < 0.05",
  "Pathway scoring" = "Adjust.P.Value < 0.05",
  "Permutations per analysis" = "100",
  "Requested resamples" = "3",
  "Sample fraction" = "0.8",
  "Resample seeds" = "20260701,20260702,20260703",
  "RNG kind" = "Mersenne-Twister | Inversion | Rejection",
  "Event-stratified resampling" = "FALSE",
  "Reuse checkpoints" = "FALSE",
  "Run full sample" = "TRUE",
  "Run stability" = "TRUE",
  "Run AI workflow" = "FALSE",
  "Model calls authorized" = "FALSE"
)
if (anyNA(settings[names(expected_settings)]) ||
    !identical(unname(settings[names(expected_settings)]),
               unname(expected_settings))) {
  stop("The example09 settings do not match the frozen PRAD contract.")
}
hash_keys <- c("Workflow script SHA-256", "Analysis fingerprint")
if (anyNA(settings[hash_keys]) ||
    any(!grepl("^[0-9a-f]{64}$", settings[hash_keys]))) {
  stop("The workflow script or analysis fingerprint is not a SHA-256 value.")
}

archive_listing <- utils::untar(analysis_archive, list = TRUE)
required_archive_members <- c(
  "InterMeta/DESCRIPTION",
  "InterMeta/inst/examples/example09_tcga_prad_complete_workflow.R"
)
if (!all(required_archive_members %in% archive_listing) ||
    any(!startsWith(archive_listing, "InterMeta/"))) {
  stop("The analysis archive does not have the expected source tree.")
}
archive_extract <- tempfile("intermeta-analysis-archive-")
dir.create(archive_extract)
on.exit(unlink(archive_extract, recursive = TRUE, force = TRUE), add = TRUE)
utils::untar(
  analysis_archive,
  files = required_archive_members,
  exdir = archive_extract
)
archive_description <- read.dcf(file.path(
  archive_extract, "InterMeta", "DESCRIPTION"
))
if (!identical(unname(archive_description[1L, "Package"]), "InterMeta") ||
    !identical(unname(archive_description[1L, "Version"]), "1.1.3")) {
  stop("The analysis archive DESCRIPTION is not InterMeta 1.1.3.")
}
archive_script <- file.path(
  archive_extract, "InterMeta", "inst", "examples",
  "example09_tcga_prad_complete_workflow.R"
)
archive_script_sha256 <- digest::digest(
  archive_script, algo = "sha256", file = TRUE
)
if (!identical(
  archive_script_sha256,
  unname(settings[["Workflow script SHA-256"]])
)) {
  stop("The archive does not contain the workflow script that ran.")
}
candidate_archive_sha256 <- digest::digest(
  analysis_archive, algo = "sha256", file = TRUE
)
candidate_archive_bytes <- unname(file.info(analysis_archive)$size)

parity <- utils::read.delim(
  parity_file, stringsAsFactors = FALSE, check.names = FALSE
)
expected_checks <- c(
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
  "bundled_modern_score_reference_parity"
)
if (!identical(names(parity), c("Check", "Observed", "Expected", "Pass")) ||
    !identical(parity$Check, expected_checks) || anyDuplicated(parity$Check) ||
    !all(toupper(as.character(parity$Pass)) == "TRUE") ||
    !identical(as.character(parity$Observed), as.character(parity$Expected))) {
  stop("The 12 full-sample parity checks did not all pass exactly.")
}

run_status <- utils::read.delim(
  stability_status_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = character()
)
run_status$Error[is.na(run_status$Error)] <- ""
required_status_columns <- c(
  "Run", "Run.ID", "Seed", "Permutation.Seed", "Status", "Samples",
  "Events", "Censored", "Query.Pairs", "Score.Pathways",
  "Elapsed.Seconds", "Error"
)
expected_seeds <- 20260701:20260703
status_ok <- all(required_status_columns %in% names(run_status)) &&
  nrow(run_status) == 3L &&
  identical(as.integer(run_status$Run), 1:3) &&
  identical(run_status$Run.ID, sprintf("resample_%03d", 1:3)) &&
  identical(as.integer(run_status$Seed), expected_seeds) &&
  all(as.integer(run_status$Permutation.Seed) == 1L) &&
  all(run_status$Status == "success") &&
  all(as.integer(run_status$Samples) == 370L) &&
  all(as.integer(run_status$Events) > 0L) &&
  all(as.integer(run_status$Censored) > 0L) &&
  all(as.integer(run_status$Events) + as.integer(run_status$Censored) == 370L) &&
  all(as.integer(run_status$Query.Pairs) > 0L) &&
  all(as.integer(run_status$Score.Pathways) > 0L) &&
  all(is.finite(as.numeric(run_status$Elapsed.Seconds))) &&
  all(as.numeric(run_status$Elapsed.Seconds) >= 0) &&
  all(!nzchar(trimws(as.character(run_status$Error))))
if (!isTRUE(status_ok)) {
  stop("The three internal resamples did not satisfy the frozen contract.")
}

stability <- readRDS(stability_result_file)
stability_expected <- list(
  analysis_profile = "TimiGP-legacy",
  n_resample = 3L,
  sample_fraction = 0.8,
  seeds = expected_seeds,
  pair_condition = "PV",
  pair_cutoff = 0.05,
  score_condition = "Adjust.P.Value",
  score_cutoff = 0.05,
  run_permutation = TRUE,
  permutation_niter = 100L,
  log = TRUE,
  GMNorm = TRUE,
  stratify_event = FALSE
)
if (!inherits(stability, "InterMetaStability") ||
    !is.list(stability$parameters) ||
    !all(vapply(names(stability_expected), function(key) {
      identical(stability$parameters[[key]], stability_expected[[key]])
    }, logical(1))) ||
    !isTRUE(all.equal(
      run_status[, required_status_columns],
      stability$run_status[, required_status_columns],
      check.attributes = FALSE,
      tolerance = 1e-12
    ))) {
  stop("The saved stability object does not match its table or contract.")
}

manifest <- utils::read.delim(
  manifest_file, stringsAsFactors = FALSE, check.names = FALSE
)
if (!identical(names(manifest), c("Relative.Path", "Bytes", "SHA256")) ||
    !nrow(manifest) || anyDuplicated(manifest$Relative.Path) ||
    !identical(manifest$Relative.Path, sort(manifest$Relative.Path)) ||
    any(!nzchar(manifest$Relative.Path)) ||
    any(grepl("\\\\", manifest$Relative.Path)) ||
    any(grepl("(^|/)\\.\\.(/|$)", manifest$Relative.Path)) ||
    any(grepl("^(/|[A-Za-z]:)", manifest$Relative.Path)) ||
    any(!grepl("^[0-9a-f]{64}$", manifest$SHA256)) ||
    any(!is.finite(as.numeric(manifest$Bytes))) ||
    any(as.numeric(manifest$Bytes) < 0)) {
  stop("The output manifest schema or path/hash values are invalid.")
}
all_output_files <- list.files(
  output_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE
)
all_output_files <- all_output_files[
  file.info(all_output_files)$isdir %in% FALSE
]
normalized_root <- normalizePath(output_dir, winslash = "/")
relative_files <- substring(
  normalizePath(all_output_files, winslash = "/"),
  nchar(normalized_root) + 2L
)
keep <- relative_files != "00_metadata/output_manifest_sha256.tsv"
all_output_files <- all_output_files[keep]
relative_files <- relative_files[keep]
file_order <- order(relative_files)
all_output_files <- all_output_files[file_order]
relative_files <- relative_files[file_order]
if (!identical(manifest$Relative.Path, relative_files)) {
  stop("The output manifest does not enumerate every output file exactly once.")
}
observed_bytes <- unname(file.info(all_output_files)$size)
observed_hashes <- vapply(
  all_output_files,
  digest::digest,
  character(1),
  algo = "sha256",
  file = TRUE
)
if (!identical(as.numeric(manifest$Bytes), as.numeric(observed_bytes)) ||
    !identical(manifest$SHA256, unname(observed_hashes))) {
  stop("At least one output file does not match its recorded SHA-256.")
}
output_manifest_sha256 <- digest::digest(
  manifest_file, algo = "sha256", file = TRUE
)

summary_lines <- readLines(summary_file, warn = FALSE, encoding = "UTF-8")
if (!any(summary_lines == "Full-sample checks passed: TRUE") ||
    !any(summary_lines == "Successful resamples: 3/3")) {
  stop("The workflow summary does not report a complete successful analysis.")
}
timing <- read_colon_file(timing_file)
elapsed <- suppressWarnings(as.numeric(timing[["Elapsed seconds"]]))
if (is.na(timing[["Run completed (UTC)"]]) ||
    !nzchar(timing[["Run completed (UTC)"]]) ||
    !is.finite(elapsed) || elapsed <= 0) {
  stop("The workflow timing record is incomplete.")
}

full_result <- readRDS(full_result_file)
expected_full_parameters <- list(
  analysis_profile = "TimiGP-legacy",
  pair_condition = "PV",
  pair_cutoff = 0.05,
  score_condition = "Adjust.P.Value",
  score_cutoff = 0.05,
  permutation_niter = 100L,
  log = TRUE,
  GMNorm = TRUE
)
if (!is.list(full_result) || !is.data.frame(full_result$enrich) ||
    !is.data.frame(full_result$score) ||
    !all(vapply(names(expected_full_parameters), function(key) {
      identical(full_result$parameters[[key]], expected_full_parameters[[key]])
    }, logical(1))) ||
    !identical(full_result$analysis_fingerprint,
               unname(settings[["Analysis fingerprint"]])) ||
    !identical(as.integer(full_result$preprocessed_dimension), c(861L, 463L)) ||
    !identical(as.integer(full_result$cox_gene_pairs), 370230L) ||
    nrow(full_result$enrich) != 2652L || nrow(full_result$score) != 50L) {
  stop("The saved full-sample result does not match the PRAD contract.")
}
selected_pairs <- utils::read.delim(
  selected_file, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(selected_pairs) != 101323L) {
  stop("The full-sample selected gene-pair count is not 101,323.")
}

bundled_environment <- new.env(parent = emptyenv())
bundled_names <- load(bundled_enrichment_file, envir = bundled_environment)
if (!identical(bundled_names, "InterMetaPRADEnrich")) {
  stop("InterMetaPRADEnrich.rda contains an unexpected object set.")
}
bundled_enrichment <- bundled_environment[["InterMetaPRADEnrich"]]
run_enrichment <- full_result$enrich
required_enrichment_columns <- c(
  "Pathway.Interaction", "Favorable.Pathway", "Unfavorable.Pathway",
  "No.Total.MGP", "No.Shared.MGP", "P.Value", "Adjust.P.Value",
  "Permutation.FDR"
)
if (!all(required_enrichment_columns %in% names(bundled_enrichment)) ||
    !all(required_enrichment_columns %in% names(run_enrichment)) ||
    nrow(bundled_enrichment) != 2652L ||
    anyDuplicated(bundled_enrichment$Pathway.Interaction) ||
    anyDuplicated(run_enrichment$Pathway.Interaction)) {
  stop("The enrichment objects cannot support an unambiguous comparison.")
}
run_match <- match(
  bundled_enrichment$Pathway.Interaction,
  run_enrichment$Pathway.Interaction
)
if (anyNA(run_match)) {
  stop("The fresh and bundled enrichment interaction keys differ.")
}
run_aligned <- run_enrichment[run_match, , drop = FALSE]
text_columns <- c(
  "Pathway.Interaction", "Favorable.Pathway", "Unfavorable.Pathway"
)
count_columns <- c("No.Total.MGP", "No.Shared.MGP")
numeric_columns <- c("P.Value", "Adjust.P.Value")
if (!all(vapply(text_columns, function(column) {
  identical(as.character(bundled_enrichment[[column]]),
            as.character(run_aligned[[column]]))
}, logical(1))) ||
    !all(vapply(count_columns, function(column) {
      identical(as.integer(bundled_enrichment[[column]]),
                as.integer(run_aligned[[column]]))
    }, logical(1))) ||
    !all(vapply(numeric_columns, function(column) {
      isTRUE(all.equal(
        bundled_enrichment[[column]], run_aligned[[column]],
        tolerance = 1e-12, check.attributes = FALSE
      ))
    }, logical(1)))) {
  stop("Fresh enrichment keys, counts, or P-values differ from the bundle.")
}
bundled_fdr <- bundled_enrichment$Permutation.FDR
fresh_fdr <- run_aligned$Permutation.FDR
if (!identical(is.na(bundled_fdr), is.na(fresh_fdr)) ||
    any(!is.finite(fresh_fdr)) || any(fresh_fdr < 0 | fresh_fdr > 1)) {
  stop("Fresh permutation FDR values are incomplete or invalid.")
}
prior_bundle_max_abs_difference <- max(
  abs(bundled_fdr - fresh_fdr), na.rm = TRUE
)
if (!is.finite(prior_bundle_max_abs_difference) ||
    prior_bundle_max_abs_difference > 1e-12) {
  stop("Fresh permutation FDR values do not reproduce the previous bundle.")
}

score_environment <- new.env(parent = emptyenv())
score_names <- load(bundled_score_file, envir = score_environment)
if (!identical(score_names, "InterMetaPRADScore") ||
    !is.data.frame(score_environment[["InterMetaPRADScore"]]) ||
    nrow(score_environment[["InterMetaPRADScore"]]) != 50L) {
  stop("InterMetaPRADScore.rda has an unexpected object or dimension.")
}

provenance_environment <- new.env(parent = emptyenv())
provenance_names <- load(provenance_file, envir = provenance_environment)
if (!identical(provenance_names, "InterMetaPRADProvenance")) {
  stop("InterMetaPRADProvenance.rda contains an unexpected object set.")
}
InterMetaPRADProvenance <-
  provenance_environment[["InterMetaPRADProvenance"]]
expected_input_names <- c(
  "conditions", "expression", "geneset", "prad_cox",
  "prad_enrichment_reference"
)
if (!identical(names(InterMetaPRADProvenance$input_sha256),
               expected_input_names) ||
    any(!grepl("^[0-9a-f]{64}$", InterMetaPRADProvenance$input_sha256))) {
  stop("The existing PRAD input fingerprints have an unexpected schema.")
}

validation_id <- paste(
  basename(dirname(dirname(output_dir))),
  basename(dirname(output_dir)),
  basename(output_dir),
  sep = "/"
)
run_provenance_text <- paste(
  readLines(run_provenance_file, warn = FALSE, encoding = "UTF-8"),
  collapse = "\n"
)
provenance_tokens <- c(
  "PCaDB", "10.1101/2021.06.29.449134",
  unname(InterMetaPRADProvenance$input_sha256)
)
if (!all(vapply(provenance_tokens, function(token) {
  grepl(token, run_provenance_text, fixed = TRUE)
}, logical(1)))) {
  stop("The run's PRAD provenance does not match the bundled input lineage.")
}

input_environment <- new.env(parent = emptyenv())
if (!identical(load(bundled_info_file, envir = input_environment),
               "InterMetaPRADInfo") ||
    !identical(load(bundled_rna_file, envir = input_environment),
               "InterMetaPRADRNA")) {
  stop("The bundled PRAD input objects have an unexpected object set.")
}
function_environment <- new.env(parent = globalenv())
sys.source(file.path("R", "InterMeta.R"), envir = function_environment)
current_pathwayset <- function_environment[["InterMetaLoadGeneset"]](
  file = geneset_file
)
current_fingerprint <- digest::digest(
  list(
    schema = "InterMeta-example09-full-v1",
    package_version = "1.1.3",
    clinical = input_environment[["InterMetaPRADInfo"]],
    expression = input_environment[["InterMetaPRADRNA"]],
    pathwayset = current_pathwayset,
    parameters = expected_full_parameters
  ),
  algo = "sha256",
  serialize = TRUE
)
if (!identical(current_fingerprint,
               unname(settings[["Analysis fingerprint"]]))) {
  stop("The current bundled PRAD inputs do not reproduce the run fingerprint.")
}

existing_fdr_provenance <- InterMetaPRADProvenance$permutation_fdr
same_reproduction <- is.list(existing_fdr_provenance$fresh_reproduction) &&
  identical(
    existing_fdr_provenance$fresh_reproduction$validation_id,
    validation_id
  )
valid_metric <- function(value) {
  is.numeric(value) && length(value) == 1L && is.finite(value) &&
    value >= 0 && value <= 1e-12
}
historical_max_abs_difference <- prior_bundle_max_abs_difference
first_publication_max_abs_difference <- prior_bundle_max_abs_difference
if (same_reproduction &&
    valid_metric(existing_fdr_provenance$generation$
      historical_reference_max_abs_difference)) {
  historical_max_abs_difference <- existing_fdr_provenance$generation$
    historical_reference_max_abs_difference
}
if (same_reproduction) {
  existing_first_metric <- existing_fdr_provenance$fresh_reproduction$
    prior_bundle_max_abs_difference_at_first_publication
  if (!valid_metric(existing_first_metric)) {
    existing_first_metric <- existing_fdr_provenance$fresh_reproduction$
      historical_reference_max_abs_difference
  }
  if (valid_metric(existing_first_metric)) {
    first_publication_max_abs_difference <- existing_first_metric
  }
}

InterMetaPRADProvenance$source$access_date <- NA_character_
InterMetaPRADProvenance$source$database_version <- NA_character_
InterMetaPRADProvenance$source$snapshot <- paste(
  "PCaDB access date, release version, and immutable snapshot were not",
  "recorded in the available source files"
)
InterMetaPRADEnrich <- bundled_enrichment
InterMetaPRADEnrich$Permutation.FDR <- fresh_fdr
InterMetaPRADProvenance$permutation_fdr <- list(
  source = paste(
    "Recomputed from bundled PCaDB-derived TCGA-PRAD expression and",
    "biochemical-recurrence inputs by the fresh InterMeta 1.1.3 example09",
    "workflow with 100 TimiGP-legacy annotation permutations; the historical",
    "enrichment result was used only to verify numerical agreement"
  ),
  iterations = 100L,
  recomputed = TRUE,
  historical_reference_used_for_validation_only = TRUE,
  generation = list(
    script = "data-raw/generate_prad_data.R",
    precomputed_cox_input = TRUE,
    analysis_profile = "TimiGP-legacy",
    niter = 100L,
    seed = 1L,
    historical_reference_max_abs_difference =
      historical_max_abs_difference
  ),
  fresh_reproduction = list(
    validation_id = validation_id,
    package_version = "1.1.3",
    precomputed_statistics_used_for_analysis = FALSE,
    reuse_checkpoints = FALSE,
    analysis_profile = "TimiGP-legacy",
    pair_condition = "PV",
    pair_cutoff = 0.05,
    score_condition = "Adjust.P.Value",
    score_cutoff = 0.05,
    niter = 100L,
    permutation_seed = 1L,
    cox_gene_pairs = 370230L,
    selected_gene_pairs = nrow(selected_pairs),
    enrichment_rows = nrow(run_enrichment),
    score_rows = nrow(full_result$score),
    parity_checks_passed = sum(
      toupper(as.character(parity$Pass)) == "TRUE"
    ),
    resamples_requested = 3L,
    resamples_successful = sum(run_status$Status == "success"),
    sample_fraction = 0.8,
    resample_seeds = expected_seeds,
    analysis_fingerprint = unname(settings[["Analysis fingerprint"]]),
    workflow_script_sha256 =
      unname(settings[["Workflow script SHA-256"]]),
    output_manifest_entries = nrow(manifest),
    output_manifest_sha256 = output_manifest_sha256,
    acceptance_candidate_archive = basename(analysis_archive),
    acceptance_candidate_archive_bytes = candidate_archive_bytes,
    acceptance_candidate_archive_sha256 = candidate_archive_sha256,
    archive_binding = paste(
      "The candidate archive contains the workflow script recorded by the",
      "run; the run did not independently record an archive or core-code hash"
    ),
    completed_at_utc = timing[["Run completed (UTC)"]],
    historical_reference_max_abs_difference =
      historical_max_abs_difference,
    prior_bundle_max_abs_difference_at_first_publication =
      first_publication_max_abs_difference
  )
)

provenance_order <- c(
  "cohort", "source_description", "source", "time_unit", "samples", "events",
  "censored", "expression_genes", "geneset", "geneset_unique_genes",
  "matched_genes", "cox_gene_pairs", "enrichment_query",
  "enrichment_adjustment", "permutation_iterations", "permutation_fdr",
  "score_condition", "input_sha256"
)
if (!setequal(names(InterMetaPRADProvenance), provenance_order)) {
  stop("The updated PRAD provenance fields have an unexpected schema.")
}
InterMetaPRADProvenance <- InterMetaPRADProvenance[provenance_order]

enrichment_temp <- tempfile(
  "InterMetaPRADEnrich-", tmpdir = dirname(bundled_enrichment_file),
  fileext = ".rda"
)
provenance_temp <- tempfile(
  "InterMetaPRADProvenance-", tmpdir = dirname(provenance_file),
  fileext = ".rda"
)
on.exit(unlink(c(enrichment_temp, provenance_temp), force = TRUE), add = TRUE)
save(
  InterMetaPRADEnrich,
  file = enrichment_temp,
  compress = "xz",
  compression_level = 9,
  version = 2
)
save(
  InterMetaPRADProvenance,
  file = provenance_temp,
  compress = "xz",
  compression_level = 9,
  version = 2
)

verify_environment <- new.env(parent = emptyenv())
if (!identical(load(enrichment_temp, envir = verify_environment),
               "InterMetaPRADEnrich") ||
    !identical(verify_environment[["InterMetaPRADEnrich"]]$Permutation.FDR,
               fresh_fdr)) {
  stop("The staged enrichment object failed verification.")
}
verify_environment <- new.env(parent = emptyenv())
if (!identical(load(provenance_temp, envir = verify_environment),
               "InterMetaPRADProvenance") ||
    !identical(
      verify_environment[["InterMetaPRADProvenance"]]$permutation_fdr$
        fresh_reproduction$validation_id,
      validation_id
    )) {
  stop("The staged provenance object failed verification.")
}

enrichment_backup <- tempfile(
  "InterMetaPRADEnrich-backup-", tmpdir = dirname(bundled_enrichment_file),
  fileext = ".rda"
)
provenance_backup <- tempfile(
  "InterMetaPRADProvenance-backup-", tmpdir = dirname(provenance_file),
  fileext = ".rda"
)
on.exit(unlink(c(enrichment_backup, provenance_backup), force = TRUE), add = TRUE)
if (!file.copy(bundled_enrichment_file, enrichment_backup, overwrite = TRUE) ||
    !file.copy(provenance_file, provenance_backup, overwrite = TRUE)) {
  stop("Could not stage package-data backups before publication.")
}
published_enrichment <- file.copy(
  enrichment_temp, bundled_enrichment_file, overwrite = TRUE
)
published_provenance <- published_enrichment && file.copy(
  provenance_temp, provenance_file, overwrite = TRUE
)
if (!published_enrichment || !published_provenance) {
  file.copy(enrichment_backup, bundled_enrichment_file, overwrite = TRUE)
  file.copy(provenance_backup, provenance_file, overwrite = TRUE)
  stop("Package-data publication failed; originals were restored.")
}

final_environment <- new.env(parent = emptyenv())
final_enrichment_names <- load(
  bundled_enrichment_file, envir = final_environment
)
final_fdr <- final_environment[["InterMetaPRADEnrich"]]$Permutation.FDR
final_environment <- new.env(parent = emptyenv())
final_provenance_names <- load(provenance_file, envir = final_environment)
final_validation_id <- final_environment[["InterMetaPRADProvenance"]]$
  permutation_fdr$fresh_reproduction$validation_id
if (!identical(final_enrichment_names, "InterMetaPRADEnrich") ||
    !identical(final_fdr, fresh_fdr) ||
    !identical(final_provenance_names, "InterMetaPRADProvenance") ||
    !identical(final_validation_id, validation_id)) {
  file.copy(enrichment_backup, bundled_enrichment_file, overwrite = TRUE)
  file.copy(provenance_backup, provenance_file, overwrite = TRUE)
  stop("Published package data failed verification; originals were restored.")
}

message(
  "Published recomputed PRAD FDR and provenance from ", validation_id,
  "; historical max |difference| = ",
  format(historical_max_abs_difference, digits = 17),
  "; manifest SHA-256 = ", output_manifest_sha256
)
}

update_prad_provenance_main()
