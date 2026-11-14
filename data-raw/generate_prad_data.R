# Generate package data from the user-curated TCGA-PRAD input and results.
# Usage from an InterMeta source-tree root (either copy is equivalent):
# Rscript data-raw/generate_prad_data.R <input-directory>
# Rscript inst/scripts/generate_prad_data.R <input-directory>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Provide the directory containing conditions.txt, exp_tumor.txt, ",
       "InterMetaPRADCox.rda, and InterMetaPRADEnrichReference.rda. ",
       "The historical Bindea2013c filenames remain accepted with warnings.")
}
if (!file.exists("DESCRIPTION") || !dir.exists("R")) {
  stop("Run this script from an InterMeta source-tree root.")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required to generate the PRAD data.")
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
condition_file <- file.path(input_dir, "conditions.txt")
expression_file <- file.path(input_dir, "exp_tumor.txt")

resolve_prad_reference <- function(canonical, historical) {
  canonical_file <- file.path(input_dir, canonical)
  if (file.exists(canonical_file)) {
    return(canonical_file)
  }
  historical_file <- file.path(input_dir, historical)
  if (file.exists(historical_file)) {
    warning(
      "Using historical PRAD input filename '", historical,
      "'; rename it to '", canonical, "'.",
      call. = FALSE
    )
    return(historical_file)
  }
  canonical_file
}

cox_file <- resolve_prad_reference(
  "InterMetaPRADCox.rda", "Bindea2013c_COX_MP_SKCM06.rda"
)
enrich_file <- resolve_prad_reference(
  "InterMetaPRADEnrichReference.rda", "Bindea2013c_enrich.rda"
)
required_files <- c(condition_file, expression_file, cox_file, enrich_file)
if (any(!file.exists(required_files))) {
  stop("Missing input file(s): ",
       paste(required_files[!file.exists(required_files)], collapse = ", "))
}

source(file.path("R", "InterMeta.R"), local = FALSE)
source(file.path("R", "Cox.R"), local = FALSE)

geneset_file <- file.path(
  "inst", "extdata", "metabolic-genesets", "geneset-all.txt"
)
pathwayset <- InterMetaLoadGeneset(file = geneset_file)

clinical_table <- as.data.frame(data.table::fread(
  condition_file,
  sep = "\t",
  data.table = FALSE,
  check.names = FALSE
))
if (ncol(clinical_table) < 3L) {
  stop("conditions.txt must contain sample ID, event, and time columns.")
}
sample_id <- as.character(clinical_table[[1]])
if (anyNA(sample_id) || any(!nzchar(sample_id)) || anyDuplicated(sample_id)) {
  stop("Clinical sample IDs must be non-empty and unique.")
}
InterMetaPRADInfo <- data.frame(
  event = suppressWarnings(as.numeric(clinical_table[[2]])),
  time = suppressWarnings(as.numeric(clinical_table[[3]])),
  row.names = sample_id,
  check.names = FALSE
)
InterMetaPRADInfo <- InterMetaCheckEvent(
  InterMetaPRADInfo,
  min_events = 5,
  min_censored = 5
)

message("Reading and filtering the full expression table.")
expression_table <- data.table::fread(
  expression_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  check.names = FALSE,
  showProgress = TRUE
)
expression_gene <- as.character(expression_table[[1]])
if (anyNA(expression_gene) || any(!nzchar(expression_gene)) ||
    anyDuplicated(expression_gene)) {
  stop("Expression gene identifiers must be non-empty and unique.")
}

common_samples <- row.names(InterMetaPRADInfo)[
  row.names(InterMetaPRADInfo) %in% names(expression_table)[-1]
]
if (length(common_samples) != nrow(InterMetaPRADInfo)) {
  stop("Not every clinical sample was found in exp_tumor.txt.")
}
matched_genes <- unique(pathwayset$Gene)
matched_genes <- matched_genes[matched_genes %in% expression_gene]
if (length(matched_genes) < 2L) {
  stop("Fewer than two geneset genes were found in exp_tumor.txt.")
}

keep_gene <- expression_gene %in% matched_genes
sample_columns <- match(common_samples, names(expression_table))
selected_table <- expression_table[
  keep_gene,
  c(1L, sample_columns),
  with = FALSE
]
selected_gene <- as.character(selected_table[[1]])
InterMetaPRADRNA <- as.matrix(selected_table[, -1L, with = FALSE])
storage.mode(InterMetaPRADRNA) <- "double"
row.names(InterMetaPRADRNA) <- selected_gene
colnames(InterMetaPRADRNA) <- common_samples
rm(expression_table, selected_table)
invisible(gc())

if (any(!is.finite(InterMetaPRADRNA))) {
  stop("The PRAD metabolic expression matrix contains non-finite values.")
}
if (!identical(colnames(InterMetaPRADRNA), row.names(InterMetaPRADInfo))) {
  stop("Clinical and expression samples are not in identical order.")
}

load_prad_reference <- function(file, canonical, historical) {
  result_environment <- new.env(parent = emptyenv())
  loaded <- load(file, envir = result_environment)
  if (canonical %in% loaded) {
    return(result_environment[[canonical]])
  }
  if (historical %in% loaded) {
    warning(
      "Using historical object name '", historical,
      "' from ", basename(file), "; rename it to '", canonical, "'.",
      call. = FALSE
    )
    return(result_environment[[historical]])
  }
  stop(
    basename(file), " must contain object '", canonical,
    "' (or historical object '", historical, "')."
  )
}

prad_cox <- load_prad_reference(
  cox_file, "InterMetaPRADCox", "Bindea2013c_COX_MP_SKCM06"
)
prad_enrich_reference <- load_prad_reference(
  enrich_file, "InterMetaPRADEnrichReference", "Bindea2013c_enrich"
)

pair_genes <- unique(unlist(strsplit(row.names(prad_cox), "_", fixed = TRUE)))
if (length(setdiff(pair_genes, row.names(InterMetaPRADRNA))) > 0L) {
  stop("The precomputed Cox result contains genes absent from PRAD RNA data.")
}
query_pairs <- row.names(prad_cox)[
  is.finite(prad_cox$PV) & prad_cox$PV < 0.05
]
background <- InterMetaBG(row.names(prad_cox))
pathway_pair <- InterMetaPathwayPair(pathwayset)

message("Recomputing compact enrichment to validate the migrated result.")
InterMetaPRADEnrich <- InterMetaEnrich(
  gene = query_pairs,
  background = background,
  pathway.pair = pathway_pair,
  p.adj = "BH",
  keep_gene_pairs = "none"
)
legacy_key <- paste(
  prad_enrich_reference$Favorable.Cell.Type,
  prad_enrich_reference$Unfavorable.Cell.Type,
  sep = "__"
)
legacy_match <- match(InterMetaPRADEnrich$Pathway.Interaction, legacy_key)
if (anyNA(legacy_match)) {
  stop("Could not match every migrated pathway interaction.")
}
if (!identical(
  InterMetaPRADEnrich$No.Total.MGP,
  as.integer(prad_enrich_reference$No.Total.IMGP[legacy_match])
) || !identical(
  InterMetaPRADEnrich$No.Shared.MGP,
  as.integer(prad_enrich_reference$No.Shared.IMGP[legacy_match])
)) {
  stop("Migrated enrichment counts differ from the legacy result.")
}
if (!isTRUE(all.equal(
  InterMetaPRADEnrich$P.Value,
  prad_enrich_reference$P.Value[legacy_match],
  tolerance = 1e-12
))) {
  stop("Migrated enrichment P-values differ from the legacy result.")
}
# Keep source-tree regeneration single-process. This avoids loading a possibly
# different installed InterMeta namespace on Windows PSOCK workers and ensures
# that the sourced code above is the code that performs every permutation.
permutation_cores <- 1L
message(paste(
  "Recomputing 100-iteration legacy permutation FDR from source inputs",
  "on one process."
))
InterMetaPRADEnrich <- InterMetaPermFDR(
  resdata = InterMetaPRADEnrich,
  pathwayset = pathwayset,
  gene = query_pairs,
  background = background,
  niter = 100L,
  core = permutation_cores,
  seed = 1L,
  analysis_profile = "TimiGP-legacy"
)
historical_fdr <- prad_enrich_reference$Permutation.FDR[legacy_match]
if (!isTRUE(all.equal(
  InterMetaPRADEnrich$Permutation.FDR,
  historical_fdr,
  tolerance = 1e-12
))) {
  stop("Recomputed permutation FDR differs from the historical reference.")
}
InterMetaPRADEnrich <- InterMetaPRADEnrich[, c(
  "Index", "Rank", "Pathway.Interaction", "Favorable.Pathway",
  "Unfavorable.Pathway", "No.Total.MGP", "No.Shared.MGP",
  "Enrichment.Ratio", "P.Value", "Adjust.P.Value", "Permutation.FDR",
  "Shared.MGP", "Total.MGP"
)]

InterMetaPRADScore <- InterMetaPathwayScore(
  InterMetaPRADEnrich,
  condition = "Adjust.P.Value",
  cutoff = 0.05
)

InterMetaPRADProvenance <- list(
  cohort = "TCGA-PRAD",
  source_description = paste(
    "PCaDB TCGA-PRAD biochemical-recurrence and normalized bulk expression",
    "data represented by conditions.txt and exp_tumor.txt"
  ),
  source = list(
    database = "PCaDB",
    cohort = "TCGA-PRAD",
    url = "http://bioinfo.jialab-ucr.org/PCaDB/",
    github = "https://github.com/rli012/PCaDB",
    doi = "10.1101/2021.06.29.449134",
    access_date = NA_character_,
    database_version = NA_character_,
    snapshot = paste(
      "PCaDB access date, release version, and immutable snapshot were not",
      "recorded in the available source files"
    )
  ),
  time_unit = "days",
  samples = nrow(InterMetaPRADInfo),
  events = sum(InterMetaPRADInfo$event == 1),
  censored = sum(InterMetaPRADInfo$event == 0),
  expression_genes = nrow(InterMetaPRADRNA),
  geneset = "all",
  geneset_unique_genes = length(unique(pathwayset$Gene)),
  matched_genes = length(matched_genes),
  cox_gene_pairs = nrow(prad_cox),
  enrichment_query = "legacy Cox PV < 0.05",
  enrichment_adjustment = "BH",
  permutation_iterations = 100L,
  permutation_fdr = list(
    source = paste(
      "Recomputed by the package-data generator from the supplied legacy Cox",
      "query with 100 TimiGP-legacy annotation permutations; the supplied",
      "historical enrichment result was used only to verify numerical",
      "agreement. A verified fresh full-workflow block is added separately",
      "by update_prad_provenance.R"
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
      enrichment_rows = nrow(InterMetaPRADEnrich),
      historical_reference_max_abs_difference = max(
        abs(InterMetaPRADEnrich$Permutation.FDR - historical_fdr),
        na.rm = TRUE
      )
    )
  ),
  score_condition = "Adjust.P.Value < 0.05",
  input_sha256 = c(
    conditions = digest::digest(
      condition_file, algo = "sha256", file = TRUE
    ),
    expression = digest::digest(
      expression_file, algo = "sha256", file = TRUE
    ),
    geneset = digest::digest(
      geneset_file, algo = "sha256", file = TRUE
    ),
    prad_cox = digest::digest(
      cox_file, algo = "sha256", file = TRUE
    ),
    prad_enrichment_reference = digest::digest(
      enrich_file, algo = "sha256", file = TRUE
    )
  )
)

dir.create("data", showWarnings = FALSE)
save(
  InterMetaPRADInfo,
  file = file.path("data", "InterMetaPRADInfo.rda"),
  compress = "xz",
  compression_level = 9,
  version = 2
)
save(
  InterMetaPRADRNA,
  file = file.path("data", "InterMetaPRADRNA.rda"),
  compress = "xz",
  compression_level = 9,
  version = 2
)
save(
  InterMetaPRADEnrich,
  file = file.path("data", "InterMetaPRADEnrich.rda"),
  compress = "xz",
  compression_level = 9,
  version = 2
)
save(
  InterMetaPRADScore,
  file = file.path("data", "InterMetaPRADScore.rda"),
  compress = "xz",
  compression_level = 9,
  version = 2
)
save(
  InterMetaPRADProvenance,
  file = file.path("data", "InterMetaPRADProvenance.rda"),
  compress = "xz",
  compression_level = 9,
  version = 2
)

message(
  "Generated TCGA-PRAD package data: ",
  nrow(InterMetaPRADRNA), " genes x ", ncol(InterMetaPRADRNA),
  " samples; ", nrow(InterMetaPRADEnrich),
  " pathway interactions; ", nrow(InterMetaPRADScore),
  " scored pathways."
)
