# Complete InterMeta server workflow for conditions.txt and exp_tumor.txt.
#
# Suggested server layout:
# project/
#   InterMeta/                 # uploaded InterMeta source package
#   run_geneset_all/
#     conditions.txt
#     exp_tumor.txt
#     run_intermeta.R
#
# Set the three paths below, then source this script in a fresh R session.

## 0. User settings ------------------------------------------------------------

project_dir <- getwd()
package_dir <- file.path(project_dir, "InterMeta")
input_dir <- file.path(project_dir, "run_geneset_all")

condition_file <- file.path(input_dir, "conditions.txt")
expression_file <- file.path(input_dir, "exp_tumor.txt")
workflow_dir <- file.path(input_dir, "InterMeta_output")
full_sample_dir <- file.path(workflow_dir, "full_sample")
stability_dir <- file.path(workflow_dir, "internal_resampling")

reinstall_package <- TRUE
run_full_analysis <- TRUE
run_stability <- TRUE
# Keep TRUE to enforce the documented order in a single invocation.
require_full_sample_before_stability <- TRUE

# TimiGP-legacy is the only analysis mode in InterMeta >= 1.1.4.
analysis_profile <- "TimiGP-legacy"
n_resample <- 3L
sample_fraction <- 0.8
seeds <- 20260701L + seq_len(n_resample) - 1L
pair_condition <- "PV"
pair_cutoff <- 0.05
# Set explicitly according to the documented expression scale.
log_transform <- TRUE
gene_median_normalization <- TRUE
stratify_event <- FALSE
run_permutation <- TRUE
permutation_niter <- 100L
score_condition <- "Adjust.P.Value"
score_cutoff <- 0.05

available_cores <- parallel::detectCores(logical = FALSE)
if (is.na(available_cores)) {
  available_cores <- 2L
}
core <- max(1L, min(8L, available_cores - 1L))
use_parallel_cox <- core > 1L
pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

## 1. Validate paths and install the uploaded package --------------------------

required_files <- c(condition_file, expression_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing input file(s): ", paste(missing_files, collapse = ", "))
}
if (!file.exists(file.path(package_dir, "DESCRIPTION"))) {
  stop("package_dir is not an InterMeta source package: ", package_dir)
}

if (reinstall_package || !requireNamespace("InterMeta", quietly = TRUE)) {
  if ("package:InterMeta" %in% search()) {
    detach("package:InterMeta", unload = TRUE, character.only = TRUE)
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_local(
    package_dir,
    dependencies = TRUE,
    upgrade = "never",
    force = TRUE
  )
}

library(InterMeta)
message("Using InterMeta version ", as.character(packageVersion("InterMeta")))
if (packageVersion("InterMeta") < "1.1.3") {
  stop("InterMeta 1.1.3 or newer is required.")
}

dir.create(workflow_dir, recursive = TRUE, showWarnings = FALSE)

## 2. Read input files ---------------------------------------------------------

read_tabular <- function(file) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(
      data.table::fread(file, sep = "\t", data.table = FALSE,
                        check.names = FALSE)
    ))
  }
  read.delim(
    file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

message("Reading clinical data: ", condition_file)
clinical_table <- read_tabular(condition_file)
clinical_ids <- clinical_table[[1]]
clinical_info <- clinical_table[-1]
row.names(clinical_info) <- clinical_ids

message("Reading expression matrix: ", expression_file)
expression_table <- read_tabular(expression_file)
expression_genes <- expression_table[[1]]
bulk_expression <- as.matrix(expression_table[-1])
storage.mode(bulk_expression) <- "double"
row.names(bulk_expression) <- expression_genes
rm(expression_table)
invisible(gc())

# geneset-all.txt is bundled in InterMeta and is not uploaded separately.
pathwayset <- InterMetaLoadGeneset("all")

common_samples <- intersect(colnames(bulk_expression),
                            row.names(clinical_info))
if (length(common_samples) < 10) {
  stop("Fewer than 10 common samples were found.")
}
bulk_expression <- bulk_expression[, common_samples, drop = FALSE]
clinical_info <- clinical_info[common_samples, , drop = FALSE]

writeLines(
  capture.output(sessionInfo()),
  file.path(workflow_dir, "sessionInfo.txt")
)

## 3. Prepare the compact pathway annotation once -----------------------------

message("Preparing compact pathway annotation")
pathway_pair <- InterMetaPathwayPair(pathwayset)
saveRDS(pathway_pair, file.path(workflow_dir, "pathway_pair.rds"))

## 4. Full-sample analysis -----------------------------------------------------

if (run_full_analysis) {
  dir.create(full_sample_dir, recursive = TRUE, showWarnings = FALSE)
  message("Full analysis 1/7: validating clinical data")
  info <- InterMetaCheckEvent(
    clinical_info,
    min_events = 5,
    min_censored = 5
  )

  message("Full analysis 2/7: preprocessing expression")
  rna <- InterMetaPreprocess(
    pathwayset = pathwayset,
    rna = bulk_expression,
    cohort = row.names(info),
    log = log_transform,
    GMNorm = gene_median_normalization,
    analysis_profile = analysis_profile
  )

  message("Full analysis 3/7: generating gene-pair scores")
  mps <- InterMetaGenePair(
    rna,
    chunk_size = 10000,
    analysis_profile = analysis_profile
  )

  message("Full analysis 4/7: running Cox models")
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
  saveRDS(cox_res, file.path(full_sample_dir, "cox_res.rds"))
  write.table(
    cbind(Metabolic.Gene.Pair = row.names(cox_res), cox_res),
    file = file.path(full_sample_dir, "cox_res.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("Full analysis 5/7: selecting prognostic gene pairs")
  selected_pair <- is.finite(cox_res[[pair_condition]]) &
    cox_res[[pair_condition]] < pair_cutoff
  query_pairs <- row.names(cox_res)[selected_pair]
  background <- InterMetaBG(row.names(cox_res))
  message(length(query_pairs), " query pairs selected using ",
          pair_condition, " < ", pair_cutoff)

  message("Full analysis 6/7: pathway-interaction enrichment")
  enrich <- InterMetaEnrich(
    gene = query_pairs,
    background = background,
    pathway.pair = pathway_pair,
    p.adj = "BH",
    keep_gene_pairs = "shared"
  )
  if (run_permutation) {
    enrich <- InterMetaPermFDR(
      resdata = enrich,
      pathwayset = pathwayset,
      gene = query_pairs,
      background = background,
      niter = permutation_niter,
      core = core,
      seed = 1,
      analysis_profile = analysis_profile
    )
  }
  saveRDS(enrich, file.path(full_sample_dir, "enrich.rds"))
  write.table(
    enrich,
    file = file.path(full_sample_dir, "enrich.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("Full analysis 7/7: scores, network, and plots")
  score <- InterMetaPathwayScore(
    enrich,
    condition = score_condition,
    cutoff = score_cutoff,
    analysis_profile = analysis_profile
  )
  write.table(
    score,
    file = file.path(full_sample_dir, "pathway_score.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  InterMetaPathwayNetwork(
    resdata = enrich,
    pathwayset = pathwayset,
    condition = score_condition,
    cutoff = score_cutoff,
    export = TRUE,
    path = full_sample_dir
  )
  dotplot <- InterMetaDotplot(
    enrich,
    condition = score_condition,
    cutoff = score_cutoff,
    top = 10
  )
  score_bar <- InterMetaPathwayScoreBar(score, top = 5)
  chord_file <- file.path(
    full_sample_dir, "pathway_interaction_chord.pdf"
  )
  grDevices::pdf(
    chord_file, width = 14, height = 14, useDingbats = FALSE
  )
  chord <- tryCatch(
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
  ggplot2::ggsave(
    file.path(full_sample_dir, "pathway_enrichment_dotplot.pdf"),
    dotplot,
    width = 8,
    height = 7,
    device = pdf_device
  )
  ggplot2::ggsave(
    file.path(full_sample_dir, "pathway_score_bar.pdf"),
    score_bar,
    width = 8,
    height = 7,
    device = pdf_device
  )

  gene_chord_candidate <- which(
    is.finite(enrich[[score_condition]]) &
      enrich[[score_condition]] < score_cutoff &
      !is.na(enrich$Shared.MGP) &
      nzchar(enrich$Shared.MGP) &
      is.finite(enrich$No.Shared.MGP) &
      enrich$No.Shared.MGP > 0 &
      enrich$No.Shared.MGP <= 50
  )
  if (length(gene_chord_candidate) > 0L) {
    gene_chord_candidate <- gene_chord_candidate[order(
      enrich[[score_condition]][gene_chord_candidate],
      -enrich$No.Shared.MGP[gene_chord_candidate],
      enrich$Index[gene_chord_candidate]
    )]
    gene_chord_file <- file.path(
      full_sample_dir, "selected_pathway_gene_pair_chord.pdf"
    )
    grDevices::pdf(gene_chord_file, width = 9, height = 9)
    gene_chord <- tryCatch(
      InterMetaGeneChord(
        enrich,
        select = enrich$Index[gene_chord_candidate[1]],
        condition = score_condition
      ),
      finally = grDevices::dev.off()
    )
  }

  message("Full-sample output: ", full_sample_dir)
}

## 5. Internal stability analysis ----------------------------------------------

if (run_stability) {
  if (require_full_sample_before_stability && !run_full_analysis) {
    stop(
      "run_stability = TRUE requires run_full_analysis = TRUE in this ",
      "workflow. Set require_full_sample_before_stability = FALSE only ",
      "when intentionally running the resampling stage by itself."
    )
  }
  if (n_resample < 2L) {
    stop(
      "n_resample must be at least 2 for an internal-resampling stability ",
      "analysis with between-run uncertainty."
    )
  }
  message("Starting internal resampling stability analysis")
  stability <- InterMetaStability(
    rna = bulk_expression,
    info = clinical_info,
    pathwayset = pathwayset,
    pathway.pair = pathway_pair,
    n_resample = n_resample,
    sample_fraction = sample_fraction,
    seeds = seeds,
    p.adj = "BH",
    pair_condition = pair_condition,
    pair_cutoff = pair_cutoff,
    enrich_p.adj = "BH",
    score_condition = score_condition,
    score_cutoff = score_cutoff,
    run_permutation = run_permutation,
    permutation_niter = permutation_niter,
    run_network = FALSE,
    log = log_transform,
    GMNorm = gene_median_normalization,
    parallel = use_parallel_cox,
    core = core,
    outdir = stability_dir,
    reuse = TRUE,
    stratify_event = stratify_event,
    min_events = 5,
    min_censored = 5,
    continue_on_error = TRUE,
    keep_details = FALSE,
    analysis_profile = analysis_profile
  )

  print(stability$run_status)
  successful_runs <- sum(stability$run_status$Status == "success")
  if (successful_runs < n_resample) {
    stop(
      successful_runs, " of ", n_resample,
      " internal resamples succeeded. Inspect run_status.txt and the ",
      "per-run error.txt files, then rerun; successful checkpoints will ",
      "be reused and failed rounds retried."
    )
  }
  stability_score_plot <- InterMetaStabilityScoreBar(
    stability,
    top = 5,
    min_selection_frequency = 1,
    min_direction_consistency = 1,
    error = "sd"
  )
  ggplot2::ggsave(
    file.path(stability_dir, "stability_pathway_score_bar.pdf"),
    stability_score_plot,
    width = 8.5,
    height = 5.5,
    device = pdf_device
  )
  stability_score_plot_data <- attr(stability_score_plot, "plot_data")
  if (!"Display.Pathway" %in% names(stability_score_plot_data)) {
    stability_score_plot_data$Display.Pathway <- as.character(
      stability_score_plot_data$Pathway
    )
  }
  stability_score_plot_data_export <- stability_score_plot_data
  text_columns <- vapply(
    stability_score_plot_data_export,
    function(column) is.character(column) || is.factor(column),
    logical(1)
  )
  stability_score_plot_data_export[text_columns] <- lapply(
    stability_score_plot_data_export[text_columns],
    function(column) gsub("[\t\r\n]+", " ", as.character(column))
  )
  stability_score_plot_data_file <- file.path(
    stability_dir,
    "stability_score_plot_data.txt"
  )
  utils::write.table(
    stability_score_plot_data_export,
    file = stability_score_plot_data_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  stability_score_plot_data_check <- utils::read.delim(
    stability_score_plot_data_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!identical(dim(stability_score_plot_data_check),
                 dim(stability_score_plot_data_export)) ||
      !identical(names(stability_score_plot_data_check),
                 names(stability_score_plot_data_export)) ||
      !identical(as.character(stability_score_plot_data_check$Pathway),
                 as.character(stability_score_plot_data_export$Pathway)) ||
      !identical(as.character(stability_score_plot_data_check$Display.Pathway),
                 as.character(stability_score_plot_data_export$Display.Pathway))) {
    stop("The stability plot-data TSV failed its round-trip validation.")
  }
  message("Stability output: ", stability_dir)
}

writeLines(
  c(
    paste0("Analysis profile: ", analysis_profile),
    paste0("Full-sample analysis enabled: ", run_full_analysis),
    paste0("Internal resampling enabled: ", run_stability),
    paste0(
      "Full-sample completion required before resampling: ",
      require_full_sample_before_stability
    ),
    paste0("Requested resamples: ", n_resample),
    paste0("Sample fraction: ", sample_fraction),
    paste0("Resampling seeds: ", paste(seeds, collapse = ", ")),
    paste0("Event stratification: ", stratify_event),
    paste0("Pair condition: ", pair_condition, " < ", pair_cutoff),
    paste0("Score condition: ", score_condition, " < ", score_cutoff),
    paste0("Permutation enabled: ", run_permutation)
  ),
  file.path(workflow_dir, "workflow_summary.txt")
)
message("InterMeta workflow finished: ", normalizePath(workflow_dir))
