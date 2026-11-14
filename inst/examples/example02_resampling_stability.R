# InterMeta internal resampling stability workflow with bundled synthetic data.

library(InterMeta)

data("InterMetaExampleInfo")
data("InterMetaExampleRNA")
data("InterMetaExamplePathwayset")

clinical_info <- InterMetaExampleInfo
bulk_expression <- InterMetaExampleRNA
pathwayset <- InterMetaExamplePathwayset
output_dir <- file.path(getwd(), "InterMeta_example_stability")

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
  stratify_event = FALSE,
  min_events = 5,
  min_censored = 5,
  core = 1,
  outdir = output_dir,
  reuse = TRUE,
  keep_details = FALSE
)

print(stability)
head(stability$summary)
stability$pairwise
head(stability$edge_summary)
stability$run_status
stability_score_plot <- InterMetaStabilityScoreBar(
  stability,
  top = 5,
  min_selection_frequency = 1,
  min_direction_consistency = 1,
  error = "sd"
)
ggplot2::ggsave(
  file.path(output_dir, "stability_pathway_score_bar.pdf"),
  stability_score_plot,
  width = 8.5,
  height = 5.5,
  device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
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
  output_dir,
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
print(stability_score_plot)
message("Example stability files: ", normalizePath(output_dir))

# If you already have score tables, summarize them without rerunning analysis:
# scores <- list(score1, score2, score3)
# stability_summary <- InterMetaSummarizeStability(scores = scores)
