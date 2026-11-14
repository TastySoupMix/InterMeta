# Recompute the original TimiGP analysis from the bundled real TCGA-PRAD inputs
# through InterMeta's explicit compatibility profile. Cox and 100 permutations
# can take substantial time; set run_permutation to FALSE for a faster run when
# only BH enrichment and pathway-score parity are needed.

library(InterMeta)

data("InterMetaPRADInfo")
data("InterMetaPRADRNA")

pathwayset <- InterMetaLoadGeneset("all")
output_dir <- file.path(getwd(), "InterMeta_PRAD_recomputed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_profile <- "TimiGP-legacy"
pair_condition <- "PV"
pair_cutoff <- 0.05
run_permutation <- TRUE
permutation_niter <- 100L

available_cores <- parallel::detectCores(logical = FALSE)
if (is.na(available_cores)) {
  available_cores <- 2L
}
core <- max(1L, min(8L, available_cores - 1L))
pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

## 1. Validate and preprocess -------------------------------------------------

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
stopifnot(identical(dim(rna), c(861L, 463L)))

## 2. Gene pairs and Cox models ----------------------------------------------

mps <- InterMetaGenePair(
  rna,
  chunk_size = 10000,
  analysis_profile = analysis_profile
)
stopifnot(identical(dim(mps), c(370230L, 463L)))
cox <- InterMetaCox(
  mps = mps,
  info = info,
  p.adj = "BH",
  parallel = core > 1L,
  core = core,
  keep_mps = FALSE,
  analysis_profile = analysis_profile
)
cox_res <- cox$cox_res
rm(mps, cox)
invisible(gc())

saveRDS(cox_res, file.path(output_dir, "cox_res.rds"))
utils::write.table(
  cbind(Metabolic.Gene.Pair = row.names(cox_res), cox_res),
  file = file.path(output_dir, "cox_res.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## 3. Pathway-interaction enrichment -----------------------------------------

query_pairs <- row.names(cox_res)[
  is.finite(cox_res[[pair_condition]]) &
    cox_res[[pair_condition]] < pair_cutoff
]
background <- InterMetaBG(row.names(cox_res))
pathway_pair <- InterMetaPathwayPair(pathwayset)

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

saveRDS(enrich, file.path(output_dir, "enrich.rds"))
utils::write.table(
  enrich,
  file = file.path(output_dir, "enrich.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## 4. Scores, network, and plots ---------------------------------------------

# Original TimiFS scoring used BH-adjusted P < 0.05 even when permutation FDR
# was also calculated and retained in the enrichment table.
score_condition <- "Adjust.P.Value"
score_cutoff <- 0.05

score <- InterMetaPathwayScore(
  enrich,
  condition = score_condition,
  cutoff = score_cutoff,
  analysis_profile = analysis_profile
)
network <- InterMetaPathwayNetwork(
  enrich,
  pathwayset = pathwayset,
  condition = score_condition,
  cutoff = score_cutoff,
  export = TRUE,
  path = output_dir
)
utils::write.table(
  score,
  file = file.path(output_dir, "pathway_score.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

dotplot <- InterMetaDotplot(
  enrich,
  condition = score_condition,
  cutoff = score_cutoff,
  top = 10
)
score_bar <- InterMetaPathwayScoreBar(score, top = 5)
chord_file <- file.path(output_dir, "pathway_interaction_chord.pdf")
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
    output_dir, "selected_pathway_gene_pair_chord.pdf"
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
ggplot2::ggsave(
  file.path(output_dir, "pathway_enrichment_dotplot.pdf"),
  dotplot,
  width = 8,
  height = 7,
  device = pdf_device
)
ggplot2::ggsave(
  file.path(output_dir, "pathway_score_bar.pdf"),
  score_bar,
  width = 8,
  height = 7,
  device = pdf_device
)

message("TCGA-PRAD recomputation finished: ", normalizePath(output_dir))
