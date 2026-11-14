# Generate the bundled synthetic example data and reference results.
# Run either source-tree copy from the InterMeta package root:
# Rscript data-raw/generate_example_data.R
# Rscript inst/scripts/generate_example_data.R

if (!file.exists("DESCRIPTION") || !dir.exists("R")) {
  stop("Run this script from an InterMeta source-tree root.")
}

source(file.path("R", "InterMeta.R"), local = FALSE)
source(file.path("R", "Cox.R"), local = FALSE)

set.seed(20260728)

pathway_genes <- list(
  "Glycolysis Gluconeogenesis" =
    c("HK2", "PFKP", "GAPDH", "PKM", "LDHA"),
  "Glutathione metabolism" =
    c("GCLC", "GCLM", "GSS", "GSR", "GSTP1"),
  "Arginine and proline metabolism" =
    c("ARG2", "OAT", "PRODH", "PYCR1", "ALDH18A1"),
  "Oxidative phosphorylation" =
    c("NDUFA1", "NDUFS1", "UQCRC1", "COX5A", "ATP5F1A"),
  "TCA cycle" =
    c("CS", "ACO2", "IDH3A", "OGDH", "SUCLG1"),
  "Fatty acid degradation" =
    c("ACADM", "ACADL", "CPT1A", "HADHA", "HADHB")
)

InterMetaExamplePathwayset <- do.call(
  rbind,
  lapply(names(pathway_genes), function(pathway) {
    data.frame(
      Pathway = pathway,
      Gene = pathway_genes[[pathway]],
      Dataset = "synthetic_example",
      stringsAsFactors = FALSE
    )
  })
)
row.names(InterMetaExamplePathwayset) <- NULL

n_sample <- 120L
sample_id <- sprintf("IM%03d", seq_len(n_sample))
latent_risk <- as.numeric(scale(stats::rnorm(n_sample)))

event_time <- stats::rexp(
  n_sample,
  rate = exp(1.35 * latent_risk) / 900
)
censor_time <- stats::rexp(n_sample, rate = 1 / 1100)

InterMetaExampleInfo <- data.frame(
  event = as.integer(event_time <= censor_time),
  time = pmax(1, round(pmin(event_time, censor_time))),
  row.names = sample_id,
  check.names = FALSE
)

pathway_loading <- c(
  "Glycolysis Gluconeogenesis" = -1.45,
  "Glutathione metabolism" = -0.95,
  "Arginine and proline metabolism" = -0.35,
  "Oxidative phosphorylation" = 0.20,
  "TCA cycle" = 0.85,
  "Fatty acid degradation" = 1.40
)

gene_id <- InterMetaExamplePathwayset$Gene
InterMetaExampleRNA <- matrix(
  NA_real_,
  nrow = length(gene_id),
  ncol = n_sample,
  dimnames = list(gene_id, sample_id)
)

for (i in seq_along(gene_id)) {
  pathway <- InterMetaExamplePathwayset$Pathway[i]
  gene_offset <- ((i - 1L) %% 5L - 2L) * 0.04
  InterMetaExampleRNA[i, ] <-
    8 + gene_offset + pathway_loading[[pathway]] * latent_risk +
    stats::rnorm(n_sample, sd = 0.32)
}
InterMetaExampleRNA <- round(InterMetaExampleRNA, 4)

checked_info <- InterMetaCheckEvent(
  InterMetaExampleInfo,
  min_events = 5,
  min_censored = 5
)
processed_rna <- InterMetaPreprocess(
  pathwayset = InterMetaExamplePathwayset,
  rna = InterMetaExampleRNA,
  cohort = row.names(checked_info),
  log = FALSE,
  GMNorm = TRUE
)
example_mps <- InterMetaGenePair(processed_rna)
example_cox <- InterMetaCox(
  mps = example_mps,
  info = checked_info,
  p.adj = "BH",
  parallel = FALSE,
  core = 1,
  keep_mps = FALSE
)
InterMetaExampleCox <- example_cox$cox_res

# Keep the bundled adjusted-result fixture as an explicit QV sensitivity
# example. The ordinary package workflow and runnable examples default to PV.
query_pairs <- row.names(InterMetaExampleCox)[
  is.finite(InterMetaExampleCox$QV) & InterMetaExampleCox$QV < 0.05
]
background <- InterMetaBG(row.names(InterMetaExampleCox))
pathway_pair <- InterMetaPathwayPair(InterMetaExamplePathwayset)

InterMetaExampleEnrich <- InterMetaEnrich(
  gene = query_pairs,
  background = background,
  pathway.pair = pathway_pair,
  p.adj = "BH",
  keep_gene_pairs = "shared"
)
InterMetaExampleEnrich <- InterMetaPermFDR(
  resdata = InterMetaExampleEnrich,
  pathwayset = InterMetaExamplePathwayset,
  gene = query_pairs,
  background = background,
  niter = 100,
  core = 1,
  seed = 1  # legacy permutations use 1:niter
)
InterMetaExampleScore <- InterMetaPathwayScore(
  InterMetaExampleEnrich,
  condition = "Adjust.P.Value",
  cutoff = 0.05
)

if (length(query_pairs) < 10L) {
  stop("Synthetic example produced too few adjusted prognostic gene pairs.")
}
if (!any(InterMetaExampleEnrich$Adjust.P.Value < 0.05, na.rm = TRUE)) {
  stop("Synthetic example produced no adjusted pathway enrichments.")
}
if (nrow(InterMetaExampleScore) < 2L) {
  stop("Synthetic example produced an incomplete pathway score table.")
}

dir.create("data", showWarnings = FALSE)
dir.create(file.path("inst", "extdata", "example"),
           recursive = TRUE, showWarnings = FALSE)

save(
  InterMetaExampleInfo,
  file = file.path("data", "InterMetaExampleInfo.rda"),
  compress = "xz",
  version = 2
)
save(
  InterMetaExampleRNA,
  file = file.path("data", "InterMetaExampleRNA.rda"),
  compress = "xz",
  version = 2
)
save(
  InterMetaExamplePathwayset,
  file = file.path("data", "InterMetaExamplePathwayset.rda"),
  compress = "xz",
  version = 2
)
save(
  InterMetaExampleCox,
  file = file.path("data", "InterMetaExampleCox.rda"),
  compress = "xz",
  version = 2
)
save(
  InterMetaExampleEnrich,
  file = file.path("data", "InterMetaExampleEnrich.rda"),
  compress = "xz",
  version = 2
)
save(
  InterMetaExampleScore,
  file = file.path("data", "InterMetaExampleScore.rda"),
  compress = "xz",
  version = 2
)

clinical_text <- data.frame(
  Sample = row.names(InterMetaExampleInfo),
  InterMetaExampleInfo,
  check.names = FALSE
)
expression_text <- data.frame(
  Gene = row.names(InterMetaExampleRNA),
  InterMetaExampleRNA,
  check.names = FALSE
)

utils::write.table(
  clinical_text,
  file = file.path("inst", "extdata", "example", "conditions.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  expression_text,
  file = file.path("inst", "extdata", "example", "exp_tumor.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  InterMetaExamplePathwayset,
  file = file.path("inst", "extdata", "example", "geneset-example.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "Generated InterMeta example data: ",
  nrow(InterMetaExampleRNA), " genes, ",
  ncol(InterMetaExampleRNA), " samples, ",
  length(query_pairs), " adjusted prognostic pairs, and ",
  sum(InterMetaExampleEnrich$Adjust.P.Value < 0.05, na.rm = TRUE),
  " adjusted pathway interactions."
)
