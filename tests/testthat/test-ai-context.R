ai_context_score <- function(n = 7L){
  data.frame(
    Index = seq_len(n * 2L),
    Rank = seq_len(n * 2L),
    Signed.Rank = seq_len(n * 2L),
    Pathway = c(paste0("U", seq_len(n)), paste0("F", seq_len(n))),
    OutDegree = seq_len(n * 2L),
    Favorable.Score = seq_len(n * 2L),
    InDegree = rev(seq_len(n * 2L)),
    Unfavorable.Score = rev(seq_len(n * 2L)),
    Difference = c(-seq_len(n), seq_len(n)),
    Group = c(rep("U", n), rep("F", n)),
    stringsAsFactors = FALSE
  )
}

ai_context_enrichment <- function(){
  data.frame(
    Index = 1:4,
    Rank = 1:4,
    Pathway.Interaction = c("F7__U7", "F6__U6", "F5__U5", "F1__U1"),
    Favorable.Pathway = c("F7", "F6", "F5", "F1"),
    Unfavorable.Pathway = c("U7", "U6", "U5", "U1"),
    No.Shared.MGP = c(2, 1, 1, 1),
    P.Value = c(0.001, 0.002, 0.003, 0.5),
    Adjust.P.Value = c(0.004, 0.008, 0.012, 0.8),
    Shared.MGP = c("A_B/C_D", "E_F", "G_H", "X_Y"),
    stringsAsFactors = FALSE
  )
}

test_that("AI context defaults to deterministic top five in both directions", {
  context <- InterMetaAIContext(
    enrichment = ai_context_enrichment(),
    score = ai_context_score(),
    disease = "prostate cancer",
    endpoint = "progression-free survival",
    time_unit = "days",
    pathway_synonyms = list(U7 = c("Zulu", "Alpha"), F7 = "F-seven")
  )

  expect_s3_class(context, "InterMetaAIContext")
  expect_named(
    context,
    c("analysis", "pathways", "interactions", "driver_gene_pairs",
      "warnings", "provenance")
  )
  expect_identical(context$analysis$top_n, 5L)
  expect_identical(context$analysis$direction, "both")
  expect_equal(sum(context$pathways$Direction == "unfavorable"), 5L)
  expect_equal(sum(context$pathways$Direction == "favorable"), 5L)
  expect_identical(
    context$pathways$Pathway[context$pathways$Direction == "unfavorable"],
    paste0("U", 7:3)
  )
  expect_identical(
    context$pathways$Pathway[context$pathways$Direction == "favorable"],
    paste0("F", 7:3)
  )
  expect_identical(
    context$pathways$Synonyms[context$pathways$Pathway == "U7"],
    "Alpha | Zulu"
  )
  expect_true(all(nchar(context$provenance$input_hash) == 64L))
})

test_that("AI context uses full analysis results and direct data take precedence", {
  analysis <- list(
    rna = matrix(1:4, nrow = 2,
                 dimnames = list(c("G1", "G2"), c("PATIENT_A", "PATIENT_B"))),
    info = data.frame(event = c(0, 1), row.names = c("PATIENT_A", "PATIENT_B")),
    enrich = ai_context_enrichment(),
    score = ai_context_score(),
    parameters = list(
      pair_condition = "PV",
      pair_cutoff = 0.05,
      score_condition = "Adjust.P.Value",
      score_cutoff = 0.05,
      analysis_profile = "TimiGP-legacy"
    )
  )
  direct_score <- ai_context_score(6L)
  direct_score$Pathway[direct_score$Pathway == "U6"] <- "DIRECT_U"

  context <- InterMetaAIContext(
    analysis = analysis,
    score = direct_score,
    disease = "prostate cancer",
    endpoint = "biochemical recurrence",
    time_unit = "months",
    direction = "unfavorable"
  )

  expect_true("DIRECT_U" %in% context$pathways$Pathway)
  expect_false("U7" %in% context$pathways$Pathway)
  expect_identical(context$analysis$score_condition, "Adjust.P.Value")
  expect_identical(context$analysis$pair_condition, "PV")
  expect_true(any(grepl("unadjusted P-values", context$warnings, fixed = TRUE)))
})

test_that("AI context extracts safe driver pairs and pathway genes", {
  stability <- list(
    summary = data.frame(
      Pathway = c(paste0("U", 1:7), paste0("F", 1:7)),
      Selection.Count = 3L,
      Selection.Frequency = 1,
      Direction.Consistency = 1,
      Mean.Difference = c(-1:-7, 1:7),
      SD.Difference = 0.1,
      stringsAsFactors = FALSE
    ),
    run_status = data.frame(Status = rep("success", 3L))
  )
  context <- InterMetaAIContext(
    enrichment = ai_context_enrichment(),
    score = ai_context_score(),
    stability = stability,
    disease = "prostate cancer",
    endpoint = "recurrence",
    time_unit = "days"
  )

  expect_named(
    context$driver_gene_pairs,
    c("Pathway", "Direction", "Gene.Pair", "Source.Gene", "Target.Gene",
      "Pathway.Interaction")
  )
  expect_setequal(context$driver_gene_pairs$Gene.Pair, c("A_B", "C_D", "E_F", "G_H"))
  expect_identical(
    context$pathways$Driver.Genes[context$pathways$Pathway == "F7"],
    "A,C"
  )
  expect_identical(
    context$pathways$Driver.Genes[context$pathways$Pathway == "U7"],
    "B,D"
  )
  expect_identical(context$pathways$DriverGenes, context$pathways$Driver.Genes)
  expect_true(any(grepl(
    "internal consistency check", context$warnings, fixed = TRUE
  )))
})

test_that("AI context preserves direction-neutral stability audit fields", {
  stability <- list(
    summary = data.frame(
      Pathway = c(paste0("U", 1:7), paste0("F", 1:7)),
      Selection.Count = 3L,
      Selection.Frequency = 1,
      Favorable.Frequency = c(rep(0, 7), rep(1, 7)),
      Unfavorable.Frequency = c(rep(1, 7), rep(0, 7)),
      Dominant.Direction = c(rep("U", 7), rep("F", 7)),
      Direction.Consistency = 1,
      Mean.Difference = c(-1:-7, 1:7),
      SD.Difference = 0.1,
      Mean.Rank = seq_len(14),
      SD.Rank = 0,
      Mean.Magnitude.Rank = rep(7:1, 2),
      SD.Magnitude.Rank = 0,
      Top5.Frequency = 1,
      Top5.Favorable.Frequency = c(rep(0, 7), rep(1, 7)),
      Top5.Unfavorable.Frequency = c(rep(1, 7), rep(0, 7)),
      Stability.Class = "High",
      Unsafe.Sample.ID = "must-not-leave-context",
      stringsAsFactors = FALSE
    ),
    run_status = data.frame(Status = rep("success", 3L))
  )

  context <- InterMetaAIContext(
    enrichment = ai_context_enrichment(),
    score = ai_context_score(),
    stability = stability,
    disease = "prostate cancer",
    endpoint = "recurrence",
    time_unit = "months",
    direction = "unfavorable"
  )

  expect_true(all(c(
    "Mean.Magnitude.Rank", "SD.Magnitude.Rank", "Top5.Frequency",
    "Top5.Favorable.Frequency", "Top5.Unfavorable.Frequency"
  ) %in% names(context$pathways)))
  expect_false("Unsafe.Sample.ID" %in% names(context$pathways))
})

test_that("AI context returns scientific limitation warnings", {
  context <- InterMetaAIContext(
    score = ai_context_score(2L),
    disease = "prostate cancer"
  )
  expect_true(any(grepl("endpoint is unknown", context$warnings, fixed = TRUE)))
  expect_true(any(grepl("time unit is unknown", context$warnings, fixed = TRUE)))
  expect_true(any(grepl("Stability results were not supplied", context$warnings,
                        fixed = TRUE)))
  expect_true(any(grepl("fewer than top_n", context$warnings, fixed = TRUE)))
})

test_that("AI context does not retain patient-level inputs or sample IDs", {
  patient_ids <- c("PATIENT_SECRET_001", "PATIENT_SECRET_002")
  analysis <- list(
    rna = matrix(1:4, nrow = 2, dimnames = list(c("A", "B"), patient_ids)),
    info = data.frame(event = c(0, 1), row.names = patient_ids),
    selected_samples = patient_ids,
    enrich = ai_context_enrichment(),
    score = ai_context_score(),
    parameters = list(score_condition = "Adjust.P.Value", score_cutoff = 0.05)
  )
  context <- InterMetaAIContext(
    analysis = analysis,
    disease = "prostate cancer",
    endpoint = "recurrence",
    time_unit = "days"
  )
  serialized <- paste(capture.output(dput(context)), collapse = "\n")
  expect_false(any(vapply(patient_ids, grepl, logical(1),
                          x = serialized, fixed = TRUE)))
  expect_false(any(c("rna", "info", "selected_samples") %in% names(context)))
})

test_that("AI context validates empty and malformed inputs", {
  expect_error(
    InterMetaAIContext(score = NULL, disease = "prostate cancer"),
    '"score" must be a non-empty data frame', fixed = TRUE
  )
  expect_error(
    InterMetaAIContext(score = data.frame(), disease = "prostate cancer"),
    '"score" must be a non-empty data frame', fixed = TRUE
  )
  malformed <- ai_context_score()
  malformed$Difference[1] <- 1
  expect_error(
    InterMetaAIContext(score = malformed, disease = "prostate cancer"),
    "inconsistent"
  )
  duplicated <- ai_context_score()
  duplicated$Pathway[2] <- duplicated$Pathway[1]
  expect_error(
    InterMetaAIContext(score = duplicated, disease = "prostate cancer"),
    "unique pathway names"
  )
  expect_error(
    InterMetaAIContext(score = ai_context_score(), disease = ""),
    '"disease" must be one non-empty character value', fixed = TRUE
  )
  expect_error(
    InterMetaAIContext(score = ai_context_score(), disease = "cancer", top_n = 0),
    '"top_n" must be one positive integer', fixed = TRUE
  )
  expect_error(
    InterMetaAIContext(
      enrichment = data.frame(Favorable.Pathway = "A"),
      score = ai_context_score(),
      disease = "cancer"
    ),
    "missing required column"
  )
  expect_error(
    InterMetaAIContext(
      score = ai_context_score(),
      disease = "cancer",
      pathway_synonyms = c("unnamed")
    ),
    "must have non-empty pathway names"
  )
})
