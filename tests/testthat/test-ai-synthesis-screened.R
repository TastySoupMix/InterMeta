synthesis_screened_fixture <- function() {
  context <- structure(list(
    analysis = list(
      disease = "prostate cancer", endpoint = "BCR", time_unit = "months",
      direction = "unfavorable", top_n = 1L
    ),
    pathways = data.frame(
      Pathway = "Oxidative phosphorylation", Direction = "unfavorable",
      Direction.Rank = 1L, Difference = -4.5, stringsAsFactors = FALSE
    ),
    interactions = data.frame(), driver_gene_pairs = data.frame(),
    warnings = character(), provenance = list(input_hash = paste(rep("a", 64L), collapse = ""))
  ), class = c("InterMetaAIContext", "list"))
  empty <- getFromNamespace(".InterMetaAIEmptyEvidence", "InterMeta")()
  empty$records <- data.frame(
    query_id = "q1", pathway = "Oxidative phosphorylation",
    direction = "unfavorable", layer = "biology", query = "query",
    PMID = "32198407", Title = "Mitochondrial pathway study",
    Abstract = paste(
      "Oxidative phosphorylation supported prostate cancer cell growth",
      "and treatment adaptation in the experimental models."
    ),
    Journal = "Example", Year = "2020", DOI = "", PublicationType = "Journal Article",
    MeSH = "Prostatic Neoplasms", term_scope = "pathway", search_terms = "OXPHOS",
    query_version = "3.0", query_warning = "", purpose = "biology_exploration",
    stringsAsFactors = FALSE
  )
  screening_provider <- InterMetaAIProvider(function(messages, schema) {
    candidate <- messages[[2L]]$content$candidates[1L, , drop = FALSE]
    list(pathways = list(list(
      pathway = candidate$pathway,
      articles = list(list(
        candidate_id = candidate$candidate_id,
        pathway = candidate$pathway,
        pmid = candidate$pmid,
        include = TRUE, exclusion_reason = "none",
        study_disease = "prostate cancer", disease_match = "matched",
        scope = "pathway",
        evidence_types = list("cell_experiment", "omics"),
        biological_processes = list("proliferation", "metabolic_adaptation"),
        effect = "promotes",
        observation = "The pathway supported growth and treatment adaptation in the studied models.",
        study_context = "Prostate cancer cell and omics models.",
        abstract_quote = candidate$abstract,
        relevance_reason = "Direct target-disease evidence."
      ))
    )))
  }, name = "screening-fixture", model = "fixed")
  audit <- InterMetaAIScreenEvidence(context, empty, screening_provider)
  evidence <- InterMetaAIScreenedEvidence(audit)
  list(context = context, evidence = evidence)
}

synthesis_screened_response <- function() list(pathways = list(list(
  pathway = "Oxidative phosphorylation",
  biological_plausibility = "strongly_plausible",
  biological_summary = "Verified functional and omics findings support a role in tumour growth and metabolic adaptation.",
  context_dependence = "The evidence is strongest in experimental models and needs broader patient validation.",
  research_hypothesis = "Test whether pathway inhibition reduces recurrence-linked tumour adaptation.",
  limitations = "The supplied evidence represents one publication and is not an independent replication set.",
  selected_candidate_ids = list("c000001")
)))

test_that("screened synthesis maps frozen evidence without a second quote task", {
  f <- synthesis_screened_fixture()
  calls <- 0L
  provider <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    expect_identical(schema$name, "intermeta_screened_pathway_synthesis")
    expect_false("Abstract" %in% names(messages[[2L]]$content$verified_findings))
    expect_true("verified_abstract_quote" %in% names(messages[[2L]]$content$verified_findings))
    synthesis_screened_response()
  }, name = "mock", model = "fixed")
  result <- InterMetaAISynthesizeScreened(f$context, f$evidence, provider)
  expect_s3_class(result, "InterMetaAISynthesisRun")
  expect_true(result$complete)
  expect_identical(calls, 1L)
  expect_true(result$pathway_status$Model.Call.Attempted)
  expect_true(result$pathway_status$Paid.May.Have.Occurred)
  expect_identical(result$interpretation$exploration_findings$PMID, "32198407")
  expect_true(result$interpretation$exploration_findings$Quote.Verified)
  expect_identical(
    result$interpretation$exploration_summary$Biological.Plausibility,
    "strongly_plausible"
  )
  expect_identical(
    result$interpretation$record_selection$screening_sha256,
    f$evidence$screening_sha256
  )
})

test_that("screened synthesis cannot be reported against changed evidence", {
  f <- synthesis_screened_fixture()
  provider <- InterMetaAIProvider(
    function(messages, schema) synthesis_screened_response(),
    name = "mock", model = "fixed"
  )
  synthesis <- InterMetaAISynthesizeScreened(
    f$context, f$evidence, provider
  )
  changed <- f$evidence
  changed$records$Title[1L] <- "changed after synthesis"

  expect_error(
    InterMetaAIReport(
      f$context, evidence = changed,
      interpretation = synthesis$interpretation,
      purpose = "biology_exploration"
    ),
    "Evidence records changed"
  )
  changed <- f$evidence
  changed$screening_sha256 <- paste(rep("0", 64L), collapse = "")
  expect_error(
    InterMetaAIReport(
      f$context, evidence = changed,
      interpretation = synthesis$interpretation,
      purpose = "biology_exploration"
    ),
    "screening audit changed"
  )
})

test_that("failed paid calls are not repeated without explicit retry", {
  f <- synthesis_screened_fixture()
  directory <- tempfile("synthesis-checkpoints-")
  calls <- 0L
  failing <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    stop("transport uncertain")
  }, name = "mock", model = "fixed")
  first <- InterMetaAISynthesizeScreened(
    f$context, f$evidence, failing, outdir = directory
  )
  expect_false(first$complete)
  expect_identical(calls, 1L)
  second <- InterMetaAISynthesizeScreened(
    f$context, f$evidence, failing, outdir = directory
  )
  expect_false(second$complete)
  expect_identical(calls, 1L)
  expect_true(second$pathway_status$Reused)
  expect_true(second$pathway_status$Paid.May.Have.Occurred)
  expect_true(nzchar(second$pathway_status$Checkpoint.SHA256))

  succeeding <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    synthesis_screened_response()
  }, name = "mock", model = "fixed")
  third <- InterMetaAISynthesizeScreened(
    f$context, f$evidence, succeeding, outdir = directory,
    retry_failed = TRUE
  )
  expect_true(third$complete)
  expect_identical(calls, 2L)
  expect_length(list.files(file.path(directory, "calls")), 2L)
})

test_that("synthesis checkpoint references survive moving their outdir", {
  f <- synthesis_screened_fixture()
  directory <- tempfile("synthesis-portable-source-")
  moved <- tempfile("synthesis-portable-moved-")
  calls <- 0L
  provider <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    synthesis_screened_response()
  }, name = "mock", model = "fixed")

  first <- InterMetaAISynthesizeScreened(
    f$context, f$evidence, provider, outdir = directory
  )
  expect_identical(calls, 1L)
  expect_identical(first$checkpoint_path_base, "outdir")
  expect_length(first$checkpoints, 1L)
  expect_identical(dirname(first$checkpoints), "calls")
  expect_false(grepl("^(?:[A-Za-z]:|/|\\\\)", first$checkpoints))
  expect_true(file.exists(file.path(directory, first$checkpoints)))

  expect_true(file.rename(directory, moved))
  moved_result <- readRDS(file.path(moved, "synthesis.rds"))
  expect_identical(moved_result$checkpoint_path_base, "outdir")
  expect_true(all(file.exists(file.path(moved, moved_result$checkpoints))))

  no_call_provider <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    stop("model must not be called after moving compatible checkpoints")
  }, name = "mock", model = "fixed")
  reused <- InterMetaAISynthesizeScreened(
    f$context, f$evidence, no_call_provider, outdir = moved
  )
  expect_true(reused$complete)
  expect_identical(calls, 1L)
  expect_true(reused$pathway_status$Reused)
  expect_false(reused$pathway_status$Model.Call.Attempted)
  expect_identical(reused$checkpoint_path_base, "outdir")
  expect_true(all(file.exists(file.path(moved, reused$checkpoints))))
})

test_that("legacy absolute aggregate checkpoint fields are ignored on reuse", {
  f <- synthesis_screened_fixture()
  directory <- tempfile("synthesis-legacy-source-")
  moved <- tempfile("synthesis-legacy-moved-")
  calls <- 0L
  provider <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    synthesis_screened_response()
  }, name = "mock", model = "fixed")

  first <- InterMetaAISynthesizeScreened(
    f$context, f$evidence, provider, outdir = directory
  )
  legacy_absolute <- normalizePath(
    file.path(directory, first$checkpoints), winslash = "/", mustWork = TRUE
  )
  aggregate <- readRDS(file.path(directory, "synthesis.rds"))
  aggregate$checkpoints <- legacy_absolute
  aggregate$checkpoint_path_base <- NULL
  saveRDS(aggregate, file.path(directory, "synthesis.rds"))

  expect_true(file.rename(directory, moved))
  stale <- readRDS(file.path(moved, "synthesis.rds"))
  expect_identical(stale$checkpoints, legacy_absolute)
  expect_false(file.exists(stale$checkpoints))

  no_call_provider <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    stop("legacy aggregate paths must not control checkpoint reuse")
  }, name = "mock", model = "fixed")
  reused <- InterMetaAISynthesizeScreened(
    f$context, f$evidence, no_call_provider, outdir = moved
  )
  expect_true(reused$complete)
  expect_identical(calls, 1L)
  expect_true(reused$pathway_status$Reused)
  expect_false(reused$pathway_status$Model.Call.Attempted)
  expect_identical(reused$checkpoint_path_base, "outdir")
  expect_false(any(grepl("^(?:[A-Za-z]:|/|\\\\)", reused$checkpoints)))
})

test_that("modified synthesis checkpoints are rejected before reuse", {
  f <- synthesis_screened_fixture()
  directory <- tempfile("synthesis-checkpoint-integrity-")
  provider <- InterMetaAIProvider(
    function(messages, schema) synthesis_screened_response(),
    name = "mock", model = "fixed"
  )
  first <- InterMetaAISynthesizeScreened(
    f$context, f$evidence, provider, outdir = directory
  )
  expect_true(first$complete)
  checkpoint_file <- list.files(
    file.path(directory, "calls"), full.names = TRUE
  )[[1L]]
  checkpoint <- readRDS(checkpoint_file)
  checkpoint$normalized_response$biological_summary <- "modified"
  saveRDS(checkpoint, checkpoint_file)

  expect_error(
    InterMetaAISynthesizeScreened(
      f$context, f$evidence, provider, outdir = directory
    ),
    "checkpoint is incompatible"
  )
})

test_that("an older locally invalid response is revalidated without a model call", {
  f <- synthesis_screened_fixture()
  directory <- tempfile("synthesis-local-revalidation-")
  calls <- 0L
  provider <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    synthesis_screened_response()
  }, name = "mock", model = "fixed")
  first <- InterMetaAISynthesizeScreened(
    f$context, f$evidence, provider, outdir = directory
  )
  expect_true(first$complete)
  checkpoint_file <- list.files(
    file.path(directory, "calls"), full.names = TRUE
  )[[1L]]
  checkpoint <- readRDS(checkpoint_file)
  checkpoint$status <- "validation_failed"
  checkpoint$error <- "old local validator rejected the response"
  checkpoint$normalized_response <- NULL
  checkpoint$validation_revision <- "older-local-validator"
  checkpoint$content_sha256 <- getFromNamespace(
    ".InterMetaAIScreenedSynthesisCheckpointContentSHA256", "InterMeta"
  )(checkpoint)
  saveRDS(checkpoint, checkpoint_file)

  no_call_provider <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    stop("model must not be called during local revalidation")
  }, name = "mock", model = "fixed")
  replayed <- InterMetaAISynthesizeScreened(
    f$context, f$evidence, no_call_provider, outdir = directory
  )

  expect_true(replayed$complete)
  expect_identical(calls, 1L)
  expect_true(replayed$pathway_status$Reused)
  expect_false(replayed$pathway_status$Model.Call.Attempted)
  expect_length(list.files(file.path(directory, "calls")), 2L)
})

test_that("synthesis rejects invented candidate identifiers", {
  f <- synthesis_screened_fixture()
  bad <- synthesis_screened_response()
  bad$pathways[[1L]]$selected_candidate_ids <- list("invented")
  result <- InterMetaAISynthesizeScreened(
    f$context, f$evidence,
    InterMetaAIProvider(function(messages, schema) bad, name = "mock", model = "fixed")
  )
  expect_false(result$complete)
  expect_match(result$pathway_status$Error, "selected_candidate_ids")
})

test_that("a hand-forged screened-evidence class cannot bypass its audit", {
  f <- synthesis_screened_fixture()
  forged <- f$evidence
  forged$screening <- NULL
  provider <- InterMetaAIProvider(
    function(messages, schema) synthesis_screened_response(),
    name = "mock", model = "fixed"
  )

  expect_error(
    InterMetaAISynthesizeScreened(f$context, forged, provider),
    "verifiable screening audit"
  )
})
