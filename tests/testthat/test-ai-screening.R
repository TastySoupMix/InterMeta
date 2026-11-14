screening_context <- function(pathways = c("P1", "P2")) {
  structure(
    list(
      analysis = list(
        disease = "prostate cancer", endpoint = "BCR", time_unit = "days",
        cohort = "test", direction = "both", top_n = length(pathways)
      ),
      pathways = data.frame(
        Pathway = pathways,
        Direction = rep("unfavorable", length(pathways)),
        Difference = seq_along(pathways), stringsAsFactors = FALSE
      ),
      interactions = data.frame(stringsAsFactors = FALSE),
      driver_gene_pairs = data.frame(stringsAsFactors = FALSE),
      warnings = character(),
      provenance = list(context_schema_version = "1.0", input_hash = "screen-test")
    ),
    class = c("InterMetaAIContext", "list")
  )
}

screening_evidence <- function() {
  abstract <- c(
    "Glycolysis promoted proliferation in prostate cancer cells and xenografts.",
    "Glycolysis promoted proliferation in prostate cancer cells and xenografts.",
    "This unrelated report did not investigate the requested metabolic pathway.",
    "This study measured metabolic adaptation in prostate cancer tissue samples.",
    "Fatty acid oxidation supported treatment resistance in prostate cancer organoids."
  )
  records <- data.frame(
    query_id = c("q1", "q2", "q3", "q4", "q5"),
    pathway = c("P1", "P1", "P1", "P1", "P2"),
    direction = "unfavorable",
    layer = c("biology", "review", "biology", "biology_process", "biology"),
    query = paste("query", 1:5),
    term_scope = "pathway", search_terms = c("P1", "P1", "P1", "P1", "P2"),
    query_version = "3.0", query_warning = "", purpose = "biology_exploration",
    PMID = c("1001", "1001", "1002", "1003", "2001"),
    Title = paste("Article", 1:5), Abstract = abstract,
    Journal = "Journal", Year = "2025", DOI = "", PublicationType = "Journal Article",
    MeSH = "Prostatic Neoplasms", stringsAsFactors = FALSE
  )
  structure(
    list(
      schema_version = "1.0", generated_at = "2026-09-08T00:00:00Z",
      search_date = "2026-09-08",
      queries = unique(records[, c(
        "query_id", "pathway", "direction", "layer", "query", "term_scope",
        "search_terms", "query_version", "query_warning", "purpose"
      )]),
      searches = data.frame(
        query_id = unique(records$query_id), pathway = records$pathway,
        api_status = "ok", truncated = FALSE, error = NA_character_,
        stringsAsFactors = FALSE
      ),
      records = records,
      response_hashes = data.frame(
        query_id = character(), stage = character(), batch = integer(),
        status_code = integer(), sha256 = character(), body_redacted = logical(),
        cached = logical(), stringsAsFactors = FALSE
      ),
      provenance = list(provider = "NCBI E-utilities", database = "PubMed")
    ),
    class = c("InterMetaPubMedEvidence", "list")
  )
}

screening_article <- function(row, include = TRUE, quote = row$abstract,
                              exclusion_reason = "none") {
  list(
    candidate_id = row$candidate_id,
    pathway = row$pathway,
    pmid = row$pmid,
    include = include,
    exclusion_reason = exclusion_reason,
    study_disease = if (include) "prostate cancer" else "other cancer",
    disease_match = if (include) "matched" else "other",
    scope = "pathway",
    evidence_types = if (include) list("cell_experiment") else list("other"),
    biological_processes = if (include) list("metabolic_adaptation") else list("other"),
    effect = if (include) "promotes" else "unclear",
    observation = if (include) "A source-linked biological observation." else
      "No relevant pathway observation was identified.",
    study_context = if (include) "Experimental tumour model." else "Unrelated context.",
    abstract_quote = if (is.na(quote)) "" else quote,
    relevance_reason = if (include) "Relevant to the target pathway and disease." else
      "The article does not inform the target pathway."
  )
}

screening_response <- function(messages, mutate_article = NULL) {
  user <- messages[[length(messages)]]$content
  candidates <- user$candidates
  articles <- lapply(seq_len(nrow(candidates)), function(i) {
    row <- candidates[i, , drop = FALSE]
    if (identical(row$pmid, "1002")) {
      value <- screening_article(
        row, include = FALSE, quote = "",
        exclusion_reason = "not_pathway_relevant"
      )
    } else {
      value <- screening_article(row)
    }
    if (is.function(mutate_article)) value <- mutate_article(value, row)
    value
  })
  list(pathways = list(list(pathway = user$pathway, articles = articles)))
}

test_that("screening deduplicates query hits and isolates a bad article", {
  context <- screening_context()
  evidence <- screening_evidence()
  original <- evidence
  observed_schema <- NULL
  provider <- InterMetaAIProvider(function(messages, schema) {
    observed_schema <<- schema
    screening_response(messages, function(value, row) {
      if (identical(row$pmid, "1003")) {
        value$abstract_quote <-
          "This invented quotation has enough characters but is absent from the abstract."
      }
      value
    })
  }, name = "fixture", model = "fixture-model")

  result <- InterMetaAIScreenEvidence(
    context, evidence, provider = provider, batch_size = 2L
  )

  expect_s3_class(result, "InterMetaAIScreening")
  expect_identical(observed_schema$name, "intermeta_literature_screening")
  expect_identical(evidence, original)
  expect_equal(nrow(result$candidates), 4L)
  expect_identical(result$candidates$Query.IDs[1L], "q1;q2")
  expect_identical(result$candidates$Layers[1L], "biology;review")
  expect_equal(nrow(result$included), 2L)
  expect_identical(result$included$PMID, c("1001", "2001"))
  expect_identical(result$excluded$PMID, "1002")
  expect_identical(result$unresolved$PMID, "1003")
  expect_identical(result$unresolved$Screen.Quote.Match.Method, "source_mismatch")
  expect_identical(result$status, "complete_with_quarantine")
  expect_true(result$provenance$complete)
  expect_identical(result$provenance$quarantine_count, 1L)
  expect_equal(nrow(result$batch_status), 3L)
  expect_identical(result$batch_status$Status, c("complete", "partial", "complete"))

  screened <- InterMetaAIScreenedEvidence(result)
  expect_s3_class(screened, "InterMetaAIScreenedEvidence")
  expect_s3_class(screened, "InterMetaPubMedEvidence")
  expect_equal(nrow(screened$candidate_records), 5L)
  expect_equal(nrow(screened$records), 3L)
  expect_identical(screened$records$PMID, c("1001", "1001", "2001"))
  expect_true(all(screened$records$Final.Status == "included"))
  expect_true(all(screened$records$Screen.Quote.Verified))
  expect_true(all(nzchar(screened$records$Candidate.ID)))
})

test_that("a whole-response failure affects only its pathway batch", {
  context <- screening_context()
  evidence <- screening_evidence()
  provider <- InterMetaAIProvider(function(messages, schema) {
    pathway <- messages[[length(messages)]]$content$pathway
    if (identical(pathway, "P1")) return(list(unexpected = TRUE))
    screening_response(messages)
  }, name = "fixture", model = "fixture-model")

  result <- InterMetaAIScreenEvidence(
    context, evidence, provider = provider, batch_size = 20L,
    continue_on_error = TRUE
  )

  expect_identical(result$batch_status$Status, c("failed", "complete"))
  expect_identical(result$included$PMID, "2001")
  expect_setequal(result$unresolved$PMID, c("1001", "1002", "1003"))
  expect_true(all(grepl(
    "response structure failure", result$unresolved$Validation.Reason,
    fixed = TRUE
  )))

  expect_error(
    InterMetaAIScreenEvidence(
      context, evidence, provider = provider, batch_size = 20L,
      continue_on_error = FALSE
    ),
    "stopped after failed batch"
  )
})

test_that("matching checkpoints resume and failed checkpoints retry explicitly", {
  context <- screening_context("P1")
  evidence <- screening_evidence()
  evidence$records <- evidence$records[evidence$records$pathway == "P1" &
    evidence$records$PMID == "1001", , drop = FALSE]
  evidence$queries <- evidence$queries[evidence$queries$pathway == "P1", , drop = FALSE]
  evidence$searches <- evidence$searches[evidence$searches$pathway == "P1", , drop = FALSE]
  checkpoint_dir <- tempfile("intermeta-screen-checkpoints-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE, force = TRUE), add = TRUE)
  calls <- 0L
  fail <- FALSE
  provider <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    if (fail) stop("fixture transport failure")
    screening_response(messages)
  }, name = "fixture", model = "fixture-model")

  first <- InterMetaAIScreenEvidence(
    context, evidence, provider, checkpoint_dir = checkpoint_dir
  )
  expect_identical(first$status, "complete")
  expect_identical(calls, 1L)
  fail <- TRUE
  resumed <- InterMetaAIScreenEvidence(
    context, evidence, provider, checkpoint_dir = checkpoint_dir,
    resume = TRUE
  )
  expect_identical(calls, 1L)
  expect_true(resumed$batch_status$From.Checkpoint)
  expect_false(resumed$batch_status$Attempted)
  expect_true(resumed$batch_status$Paid.May.Have.Occurred)
  expect_true(nzchar(resumed$batch_status$Checkpoint.File.SHA256))
  expect_identical(resumed$status, "complete")

  failed_dir <- tempfile("intermeta-screen-failed-")
  on.exit(unlink(failed_dir, recursive = TRUE, force = TRUE), add = TRUE)
  calls <- 0L
  failed <- InterMetaAIScreenEvidence(
    context, evidence, provider, checkpoint_dir = failed_dir,
    continue_on_error = TRUE
  )
  expect_identical(failed$status, "failed")
  expect_identical(calls, 1L)
  cached_failure <- InterMetaAIScreenEvidence(
    context, evidence, provider, checkpoint_dir = failed_dir,
    resume = TRUE, retry_failed = FALSE
  )
  expect_identical(calls, 1L)
  expect_true(cached_failure$batch_status$From.Checkpoint)
  fail <- FALSE
  retried <- InterMetaAIScreenEvidence(
    context, evidence, provider, checkpoint_dir = failed_dir,
    resume = TRUE, retry_failed = TRUE
  )
  expect_identical(calls, 2L)
  expect_identical(retried$status, "complete")
  expect_false(retried$batch_status$From.Checkpoint)
})

test_that("older screening responses are revalidated without a model call", {
  context <- screening_context("P1")
  evidence <- screening_evidence()
  evidence$records <- evidence$records[evidence$records$pathway == "P1" &
    evidence$records$PMID == "1001", , drop = FALSE]
  evidence$queries <- evidence$queries[evidence$queries$pathway == "P1", , drop = FALSE]
  evidence$searches <- evidence$searches[evidence$searches$pathway == "P1", , drop = FALSE]
  checkpoint_dir <- tempfile("intermeta-screen-local-revalidation-")
  calls <- 0L
  provider <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    screening_response(messages)
  }, name = "fixture", model = "fixture-model")
  first <- InterMetaAIScreenEvidence(
    context, evidence, provider, checkpoint_dir = checkpoint_dir
  )
  expect_identical(first$status, "complete")
  checkpoint_file <- list.files(checkpoint_dir, full.names = TRUE)[[1L]]
  checkpoint <- readRDS(checkpoint_file)
  checkpoint$status <- "partial"
  checkpoint$decisions <- checkpoint$decisions[0, , drop = FALSE]
  checkpoint$validation <- checkpoint$validation[0, , drop = FALSE]
  checkpoint$validation_revision <- "older-local-validator"
  checkpoint$content_sha256 <- getFromNamespace(
    ".InterMetaAIScreeningCheckpointContentSHA256", "InterMeta"
  )(checkpoint)
  saveRDS(checkpoint, checkpoint_file)

  no_call_provider <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    stop("model must not be called during local revalidation")
  }, name = "fixture", model = "fixture-model")
  replayed <- InterMetaAIScreenEvidence(
    context, evidence, no_call_provider, checkpoint_dir = checkpoint_dir
  )

  expect_identical(replayed$status, "complete")
  expect_identical(calls, 1L)
  expect_true(replayed$batch_status$From.Checkpoint)
  expect_false(replayed$batch_status$Attempted)
  expect_length(list.files(checkpoint_dir), 2L)
})

test_that("invalid identifiers are unresolved and missing abstracts are rule-excluded", {
  context <- screening_context("P1")
  evidence <- screening_evidence()
  records <- evidence$records[c(1L, 3L), , drop = FALSE]
  records$pathway <- "P1"
  records$PMID <- c("not-a-pmid", "9001")
  records$Abstract <- c("A usable abstract that cannot repair an invalid identifier.", NA_character_)
  evidence$records <- records
  calls <- 0L
  provider <- InterMetaAIProvider(function(messages, schema) {
    calls <<- calls + 1L
    stop("provider must not be called")
  })

  result <- InterMetaAIScreenEvidence(context, evidence, provider)

  expect_identical(calls, 0L)
  expect_equal(nrow(result$batch_status), 0L)
  expect_identical(result$unresolved$PMID, "not-a-pmid")
  expect_identical(result$excluded$PMID, "9001")
  expect_identical(
    result$excluded$Screen.Exclusion.Reason, "insufficient_source_text"
  )
  expect_identical(result$status, "complete_with_quarantine")
})

test_that("screened conversion detects source mutation", {
  context <- screening_context("P1")
  evidence <- screening_evidence()
  evidence$records <- evidence$records[1L, , drop = FALSE]
  provider <- InterMetaAIProvider(function(messages, schema) {
    screening_response(messages)
  })
  result <- InterMetaAIScreenEvidence(context, evidence, provider)
  expect_identical(result$status, "complete")
  result$source_evidence$records$Title <- "changed after screening"
  expect_error(
    InterMetaAIScreenedEvidence(result), "source evidence changed"
  )
})

test_that("quotes are grounded only against provider-visible abstract text", {
  context <- screening_context("P1")
  evidence <- screening_evidence()
  evidence$records <- evidence$records[
    evidence$records$pathway == "P1" & evidence$records$PMID == "1001",
    , drop = FALSE
  ]
  evidence$queries <- evidence$queries[evidence$queries$pathway == "P1", , drop = FALSE]
  evidence$searches <- evidence$searches[evidence$searches$pathway == "P1", , drop = FALSE]
  hidden <- "This provider-hidden suffix must never qualify as a grounded quotation."
  evidence$records$Abstract <- paste0(
    paste(rep("A", 4000L), collapse = ""), hidden
  )
  provider <- InterMetaAIProvider(function(messages, schema) {
    response <- screening_response(messages)
    response$pathways[[1L]]$articles[[1L]]$abstract_quote <- hidden
    response
  }, name = "fixture", model = "fixture-model")

  result <- InterMetaAIScreenEvidence(context, evidence, provider)

  expect_identical(result$status, "complete_with_quarantine")
  expect_identical(
    result$unresolved$Screen.Quote.Match.Method, "source_mismatch"
  )
  expect_identical(result$provenance$provider_abstract_max_chars, 4000L)
})

test_that("an empty evidence bundle is complete without contacting DeepSeek", {
  context <- screening_context("P1")
  evidence <- screening_evidence()
  evidence$records <- evidence$records[FALSE, , drop = FALSE]
  result <- InterMetaAIScreenEvidence(context, evidence, provider = NULL)
  expect_identical(result$status, "complete")
  expect_equal(nrow(result$candidates), 0L)
  expect_equal(nrow(result$batch_status), 0L)
  screened <- InterMetaAIScreenedEvidence(result)
  expect_equal(nrow(screened$records), 0L)
  expect_true(all(c("Candidate.ID", "Screen.Scope", "Final.Status") %in%
                    names(screened$records)))
})
