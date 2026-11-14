workflow_fixture <- function() {
  context <- structure(list(
    analysis = list(
      disease = "prostate cancer", endpoint = "BCR", time_unit = "months",
      cohort = "public-test", direction = "unfavorable", top_n = 1L
    ),
    pathways = data.frame(
      Index = 1L, Rank = 1L, Signed.Rank = 1L,
      Pathway = "Oxidative phosphorylation", Difference = -5,
      Direction = "unfavorable", Direction.Rank = 1L,
      stringsAsFactors = FALSE
    ),
    interactions = data.frame(), driver_gene_pairs = data.frame(),
    warnings = character(), provenance = list(input_hash = paste(rep("a", 64), collapse = ""))
  ), class = c("InterMetaAIContext", "list"))
  query <- data.frame(
    query_id = "q0001", pathway = "Oxidative phosphorylation",
    direction = "unfavorable", layer = "biology",
    disease = "prostate cancer", endpoint = "BCR", driver_genes = "",
    query = "prostate cancer AND oxidative phosphorylation",
    term_scope = "pathway", search_terms = "oxidative phosphorylation",
    query_version = "3.0", query_warning = "", purpose = "biology_exploration",
    stringsAsFactors = FALSE
  )
  evidence <- getFromNamespace(".InterMetaAIEmptyEvidence", "InterMeta")()
  evidence$schema_version <- "3.0"
  evidence$queries <- query
  evidence$searches <- data.frame(
    query_id = "q0001", pathway = "Oxidative phosphorylation",
    direction = "unfavorable", layer = "biology",
    query = query$query, searched_at = "2026-09-08T00:00:00Z",
    retmax = 10L, total_count = 1L, returned_pmids = 1L,
    records_fetched = 1L, truncated = FALSE, api_status = "ok",
    error = "", cached = TRUE, stringsAsFactors = FALSE
  )
  evidence$records <- data.frame(
    query_id = "q0001", pathway = "Oxidative phosphorylation",
    direction = "unfavorable", layer = "biology", query = query$query,
    PMID = "32198407", Title = "Mitochondrial pathway study",
    Abstract = paste(
      "Oxidative phosphorylation supported prostate cancer cell growth",
      "and treatment adaptation in the experimental models."
    ),
    Journal = "Example", Year = "2020", DOI = "",
    PublicationType = "Journal Article", MeSH = "Prostatic Neoplasms",
    term_scope = "pathway", search_terms = "oxidative phosphorylation",
    query_version = "3.0", query_warning = "", purpose = "biology_exploration",
    stringsAsFactors = FALSE
  )
  list(context = context, queries = query, evidence = evidence)
}

workflow_provider <- function(counter) InterMetaAIProvider(function(messages, schema) {
  counter$n <- counter$n + 1L
  if (identical(schema$name, "intermeta_literature_screening")) {
    article <- messages[[2L]]$content$candidates[1L, , drop = FALSE]
    return(list(pathways = list(list(
      pathway = messages[[2L]]$content$pathway,
      articles = list(list(
        candidate_id = article$candidate_id,
        pathway = article$pathway,
        pmid = article$pmid,
        include = TRUE, exclusion_reason = "none",
        study_disease = "prostate cancer", disease_match = "matched",
        scope = "pathway", evidence_types = list("cell_experiment", "omics"),
        biological_processes = list("proliferation", "metabolic_adaptation"),
        effect = "promotes",
        observation = "The pathway supported tumour growth and adaptation in the studied models.",
        study_context = "Prostate cancer cell and omics models.",
        abstract_quote = article$abstract,
        relevance_reason = "Direct target-disease pathway evidence."
      ))
    ))))
  }
  if (identical(schema$name, "intermeta_screened_pathway_synthesis")) {
    candidate <- messages[[2L]]$content$verified_findings$candidate_id[1L]
    return(list(pathways = list(list(
      pathway = messages[[2L]]$content$pathway$Pathway[1L],
      biological_plausibility = "strongly_plausible",
      biological_summary = "Verified functional and omics evidence supports a role in tumour growth and metabolic adaptation.",
      context_dependence = "Evidence is strongest in experimental models.",
      research_hypothesis = "Test whether pathway inhibition limits tumour adaptation in independent models.",
      limitations = "The frozen bundle contains one publication.",
      selected_candidate_ids = list(candidate)
    ))))
  }
  stop("unexpected schema")
}, name = "workflow-mock", model = "fixed-v1")

workflow_figure_paths <- function(figure) {
  paths <- character()
  visit <- function(value) {
    if (!is.list(value)) return(invisible(NULL))
    if ("output" %in% names(value) && is.list(value$output)) {
      candidate <- unlist(value$output[c("pdf", "png")], use.names = FALSE)
      paths <<- c(paths, as.character(candidate[!is.na(candidate)]))
    }
    nested <- which(vapply(value, is.list, logical(1L)))
    value_names <- names(value)
    if (!is.null(value_names)) {
      nested_names <- value_names[nested]
      nested <- nested[is.na(nested_names) | nested_names != "output"]
    }
    for (index in nested) visit(value[[index]])
    invisible(NULL)
  }
  visit(figure)
  paths[nzchar(paths)]
}

workflow_mock_figure <- function(report, pdf_file, png_file, render, ...) {
  writeLines("portable test PDF", pdf_file, useBytes = TRUE)
  writeLines("portable test PNG", png_file, useBytes = TRUE)
  list(
    output = list(
      pdf = normalizePath(pdf_file, winslash = "/", mustWork = TRUE),
      png = normalizePath(png_file, winslash = "/", mustWork = TRUE)
    ),
    draw = function() NULL
  )
}

workflow_absolute_path <- function(path) {
  grepl("^(?:[A-Za-z]:[/\\\\]|[/\\\\])", path, perl = TRUE)
}

test_that("workflow pauses before paid calls and resumes without repeating them", {
  fixture <- workflow_fixture()
  run_dir <- tempfile("intermeta-workflow-")
  replay_dir <- tempfile("intermeta-replay-")
  cache_dir <- tempfile("intermeta-pubmed-cache-")
  on.exit(unlink(c(run_dir, replay_dir, cache_dir), recursive = TRUE, force = TRUE), add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  provider <- workflow_provider(counter)

  paused <- InterMetaAIWorkflow(
    fixture$context, run_dir, evidence = fixture$evidence,
    queries = fixture$queries, provider = provider,
    allow_model_calls = FALSE, pubmed_cache_dir = cache_dir,
    render = FALSE
  )
  expect_identical(paused$status, "paused")
  expect_identical(counter$n, 0L)
  expect_identical(InterMetaAIStatus(run_dir)$stages$Status[4L], "paused")

  complete <- InterMetaAIResume(
    run_dir, provider = provider, allow_model_calls = TRUE
  )
  expect_identical(complete$status, "complete")
  expect_identical(counter$n, 2L)
  expect_identical(complete$report$purpose, "biology_exploration")
  expect_identical(
    complete$report$exploration_summary$Biological.Plausibility,
    "strongly_plausible"
  )
  status <- InterMetaAIStatus(run_dir)
  expect_true(status$config_sha256_ok)
  expect_true(all(status$stages$Status %in% c("complete", "skipped")))
  expect_true(all(status$stages$Integrity.OK %in% TRUE))

  again <- InterMetaAIResume(run_dir, provider = provider)
  expect_identical(again$status, "complete")
  expect_identical(counter$n, 2L)

  replay <- InterMetaAIReplay(run_dir, replay_dir, render = FALSE)
  expect_s3_class(replay, "InterMetaAIWorkflowReplay")
  expect_identical(counter$n, 2L)
  expect_true(file.exists(file.path(replay_dir, "replay.rds")))
  expect_true(file.exists(file.path(replay_dir, "report", "report.rds")))
})

test_that("a bounded PubMed snapshot is complete but explicitly truncated", {
  fixture <- workflow_fixture()
  fixture$evidence$searches$total_count <- 25L
  fixture$evidence$searches$retmax <- 10L
  fixture$evidence$searches$truncated <- TRUE
  run_dir <- tempfile("intermeta-workflow-truncated-")
  cache_dir <- tempfile("intermeta-pubmed-cache-")
  on.exit(unlink(c(run_dir, cache_dir), recursive = TRUE, force = TRUE),
          add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  provider <- workflow_provider(counter)

  result <- InterMetaAIWorkflow(
    fixture$context, run_dir, evidence = fixture$evidence,
    queries = fixture$queries, provider = provider,
    allow_model_calls = TRUE, pubmed_cache_dir = cache_dir,
    render = FALSE
  )

  expect_identical(result$status, "complete")
  expect_identical(counter$n, 2L)
  pubmed <- result$stages[result$stages$Stage == "pubmed", , drop = FALSE]
  expect_identical(pubmed$Status, "complete_with_truncation")
  expect_match(pubmed$Message, "bounded evidence snapshot")
  expect_match(pubmed$Message, "not an exhaustive literature search")
  status <- InterMetaAIStatus(run_dir)
  expect_true(status$stages$Integrity.OK[status$stages$Stage == "pubmed"])

  resumed <- InterMetaAIResume(run_dir, provider = provider)
  expect_identical(resumed$status, "complete")
  expect_identical(counter$n, 2L)
  expect_identical(
    resumed$stages$Status[resumed$stages$Stage == "pubmed"],
    "complete_with_truncation"
  )
})

test_that("workflow and replay persist portable figure paths", {
  fixture <- workflow_fixture()
  run_dir <- tempfile("intermeta-workflow-portable-figure-")
  replay_dir <- tempfile("intermeta-replay-portable-figure-")
  cache_dir <- tempfile("intermeta-pubmed-cache-")
  on.exit(unlink(c(run_dir, replay_dir, cache_dir), recursive = TRUE,
                 force = TRUE), add = TRUE)
  old_dependencies <- getOption("InterMeta.ai_workflow_dependencies", NULL)
  dependencies <- if (is.list(old_dependencies)) old_dependencies else list()
  dependencies$InterMetaAIExplorationFigure <- workflow_mock_figure
  options(InterMeta.ai_workflow_dependencies = dependencies)
  on.exit(options(InterMeta.ai_workflow_dependencies = old_dependencies),
          add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  provider <- workflow_provider(counter)

  result <- InterMetaAIWorkflow(
    fixture$context, run_dir, evidence = fixture$evidence,
    queries = fixture$queries, provider = provider,
    allow_model_calls = TRUE, pubmed_cache_dir = cache_dir,
    render = TRUE
  )

  stored <- readRDS(file.path(run_dir, "07_figure", "figure.rds"))
  stored_paths <- workflow_figure_paths(stored)
  returned_paths <- workflow_figure_paths(result$figure)
  expect_length(stored_paths, 2L)
  expect_false(any(vapply(stored_paths, workflow_absolute_path, logical(1L))))
  expect_true(all(startsWith(stored_paths, "07_figure/")))
  expect_true(all(vapply(returned_paths, workflow_absolute_path, logical(1L))))
  expect_true(all(file.exists(returned_paths)))

  replay <- InterMetaAIReplay(run_dir, replay_dir, render = TRUE)
  saved_replay <- readRDS(file.path(replay_dir, "replay.rds"))
  replay_stored_paths <- workflow_figure_paths(saved_replay$figure)
  replay_returned_paths <- workflow_figure_paths(replay$figure)
  expect_length(replay_stored_paths, 2L)
  expect_false(any(vapply(
    replay_stored_paths, workflow_absolute_path, logical(1L)
  )))
  expect_true(all(startsWith(replay_stored_paths, "figures/")))
  expect_true(all(vapply(
    replay_returned_paths, workflow_absolute_path, logical(1L)
  )))
  expect_true(all(file.exists(replay_returned_paths)))
})

test_that("resume migrates legacy absolute figure checkpoint paths", {
  fixture <- workflow_fixture()
  run_dir <- tempfile("intermeta-workflow-legacy-figure-")
  cache_dir <- tempfile("intermeta-pubmed-cache-")
  on.exit(unlink(c(run_dir, cache_dir), recursive = TRUE, force = TRUE),
          add = TRUE)
  old_dependencies <- getOption("InterMeta.ai_workflow_dependencies", NULL)
  dependencies <- if (is.list(old_dependencies)) old_dependencies else list()
  dependencies$InterMetaAIExplorationFigure <- workflow_mock_figure
  options(InterMeta.ai_workflow_dependencies = dependencies)
  on.exit(options(InterMeta.ai_workflow_dependencies = old_dependencies),
          add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  provider <- workflow_provider(counter)

  InterMetaAIWorkflow(
    fixture$context, run_dir, evidence = fixture$evidence,
    queries = fixture$queries, provider = provider,
    allow_model_calls = TRUE, pubmed_cache_dir = cache_dir,
    render = TRUE
  )
  figure_file <- file.path(run_dir, "07_figure", "figure.rds")
  legacy <- readRDS(figure_file)
  legacy <- getFromNamespace(
    ".InterMetaAIWorkflowMapFigurePaths", "InterMeta"
  )(legacy, function(path) file.path("Z:/retired-server/job", path))
  saveRDS(legacy, figure_file, version = 3L)
  hash <- digest::digest(figure_file, algo = "sha256", file = TRUE)
  marker_file <- file.path(run_dir, "07_figure", "_SUCCESS.rds")
  marker <- readRDS(marker_file)
  marker$output_sha256 <- hash
  saveRDS(marker, marker_file, version = 3L)
  state_file <- file.path(run_dir, "workflow_state.rds")
  state <- readRDS(state_file)
  figure_index <- match("figure", state$stages$Stage)
  state$stages$Output.SHA256[figure_index] <- hash
  saveRDS(state, state_file, version = 3L)

  resumed <- InterMetaAIResume(run_dir, provider = provider)

  expect_identical(resumed$status, "complete")
  expect_identical(counter$n, 2L)
  migrated_paths <- workflow_figure_paths(readRDS(figure_file))
  returned_paths <- workflow_figure_paths(resumed$figure)
  expect_false(any(vapply(
    migrated_paths, workflow_absolute_path, logical(1L)
  )))
  expect_true(all(startsWith(migrated_paths, "07_figure/")))
  expect_true(all(file.exists(returned_paths)))
  expect_true(all(InterMetaAIStatus(run_dir)$stages$Integrity.OK %in% TRUE))
})

test_that("workflow configuration never serializes provider functions or secrets", {
  fixture <- workflow_fixture()
  run_dir <- tempfile("intermeta-workflow-config-")
  cache_dir <- tempfile("intermeta-pubmed-cache-")
  on.exit(unlink(c(run_dir, cache_dir), recursive = TRUE, force = TRUE), add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  provider <- workflow_provider(counter)
  result <- InterMetaAIWorkflow(
    fixture$context, run_dir, evidence = fixture$evidence,
    queries = fixture$queries, provider = provider,
    pubmed_cache_dir = cache_dir, render = FALSE
  )
  config <- readRDS(file.path(run_dir, "workflow_config.rds"))
  expect_identical(result$status, "paused")
  expect_false(any(vapply(config, is.function, logical(1))))
  expect_identical(config$provider$name, "workflow-mock")
  expect_false("generate" %in% names(config$provider))
})

test_that("two-stage workflow rejects the legacy prognosis-validation purpose", {
  fixture <- workflow_fixture()
  expect_error(
    InterMetaAIWorkflow(
      fixture$context, tempfile("intermeta-invalid-purpose-"),
      evidence = fixture$evidence, queries = fixture$queries,
      purpose = "prognostic_validation", render = FALSE
    ),
    "biology_exploration"
  )
})

test_that("workflow credential parameters accept environment names only", {
  fixture <- workflow_fixture()
  run_dir <- tempfile("intermeta-invalid-env-")

  expect_error(
    InterMetaAIWorkflow(
      fixture$context, run_dir, evidence = fixture$evidence,
      queries = fixture$queries, api_key_env = "not an env name",
      render = FALSE
    ),
    "environment-variable name"
  )
  expect_false(file.exists(run_dir) || dir.exists(run_dir))
})

test_that("one invalid screening row is quarantined without blocking outputs", {
  fixture <- workflow_fixture()
  run_dir <- tempfile("intermeta-workflow-quarantine-")
  cache_dir <- tempfile("intermeta-pubmed-cache-")
  on.exit(unlink(c(run_dir, cache_dir), recursive = TRUE, force = TRUE),
          add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  provider <- InterMetaAIProvider(function(messages, schema) {
    counter$n <- counter$n + 1L
    if (!identical(schema$name, "intermeta_literature_screening")) {
      stop("synthesis must not be called without verified evidence")
    }
    candidate <- messages[[2L]]$content$candidates[1L, , drop = FALSE]
    list(pathways = list(list(
      pathway = candidate$pathway,
      articles = list(list(
        candidate_id = candidate$candidate_id,
        pathway = candidate$pathway, pmid = candidate$pmid,
        include = TRUE, exclusion_reason = "none",
        study_disease = "prostate cancer", disease_match = "matched",
        scope = "pathway", evidence_types = list("omics"),
        biological_processes = list("metabolic_adaptation"),
        effect = "unclear", observation = "A proposed observation.",
        study_context = "Example context.",
        abstract_quote = "This quotation is not present in the supplied abstract.",
        relevance_reason = "Proposed relevance."
      ))
    )))
  }, name = "quarantine-mock", model = "fixed-v1")

  result <- InterMetaAIWorkflow(
    fixture$context, run_dir, evidence = fixture$evidence,
    queries = fixture$queries, provider = provider,
    allow_model_calls = TRUE, pubmed_cache_dir = cache_dir,
    render = FALSE
  )

  expect_identical(result$status, "complete")
  expect_identical(counter$n, 1L)
  status <- InterMetaAIStatus(run_dir)
  screen_row <- status$stages$Stage == "screening"
  expect_identical(
    status$stages$Status[screen_row], "complete_with_quarantine"
  )
  expect_identical(
    result$report$exploration_summary$Biological.Plausibility,
    "not_assessable"
  )
})

test_that("an empty retrieved evidence set completes with model calls disabled", {
  fixture <- workflow_fixture()
  fixture$evidence$records <- fixture$evidence$records[0, , drop = FALSE]
  fixture$evidence$searches$total_count <- 0L
  fixture$evidence$searches$returned_pmids <- 0L
  fixture$evidence$searches$records_fetched <- 0L
  run_dir <- tempfile("intermeta-workflow-empty-")
  cache_dir <- tempfile("intermeta-pubmed-cache-")
  on.exit(unlink(c(run_dir, cache_dir), recursive = TRUE, force = TRUE),
          add = TRUE)
  provider <- InterMetaAIProvider(function(messages, schema) {
    stop("no model call should be made")
  }, name = "no-call-mock", model = "fixed-v1")

  result <- InterMetaAIWorkflow(
    fixture$context, run_dir, evidence = fixture$evidence,
    queries = fixture$queries, provider = provider,
    allow_model_calls = FALSE, pubmed_cache_dir = cache_dir,
    render = FALSE
  )

  expect_identical(result$status, "complete")
  expect_true(all(result$stages$Status %in% c("complete", "skipped")))
})

test_that("workflow resumes an older saved synthesis response with zero model calls", {
  fixture <- workflow_fixture()
  run_dir <- tempfile("intermeta-workflow-local-replay-")
  cache_dir <- tempfile("intermeta-pubmed-cache-")
  on.exit(unlink(c(run_dir, cache_dir), recursive = TRUE, force = TRUE),
          add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  base_provider <- workflow_provider(counter)
  failing_provider <- InterMetaAIProvider(function(messages, schema) {
    value <- base_provider$generate(messages, schema)
    if (identical(schema$name, "intermeta_screened_pathway_synthesis")) {
      value$pathways[[1L]]$selected_candidate_ids <- list("invented")
    }
    value
  }, name = "workflow-mock", model = "fixed-v1")
  failed <- InterMetaAIWorkflow(
    fixture$context, run_dir, evidence = fixture$evidence,
    queries = fixture$queries, provider = failing_provider,
    allow_model_calls = TRUE, pubmed_cache_dir = cache_dir,
    render = FALSE
  )
  expect_identical(failed$status, "failed")
  expect_identical(counter$n, 2L)

  checkpoint_file <- list.files(
    file.path(run_dir, "05_synthesis", "checkpoints", "calls"),
    full.names = TRUE
  )[[1L]]
  checkpoint <- readRDS(checkpoint_file)
  checkpoint$raw_parsed_response$pathways[[1L]]$selected_candidate_ids <-
    list("c000001")
  checkpoint$validation_revision <- "older-local-validator"
  checkpoint$content_sha256 <- getFromNamespace(
    ".InterMetaAIScreenedSynthesisCheckpointContentSHA256", "InterMeta"
  )(checkpoint)
  saveRDS(checkpoint, checkpoint_file)

  no_call_provider <- InterMetaAIProvider(function(messages, schema) {
    counter$n <- counter$n + 1L
    stop("no model call is authorized")
  }, name = "workflow-mock", model = "fixed-v1")
  resumed <- InterMetaAIResume(
    run_dir, provider = no_call_provider, allow_model_calls = FALSE
  )

  expect_identical(resumed$status, "complete")
  expect_identical(counter$n, 2L)
  expect_true(file.exists(file.path(
    run_dir, "06_report", "outputs", "attempt-001", "report.rds"
  )))
})

test_that("a failed PubMed stage follows its portable cache after moving", {
  fixture <- workflow_fixture()
  origin <- tempfile("intermeta-portable-cache-origin-")
  moved <- tempfile("intermeta-portable-cache-moved-")
  dir.create(origin)
  on.exit(unlink(c(origin, moved), recursive = TRUE, force = TRUE), add = TRUE)
  run_dir <- file.path(origin, "run")
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  calls$cache <- character()
  pubmed <- function(context, queries, cache_dir, ...) {
    calls$n <- calls$n + 1L
    calls$cache <- c(calls$cache, normalizePath(
      cache_dir, winslash = "/", mustWork = FALSE
    ))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    value <- fixture$evidence
    if (calls$n == 1L) {
      value$searches$api_status <- "error"
      value$searches$error <- "simulated interrupted retrieval"
    }
    value
  }
  old_dependencies <- getOption("InterMeta.ai_workflow_dependencies", NULL)
  dependencies <- if (is.list(old_dependencies)) old_dependencies else list()
  dependencies$InterMetaPubMedSearch <- pubmed
  options(InterMeta.ai_workflow_dependencies = dependencies)
  on.exit(options(InterMeta.ai_workflow_dependencies = old_dependencies),
          add = TRUE)

  expect_error(
    InterMetaAIWorkflow(
      fixture$context, run_dir, queries = fixture$queries, render = FALSE
    ),
    "PubMed retrieval is incomplete"
  )
  config <- readRDS(file.path(run_dir, "workflow_config.rds"))
  expect_identical(config$pubmed_settings$cache_dir, "pubmed_cache")
  expect_identical(config$pubmed_settings$cache_anchor, "run_parent")
  expect_identical(calls$n, 1L)
  expect_true(file.rename(origin, moved))
  moved_run <- file.path(moved, "run")

  result <- InterMetaAIResume(moved_run)

  expect_identical(result$status, "paused")
  expect_identical(calls$n, 2L)
  expect_identical(
    calls$cache[[2L]],
    normalizePath(file.path(moved, "pubmed_cache"), winslash = "/",
                  mustWork = TRUE)
  )
  expect_true(InterMetaAIStatus(moved_run)$config_sha256_ok)
})

test_that("legacy absolute cache paths map safely or require an override", {
  run_dir <- tempfile("intermeta-legacy-cache-run-")
  dir.create(run_dir)
  existing_cache <- tempfile("intermeta-existing-legacy-cache-")
  dir.create(existing_cache)
  on.exit(unlink(c(run_dir, existing_cache), recursive = TRUE, force = TRUE),
          add = TRUE)
  resolve <- getFromNamespace(
    ".InterMetaAIWorkflowResolveCacheDir", "InterMeta"
  )
  cache_spec <- getFromNamespace(
    ".InterMetaAIWorkflowCacheSpec", "InterMeta"
  )
  expect_error(cache_spec("C:drive-relative-cache", FALSE), "relative path")
  expect_error(
    resolve(
      list(pubmed_settings = list(
        cache_dir = "C:drive-relative-cache", cache_anchor = "run_parent"
      )),
      normalizePath(run_dir, winslash = "/")
    ),
    "relative path"
  )
  existing_legacy <- list(pubmed_settings = list(
    cache_dir = normalizePath(existing_cache, winslash = "/", mustWork = TRUE)
  ))
  expect_identical(
    resolve(existing_legacy, normalizePath(run_dir, winslash = "/")),
    normalizePath(existing_cache, winslash = "/", mustWork = TRUE)
  )
  legacy_default <- list(pubmed_settings = list(
    cache_dir = file.path("Z:/retired-server/project", "pubmed_cache")
  ))
  expect_identical(
    resolve(legacy_default, normalizePath(run_dir, winslash = "/")),
    normalizePath(file.path(dirname(run_dir), "pubmed_cache"),
                  winslash = "/", mustWork = FALSE)
  )

  legacy_custom <- list(pubmed_settings = list(
    cache_dir = file.path("Z:/retired-server/project", "custom-cache")
  ))
  expect_error(
    resolve(legacy_custom, normalizePath(run_dir, winslash = "/")),
    "pubmed_cache_dir"
  )
  expect_identical(
    resolve(
      legacy_custom, normalizePath(run_dir, winslash = "/"),
      override = "recovered-cache"
    ),
    normalizePath(file.path(dirname(run_dir), "recovered-cache"),
                  winslash = "/", mustWork = FALSE)
  )
})

test_that("completed resume and replay neither access PubMed nor lose portability", {
  fixture <- workflow_fixture()
  origin <- tempfile("intermeta-complete-origin-")
  moved <- tempfile("intermeta-complete-moved-")
  moved_again <- tempfile("intermeta-complete-moved-again-")
  dir.create(origin)
  on.exit(unlink(c(origin, moved, moved_again), recursive = TRUE, force = TRUE),
          add = TRUE)
  run_dir <- file.path(origin, "run")
  pubmed_calls <- new.env(parent = emptyenv())
  pubmed_calls$n <- 0L
  old_dependencies <- getOption("InterMeta.ai_workflow_dependencies", NULL)
  dependencies <- if (is.list(old_dependencies)) old_dependencies else list()
  dependencies$InterMetaPubMedSearch <- function(...) {
    pubmed_calls$n <- pubmed_calls$n + 1L
    stop("completed workflows must not access PubMed")
  }
  options(InterMeta.ai_workflow_dependencies = dependencies)
  on.exit(options(InterMeta.ai_workflow_dependencies = old_dependencies),
          add = TRUE)
  model_calls <- new.env(parent = emptyenv())
  model_calls$n <- 0L
  provider <- workflow_provider(model_calls)

  result <- InterMetaAIWorkflow(
    fixture$context, run_dir, evidence = fixture$evidence,
    queries = fixture$queries, provider = provider,
    allow_model_calls = TRUE, render = FALSE
  )
  expect_identical(result$status, "complete")
  expect_identical(pubmed_calls$n, 0L)
  expect_identical(model_calls$n, 2L)
  synthesis_file <- file.path(run_dir, "05_synthesis", "synthesis.rds")
  workflow_synthesis <- readRDS(synthesis_file)
  expect_identical(
    workflow_synthesis$checkpoint_path_base,
    "synthesis_artifact_dir"
  )
  expect_true(length(workflow_synthesis$checkpoints) > 0L)
  expect_true(all(startsWith(
    workflow_synthesis$checkpoints, "checkpoints/calls/"
  )))
  expect_true(all(file.exists(file.path(
    dirname(synthesis_file), workflow_synthesis$checkpoints
  ))))
  engine_synthesis <- readRDS(file.path(
    run_dir, "05_synthesis", "checkpoints", "synthesis.rds"
  ))
  expect_identical(engine_synthesis$checkpoint_path_base, "outdir")
  expect_true(all(startsWith(engine_synthesis$checkpoints, "calls/")))
  expect_true(file.rename(origin, moved))
  moved_run <- file.path(moved, "run")

  resumed <- InterMetaAIResume(
    moved_run, provider = provider,
    pubmed_cache_dir = "an-unused-operational-override"
  )
  replay_dir <- file.path(moved_run, "replay")
  replayed <- InterMetaAIReplay(moved_run, replay_dir, render = FALSE)

  expect_identical(resumed$status, "complete")
  expect_identical(pubmed_calls$n, 0L)
  expect_identical(model_calls$n, 2L)
  saved <- readRDS(file.path(replay_dir, "replay.rds"))
  expect_type(saved$source_run, "list")
  expect_identical(saved$source_run$kind, "relative_path")
  expect_identical(saved$source_run$anchor, "replay_dir")
  expect_identical(saved$source_run$path, "..")
  expect_false(any(vapply(
    unlist(saved$source_run, use.names = FALSE),
    workflow_absolute_path, logical(1L)
  )))
  expect_identical(
    replayed$source_run,
    normalizePath(moved_run, winslash = "/", mustWork = TRUE)
  )

  expect_true(file.rename(moved, moved_again))
  final_run <- file.path(moved_again, "run")
  final_replay <- file.path(final_run, "replay")
  saved_after_move <- readRDS(file.path(final_replay, "replay.rds"))
  resolve_source <- getFromNamespace(
    ".InterMetaAIWorkflowReplaySourceForSession", "InterMeta"
  )
  expect_identical(
    resolve_source(saved_after_move$source_run, final_replay),
    normalizePath(final_run, winslash = "/", mustWork = TRUE)
  )
  expect_identical(
    resolve_source(
      normalizePath(final_run, winslash = "/", mustWork = TRUE),
      final_replay
    ),
    normalizePath(final_run, winslash = "/", mustWork = TRUE)
  )

  config_file <- file.path(final_run, "workflow_config.rds")
  tampered_config <- readRDS(config_file)
  tampered_config$language <- if (identical(tampered_config$language, "en")) {
    "zh"
  } else {
    "en"
  }
  saveRDS(tampered_config, config_file)
  rejected_replay <- file.path(final_run, "replay-after-config-tamper")
  expect_error(
    InterMetaAIReplay(final_run, rejected_replay, render = FALSE),
    "configuration failed its SHA256 check"
  )
  expect_false(file.exists(rejected_replay) || dir.exists(rejected_replay))
})
