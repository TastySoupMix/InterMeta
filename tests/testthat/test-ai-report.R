# Existing clinical-mode regression fixtures explicitly retain their original goal.
clinical_ai_interpret <- function(..., max_records_per_pathway = 20L, purpose = "prognostic_validation")
  InterMeta::InterMetaAIInterpret(..., max_records_per_pathway = max_records_per_pathway, purpose = purpose)
clinical_ai_report <- function(..., max_records_per_pathway = 20L, purpose = "prognostic_validation")
  InterMeta::InterMetaAIReport(..., max_records_per_pathway = max_records_per_pathway, purpose = purpose)

ai_function <- function(name) {
  if (name == "InterMetaAIInterpret") return(clinical_ai_interpret)
  if (name == "InterMetaAIReport") return(clinical_ai_report)
  getFromNamespace(name, "InterMeta")
}

make_ai_context <- function() {
  pathways <- rbind(
    data.frame(
      Pathway = paste0("U", 1:6), Direction = "Unfavorable",
      Direction.Rank = 1:6, Group = "U", Difference = -(6:1),
      Synonyms = paste0("Unfavorable ", 1:6),
      Patient.ID = paste0("PATIENT-U", 1:6),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Pathway = paste0("F", 1:6), Direction = "Favorable",
      Direction.Rank = 1:6, Group = "F", Difference = 6:1,
      Synonyms = paste0("Favorable ", 1:6),
      Patient.ID = paste0("PATIENT-F", 1:6),
      stringsAsFactors = FALSE
    )
  )
  structure(
    list(
      analysis = list(
        disease = "prostate adenocarcinoma", endpoint = "overall survival",
        time_unit = "days", cohort = "toy", direction = "both", top_n = 6L,
        score_condition = "Adjust.P.Value", score_cutoff = 0.05,
        pair_condition = "PV", pair_cutoff = 0.05,
        analysis_profile = "InterMeta",
        patient_data = "PATIENT-SECRET"
      ),
      pathways = pathways,
      interactions = data.frame(
        Pathway1 = "U1", Pathway2 = "F1", Adjust.P.Value = 0.01,
        stringsAsFactors = FALSE
      ),
      driver_gene_pairs = data.frame(
        Pathway.Interaction = "U1--F1", Favorable.Pathway = "F1",
        Unfavorable.Pathway = "U1", Gene.Pair = "GENE1|GENE2",
        Source = "Shared.MGP", stringsAsFactors = FALSE
      ),
      warnings = "Only three exploratory resamples were available.",
      provenance = list(
        context_schema_version = "1.0", input_hash = "abc",
        enrichment_hash = "def", score_hash = "ghi",
        stability_hash = "unknown", pathway_synonyms_hash = "jkl"
      ),
      patient_data = data.frame(
        Patient.ID = "PATIENT-SECRET", Time = 42,
        stringsAsFactors = FALSE
      ),
      rna = matrix(99, nrow = 1L)
    ),
    class = c("InterMetaAIContext", "list")
  )
}

make_ai_evidence <- function() {
  structure(
    list(
      schema_version = "1.0",
      generated_at = "2026-09-02T00:00:00Z",
      search_date = "2026-09-02",
      queries = data.frame(
        query_id = "q1", pathway = "U1", direction = "Unfavorable",
        layer = "disease_pathway", query = "prostate cancer AND U1",
        stringsAsFactors = FALSE
      ),
      searches = data.frame(
        query_id = "q1", pathway = "U1", direction = "Unfavorable",
        layer = "disease_pathway", query = "prostate cancer AND U1",
        total_count = 1L, returned_pmids = 1L, records_fetched = 1L,
        retmax = 10L, truncated = FALSE,
        api_status = "success", error = "", stringsAsFactors = FALSE
      ),
      records = data.frame(
        query_id = "q1", pathway = "U1", direction = "Unfavorable",
        layer = "disease_pathway", query = "prostate cancer AND U1",
        PMID = "12345678", Title = "Evidence for U1", Abstract = "Abstract",
        Journal = "Toy Journal", Year = "2025", DOI = "10.1/toy",
        PublicationType = "Journal Article", MeSH = "Prostatic Neoplasms",
        stringsAsFactors = FALSE
      ),
      response_hashes = data.frame(
        query_id = c("q1", "q1"), stage = c("esearch", "efetch"),
        batch = c(0L, 1L), status_code = c(200L, 200L),
        sha256 = c(paste(rep("a", 64L), collapse = ""),
                   paste(rep("b", 64L), collapse = "")),
        cached = c(FALSE, FALSE), stringsAsFactors = FALSE
      ),
      provenance = list(
        provider = "NCBI E-utilities", database = "PubMed", tool = "InterMeta",
        email_supplied = FALSE, api_key_used = FALSE,
        cache_enabled = FALSE, refresh = FALSE
      )
    ),
    class = c("InterMetaPubMedEvidence", "list")
  )
}

make_fake_provider <- function(pmid = "12345678") {
  provider_constructor <- ai_function("InterMetaAIProvider")
  provider_constructor(
    generate = function(messages, schema) {
      expect_type(messages, "list")
      expect_identical(schema$name, "intermeta_pathway_interpretation")
      list(
        pathways = list(
          list(
            pathway = "U1",
            literature_consistency = "supported",
            summary = "The supplied record is consistent with U1.",
            claims = list(list(
              claim = "The supplied publication discusses U1.",
              pmids = pmid
            ))
          )
        )
      )
    },
    name = "fake", model = "fake-1"
  )
}

test_that("AI provider is callable metadata without a key field", {
  provider_constructor <- ai_function("InterMetaAIProvider")
  expect_error(provider_constructor("not a function"), "generate must be")
  provider <- make_fake_provider()
  expect_s3_class(provider, "InterMetaAIProvider")
  expect_identical(names(provider), c("generate", "name", "model"))
  expect_false(any(grepl("key|token|secret", names(provider), ignore.case = TRUE)))
})

test_that("fake provider produces citation-constrained interpretation", {
  interpret <- ai_function("InterMetaAIInterpret")
  result <- interpret(
    make_ai_context(), make_ai_evidence(), make_fake_provider(), language = "en"
  )
  expect_s3_class(result, "InterMetaAIInterpretation")
  expect_identical(result$pathways$Literature.Consistency, "insufficient_evidence")
  expect_identical(result$pathways$Model.Literature.Consistency, "supported")
  expect_identical(result$assessments$Assessment.Status, "unassessed")
  expect_identical(result$claims$PMID, "12345678")
  expect_identical(result$provider$name, "fake")
})

test_that("provider payload contains aggregate context only", {
  interpret <- ai_function("InterMetaAIInterpret")
  provider_constructor <- ai_function("InterMetaAIProvider")
  observed <- new.env(parent = emptyenv())
  provider <- provider_constructor(function(messages, schema) {
    observed$messages <- messages
    list(pathways = list(list(
      pathway = "U1", literature_consistency = "insufficient_evidence",
      summary = "No evidence conclusion is asserted.", claims = list()
    )))
  })
  interpret(make_ai_context(), make_ai_evidence(), provider)
  payload_context <- observed$messages[[2L]]$content$context
  expect_identical(
    names(payload_context),
    c(
      "analysis", "pathways", "interactions", "driver_gene_pairs",
      "warnings", "provenance"
    )
  )
  all_names <- unique(unlist(lapply(payload_context, function(x) {
    if (is.data.frame(x) || is.list(x)) names(x) else character()
  }), use.names = FALSE))
  expect_false(any(grepl(
    "patient|sample|rna|expression|selected.samples",
    all_names, ignore.case = TRUE
  )))
  expect_false(grepl(
    "PATIENT-SECRET",
    paste(capture.output(str(payload_context)), collapse = "\n"),
    fixed = TRUE
  ))
})

test_that("provider payload is deterministically bounded per pathway", {
  interpret <- ai_function("InterMetaAIInterpret")
  evidence <- make_ai_evidence()
  evidence$records <- do.call(rbind, lapply(seq_len(25L), function(i) {
    row <- evidence$records
    row$PMID <- sprintf("%08d", i)
    row$Title <- paste("Record", i)
    row
  }))
  # Duplicate the first PMID through another query; it must not consume a slot.
  evidence$records <- rbind(evidence$records, evidence$records[1L, ])
  observed <- new.env(parent = emptyenv())
  provider_constructor <- ai_function("InterMetaAIProvider")
  provider <- provider_constructor(function(messages, schema) {
    sent <- messages[[2L]]$content$evidence$records
    observed$pmids <- sent$PMID
    list(pathways = list(list(
      pathway = "U1", literature_consistency = "supported",
      summary = "Bounded evidence.",
      claims = list(list(claim = "First record.", pmids = sent$PMID[1L]))
    )))
  })
  result <- interpret(
    make_ai_context(), evidence, provider, max_records_per_pathway = 20L
  )
  expect_length(observed$pmids, 20L)
  expect_identical(observed$pmids, sprintf("%08d", 1:20))
  expect_true(result$record_selection$truncated)
  expect_identical(result$record_selection$sent_record_rows, 20L)
  expect_true(any(grepl("truncated", result$warnings, fixed = TRUE)))
})

test_that("invented PMID and invalid consistency are rejected", {
  interpret <- ai_function("InterMetaAIInterpret")
  expect_error(
    interpret(
      make_ai_context(), make_ai_evidence(),
      make_fake_provider("99999999"), language = "en"
    ),
    "absent from pathway-specific evidence\\$records"
  )

  provider_constructor <- ai_function("InterMetaAIProvider")
  invalid_provider <- provider_constructor(function(messages, schema) {
    list(pathways = list(list(
      pathway = "U1", literature_consistency = "proven",
      summary = "Overstated", claims = list()
    )))
  })
  expect_error(
    interpret(
      make_ai_context(), make_ai_evidence(), invalid_provider, language = "en"
    ),
    "Invalid literature_consistency"
  )
})

test_that("citations cannot borrow a PMID retrieved for another pathway", {
  interpret <- ai_function("InterMetaAIInterpret")
  provider_constructor <- ai_function("InterMetaAIProvider")
  provider <- provider_constructor(function(messages, schema) {
    list(pathways = list(list(
      pathway = "U2", literature_consistency = "supported",
      summary = "A record for U1 cannot support U2.",
      claims = list(list(claim = "Wrong pathway citation.", pmids = "12345678"))
    )))
  })
  expect_error(
    interpret(make_ai_context(), make_ai_evidence(), provider),
    "pathway-specific evidence\\$records: U2:12345678"
  )
})

test_that("non-insufficient classifications require pathway-specific evidence and citations", {
  interpret <- ai_function("InterMetaAIInterpret")
  provider_constructor <- ai_function("InterMetaAIProvider")
  for (classification in c("supported", "mixed", "potentially_contradicted")) {
    provider <- provider_constructor(function(messages, schema) {
      list(pathways = list(list(
        pathway = "U1", literature_consistency = classification,
        summary = "Uncited classification.", claims = list()
      )))
    })
    expect_error(
      interpret(make_ai_context(), make_ai_evidence(), provider),
      "requires at least one pathway-specific cited claim"
    )
    no_records <- make_ai_evidence()
    no_records$records <- no_records$records[0L, , drop = FALSE]
    expect_error(
      interpret(make_ai_context(), no_records, provider),
      "has no pathway-specific evidence records"
    )
  }
})

test_that("empty claim text and uncited claims are rejected", {
  interpret <- ai_function("InterMetaAIInterpret")
  provider_constructor <- ai_function("InterMetaAIProvider")
  provider_with_claim <- function(claim) {
    provider_constructor(function(messages, schema) {
      list(pathways = list(list(
        pathway = "U1", literature_consistency = "supported",
        summary = "Invalid claim.", claims = list(claim)
      )))
    })
  }
  expect_error(
    interpret(make_ai_context(), make_ai_evidence(),
              provider_with_claim(list(claim = "  ", pmids = "12345678"))),
    "non-empty text"
  )
  expect_error(
    interpret(make_ai_context(), make_ai_evidence(),
              provider_with_claim(list(claim = "Uncited.", pmids = character()))),
    "must cite at least one PMID"
  )
})

test_that("model evidence uses column allowlists and excludes error details", {
  interpret <- ai_function("InterMetaAIInterpret")
  provider_constructor <- ai_function("InterMetaAIProvider")
  evidence <- make_ai_evidence()
  for (table_name in c("queries", "searches", "records", "response_hashes")) {
    evidence[[table_name]]$Patient.ID <- "EVIDENCE-PATIENT-SECRET"
    evidence[[table_name]]$api_key <- "EVIDENCE-API-SECRET"
    evidence[[table_name]]$error <- "EVIDENCE-ERROR-SECRET"
  }
  evidence$api_key <- "TOP-LEVEL-SECRET"
  evidence$provenance$api_key <- "PROVENANCE-API-SECRET"
  evidence$provenance$Patient.ID <- "PROVENANCE-PATIENT-SECRET"
  observed <- new.env(parent = emptyenv())
  provider <- provider_constructor(function(messages, schema) {
    observed$payload <- messages[[2L]]$content$evidence
    list(pathways = list(list(
      pathway = "U1", literature_consistency = "supported",
      summary = "Allowed evidence only.",
      claims = list(list(claim = "U1 evidence.", pmids = "12345678"))
    )))
  })
  result <- interpret(make_ai_context(), evidence, provider)
  payload <- observed$payload
  for (table_name in c("queries", "searches", "records")) {
    expect_false(any(c("Patient.ID", "api_key", "error") %in%
                       names(payload[[table_name]])))
  }
  rendered <- paste(capture.output(str(payload)), collapse = "\n")
  expect_false(grepl("SECRET", rendered, fixed = TRUE))
  expect_identical(result$claims$PMID, "12345678")
})

test_that("model payload bounds long abstracts and aggregate context tables", {
  interpret <- ai_function("InterMetaAIInterpret")
  provider_constructor <- ai_function("InterMetaAIProvider")
  context <- make_ai_context()
  context$interactions <- context$interactions[rep(1L, 250L), , drop = FALSE]
  context$interactions$Shared.MGP <- paste(rep("LARGE-CONCATENATED-PAIR", 100L),
                                         collapse = "|")
  context$driver_gene_pairs <- context$driver_gene_pairs[
    rep(1L, 250L), , drop = FALSE
  ]
  context$pathways$Synonyms[1L] <- paste(rep("S", 2500L), collapse = "")
  evidence <- make_ai_evidence()
  evidence$records$Abstract <- paste(rep("A", 5000L), collapse = "")
  observed <- new.env(parent = emptyenv())
  provider <- provider_constructor(function(messages, schema) {
    observed$payload <- messages[[2L]]$content
    list(pathways = list(list(
      pathway = "U1", literature_consistency = "supported",
      summary = "Bounded payload.",
      claims = list(list(claim = "U1 evidence.", pmids = "12345678"))
    )))
  })
  result <- interpret(context, evidence, provider)
  expect_identical(nchar(observed$payload$evidence$records$Abstract), 4000L)
  expect_equal(nrow(observed$payload$context$interactions), 200L)
  expect_equal(nrow(observed$payload$context$driver_gene_pairs), 200L)
  expect_false("Shared.MGP" %in% names(observed$payload$context$interactions))
  expect_identical(nchar(observed$payload$context$pathways$Synonyms[1L]), 2000L)
  expect_identical(result$context_selection$text_cells_truncated, 1L)
  expect_identical(result$record_selection$abstracts_truncated, 1L)
  expect_identical(result$record_selection$abstract_characters_before, 5000)
  expect_identical(result$record_selection$abstract_characters_sent, 4000)
  expect_true(result$context_selection$interactions_truncated)
  expect_true(result$context_selection$driver_gene_pairs_truncated)
  expect_true(any(grepl("abstract(s) were truncated", result$warnings,
                       fixed = TRUE)))
  expect_true(any(grepl("interactions were deterministically truncated",
                       result$warnings, fixed = TRUE)))
  expect_identical(nchar(evidence$records$Abstract), 5000L)
  expect_equal(nrow(context$interactions), 250L)
  expect_identical(nchar(context$pathways$Synonyms[1L]), 2500L)
})

test_that("aggregate abstract budget applies across multiple pathways", {
  interpret <- ai_function("InterMetaAIInterpret")
  provider_constructor <- ai_function("InterMetaAIProvider")
  evidence <- make_ai_evidence()
  evidence$records <- do.call(rbind, lapply(1:3, function(pathway_index) {
    rows <- evidence$records[rep(1L, 20L), , drop = FALSE]
    rows$pathway <- paste0("U", pathway_index)
    rows$PMID <- as.character(pathway_index * 100L + 1:20)
    rows$Abstract <- paste(rep("A", 5000L), collapse = "")
    rows
  }))
  observed <- new.env(parent = emptyenv())
  provider <- provider_constructor(function(messages, schema) {
    observed$records <- messages[[2L]]$content$evidence$records
    list(pathways = list(list(
      pathway = "U1", literature_consistency = "supported",
      summary = "Bounded multi-pathway evidence.",
      claims = list(list(claim = "U1 evidence.", pmids = "101"))
    )))
  })
  result <- interpret(make_ai_context(), evidence, provider)
  lengths <- nchar(observed$records$Abstract)
  expect_equal(nrow(observed$records), 60L)
  expect_true(all(lengths == floor(200000L / 60L)))
  expect_lte(sum(lengths), 200000L)
  expect_identical(result$record_selection$abstracts_truncated, 60L)
  expect_identical(result$record_selection$abstract_budget_chars, 200000L)
  expect_equal(result$record_selection$abstract_characters_before, 300000)
  expect_equal(result$record_selection$abstract_characters_sent, sum(lengths))
  expect_true(all(nchar(evidence$records$Abstract) == 5000L))
})

test_that("round-robin selection represents every evidence layer before repeats", {
  interpret <- ai_function("InterMetaAIInterpret")
  provider_constructor <- ai_function("InterMetaAIProvider")
  evidence <- make_ai_evidence()
  layers <- c("review", "prognosis", "gene_bridge", "biology")
  evidence$records <- do.call(rbind, lapply(seq_along(layers), function(i) {
    rows <- evidence$records[rep(1L, 3L), , drop = FALSE]
    rows$layer <- layers[i]
    rows$PMID <- as.character(i * 100L + 1:3)
    rows
  }))
  observed <- new.env(parent = emptyenv())
  provider <- provider_constructor(function(messages, schema) {
    observed$records <- messages[[2L]]$content$evidence$records
    list(pathways = list(list(
      pathway = "U1", literature_consistency = "supported",
      summary = "Evidence from all layers.",
      claims = list(list(claim = "Selected record.",
                         pmids = observed$records$PMID[1L]))
    )))
  })
  result <- interpret(make_ai_context(), evidence, provider,
                      max_records_per_pathway = 5L)
  expect_identical(observed$records$layer,
                   c("prognosis", "biology", "review", "gene_bridge", "prognosis"))
  expect_identical(observed$records$PMID, c("201", "401", "101", "301", "202"))
  expect_identical(result$record_selection$selection_method,
                   "pathway_order_layer_round_robin_v1")
  per_layer <- result$record_selection$per_layer
  expect_equal(per_layer$Sent[match(c("prognosis", "biology", "review", "gene_bridge"),
                                   per_layer$Layer)], c(2, 1, 1, 1))
})

test_that("report retains evidence schema and retrieval audit metadata", {
  report_builder <- ai_function("InterMetaAIReport")
  evidence <- make_ai_evidence()
  report <- report_builder(make_ai_context(), evidence = evidence)
  expect_identical(report$evidence$schema_version, evidence$schema_version)
  expect_identical(report$evidence$generated_at, evidence$generated_at)
  expect_identical(report$evidence$search_date, evidence$search_date)
  expect_identical(report$evidence$response_hashes, evidence$response_hashes)
  expect_identical(report$evidence$provenance, evidence$provenance)
})

test_that("truncated PubMed searches surface in report warnings", {
  report_builder <- ai_function("InterMetaAIReport")
  evidence <- make_ai_evidence()
  evidence$searches$truncated <- TRUE
  evidence$searches$total_count <- 100L
  report <- report_builder(make_ai_context(), evidence = evidence)
  expect_true(any(grepl("truncat", report$warnings, ignore.case = TRUE)))
})

test_that("report has deterministic Top 5 tables and works without a provider", {
  report_builder <- ai_function("InterMetaAIReport")
  context <- make_ai_context()
  report <- report_builder(context, evidence = make_ai_evidence(), language = "zh")
  expect_s3_class(report, "InterMetaAIReport")
  expect_equal(nrow(report$top_pathways), 10L)
  expect_setequal(report$top_pathways$Pathway, c(paste0("U", 1:5), paste0("F", 1:5)))
  expect_equal(
    report$top_pathways$Difference,
    context$pathways$Difference[context$pathways$Direction.Rank <= 5]
  )
  expect_equal(nrow(report$interpretation$claims), 0L)
  expect_true(any(grepl("evidence skeleton", report$warnings, fixed = TRUE)))
  expect_length(report$notices, 3L)
})

test_that("report publishes a complete immutable artifact directory", {
  report_builder <- ai_function("InterMetaAIReport")
  root <- tempfile("intermeta-ai-report-test-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  outdir <- file.path(root, "report")

  report <- report_builder(
    make_ai_context(), evidence = make_ai_evidence(),
    provider = make_fake_provider(), language = "en", outdir = outdir
  )
  expected <- c(
    "report.md", "report.rds", "evidence.tsv", "claims.tsv", "assessments.tsv", "queries.tsv",
    "pathway_synthesis.tsv", "biology_findings.tsv", "exploration_summary.tsv", "exploration_findings.tsv",
    "searches.tsv", "response_hashes.tsv", "warnings.txt", "provenance.rds",
    "output_manifest.txt",
    "report_complete.rds"
  )
  expect_setequal(list.files(outdir), expected)
  expect_length(list.files(outdir), 16L)
  completion <- readRDS(file.path(outdir, "report_complete.rds"))
  expect_identical(completion$status, "complete")

  manifest <- utils::read.delim(
    file.path(outdir, "output_manifest.txt"), stringsAsFactors = FALSE
  )
  expect_setequal(
    manifest$File,
    c(
      "report.md", "report.rds", "evidence.tsv", "claims.tsv", "assessments.tsv", "queries.tsv",
      "pathway_synthesis.tsv", "biology_findings.tsv", "exploration_summary.tsv", "exploration_findings.tsv",
      "searches.tsv", "response_hashes.tsv", "warnings.txt", "provenance.rds",
      "report_complete.rds"
    )
  )
  hashes <- vapply(
    file.path(outdir, manifest$File),
    digest::digest, character(1), algo = "sha256", file = TRUE
  )
  expect_identical(unname(hashes), manifest$SHA256)
  expect_equal(as.numeric(file.info(file.path(outdir, manifest$File))$size),
               manifest$Bytes)
  expect_setequal(completion$content_files,
                  setdiff(expected, c("output_manifest.txt", "report_complete.rds")))
  stored_searches <- utils::read.delim(
    file.path(outdir, "searches.tsv"), stringsAsFactors = FALSE,
    colClasses = c(error = "character")
  )
  expected_searches <- make_ai_evidence()$searches
  expect_setequal(names(stored_searches), names(expected_searches))
  expect_equal(stored_searches[, names(expected_searches)], expected_searches)
  expect_equal(
    utils::read.delim(file.path(outdir, "response_hashes.tsv"),
                     stringsAsFactors = FALSE),
    make_ai_evidence()$response_hashes
  )
  expect_error(
    report_builder(
      make_ai_context(), evidence = make_ai_evidence(), outdir = outdir
    ),
    "already exists"
  )
  expect_s3_class(report, "InterMetaAIReport")
})

test_that("patient-level objects never enter reports or provenance", {
  report_builder <- ai_function("InterMetaAIReport")
  context <- make_ai_context()
  report <- report_builder(context, evidence = make_ai_evidence())
  rendered <- paste(capture.output(str(report)), collapse = "\n")
  expect_false(grepl("PATIENT-SECRET", rendered, fixed = TRUE))
  expect_false(any(grepl("patient|sample|rna|expression", names(report),
                         ignore.case = TRUE)))
  expect_false("Patient.ID" %in% names(report$top_pathways))

  root <- tempfile("intermeta-ai-privacy-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  outdir <- file.path(root, "report")
  report_builder(context, evidence = make_ai_evidence(), outdir = outdir)
  stored_report <- readRDS(file.path(outdir, "report.rds"))
  stored_provenance <- readRDS(file.path(outdir, "provenance.rds"))
  stored <- paste(
    capture.output(str(stored_report)),
    capture.output(str(stored_provenance)),
    collapse = "\n"
  )
  expect_false(grepl("PATIENT-SECRET", stored, fixed = TRUE))
})

test_that("free text cannot bypass structured citation validation", {
  context <- make_ai_context()
  evidence <- make_ai_evidence()
  evidence$records <- evidence$records[FALSE, , drop = FALSE]
  provider <- InterMetaAIProvider(function(messages, schema) {
    list(pathways = list(list(pathway = "U1",
      literature_consistency = "insufficient_evidence",
      summary = "PMID: 99999999 shows U1 causes this cancer.", claims = list())))
  })
  expect_error(clinical_ai_interpret(context, evidence, provider), "Inline citation")

  provider <- InterMetaAIProvider(function(messages, schema) {
    list(pathways = list(list(pathway = "U1",
      literature_consistency = "insufficient_evidence",
      summary = "U1 causes this cancer.", claims = list())))
  })
  result <- clinical_ai_interpret(context, evidence, provider)
  expect_match(result$pathways$Summary, "evidence is insufficient")
  expect_false(grepl("causes", result$pathways$Summary))
})

test_that("safe provider call metadata is retained without arbitrary fields", {
  context <- make_ai_context()
  evidence <- make_ai_evidence()
  provider <- make_fake_provider()
  original_generate <- provider$generate
  provider$generate <- function(messages, schema) {
    result <- original_generate(messages, schema)
    attr(result, "intermeta_provider_metadata") <- list(
      endpoint = "https://api.deepseek.com/chat/completions",
      requested_model = "deepseek-v4-pro", response_model = "deepseek-v4-pro",
      response_id = "test-response", finish_reason = "stop", max_tokens = 16384L,
       request_sha256 = paste(rep("a", 64L), collapse = ""),
       response_sha256 = paste(rep("b", 64L), collapse = ""),
       replay = TRUE,
       replay_candidate_sha256 = paste(rep("c", 64L), collapse = ""),
       replay_candidate_format_version = 1L,
       replay_candidate_validation_revision = "legacy_missing",
      usage = list(prompt_tokens = 12L, completion_tokens = 34L,
                   total_tokens = 46L, secret = "DO_NOT_STORE"),
      api_key = "DO_NOT_STORE", raw_response = "DO_NOT_STORE",
      reasoning_content = "DO_NOT_STORE"
    )
    result
  }
  interpretation <- clinical_ai_interpret(context, evidence, provider)
  report <- clinical_ai_report(context, evidence, interpretation = interpretation)
  expect_identical(report$provenance$provider_call$response_id, "test-response")
  expect_identical(report$provenance$provider_call$usage$total_tokens, 46L)
  expect_true(report$provenance$provider_call$replay)
  expect_identical(report$provenance$provider_call$replay_candidate_sha256,
    paste(rep("c", 64L), collapse = ""))
  expect_identical(report$provenance$provider_call$replay_candidate_format_version, 1L)
  expect_identical(report$provenance$provider_call$replay_candidate_validation_revision,
    "legacy_missing")
  expect_identical(report$interpretation$provider_call, interpretation$provider_call)
  expect_false(grepl("DO_NOT_STORE", paste(capture.output(dput(report)), collapse = "")))
})
