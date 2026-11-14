explore_fn <- function(name) getFromNamespace(name, "InterMeta")

explore_fixture <- function() {
  context <- structure(list(analysis = list(disease = "prostate adenocarcinoma", endpoint = "BCR", time_unit = "months"),
    pathways = data.frame(Pathway = c("P1", "P2"), Direction = "Unfavorable", Direction.Rank = 1:2,
      Difference = c(-7.25, -3.5)), interactions = data.frame(), driver_gene_pairs = data.frame(),
    warnings = character(), provenance = list()), class = c("InterMetaAIContext", "list"))
  evidence <- explore_fn(".InterMetaAIEmptyEvidence")()
  evidence$records <- data.frame(pathway = c("P1", "P1", "P2"), PMID = c("12345678", "12345679", "12345680"),
    layer = c("biology", "biology_process", "biology"), PublicationType = "Journal Article",
    Title = c("An experimental study", "A null observation", "Other pathway"),
    Abstract = c("Pathway inhibition reduced growth in prostate cancer cells. Tumour profiles differed across treatment contexts.",
      "Patient tissue showed no association between the component and tumour stage.", "An unrelated finding."))
  finding <- list(finding_id = "f1", pmid = "12345678", study_disease = "prostate cancer",
    disease_match = "matched", scope = "component", evidence_types = c("cell_experiment", "omics"),
    biological_processes = c("proliferation", "treatment_response"),
    observation = "Component inhibition reduced growth in the studied cells.", effect = "suppresses",
    study_context = "Prostate cancer cell experiment and expression profiling.",
    abstract_quote = "Pathway inhibition reduced growth in prostate cancer cells.", cohort_id = "unknown")
  entry <- list(pathway = "P1", exploration = list(biological_summary = "The component participates in context-dependent growth and treatment adaptation.",
    context_dependence = "The experimental setting differs from patient tissue observations.",
    research_hypothesis = "Test whether component perturbation alters treatment sensitivity.", findings = list(finding)))
  list(context = context, evidence = evidence, raw = list(pathways = list(entry)))
}

explore_run <- function(fixture, ...) InterMeta::InterMetaAIInterpret(fixture$context, fixture$evidence,
  InterMeta::InterMetaAIProvider(function(messages, schema) fixture$raw), ...)

test_that("Chinese reports do not lose exact target-disease aliases", {
  f <- explore_fixture()
  for (label in c("前列腺癌", "前列腺腺癌")) {
    f$raw$pathways[[1L]]$exploration$findings[[1L]]$study_disease <- label
    out <- explore_run(f, language = "zh")
    expect_true(out$exploration_findings$Disease.Match.Verified)
    expect_identical(out$exploration_summary$Exploration.Status[1], "provided")
  }
  for (label in c("肺癌", "前列腺癌和膀胱癌", "癌症")) {
    f$raw$pathways[[1L]]$exploration$findings[[1L]]$study_disease <- label
    expect_false(explore_run(f, language = "zh")$exploration_findings$Disease.Match.Verified)
  }
  key <- explore_fn(".InterMetaAIExplorationDiseaseKey")
  expect_true(nzchar(key("肝癌")))
  expect_false(identical(key("乳腺癌 HER2"), key("胃癌 HER2")))
  expect_identical(key("  肝癌  "), key("肝癌"))
  expect_identical(key("未知"), "")
})

test_that("fresh interpretation uses exploration schema in one provider request", {
  f <- explore_fixture()
  calls <- 0L
  result <- InterMeta::InterMetaAIInterpret(f$context, f$evidence,
    InterMeta::InterMetaAIProvider(function(messages, schema) {
      calls <<- calls + 1L
      expect_identical(schema$properties$pathways$items$required, c("pathway", "exploration"))
      expect_false("assessments" %in% names(schema$properties$pathways$items$properties))
      expect_identical(
        schema$properties$pathways$items$properties$exploration$properties$findings$maxItems,
        5L
      )
      quote_schema <- schema$properties$pathways$items$properties$exploration$properties$findings$
        items$properties$abstract_quote
      expect_match(quote_schema$description, "verbatim contiguous excerpt", fixed = TRUE)
      expect_match(quote_schema$description, "original language", fixed = TRUE)
      expect_match(messages[[1L]]$content, "NOT a literature admission criterion", fixed = TRUE)
      expect_match(messages[[1L]]$content, "no more than 5 non-redundant", fixed = TRUE)
      expect_match(messages[[1L]]$content, "original source language", fixed = TRUE)
      expect_match(messages[[1L]]$content, "Never translate, paraphrase", fixed = TRUE)
      expect_match(messages[[1L]]$content, "shorten with ellipses", fixed = TRUE)
      expect_identical(messages[[2L]]$content$purpose, "biology_exploration")
      f$raw
    }))
  expect_identical(calls, 1L)
  expect_identical(result$schema_version, "4.3.0")
  expect_identical(result$purpose, "biology_exploration")
  expect_identical(result$record_selection$purpose, "biology_exploration")
  expect_identical(result$record_selection$max_records_per_pathway, 30L)
  expect_identical(result$exploration_summary$Exploration.Status, c("provided", "not_assessed"))
  expect_identical(result$exploration_summary$Quote.Grounding,
    c("complete", "no_findings"))
  expect_identical(result$exploration_summary$Evidence.Eligibility,
    c("complete", "no_findings"))
  expect_identical(result$exploration_summary$Difference, f$context$pathways$Difference)
  expect_equal(nrow(result$assessments), 0)
})

test_that("fresh provider output enforces the five-finding boundary locally", {
  f <- explore_fixture()
  item <- f$raw$pathways[[1L]]$exploration$findings[[1L]]
  f$raw$pathways[[1L]]$exploration$findings <- lapply(seq_len(6L), function(i) {
    out <- item
    out$finding_id <- paste0("f", i)
    out
  })
  expect_error(explore_run(f), "exceeds 5 findings", fixed = TRUE)

  normalized <- explore_run(explore_fixture())
  # Normalized stored interpretations follow their original recorded selection;
  # the fresh-provider boundary must not relabel them as raw schema responses.
  expect_s3_class(InterMeta::InterMetaAIReport(
    explore_fixture()$context, explore_fixture()$evidence, normalized
  ), "InterMetaAIReport")

  f <- explore_fixture()
  f$raw <- list(exploration_summary = data.frame(), exploration_findings = data.frame())
  expect_error(explore_run(f), "must use the pathways schema", fixed = TRUE)
})

test_that("fresh provider keeps the supported single-finding object form", {
  f <- explore_fixture()
  f$raw$pathways[[1L]]$exploration$findings <-
    f$raw$pathways[[1L]]$exploration$findings[[1L]]
  result <- explore_run(f)
  expect_identical(nrow(result$exploration_findings), 1L)
  expect_identical(result$exploration_findings$Finding.ID, "f1")
  expect_identical(result$exploration_findings$PMID, "12345678")
})

test_that("report and installed metadata expose the same AI revisions", {
  revision <- explore_fn(".InterMetaAIContractRevision")
  validation_revision <- explore_fn(".InterMetaAIValidationRevision")
  description <- read.dcf(file.path(getNamespaceInfo(asNamespace("InterMeta"), "path"),
    "DESCRIPTION"))
  expect_identical(unname(description[1L, "Config/InterMeta/AIContractRevision"]),
    revision)
  expect_identical(unname(description[1L, "Config/InterMeta/AIValidationRevision"]),
    validation_revision)
  f <- explore_fixture()
  report <- InterMeta::InterMetaAIReport(f$context, f$evidence, explore_run(f))
  expect_identical(report$schema_version, "4.3.0")
  expect_identical(report$provenance$ai_contract_revision, revision)
  expect_identical(report$provenance$ai_validation_revision, validation_revision)
})

test_that("nonexclusive channels count component experiments without an endpoint gate", {
  f <- explore_fixture()
  result <- explore_run(f)
  row <- result$exploration_summary[1L, ]
  expect_identical(row$Functional.PMIDs, "12345678")
  expect_identical(row$Omics.PMIDs, "12345678")
  expect_identical(row$Patient.PMIDs, "")
  expect_identical(row$Unique.PMID.Count, 1L)
  expect_true(result$exploration_findings$Eligible.For.Summary)
  expect_identical(result$exploration_findings$Quote.Match.Method, "whitespace_exact")
  expect_identical(result$exploration_findings$Scope, "component")
  f$context$analysis$endpoint <- "unknown"
  expect_identical(explore_run(f)$exploration_summary, result$exploration_summary)
})

test_that("multiple findings per publication retain null observations without double counting", {
  f <- explore_fixture()
  second <- f$raw$pathways[[1L]]$exploration$findings[[1L]]
  second$finding_id <- "f2"
  second$evidence_types <- "bioinformatics"
  second$effect <- "no_association"
  second$observation <- "A distinct null finding was observed in another context."
  f$raw$pathways[[1L]]$exploration$findings[[2L]] <- second
  result <- explore_run(f)
  expect_equal(nrow(result$exploration_findings), 2L)
  expect_identical(result$exploration_summary$Unique.PMID.Count[1L], 1L)
  expect_identical(result$exploration_summary$Omics.PMIDs[1L], "12345678")
  expect_match(result$exploration_summary$Source.Linked.Observations[1L], "null finding")
  expect_true("no_association" %in% result$exploration_findings$Effect)
  f$raw$pathways[[1L]]$exploration$findings[[2L]]$finding_id <- "f1"
  expect_error(explore_run(f), "Duplicate exploration finding_id")
})

test_that("patient observations and reviews remain distinct from clinical interventions", {
  f <- explore_fixture()
  item <- f$raw$pathways[[1L]]$exploration$findings[[1L]]
  item$evidence_types <- "human_tissue"
  item$pmid <- "12345679"
  item$abstract_quote <- f$evidence$records$Abstract[2L]
  item$effect <- "no_association"
  f$raw$pathways[[1L]]$exploration$findings <- list(item)
  expect_identical(explore_run(f)$exploration_summary$Patient.PMIDs[1L], "12345679")
  f$raw$pathways[[1L]]$exploration$findings[[1L]]$evidence_types <- "review"
  result <- explore_run(f)
  expect_identical(result$exploration_summary$Review.PMIDs[1L], "12345679")
  expect_identical(result$exploration_summary$Unique.PMID.Count[1L], 0L)
  expect_identical(result$exploration_summary$Functional.PMIDs[1L], "")
  f$raw$pathways[[1L]]$exploration$findings[[1L]]$evidence_types <- c("review", "cell_experiment")
  expect_error(explore_run(f), "Review findings cannot also count")
})

test_that("pathway membership fails closed while quote failures are quarantined", {
  f <- explore_fixture()
  f$raw$pathways[[1L]]$exploration$findings[[1L]]$pmid <- "99999999"
  expect_error(explore_run(f), "pathway/PMID absent")
  f$raw$pathways[[1L]]$exploration$findings[[1L]]$pmid <- "12345680"
  expect_error(explore_run(f), "pathway/PMID absent")
  f <- explore_fixture()
  invalid_quote <- "An experimental study title that is not abstract text."
  f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <- invalid_quote
  result <- explore_run(f)
  expect_identical(result$exploration_findings$Abstract.Quote, invalid_quote)
  expect_false(result$exploration_findings$Quote.Verified)
  expect_false(result$exploration_findings$Eligible.For.Summary)
  expect_identical(result$exploration_findings$Quote.Match.Method, "source_mismatch")
  expect_identical(result$exploration_summary$Exploration.Status[1L], "ungrounded_findings")
  expect_identical(result$exploration_summary$Quote.Grounding[1L], "none")
  expect_identical(result$exploration_summary$Evidence.Eligibility[1L], "none")
  expect_identical(result$exploration_summary$Unique.PMID.Count[1L], 0L)
  mismatch_warning <- result$warnings[grepl("were quarantined", result$warnings, fixed = TRUE)]
  expect_length(mismatch_warning, 1L)
  expect_match(mismatch_warning, 'pathway="P1"', fixed = TRUE)
  expect_match(mismatch_warning, 'finding="f1"', fixed = TRUE)
  expect_match(mismatch_warning, "PMID=12345678", fixed = TRUE)
  expect_false(grepl(invalid_quote, mismatch_warning, fixed = TRUE))
  f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <- ""
  result <- explore_run(f)
  expect_identical(result$exploration_summary$Exploration.Status[1L], "ungrounded_findings")
  expect_identical(result$exploration_summary$Unique.PMID.Count[1L], 0L)
  expect_false(result$exploration_findings$Quote.Verified)
  expect_identical(result$exploration_findings$Quote.Match.Method, "empty")

  for (uninformative in c("a", ".", "-- --")) {
    f <- explore_fixture()
    f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <- uninformative
    result <- explore_run(f)
    expect_false(result$exploration_findings$Quote.Verified)
    expect_false(result$exploration_findings$Eligible.For.Summary)
    expect_identical(result$exploration_findings$Quote.Match.Method,
      "insufficient_information")
    expect_match(result$warnings[grepl("were quarantined", result$warnings, fixed = TRUE)],
      "reason=insufficient_information", fixed = TRUE)
  }

  for (size in c(19L, 20L)) {
    f <- explore_fixture()
    boundary <- paste(rep("a", size), collapse = "")
    f$evidence$records$Abstract[1L] <- boundary
    f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <- boundary
    result <- explore_run(f)
    expect_identical(result$exploration_findings$Quote.Verified, size == 20L)
    expect_identical(result$exploration_findings$Quote.Match.Method,
      if (size == 20L) "whitespace_exact" else "insufficient_information")
  }

  f <- explore_fixture()
  f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <- "\u00a0\u2009\u00ad"
  result <- explore_run(f)
  expect_identical(result$exploration_findings$Abstract.Quote, "")
  expect_identical(result$exploration_findings$Quote.Match.Method, "empty")
})

test_that("overlength quotes are preserved and quarantined per row", {
  f <- explore_fixture()
  overlength <- paste(rep("z", 1001L), collapse = "")
  # Even an otherwise exact source excerpt cannot bypass the provider contract
  # bound; retaining it makes the failed row auditable without losing the report.
  f$evidence$records$Abstract[1L] <- overlength
  f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <- overlength
  result <- explore_run(f)
  expect_identical(result$exploration_findings$Abstract.Quote, overlength)
  expect_false(result$exploration_findings$Quote.Verified)
  expect_false(result$exploration_findings$Eligible.For.Summary)
  expect_identical(result$exploration_findings$Quote.Match.Method, "overlength")
  expect_identical(result$exploration_summary$Quote.Grounding[1L], "none")
  expect_identical(result$exploration_summary$Evidence.Eligibility[1L], "none")
  expect_match(result$warnings[grepl("were quarantined", result$warnings, fixed = TRUE)],
    "reason=overlength", fixed = TRUE)
})

test_that("explicit Unicode presentation variants retain lexical source matching", {
  whitespace <- explore_fixture()
  whitespace$evidence$records$Abstract[1L] <-
    "Pathway inhibition\n reduced\t growth in prostate cancer cells."
  whitespace_result <- explore_run(whitespace)
  expect_true(whitespace_result$exploration_findings$Quote.Verified)
  expect_identical(whitespace_result$exploration_findings$Quote.Match.Method,
    "whitespace_exact")

  cases <- list(
    unicode_spaces = c("alpha\u00a0beta\u2009gamma\u202fdelta", "alpha beta gamma delta"),
    curly_single_quotes = c("the pathway\u2019s role", "the pathway's role"),
    curly_double_quotes = c("\u201cpathways\u201d", '"pathways"'),
    hyphen_and_dashes = c("alpha\u2010beta\u2011gamma\u2012delta\u2013epsilon\u2014zeta",
      "alpha-beta-gamma-delta-epsilon-zeta"),
    soft_hyphen = c("meta\u00adbolic", "metabolic"),
    ellipsis = c("growth\u2026continued", "growth...continued")
  )
  for (label in names(cases)) {
    f <- explore_fixture()
    f$evidence$records$Abstract[1L] <- paste(
      "A sufficiently informative prefix", cases[[label]][1L], "and source suffix.")
    f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <- paste(
      "A sufficiently informative prefix", cases[[label]][2L], "and source suffix.")
    result <- explore_run(f)
    expect_true(result$exploration_findings$Quote.Verified, info = label)
    expect_identical(result$exploration_findings$Quote.Match.Method,
      "typography_equivalent", info = label)
  }
})

test_that("semantic and scientific-symbol variants are quarantined rather than fuzzily matched", {
  failures <- list(
    case_change = c("Pathway inhibition reduced growth.", "pathway inhibition reduced growth."),
    paraphrase = c("Pathway inhibition reduced growth.", "Inhibiting the pathway decreased growth."),
    greek_name = c("beta oxidation differed.", "\u03b2 oxidation differed."),
    mathematical_minus = c("The coefficient was \u22121.", "The coefficient was -1."),
    inequality = c("The result was P < 0.05.", "The result was P \u2264 0.05."),
    composed_accent = c("The caf\u00e9 signal differed.", "The cafe\u0301 signal differed."),
    zero_width = c("Metabolic adaptation differed.", "Metabolic\u200badaptation differed."),
    ellipsis_wildcard = c("Pathway inhibition strongly reduced growth.",
      "Pathway inhibition\u2026reduced growth.")
  )
  for (label in names(failures)) {
    f <- explore_fixture()
    f$evidence$records$Abstract[1L] <- paste(
      "A sufficiently informative source says", failures[[label]][1L], "in this experiment.")
    f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <- paste(
      "A sufficiently informative source says", failures[[label]][2L], "in this experiment.")
    result <- explore_run(f)
    expect_false(result$exploration_findings$Quote.Verified, info = label)
    expect_false(result$exploration_findings$Eligible.For.Summary, info = label)
    expect_identical(result$exploration_findings$Quote.Match.Method,
      "source_mismatch", info = label)
  }
})

test_that("non-contiguous and clipped source text remains quarantined", {
  cases <- list(
    skipped_sentence = c(
      "First exact source sentence. A middle source sentence must not be skipped. Final exact source sentence.",
      "First exact source sentence. Final exact source sentence."
    ),
    clipped_clause = c(
      "The effect was present in one cell line, but not in the comparison subline.",
      "The effect was present in one cell line."
    ),
    clipped_parenthesis = c(
      "The effect was observed in parental control cells (PC-3 and DU-145).",
      "The effect was observed in parental control cells."
    )
  )
  for (label in names(cases)) {
    f <- explore_fixture()
    f$evidence$records$Abstract[1L] <- cases[[label]][1L]
    f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <- cases[[label]][2L]
    result <- explore_run(f)
    expect_false(result$exploration_findings$Quote.Verified, info = label)
    expect_identical(result$exploration_findings$Quote.Match.Method,
      "source_mismatch", info = label)
  }
})

test_that("quote grounding is restricted to the provider-visible truncated abstract", {
  f <- explore_fixture()
  f$evidence$records$Abstract[1L] <- paste0(
    "Verified prefix phrase. ", paste(rep("x", 4100L), collapse = ""),
    " Hidden source suffix after the provider limit."
  )
  f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <-
    "Hidden source suffix after the provider limit."
  result <- explore_run(f)
  expect_false(result$exploration_findings$Quote.Verified)
  expect_identical(result$exploration_findings$Quote.Match.Method, "source_mismatch")

  f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <-
    "Verified prefix phrase."
  result <- explore_run(f)
  expect_true(result$exploration_findings$Quote.Verified)
  expect_identical(result$exploration_findings$Quote.Match.Method, "whitespace_exact")
})

test_that("mixed quote grounding keeps only verified findings in evidence summaries", {
  f <- explore_fixture()
  f$raw$pathways[[1L]]$exploration$biological_summary <-
    "MODEL-LEVEL SUMMARY THAT MAY USE THE REJECTED FINDING."
  f$raw$pathways[[1L]]$exploration$context_dependence <-
    "MODEL-LEVEL CONTEXT THAT MAY USE THE REJECTED FINDING."
  f$raw$pathways[[1L]]$exploration$research_hypothesis <-
    "MODEL-LEVEL HYPOTHESIS THAT MAY USE THE REJECTED FINDING."
  rejected <- f$raw$pathways[[1L]]$exploration$findings[[1L]]
  rejected$finding_id <- "f2"
  rejected$pmid <- "12345679"
  rejected$evidence_types <- "human_tissue"
  rejected$biological_processes <- "disease_development"
  rejected$observation <- "This rejected observation must not enter source-linked summaries."
  rejected$abstract_quote <- "Patient tissue showed an altered association."
  f$raw$pathways[[1L]]$exploration$findings[[2L]] <- rejected
  result <- explore_run(f)
  rows <- result$exploration_findings
  expect_identical(rows$Quote.Match.Method, c("whitespace_exact", "source_mismatch"))
  expect_identical(rows$Eligible.For.Summary, c(TRUE, FALSE))
  expect_identical(result$exploration_summary$Exploration.Status[1L], "provided")
  expect_identical(result$exploration_summary$Quote.Grounding[1L], "partial")
  expect_identical(result$exploration_summary$Evidence.Eligibility[1L], "partial")
  expect_identical(result$exploration_summary$Key.PMIDs[1L], "12345678")
  expect_identical(result$exploration_summary$Unique.PMID.Count[1L], 1L)
  expect_false(grepl("rejected observation",
    result$exploration_summary$Source.Linked.Observations[1L], fixed = TRUE))
  expect_match(result$exploration_summary$Biological.Summary[1L],
    "Component inhibition reduced growth", fixed = TRUE)
  expect_false(grepl("MODEL-LEVEL", result$exploration_summary$Biological.Summary[1L],
    fixed = TRUE))
  expect_match(result$exploration_summary$Context.Dependence[1L], "Suppressed", fixed = TRUE)
  expect_match(result$exploration_summary$Research.Hypothesis[1L], "Suppressed", fixed = TRUE)
  expect_identical(result$exploration_summary$Model.Biological.Summary[1L],
    "MODEL-LEVEL SUMMARY THAT MAY USE THE REJECTED FINDING.")
  expect_identical(result$exploration_summary$Model.Context.Dependence[1L],
    "MODEL-LEVEL CONTEXT THAT MAY USE THE REJECTED FINDING.")
  expect_identical(result$exploration_summary$Model.Research.Hypothesis[1L],
    "MODEL-LEVEL HYPOTHESIS THAT MAY USE THE REJECTED FINDING.")
})

test_that("exact quotes with mismatched disease are separately excluded from target evidence", {
  f <- explore_fixture()
  f$raw$pathways[[1L]]$exploration$findings[[1L]]$study_disease <- "breast cancer"
  result <- explore_run(f)
  expect_true(result$exploration_findings$Quote.Verified)
  expect_false(result$exploration_findings$Disease.Match.Verified)
  expect_false(result$exploration_findings$Eligible.For.Summary)
  expect_identical(result$exploration_summary$Quote.Grounding[1L], "complete")
  expect_identical(result$exploration_summary$Evidence.Eligibility[1L], "none")
  expect_identical(result$exploration_summary$Key.PMIDs[1L], "")
  expect_identical(result$exploration_summary$Model.Biological.Summary[1L],
    f$raw$pathways[[1L]]$exploration$biological_summary)
  expect_false(identical(result$exploration_summary$Biological.Summary[1L],
    result$exploration_summary$Model.Biological.Summary[1L]))
  expect_identical(result$exploration_summary$Context.Dependence[1L], "")
  expect_identical(result$exploration_summary$Research.Hypothesis[1L], "")
})

test_that("missing assessment, empty relevant set and failed retrieval are distinguished", {
  f <- explore_fixture()
  f$raw$pathways[[1L]]$exploration$findings <- list()
  result <- explore_run(f)
  expect_identical(result$exploration_summary$Exploration.Status, c("no_relevant_findings", "not_assessed"))
  f$evidence$searches <- data.frame(pathway = "P1", api_status = "error")
  result <- explore_run(f)
  expect_identical(result$exploration_summary$Exploration.Status[1L], "input_failed")
  expect_match(result$exploration_summary$Biological.Summary[1L], "no literature-absence")
  f$raw <- NULL
  expect_error(explore_run(f), "pathways schema")
})

test_that("invalid enums, inline citations and sample-level cohort IDs are rejected", {
  for (case in list(c("effect", "supported"), c("scope", "direct"), c("evidence_types", "preclinical"),
    c("biological_processes", "BCR"), c("cohort_id", "TCGA-AA-1234-01"))) {
    f <- explore_fixture()
    f$raw$pathways[[1L]]$exploration$findings[[1L]][[case[1L]]] <- case[2L]
    expect_error(explore_run(f))
  }
  f <- explore_fixture()
  f$raw$pathways[[1L]]$exploration$biological_summary <- "See PMID 12345678 for details."
  expect_error(explore_run(f), "Inline citation")
  f$raw$pathways[[1L]]$exploration$biological_summary <- paste(rep("x", 421L), collapse = "")
  expect_error(explore_run(f), "exceeds 420")
})

test_that("exploration report round trips recompute source-grounded derived fields", {
  f <- explore_fixture()
  attr(f$raw, "intermeta_assessment_source") <- "assistant_abstract_audit"
  result <- explore_run(f)
  report <- InterMeta::InterMetaAIReport(f$context, f$evidence, result)
  again <- InterMeta::InterMetaAIReport(f$context, f$evidence, report$interpretation)
  expect_identical(again$exploration_summary, report$exploration_summary)
  expect_identical(again$exploration_findings, report$exploration_findings)
  expect_identical(report$exploration_summary$Assessment.Source, c("assistant_abstract_audit", "not_assessed"))
  result$exploration_findings$Functional <- FALSE
  result$exploration_findings$Quote.Match.Method <- "typography_equivalent"
  result$exploration_summary$Functional.PMIDs <- "99999999"
  fixed <- InterMeta::InterMetaAIReport(f$context, f$evidence, result)
  expect_true(fixed$exploration_findings$Functional)
  expect_identical(fixed$exploration_findings$Quote.Match.Method, "whitespace_exact")
  expect_identical(fixed$exploration_summary$Functional.PMIDs[1L], "12345678")

  f$raw$pathways[[1L]]$exploration$findings[[1L]]$abstract_quote <-
    "A model sentence absent from the supplied Abstract."
  rejected <- explore_run(f)
  rejected$exploration_findings$Quote.Verified <- TRUE
  rejected$exploration_findings$Quote.Match.Method <- "whitespace_exact"
  rejected$exploration_findings$Eligible.For.Summary <- TRUE
  rechecked <- InterMeta::InterMetaAIReport(f$context, f$evidence, rejected)
  expect_false(rechecked$exploration_findings$Quote.Verified)
  expect_false(rechecked$exploration_findings$Eligible.For.Summary)
  expect_identical(rechecked$exploration_findings$Quote.Match.Method,
    "source_mismatch")
})

test_that("purpose and provider-visible payload cannot change on revalidation", {
  f <- explore_fixture()
  result <- explore_run(f)
  expect_error(InterMeta::InterMetaAIReport(f$context, f$evidence, result,
    purpose = "prognostic_validation"), "purpose changed")
  result$purpose <- "prognostic_validation"
  expect_error(InterMeta::InterMetaAIReport(f$context, f$evidence, result), "selection purpose")
  result <- explore_run(f)
  f$evidence$records$Abstract[1L] <- "Changed source."
  expect_error(InterMeta::InterMetaAIReport(f$context, f$evidence, result), "Evidence records changed")
})

test_that("legacy interpretation without purpose keeps prognostic validation", {
  f <- explore_fixture()
  legacy <- list(pathways = list(list(pathway = "P1", literature_consistency = "insufficient_evidence",
    summary = "No clinical conclusion.", claims = list(), assessments = list())))
  report <- InterMeta::InterMetaAIReport(f$context, f$evidence, legacy)
  expect_identical(report$purpose, "prognostic_validation")
  expect_true(all(report$exploration_summary$Exploration.Status == "not_assessed"))
  expect_error(InterMeta::InterMetaAIReport(f$context, f$evidence, legacy,
    purpose = "biology_exploration"), "legacy interpretation")
  expect_identical(InterMeta::InterMetaAIReport(f$context)$purpose, "biology_exploration")
})

test_that("legacy fingerprints use the original selection and abstract allocation", {
  f <- explore_fixture()
  safe <- explore_fn(".InterMetaAISafeContext")(f$context)
  evidence <- explore_fn(".InterMetaAISafeEvidence")(f$evidence)
  selected <- explore_fn(".InterMetaAISelectRecords")(safe, evidence, 20L)
  prompt <- explore_fn(".InterMetaAIPromptEvidence")(evidence, selected)
  legacy <- list(pathways = list(list(pathway = "P1", literature_consistency = "insufficient_evidence",
    summary = "No clinical conclusion.", claims = list(), assessments = list())), record_selection = prompt$provenance)
  report <- InterMeta::InterMetaAIReport(f$context, f$evidence, legacy)
  expect_identical(report$interpretation$record_selection, prompt$provenance)
  expect_error(InterMeta::InterMetaAIReport(f$context, f$evidence, legacy,
    purpose = "biology_exploration"), "legacy interpretation")
})

test_that("waterfill keeps full short abstracts unless actual aggregate cap is exceeded", {
  f <- explore_fixture()
  f$evidence$records <- f$evidence$records[rep(1L, 100L), ]
  f$evidence$records$PMID <- as.character(20000000L + seq_len(100L))
  f$evidence$records$Abstract <- c(paste(rep("x", 3500L), collapse = ""), rep("short", 99L))
  selected <- explore_fn(".InterMetaAISelectRecords")(f$context, f$evidence, 100L, "biology_exploration")
  prompt <- explore_fn(".InterMetaAIPromptEvidence")(f$evidence, selected, "biology_exploration")
  expect_identical(prompt$provenance$effective_abstract_max_chars, 4000L)
  expect_identical(prompt$provenance$abstracts_truncated, 0L)
  f$evidence$records$Abstract <- rep(paste(rep("x", 4000L), collapse = ""), 100L)
  selected <- explore_fn(".InterMetaAISelectRecords")(f$context, f$evidence, 100L, "biology_exploration")
  prompt <- explore_fn(".InterMetaAIPromptEvidence")(f$evidence, selected, "biology_exploration")
  expect_identical(prompt$provenance$effective_abstract_max_chars, 2000L)
  expect_equal(prompt$provenance$abstract_characters_sent, 200000)
})

test_that("selection type hints diversify input without assigning evidence grades", {
  f <- explore_fixture()
  records <- f$evidence$records[rep(1L, 4L), ]
  records$PMID <- as.character(30000001:30000004)
  records$Title <- c("Cell lines knockdown", "Cell lines knockdown", "TCGA profiles", "Clinical cohort")
  records$Abstract <- "Study abstract."
  f$evidence$records <- records
  selected <- explore_fn(".InterMetaAISelectRecords")(f$context, f$evidence, 3L, "biology_exploration")
  expect_identical(selected$records$PMID, as.character(c(30000001, 30000003, 30000004)))
  prompt <- explore_fn(".InterMetaAIPromptEvidence")(f$evidence, selected, "biology_exploration")
  expect_match(prompt$provenance$selection_note, "not validated study classifications")
  expect_identical(prompt$provenance$selected_type_hints$Type.Hint,
    c("functional_hint", "omics_hint", "patient_hint"))
})

test_that("exploration artifacts have a complete sixteen-file manifest and no BCR grading section", {
  f <- explore_fixture()
  target <- tempfile("exploration-report-")
  on.exit(unlink(target, recursive = TRUE), add = TRUE)
  report <- InterMeta::InterMetaAIReport(f$context, f$evidence, explore_run(f), outdir = target)
  expect_length(list.files(target), 16L)
  expect_true(all(c("exploration_summary.tsv", "exploration_findings.tsv") %in% list.files(target)))
  text <- paste(readLines(file.path(target, "report.md")), collapse = "\n")
  expect_match(text, "Nonexclusive evidence channels")
  expect_false(grepl("Endpoint-specific directional assessment", text, fixed = TRUE))
  manifest <- utils::read.delim(file.path(target, "output_manifest.txt"))
  expect_true(all(file.exists(file.path(target, manifest$File))))
  expect_identical(report$provenance$purpose, "biology_exploration")
})

test_that("exploration reports label prognostic compatibility tables without fake evidence rows", {
  f <- explore_fixture()
  target <- tempfile("exploration-semantics-")
  on.exit(unlink(target, recursive = TRUE), add = TRUE)
  report <- InterMeta::InterMetaAIReport(
    f$context, f$evidence, explore_run(f), outdir = target
  )

  semantics <- report$table_semantics
  legacy <- semantics$Table %in% c(
    "claims.tsv", "assessments.tsv", "pathway_synthesis.tsv",
    "biology_findings.tsv"
  )
  expect_true(all(semantics$Role[legacy] == "legacy_prognostic_compatibility"))
  expect_true(all(semantics$Applicability[legacy] == "not_applicable"))
  expect_false(any(semantics$Authoritative[legacy]))
  expect_identical(
    semantics$Table[semantics$Authoritative],
    c("exploration_summary.tsv", "exploration_findings.tsv")
  )
  expect_identical(report$provenance$table_semantics, semantics)
  expect_identical(report$interpretation$table_semantics, semantics)

  # Empty legacy tables remain empty: explicit semantics must not be encoded as
  # fabricated publications, claims, or PMIDs.
  expect_identical(nrow(report$interpretation$claims), 0L)
  expect_identical(nrow(report$interpretation$assessments), 0L)
  expect_identical(nrow(report$biology_findings), 0L)
  expect_identical(
    names(report$interpretation$claims),
    names(explore_fn(".InterMetaAIEmptyClaims")())
  )
  expect_identical(
    names(report$interpretation$assessments),
    names(explore_fn(".InterMetaAIEmptyAssessments")(include_derived = TRUE))
  )
  expect_identical(
    names(report$biology_findings),
    c(
      names(explore_fn(".InterMetaAIEmptyBiologyFindings")()),
      "Quote.Verified", "Disease.Match.Verified", "Evidence.Type",
      "Null.Finding", "Eligible.For.Mechanism", "Assessment.Source"
    )
  )

  synthesis <- report$pathway_synthesis
  expect_true(all(synthesis$Synthesis.Status == "not_applicable"))
  expect_true(all(synthesis$Biological.Explanation == "not_applicable"))
  expect_true(all(synthesis$Literature.Consistency == "not_applicable"))
  expect_true(all(synthesis$Assessment.Source == "legacy_prognostic_compatibility"))
  expect_false(any(grepl("not_assessed|insufficient_evidence", unlist(synthesis))))

  published <- utils::read.delim(
    file.path(target, "pathway_synthesis.tsv"), check.names = FALSE
  )
  expect_true(all(published$Synthesis.Status == "not_applicable"))
  expect_true(all(published$Assessment.Source == "legacy_prognostic_compatibility"))
  expect_identical(nrow(utils::read.delim(file.path(target, "claims.tsv"))), 0L)
  expect_identical(nrow(utils::read.delim(file.path(target, "assessments.tsv"))), 0L)
  expect_identical(nrow(utils::read.delim(file.path(target, "biology_findings.tsv"))), 0L)

  text <- paste(readLines(file.path(target, "report.md")), collapse = "\n")
  expect_match(text, "legacy prognostic-output compatibility", fixed = TRUE)
  expect_match(text, "empty compatibility tables", fixed = TRUE)
  expect_match(text, "authoritative exploration outputs", fixed = TRUE)
})

test_that("prognostic validation keeps its established compatibility-table semantics", {
  f <- explore_fixture()
  legacy <- list(pathways = list(list(
    pathway = "P1", literature_consistency = "insufficient_evidence",
    summary = "No clinical conclusion.", claims = list(), assessments = list()
  )))
  report <- InterMeta::InterMetaAIReport(f$context, f$evidence, legacy)
  expect_identical(report$purpose, "prognostic_validation")
  expect_null(report$table_semantics)
  expect_null(report$interpretation$table_semantics)
  expect_true(all(report$pathway_synthesis$Synthesis.Status == "not_assessed"))
  expect_true(all(report$pathway_synthesis$Literature.Consistency == "insufficient_evidence"))
})
