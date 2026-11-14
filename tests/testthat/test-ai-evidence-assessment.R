# Existing clinical-mode regression fixtures explicitly retain their original goal.
clinical_ai_interpret <- function(..., max_records_per_pathway = 20L, purpose = "prognostic_validation")
  InterMeta::InterMetaAIInterpret(..., max_records_per_pathway = max_records_per_pathway, purpose = purpose)
clinical_ai_report <- function(..., max_records_per_pathway = 20L, purpose = "prognostic_validation")
  InterMeta::InterMetaAIReport(..., max_records_per_pathway = max_records_per_pathway, purpose = purpose)

assessment_fixture <- function(
    pmid = "45678901", evidence_type = "human_prognosis",
    study_endpoint = "biochemical recurrence", endpoint_match = "matched",
    direction = "same_direction", directness = "direct",
    quote = "Higher pathway score was associated with shorter biochemical recurrence-free survival.",
    grade = "supported", endpoint = "BCR (biochemical recurrence)") {
  context <- structure(list(
    analysis = list(disease = "prostate cancer", endpoint = endpoint,
                    time_unit = "months"),
    pathways = data.frame(Pathway = "Lysine degradation", Direction = "Unfavorable",
                          Direction.Rank = 1L, stringsAsFactors = FALSE),
    interactions = data.frame(), driver_gene_pairs = data.frame(),
    warnings = character(), provenance = list()
  ), class = c("InterMetaAIContext", "list"))
  evidence <- getFromNamespace(".InterMetaAIEmptyEvidence", "InterMeta")()
  evidence$records <- data.frame(
    pathway = "Lysine degradation", PMID = pmid, layer = "prognosis",
    Title = "An illustrative pathway study", Abstract = quote,
    stringsAsFactors = FALSE
  )
  entry <- list(
    pathway = "Lysine degradation", literature_consistency = grade,
    summary = "The model's raw directional interpretation.",
    claims = list(list(claim = "The study discusses this metabolic pathway.",
                       pmids = pmid)),
    assessments = list(list(
      pmid = pmid, evidence_type = evidence_type, study_endpoint = study_endpoint,
      endpoint_match = endpoint_match, direction = direction,
      directness = directness, abstract_quote = quote
    ))
  )
  list(context = context, evidence = evidence, raw = list(pathways = list(entry)))
}

interpret_assessment_fixture <- function(fixture) {
  clinical_ai_interpret(
    fixture$context, fixture$evidence,
    InterMetaAIProvider(function(messages, schema) fixture$raw)
  )
}

test_that("schema 3 preserves distinct clinical assessments for every cited PMID", {
  schema <- getFromNamespace(".InterMetaAIInterpretationSchema", "InterMeta")()
  item <- schema$properties$pathways$items
  expect_true("assessments" %in% item$required)
  expect_setequal(item$properties$assessments$items$required, c(
    "pmid", "evidence_type", "study_endpoint", "endpoint_match", "direction",
    "directness", "abstract_quote"
  ))
  result <- interpret_assessment_fixture(assessment_fixture())
  expect_identical(result$schema_version, "4.3.0")
})

test_that("PMID 25683757 null incidence evidence cannot contradict BCR prognosis", {
  fixture <- assessment_fixture(
    pmid = "25683757", evidence_type = "human_incidence",
    study_endpoint = "prostate cancer incidence", endpoint_match = "mismatched",
    direction = "no_association", grade = "potentially_contradicted",
    quote = "Lysine concentrations were not associated with prostate cancer risk."
  )
  result <- interpret_assessment_fixture(fixture)
  expect_identical(result$pathways$Literature.Consistency, "insufficient_evidence")
  expect_identical(result$pathways$Biological.Relevance, "contextual_evidence")
  expect_identical(result$pathways$Model.Literature.Consistency, "potentially_contradicted")
  expect_identical(result$pathways$Model.Summary, fixture$raw$pathways[[1L]]$summary)
  expect_false(result$assessments$Eligible.For.Consistency)
  expect_match(result$pathways$Summary, "insufficient to support or contradict")
})

test_that("preclinical mTOR and tissue evidence give context, not directional support", {
  for (type in c("preclinical", "human_tissue", "review", "other", "unknown")) {
    fixture <- assessment_fixture(
      evidence_type = type, directness = "indirect", study_endpoint = "mTOR activity",
      endpoint_match = "unknown", quote = "Lysine affected mTOR activity in cultured cells."
    )
    result <- interpret_assessment_fixture(fixture)
    expect_identical(result$pathways$Literature.Consistency, "insufficient_evidence")
    expect_identical(result$pathways$Biological.Relevance, "contextual_evidence")
  }
  fixture <- assessment_fixture(directness = "indirect")
  expect_identical(interpret_assessment_fixture(fixture)$pathways$Literature.Consistency,
                   "insufficient_evidence")
})

test_that("only direct matched human BCR directional findings contribute to grades", {
  result <- interpret_assessment_fixture(assessment_fixture())
  expect_identical(result$pathways$Literature.Consistency, "supported")
  expect_true(result$assessments$Quote.Verified)
  expect_true(result$assessments$Eligible.For.Consistency)
  expect_identical(interpret_assessment_fixture(assessment_fixture(
    direction = "opposite_direction"))$pathways$Literature.Consistency,
    "potentially_contradicted")
  null <- interpret_assessment_fixture(assessment_fixture(direction = "no_association"))
  expect_true(null$assessments$Eligible.For.Consistency)
  expect_identical(null$pathways$Literature.Consistency, "insufficient_evidence")
  not_assessed <- interpret_assessment_fixture(assessment_fixture(direction = "not_assessed"))
  expect_identical(not_assessed$pathways$Literature.Consistency, "insufficient_evidence")
})

test_that("OS and risk labels cannot claim an endpoint match to BCR", {
  for (study_endpoint in c("OS", "overall survival", "prostate cancer risk", "unknown")) {
    fixture <- assessment_fixture(study_endpoint = study_endpoint, endpoint_match = "matched")
    result <- interpret_assessment_fixture(fixture)
    expect_false(result$assessments$Eligible.For.Consistency)
    expect_identical(result$pathways$Literature.Consistency, "insufficient_evidence")
  }
  for (study_endpoint in c("BCR", "biochemical relapse", "biochemical recurrence")) {
    result <- interpret_assessment_fixture(assessment_fixture(study_endpoint = study_endpoint))
    expect_identical(result$pathways$Literature.Consistency, "supported")
  }
  result <- interpret_assessment_fixture(assessment_fixture(endpoint_match = "unknown"))
  expect_identical(result$pathways$Literature.Consistency, "insufficient_evidence")
})

test_that("missing target endpoints block all endpoint-specific grades", {
  for (endpoint in list(NULL, NA_character_, "", "unknown", "not specified")) {
    fixture <- assessment_fixture(endpoint = endpoint)
    result <- interpret_assessment_fixture(fixture)
    expect_false(result$assessments$Eligible.For.Consistency)
    expect_identical(result$pathways$Literature.Consistency, "insufficient_evidence")
  }
})

test_that("biochemical recurrence-free survival is BCR without a spurious RFS endpoint", {
  endpoint_key <- getFromNamespace(".InterMetaAIEndpointKey", "InterMeta")
  for (endpoint in c("biochemical recurrence-free survival",
                     "biochemical recurrence free survival",
                     "biochemical relapse-free survival")) {
    expect_identical(endpoint_key(endpoint), "bcr")
    result <- interpret_assessment_fixture(assessment_fixture(study_endpoint = endpoint))
    expect_identical(result$pathways$Literature.Consistency, "supported")
    reverse <- interpret_assessment_fixture(assessment_fixture(
      study_endpoint = "BCR", endpoint = endpoint))
    expect_identical(reverse$pathways$Literature.Consistency, "supported")
  }
  for (endpoint in c("biochemical recurrence-free survival and OS",
                     "BCR and overall survival", "BCR and RFS",
                     "biochemical recurrence-free survival and recurrence-free survival")) {
    expect_match(endpoint_key(endpoint), "|", fixed = TRUE)
    result <- interpret_assessment_fixture(assessment_fixture(study_endpoint = endpoint))
    expect_identical(result$pathways$Literature.Consistency, "insufficient_evidence")
    expect_false(result$assessments$Eligible.For.Consistency)
  }
})

test_that("legacy provider outputs are unassessed and audit-preserving on revalidation", {
  fixture <- assessment_fixture()
  fixture$raw$pathways[[1L]]$assessments <- NULL
  result <- interpret_assessment_fixture(fixture)
  expect_identical(result$assessments$Assessment.Status, "unassessed")
  expect_identical(result$pathways$Literature.Consistency, "insufficient_evidence")
  expect_true(any(grepl("Legacy citations", result$warnings)))
  report <- clinical_ai_report(fixture$context, fixture$evidence, interpretation = result)
  expect_identical(report$interpretation$pathways, result$pathways)
  expect_identical(report$interpretation$assessments, result$assessments)
  # Old normalized objects may be passed without their S3 class.
  normalized <- unclass(result)
  normalized$assessments <- NULL
  report2 <- clinical_ai_report(fixture$context, fixture$evidence, interpretation = normalized)
  expect_identical(report2$interpretation$pathways, result$pathways)
  expect_identical(report2$interpretation$pathways$Model.Literature.Consistency, "supported")
})

test_that("normalized structured interpretation is idempotent", {
  fixture <- assessment_fixture()
  result <- interpret_assessment_fixture(fixture)
  report <- clinical_ai_report(fixture$context, fixture$evidence, interpretation = result)
  expect_identical(report$interpretation$pathways, result$pathways)
  expect_identical(report$interpretation$assessments, result$assessments)
  report2 <- clinical_ai_report(fixture$context, fixture$evidence,
                              interpretation = report$interpretation)
  expect_identical(report2$interpretation$pathways, result$pathways)
  expect_identical(report2$interpretation$assessments, result$assessments)
})

test_that("quotes must be exact abstract excerpts, with whitespace-only normalization", {
  fixture <- assessment_fixture(quote = "Biochemical recurrence was associated with the pathway.")
  fixture$evidence$records$Abstract <- "Results: Biochemical\n recurrence was\tassociated with the pathway. Conclusion."
  expect_true(interpret_assessment_fixture(fixture)$assessments$Quote.Verified)
  fixture$raw$pathways[[1L]]$assessments[[1L]]$abstract_quote <- "Invented abstract result."
  fixture$evidence$records$Title <- "Invented abstract result."
  expect_error(interpret_assessment_fixture(fixture), "not an exact excerpt.*Abstract")
  fixture$raw$pathways[[1L]]$assessments[[1L]]$abstract_quote <-
    "biochemical recurrence was associated with the pathway."
  expect_error(interpret_assessment_fixture(fixture), "not an exact excerpt")
  fixture$raw$pathways[[1L]]$assessments[[1L]]$abstract_quote <- ""
  result <- interpret_assessment_fixture(fixture)
  expect_false(result$assessments$Quote.Verified)
  expect_identical(result$pathways$Literature.Consistency, "insufficient_evidence")
})

test_that("assessments cannot cite another pathway, uncited PMID, or duplicate a source", {
  fixture <- assessment_fixture()
  fixture$raw$pathways[[1L]]$assessments[[1L]]$pmid <- "99999999"
  expect_error(interpret_assessment_fixture(fixture), "must be cited for the same pathway")
  fixture <- assessment_fixture()
  fixture$raw$pathways[[1L]]$assessments <- rep(fixture$raw$pathways[[1L]]$assessments, 2L)
  expect_error(interpret_assessment_fixture(fixture), "Duplicate pathway/PMID")
  fixture <- assessment_fixture()
  fixture$evidence$records$pathway <- "Another pathway"
  expect_error(interpret_assessment_fixture(fixture), "pathway-specific evidence")
})

test_that("invalid enums, empty endpoint and malformed assessment identifiers fail", {
  for (field in c("evidence_type", "endpoint_match", "direction", "directness")) {
    fixture <- assessment_fixture()
    fixture$raw$pathways[[1L]]$assessments[[1L]][[field]] <- "invented"
    expect_error(interpret_assessment_fixture(fixture), paste("Invalid assessment", field))
  }
  fixture <- assessment_fixture(study_endpoint = " ")
  expect_error(interpret_assessment_fixture(fixture), "non-empty character")
  fixture <- assessment_fixture()
  fixture$raw$pathways[[1L]]$assessments[[1L]]$pmid <- "45678901,45678902"
  expect_error(interpret_assessment_fixture(fixture), "exactly one numeric PMID")
})

test_that("multiple cited PMIDs require individual assessments and opposing findings yield mixed", {
  fixture <- assessment_fixture()
  fixture$evidence$records <- fixture$evidence$records[c(1L, 1L), , drop = FALSE]
  fixture$evidence$records$PMID[2L] <- "45678902"
  fixture$raw$pathways[[1L]]$claims[[1L]]$pmids <- c("45678901", "45678902")
  expect_error(interpret_assessment_fixture(fixture), "cover every cited PMID")
  second <- fixture$raw$pathways[[1L]]$assessments[[1L]]
  second$pmid <- "45678902"
  second$direction <- "opposite_direction"
  fixture$raw$pathways[[1L]]$assessments[[2L]] <- second
  result <- interpret_assessment_fixture(fixture)
  expect_equal(nrow(result$assessments), 2L)
  expect_identical(result$pathways$Literature.Consistency, "mixed")
  # Duplicated retrieval rows are not repeated votes or independent assessments.
  fixture$evidence$records <- fixture$evidence$records[c(1L, 1L, 2L), , drop = FALSE]
  expect_equal(nrow(interpret_assessment_fixture(fixture)$assessments), 2L)
})

test_that("quotes outside truncated provider abstracts are rejected", {
  fixture <- assessment_fixture(quote = "The relevant endpoint result appears only at the end.")
  fixture$evidence$records$Abstract <- paste0(paste(rep("x", 4100L), collapse = ""),
                                            fixture$evidence$records$Abstract)
  expect_error(interpret_assessment_fixture(fixture), "not an exact excerpt")
})

test_that("component retrieval metadata is available to providers and report previews", {
  fixture <- assessment_fixture()
  fixture$evidence$records$layer <- "component_prognosis"
  fixture$evidence$records$term_scope <- "component"
  fixture$evidence$records$search_terms <- "lysine"
  fixture$evidence$records$query_version <- "2.0.0"
  fixture$evidence$records$query_warning <- "Component evidence is indirect."
  seen <- new.env(parent = emptyenv())
  provider <- InterMetaAIProvider(function(messages, schema) {
    seen$system <- messages[[1L]]$content
    seen$records <- messages[[2L]]$content$evidence$records
    fixture$raw
  })
  result <- clinical_ai_interpret(fixture$context, fixture$evidence, provider)
  expect_identical(seen$records$term_scope, "component")
  expect_identical(seen$records$search_terms, "lysine")
  expect_identical(seen$records$query_version, "2.0.0")
  expect_match(seen$system, "No association")
  expect_match(seen$system, "Component or driver-gene retrieval")
  report <- clinical_ai_report(fixture$context, fixture$evidence, interpretation = result)
  markdown <- paste(getFromNamespace(".InterMetaAIRenderMarkdown", "InterMeta")(report), collapse = "\n")
  expect_match(markdown, "Endpoint-specific directional assessment", fixed = TRUE)
  expect_match(markdown, "Biological context", fixed = TRUE)
  expect_match(markdown, "human review", fixed = TRUE)
  expect_match(markdown, "term_scope", fixed = TRUE)
})

test_that("report revalidation cannot borrow records excluded from the provider payload", {
  fixture <- assessment_fixture()
  fixture$evidence$records <- fixture$evidence$records[c(1L, 1L), , drop = FALSE]
  fixture$evidence$records$PMID[2L] <- "45678902"
  result <- clinical_ai_interpret(
    fixture$context, fixture$evidence,
    InterMetaAIProvider(function(messages, schema) fixture$raw),
    max_records_per_pathway = 1L
  )
  result$claims$PMID <- "45678902"
  result$assessments$PMID <- "45678902"
  expect_error(clinical_ai_report(fixture$context, fixture$evidence,
                                interpretation = result), "pathway-specific evidence")
})

test_that("report revalidation retains the actual abstract length boundary", {
  fixture <- assessment_fixture(quote = "BCR was associated with the pathway.")
  fixture$evidence$records$Abstract <- paste0(fixture$evidence$records$Abstract,
    paste(rep("x", 4100L), collapse = ""), "An excluded final result.")
  result <- interpret_assessment_fixture(fixture)
  result$assessments$Abstract.Quote <- "An excluded final result."
  expect_error(clinical_ai_report(fixture$context, fixture$evidence,
                                interpretation = result), "not an exact excerpt")
})

test_that("selection boundaries survive classless report round trips", {
  fixture <- assessment_fixture()
  fixture$evidence$records <- fixture$evidence$records[c(1L, 1L), , drop = FALSE]
  fixture$evidence$records$PMID[2L] <- "45678902"
  result <- clinical_ai_interpret(
    fixture$context, fixture$evidence,
    InterMetaAIProvider(function(messages, schema) fixture$raw),
    max_records_per_pathway = 1L
  )
  first <- clinical_ai_report(fixture$context, fixture$evidence, interpretation = result)
  second <- clinical_ai_report(fixture$context, fixture$evidence,
                              interpretation = first$interpretation)
  expect_identical(first$interpretation$record_selection, result$record_selection)
  expect_identical(second$interpretation$record_selection, result$record_selection)
  expect_identical(second$provenance$interpretation_record_selection,
                   result$record_selection)
  expect_match(second$interpretation$record_selection$source_records_sha256,
                "^[a-f0-9]{64}$")
  second$interpretation$claims$PMID <- "45678902"
  second$interpretation$assessments$PMID <- "45678902"
  expect_error(clinical_ai_report(fixture$context, fixture$evidence,
    interpretation = second$interpretation), "pathway-specific evidence")
})

test_that("source fingerprints reject reordered records before payload reselection", {
  fixture <- assessment_fixture()
  fixture$evidence$records <- fixture$evidence$records[c(1L, 1L), , drop = FALSE]
  fixture$evidence$records$PMID[2L] <- "45678902"
  result <- clinical_ai_interpret(
    fixture$context, fixture$evidence,
    InterMetaAIProvider(function(messages, schema) fixture$raw),
    max_records_per_pathway = 1L
  )
  original <- clinical_ai_report(fixture$context, fixture$evidence,
                                interpretation = result)
  expect_identical(original$interpretation$pathways$Literature.Consistency, "supported")
  reordered <- fixture$evidence
  reordered$records <- reordered$records[2:1, , drop = FALSE]
  expect_error(clinical_ai_report(fixture$context, reordered, interpretation = result),
                "Evidence records changed since the provider call")
  changed <- result
  changed$claims$PMID <- "45678902"
  changed$assessments$PMID <- "45678902"
  expect_error(clinical_ai_report(fixture$context, reordered, interpretation = changed),
                "Evidence records changed since the provider call")
  classless <- original$interpretation
  expect_error(clinical_ai_report(fixture$context, reordered, interpretation = classless),
                "Evidence records changed since the provider call")
})

test_that("source fingerprints detect changed and truncated abstracts", {
  fixture <- assessment_fixture()
  result <- interpret_assessment_fixture(fixture)
  for (abstract in c(paste(fixture$evidence$records$Abstract, "New final result."),
                     substr(fixture$evidence$records$Abstract, 1L, 20L))) {
    changed <- fixture$evidence
    changed$records$Abstract <- abstract
    expect_error(clinical_ai_report(fixture$context, changed, interpretation = result),
                  "Evidence records changed since the provider call")
  }
})

test_that("fingerprints use safe evidence fields only and remain stable on safe round trips", {
  fixture <- assessment_fixture()
  fixture$evidence$records$Patient.ID <- "PRIVATE-SOURCE-IDENTIFIER"
  fixture$evidence$records$api_key <- "PRIVATE-SOURCE-KEY"
  result <- interpret_assessment_fixture(fixture)
  fingerprint <- result$record_selection$source_records_sha256
  expect_match(fingerprint, "^[a-f0-9]{64}$")
  fixture$evidence$records$Patient.ID <- "A-DIFFERENT-PRIVATE-IDENTIFIER"
  fixture$evidence$records$api_key <- NULL
  report <- clinical_ai_report(fixture$context, fixture$evidence, interpretation = result)
  expect_identical(report$interpretation$record_selection$source_records_sha256, fingerprint)
  expect_false(grepl("PRIVATE", paste(capture.output(dput(report)), collapse = "")))
  report2 <- clinical_ai_report(fixture$context, fixture$evidence,
                               interpretation = report$interpretation)
  expect_identical(report2$provenance$interpretation_record_selection$source_records_sha256,
                   fingerprint)
})

test_that("legacy selection metadata without a source fingerprint gives an explicit warning", {
  fixture <- assessment_fixture()
  fixture$raw$pathways[[1L]]$assessments <- NULL
  result <- interpret_assessment_fixture(fixture)
  result$schema_version <- "1.0.0"
  result$record_selection$source_records_sha256 <- NULL
  report <- clinical_ai_report(fixture$context, fixture$evidence, interpretation = result)
  expect_true(any(grepl("no source-record fingerprint", report$warnings, fixed = TRUE)))
  expect_true(any(grepl("original provider evidence cannot be verified",
                         report$warnings, fixed = TRUE)))
  expect_identical(report$interpretation$pathways$Literature.Consistency,
                   "insufficient_evidence")
  expect_null(report$interpretation$record_selection$source_records_sha256)
})

test_that("Chinese reports have Chinese safe summaries without replacing raw model audit", {
  fixture <- assessment_fixture(evidence_type = "preclinical", directness = "indirect")
  report <- clinical_ai_report(fixture$context, fixture$evidence,
    provider = InterMetaAIProvider(function(messages, schema) fixture$raw), language = "zh")
  expect_match(report$interpretation$pathways$Summary, "\u8bc1\u636e\u4e0d\u8db3")
  expect_identical(report$interpretation$pathways$Model.Summary,
                   fixture$raw$pathways[[1L]]$summary)
  expect_identical(report$interpretation$pathways$Literature.Consistency,
                   "insufficient_evidence")
})
