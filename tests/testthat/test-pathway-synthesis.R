# Existing clinical-mode regression fixtures explicitly retain their original goal.
clinical_ai_interpret <- function(..., max_records_per_pathway = 20L, purpose = "prognostic_validation")
  InterMeta::InterMetaAIInterpret(..., max_records_per_pathway = max_records_per_pathway, purpose = purpose)
clinical_ai_report <- function(..., max_records_per_pathway = 20L, purpose = "prognostic_validation")
  InterMeta::InterMetaAIReport(..., max_records_per_pathway = max_records_per_pathway, purpose = purpose)

synthesis_fn <- function(name) getFromNamespace(name, "InterMeta")

synthesis_fixture <- function() {
  context <- structure(list(
    analysis = list(disease = "prostate adenocarcinoma", endpoint = "BCR", time_unit = "months"),
    pathways = data.frame(Pathway = c("P1", "P2"), Direction = "Unfavorable",
      Direction.Rank = 1:2, Difference = c(-7.25, -3.5), stringsAsFactors = FALSE),
    interactions = data.frame(), driver_gene_pairs = data.frame(), warnings = character(),
    provenance = list()), class = c("InterMetaAIContext", "list"))
  evidence <- synthesis_fn(".InterMetaAIEmptyEvidence")()
  quote <- "Experimental pathway inhibition reduced prostate tumor growth."
  evidence$records <- data.frame(pathway = "P1", PMID = "12345678", layer = "biology",
    Title = "An experimental study", Abstract = quote, stringsAsFactors = FALSE)
  finding <- list(pmid = "12345678", study_disease = "prostate cancer",
    disease_match = "matched", scope = "pathway", finding_type = "mechanistic",
    study_design = "experimental", relation = "supports",
    claim = "Pathway perturbation changed tumor growth in the studied prostate model.",
    abstract_quote = quote)
  entry <- list(pathway = "P1", literature_consistency = "insufficient_evidence",
    summary = "Raw endpoint interpretation.",
    claims = list(list(claim = "The study examines pathway perturbation.", pmids = "12345678")),
    assessments = list(list(pmid = "12345678", evidence_type = "preclinical",
      study_endpoint = "tumor growth", endpoint_match = "mismatched", direction = "not_assessed",
      directness = "indirect", abstract_quote = quote)),
    biology_synthesis = list(mechanism_summary = "A proposed cross-study biological mechanism.",
      limitations = "Experiments used a specific model and do not establish clinical utility.",
      findings = list(finding)))
  list(context = context, evidence = evidence, raw = list(pathways = list(entry)))
}

synthesis_interpret <- function(fixture) {
  clinical_ai_interpret(fixture$context, fixture$evidence,
    InterMetaAIProvider(function(messages, schema) fixture$raw))
}

synthesis_second <- function(fixture, relation = "opposes", scope = "pathway") {
  row <- fixture$evidence$records
  row$PMID <- "12345679"
  row$Abstract <- "A comparable experiment found a different pathway-mediated growth response."
  fixture$evidence$records <- rbind(fixture$evidence$records, row)
  fixture$raw$pathways[[1L]]$claims[[1L]]$pmids <- c("12345678", "12345679")
  second <- fixture$raw$pathways[[1L]]$assessments[[1L]]
  second$pmid <- "12345679"
  second$abstract_quote <- row$Abstract
  fixture$raw$pathways[[1L]]$assessments[[2L]] <- second
  finding <- fixture$raw$pathways[[1L]]$biology_synthesis$findings[[1L]]
  finding$pmid <- "12345679"
  finding$scope <- scope
  finding$relation <- relation
  finding$claim <- "A different perturbation context showed an opposing mechanism-related growth response."
  finding$abstract_quote <- row$Abstract
  fixture$raw$pathways[[1L]]$biology_synthesis$findings[[2L]] <- finding
  fixture
}

test_that("schema 3 adds complete synthesis to the same provider response", {
  schema <- synthesis_fn(".InterMetaAIInterpretationSchema")()
  path <- schema$properties$pathways$items
  expect_true("biology_synthesis" %in% path$required)
  expect_setequal(path$properties$biology_synthesis$required,
                   c("mechanism_summary", "limitations", "findings"))
  expect_true("study_design" %in% path$properties$biology_synthesis$properties$findings$items$required)
  fixture <- synthesis_fixture()
  calls <- 0L
  seen <- new.env(parent = emptyenv())
  result <- clinical_ai_interpret(fixture$context, fixture$evidence,
    InterMetaAIProvider(function(messages, schema) {
      calls <<- calls + 1L
      seen$prompt <- messages[[1L]]$content
      fixture$raw
    }))
  expect_identical(calls, 1L)
  expect_identical(result$schema_version, "4.3.0")
  expect_match(seen$prompt, "not pathway activation, flux", fixed = TRUE)
  expect_match(seen$prompt, "25-35 word claim", fixed = TRUE)
})

test_that("pathway mechanism evidence is distinct from BCR direction and signed score", {
  fixture <- synthesis_fixture()
  result <- synthesis_interpret(fixture)
  table <- result$pathway_synthesis
  expect_identical(table$Pathway, fixture$context$pathways$Pathway)
  expect_identical(table$Difference, c(-7.25, -3.5))
  expect_identical(table$Biological.Explanation, c("pathway_mechanism", "not_assessed"))
  expect_true(all(table$Literature.Consistency == "insufficient_evidence"))
  expect_true(result$biology_findings$Eligible.For.Mechanism)
  expect_true(result$biology_findings$Quote.Verified)
  expect_true(result$biology_findings$Disease.Match.Verified)
  expect_identical(table$Mechanism.Summary[1L],
    fixture$raw$pathways[[1L]]$biology_synthesis$findings[[1L]]$claim)
  expect_identical(table$Model.Mechanism.Summary[1L],
    fixture$raw$pathways[[1L]]$biology_synthesis$mechanism_summary)
  expect_identical(table$Key.PMIDs, c("12345678", ""))
})

test_that("component signatures, associations, and reviews cannot become pathway mechanism", {
  variants <- list(scope = "component", finding_type = "associational", study_design = "observational")
  for (field in names(variants)) {
    fixture <- synthesis_fixture()
    fixture$raw$pathways[[1L]]$biology_synthesis$findings[[1L]][[field]] <- variants[[field]]
    result <- synthesis_interpret(fixture)
    expect_identical(result$pathway_synthesis$Biological.Explanation[1L], "indirect_context")
    expect_false(result$biology_findings$Eligible.For.Mechanism)
  }
  fixture <- synthesis_fixture()
  fixture$raw$pathways[[1L]]$assessments[[1L]]$evidence_type <- "review"
  result <- synthesis_interpret(fixture)
  expect_identical(result$pathway_synthesis$Biological.Explanation[1L], "indirect_context")
  expect_false(result$biology_findings$Eligible.For.Mechanism)
})

test_that("only eligible comparable mechanism findings can yield the mixed biological label", {
  fixture <- synthesis_second(synthesis_fixture())
  result <- synthesis_interpret(fixture)
  expect_identical(result$pathway_synthesis$Biological.Explanation[1L], "conflicting_evidence")
  expect_match(result$pathway_synthesis$Mechanism.Summary[1L], "Contrasting finding")
  expect_identical(result$pathway_synthesis$Key.PMIDs[1L], "12345678;12345679")
  expect_match(result$pathway_synthesis$Limitations[1L], "Context-dependent or mixed")
  expect_identical(result$pathways$Literature.Consistency, "insufficient_evidence")
  fixture <- synthesis_second(synthesis_fixture(), scope = "component")
  expect_identical(synthesis_interpret(fixture)$pathway_synthesis$Biological.Explanation[1L], "pathway_mechanism")
  fixture <- synthesis_second(synthesis_fixture())
  fixture$raw$pathways[[1L]]$biology_synthesis$findings[[2L]]$finding_type <- "associational"
  expect_identical(synthesis_interpret(fixture)$pathway_synthesis$Biological.Explanation[1L], "pathway_mechanism")
})

test_that("null findings and other-disease findings never fabricate biological opposition", {
  fixture <- synthesis_second(synthesis_fixture(), relation = "neutral")
  fixture$raw$pathways[[1L]]$assessments[[2L]]$direction <- "no_association"
  result <- synthesis_interpret(fixture)
  expect_identical(result$pathway_synthesis$Biological.Explanation[1L], "pathway_mechanism")
  expect_true(result$biology_findings$Null.Finding[2L])
  expect_false(result$biology_findings$Eligible.For.Mechanism[2L])
  expect_match(result$pathway_synthesis$Mechanism.Summary[1L], "Neutral finding")
  fixture <- synthesis_second(synthesis_fixture())
  fixture$raw$pathways[[1L]]$biology_synthesis$findings[[2L]]$study_disease <- "breast cancer"
  expect_identical(synthesis_interpret(fixture)$pathway_synthesis$Biological.Explanation[1L], "pathway_mechanism")
})

test_that("disease claims must be explicit and lexically compatible with the target", {
  for (disease in c("breast cancer", "unknown", "cancer", "benign prostate hyperplasia",
                     "prostate cells", "prostate cancer and breast cancer",
                     "prostate and bladder cancer", "urinary cancers including prostate and bladder cancer",
                     "prostate cancer or melanoma", "prostate cancer; glioblastoma")) {
    fixture <- synthesis_fixture()
    fixture$raw$pathways[[1L]]$biology_synthesis$findings[[1L]]$study_disease <- disease
    result <- synthesis_interpret(fixture)
    expect_identical(result$pathway_synthesis$Biological.Explanation[1L], "insufficient_evidence")
    expect_false(result$biology_findings$Disease.Match.Verified)
    expect_false(grepl("Pathway perturbation changed", result$pathway_synthesis$Mechanism.Summary[1L]))
  }
  fixture <- synthesis_fixture()
  fixture$raw$pathways[[1L]]$biology_synthesis$findings[[1L]]$disease_match <- "unknown"
  expect_identical(synthesis_interpret(fixture)$pathway_synthesis$Biological.Explanation[1L], "insufficient_evidence")
})

test_that("null clinical association cannot suppress independent experimental biology", {
  fixture <- synthesis_fixture()
  fixture$raw$pathways[[1L]]$assessments[[1L]]$direction <- "no_association"
  result <- synthesis_interpret(fixture)
  expect_identical(result$pathway_synthesis$Biological.Explanation[1L], "pathway_mechanism")
  expect_false(result$biology_findings$Null.Finding[1L])
  expect_true(result$biology_findings$Eligible.For.Mechanism[1L])
  expect_identical(result$pathway_synthesis$Literature.Consistency[1L], "insufficient_evidence")
})

test_that("revalidation cannot expand or falsify the actual provider selection", {
  fixture <- synthesis_second(synthesis_fixture(), relation = "supports")
  one <- fixture
  one$raw$pathways[[1L]]$claims[[1L]]$pmids <- "12345678"
  one$raw$pathways[[1L]]$assessments <- one$raw$pathways[[1L]]$assessments[1L]
  one$raw$pathways[[1L]]$biology_synthesis$findings <-
    one$raw$pathways[[1L]]$biology_synthesis$findings[1L]
  result <- clinical_ai_interpret(one$context, one$evidence,
    InterMetaAIProvider(function(messages, schema) one$raw), max_records_per_pathway = 1L)
  expect_match(result$record_selection$provider_visible_sha256, "^[a-f0-9]{64}$")
  modified <- fixture$raw
  modified$record_selection <- result$record_selection
  modified$context_selection <- result$context_selection
  modified$record_selection$max_records_per_pathway <- 2L
  expect_error(clinical_ai_report(fixture$context, fixture$evidence, interpretation = modified),
    "Provider-visible evidence or selection metadata changed")
  modified <- result
  modified$record_selection$sent_record_rows <- 2L
  expect_error(clinical_ai_report(fixture$context, fixture$evidence, interpretation = modified),
    "selection metadata changed")
  # Even internally consistent replacement metadata cannot retain the original
  # payload hash after expanding the selected records.
  selected <- synthesis_fn(".InterMetaAISelectRecords")(
    synthesis_fn(".InterMetaAISafeContext")(fixture$context),
    synthesis_fn(".InterMetaAISafeEvidence")(fixture$evidence), 2L)
  replacement <- synthesis_fn(".InterMetaAIPromptEvidence")(
    synthesis_fn(".InterMetaAISafeEvidence")(fixture$evidence), selected)$provenance
  replacement$provider_visible_sha256 <- result$record_selection$provider_visible_sha256
  modified <- fixture$raw
  modified$record_selection <- replacement
  modified$context_selection <- result$context_selection
  expect_error(clinical_ai_report(fixture$context, fixture$evidence, interpretation = modified),
    "Provider-visible evidence or selection metadata changed")
  old <- result
  old$record_selection$provider_visible_sha256 <- NULL
  report <- clinical_ai_report(fixture$context, fixture$evidence, interpretation = old)
  expect_true(any(grepl("no provider-visible payload fingerprint", report$warnings)))
})

test_that("missing synthesis is unassessed, not an upgrade from having claims", {
  fixture <- synthesis_fixture()
  fixture$raw$pathways[[1L]]$biology_synthesis <- NULL
  result <- synthesis_interpret(fixture)
  expect_identical(result$pathways$Biological.Relevance, "contextual_evidence")
  expect_true(all(result$pathway_synthesis$Biological.Explanation == "not_assessed"))
  expect_equal(nrow(result$biology_findings), 0L)
  fixture <- synthesis_fixture()
  fixture$raw$pathways[[1L]]$biology_synthesis$findings <- list()
  result <- synthesis_interpret(fixture)
  expect_identical(result$pathway_synthesis$Biological.Explanation[1L], "insufficient_evidence")
  expect_false(grepl("proposed cross-study", result$pathway_synthesis$Mechanism.Summary[1L]))
})

test_that("biology findings require existing pathway citations and abstract grounding", {
  fixture <- synthesis_fixture()
  fixture$raw$pathways[[1L]]$biology_synthesis$findings[[1L]]$pmid <- "99999999"
  expect_error(synthesis_interpret(fixture), "must be cited and assessed")
  fixture <- synthesis_fixture()
  fixture$raw$pathways[[1L]]$biology_synthesis$findings[[1L]]$abstract_quote <- "Invented mechanism."
  fixture$evidence$records$Title <- "Invented mechanism."
  expect_error(synthesis_interpret(fixture), "not an exact excerpt")
  fixture$raw$pathways[[1L]]$biology_synthesis$findings[[1L]]$abstract_quote <- ""
  expect_identical(synthesis_interpret(fixture)$pathway_synthesis$Biological.Explanation[1L], "insufficient_evidence")
  fixture <- synthesis_fixture()
  fixture$raw$pathways[[1L]]$biology_synthesis$findings <- rep(fixture$raw$pathways[[1L]]$biology_synthesis$findings, 2L)
  expect_error(synthesis_interpret(fixture), "Duplicate pathway/PMID")
  fixture <- synthesis_fixture()
  fixture$raw$pathways[[1L]]$assessments <- NULL
  expect_error(synthesis_interpret(fixture), "provided structured evidence assessment")
})

test_that("biology narrative cannot bypass inline citation and enum validation", {
  for (field in c("mechanism_summary", "limitations")) {
    fixture <- synthesis_fixture()
    fixture$raw$pathways[[1L]]$biology_synthesis[[field]] <- "PMID 99999999 proves the mechanism."
    expect_error(synthesis_interpret(fixture), "Inline citation")
  }
  fixture <- synthesis_fixture()
  fixture$raw$pathways[[1L]]$biology_synthesis$findings[[1L]]$claim <- "See https://evil.example."
  expect_error(synthesis_interpret(fixture), "Inline citation")
  for (field in names(synthesis_fn(".InterMetaAIBiologyEnums")())) {
    fixture <- synthesis_fixture()
    fixture$raw$pathways[[1L]]$biology_synthesis$findings[[1L]][[field]] <- "invented"
    expect_error(synthesis_interpret(fixture), "Invalid biology finding")
  }
})

test_that("synthesis, source attribution and evidence boundaries survive repeated reports", {
  fixture <- synthesis_fixture()
  attr(fixture$raw, "intermeta_assessment_source") <- "assistant_abstract_audit"
  result <- synthesis_interpret(fixture)
  report <- clinical_ai_report(fixture$context, fixture$evidence, interpretation = result)
  second <- clinical_ai_report(fixture$context, fixture$evidence, interpretation = report$interpretation)
  expect_identical(second$pathway_synthesis, result$pathway_synthesis)
  expect_identical(second$biology_findings, result$biology_findings)
  expect_identical(second$pathway_synthesis$Assessment.Source, c("assistant_abstract_audit", "not_assessed"))
  expect_identical(second$interpretation$record_selection, result$record_selection)
  expect_identical(second$interpretation$context_selection, result$context_selection)
  changed <- fixture$context
  changed$analysis$disease <- "breast cancer"
  expect_error(clinical_ai_report(changed, fixture$evidence, interpretation = second$interpretation), "Analysis context changed")
  changed <- fixture$context
  changed$pathways$Difference[1L] <- -20
  expect_error(clinical_ai_report(changed, fixture$evidence, interpretation = second$interpretation), "Analysis context changed")
  changed_evidence <- fixture$evidence
  changed_evidence$records$Abstract <- "Different findings."
  expect_error(clinical_ai_report(fixture$context, changed_evidence, interpretation = second$interpretation), "Evidence records changed")
})

test_that("reports lead with pathway synthesis and preserve complete proposed mechanisms for audit", {
  fixture <- synthesis_fixture()
  report <- clinical_ai_report(fixture$context, fixture$evidence, interpretation = fixture$raw)
  markdown <- paste(synthesis_fn(".InterMetaAIRenderMarkdown")(report), collapse = "\n")
  expect_match(markdown, "Pathway-level biological synthesis", fixed = TRUE)
  expect_match(markdown, "proposed cross-study biological mechanism", fixed = TRUE)
  expect_match(markdown, "not independently verified", fixed = TRUE)
  expect_match(markdown, "Supporting biological findings (supplement)", fixed = TRUE)
  expect_lt(regexpr("Pathway-level biological synthesis", markdown)[1L],
             regexpr("Endpoint-specific directional assessment", markdown)[1L])
  skeleton <- clinical_ai_report(fixture$context, fixture$evidence)
  expect_equal(nrow(skeleton$pathway_synthesis), 2L)
  expect_true(all(skeleton$pathway_synthesis$Biological.Explanation == "not_assessed"))
})
