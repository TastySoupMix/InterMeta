legacy_pubmed_queries <- function(...) InterMeta::InterMetaBuildPubMedQueries(
  ..., purpose = "prognostic_validation")

pubmed_test_context <- function(with_genes = TRUE) {
  pathways <- data.frame(
    Pathway = c("Glycolysis", "Fatty acid oxidation"),
    Direction = c("unfavorable", "favorable"),
    stringsAsFactors = FALSE
  )
  if (with_genes) {
    pathways$Driver.Genes <- I(list(c("HK2", "LDHA", "PKM"), character()))
  }
  list(
    analysis = list(disease = "prostate cancer", endpoint = "overall survival"),
    pathways = pathways
  )
}

pubmed_test_esearch <- paste0(
  '{"esearchresult":{"count":"2","retmax":"1",',
  '"retstart":"0","idlist":["12345"]}}'
)

pubmed_test_zero_esearch <- paste0(
  '{"esearchresult":{"count":"0","retmax":"0",',
  '"retstart":"0","idlist":[]}}'
)

pubmed_test_efetch <- paste0(
  '<?xml version="1.0"?><PubmedArticleSet><PubmedArticle>',
  '<MedlineCitation><PMID>12345</PMID><Article>',
  '<Journal><JournalIssue><PubDate><Year>2024</Year></PubDate></JournalIssue>',
  '<Title>Journal of Metabolism</Title></Journal>',
  '<ArticleTitle>Glycolysis and prostate cancer survival</ArticleTitle>',
  '<Abstract><AbstractText Label="BACKGROUND">Background text.</AbstractText>',
  '<AbstractText Label="RESULTS">Result text.</AbstractText></Abstract>',
  '<PublicationTypeList><PublicationType>Journal Article</PublicationType>',
  '<PublicationType>Review</PublicationType></PublicationTypeList>',
  '</Article><MeshHeadingList><MeshHeading>',
  '<DescriptorName>Prostatic Neoplasms</DescriptorName>',
  '</MeshHeading></MeshHeadingList></MedlineCitation>',
  '<PubmedData><ArticleIdList>',
  '<ArticleId IdType="doi">10.1000/intermeta.1</ArticleId>',
  '</ArticleIdList></PubmedData></PubmedArticle></PubmedArticleSet>'
)

test_that("PubMed queries are layered, field-tagged, and injection-safe", {
  context <- pubmed_test_context()
  context$analysis$disease <- 'prostate cancer" OR all[sb]'
  queries <- legacy_pubmed_queries(
    context,
    disease_mesh = "Prostatic Neoplasms",
    pathway_synonyms = list(Glycolysis = c("glycolytic pathway", "Embden Meyerhof")),
    max_genes = 2L
  )

  expect_equal(nrow(queries), 7L)
  expect_equal(
    queries$layer[queries$pathway == "Glycolysis"],
    c("prognosis", "biology", "review", "gene_bridge")
  )
  expect_equal(
    queries$layer[queries$pathway == "Fatty acid oxidation"],
    c("prognosis", "biology", "review")
  )
  expect_true(all(grepl('NOT "Retracted Publication"\\[pt\\]', queries$query)))
  expect_true(all(grepl('"Prostatic Neoplasms"\\[mh\\]', queries$query)))
  expect_true(all(grepl('\\[tiab\\]', queries$query)))
  expect_false(any(grepl("all\\[sb\\]", queries$query, fixed = FALSE)))
  gene_query <- queries$query[queries$layer == "gene_bridge"]
  expect_match(gene_query, '"HK2"\\[tiab\\]')
  expect_match(gene_query, '"LDHA"\\[tiab\\]')
  expect_false(grepl('"PKM"\\[tiab\\]', gene_query))
})

test_that("unknown endpoints are omitted and BCR endpoints use explicit terms", {
  context <- pubmed_test_context(with_genes = FALSE)
  for (endpoint in c("unknown", " UNKNOWN ")) {
    context$analysis$endpoint <- endpoint
    queries <- legacy_pubmed_queries(context)
    expect_false(any(grepl('"unknown"[tiab]', tolower(queries$query), fixed = TRUE)))
    prognosis <- queries[queries$layer == "prognosis", , drop = FALSE]
    expect_true(all(grepl('"Prognosis"[mh]', prognosis$query, fixed = TRUE)))
    expect_true(all(grepl("survival[tiab]", prognosis$query, fixed = TRUE)))
    expect_true(all(grepl("recurrence[tiab]", prognosis$query, fixed = TRUE)))
    expect_true(all(grepl("Endpoint 'unknown' was omitted", queries$query_warning,
                          fixed = TRUE)))
    expect_identical(queries$endpoint, rep(trimws(endpoint), nrow(queries)))
  }
  for (endpoint in c("BCR", "biochemical recurrence", "BCR (biochemical recurrence)")) {
    context$analysis$endpoint <- endpoint
    queries <- legacy_pubmed_queries(context)
    prognosis <- queries$query[queries$layer == "prognosis"]
    for (term in c("BCR", "biochemical recurrence", "biochemical relapse")) {
      expect_true(all(grepl(paste0('"', term, '"[tiab]'), prognosis, fixed = TRUE)))
    }
    expect_true(all(grepl("prognos*[tiab]", prognosis, fixed = TRUE)))
    expect_true(all(grepl("survival[tiab]", prognosis, fixed = TRUE)))
  }
})

test_that("v2 dictionary preserves the Top5 labels and separates components", {
  pathways <- c(
    "Oxidative phosphorylation", "Taurine and hypotaurine metabolism",
    "Inositol phosphate metabolism", "Glycerophospholipid metabolism",
    "Lysine degradation"
  )
  context <- list(
    analysis = list(disease = "prostate cancer", endpoint = "BCR"),
    pathways = data.frame(Pathway = pathways, Direction = "unfavorable")
  )
  original <- context
  queries <- legacy_pubmed_queries(context)
  expect_identical(context, original)
  expect_identical(unique(queries$pathway), pathways)
  expect_true(all(queries$direction == "unfavorable"))
  expect_equal(nrow(queries), 27L)
  expect_true(all(queries$query_version == "2.1"))
  expect_equal(sum(queries$term_scope == "pathway"), 15L)
  expect_equal(sum(queries$term_scope == "component"), 12L)
  expect_true(all(grepl('NOT "Retracted Publication"[pt]', queries$query, fixed = TRUE)))

  expansions <- list(
    "Oxidative phosphorylation" = "OXPHOS",
    "Lysine degradation" = "lysine catabolism"
  )
  components <- list(
    "Taurine and hypotaurine metabolism" = c("taurine", "hypotaurine"),
    "Inositol phosphate metabolism" = c("inositol phosphate", "inositol phosphates"),
    "Glycerophospholipid metabolism" = c("glycerophospholipid", "glycerophospholipids"),
    "Lysine degradation" = c(
      "saccharopine", "aminoadipate", "aminoadipic acid", "pipecolate"
    )
  )
  for (pathway in pathways) {
    pathway_rows <- queries[queries$pathway == pathway &
                              queries$term_scope == "pathway", , drop = FALSE]
    expect_identical(pathway_rows$layer, c("prognosis", "biology", "review"))
    expect_identical(pathway_rows$search_terms,
                     rep(paste(c(pathway, expansions[[pathway]]), collapse = ";"), 3L))
    expect_true(all(grepl(paste0('"', pathway, '"[tiab]'),
                         pathway_rows$query, fixed = TRUE)))
    if (pathway %in% names(expansions)) {
      expect_true(all(grepl(paste0('"', expansions[[pathway]], '"[tiab]'),
                           pathway_rows$query, fixed = TRUE)))
      expect_true(all(grepl("not proof", pathway_rows$query_warning, fixed = TRUE)))
    }
    component_rows <- queries[queries$pathway == pathway &
                                queries$term_scope == "component", , drop = FALSE]
    if (!pathway %in% names(components)) {
      expect_equal(nrow(component_rows), 0L)
      next
    }
    expect_identical(component_rows$layer,
                     c("component_prognosis", "component_biology", "component_review"))
    expect_identical(component_rows$search_terms,
                     rep(paste(components[[pathway]], collapse = ";"), 3L))
    expect_false(any(grepl(paste0('"', pathway, '"[tiab]'),
                          component_rows$query, fixed = TRUE)))
    expect_true(all(grepl("not whole-pathway evidence", component_rows$query_warning,
                          fixed = TRUE)))
    for (component in components[[pathway]]) {
      clause <- paste0('"', component, '"[tiab]')
      expect_true(all(grepl(clause, component_rows$query, fixed = TRUE)))
      expect_false(any(grepl(clause, pathway_rows$query, fixed = TRUE)))
    }
    prognosis <- component_rows$query[component_rows$layer == "component_prognosis"]
    expect_match(prognosis, '"biochemical recurrence"[tiab]', fixed = TRUE)
  }
  without_components <- legacy_pubmed_queries(context, include_components = FALSE)
  expect_equal(nrow(without_components), 15L)
  expect_true(all(without_components$term_scope == "pathway"))
  expect_identical(without_components$query, queries$query[queries$term_scope == "pathway"])
  expect_identical(queries, legacy_pubmed_queries(context))
})

test_that("lysine degradation avoids the nonspecific bare lysine component query", {
  context <- list(
    analysis = list(disease = "prostate cancer", endpoint = "BCR"),
    pathways = data.frame(
      Pathway = "Lysine degradation", Direction = "unfavorable",
      stringsAsFactors = FALSE
    )
  )
  queries <- InterMetaBuildPubMedQueries(
    context, purpose = "biology_exploration"
  )
  component <- queries[queries$term_scope == "component", , drop = FALSE]

  expect_equal(nrow(component), 3L)
  expect_false(any(grepl('("lysine"\\[tiab\\])', component$query)))
  expect_true(all(grepl('"saccharopine"[tiab]', component$query, fixed = TRUE)))
  expect_true(all(grepl('"aminoadipate"[tiab]', component$query, fixed = TRUE)))
  expect_true(all(grepl('"pipecolate"[tiab]', component$query, fixed = TRUE)))
})

test_that("component mappings are opt-in extensions for unlisted pathways and sanitized", {
  context <- pubmed_test_context(with_genes = FALSE)
  original_queries <- legacy_pubmed_queries(context)
  expect_equal(nrow(original_queries), 6L)
  expect_true(all(original_queries$term_scope == "pathway"))
  expect_identical(original_queries$search_terms, rep(context$pathways$Pathway, each = 3L))
  custom <- data.frame(
    Pathway = c("Glycolysis", "Glycolysis"),
    Component = c('unsafe" OR all[sb]; glucose', "glucose")
  )
  queries <- legacy_pubmed_queries(context, pathway_components = custom)
  component_rows <- queries[queries$term_scope == "component", , drop = FALSE]
  expect_equal(nrow(component_rows), 3L)
  expect_true(all(component_rows$pathway == "Glycolysis"))
  expect_true(all(component_rows$search_terms == "unsafe OR all sb;glucose"))
  expect_true(all(grepl('"unsafe OR all sb"[tiab]', component_rows$query, fixed = TRUE)))
  expect_false(any(grepl("all[sb]", queries$query, fixed = TRUE)))
  expect_identical(queries$query[queries$term_scope == "pathway"], original_queries$query)
  expect_identical(
    legacy_pubmed_queries(context, pathway_components = custom,
                               include_components = FALSE),
    original_queries
  )
  expect_identical(
    legacy_pubmed_queries(context, NULL, NULL, 0L),
    legacy_pubmed_queries(context, max_genes = 0L)
  )
  for (invalid in list(NA, NULL, c(TRUE, FALSE), 1, "yes")) {
    expect_error(legacy_pubmed_queries(context, include_components = invalid),
                 "include_components")
  }
})

test_that("scope and version metadata survive retrieval and legacy query tables", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("xml2")
  context <- list(
    analysis = list(disease = "prostate cancer", endpoint = "BCR"),
    pathways = data.frame(Pathway = "Lysine degradation", Direction = "unfavorable")
  )
  queries <- legacy_pubmed_queries(context)
  requested <- character()
  request <- function(endpoint, params) {
    if (grepl("esearch", endpoint, fixed = TRUE)) {
      requested <<- c(requested, params$term)
      pubmed_test_esearch
    } else pubmed_test_efetch
  }
  evidence <- InterMetaPubMedSearch(
    context, queries = queries, api_key = "", request_fun = request,
    sleep_fun = function(seconds) NULL
  )
  columns <- c("term_scope", "search_terms", "query_version", "query_warning")
  expect_identical(requested, queries$query)
  for (table in c("queries", "searches", "records")) {
    expect_identical(evidence[[table]][, columns], queries[, columns])
  }
  legacy <- queries[, c("query_id", "pathway", "direction", "layer", "query")]
  legacy_evidence <- InterMetaPubMedSearch(
    context, queries = legacy, api_key = "", request_fun = request,
    sleep_fun = function(seconds) NULL
  )
  expect_identical(legacy_evidence$queries$term_scope, queries$term_scope)
  for (table in c("queries", "searches", "records")) {
    expect_true(all(legacy_evidence[[table]]$query_version == "legacy"))
    expect_true(all(is.na(legacy_evidence[[table]]$search_terms)))
  }
  empty <- InterMetaPubMedSearch(
    context, queries = queries, api_key = "",
    request_fun = function(endpoint, params) pubmed_test_zero_esearch,
    sleep_fun = function(seconds) NULL
  )
  expect_true(all(columns %in% names(empty$records)))
  expect_equal(nrow(empty$records), 0L)
  expect_identical(empty$searches[, columns], queries[, columns])
  queries$term_scope[[1L]] <- "whole organism"
  expect_error(InterMetaPubMedSearch(context, queries = queries, api_key = "",
                                    request_fun = request), "term_scope")
})

test_that("expanded v2 query text separates legacy cache keys without altering caches", {
  context <- list(
    analysis = list(disease = "prostate cancer", endpoint = "unknown"),
    pathways = data.frame(Pathway = "Oxidative phosphorylation", Direction = "unfavorable")
  )
  query <- legacy_pubmed_queries(context)$query[[1L]]
  legacy_query <- sub(' OR "OXPHOS"[tiab]', "", query, fixed = TRUE)
  expect_false(identical(.InterMetaPubMedCacheKey(query, 20L),
                         .InterMetaPubMedCacheKey(legacy_query, 20L)))
})

test_that("context pathway synonyms flow into default searches safely", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("xml2")
  context <- pubmed_test_context(with_genes = FALSE)
  context$pathways$Synonyms <- c(
    'glycolytic pathway | Embden-Meyerhof pathway; unsafe" OR all[sb]',
    "beta oxidation | fatty-acid catabolism"
  )
  requested_terms <- character()
  mock_zero <- function(endpoint, params) {
    requested_terms <<- c(requested_terms, params$term)
    list(status_code = 200L, body = pubmed_test_zero_esearch)
  }

  evidence <- InterMetaPubMedSearch(
    context,
    request_fun = mock_zero,
    sleep_fun = function(seconds) NULL
  )

  expect_equal(nrow(evidence$queries), 6L)
  glycolysis_queries <- evidence$queries$query[
    evidence$queries$pathway == "Glycolysis"
  ]
  expect_true(all(grepl('"glycolytic pathway"\\[tiab\\]', glycolysis_queries)))
  expect_true(all(grepl('"Embden-Meyerhof pathway"\\[tiab\\]', glycolysis_queries)))
  expect_true(all(grepl('"unsafe OR all sb"\\[tiab\\]', glycolysis_queries)))
  expect_false(any(grepl("all\\[sb\\]", glycolysis_queries)))
  expect_true(all(glycolysis_queries %in% requested_terms))
  expect_true(all(evidence$searches$api_status == "ok"))
  expect_true(all(evidence$searches$total_count == 0))
})

test_that("PubMed search parses traceable records without network access", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("xml2")
  context <- pubmed_test_context(with_genes = FALSE)
  queries <- legacy_pubmed_queries(context)
  queries <- queries[queries$pathway == "Glycolysis", , drop = FALSE]
  calls <- list()
  mock_request <- function(endpoint, params) {
    calls[[length(calls) + 1L]] <<- list(endpoint = endpoint, params = params)
    if (grepl("esearch", endpoint, fixed = TRUE)) {
      return(list(status_code = 200L, body = pubmed_test_esearch))
    }
    list(status_code = 200L, body = pubmed_test_efetch)
  }

  evidence <- InterMetaPubMedSearch(
    context,
    queries = queries,
    retmax = 1L,
    email = "maintainer@example.org",
    request_fun = mock_request,
    sleep_fun = function(seconds) invisible(seconds)
  )

  expect_s3_class(evidence, "InterMetaPubMedEvidence")
  expect_equal(nrow(evidence$searches), 3L)
  expect_true(all(evidence$searches$api_status == "ok"))
  expect_true(all(evidence$searches$total_count == 2))
  expect_true(all(evidence$searches$returned_pmids == 1L))
  expect_true(all(evidence$searches$records_fetched == 1L))
  expect_true(all(evidence$searches$truncated))
  expect_equal(evidence$provenance$retmax, 1L)
  expect_equal(evidence$provenance$queries_truncated, 3L)
  expect_true(evidence$provenance$truncation_present)
  expect_equal(nrow(evidence$records), 3L)
  expect_true(all(evidence$records$PMID == "12345"))
  expect_true(all(evidence$records$DOI == "10.1000/intermeta.1"))
  expect_true(all(evidence$records$Year == "2024"))
  expect_match(evidence$records$Abstract[[1L]], "BACKGROUND: Background text")
  expect_match(evidence$records$PublicationType[[1L]], "Review")
  expect_match(evidence$records$MeSH[[1L]], "Prostatic Neoplasms")
  expect_equal(nrow(evidence$response_hashes), 6L)
  expect_true(all(nchar(evidence$response_hashes$sha256) == 64L))
  expect_false(evidence$provenance$api_key_used)
  expect_true(all(vapply(calls, function(x) {
    identical(x$params$email, "maintainer@example.org") &&
      is.null(x$params$api_key)
  }, logical(1L))))
})

test_that("PubMed cache retains raw-response hashes and avoids repeat requests", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("xml2")
  context <- pubmed_test_context(with_genes = FALSE)
  queries <- legacy_pubmed_queries(context)
  queries <- queries[1L, , drop = FALSE]
  cache_dir <- tempfile("intermeta_pubmed_cache_")
  dir.create(cache_dir)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  call_count <- 0L
  mock_request <- function(endpoint, params) {
    call_count <<- call_count + 1L
    if (grepl("esearch", endpoint, fixed = TRUE)) {
      return(list(status_code = 200L, body = pubmed_test_esearch))
    }
    list(status_code = 200L, body = pubmed_test_efetch)
  }

  first <- InterMetaPubMedSearch(
    context, queries = queries, retmax = 1L, cache_dir = cache_dir,
    request_fun = mock_request, sleep_fun = function(seconds) NULL
  )
  expect_equal(call_count, 2L)
  cache_file <- list.files(cache_dir, pattern = "^pubmed_.*\\.rds$", full.names = TRUE)
  expect_length(cache_file, 1L)
  cache_value <- readRDS(cache_file)
  expect_identical(cache_value$query, queries$query[[1L]])
  expect_true(all(vapply(cache_value$raw_responses, function(x) {
    identical(
      digest::digest(charToRaw(enc2utf8(x$body)), algo = "sha256", serialize = FALSE),
      x$sha256
    )
  }, logical(1L))))

  second <- InterMetaPubMedSearch(
    context, queries = queries, retmax = 1L, cache_dir = cache_dir,
    request_fun = function(endpoint, params) stop("network should not be called"),
    sleep_fun = function(seconds) stop("rate limiter should not run for cache")
  )
  expect_equal(call_count, 2L)
  expect_false(first$searches$cached)
  expect_true(second$searches$cached)
  expect_equal(second$records$PMID, first$records$PMID)
  expect_true(all(second$response_hashes$cached))
  expect_true(first$provenance$cache_enabled)
  expect_null(first$provenance$cache_dir)
  expect_equal(nchar(first$provenance$cache_path_sha256), 64L)
  expect_false(grepl(normalizePath(cache_dir, winslash = "/"),
                     paste(capture.output(str(first)), collapse = "\n"),
                     fixed = TRUE))
})

test_that("API failures remain unknown rather than becoming zero evidence", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("xml2")
  context <- pubmed_test_context(with_genes = FALSE)
  queries <- legacy_pubmed_queries(context)[1L, , drop = FALSE]

  expect_warning(
    evidence <- InterMetaPubMedSearch(
      context,
      queries = queries,
      request_fun = function(endpoint, params) {
        list(status_code = 503L, body = "service unavailable")
      },
      sleep_fun = function(seconds) NULL
    ),
    "must not be interpreted as zero evidence"
  )
  expect_equal(evidence$searches$api_status, "error")
  expect_true(is.na(evidence$searches$total_count))
  expect_true(is.na(evidence$searches$returned_pmids))
  expect_true(is.na(evidence$searches$records_fetched))
  expect_match(evidence$searches$error, "HTTP status 503")
  expect_equal(nrow(evidence$records), 0L)
})

test_that("NCBI request delays reflect API-key limits", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("xml2")
  context <- pubmed_test_context(with_genes = FALSE)
  queries <- legacy_pubmed_queries(context)[1L, , drop = FALSE]
  mock_zero <- function(endpoint, params) {
    list(status_code = 200L, body = pubmed_test_zero_esearch)
  }
  delays <- numeric()

  no_key <- InterMetaPubMedSearch(
    context, queries = queries, api_key = "", request_fun = mock_zero,
    sleep_fun = function(seconds) delays <<- c(delays, seconds)
  )
  expect_equal(delays, 0.34)
  expect_equal(no_key$searches$retmax, 20L)
  expect_equal(no_key$provenance$retmax, 20L)
  expect_equal(no_key$searches$total_count, 0)
  expect_false(no_key$searches$truncated)
  expect_false(no_key$provenance$truncation_present)

  delays <- numeric()
  with_key <- InterMetaPubMedSearch(
    context, queries = queries, api_key = "secret-not-returned",
    request_fun = mock_zero,
    sleep_fun = function(seconds) delays <<- c(delays, seconds)
  )
  expect_equal(delays, 0.11)
  expect_true(with_key$provenance$api_key_used)
  expect_false("api_key" %in% names(with_key$provenance))
})

test_that("JSON/XML parsing errors cannot disclose an echoed API key", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("xml2")
  context <- pubmed_test_context(with_genes = FALSE)
  queries <- legacy_pubmed_queries(context)[1L, , drop = FALSE]
  secret <- "PRIVATE-NCBI-API-KEY-123"
  bodies <- list(
    malformed_json = paste0('{"esearchresult":', secret),
    json_error = paste0('{"esearchresult":{"error":"', secret, '"}}'),
    malformed_xml = paste0("<PubmedArticleSet><", secret),
    xml_error = paste0("<PubmedArticleSet><ERROR>", secret,
                       "</ERROR></PubmedArticleSet>")
  )
  for (case in names(bodies)) {
    body <- bodies[[case]]
    use_fetch <- grepl("xml", case)
    mock_request <- function(endpoint, params) {
      if (use_fetch && grepl("esearch", endpoint, fixed = TRUE)) {
        return(list(status_code = 200L, body = pubmed_test_esearch))
      }
      list(status_code = 200L, body = body)
    }
    expect_warning(
      evidence <- InterMetaPubMedSearch(
        context, queries = queries, api_key = secret,
        request_fun = mock_request, sleep_fun = function(seconds) NULL
      ),
      "must not be interpreted as zero evidence"
    )
    expect_identical(evidence$searches$api_status, "error", info = case)
    expect_false(any(grepl(secret, evidence$searches$error, fixed = TRUE)),
                 info = case)
    expect_false(grepl(secret, paste(capture.output(dput(evidence)),
                                    collapse = "\n"), fixed = TRUE), info = case)
    expect_true(any(evidence$response_hashes$body_redacted), info = case)
  }
})

test_that("successful response bodies are redacted before records and cache", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("xml2")
  context <- pubmed_test_context(with_genes = FALSE)
  queries <- legacy_pubmed_queries(context)[1L, , drop = FALSE]
  secret <- "PRIVATE-NCBI-API-KEY-123"
  cache_dir <- tempfile("intermeta_pubmed_redaction_")
  dir.create(cache_dir)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  echo_response <- sub("Glycolysis and prostate cancer survival", secret,
                       pubmed_test_efetch, fixed = TRUE)
  evidence <- InterMetaPubMedSearch(
    context, queries = queries, retmax = 1L, cache_dir = cache_dir,
    api_key = secret,
    request_fun = function(endpoint, params) {
      list(status_code = 200L, body = if (grepl("esearch", endpoint)) {
        pubmed_test_esearch
      } else echo_response)
    },
    sleep_fun = function(seconds) NULL
  )
  cache_file <- list.files(cache_dir, full.names = TRUE)
  cache_value <- readRDS(cache_file[[1L]])
  expect_identical(evidence$records$Title, "[REDACTED]")
  expect_true(any(evidence$response_hashes$body_redacted))
  expect_identical(evidence$provenance$response_hash_scope,
                   "stored_body_after_secret_redaction")
  expect_false(grepl(secret, paste(capture.output(dput(cache_value)),
                                  collapse = "\n"), fixed = TRUE))
  expect_true(.InterMetaPubMedValidCache(
    cache_value, cache_value$key, queries$query[[1L]], 1L,
    request_source = "custom_request_fun"
  ))
})

test_that("custom transport caches never populate official NCBI requests", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("xml2")
  skip_if_not_installed("httr2")
  context <- pubmed_test_context(with_genes = FALSE)
  queries <- legacy_pubmed_queries(context)[1L, , drop = FALSE]
  cache_dir <- tempfile("intermeta_pubmed_sources_")
  dir.create(cache_dir)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  custom <- InterMetaPubMedSearch(
    context, queries = queries, retmax = 1L, cache_dir = cache_dir,
    request_fun = function(endpoint, params) pubmed_test_zero_esearch,
    sleep_fun = function(seconds) NULL
  )
  official_calls <- 0L
  local_mocked_bindings(
    .InterMetaPubMedDefaultRequest = function(endpoint, params) {
      official_calls <<- official_calls + 1L
      if (grepl("esearch", endpoint)) pubmed_test_esearch else pubmed_test_efetch
    },
    .package = "InterMeta"
  )
  official <- InterMetaPubMedSearch(
    context, queries = queries, retmax = 1L, cache_dir = cache_dir,
    sleep_fun = function(seconds) NULL
  )
  expect_equal(official_calls, 2L)
  expect_false(official$searches$cached)
  expect_equal(nrow(official$records), 1L)
  expect_false(identical(custom$searches$cache_key, official$searches$cache_key))
  expect_identical(custom$provenance$request_source, "custom_request_fun")
  expect_identical(official$provenance$request_source, "ncbi_eutilities")
  expect_length(list.files(cache_dir, pattern = "\\.rds$"), 2L)
})

test_that("cache requires canonical integrity and agreement with raw responses", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("xml2")
  context <- pubmed_test_context(with_genes = FALSE)
  queries <- legacy_pubmed_queries(context)[1L, , drop = FALSE]
  cache_dir <- tempfile("intermeta_pubmed_integrity_")
  dir.create(cache_dir)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  calls <- 0L
  mock_request <- function(endpoint, params) {
    calls <<- calls + 1L
    if (grepl("esearch", endpoint)) pubmed_test_esearch else pubmed_test_efetch
  }
  invisible(InterMetaPubMedSearch(
    context, queries = queries, retmax = 1L, cache_dir = cache_dir,
    request_fun = mock_request, sleep_fun = function(seconds) NULL
  ))
  cache_file <- list.files(cache_dir, full.names = TRUE)[[1L]]
  original <- readRDS(cache_file)
  valid <- function(x) .InterMetaPubMedValidCache(
    x, original$key, queries$query[[1L]], 1L, "custom_request_fun"
  )
  expect_true(valid(original))
  expect_identical(original$cache_version, 2L)
  expect_equal(nchar(original$ids_sha256), 64L)
  expect_equal(nchar(original$records_sha256), 64L)

  tampered <- original
  tampered$records$Title <- "Invented article"
  expect_false(valid(tampered))
  tampered$records_sha256 <- .InterMetaPubMedHashRecords(tampered$records)
  expect_false(valid(tampered))
  saveRDS(tampered, cache_file)
  expect_warning(
    refreshed <- InterMetaPubMedSearch(
      context, queries = queries, retmax = 1L, cache_dir = cache_dir,
      request_fun = mock_request, sleep_fun = function(seconds) NULL
    ),
    "Ignoring an invalid PubMed cache entry"
  )
  expect_equal(calls, 4L)
  expect_false(refreshed$searches$cached)
  expect_identical(refreshed$records$Title,
                   "Glycolysis and prostate cancer survival")

  tampered <- original
  tampered$ids <- "54321"
  tampered$ids_sha256 <- .InterMetaPubMedHashIDs(tampered$ids)
  expect_false(valid(tampered))
  tampered <- original
  tampered$records$Title <- NULL
  expect_false(valid(tampered))
  tampered <- original
  tampered$records$PMID <- 12345
  expect_false(valid(tampered))
  tampered <- original
  tampered$raw_responses[[1L]]$body <- pubmed_test_zero_esearch
  expect_false(valid(tampered))
  tampered <- original
  tampered$raw_responses[[1L]]$stage <- NA_character_
  expect_false(valid(tampered))
  tampered <- original
  tampered$request_source <- "ncbi_eutilities"
  expect_false(valid(tampered))
})

test_that("default HTTP requests have finite timeout and transient retries", {
  skip_if_not_installed("httr2")
  request <- .InterMetaPubMedDefaultRequestObject(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    list(db = "pubmed", term = "test")
  )
  expect_equal(request$options$timeout_ms, 30000)
  expect_equal(request$policies$retry_max_tries, 3L)
  expect_equal(request$policies$retry_max_wait, 60)
  transient <- request$policies$retry_is_transient
  expect_true(transient(httr2::response(429L)))
  expect_true(transient(httr2::response(500L)))
  expect_true(transient(httr2::response(503L)))
  expect_false(transient(httr2::response(200L)))
  expect_false(transient(httr2::response(400L)))
  expect_false(transient(httr2::response(404L)))
})
