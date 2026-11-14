explore_query_context <- function() list(
  analysis = list(disease = "prostate cancer", endpoint = "BCR"),
  pathways = data.frame(Pathway = "Oxidative phosphorylation", Direction = "unfavorable"))

paged_xml <- function(ids) paste0("<PubmedArticleSet>", paste0(vapply(ids, function(id)
  paste0('<PubmedArticle><MedlineCitation><PMID>', id,
    '</PMID><Article><ArticleTitle>Respiration in prostate cancer</ArticleTitle>',
    '<Abstract><AbstractText>Prostate cancer cells use respiration.</AbstractText></Abstract>',
    '</Article></MedlineCitation></PubmedArticle>'), character(1)), collapse = ""),
  "</PubmedArticleSet>")

paged_mock <- function(count = 5L, mode = "ok") {
  seen <- new.env(parent = emptyenv()); seen$params <- list()
  fun <- function(endpoint, params) {
    url <- endpoint
    query <- params
    seen$params[[length(seen$params) + 1L]] <- query
    if (grepl("esearch", url)) {
      if (mode == "error") stop("Unavailable")
      offset <- query$retstart
      ids <- if (count == 0L) character() else as.character(seq_len(count))
      page <- utils::head(ids[seq_along(ids) > offset], query$retmax)
      if (offset > 0L && mode == "duplicate") page[1] <- "1"
      if (offset > 0L && mode == "short") page <- character()
      n <- if (offset > 0L && mode == "changed") count + 1L else count
      body <- jsonlite::toJSON(list(esearchresult = list(count = as.character(n), idlist = page)), auto_unbox = TRUE)
    } else {
      ids <- strsplit(query$id, ",", fixed = TRUE)[[1]]
      if (mode == "missing") ids <- utils::head(ids, -1L)
      body <- paged_xml(ids)
    }
    list(status_code = 200L, body = body)
  }
  list(request = fun, seen = seen)
}

test_that("exploration queries remove endpoint gate but preserve disease and scope", {
  ctx <- explore_query_context()
  q <- InterMetaBuildPubMedQueries(ctx)
  expect_true(all(q$purpose == "biology_exploration"))
  expect_true(all(q$query_version == "3.1"))
  expect_identical(q$layer, c("biology", "biology_process", "review"))
  expect_false(any(grepl("BCR|survival|recurrence|Prognosis", q$query)))
  expect_true(all(grepl('"Prostatic Neoplasms"[mh]', q$query, fixed = TRUE)))
  ctx$analysis$endpoint <- "unknown"
  expect_identical(q$query, InterMetaBuildPubMedQueries(ctx)$query)
  old <- InterMetaBuildPubMedQueries(ctx, purpose = "prognostic_validation")
  expect_true(all(old$query_version == "2.1"))
  expect_false("purpose" %in% names(old))
  expect_identical(old$layer, c("prognosis", "biology", "review"))
  ctx$pathways$Pathway <- "Taurine and hypotaurine metabolism"
  component <- InterMetaBuildPubMedQueries(ctx)
  expect_setequal(component$layer, c("biology", "biology_process", "review",
    "component_biology", "component_process", "component_review"))
  expect_true(all(component$term_scope[startsWith(component$layer, "component")] == "component"))
})

test_that("paged retrieval records sort, pages, order, hashes and validated cache", {
  mock <- paged_mock(); ctx <- explore_query_context()
  q <- InterMetaBuildPubMedQueries(ctx)[1L, , drop = FALSE]
  cache <- tempfile(); dir.create(cache)
  withr::defer(unlink(cache, recursive = TRUE))
  run <- function(...) InterMetaPubMedSearch(ctx, queries = q, retmax = 5L, page_size = 2L,
    sort = "relevance", request_fun = mock$request, sleep_fun = function(x) NULL,
    api_key = "", cache_dir = cache, ...)
  a <- run()
  expect_identical(a$searches$api_status, "ok")
  expect_identical(a$records$PMID, as.character(1:5))
  expect_equal(a$searches$pages_retrieved, 3L)
  safe <- InterMeta:::.InterMetaAISafeEvidence(a)
  expect_identical(safe$searches$pages_retrieved, a$searches$pages_retrieved)
  selected <- InterMeta:::.InterMetaAISelectRecords(ctx, safe, 30L, purpose = "biology_exploration")
  payload <- InterMeta:::.InterMetaAIPromptEvidence(safe, selected, purpose = "biology_exploration")
  expect_identical(payload$payload$searches$pages_retrieved, a$searches$pages_retrieved)
  expect_equal(vapply(mock$seen$params[1:3], `[[`, integer(1), "retstart"), c(0L, 2L, 4L))
  expect_equal(a$searches$sort, "relevance")
  expect_true(all(a$records$purpose == "biology_exploration"))
  expect_false(a$searches$truncated)
  calls <- length(mock$seen$params)
  b <- run()
  expect_true(b$searches$cached)
  expect_identical(length(mock$seen$params), calls)
  path <- list.files(cache, full.names = TRUE)[1]
  tampered <- readRDS(path); tampered$query <- "tampered"; saveRDS(tampered, path)
  expect_warning(c <- run(), "invalid paged")
  expect_false(c$searches$cached)
})

test_that("paged failures, successful zero hits and partial results remain distinct", {
  ctx <- explore_query_context(); q <- InterMetaBuildPubMedQueries(ctx)[1L, , drop = FALSE]
  for (mode in c("error", "duplicate", "short", "changed", "missing")) {
    mock <- paged_mock(mode = mode)
    expect_warning(ev <- InterMetaPubMedSearch(ctx, queries = q, retmax = 5L,
      page_size = 2L, request_fun = mock$request, sleep_fun = function(x) NULL, api_key = ""), "incomplete or failed")
    expect_false(ev$searches$api_status == "ok")
    if (mode == "error") {
      expect_true(is.na(ev$searches$total_count))
      expect_true(is.na(ev$searches$records_fetched))
    } else expect_gt(nrow(ev$records), 0L)
  }
  zero <- paged_mock(count = 0L)
  ev <- InterMetaPubMedSearch(ctx, queries = q, retmax = 5L,
    page_size = 2L, request_fun = zero$request, sleep_fun = function(x) NULL, api_key = "")
  expect_identical(ev$searches$api_status, "ok")
  expect_equal(ev$searches$total_count, 0)
  expect_equal(ev$searches$records_fetched, 0L)
  expect_equal(nrow(ev$records), 0L)
})

test_that("paged caps and explicit purpose mismatches fail before transport", {
  ctx <- explore_query_context(); q <- InterMetaBuildPubMedQueries(ctx)
  blocked <- function(...) stop("SHOULD_NOT_REQUEST")
  expect_error(InterMetaPubMedSearch(ctx, q, purpose = "prognostic_validation", request_fun = blocked), "purpose")
  expect_error(InterMetaPubMedSearch(ctx, q, retmax = 10001L, page_size = 50L, request_fun = blocked), "10000")
  expect_error(InterMetaPubMedSearch(ctx, q, page_size = 0, request_fun = blocked), "page_size")
  expect_error(InterMetaPubMedSearch(ctx, q, sort = "bogus", request_fun = blocked), "arg")
  bad <- q; bad$purpose[1] <- NA_character_
  expect_error(InterMetaPubMedSearch(ctx, bad, request_fun = blocked), "purpose")
})
