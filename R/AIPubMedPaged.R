# Paged public-literature retrieval. Legacy v2 caches remain untouched.
.InterMetaPubMedPagedCacheValid <- function(x, key, retmax, page_size, sort, query, request_source) {
  tryCatch({
    if (!is.list(x) || !identical(x$cache_version, 3L) || !identical(x$key, key) ||
        !identical(x$retmax, retmax) || !identical(x$page_size, page_size) ||
        !identical(x$sort, sort) || !identical(x$query, query) ||
        !identical(x$request_source, request_source) || !identical(x$api_status, "ok") ||
        !is.character(x$ids) || anyDuplicated(x$ids) || anyNA(x$ids) ||
        !is.list(x$raw_responses) || !length(x$raw_responses) ||
        !is.data.frame(x$records) || anyDuplicated(x$records$PMID) ||
        !identical(names(x$records), names(.InterMetaPubMedEmptyCoreRecords())) ||
        !identical(x$ids_sha256, .InterMetaPubMedHashIDs(x$ids)) ||
        !identical(x$records_sha256, .InterMetaPubMedHashRecords(x$records))) return(FALSE)
    raw_ok <- vapply(x$raw_responses, function(r) {
      is.list(r) && is.character(r$body) && length(r$body) == 1L && !is.na(r$body) &&
        identical(.InterMetaPubMedHashText(r$body), r$sha256) &&
        r$stage %in% c("esearch", "efetch") && length(r$status_code) == 1L &&
        is.finite(r$status_code) && r$status_code >= 200L && r$status_code < 300L &&
        is.logical(r$body_redacted) && length(r$body_redacted) == 1L && !is.na(r$body_redacted)
    }, logical(1))
    if (!all(raw_ok)) return(FALSE)
    search <- x$raw_responses[vapply(x$raw_responses, function(r) r$stage == "esearch", logical(1))]
    if (!length(search)) return(FALSE)
    ids <- character(); offset <- 0L
    for (i in seq_along(search)) {
      r <- search[[i]]
      expected <- as.integer(min(page_size, retmax - offset))
      if (!identical(r$retstart, offset) || !identical(r$requested, expected) ||
          !identical(r$batch, as.integer(i))) return(FALSE)
      p <- .InterMetaPubMedParseESearch(r$body)
      if (!identical(p$count, x$total_count) ||
          length(p$ids) != min(expected, max(0, p$count - offset))) return(FALSE)
      ids <- c(ids, p$ids)
      offset <- offset + expected
    }
    if (anyDuplicated(ids) || !identical(ids, x$ids) ||
        length(ids) != min(retmax, x$total_count)) return(FALSE)
    fetch <- x$raw_responses[vapply(x$raw_responses, function(r) r$stage == "efetch", logical(1))]
    rec <- if (length(fetch)) do.call(rbind, lapply(fetch, function(r)
      .InterMetaPubMedParseEFetch(r$body))) else .InterMetaPubMedEmptyCoreRecords()
    if (any(!rec$PMID %in% ids)) return(FALSE)
    rec <- rec[!duplicated(rec$PMID), , drop = FALSE]
    rec <- rec[match(ids, rec$PMID), , drop = FALSE]
    if (anyNA(rec$PMID)) return(FALSE)
    identical(.InterMetaPubMedHashRecords(rec), x$records_sha256)
  }, error = function(e) FALSE)
}

.InterMetaPubMedPagedSearch <- function(queries, retmax, page_size, sort,
    email, tool, api_key, cache_dir, refresh, request_fun, request_source,
    sleep_fun, delay) {
  if (retmax > 10000L) stop("Paged PubMed retrieval supports retmax <= 10000 per query.", call. = FALSE)
  if (is.null(page_size)) page_size <- min(100L, retmax)
  if (!is.numeric(page_size) || length(page_size) != 1L || is.na(page_size) ||
      !is.finite(page_size) || page_size < 1 || page_size > 10000 || page_size != floor(page_size))
    stop("page_size must be a positive integer <= 10000.", call. = FALSE)
  page_size <- as.integer(page_size)
  if (is.null(sort)) sort <- "relevance"
  sort <- match.arg(sort, c("relevance", "pub_date"))
  urls <- c(esearch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            efetch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi")
  search_rows <- record_rows <- vector("list", nrow(queries))
  hashes <- list()
  request <- function(stage, params) {
    params$tool <- tool
    if (!is.null(email)) params$email <- email
    if (nzchar(api_key)) params$api_key <- api_key
    .InterMetaPubMedRequest(request_fun, urls[[stage]], params, delay, sleep_fun, api_key)
  }
  raw_row <- function(stage, response, batch, retstart = NULL, requested = NULL) list(
    stage = stage, batch = as.integer(batch), retstart = retstart, requested = requested,
    status_code = response$status_code, body = response$body,
    body_redacted = response$body_redacted, sha256 = .InterMetaPubMedHashText(response$body))
  for (i in seq_len(nrow(queries))) {
    q <- queries[i, , drop = FALSE]
    key <- digest::digest(list(cache_version = 3L, query = q$query, retmax = retmax,
      page_size = page_size, sort = sort, request_source = request_source), algo = "sha256")
    path <- if (is.null(cache_dir)) NULL else file.path(cache_dir, paste0("pubmed_paged_", key, ".rds"))
    cached <- NULL
    if (!refresh && !is.null(path) && file.exists(path)) {
      cached <- tryCatch(readRDS(path), error = function(e) NULL)
      if (!.InterMetaPubMedPagedCacheValid(cached, key, retmax, page_size, sort, q$query, request_source)) {
        warning("Ignoring invalid paged PubMed cache for query ", q$query_id, call. = FALSE)
        cached <- NULL
      }
    }
    used_cache <- !is.null(cached)
    if (used_cache) {
      value <- cached
    } else {
      value <- list(cache_version = 3L, key = key, retmax = retmax, page_size = page_size,
        sort = sort, query = q$query, request_source = request_source, searched_at = .InterMetaPubMedNow(),
        total_count = NA_real_, ids = character(), records = .InterMetaPubMedEmptyCoreRecords(),
        api_status = "ok", error = NA_character_, raw_responses = list())
      offset <- 0L; page <- 0L
      repeat {
        page <- page + 1L
        requested <- as.integer(min(page_size, retmax - offset))
        got <- tryCatch(request("esearch", list(db = "pubmed", term = q$query,
          retmode = "json", retmax = requested, retstart = offset, sort = sort)), error = function(e) e)
        if (inherits(got, "error")) {
          value$api_status <- if (length(value$ids)) "partial_error" else "error"
          value$error <- paste("ESearch page failed:", conditionMessage(got)); break
        }
        value$raw_responses[[length(value$raw_responses) + 1L]] <-
          raw_row("esearch", got, page, offset, requested)
        parsed <- tryCatch(.InterMetaPubMedParseESearch(got$body), error = function(e) e)
        if (inherits(parsed, "error")) {
          value$api_status <- if (length(value$ids)) "partial_error" else "error"
          value$error <- paste("ESearch page parsing failed:", conditionMessage(parsed)); break
        }
        if (is.na(value$total_count)) value$total_count <- parsed$count
        if (!identical(parsed$count, value$total_count) ||
            length(parsed$ids) != min(requested, max(0, parsed$count - offset)) ||
            anyDuplicated(c(value$ids, parsed$ids))) {
          value$api_status <- "partial_error"
          value$error <- "ESearch pages were incomplete, duplicated, or changed during pagination."; break
        }
        value$ids <- c(value$ids, parsed$ids)
        offset <- offset + requested
        if (length(value$ids) >= min(retmax, value$total_count)) break
      }
      if (length(value$ids)) {
        batches <- split(value$ids, ceiling(seq_along(value$ids) / 200L))
        fetched <- list()
        for (j in seq_along(batches)) {
          got <- tryCatch(request("efetch", list(db = "pubmed", id = paste(batches[[j]], collapse = ","),
            retmode = "xml")), error = function(e) e)
          if (inherits(got, "error")) {
            value$api_status <- "partial_error"; value$error <- paste("EFetch failed:", conditionMessage(got)); next
          }
          value$raw_responses[[length(value$raw_responses) + 1L]] <- raw_row("efetch", got, j)
          parsed <- tryCatch(.InterMetaPubMedParseEFetch(got$body), error = function(e) e)
          if (inherits(parsed, "error")) {
            value$api_status <- "partial_error"; value$error <- paste("EFetch parsing failed:", conditionMessage(parsed)); next
          }
          fetched[[length(fetched) + 1L]] <- parsed[parsed$PMID %in% batches[[j]], , drop = FALSE]
        }
        if (length(fetched)) {
          rec <- do.call(rbind, fetched); rec <- rec[!duplicated(rec$PMID), , drop = FALSE]
          positions <- match(value$ids, rec$PMID)
          value$records <- rec[positions[!is.na(positions)], , drop = FALSE]
          rownames(value$records) <- NULL
        }
        if (nrow(value$records) < length(value$ids) && value$api_status == "ok") {
          value$api_status <- "partial"; value$error <- "EFetch metadata was missing for some requested PMIDs."
        }
      }
      value$error <- .InterMetaPubMedRedact(value$error, api_key)
      value$ids_sha256 <- .InterMetaPubMedHashIDs(value$ids)
      value$records_sha256 <- .InterMetaPubMedHashRecords(value$records)
      if (!is.null(path) && value$api_status == "ok") .InterMetaPubMedWriteCache(value, path)
    }
    known <- !is.na(value$total_count)
    s <- .InterMetaPubMedSearchRow(q, value$searched_at, retmax, key, used_cache,
      value$total_count, if (known) length(value$ids) else NA_integer_,
      if (known) nrow(value$records) else NA_integer_,
      if (is.na(value$total_count)) NA else value$total_count > length(value$ids),
      value$api_status, value$error)
    s$page_size <- page_size; s$sort <- sort
    s$pages_retrieved <- sum(vapply(value$raw_responses, function(r) r$stage == "esearch", logical(1)))
    if ("purpose" %in% names(q)) s$purpose <- q$purpose
    search_rows[[i]] <- s
    rec <- .InterMetaPubMedAttachQuery(value$records, q)
    if ("purpose" %in% names(q)) rec$purpose <- rep(q$purpose, nrow(rec))
    record_rows[[i]] <- rec
    for (r in value$raw_responses) hashes[[length(hashes) + 1L]] <- data.frame(
      query_id = q$query_id, stage = r$stage, batch = r$batch, status_code = r$status_code,
      sha256 = r$sha256, body_redacted = r$body_redacted, cached = used_cache,
      stringsAsFactors = FALSE)
  }
  searches <- do.call(rbind, search_rows); rownames(searches) <- NULL
  records <- do.call(rbind, record_rows); rownames(records) <- NULL
  response_hashes <- if (length(hashes)) do.call(rbind, hashes) else data.frame(
    query_id = character(), stage = character(), batch = integer(), status_code = integer(),
    sha256 = character(), body_redacted = logical(), cached = logical())
  if (any(searches$api_status != "ok")) warning("Some paged PubMed queries were incomplete or failed; inspect searches before interpretation.", call. = FALSE)
  structure(list(schema_version = "1.0", generated_at = .InterMetaPubMedNow(),
    search_date = format(Sys.Date(), "%Y-%m-%d"), queries = queries, searches = searches,
    records = records, response_hashes = response_hashes, provenance = list(
      provider = "NCBI E-utilities", request_source = request_source, database = "PubMed",
      tool = tool, email_supplied = !is.null(email), api_key_used = nzchar(api_key),
      cache_enabled = !is.null(cache_dir), refresh = refresh, rate_limit_seconds = delay,
      retmax = retmax, page_size = page_size, sort = sort, cache_version = 3L,
      queries_truncated = searches$query_id[!is.na(searches$truncated) & searches$truncated],
      truncation_present = any(searches$truncated, na.rm = TRUE),
      response_hash_scope = "stored_body_after_secret_redaction")), class = c("InterMetaPubMedEvidence", "list"))
}
