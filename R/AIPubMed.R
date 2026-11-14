# PubMed evidence retrieval for InterMeta AI reports.
#
# This module deliberately keeps literature retrieval separate from model-based
# interpretation. Query construction is deterministic, every record retains the
# query that retrieved it, and transport failures are represented as failures
# rather than as searches with zero results.

.InterMetaPubMedScalarString <- function(x, name, allow_empty = FALSE) {
  if (length(x) != 1L || is.na(x) || !is.character(x)) {
    stop(sprintf("`%s` must be a single non-missing character string.", name),
         call. = FALSE)
  }
  x <- trimws(x)
  if (!allow_empty && !nzchar(x)) {
    stop(sprintf("`%s` must not be empty.", name), call. = FALSE)
  }
  x
}

.InterMetaPubMedSafeText <- function(x) {
  x <- enc2utf8(as.character(x))
  x[is.na(x)] <- ""
  # PubMed has no general-purpose escaping syntax for field-tag delimiters.
  # Removing syntax characters before quoting prevents a user-provided disease,
  # pathway, or gene label from injecting a second query clause.
  x <- gsub("[\\x00-\\x1f\\x7f]", " ", x, perl = TRUE)
  x <- gsub("[\\\"\\[\\]\\{\\}]", " ", x, perl = TRUE)
  x <- gsub("[()]", " ", x, perl = TRUE)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  trimws(x)
}

.InterMetaPubMedRedact <- function(x, secrets = character()) {
  x <- as.character(x)
  secrets <- unique(as.character(secrets))
  secrets <- secrets[!is.na(secrets) & nzchar(secrets)]
  for (secret in secrets) {
    x <- gsub(secret, "[REDACTED]", x, fixed = TRUE)
  }
  x
}

.InterMetaPubMedTaggedPhrase <- function(x, tag) {
  x <- .InterMetaPubMedSafeText(x)
  x <- x[nzchar(x)]
  if (!length(x)) {
    return(character())
  }
  sprintf('"%s"[%s]', x, tag)
}

.InterMetaPubMedOrClause <- function(terms) {
  terms <- unique(terms[nzchar(terms)])
  if (!length(terms)) {
    return("")
  }
  paste0("(", paste(terms, collapse = " OR "), ")")
}

.InterMetaPubMedLookupTerms <- function(mapping, key, value_names = NULL) {
  if (is.null(mapping)) {
    return(character())
  }

  if (is.data.frame(mapping)) {
    if (!nrow(mapping)) {
      return(character())
    }
    key_col <- intersect(c("Pathway", "pathway", "Disease", "disease"),
                         names(mapping))
    value_col <- intersect(c(value_names, "Synonym", "synonym", "Term", "term",
                             "Gene", "gene"), names(mapping))
    if (!length(key_col) || !length(value_col)) {
      stop("A term mapping data frame must contain a pathway/disease column and a term column.",
           call. = FALSE)
    }
    keep <- !is.na(mapping[[key_col[[1L]]]]) &
      as.character(mapping[[key_col[[1L]]]]) == key
    return(as.character(mapping[[value_col[[1L]]]][keep]))
  }

  if (is.list(mapping)) {
    if (!is.null(names(mapping)) && key %in% names(mapping)) {
      return(as.character(unlist(mapping[[key]], use.names = FALSE)))
    }
    if (is.null(names(mapping)) && length(mapping) == 1L) {
      return(as.character(unlist(mapping[[1L]], use.names = FALSE)))
    }
    return(character())
  }

  if (is.character(mapping)) {
    if (!is.null(names(mapping))) {
      return(as.character(mapping[names(mapping) == key]))
    }
    return(as.character(mapping))
  }

  stop("Term mappings must be character vectors, named lists, or data frames.",
       call. = FALSE)
}

.InterMetaPubMedSplitGenes <- function(x) {
  if (is.null(x) || !length(x)) {
    return(character())
  }
  x <- as.character(unlist(x, use.names = FALSE))
  x <- unlist(strsplit(x, "[,;|/]", perl = TRUE), use.names = FALSE)
  x <- trimws(x)
  unique(x[!is.na(x) & nzchar(x)])
}

.InterMetaPubMedSplitSynonyms <- function(x) {
  if (is.null(x) || !length(x)) {
    return(character())
  }
  x <- as.character(unlist(x, use.names = FALSE))
  # InterMetaAIContext serializes pathway synonyms with ` | `. Semicolons are
  # also accepted for hand-authored contexts; commas are intentionally retained
  # because they can be part of a legitimate pathway name.
  x <- unlist(strsplit(x, "\\s*\\|\\s*|\\s*;\\s*", perl = TRUE),
              use.names = FALSE)
  x <- trimws(x)
  unique(x[!is.na(x) & nzchar(x)])
}

.InterMetaPubMedContextSynonyms <- function(context, pathway, direction) {
  pathway_table <- context$pathways
  synonym_cols <- intersect(
    c("Synonyms", "synonyms", "Pathway.Synonyms", "pathway_synonyms"),
    names(pathway_table)
  )
  if (!length(synonym_cols)) {
    return(character())
  }
  keep <- !is.na(pathway_table$Pathway) &
    as.character(pathway_table$Pathway) == pathway
  if ("Direction" %in% names(pathway_table)) {
    keep <- keep & !is.na(pathway_table$Direction) &
      as.character(pathway_table$Direction) == direction
  }
  column <- pathway_table[[synonym_cols[[1L]]]]
  values <- if (is.list(column)) column[keep] else column[keep]
  .InterMetaPubMedSplitSynonyms(values)
}

.InterMetaPubMedBuiltInTerms <- function(pathway) {
  # This small dictionary expands retrieval, not the biological definition of
  # an enriched pathway. In particular, component names must never enter the
  # whole-pathway synonym clause or silently replace a ranked pathway label.
  key <- tolower(trimws(pathway))
  synonyms <- switch(key,
    "oxidative phosphorylation" = "OXPHOS",
    "lysine degradation" = "lysine catabolism",
    character()
  )
  components <- switch(key,
    "taurine and hypotaurine metabolism" = c("taurine", "hypotaurine"),
    "inositol phosphate metabolism" = c("inositol phosphate", "inositol phosphates"),
    "glycerophospholipid metabolism" = c("glycerophospholipid", "glycerophospholipids"),
    # A bare "lysine" query retrieves protein modification, drug-conjugate and
    # nutrition papers that are not evidence about lysine degradation.  Use
    # pathway-specific intermediates instead; screening still decides whether
    # an individual record is relevant to the target pathway and disease.
    "lysine degradation" = c(
      "saccharopine", "aminoadipate", "aminoadipic acid", "pipecolate"
    ),
    character()
  )
  list(synonyms = synonyms, components = components)
}

.InterMetaPubMedDriverGenes <- function(context, pathways, max_genes) {
  out <- stats::setNames(vector("list", length(pathways)), pathways)
  pathway_table <- context$pathways
  candidate_cols <- intersect(
    c("Driver.Genes", "DriverGenes", "driver_genes", "Genes", "genes",
      "TopGenes", "top_genes"),
    names(pathway_table)
  )

  if (length(candidate_cols)) {
    col <- pathway_table[[candidate_cols[[1L]]]]
    for (pathway in pathways) {
      idx <- which(as.character(pathway_table$Pathway) == pathway)
      values <- if (is.list(col)) col[idx] else col[idx]
      out[[pathway]] <- c(out[[pathway]], .InterMetaPubMedSplitGenes(values))
    }
  }

  external <- context$driver_genes
  if (is.null(external)) {
    external <- context$pathway_drivers
  }
  if (!is.null(external)) {
    for (pathway in pathways) {
      out[[pathway]] <- c(
        out[[pathway]],
        .InterMetaPubMedLookupTerms(
          external, pathway,
          value_names = c("DriverGene", "driver_gene", "Symbol", "symbol")
        )
      )
    }
  }

  lapply(out, function(x) {
    x <- .InterMetaPubMedSplitGenes(x)
    if (!length(x) || max_genes == 0L) {
      return(character())
    }
    utils::head(x, max_genes)
  })
}

#' Build reproducible PubMed queries for an InterMeta AI context
#'
#' @param context A list produced by `InterMetaAIContext()`. It must contain a
#'   `pathways` data frame with `Pathway` and `Direction`, and an `analysis`
#'   list with `disease` and `endpoint`.
#' @param disease_mesh Optional MeSH disease term or disease-to-term mapping.
#' @param pathway_synonyms Optional pathway-to-synonym mapping supplied as a
#'   named list, named character vector, or data frame. These terms supplement a
#'   conservative built-in retrieval dictionary; they do not establish strict
#'   biological equivalence.
#' @param max_genes Maximum number of driver genes used in an optional
#'   `gene_bridge` query for each pathway.
#' @param pathway_components Optional pathway-to-component mapping, in the same
#'   forms as `pathway_synonyms`. Component terms are searched independently in
#'   `component_prognosis`, `component_biology`, and `component_review` layers,
#'   never as whole-pathway synonyms. Built-in component candidates are included
#'   for selected metabolic pathways. Generic single-amino-acid names are
#'   avoided when they would produce a predominantly unrelated result set.
#' @param include_components Whether to include component-level queries.
#' @param purpose Report goal: exploratory tumor biology (default), or the
#'   legacy endpoint-specific prognostic validation strategy.
#'
#' @return A data frame containing deterministic PubMed queries. `pathway`
#'   retains the original pathway label. `term_scope` distinguishes `pathway`
#'   from `component`; `search_terms` contains the semicolon-delimited, sanitized
#'   pathway or component terms actually queried. `query_version` is `"2.1"`
#'   for endpoint-oriented retrieval and `"3.1"` for biology exploration;
#'   `query_warning` records interpretation limits. The placeholder endpoint
#'   `"unknown"` is omitted from search terms, while generic prognosis terms are
#'   retained. BCR endpoints also include biochemical recurrence and relapse.
#' @export
InterMetaBuildPubMedQueries <- function(context, disease_mesh = NULL,
                                         pathway_synonyms = NULL,
                                         max_genes = 5L,
                                         pathway_components = NULL,
                                         include_components = TRUE,
                                         purpose = c("biology_exploration", "prognostic_validation")) {
  purpose <- match.arg(purpose)
  if (!is.list(context) || !is.data.frame(context$pathways) ||
      !is.list(context$analysis)) {
    stop("`context` must contain a `pathways` data frame and an `analysis` list.",
         call. = FALSE)
  }
  missing_pathway_cols <- setdiff(c("Pathway", "Direction"),
                                  names(context$pathways))
  if (length(missing_pathway_cols)) {
    stop(sprintf("`context$pathways` is missing: %s.",
                 paste(missing_pathway_cols, collapse = ", ")), call. = FALSE)
  }
  disease <- .InterMetaPubMedScalarString(context$analysis$disease,
                                           "context$analysis$disease")
  endpoint <- .InterMetaPubMedScalarString(context$analysis$endpoint,
                                            "context$analysis$endpoint")
  if (length(max_genes) != 1L || is.na(max_genes) || !is.numeric(max_genes) ||
      !is.finite(max_genes) || max_genes < 0 ||
      max_genes > .Machine$integer.max || max_genes != floor(max_genes)) {
    stop("`max_genes` must be one non-negative integer.", call. = FALSE)
  }
  max_genes <- as.integer(max_genes)
  if (length(include_components) != 1L || is.na(include_components) ||
      !is.logical(include_components)) {
    stop("`include_components` must be TRUE or FALSE.", call. = FALSE)
  }

  pathway_table <- context$pathways[, c("Pathway", "Direction"), drop = FALSE]
  pathway_table$Pathway <- as.character(pathway_table$Pathway)
  pathway_table$Direction <- as.character(pathway_table$Direction)
  valid <- !is.na(pathway_table$Pathway) & nzchar(trimws(pathway_table$Pathway)) &
    !is.na(pathway_table$Direction) & nzchar(trimws(pathway_table$Direction))
  pathway_table <- unique(pathway_table[valid, , drop = FALSE])
  if (!nrow(pathway_table)) {
    stop("`context$pathways` contains no usable pathway rows.", call. = FALSE)
  }

  mesh_terms <- .InterMetaPubMedLookupTerms(
    disease_mesh, disease, value_names = c("MeSH", "mesh", "Mesh")
  )
  if (!length(mesh_terms)) {
    mesh_terms <- disease
  }
  disease_clause <- .InterMetaPubMedOrClause(c(
    .InterMetaPubMedTaggedPhrase(mesh_terms, "mh"),
    .InterMetaPubMedTaggedPhrase(disease, "tiab")
  ))
  if (purpose == "biology_exploration" && is.null(disease_mesh) &&
      tolower(disease) %in% c("prostate cancer", "prostate adenocarcinoma", "prad")) {
    disease_clause <- .InterMetaPubMedOrClause(c(
      '"Prostatic Neoplasms"[mh]', '"prostate cancer"[tiab]',
      '"prostatic cancer"[tiab]', '"prostate adenocarcinoma"[tiab]'))
  }
  endpoint_unknown <- tolower(endpoint) == "unknown"
  endpoint_terms <- if (endpoint_unknown) character() else endpoint
  if (!endpoint_unknown &&
      grepl("\\bbcr\\b|biochemical (recurrence|relapse)", endpoint,
            ignore.case = TRUE, perl = TRUE)) {
    endpoint_terms <- unique(c(endpoint_terms, "biochemical recurrence",
                               "biochemical relapse", "BCR"))
  }
  endpoint_clause <- .InterMetaPubMedOrClause(c(
    .InterMetaPubMedTaggedPhrase(endpoint_terms, "tiab"),
    '"Prognosis"[mh]', "prognos*[tiab]", "survival[tiab]",
    "recurrence[tiab]"
  ))
  review_clause <- .InterMetaPubMedOrClause(c(
    '"Review"[pt]', '"Systematic Review"[pt]', '"Meta-Analysis"[pt]',
    "review[tiab]"
  ))
  process_clause <- .InterMetaPubMedOrClause(c(
    "proliferat*[tiab]", "invasion[tiab]", "metasta*[tiab]",
    "resistan*[tiab]", "apoptos*[tiab]", "ferroptos*[tiab]",
    "microenvironment[tiab]", "immune[tiab]", "metabolic[tiab]",
    "progression[tiab]", "tumorigen*[tiab]"))
  exclusion <- 'NOT "Retracted Publication"[pt]'
  genes <- .InterMetaPubMedDriverGenes(
    context, unique(pathway_table$Pathway), max_genes
  )

  rows <- vector("list", nrow(pathway_table) * 7L)
  cursor <- 0L
  for (i in seq_len(nrow(pathway_table))) {
    pathway <- pathway_table$Pathway[[i]]
    direction <- pathway_table$Direction[[i]]
    built_in <- .InterMetaPubMedBuiltInTerms(pathway)
    synonyms <- c(
      pathway,
      built_in$synonyms,
      .InterMetaPubMedContextSynonyms(context, pathway, direction),
      .InterMetaPubMedLookupTerms(
        pathway_synonyms, pathway,
        value_names = c("PathwaySynonym", "pathway_synonym")
      )
    )
    synonyms <- unique(.InterMetaPubMedSafeText(
      .InterMetaPubMedSplitSynonyms(synonyms)
    ))
    synonyms <- synonyms[nzchar(synonyms)]
    components <- if (include_components) {
      unique(.InterMetaPubMedSafeText(.InterMetaPubMedSplitSynonyms(c(
        built_in$components,
        .InterMetaPubMedLookupTerms(
          pathway_components, pathway,
          value_names = c("PathwayComponent", "pathway_component", "Component", "component")
        )
      ))))
    } else character()
    components <- components[nzchar(components)]
    pathway_clause <- .InterMetaPubMedOrClause(
      .InterMetaPubMedTaggedPhrase(unique(synonyms), "tiab")
    )
    base <- paste(disease_clause, "AND", pathway_clause)
    layer_queries <- c(
      prognosis = paste(base, "AND", endpoint_clause, exclusion),
      biology = paste(base, exclusion),
      review = paste(base, "AND", review_clause, exclusion)
    )
    if (purpose == "biology_exploration") {
      layer_queries <- c(biology = paste(base, exclusion),
        biology_process = paste(base, "AND", process_clause, exclusion),
        review = paste(base, "AND", review_clause, exclusion))
    }
    driver_genes <- genes[[pathway]]
    if (length(driver_genes)) {
      gene_clause <- .InterMetaPubMedOrClause(
        .InterMetaPubMedTaggedPhrase(driver_genes, "tiab")
      )
      layer_queries <- c(
        layer_queries,
        gene_bridge = paste(base, "AND", gene_clause, exclusion)
      )
    }
    if (length(components)) {
      component_clause <- .InterMetaPubMedOrClause(
        .InterMetaPubMedTaggedPhrase(components, "tiab")
      )
      component_base <- paste(disease_clause, "AND", component_clause)
      component_queries <- if (purpose == "biology_exploration") c(
        component_biology = paste(component_base, exclusion),
        component_process = paste(component_base, "AND", process_clause, exclusion),
        component_review = paste(component_base, "AND", review_clause, exclusion)) else c(
        component_prognosis = paste(component_base, "AND", endpoint_clause, exclusion),
        component_biology = paste(component_base, exclusion),
        component_review = paste(component_base, "AND", review_clause, exclusion))
      layer_queries <- c(layer_queries, component_queries)
    }

    for (layer in names(layer_queries)) {
      cursor <- cursor + 1L
      is_component <- startsWith(layer, "component_")
      query_warning <- if (is_component) {
        paste("Component terms are broader retrieval candidates, not synonyms",
              "for the complete pathway; component hits are not whole-pathway evidence.")
      } else if (length(synonyms) > 1L) {
        paste("Pathway search expansions are retrieval aids, not proof of",
              "strict pathway equivalence.")
      } else ""
      if (endpoint_unknown && purpose == "prognostic_validation") {
        query_warning <- trimws(paste(query_warning,
          "Endpoint 'unknown' was omitted; generic prognosis terms were retained."))
      }
      rows[[cursor]] <- data.frame(
        query_id = sprintf("q%03d_%s", cursor, layer),
        pathway = pathway,
        direction = direction,
        layer = layer,
        term_scope = if (is_component) "component" else "pathway",
        search_terms = paste(if (is_component) components else synonyms, collapse = ";"),
        query_version = if (purpose == "biology_exploration") "3.1" else "2.1",
        query_warning = query_warning,
        disease = disease,
        endpoint = endpoint,
        driver_genes = paste(driver_genes, collapse = ";"),
        query = unname(layer_queries[[layer]]),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows[seq_len(cursor)])
  if (purpose == "biology_exploration") out$purpose <- purpose
  rownames(out) <- NULL
  out
}

.InterMetaPubMedRequire <- function(request_fun) {
  packages <- c("jsonlite", "xml2")
  if (is.null(request_fun)) {
    packages <- c("httr2", packages)
  }
  missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing)) {
    stop(
      sprintf(
        "InterMetaPubMedSearch() requires the optional package%s %s. Install %s before searching PubMed.",
        if (length(missing) == 1L) "" else "s",
        paste(sprintf("'%s'", missing), collapse = ", "),
        paste(sprintf("install.packages('%s')", missing), collapse = "; ")
      ),
      call. = FALSE
    )
  }
}

.InterMetaPubMedDefaultRequestObject <- function(endpoint, params) {
  request <- httr2::request(endpoint)
  request <- do.call(httr2::req_url_query, c(list(request), params))
  request <- httr2::req_user_agent(request, "InterMeta PubMed evidence retrieval")
  request <- httr2::req_timeout(request, seconds = 30)
  request <- httr2::req_retry(
    request, max_tries = 3L, max_seconds = 60,
    is_transient = function(response) {
      status <- httr2::resp_status(response)
      status == 429L || (status >= 500L && status < 600L)
    },
    backoff = function(tries) min(8, 2^(tries - 1L)),
    after = function(response) {
      seconds <- suppressWarnings(as.numeric(
        httr2::resp_header(response, "retry-after", default = NA_character_)
      ))
      if (length(seconds) != 1L || is.na(seconds) || !is.finite(seconds)) {
        return(NA_real_)
      }
      min(10, max(0.34, seconds))
    }
  )
  request <- httr2::req_error(request, is_error = function(resp) FALSE)
  request
}

.InterMetaPubMedDefaultRequest <- function(endpoint, params) {
  request <- .InterMetaPubMedDefaultRequestObject(endpoint, params)
  response <- httr2::req_perform(request)
  list(
    status_code = httr2::resp_status(response),
    body = httr2::resp_body_raw(response)
  )
}

.InterMetaPubMedNormalizeResponse <- function(response) {
  if (inherits(response, "httr2_response")) {
    return(list(
      status_code = httr2::resp_status(response),
      body = rawToChar(httr2::resp_body_raw(response))
    ))
  }
  if (is.raw(response)) {
    return(list(status_code = 200L, body = rawToChar(response)))
  }
  if (is.character(response) && length(response) == 1L) {
    return(list(status_code = 200L, body = response))
  }
  if (!is.list(response)) {
    stop("`request_fun` must return text, raw bytes, an httr2 response, or a response list.",
         call. = FALSE)
  }
  status <- response$status_code
  if (is.null(status)) status <- response$status
  if (is.null(status)) status <- 200L
  body <- response$body
  if (is.null(body)) body <- response$content
  if (is.null(body)) body <- response$text
  if (is.null(body)) body <- ""
  if (is.raw(body)) body <- rawToChar(body)
  if (!is.character(body) || length(body) != 1L) {
    stop("The response body returned by `request_fun` must be one character string or raw vector.",
         call. = FALSE)
  }
  list(status_code = as.integer(status[[1L]]), body = enc2utf8(body))
}

.InterMetaPubMedRequest <- function(request_fun, endpoint, params, delay,
                                     sleep_fun, secrets = character()) {
  sleep_fun(delay)
  response <- do.call(request_fun, list(endpoint = endpoint, params = params))
  response <- .InterMetaPubMedNormalizeResponse(response)
  if (is.na(response$status_code) || response$status_code < 200L ||
      response$status_code >= 300L) {
    stop(sprintf("NCBI returned HTTP status %s.", response$status_code),
         call. = FALSE)
  }
  original_body <- response$body
  response$body <- .InterMetaPubMedRedact(response$body, secrets)
  response$body_redacted <- !identical(original_body, response$body)
  response
}

.InterMetaPubMedParseESearch <- function(text) {
  parsed <- tryCatch(
    jsonlite::fromJSON(text, simplifyVector = TRUE),
    error = function(e) stop(sprintf("Invalid ESearch JSON: %s", conditionMessage(e)),
                             call. = FALSE)
  )
  result <- parsed$esearchresult
  if (is.null(result)) {
    stop("ESearch JSON does not contain `esearchresult`.", call. = FALSE)
  }
  if (!is.null(result$error) && length(result$error) && nzchar(result$error[[1L]])) {
    stop(sprintf("NCBI ESearch error: %s", result$error[[1L]]), call. = FALSE)
  }
  count <- suppressWarnings(as.numeric(result$count[[1L]]))
  if (length(count) != 1L || is.na(count) || !is.finite(count) || count < 0 ||
      count != floor(count)) {
    stop("ESearch returned an invalid result count.", call. = FALSE)
  }
  ids <- as.character(unlist(result$idlist, use.names = FALSE))
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (any(!grepl("^[0-9]+$", ids)) || length(ids) > count) {
    stop("ESearch returned invalid PubMed identifiers.", call. = FALSE)
  }
  list(count = count, ids = ids)
}

.InterMetaPubMedNodeText <- function(node, paths) {
  for (path in paths) {
    found <- xml2::xml_find_first(node, path)
    if (!inherits(found, "xml_missing")) {
      value <- trimws(xml2::xml_text(found))
      if (nzchar(value)) return(value)
    }
  }
  NA_character_
}

.InterMetaPubMedManyText <- function(node, path) {
  found <- xml2::xml_find_all(node, path)
  if (!length(found)) return(character())
  values <- trimws(xml2::xml_text(found))
  unique(values[nzchar(values)])
}

.InterMetaPubMedParseEFetch <- function(text) {
  doc <- tryCatch(
    xml2::read_xml(text),
    error = function(e) stop(sprintf("Invalid EFetch XML: %s", conditionMessage(e)),
                             call. = FALSE)
  )
  api_error <- xml2::xml_find_first(doc, ".//ERROR | .//Error")
  if (!inherits(api_error, "xml_missing")) {
    stop(sprintf("NCBI EFetch error: %s", trimws(xml2::xml_text(api_error))),
         call. = FALSE)
  }
  articles <- xml2::xml_find_all(doc, ".//PubmedArticle | .//PubmedBookArticle")
  if (!length(articles)) {
    return(.InterMetaPubMedEmptyCoreRecords())
  }

  rows <- lapply(articles, function(article) {
    abstract_nodes <- xml2::xml_find_all(
      article,
      ".//Article/Abstract/AbstractText | .//BookDocument/Abstract/AbstractText"
    )
    abstract <- character()
    if (length(abstract_nodes)) {
      abstract <- vapply(abstract_nodes, function(node) {
        value <- trimws(xml2::xml_text(node))
        label <- xml2::xml_attr(node, "Label")
        if (!is.na(label) && nzchar(label)) paste0(label, ": ", value) else value
      }, character(1L))
    }
    year <- .InterMetaPubMedNodeText(article, c(
      ".//Article/Journal/JournalIssue/PubDate/Year",
      ".//PubmedData/History/PubMedPubDate[@PubStatus='pubmed']/Year",
      ".//BookDocument/Book/PubDate/Year"
    ))
    if (is.na(year)) {
      medline_date <- .InterMetaPubMedNodeText(
        article, c(".//Article/Journal/JournalIssue/PubDate/MedlineDate")
      )
      if (!is.na(medline_date)) {
        matched <- regmatches(medline_date, regexpr("[12][0-9]{3}", medline_date))
        if (length(matched) && nzchar(matched)) year <- matched
      }
    }
    doi <- .InterMetaPubMedNodeText(article, c(
      ".//PubmedData/ArticleIdList/ArticleId[@IdType='doi']",
      ".//Article/ELocationID[@EIdType='doi']"
    ))
    data.frame(
      PMID = .InterMetaPubMedNodeText(article, c(
        ".//MedlineCitation/PMID", ".//BookDocument/PMID"
      )),
      Title = .InterMetaPubMedNodeText(article, c(
        ".//Article/ArticleTitle", ".//BookDocument/ArticleTitle"
      )),
      Abstract = if (length(abstract)) paste(abstract, collapse = " ") else NA_character_,
      Journal = .InterMetaPubMedNodeText(article, c(
        ".//Article/Journal/Title", ".//MedlineJournalInfo/MedlineTA",
        ".//BookDocument/Book/BookTitle"
      )),
      Year = year,
      DOI = doi,
      PublicationType = paste(.InterMetaPubMedManyText(
        article, ".//PublicationTypeList/PublicationType"
      ), collapse = ";"),
      MeSH = paste(.InterMetaPubMedManyText(
        article, ".//MeshHeadingList/MeshHeading/DescriptorName"
      ), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[!is.na(out$PMID) & nzchar(out$PMID), , drop = FALSE]
  out <- out[!duplicated(out$PMID), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.InterMetaPubMedEmptyCoreRecords <- function() {
  data.frame(
    PMID = character(), Title = character(), Abstract = character(),
    Journal = character(), Year = character(), DOI = character(),
    PublicationType = character(), MeSH = character(),
    stringsAsFactors = FALSE
  )
}

.InterMetaPubMedEmptyRecords <- function() {
  data.frame(
    query_id = character(), pathway = character(), direction = character(),
    layer = character(), term_scope = character(), search_terms = character(),
    query_version = character(), query_warning = character(),
    query = character(), PMID = character(),
    Title = character(), Abstract = character(), Journal = character(),
    Year = character(), DOI = character(), PublicationType = character(),
    MeSH = character(), stringsAsFactors = FALSE
  )
}

.InterMetaPubMedHashText <- function(text) {
  digest::digest(charToRaw(enc2utf8(text)), algo = "sha256", serialize = FALSE)
}

.InterMetaPubMedCacheKey <- function(query, retmax,
                                     request_source = "ncbi_eutilities") {
  digest::digest(
    list(cache_version = 2L, database = "pubmed", query = query,
         retmax = as.integer(retmax), request_source = request_source),
    algo = "sha256", serialize = TRUE
  )
}

.InterMetaPubMedHashIDs <- function(ids) {
  .InterMetaPubMedHashText(as.character(jsonlite::toJSON(
    unname(enc2utf8(ids)), auto_unbox = FALSE, na = "null"
  )))
}

.InterMetaPubMedHashRecords <- function(records) {
  columns <- names(.InterMetaPubMedEmptyCoreRecords())
  canonical <- lapply(records[, columns, drop = FALSE], function(x) {
    unname(enc2utf8(as.character(x)))
  })
  .InterMetaPubMedHashText(as.character(jsonlite::toJSON(
    canonical, auto_unbox = FALSE, na = "null"
  )))
}

.InterMetaPubMedValidCache <- function(x, key, query, retmax,
                                       request_source = "ncbi_eutilities") {
  if (!is.list(x) || !identical(x$cache_version, 2L) ||
      !identical(x$key, key) || !identical(x$query, query) ||
      !identical(x$retmax, as.integer(retmax)) ||
      !identical(x$request_source, request_source) ||
      !is.list(x$raw_responses) || !length(x$raw_responses) ||
      !is.data.frame(x$records) ||
      !identical(names(x$records), names(.InterMetaPubMedEmptyCoreRecords())) ||
      !all(vapply(x$records, is.character, logical(1L))) ||
      !is.character(x$ids) || anyNA(x$ids) || anyDuplicated(x$ids) ||
      any(!grepl("^[0-9]+$", x$ids)) || length(x$ids) > retmax ||
      anyNA(x$records$PMID) || anyDuplicated(x$records$PMID) ||
      any(!x$records$PMID %in% x$ids) ||
      length(x$total_count) != 1L || !is.numeric(x$total_count) ||
      is.na(x$total_count) || !is.finite(x$total_count) ||
      x$total_count < length(x$ids) ||
      !is.character(x$searched_at) || length(x$searched_at) != 1L ||
      is.na(x$searched_at) || !nzchar(x$searched_at) ||
      !is.character(x$error) || length(x$error) != 1L ||
      !identical(x$ids_sha256, .InterMetaPubMedHashIDs(x$ids)) ||
      !identical(x$records_sha256, .InterMetaPubMedHashRecords(x$records))) {
    return(FALSE)
  }
  valid_raw <- all(vapply(x$raw_responses, function(response) {
    is.list(response) && is.character(response$body) &&
      length(response$body) == 1L && !is.na(response$body) &&
      is.character(response$sha256) && length(response$sha256) == 1L &&
      is.character(response$stage) && length(response$stage) == 1L &&
      !is.na(response$stage) &&
      response$stage %in% c("esearch", "efetch") &&
      is.numeric(response$batch) && length(response$batch) == 1L &&
      !is.na(response$batch) && is.finite(response$batch) &&
      response$batch >= 1 && response$batch == floor(response$batch) &&
      is.numeric(response$status_code) && length(response$status_code) == 1L &&
      !is.na(response$status_code) && response$status_code >= 200L &&
      response$status_code < 300L &&
      is.logical(response$body_redacted) && length(response$body_redacted) == 1L &&
      !is.na(response$body_redacted) &&
      identical(.InterMetaPubMedHashText(response$body), response$sha256)
  }, logical(1L)))
  if (!valid_raw) return(FALSE)
  # Raw responses remain the source of truth. Re-parsing prevents a modified
  # canonical table (even with an updated table hash) from masquerading as the
  # actual E-utilities response.
  tryCatch({
    stages <- vapply(x$raw_responses, `[[`, character(1L), "stage")
    if (!identical(stages[[1L]], "esearch") || sum(stages == "esearch") != 1L) {
      return(FALSE)
    }
    parsed_search <- .InterMetaPubMedParseESearch(x$raw_responses[[1L]]$body)
    if (!identical(x$total_count, parsed_search$count) ||
        !identical(x$ids, utils::head(parsed_search$ids, retmax))) {
      return(FALSE)
    }
    fetched <- lapply(x$raw_responses[stages == "efetch"], function(response) {
      .InterMetaPubMedParseEFetch(response$body)
    })
    reconstructed <- if (length(fetched)) do.call(rbind, fetched) else
      .InterMetaPubMedEmptyCoreRecords()
    reconstructed <- reconstructed[reconstructed$PMID %in% x$ids, , drop = FALSE]
    reconstructed <- reconstructed[!duplicated(reconstructed$PMID), , drop = FALSE]
    expected_status <- if (nrow(reconstructed) < length(x$ids)) "partial" else "ok"
    identical(x$api_status, expected_status) &&
      identical(.InterMetaPubMedHashRecords(reconstructed), x$records_sha256)
  }, error = function(e) FALSE)
}

.InterMetaPubMedWriteCache <- function(x, path) {
  tmp <- tempfile("intermeta_pubmed_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, version = 3L)
  if (!file.copy(tmp, path, overwrite = TRUE)) {
    warning(sprintf("Could not write PubMed cache file '%s'.", path),
            call. = FALSE)
  }
  invisible(path)
}

.InterMetaPubMedNow <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

.InterMetaPubMedSearchRow <- function(query_row, searched_at, retmax, key,
                                      cached, total_count = NA_real_,
                                      returned_pmids = NA_integer_,
                                      records_fetched = NA_integer_,
                                      truncated = NA,
                                      api_status = "error", error = NA_character_) {
  out <- data.frame(
    query_id = query_row$query_id,
    pathway = query_row$pathway,
    direction = query_row$direction,
    layer = query_row$layer,
    term_scope = query_row$term_scope,
    search_terms = query_row$search_terms,
    query_version = query_row$query_version,
    query_warning = query_row$query_warning,
    query = query_row$query,
    searched_at = searched_at,
    retmax = as.integer(retmax),
    total_count = as.numeric(total_count),
    returned_pmids = as.integer(returned_pmids),
    records_fetched = as.integer(records_fetched),
    truncated = as.logical(truncated),
    api_status = api_status,
    error = error,
    cache_key = key,
    cached = cached,
    stringsAsFactors = FALSE
  )
  if ("purpose" %in% names(query_row)) out$purpose <- query_row$purpose
  out
}

.InterMetaPubMedAttachQuery <- function(records, query_row) {
  if (!nrow(records)) {
    out <- .InterMetaPubMedEmptyRecords()
    if ("purpose" %in% names(query_row)) out$purpose <- character()
    return(out)
  }
  out <- data.frame(
    query_id = rep(query_row$query_id, nrow(records)),
    pathway = rep(query_row$pathway, nrow(records)),
    direction = rep(query_row$direction, nrow(records)),
    layer = rep(query_row$layer, nrow(records)),
    term_scope = rep(query_row$term_scope, nrow(records)),
    search_terms = rep(query_row$search_terms, nrow(records)),
    query_version = rep(query_row$query_version, nrow(records)),
    query_warning = rep(query_row$query_warning, nrow(records)),
    query = rep(query_row$query, nrow(records)),
    records,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if ("purpose" %in% names(query_row)) out$purpose <- rep(query_row$purpose, nrow(records))
  out
}

#' Retrieve traceable PubMed evidence for an InterMeta AI context
#'
#' @param context A list produced by `InterMetaAIContext()`.
#' @param queries Optional query table produced by
#'   `InterMetaBuildPubMedQueries()`. Query-scope metadata is retained in search
#'   and article tables. Legacy tables without this metadata remain accepted;
#'   their version is marked `"legacy"` and unavailable terms remain missing.
#' @param retmax Maximum PubMed identifiers returned per query.
#' @param email Optional contact email sent to NCBI.
#' @param tool Tool name sent to NCBI.
#' @param api_key NCBI API key. The key is never stored in the result or cache.
#' @param cache_dir Optional directory for query-level raw-response caches.
#' @param refresh Ignore valid cache entries and retrieve fresh records.
#' @param request_fun Optional function with arguments `endpoint` and `params`.
#'   It is intended for controlled testing; by default `httr2` performs explicit
#'   network requests to NCBI E-utilities.
#' @param sleep_fun Function used for NCBI rate limiting.
#' @param page_size Optional ESearch page size. Enables fingerprinted paged
#'   retrieval; retmax remains the total per-query limit (at most 10,000).
#' @param sort Optional PubMed sort, relevance or pub_date. An explicit sort
#'   enables the paged retriever and is included in its cache identity.
#' @param purpose Query goal when queries are built automatically. If omitted,
#'   new queries use biology_exploration; supplied query metadata is preserved.
#'
#' @return An object of class `InterMetaPubMedEvidence` containing query,
#'   retrieval-status, article, response-hash, and provenance tables.
#' @export
InterMetaPubMedSearch <- function(context, queries = NULL, retmax = 20L,
                                   email = NULL, tool = "InterMeta",
                                   api_key = Sys.getenv("NCBI_API_KEY"),
                                   cache_dir = NULL, refresh = FALSE,
                                   request_fun = NULL,
                                   sleep_fun = Sys.sleep,
                                   page_size = NULL, sort = NULL, purpose = NULL) {
  .InterMetaPubMedRequire(request_fun)
  if (is.null(queries)) {
    queries <- InterMetaBuildPubMedQueries(context,
      purpose = if (is.null(purpose)) "biology_exploration" else purpose)
  } else if (!is.null(purpose)) {
    purpose <- match.arg(purpose, c("biology_exploration", "prognostic_validation"))
    known_purpose <- if ("purpose" %in% names(queries)) unique(queries$purpose) else "prognostic_validation"
    if (!identical(as.character(known_purpose), purpose))
      stop("Supplied query purpose does not match the requested purpose.", call. = FALSE)
  }
  required_query_cols <- c("query_id", "pathway", "direction", "layer", "query")
  if (!is.data.frame(queries) ||
      length(setdiff(required_query_cols, names(queries)))) {
    stop(sprintf("`queries` must contain columns: %s.",
                 paste(required_query_cols, collapse = ", ")), call. = FALSE)
  }
  if (!nrow(queries)) stop("`queries` must contain at least one query.", call. = FALSE)
  if (anyDuplicated(queries$query_id)) {
    stop("`queries$query_id` values must be unique.", call. = FALSE)
  }
  if (anyNA(queries$query) || any(!nzchar(trimws(queries$query)))) {
    stop("Every PubMed query must be a non-empty string.", call. = FALSE)
  }
  if ("purpose" %in% names(queries) && (anyNA(queries$purpose) ||
      length(unique(queries$purpose)) != 1L ||
      !all(queries$purpose %in% c("biology_exploration", "prognostic_validation")))) {
    stop("Query purpose must be one consistent supported goal.", call. = FALSE)
  }
  # Accept pre-v2 user-supplied query tables without inventing search terms or
  # labelling their query text as generated by the current dictionary.
  if (!"term_scope" %in% names(queries)) {
    queries$term_scope <- ifelse(startsWith(as.character(queries$layer), "component_"),
                                 "component", "pathway")
  }
  queries$term_scope <- as.character(queries$term_scope)
  if (anyNA(queries$term_scope) ||
      any(!queries$term_scope %in% c("pathway", "component"))) {
    stop("`queries$term_scope` must contain only 'pathway' or 'component'.", call. = FALSE)
  }
  if (!"search_terms" %in% names(queries)) queries$search_terms <- NA_character_
  if (!"query_version" %in% names(queries)) queries$query_version <- "legacy"
  if (!"query_warning" %in% names(queries)) {
    queries$query_warning <- paste("Legacy query metadata is incomplete;",
                                    "retrieval terms do not establish pathway equivalence.")
  }
  for (column in c("search_terms", "query_version", "query_warning")) {
    queries[[column]] <- as.character(queries[[column]])
  }
  if (length(retmax) != 1L || is.na(retmax) || !is.numeric(retmax) ||
      !is.finite(retmax) || retmax < 1 || retmax > .Machine$integer.max ||
      retmax != floor(retmax)) {
    stop("`retmax` must be one positive integer.", call. = FALSE)
  }
  retmax <- as.integer(retmax)
  tool <- .InterMetaPubMedScalarString(tool, "tool")
  if (!is.null(email)) email <- .InterMetaPubMedScalarString(email, "email")
  if (length(api_key) > 1L || anyNA(api_key)) {
    stop("`api_key` must be one character string.", call. = FALSE)
  }
  api_key <- if (length(api_key)) trimws(as.character(api_key)) else ""
  if (length(refresh) != 1L || is.na(refresh) || !is.logical(refresh)) {
    stop("`refresh` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.function(sleep_fun)) stop("`sleep_fun` must be a function.", call. = FALSE)
  request_source <- if (is.null(request_fun)) "ncbi_eutilities" else
    "custom_request_fun"
  if (is.null(request_fun)) request_fun <- .InterMetaPubMedDefaultRequest
  if (!is.function(request_fun)) stop("`request_fun` must be a function.", call. = FALSE)

  if (!is.null(cache_dir)) {
    cache_dir <- .InterMetaPubMedScalarString(cache_dir, "cache_dir")
    if (!dir.exists(cache_dir) &&
        !dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop(sprintf("Could not create PubMed cache directory '%s'.", cache_dir),
           call. = FALSE)
    }
    cache_dir <- normalizePath(cache_dir, winslash = "/", mustWork = TRUE)
  }

  delay <- if (nzchar(api_key)) 0.11 else 0.34
  if (!is.null(page_size) || !is.null(sort)) {
    return(.InterMetaPubMedPagedSearch(queries, retmax, page_size, sort,
      email, tool, api_key, cache_dir, refresh, request_fun, request_source,
      sleep_fun, delay))
  }
  esearch_endpoint <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  )
  efetch_endpoint <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  )
  search_rows <- vector("list", nrow(queries))
  record_rows <- vector("list", nrow(queries))
  hash_rows <- list()
  hash_cursor <- 0L

  for (i in seq_len(nrow(queries))) {
    query_row <- queries[i, , drop = FALSE]
    key <- .InterMetaPubMedCacheKey(query_row$query, retmax, request_source)
    cache_path <- if (is.null(cache_dir)) NULL else
      file.path(cache_dir, paste0("pubmed_", key, ".rds"))
    cached_value <- NULL
    if (!refresh && !is.null(cache_path) && file.exists(cache_path)) {
      cached_value <- tryCatch(readRDS(cache_path), error = function(e) NULL)
      if (!.InterMetaPubMedValidCache(cached_value, key, query_row$query, retmax,
                                    request_source)) {
        cached_value <- NULL
        warning(sprintf("Ignoring an invalid PubMed cache entry for query '%s'.",
                        query_row$query_id), call. = FALSE)
      }
    }

    if (!is.null(cached_value)) {
      search_rows[[i]] <- .InterMetaPubMedSearchRow(
        query_row = query_row,
        searched_at = cached_value$searched_at,
        retmax = retmax,
        key = key,
        cached = TRUE,
        total_count = cached_value$total_count,
        returned_pmids = length(cached_value$ids),
        records_fetched = nrow(cached_value$records),
        truncated = cached_value$total_count > length(cached_value$ids),
        api_status = cached_value$api_status,
        error = cached_value$error
      )
      record_rows[[i]] <- .InterMetaPubMedAttachQuery(
        cached_value$records, query_row
      )
      for (j in seq_along(cached_value$raw_responses)) {
        raw_response <- cached_value$raw_responses[[j]]
        hash_cursor <- hash_cursor + 1L
        hash_rows[[hash_cursor]] <- data.frame(
          query_id = query_row$query_id,
          stage = raw_response$stage,
          batch = raw_response$batch,
          status_code = raw_response$status_code,
          sha256 = raw_response$sha256,
          body_redacted = raw_response$body_redacted,
          cached = TRUE,
          stringsAsFactors = FALSE
        )
      }
      next
    }

    searched_at <- .InterMetaPubMedNow()
    raw_responses <- list()
    status <- "ok"
    error_message <- NA_character_
    total_count <- NA_real_
    ids <- character()
    core_records <- .InterMetaPubMedEmptyCoreRecords()

    esearch_params <- list(
      db = "pubmed", term = query_row$query, retmode = "json",
      retmax = retmax, retstart = 0L, tool = tool
    )
    if (!is.null(email)) esearch_params$email <- email
    if (nzchar(api_key)) esearch_params$api_key <- api_key
    esearch <- tryCatch(
      .InterMetaPubMedRequest(
        request_fun, esearch_endpoint, esearch_params, delay, sleep_fun, api_key
      ),
      error = function(e) e
    )
    if (inherits(esearch, "error")) {
      status <- "error"
      error_message <- .InterMetaPubMedRedact(
        paste("ESearch failed:", conditionMessage(esearch)), api_key
      )
    } else {
      raw_responses[[1L]] <- list(
        stage = "esearch", batch = 1L, status_code = esearch$status_code,
        body = esearch$body,
        body_redacted = esearch$body_redacted,
        sha256 = .InterMetaPubMedHashText(esearch$body)
      )
      parsed_search <- tryCatch(
        .InterMetaPubMedParseESearch(esearch$body), error = function(e) e
      )
      if (inherits(parsed_search, "error")) {
        status <- "error"
        error_message <- .InterMetaPubMedRedact(
          paste("ESearch parsing failed:", conditionMessage(parsed_search)),
          api_key
        )
      } else {
        total_count <- parsed_search$count
        ids <- utils::head(parsed_search$ids, retmax)
      }
    }

    if (identical(status, "ok") && length(ids)) {
      batches <- split(ids, ceiling(seq_along(ids) / 200L))
      fetched <- vector("list", length(batches))
      fetch_errors <- character()
      for (batch_index in seq_along(batches)) {
        efetch_params <- list(
          db = "pubmed", id = paste(batches[[batch_index]], collapse = ","),
          retmode = "xml", tool = tool
        )
        if (!is.null(email)) efetch_params$email <- email
        if (nzchar(api_key)) efetch_params$api_key <- api_key
        fetched_response <- tryCatch(
          .InterMetaPubMedRequest(
            request_fun, efetch_endpoint, efetch_params, delay, sleep_fun, api_key
          ),
          error = function(e) e
        )
        if (inherits(fetched_response, "error")) {
          fetch_errors <- c(
            fetch_errors,
            .InterMetaPubMedRedact(
              sprintf("EFetch batch %d failed: %s", batch_index,
                      conditionMessage(fetched_response)),
              api_key
            )
          )
          next
        }
        raw_index <- length(raw_responses) + 1L
        raw_responses[[raw_index]] <- list(
          stage = "efetch", batch = as.integer(batch_index),
          status_code = fetched_response$status_code,
          body = fetched_response$body,
          body_redacted = fetched_response$body_redacted,
          sha256 = .InterMetaPubMedHashText(fetched_response$body)
        )
        parsed_fetch <- tryCatch(
          .InterMetaPubMedParseEFetch(fetched_response$body),
          error = function(e) e
        )
        if (inherits(parsed_fetch, "error")) {
          fetch_errors <- c(
            fetch_errors,
            .InterMetaPubMedRedact(
              sprintf("EFetch batch %d parsing failed: %s", batch_index,
                      conditionMessage(parsed_fetch)), api_key
            )
          )
        } else {
          fetched[[batch_index]] <- parsed_fetch[
            parsed_fetch$PMID %in% batches[[batch_index]], , drop = FALSE
          ]
        }
      }
      fetched <- fetched[vapply(fetched, is.data.frame, logical(1L))]
      if (length(fetched)) {
        core_records <- do.call(rbind, fetched)
        core_records <- core_records[!duplicated(core_records$PMID), , drop = FALSE]
        rownames(core_records) <- NULL
      }
      if (length(fetch_errors)) {
        status <- if (nrow(core_records)) "partial_error" else "error"
        error_message <- paste(fetch_errors, collapse = " | ")
      } else if (nrow(core_records) < length(ids)) {
        status <- "partial"
        error_message <- sprintf(
          "EFetch returned metadata for %d of %d requested PMID(s).",
          nrow(core_records), length(ids)
        )
      }
    }

    search_rows[[i]] <- .InterMetaPubMedSearchRow(
      query_row = query_row,
      searched_at = searched_at,
      retmax = retmax,
      key = key,
      cached = FALSE,
      total_count = total_count,
      returned_pmids = if (is.na(total_count)) NA_integer_ else length(ids),
      records_fetched = if (is.na(total_count)) NA_integer_ else nrow(core_records),
      truncated = if (is.na(total_count)) NA else total_count > length(ids),
      api_status = status,
      error = error_message
    )
    record_rows[[i]] <- .InterMetaPubMedAttachQuery(core_records, query_row)

    for (j in seq_along(raw_responses)) {
      raw_response <- raw_responses[[j]]
      hash_cursor <- hash_cursor + 1L
      hash_rows[[hash_cursor]] <- data.frame(
        query_id = query_row$query_id,
        stage = raw_response$stage,
        batch = raw_response$batch,
        status_code = raw_response$status_code,
        sha256 = raw_response$sha256,
        body_redacted = raw_response$body_redacted,
        cached = FALSE,
        stringsAsFactors = FALSE
      )
    }

    if (!is.null(cache_path) && status %in% c("ok", "partial")) {
      cache_value <- list(
        cache_version = 2L,
        key = key,
        query = query_row$query,
        retmax = retmax,
        request_source = request_source,
        searched_at = searched_at,
        total_count = total_count,
        ids = ids,
        ids_sha256 = .InterMetaPubMedHashIDs(ids),
        records = core_records,
        records_sha256 = .InterMetaPubMedHashRecords(core_records),
        api_status = status,
        error = error_message,
        raw_responses = raw_responses
      )
      .InterMetaPubMedWriteCache(cache_value, cache_path)
    }
  }

  searches <- do.call(rbind, search_rows)
  rownames(searches) <- NULL
  nonempty_records <- record_rows[vapply(record_rows, nrow, integer(1L)) > 0L]
  records <- if (length(nonempty_records)) {
    out <- do.call(rbind, nonempty_records)
    rownames(out) <- NULL
    out
  } else {
    .InterMetaPubMedEmptyRecords()
  }
  response_hashes <- if (length(hash_rows)) {
    out <- do.call(rbind, hash_rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame(
      query_id = character(), stage = character(), batch = integer(),
      status_code = integer(), sha256 = character(), body_redacted = logical(),
      cached = logical(),
      stringsAsFactors = FALSE
    )
  }
  failure_statuses <- c("error", "partial_error")
  failed <- searches$api_status %in% failure_statuses
  if (any(failed)) {
    warning(
      sprintf(
        "PubMed retrieval failed for %d of %d queries; failed searches retain NA/error status and must not be interpreted as zero evidence.",
        sum(failed), nrow(searches)
      ),
      call. = FALSE
    )
  }

  structure(
    list(
      schema_version = "1.0",
      generated_at = .InterMetaPubMedNow(),
      search_date = format(Sys.Date(), "%Y-%m-%d"),
      queries = queries,
      searches = searches,
      records = records,
      response_hashes = response_hashes,
      provenance = list(
        provider = "NCBI E-utilities",
        request_source = request_source,
        database = "PubMed",
        tool = tool,
        email_supplied = !is.null(email),
        api_key_used = nzchar(api_key),
        cache_enabled = !is.null(cache_dir),
        cache_path_sha256 = if (is.null(cache_dir)) NULL else
          .InterMetaPubMedHashText(cache_dir),
        response_hash_scope = "stored_body_after_secret_redaction",
        refresh = refresh,
        rate_limit_seconds = delay,
        retmax = retmax,
        queries_truncated = sum(searches$truncated %in% TRUE, na.rm = TRUE),
        truncation_present = any(searches$truncated %in% TRUE, na.rm = TRUE)
      )
    ),
    class = c("InterMetaPubMedEvidence", "list")
  )
}
