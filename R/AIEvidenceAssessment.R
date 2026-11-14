# Endpoint-aware evidence grading. Model classifications remain an audit trail;
# final grades use only eligible, source-grounded per-publication assessments.

.InterMetaAIAssessmentEnums <- function() {
  list(
    evidence_type = c("human_prognosis", "human_incidence", "human_tissue",
                      "preclinical", "review", "other", "unknown"),
    endpoint_match = c("matched", "mismatched", "unknown"),
    direction = c("same_direction", "opposite_direction", "no_association",
                  "not_assessed"),
    directness = c("direct", "indirect")
  )
}

.InterMetaAIAssessmentSchema <- function() {
  enums <- .InterMetaAIAssessmentEnums()
  list(
    type = "array",
    description = "Exactly one assessment per distinct cited PMID for this pathway.",
    items = list(
      type = "object", additionalProperties = FALSE,
      required = c("pmid", "evidence_type", "study_endpoint", "endpoint_match",
                   "direction", "directness", "abstract_quote"),
      properties = list(
        pmid = list(type = "string"),
        evidence_type = list(type = "string", enum = enums$evidence_type),
        study_endpoint = list(type = "string", minLength = 1L,
                              description = "Actual study endpoint; use unknown if unavailable."),
        endpoint_match = list(type = "string", enum = enums$endpoint_match),
        direction = list(type = "string", enum = enums$direction),
        directness = list(type = "string", enum = enums$directness),
        abstract_quote = list(
          type = "string", maxLength = 1000L,
          description = paste(
            "Short exact excerpt from the supplied Abstract, not Title.",
            "Use an empty string when no relevant excerpt is available;",
            "such evidence cannot support a directional grade."
          )
        )
      )
    )
  )
}

.InterMetaAIAssessmentPrompt <- function() {
  paste(
    "For every distinct cited PMID in each pathway, supply one assessment with",
    "evidence_type, actual study_endpoint, endpoint_match, direction, directness,",
    "and a short exact abstract_quote from its supplied Abstract (not its Title).",
    "A missing abstract or relevant passage requires an empty quote and cannot",
    "support a directional grade. Compare the actual study endpoint with",
    "context.analysis.endpoint; unknown target endpoints cannot be matched.",
    "BCR means biochemical recurrence or biochemical relapse, not cancer risk,",
    "incidence, overall survival, or generic prognosis. A matched BCR study_endpoint",
    "must explicitly name BCR, biochemical recurrence, or biochemical relapse.",
    "Human incidence or cancer risk is not human prognosis. No association,",
    "a null result, or lack of significance is never opposite_direction.",
    "Use same_direction or opposite_direction relative to the pathway's",
    "favorable/unfavorable prognostic direction only for actual directional findings.",
    "Only direct human_prognosis evidence with a matched endpoint and an exact",
    "abstract quote may contribute to literature_consistency: same_direction",
    "supports; opposite_direction potentially contradicts; both yield mixed.",
    "Otherwise use insufficient_evidence. Tissue differences, cell/animal models,",
    "reviews, lysine abundance, mTOR links and other indirect biology provide",
    "context only, not validation of the pathway's prognostic direction.",
    "Component or driver-gene retrieval is not evidence for a whole pathway's",
    "activity, flux, direction, or the InterMeta gene-pair network. Respect term_scope.",
    "Directness direct also requires a comparable disease/population and exposure:",
    "a gene signature enriched for a pathway (including a BCR lipid signature",
    "enriched for glycerophospholipid metabolism) is not direct whole-pathway validation.",
    "External TCGA publications may reuse or overlap the analysis cohort; they must",
    "not be described as independent validation without demonstrated cohort independence.",
    "Repeated retrieval of the same publication and reviews of the same studies",
    "are not independent replication. Citation membership and excerpt matching",
    "only establish grounding, not the correctness of semantic classifications;",
    "human review is required. Preserve exploratory, non-causal language."
  )
}

.InterMetaAIEmptyAssessments <- function(include_derived = FALSE) {
  out <- data.frame(
    Pathway = character(), PMID = character(), Evidence.Type = character(),
    Study.Endpoint = character(), Endpoint.Match = character(),
    Direction = character(), Directness = character(), Abstract.Quote = character(),
    Assessment.Status = character(), stringsAsFactors = FALSE
  )
  if (include_derived) {
    out$Quote.Verified <- logical()
    out$Eligible.For.Consistency <- logical()
  }
  out
}

.InterMetaAIAssessmentRows <- function(assessments, pathway) {
  if (is.null(assessments) || !length(assessments)) {
    return(.InterMetaAIEmptyAssessments())
  }
  if (is.data.frame(assessments)) {
    assessments <- lapply(seq_len(nrow(assessments)), function(i) {
      as.list(assessments[i, , drop = FALSE])
    })
  } else if (!is.list(assessments)) {
    stop("Provider assessments must be a list or data frame.", call. = FALSE)
  } else if ("pmid" %in% .InterMetaAINormalizeName(names(assessments))) {
    assessments <- list(assessments)
  }
  if (!length(assessments)) return(.InterMetaAIEmptyAssessments())
  rows <- lapply(assessments, function(item) {
    if (!is.list(item)) {
      stop("Each provider assessment must be a list.", call. = FALSE)
    }
    value <- function(field) {
      .InterMetaAIScalarString(.InterMetaAIField(item, field),
                              paste("assessment", field))
    }
    pmid <- value("pmid")
    if (!grepl("^[0-9]+$", pmid)) {
      stop("Each assessment must contain exactly one numeric PMID.", call. = FALSE)
    }
    quote <- .InterMetaAIField(item, "abstract_quote")
    if (!is.character(quote) || length(quote) != 1L || is.na(quote) ||
        nchar(quote) > 1000L) {
      stop("assessment abstract_quote must be one string of at most 1000 characters.",
           call. = FALSE)
    }
    status <- .InterMetaAIField(item, "Assessment.Status", "provided")
    if (!identical(status, "unassessed")) status <- "provided"
    data.frame(
      Pathway = pathway, PMID = pmid, Evidence.Type = value("evidence_type"),
      Study.Endpoint = value("study_endpoint"), Endpoint.Match = value("endpoint_match"),
      Direction = value("direction"), Directness = value("directness"),
      Abstract.Quote = quote, Assessment.Status = status, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, c(rows, list(make.row.names = FALSE)))
}

.InterMetaAIWhitespace <- function(x) {
  trimws(gsub("[[:space:]]+", " ", x))
}

.InterMetaAIEndpointKey <- function(x) {
  if (is.null(x) || length(x) != 1L || is.na(x)) return("")
  x <- tolower(.InterMetaAIWhitespace(as.character(x)))
  if (!nzchar(x) || x %in% c("unknown", "unspecified", "not specified", "not reported",
                             "not provided", "missing", "none", "na", "n/a",
                             "not available", "unavailable", "unknown endpoint",
                             "unspecified endpoint", "endpoint unknown")) return("")
  # This is one biochemical endpoint, not separate BCR and generic RFS.
  # Replace only the complete phrase so a separately named OS/RFS remains
  # detectable and does not silently become a single-endpoint match.
  x <- gsub(
    "\\bbiochemical[ -]+(?:recurrence|relapse)[ -]+free[ -]+survival\\b",
    "bcr", x, perl = TRUE
  )
  aliases <- c(
    bcr = "\\bbcr\\b|biochemical[ -]+recurrence|biochemical[ -]+relapse",
    os = "\\bos\\b|overall[ -]+survival",
    pfs = "\\bpfs\\b|progression[ -]+free[ -]+survival",
    dfs = "\\bdfs\\b|disease[ -]+free[ -]+survival",
    rfs = "\\brfs\\b|recurrence[ -]+free[ -]+survival|relapse[ -]+free[ -]+survival",
    dss = "\\bdss\\b|disease[ -]+specific[ -]+survival",
    css = "\\bcss\\b|cancer[ -]+specific[ -]+survival"
  )
  matches <- names(aliases)[vapply(aliases, grepl, logical(1), x = x, perl = TRUE)]
  if (length(matches)) return(paste(sort(matches), collapse = "|"))
  .InterMetaAINormalizeName(x)
}

.InterMetaAIGradeAssessments <- function(pathways, claims, assessments,
                                         assessed_pathways, context, evidence,
                                         language = "en") {
  enums <- .InterMetaAIAssessmentEnums()
  for (field in names(enums)) {
    values <- .InterMetaAIField(assessments, field)
    bad <- setdiff(unique(values), enums[[field]])
    if (length(bad)) {
      stop("Invalid assessment ", field, " value(s): ", paste(bad, collapse = ", "),
           ".", call. = FALSE)
    }
  }
  keys <- paste(assessments$Pathway, assessments$PMID, sep = "\r")
  if (anyDuplicated(keys)) {
    stop("Duplicate pathway/PMID assessments are not allowed.", call. = FALSE)
  }
  cited <- unique(claims[, c("Pathway", "PMID"), drop = FALSE])
  cited_keys <- paste(cited$Pathway, cited$PMID, sep = "\r")
  if (any(!keys %in% cited_keys)) {
    stop("Every assessment PMID must be cited for the same pathway.", call. = FALSE)
  }
  missing <- !cited_keys %in% keys
  if (any(missing & cited$Pathway %in% assessed_pathways)) {
    stop("Structured assessments must cover every cited PMID for each pathway.",
         call. = FALSE)
  }
  if (any(missing)) {
    unassessed <- data.frame(
      cited[missing, , drop = FALSE], Evidence.Type = "unknown",
      Study.Endpoint = "unknown", Endpoint.Match = "unknown",
      Direction = "not_assessed", Directness = "indirect", Abstract.Quote = "",
      Assessment.Status = "unassessed", stringsAsFactors = FALSE
    )
    assessments <- rbind(assessments, unassessed)
  }
  rownames(assessments) <- NULL
  # Only Abstract text actually made available to the provider is eligible.
  # Case and punctuation are unchanged; normalization collapses whitespace only.
  assessments$Quote.Verified <- rep(FALSE, nrow(assessments))
  record_pathways <- .InterMetaAIField(evidence$records, "pathway", character())
  for (i in seq_len(nrow(assessments))) {
    quote <- .InterMetaAIWhitespace(assessments$Abstract.Quote[i])
    if (!nzchar(quote)) next
    rows <- which(record_pathways == assessments$Pathway[i] &
                    evidence$records$PMID == assessments$PMID[i])
    abstracts <- if ("Abstract" %in% names(evidence$records)) {
      .InterMetaAIWhitespace(as.character(evidence$records$Abstract[rows]))
    } else character()
    matched <- any(!is.na(abstracts) & vapply(abstracts, function(abstract) {
      !is.na(abstract) && grepl(quote, abstract, fixed = TRUE)
    }, logical(1)))
    if (!matched) {
      stop("assessment abstract_quote is not an exact excerpt of the supplied Abstract for ",
           assessments$Pathway[i], ":", assessments$PMID[i], ".", call. = FALSE)
    }
    assessments$Quote.Verified[i] <- TRUE
  }
  target <- .InterMetaAIField(context$analysis, "endpoint")
  target_key <- .InterMetaAIEndpointKey(target)
  study_keys <- vapply(assessments$Study.Endpoint, .InterMetaAIEndpointKey, character(1))
  endpoint_matches <- nzchar(target_key) & nzchar(study_keys) & study_keys == target_key
  assessments$Eligible.For.Consistency <-
    assessments$Assessment.Status == "provided" &
    assessments$Evidence.Type == "human_prognosis" &
    assessments$Endpoint.Match == "matched" & endpoint_matches &
    assessments$Directness == "direct" & assessments$Quote.Verified

  pathways$Biological.Relevance <- ifelse(
    pathways$Pathway %in% cited$Pathway, "contextual_evidence", "insufficient_evidence")
  endpoint_label <- if (nzchar(target_key)) as.character(target) else "an unspecified endpoint"
  for (i in seq_len(nrow(pathways))) {
    subset <- assessments[assessments$Pathway == pathways$Pathway[i] &
                            assessments$Eligible.For.Consistency, , drop = FALSE]
    same <- any(subset$Direction == "same_direction")
    opposite <- any(subset$Direction == "opposite_direction")
    grade <- if (same && opposite) "mixed" else if (same) "supported" else if (opposite) {
      "potentially_contradicted"
    } else "insufficient_evidence"
    pathways$Literature.Consistency[i] <- grade
    has_claim <- pathways$Pathway[i] %in% cited$Pathway
    pathways$Summary[i] <- if (!has_claim) {
      paste("No cited pathway interpretation is available; the evidence is insufficient",
            "for a biological conclusion or endpoint-specific directional assessment.")
    } else if (grade == "insufficient_evidence") {
      paste0("Cited records provide biological context only; the evidence is insufficient ",
             "to support or contradict the pathway's prognostic direction for ",
             endpoint_label, ". Null, incomparable, indirect, or unassessed findings ",
             "do not establish an opposite prognostic direction.")
    } else {
      directional_text <- switch(grade,
        supported = "the same prognostic direction",
        potentially_contradicted = "an opposite prognostic direction",
        mixed = "both same and opposite prognostic directions")
      paste0("Eligible source assessments indicate ", directional_text, " for ",
             endpoint_label, "; this is exploratory evidence requiring human review, ",
             "not proof of pathway activity, causality, or independent validation.")
    }
    if (identical(language, "zh")) {
      target_text <- if (nzchar(target_key)) as.character(target) else "\u672a\u6307\u5b9a\u7684\u7ec8\u70b9"
      pathways$Summary[i] <- if (!has_claim) {
        "\u76ee\u524d\u6ca1\u6709\u5e26\u5f15\u6587\u7684\u901a\u8def\u89e3\u8bfb\uff1b\u8bc1\u636e\u4e0d\u8db3\u4ee5\u5f62\u6210\u751f\u7269\u5b66\u7ed3\u8bba\u6216\u7ec8\u70b9\u7279\u5f02\u7684\u9884\u540e\u65b9\u5411\u5224\u65ad\u3002"
      } else if (grade == "insufficient_evidence") {
        paste0("\u6240\u5f15\u6587\u732e\u4ec5\u63d0\u4f9b\u5f85\u5ba1\u6838\u7684\u751f\u7269\u5b66\u80cc\u666f\uff0c\u8bc1\u636e\u4e0d\u8db3\u4ee5\u652f\u6301\u6216\u53cd\u9a73\u8be5\u901a\u8def\u9488\u5bf9",
               target_text, "\u7684\u9884\u540e\u65b9\u5411\u3002\u65e0\u5173\u8054\u3001\u4e0d\u53ef\u6bd4\u3001\u95f4\u63a5\u6216\u672a\u7ecf\u7ed3\u6784\u5316\u8bc4\u4f30\u7684\u7ed3\u679c\uff0c",
               "\u4e0d\u80fd\u89c6\u4e3a\u76f8\u53cd\u9884\u540e\u65b9\u5411\u7684\u8bc1\u636e\u3002")
      } else {
        directional_text <- switch(grade, supported = "\u76f8\u540c\u7684\u9884\u540e\u65b9\u5411",
          potentially_contradicted = "\u76f8\u53cd\u7684\u9884\u540e\u65b9\u5411", mixed = "\u76f8\u540c\u53ca\u76f8\u53cd\u7684\u9884\u540e\u65b9\u5411")
        paste0("\u7b26\u5408\u6761\u4ef6\u7684\u6587\u732e\u8bc4\u4f30\u63d0\u793a\u9488\u5bf9", target_text, "\u5b58\u5728", directional_text,
               "\uff1b\u8fd9\u4ecd\u662f\u9700\u8981\u4eba\u5de5\u590d\u6838\u7684\u63a2\u7d22\u6027\u8bc1\u636e\uff0c\u4e0d\u80fd\u8bc1\u660e\u901a\u8def\u6d3b\u6027\u3001\u56e0\u679c\u5173\u7cfb\u6216\u72ec\u7acb\u9a8c\u8bc1\u3002")
      }
    }
  }
  warnings <- character()
  if (any(assessments$Assessment.Status == "unassessed")) {
    warnings <- c(warnings, paste(
      "Legacy citations without structured assessments were retained as unassessed",
      "context only and cannot support or contradict endpoint-specific prognosis."
    ))
  }
  if (!nzchar(target_key) && nrow(pathways)) {
    warnings <- c(warnings, "The target endpoint is unspecified; endpoint-specific directional evidence cannot be graded.")
  }
  if (any(pathways$Model.Literature.Consistency != pathways$Literature.Consistency)) {
    warnings <- c(warnings, paste(
      "Rule-derived endpoint grades differ from model classifications; raw",
      "Model.Literature.Consistency and Model.Summary are retained for audit only."
    ))
  }
  list(pathways = pathways, assessments = assessments, warnings = warnings)
}
