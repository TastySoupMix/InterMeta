# Pathway-level biological synthesis is distinct from endpoint direction.
# Model classifications are source-grounded, not automatically proven true.

.InterMetaAISynthesisSchema <- function() {
  enums <- .InterMetaAIBiologyEnums()
  fields <- lapply(enums, function(x) list(type = "string", enum = x))
  finding <- list(type = "object", additionalProperties = FALSE,
    required = c("pmid", "study_disease", names(enums), "claim", "abstract_quote"),
    properties = c(list(pmid = list(type = "string"),
      study_disease = list(type = "string", minLength = 1L)), fields,
      list(claim = list(type = "string", minLength = 1L),
           abstract_quote = list(type = "string", maxLength = 1000L))))
  list(type = "object", additionalProperties = FALSE,
    required = c("mechanism_summary", "limitations", "findings"),
    properties = list(
      mechanism_summary = list(type = "string"), limitations = list(type = "string"),
      findings = list(type = "array", items = finding)))
}

.InterMetaAIBiologyEnums <- function() list(
  disease_match = c("matched", "other", "unknown"),
  scope = c("pathway", "component", "unclear"),
  finding_type = c("mechanistic", "associational", "background", "unknown"),
  study_design = c("experimental", "observational", "review", "unknown"),
  relation = c("supports", "opposes", "neutral", "unclear")
)

.InterMetaAISynthesisPrompt <- function() paste(
  "Also return biology_synthesis for each pathway: a concise cross-publication",
  "mechanism_summary, limitations, and one finding per selected cited PMID.",
  "This is biological explanation of pathway involvement, NOT BCR direction validation.",
  "InterMeta Difference is normalized outgoing minus incoming gene-pair network degree;",
  "it is not pathway activation, flux, or a measured biological exposure.",
  "Every biology finding PMID must already occur in claims.pmids and assessments.",
  "Use actual study_disease, disease_match, scope, finding_type, study_design, relation,",
  "a concise 25-35 word claim, and a short exact supplied Abstract excerpt.",
  "Use an empty quote if no relevant Abstract excerpt exists. Never quote Title instead.",
  "disease_match matched requires the target disease, not another cancer or cancer in general.",
  "scope pathway requires evidence about the named whole pathway; single genes, signatures,",
  "enrichment, metabolites, lysine abundance and mTOR associations are component or unclear.",
  "finding_type mechanistic plus study_design experimental requires an actual perturbation",
  "or mechanistic experiment, not expression differences or outcome associations.",
  "The particular claim must demonstrate a mechanism-mediated effect; a paper containing",
  "other experiments does not make every association within it mechanistic.",
  "Reviews are study_design review and background; do not relabel a reviewed experiment",
  "as an independent experimental publication. Reused cohorts are not independent validation.",
  "relation supports/opposes refers to the proposed biological mechanism, not the BCR",
  "or signed network direction. Null/no-association findings are neutral, never opposes.",
  "Use opposes only for an explicitly comparable proposed mechanism and study conditions;",
  "different doses, perturbations, cell types or disease systems alone are not opposition.",
  "Discuss genuine comparable counter-evidence, component-to-pathway mapping gaps,",
  "study systems and limitations; incomparable studies are not opposing evidence.",
  "No inline PMIDs, URLs or DOI identifiers in any narrative field; use structured pmid fields.",
  "These are model-assisted classifications requiring human review, not semantic proof."
)

.InterMetaAIEmptySynthesis <- function() data.frame(
  Pathway = character(), Synthesis.Status = character(),
  Model.Mechanism.Summary = character(), Model.Limitations = character(),
  Assessment.Source = character(),
  stringsAsFactors = FALSE)

.InterMetaAIEmptyBiologyFindings <- function() data.frame(
  Pathway = character(), PMID = character(), Study.Disease = character(),
  Disease.Match = character(), Scope = character(), Finding.Type = character(),
  Study.Design = character(), Relation = character(), Claim = character(),
  Abstract.Quote = character(), stringsAsFactors = FALSE)

.InterMetaAISynthesisText <- function(x, name, empty = TRUE) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      (!empty && !nzchar(trimws(x)))) stop(name, " must be one character value.", call. = FALSE)
  trimws(x)
}

.InterMetaAIBiologyRows <- function(items, pathway) {
  if (is.null(items)) stop("biology_synthesis findings must be an array.", call. = FALSE)
  if (is.data.frame(items)) items <- lapply(seq_len(nrow(items)), function(i) as.list(items[i, , drop = FALSE]))
  if (!is.list(items)) stop("biology_synthesis findings must be an array.", call. = FALSE)
  if (!length(items)) return(.InterMetaAIEmptyBiologyFindings())
  if ("pmid" %in% .InterMetaAINormalizeName(names(items))) items <- list(items)
  rows <- lapply(items, function(item) {
    if (!is.list(item)) stop("Each biology finding must be an object.", call. = FALSE)
    value <- function(name, empty = FALSE) .InterMetaAISynthesisText(
      .InterMetaAIField(item, name), paste("biology finding", name), empty)
    pmid <- value("pmid")
    if (!grepl("^[0-9]+$", pmid)) stop("Each biology finding requires one numeric PMID.", call. = FALSE)
    quote <- value("abstract_quote", TRUE)
    if (nchar(quote) > 1000L) stop("biology abstract_quote exceeds 1000 characters.", call. = FALSE)
    data.frame(Pathway = pathway, PMID = pmid, Study.Disease = value("study_disease"),
      Disease.Match = value("disease_match"), Scope = value("scope"),
      Finding.Type = value("finding_type"), Study.Design = value("study_design"),
      Relation = value("relation"), Claim = value("claim"), Abstract.Quote = quote,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, c(rows, list(make.row.names = FALSE)))
}

.InterMetaAINormalizeSynthesis <- function(raw, entries, pathways) {
  rows <- list()
  findings <- .InterMetaAIEmptyBiologyFindings()
  source <- attr(raw, "intermeta_assessment_source", exact = TRUE)
  if (!identical(source, "assistant_abstract_audit")) source <- "model_assisted"
  for (entry in entries) {
    if (!"biologysynthesis" %in% .InterMetaAINormalizeName(names(entry))) next
    pathway <- .InterMetaAIScalarString(.InterMetaAIField(entry, c("pathway", "pathway_name")), "pathway")
    biology <- .InterMetaAIField(entry, "biology_synthesis")
    if (!is.list(biology)) stop("biology_synthesis must be an object.", call. = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(Pathway = pathway, Synthesis.Status = "provided",
      Model.Mechanism.Summary = .InterMetaAISynthesisText(biology$mechanism_summary, "mechanism_summary"),
      Model.Limitations = .InterMetaAISynthesisText(biology$limitations, "limitations"),
      Assessment.Source = source, stringsAsFactors = FALSE)
    findings <- rbind(findings, .InterMetaAIBiologyRows(biology$findings, pathway))
  }
  if (is.list(raw) && !is.null(raw$pathway_synthesis)) {
    table <- raw$pathway_synthesis
    if (!is.data.frame(table)) stop("pathway_synthesis must be a data frame.", call. = FALSE)
    for (i in seq_len(nrow(table))) {
      row <- as.list(table[i, , drop = FALSE])
      rows[[length(rows) + 1L]] <- data.frame(
        Pathway = .InterMetaAIScalarString(row$Pathway, "synthesis Pathway"),
        Synthesis.Status = .InterMetaAIScalarString(row$Synthesis.Status, "Synthesis.Status"),
        Model.Mechanism.Summary = .InterMetaAISynthesisText(row$Model.Mechanism.Summary, "Model.Mechanism.Summary"),
        Model.Limitations = .InterMetaAISynthesisText(row$Model.Limitations, "Model.Limitations"),
        Assessment.Source = if (!is.null(row$Assessment.Source) && row$Assessment.Source %in%
          c("model_assisted", "assistant_abstract_audit", "not_assessed")) row$Assessment.Source else source,
        stringsAsFactors = FALSE)
    }
  }
  if (is.list(raw) && !is.null(raw$biology_findings)) {
    table <- raw$biology_findings
    if (!is.data.frame(table)) stop("biology_findings must be a data frame.", call. = FALSE)
    for (i in seq_len(nrow(table))) {
      row <- as.list(table[i, , drop = FALSE])
      findings <- rbind(findings, .InterMetaAIBiologyRows(row,
        .InterMetaAIScalarString(row$Pathway, "finding Pathway")))
    }
  }
  synthesis <- if (length(rows)) do.call(rbind, c(rows, list(make.row.names = FALSE))) else .InterMetaAIEmptySynthesis()
  rownames(findings) <- NULL
  list(pathway_synthesis = synthesis, biology_findings = findings)
}

.InterMetaAIDiseaseKey <- function(x) {
  if (is.null(x) || length(x) != 1L || is.na(x)) return("")
  x <- tolower(.InterMetaAIWhitespace(as.character(x)))
  if (!nzchar(x) || x %in% c("unknown", "not specified", "unspecified", "cancer",
                             "tumor", "tumour", "none", "na", "n/a")) return("")
  # Lexical aliases only; biological population comparability remains a
  # model-assisted assessment, not a disease ontology or semantic verifier.
  # Match a complete single-disease phrase, never just one known organ within
  # a longer mixed-population label (e.g. prostate AND bladder cancer).
  single <- sub("^(?:(?:human|primary|advanced|metastatic|castration[- ]resistant|castration[- ]sensitive|hormone[- ]sensitive) +)+",
    "", x, perl = TRUE)
  malignancy <- "(?:cancer|carcinoma|adenocarcinoma|neoplasms?|tumou?rs?)"
  aliases <- c(prostate = paste0("^(?:(?:prostate|prostatic) ", malignancy, "|prad)$"),
    breast = paste0("^breast ", malignancy, "$"),
    lung = paste0("^lung ", malignancy, "$"),
    colorectal = paste0("^(?:colorectal|colon|rectal) ", malignancy, "$"))
  hits <- names(aliases)[vapply(aliases, grepl, logical(1), x = single, perl = TRUE)]
  if (length(hits) == 1L) return(hits)
  .InterMetaAINormalizeName(x)
}

.InterMetaAIValidateSynthesis <- function(normalized, context, evidence, graded,
                                         claims, language = "en") {
  synthesis <- normalized$pathway_synthesis
  findings <- normalized$biology_findings
  known <- as.character(context$pathways$Pathway)
  if (anyDuplicated(synthesis$Pathway) || any(!synthesis$Pathway %in% known)) {
    stop("Duplicate or unknown pathway in pathway synthesis.", call. = FALSE)
  }
  if (any(!synthesis$Synthesis.Status %in% c("provided", "not_assessed"))) {
    stop("Invalid Synthesis.Status.", call. = FALSE)
  }
  enums <- .InterMetaAIBiologyEnums()
  for (field in names(enums)) {
    if (any(!.InterMetaAIField(findings, field) %in% enums[[field]])) {
      stop("Invalid biology finding ", field, ".", call. = FALSE)
    }
  }
  narrative <- c(synthesis$Model.Mechanism.Summary, synthesis$Model.Limitations, findings$Claim)
  if (any(grepl("\\bPMIDs?\\b|https?://|www[.]|\\b10[.][0-9]{4,9}/|\\[[[:space:]]*[0-9]{4,9}[[:space:]]*\\]",
                narrative, ignore.case = TRUE, perl = TRUE))) {
    stop("Inline citation identifiers or URLs are not allowed in biology narrative; use structured PMID fields.", call. = FALSE)
  }
  pairs <- paste(findings$Pathway, findings$PMID, sep = "\r")
  if (anyDuplicated(pairs)) stop("Duplicate pathway/PMID biology findings are not allowed.", call. = FALSE)
  claim_pairs <- paste(claims$Pathway, claims$PMID, sep = "\r")
  assessment_pairs <- paste(graded$assessments$Pathway, graded$assessments$PMID, sep = "\r")
  if (any(!pairs %in% claim_pairs) || any(!pairs %in% assessment_pairs)) {
    stop("Each biology finding PMID must be cited and assessed for the same pathway.", call. = FALSE)
  }
  assessment_index <- match(pairs, assessment_pairs)
  if (any(graded$assessments$Assessment.Status[assessment_index] != "provided")) {
    stop("Each biology finding PMID requires a provided structured evidence assessment; legacy unassessed placeholders are not sufficient.", call. = FALSE)
  }
  provided <- synthesis$Pathway[synthesis$Synthesis.Status == "provided"]
  if (any(!findings$Pathway %in% provided)) stop("Biology findings require a provided synthesis for that pathway.", call. = FALSE)
  records <- evidence$records
  record_paths <- .InterMetaAIField(records, "pathway", character())
  findings$Quote.Verified <- rep(FALSE, nrow(findings))
  for (i in seq_len(nrow(findings))) {
    quote <- .InterMetaAIWhitespace(findings$Abstract.Quote[i])
    if (!nzchar(quote)) next
    hit <- which(record_paths == findings$Pathway[i] & records$PMID == findings$PMID[i])
    abstracts <- if ("Abstract" %in% names(records)) records$Abstract[hit] else character()
    matched <- any(vapply(abstracts, function(x) !is.na(x) &&
      grepl(quote, .InterMetaAIWhitespace(x), fixed = TRUE), logical(1)))
    if (!matched) stop("biology abstract_quote is not an exact excerpt of the supplied Abstract.", call. = FALSE)
    findings$Quote.Verified[i] <- TRUE
  }
  disease_key <- .InterMetaAIDiseaseKey(context$analysis$disease)
  keys <- vapply(findings$Study.Disease, .InterMetaAIDiseaseKey, character(1))
  findings$Disease.Match.Verified <- nzchar(disease_key) & nzchar(keys) & keys == disease_key &
    findings$Disease.Match == "matched"
  matched_assessments <- assessment_index
  findings$Evidence.Type <- graded$assessments$Evidence.Type[matched_assessments]
  # Biology and clinical outcome are separate dimensions. A study may show a
  # mechanism while its clinical endpoint association is null (or vice versa).
  # Only the biological relation marks this finding as neutral/null.
  findings$Null.Finding <- findings$Relation == "neutral"
  grounded <- findings$Quote.Verified & findings$Disease.Match.Verified
  findings$Eligible.For.Mechanism <- grounded & findings$Scope == "pathway" &
    findings$Finding.Type == "mechanistic" & findings$Study.Design == "experimental" &
    findings$Evidence.Type %in% c("human_prognosis", "human_tissue", "preclinical") &
    !findings$Null.Finding & findings$Relation %in% c("supports", "opposes")
  # Emit one row per selected context pathway, including missing model outputs.
  ordered <- .InterMetaAIEmptySynthesis()
  for (pathway in known) {
    hit <- match(pathway, synthesis$Pathway)
    row <- if (!is.na(hit)) synthesis[hit, , drop = FALSE] else data.frame(
      Pathway = pathway, Synthesis.Status = "not_assessed", Model.Mechanism.Summary = "",
      Model.Limitations = "", Assessment.Source = "not_assessed", stringsAsFactors = FALSE)
    ordered <- rbind(ordered, row)
  }
  synthesis <- ordered
  synthesis$Biological.Explanation <- rep("not_assessed", nrow(synthesis))
  synthesis$Mechanism.Summary <- rep("", nrow(synthesis))
  synthesis$Limitations <- rep("", nrow(synthesis))
  synthesis$Key.PMIDs <- rep("", nrow(synthesis))
  for (i in seq_len(nrow(synthesis))) {
    indices <- which(findings$Pathway == synthesis$Pathway[i] & grounded)
    mechanism <- indices[findings$Eligible.For.Mechanism[indices]]
    same <- mechanism[findings$Relation[mechanism] == "supports"]
    opposite <- mechanism[findings$Relation[mechanism] == "opposes"]
    grade <- if (synthesis$Synthesis.Status[i] == "not_assessed") "not_assessed" else
      if (!length(indices)) "insufficient_evidence" else
      if (length(same) && length(opposite)) "conflicting_evidence" else
      if (length(same)) "pathway_mechanism" else "indirect_context"
    synthesis$Biological.Explanation[i] <- grade
    chosen <- if (grade == "conflicting_evidence") c(same[1L], opposite[1L]) else
      head(unique(c(mechanism, indices)), 2L)
    if (length(chosen)) {
      texts <- vapply(chosen, function(j) {
        prefix <- if (findings$Null.Finding[j] || findings$Relation[j] == "neutral") {
          if (language == "zh") "\u4e2d\u6027\u6216\u65e0\u5173\u8054\u7ed3\u679c\uff1a" else "Neutral finding: "
        } else if (findings$Relation[j] == "opposes") {
          if (language == "zh") "\u4e0d\u540c\u65b9\u5411\u7684\u673a\u5236\u7ebf\u7d22\uff1a" else "Contrasting finding: "
        } else ""
        paste0(prefix, findings$Claim[j])
      }, character(1))
      synthesis$Mechanism.Summary[i] <- paste(texts, collapse = " ")
      synthesis$Key.PMIDs[i] <- paste(findings$PMID[chosen], collapse = ";")
    } else {
      synthesis$Mechanism.Summary[i] <- if (grade == "not_assessed") {
        if (language == "zh") "\u672c\u6b21\u8f93\u51fa\u5c1a\u672a\u8bc4\u4f30\u901a\u8def\u5c42\u9762\u7684\u751f\u7269\u5b66\u7efc\u5408\u89e3\u91ca\u3002" else
          "Pathway-level biological synthesis was not assessed in this provider output."
      } else if (language == "zh") "\u5c1a\u65e0\u540c\u65f6\u6ee1\u8db3\u76ee\u6807\u75c5\u79cd\u5339\u914d\u4e0e\u6458\u8981\u51fa\u5904\u6838\u9a8c\u7684\u751f\u7269\u5b66\u53d1\u73b0\uff0c\u65e0\u6cd5\u5f62\u6210\u901a\u8def\u7efc\u5408\u89e3\u91ca\u3002" else
        "No source-grounded finding with a verified target-disease match is available for biological synthesis."
    }
    synthesis$Limitations[i] <- trimws(paste(
      if (grade == "conflicting_evidence") "Context-dependent or mixed biological findings are not automatically a direct experimental contradiction." else "",
      if (grade == "indirect_context") "Evidence is indirect, component-level, observational, background, null, or counter-evidence only." else "",
      "Source-linked, model-assisted interpretation requires human review; disease and mechanism semantics are not automatically verified.",
      "Neither this label nor signed Difference establishes pathway activity, flux, causality or independent endpoint-specific prognostic validation."))
    if (language == "zh") synthesis$Limitations[i] <- paste0(
      if (grade == "conflicting_evidence") "\u4e0d\u540c\u8bed\u5883\u4e0b\u7684\u673a\u5236\u7ebf\u7d22\u53ef\u80fd\u5e76\u4e0d\u4e00\u81f4\uff0c\u4e0d\u80fd\u81ea\u52a8\u8ba4\u5b9a\u5b9e\u9a8c\u76f8\u4e92\u53cd\u9a73\u3002" else "",
      if (grade == "indirect_context") "\u76ee\u524d\u4ec5\u6709\u7ec4\u5206\u3001\u5173\u8054\u3001\u80cc\u666f\u3001\u4e2d\u6027\u6216\u5176\u4ed6\u95f4\u63a5\u7ebf\u7d22\u3002" else "",
      "\u7ed3\u8bba\u4e3a\u53ef\u8ffd\u6eaf\u7684\u6a21\u578b\u8f85\u52a9\u7efc\u5408\uff0c\u75c5\u79cd\u5339\u914d\u53ca\u673a\u5236\u542b\u4e49\u4ecd\u9700\u4eba\u5de5\u590d\u6838\uff1b\u8bc4\u5206\u4e0d\u4ee3\u8868\u901a\u8def\u6d3b\u6027\u3001\u901a\u91cf\u3001\u56e0\u679c\u5173\u7cfb\u6216\u72ec\u7acb\u7684\u7ec8\u70b9\u9884\u540e\u9a8c\u8bc1\u3002")
  }
  synthesis$Literature.Consistency <- graded$pathways$Literature.Consistency[
    match(synthesis$Pathway, graded$pathways$Pathway)]
  synthesis$Literature.Consistency[is.na(synthesis$Literature.Consistency)] <- "insufficient_evidence"
  for (field in intersect(c("Direction", "Direction.Rank", "Difference"), names(context$pathways))) {
    synthesis[[field]] <- context$pathways[[field]][match(synthesis$Pathway, context$pathways$Pathway)]
  }
  synthesis <- synthesis[, c("Pathway", intersect(c("Direction", "Direction.Rank", "Difference"), names(synthesis)),
    "Biological.Explanation", "Literature.Consistency", "Mechanism.Summary", "Limitations",
    "Key.PMIDs", "Model.Mechanism.Summary", "Model.Limitations", "Synthesis.Status", "Assessment.Source"), drop = FALSE]
  findings$Assessment.Source <- synthesis$Assessment.Source[match(findings$Pathway, synthesis$Pathway)]
  rownames(synthesis) <- NULL
  list(pathway_synthesis = synthesis, biology_findings = findings)
}

.InterMetaAISourceContextSHA256 <- function(safe_context) {
  digest::digest(safe_context, algo = "sha256", serialize = TRUE)
}
