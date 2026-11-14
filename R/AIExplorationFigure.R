#' Render the fixed InterMeta AI exploration figure
#'
#' Render the publication-style summary of a validated tumour-biology
#' exploration report.  Rendering is entirely local: this function performs no
#' PubMed retrieval and makes no model request.  Evidence-channel counts are
#' non-exclusive counts of unique PMID records, not evidence grades.
#'
#' @param report A validated `InterMetaAIReport` whose purpose is
#'   `"biology_exploration"`.  Its `top_pathways` must contain one to five
#'   pathways from one prognostic direction.
#' @param pdf_file Optional path to a new PDF file.
#' @param png_file Optional path to a new PNG file.
#' @param render Whether to render files.  Use `FALSE` to prepare and inspect
#'   the figure data without opening a graphics device.
#' @param font_family `NULL`, or one installed font-family name.  ASCII-only
#'   output uses the ordinary `sans` device font.  Non-ASCII output and custom
#'   fonts require Cairo, the optional `systemfonts` package, and complete glyph
#'   coverage in both regular and bold faces.
#' @param text_mode `"compact_english"` (the publication default) or `"full"`.
#'   Compact mode leaves the complete narrative in the report and displays a
#'   short source-count summary.
#'
#' @return Invisibly, a list containing the validated display tables, caption,
#'   dimensions, font audit, evidence-source audit, output paths, and a local
#'   drawing function.
#' @export
InterMetaAIExplorationFigure <- function(
    report, pdf_file = NULL, png_file = NULL, render = TRUE,
    font_family = NULL,
    text_mode = c("compact_english", "full")) {
  fail <- function(...) stop(..., call. = FALSE)
  scalar <- function(x, default = "unknown") {
    if (is.null(x) || !length(x) || is.na(x[1L]) ||
        !nzchar(trimws(as.character(x[1L])))) {
      return(default)
    }
    trimws(as.character(x[1L]))
  }
  columns <- function(x, required, label) {
    if (!is.data.frame(x) || !all(required %in% names(x))) {
      fail(label, " must contain: ", paste(required, collapse = ", "), ".")
    }
  }
  wrap <- function(x, width) {
    text <- gsub("[[:space:]]+", " ", enc2utf8(x))
    hard_wrap <- function(line) {
      if (nchar(line, type = "width") <= width) return(line)
      chars <- strsplit(line, "", fixed = TRUE)[[1L]]
      sizes <- nchar(chars, type = "width")
      chunks <- character()
      start <- 1L
      used <- 0L
      for (j in seq_along(chars)) {
        if (used + sizes[j] > width && j > start) {
          chunks <- c(chunks, paste(chars[start:(j - 1L)], collapse = ""))
          start <- j
          used <- 0L
        }
        used <- used + sizes[j]
      }
      c(chunks, paste(chars[start:length(chars)], collapse = ""))
    }
    lines <- strwrap(text, width = width)
    if (any(utf8ToInt(text) > 127L)) {
      tokens <- regmatches(
        text,
        gregexpr("[A-Za-z0-9_./-]+| +|.", text, perl = TRUE)
      )[[1L]]
      tokens <- unlist(lapply(tokens, hard_wrap), use.names = FALSE)
      lines <- character()
      current <- ""
      used <- 0L
      for (token in tokens) {
        size <- nchar(token, type = "width")
        if (used + size > width && nzchar(current)) {
          lines <- c(lines, trimws(current))
          current <- ""
          used <- 0L
        }
        if (!nzchar(current) && !nzchar(trimws(token))) next
        current <- paste0(current, token)
        used <- used + size
      }
      if (nzchar(current)) lines <- c(lines, trimws(current))
    }
    paste(unlist(lapply(lines, hard_wrap), use.names = FALSE), collapse = "\n")
  }
  pmids <- function(x) {
    if (length(x) != 1L || is.na(x)) {
      fail("PMID fields must be non-missing scalar strings.")
    }
    ids <- unique(strsplit(trimws(as.character(x)), "[;,[:space:]]+")[[1L]])
    ids <- ids[nzchar(ids)]
    if (any(!grepl("^[0-9]+$", ids))) {
      fail("PMID fields must contain numeric identifiers only.")
    }
    ids
  }
  output_path <- function(x, label, extension) {
    if (!is.character(x) || length(x) != 1L || is.na(x) ||
        !nzchar(trimws(x)) || grepl("[[:cntrl:]]", x)) {
      fail(label, " must be one new filename.")
    }
    if (!grepl(paste0("\\.", extension, "$"), x, ignore.case = TRUE)) {
      fail(label, " must end in .", extension, ".")
    }
    if (file.exists(x) || dir.exists(x)) fail("Refusing to overwrite: ", x)
    if (!dir.exists(dirname(x))) {
      fail("Output directory does not exist: ", dirname(x))
    }
    file.path(
      normalizePath(dirname(x), winslash = "/", mustWork = TRUE), basename(x)
    )
  }

  if (!is.logical(render) || length(render) != 1L || is.na(render)) {
    fail("render must be TRUE or FALSE.")
  }
  if (!is.character(text_mode) || !length(text_mode) || anyNA(text_mode)) {
    fail("text_mode must be compact_english or full.")
  }
  text_mode <- tryCatch(match.arg(text_mode), error = function(e) {
    fail("text_mode must be compact_english or full.")
  })
  if (!is.null(font_family) &&
      (!is.character(font_family) || length(font_family) != 1L ||
       is.na(font_family) || !nzchar(trimws(font_family)) ||
       grepl("[[:cntrl:]]", font_family))) {
    fail("font_family must be NULL or one installed font-family name.")
  }
  if (!inherits(report, "InterMetaAIReport")) {
    fail("report must be a validated InterMetaAIReport.")
  }
  if (!identical(report$purpose, "biology_exploration")) {
    fail("report must have purpose = \"biology_exploration\".")
  }
  schema_version <- scalar(report$schema_version, default = "")
  schema_major <- suppressWarnings(as.integer(
    strsplit(schema_version, ".", fixed = TRUE)[[1L]][1L]
  ))
  if (!nzchar(schema_version) || is.na(schema_major) || schema_major < 4L) {
    fail("report must use a supported InterMeta AI report schema (major version 4 or later).")
  }

  required <- c(
    "Pathway", "Direction", "Direction.Rank", "Difference",
    "Biological.Summary", "Biological.Processes", "Context.Dependence",
    "Research.Hypothesis", "Functional.PMIDs", "Omics.PMIDs",
    "Patient.PMIDs", "Review.PMIDs", "Key.PMIDs", "Unique.PMID.Count",
    "Assessment.Source", "Exploration.Status"
  )
  summary <- report$exploration_summary
  columns(summary, required, "report$exploration_summary (report schema >= 4)")
  plausibility_levels <- c(
    "strongly_plausible", "plausible", "mixed_context_dependent",
    "limited_evidence", "not_assessable"
  )
  if (!"Biological.Plausibility" %in% names(summary)) {
    # Reports created before the two-stage 1.1.0 workflow remain drawable, but
    # the renderer must not infer a plausibility result that was never supplied.
    summary$Biological.Plausibility <- rep("not_assessable", nrow(summary))
  } else {
    summary$Biological.Plausibility <- as.character(
      summary$Biological.Plausibility
    )
    if (anyNA(summary$Biological.Plausibility) ||
        any(!summary$Biological.Plausibility %in% plausibility_levels)) {
      fail(
        "Biological.Plausibility must use strongly_plausible, plausible, ",
        "mixed_context_dependent, limited_evidence or not_assessable."
      )
    }
  }
  original <- report$top_pathways
  columns(
    original,
    c("Pathway", "Direction", "Direction.Rank", "Difference"),
    "report$top_pathways"
  )
  if (!nrow(original) || nrow(original) > 5L) {
    fail("Use 1 to 5 selected pathways from one direction for the main figure.")
  }
  if (anyNA(original$Pathway) || anyDuplicated(original$Pathway) ||
      any(!nzchar(trimws(as.character(original$Pathway)))) ||
      anyNA(summary$Pathway) || anyDuplicated(summary$Pathway)) {
    fail("Pathway identifiers must be unique and non-empty.")
  }
  original$Pathway <- as.character(original$Pathway)
  summary$Pathway <- as.character(summary$Pathway)
  if (anyNA(original$Direction) || length(unique(original$Direction)) != 1L ||
      !all(original$Direction %in% c("unfavorable", "favorable"))) {
    fail("Use one unfavorable or favorable direction.")
  }
  if (!is.numeric(original$Difference) || any(!is.finite(original$Difference)) ||
      !is.numeric(summary$Difference) || any(!is.finite(summary$Difference))) {
    fail("Difference must be finite numeric values.")
  }
  direction <- as.character(original$Direction[1L])
  if ((direction == "unfavorable" && any(original$Difference >= 0)) ||
      (direction == "favorable" && any(original$Difference <= 0))) {
    fail("Direction is inconsistent with the sign of Difference.")
  }
  if (!is.numeric(original$Direction.Rank) || anyNA(original$Direction.Rank) ||
      any(original$Direction.Rank < 1 |
          original$Direction.Rank != floor(original$Direction.Rank)) ||
      anyDuplicated(original$Direction.Rank) ||
      is.unsorted(original$Direction.Rank, strictly = TRUE)) {
    fail("Direction.Rank must contain unique increasing positive integers.")
  }
  index <- match(original$Pathway, summary$Pathway)
  if (anyNA(index)) fail("A selected pathway is missing from exploration_summary.")
  state_fields <- intersect(
    c("Quote.Grounding", "Evidence.Eligibility"), names(summary)
  )
  rows <- summary[
    index,
    c(required, "Biological.Plausibility", state_fields),
    drop = FALSE
  ]
  rownames(rows) <- NULL
  if (!is.numeric(rows$Direction.Rank) || anyNA(rows$Direction.Rank) ||
      any(rows$Difference != original$Difference) ||
      !identical(as.character(rows$Direction), as.character(original$Direction)) ||
      any(rows$Direction.Rank != original$Direction.Rank)) {
    fail("Exploration scores, direction or rank disagree with original top_pathways.")
  }
  text_fields <- setdiff(
    required, c("Difference", "Direction.Rank", "Unique.PMID.Count")
  )
  for (field in text_fields) {
    if (anyNA(rows[[field]])) {
      fail(field, " cannot contain NA; preserve explicit unassessed states.")
    }
    rows[[field]] <- as.character(rows[[field]])
  }
  statuses <- c(
    "provided", "not_assessed", "no_relevant_findings",
    "ungrounded_findings", "input_failed"
  )
  if (any(!rows$Exploration.Status %in% statuses)) {
    fail("Unknown Exploration.Status.")
  }
  grounding_states <- c("complete", "partial", "none", "no_findings")
  for (field in state_fields) {
    if (!is.character(rows[[field]]) || anyNA(rows[[field]]) ||
        any(!rows[[field]] %in% grounding_states)) {
      fail(field, " must use complete, partial, none or no_findings.")
    }
  }

  records <- report$evidence$records
  columns(records, c("pathway", "PMID"), "report$evidence$records")
  if (anyNA(records$pathway) || anyNA(records$PMID)) {
    fail("Evidence record pathway and PMID values cannot be missing.")
  }
  record_keys <- paste(
    as.character(records$pathway), as.character(records$PMID), sep = "\r"
  )
  findings <- report$exploration_findings
  columns(
    findings,
    c(
      "Pathway", "Finding.ID", "PMID", "Scope", "Eligible.For.Summary",
      "Functional", "Omics", "Patient", "Review"
    ),
    "report$exploration_findings"
  )
  for (field in c(
    "Eligible.For.Summary", "Functional", "Omics", "Patient", "Review"
  )) {
    if (!is.logical(findings[[field]]) || anyNA(findings[[field]])) {
      fail("Exploration finding flags must be explicit logical values: ", field)
    }
  }
  quote_grounding_available <- "Quote.Verified" %in% names(findings)
  if (quote_grounding_available &&
      (!is.logical(findings$Quote.Verified) || anyNA(findings$Quote.Verified))) {
    fail("Quote.Verified must contain explicit logical values when supplied.")
  }
  if (quote_grounding_available &&
      any(findings$Eligible.For.Summary & !findings$Quote.Verified)) {
    fail(paste(
      "Eligible.For.Summary cannot include a finding whose strict quote is",
      "not verified."
    ))
  }
  if (anyNA(findings$Pathway) || anyNA(findings$Finding.ID) ||
      anyNA(findings$PMID) || anyNA(findings$Scope) ||
      any(!nzchar(trimws(as.character(findings$Finding.ID)))) ||
      any(!grepl("^[0-9]+$", as.character(findings$PMID))) ||
      any(!findings$Scope %in% c("pathway", "component", "unclear")) ||
      anyDuplicated(paste(findings$Pathway, findings$Finding.ID, sep = "\r"))) {
    fail(paste(
      "Exploration findings require unique non-empty finding identifiers,",
      "numeric PMIDs and explicit scope."
    ))
  }
  findings$Pathway <- as.character(findings$Pathway)
  findings$PMID <- as.character(findings$PMID)

  channels <- c(
    Functional = "Functional.PMIDs", Omics = "Omics.PMIDs",
    Patient = "Patient.PMIDs", Review = "Review.PMIDs"
  )
  channel_counts <- data.frame(
    Pathway = rows$Pathway, Functional = NA_integer_, Omics = NA_integer_,
    Patient = NA_integer_, Review = NA_integer_, Unique = NA_integer_,
    stringsAsFactors = FALSE
  )
  channel_sets <- vector("list", nrow(rows))
  channel_scopes <- matrix(
    "", nrow = nrow(rows), ncol = 3,
    dimnames = list(NULL, c("Functional", "Omics", "Patient"))
  )
  references <- character(nrow(rows))
  grounding_counts <- data.frame(
    Pathway = rows$Pathway,
    Total.Findings = rep(NA_integer_, nrow(rows)),
    Quote.Verified = rep(NA_integer_, nrow(rows)),
    Eligible.For.Summary = rep(NA_integer_, nrow(rows)),
    Quote.Quarantined = rep(NA_integer_, nrow(rows)),
    Disease.Excluded = rep(NA_integer_, nrow(rows)),
    Status = rep("unavailable", nrow(rows)),
    Eligibility.Status = rep("unavailable", nrow(rows)),
    Display = rep("", nrow(rows)), stringsAsFactors = FALSE
  )
  scope_labels <- c(pathway = "P", component = "C", unclear = "?")

  for (i in seq_len(nrow(rows))) {
    sets <- lapply(channels, function(field) pmids(rows[[field]][i]))
    keys <- pmids(rows$Key.PMIDs[i])
    pathway_findings <- findings[
      findings$Pathway == rows$Pathway[i], , drop = FALSE
    ]
    eligible <- pathway_findings[
      pathway_findings$Eligible.For.Summary, , drop = FALSE
    ]
    if (quote_grounding_available) {
      total <- nrow(pathway_findings)
      verified <- sum(pathway_findings$Quote.Verified)
      eligible_count <- sum(pathway_findings$Eligible.For.Summary)
      quote_quarantined <- total - verified
      disease_excluded <- verified - eligible_count
      grounding_status <- if (!total) {
        "no_findings"
      } else if (verified == total) {
        "complete"
      } else if (verified > 0L) {
        "partial"
      } else {
        "none"
      }
      eligibility_status <- if (!total) {
        "no_findings"
      } else if (eligible_count == total) {
        "complete"
      } else if (eligible_count > 0L) {
        "partial"
      } else {
        "none"
      }
      grounding_counts[i, c(
        "Total.Findings", "Quote.Verified", "Eligible.For.Summary",
        "Quote.Quarantined", "Disease.Excluded"
      )] <- c(
        total, verified, eligible_count, quote_quarantined, disease_excluded
      )
      grounding_counts$Status[i] <- grounding_status
      grounding_counts$Eligibility.Status[i] <- eligibility_status
      grounding_counts$Display[i] <- if (!total) {
        "Source grounding: no findings"
      } else {
        sprintf(
          paste0(
            "Source grounding: %d/%d quotes verified; ",
            "quote-quarantined: %d; disease-excluded: %d"
          ),
          verified, total, quote_quarantined, disease_excluded
        )
      }
      if ("Quote.Grounding" %in% names(rows) &&
          !identical(rows$Quote.Grounding[i], grounding_status)) {
        fail(paste(
          "Quote.Grounding disagrees with strict finding-level",
          "Quote.Verified flags."
        ))
      }
      if ("Evidence.Eligibility" %in% names(rows) &&
          !identical(rows$Evidence.Eligibility[i], eligibility_status)) {
        fail(paste(
          "Evidence.Eligibility disagrees with finding-level",
          "Eligible.For.Summary flags."
        ))
      }
    } else {
      if ("Quote.Grounding" %in% names(rows)) {
        fail("Quote.Grounding requires finding-level Quote.Verified flags.")
      }
      if ("Evidence.Eligibility" %in% names(rows)) {
        total <- nrow(pathway_findings)
        eligible_count <- sum(pathway_findings$Eligible.For.Summary)
        eligibility_status <- if (!total) {
          "no_findings"
        } else if (eligible_count == total) {
          "complete"
        } else if (eligible_count > 0L) {
          "partial"
        } else {
          "none"
        }
        if (!identical(rows$Evidence.Eligibility[i], eligibility_status)) {
          fail(paste(
            "Evidence.Eligibility disagrees with finding-level",
            "Eligible.For.Summary flags."
          ))
        }
      }
    }
    for (name in names(sets)) {
      actual <- unique(as.character(eligible$PMID[eligible[[name]]]))
      if (!setequal(sets[[name]], actual)) {
        fail("Exploration channel PMID sets disagree with eligible findings: ", name)
      }
    }
    if (!setequal(keys, unique(as.character(eligible$PMID)))) {
      fail("Key.PMIDs must match the eligible exploration findings.")
    }
    all_ids <- unique(c(unlist(sets, use.names = FALSE), keys))
    if (length(all_ids) && any(
      !paste(rows$Pathway[i], all_ids, sep = "\r") %in% record_keys
    )) {
      fail("An exploration PMID is absent from retrieved evidence for its pathway.")
    }
    if (rows$Exploration.Status[i] %in% c("not_assessed", "input_failed")) {
      if (length(all_ids)) {
        fail("Unassessed exploration cannot contain channel or key PMIDs.")
      }
      references[i] <- if (rows$Exploration.Status[i] == "input_failed") {
        "Key evidence: input/search failed"
      } else {
        "Key evidence: not assessed"
      }
    } else {
      if (!nzchar(trimws(rows$Biological.Summary[i]))) {
        fail("An assessed pathway requires a meaningful Biological.Summary.")
      }
      for (name in names(sets)) {
        channel_counts[[name]][i] <- length(sets[[name]])
      }
      channel_counts$Unique[i] <- length(unique(
        as.character(eligible$PMID[!eligible$Review])
      ))
      if (!is.numeric(rows$Unique.PMID.Count) ||
          !is.finite(rows$Unique.PMID.Count[i]) ||
          rows$Unique.PMID.Count[i] != channel_counts$Unique[i]) {
        fail(paste(
          "Unique.PMID.Count disagrees with unique eligible non-review",
          "PMID records."
        ))
      }
      references[i] <- paste0(
        if (length(keys)) {
          paste("Key PMID:", paste(keys, collapse = "; "))
        } else {
          "Key PMID: none in the assessed selection"
        },
        "  |  Non-review sources: ", channel_counts$Unique[i],
        if (channel_counts$Review[i] > 0) {
          paste0("; reviews: ", channel_counts$Review[i], " (background)")
        } else {
          ""
        }
      )
    }
    if (rows$Exploration.Status[i] == "provided" && !nrow(eligible)) {
      fail("Provided exploration requires at least one eligible finding.")
    }
    if (rows$Exploration.Status[i] != "provided" && nrow(eligible)) {
      fail("Non-provided exploration cannot contain eligible findings.")
    }
    for (name in colnames(channel_scopes)) {
      present <- unique(as.character(eligible$Scope[eligible[[name]]]))
      ordered <- names(scope_labels)[names(scope_labels) %in% present]
      channel_scopes[i, name] <- paste(scope_labels[ordered], collapse = "/")
    }
    channel_sets[[i]] <- sets
  }

  unassessed <- rows$Exploration.Status == "not_assessed"
  input_failed <- rows$Exploration.Status == "input_failed"
  no_relevant <- rows$Exploration.Status == "no_relevant_findings"
  ungrounded <- rows$Exploration.Status == "ungrounded_findings"
  if (identical(text_mode, "full")) {
    rows$Summary.Display <- rows$Biological.Summary
  } else {
    publications <- channel_counts$Unique
    rows$Summary.Display <- sprintf(
      paste0(
        "Source-linked evidence: %d unique non-review publication%s; ",
        "see the accompanying report for the full synthesis."
      ),
      publications, ifelse(publications == 1L, "", "s")
    )
    rows$Summary.Display[no_relevant] <- paste(
      "Source-linked evidence: no eligible findings in the assessed selection;",
      "see the accompanying report."
    )
    rows$Summary.Display[ungrounded] <- paste(
      "Source-linked evidence: no eligible grounded findings in the assessed",
      "selection; see the accompanying report."
    )
  }
  rows$Summary.Display[unassessed] <-
    "Tumor-biology exploration has not been assessed."
  rows$Summary.Display[input_failed] <- paste(
    "The literature input/search failed; biological exploration is not",
    "assessable."
  )
  process_labels <- c(
    proliferation = "Proliferation",
    invasion_metastasis = "Invasion / metastasis",
    treatment_response = "Treatment response", cell_death = "Cell death",
    immune_microenvironment = "Immune microenvironment",
    metabolic_adaptation = "Metabolic adaptation",
    disease_development = "Disease development", angiogenesis = "Angiogenesis",
    stromal_microenvironment = "Stromal crosstalk", other = "Other process"
  )
  display_processes <- vapply(rows$Biological.Processes, function(value) {
    terms <- trimws(strsplit(value, ";", fixed = TRUE)[[1L]])
    terms <- terms[nzchar(terms)]
    known <- terms %in% names(process_labels)
    terms[known] <- unname(process_labels[terms[known]])
    paste(terms, collapse = "; ")
  }, character(1L))
  rows$Processes.Display <- ifelse(
    nzchar(trimws(rows$Biological.Processes)),
    paste("Processes:", display_processes), "Processes: not assigned"
  )
  rows$Processes.Display[
    rows$Exploration.Status %in% c("not_assessed", "input_failed")
  ] <- "Processes: not assessed"
  plausibility_labels <- c(
    strongly_plausible = "Strongly plausible",
    plausible = "Plausible",
    mixed_context_dependent = "Mixed / context-dependent",
    limited_evidence = "Limited evidence",
    not_assessable = "Not assessable"
  )
  rows$Plausibility.Display <- paste0(
    "Biological plausibility (AI synthesis): ",
    unname(plausibility_labels[rows$Biological.Plausibility])
  )
  rows$Context.Display <- if (identical(text_mode, "full")) {
    ifelse(
      nzchar(trimws(rows$Context.Dependence)),
      paste("Context:", rows$Context.Dependence), ""
    )
  } else {
    rep("", nrow(rows))
  }
  rows$Context.Display[unassessed | input_failed] <- ""
  rows$Key.Evidence.Display <- references

  sources <- unique(rows$Assessment.Source[nzchar(rows$Assessment.Source)])
  if (!length(sources)) sources <- "unspecified"
  provider <- report$interpretation$provider
  provider_name <- if (is.list(provider)) scalar(provider$name) else scalar(provider)
  provider_model <- if (is.list(provider)) scalar(provider$model) else "unknown"
  disease <- scalar(report$analysis$disease)
  endpoint <- scalar(report$analysis$endpoint)
  time_unit <- scalar(report$analysis$time_unit)
  evidence_provenance <- report$evidence$provenance
  if (!is.list(evidence_provenance)) evidence_provenance <- list()
  evidence_database <- scalar(evidence_provenance$database)
  source_kind <- scalar(evidence_provenance$source_kind)
  synthetic_fixture <- tolower(source_kind) %in% c(
    "synthetic", "synthetic_fixture", "fixture"
  ) || tolower(evidence_database) %in% c(
    "synthetic", "synthetic_not_pubmed", "synthetic_fixture"
  )
  synthetic_banner <- if (synthetic_fixture) {
    paste(
      "INTERMETA PATHWAY SCORES | SYNTHETIC AI-EVIDENCE FIXTURE |",
      "NO PUBMED / MODEL EVIDENCE"
    )
  } else {
    ""
  }
  if (synthetic_fixture) {
    rows$Summary.Display <- sprintf(
      paste0(
        "Synthetic AI-evidence fixture: %d generated test record%s; ",
        "no literature evidence was used."
      ),
      channel_counts$Unique,
      ifelse(channel_counts$Unique == 1L, "", "s")
    )
    rows$Summary.Display[unassessed] <- "Synthetic exploration was not assessed."
    rows$Summary.Display[input_failed] <- "Synthetic exploration input failed."
    references <- gsub(
      "Key PMID:", "Generated record IDs:", references, fixed = TRUE
    )
    references <- gsub(
      "Non-review sources:", "Generated records:", references, fixed = TRUE
    )
    grounding_counts$Display <- sub(
      "^Source grounding:", "Fixture grounding:", grounding_counts$Display
    )
    rows$Key.Evidence.Display <- references
  }

  displayed_text <- enc2utf8(c(
    rows$Pathway, rows$Plausibility.Display, rows$Summary.Display,
    rows$Processes.Display,
    rows$Context.Display, references, grounding_counts$Display, disease,
    synthetic_banner
  ))
  non_ascii <- any(vapply(
    displayed_text,
    function(x) any(utf8ToInt(x) > 127L), logical(1L)
  ))
  drawing_font <- if (is.null(font_family)) "sans" else trimws(font_family)
  font_checked <- !non_ascii && is.null(font_family)
  verify_font <- function() {
    if (!capabilities("cairo")) {
      fail(paste(
        "Non-ASCII/Chinese rendering or custom fonts require R Cairo support.",
        "No output was created."
      ))
    }
    if (!requireNamespace("systemfonts", quietly = TRUE)) {
      fail(paste(
        "Non-ASCII text or an explicit font_family requires optional",
        "systemfonts and a local font with complete glyph coverage.",
        "Nothing is downloaded. No output was created."
      ))
    }
    families <- unique(systemfonts::system_fonts()$family)
    candidates <- if (!is.null(font_family)) {
      exact <- families[tolower(families) == tolower(trimws(font_family))]
      if (!length(exact)) {
        fail(
          "Requested font_family is not installed: ", font_family,
          ". Choose an installed family from systemfonts::system_fonts(). ",
          "No output was created."
        )
      }
      exact[1L]
    } else {
      preferred <- c(
        "Noto Sans CJK SC", "Noto Sans SC", "Source Han Sans SC",
        "Microsoft YaHei", "Microsoft YaHei UI", "WenQuanYi Micro Hei",
        "SimHei", "SimSun"
      )
      found <- match(tolower(preferred), tolower(families))
      unique(c(
        families[found[!is.na(found)]],
        families[grepl("CJK|Han Sans|Hei|Noto Sans", families, ignore.case = TRUE)]
      ))
    }
    glyphs <- unique(unlist(strsplit(
      paste(displayed_text, collapse = ""), "", fixed = TRUE
    )))
    glyphs <- glyphs[!grepl("[[:space:]]", glyphs)]
    for (family in candidates) {
      okay <- tryCatch(all(vapply(c(FALSE, TRUE), function(bold) {
        info <- systemfonts::glyph_info(glyphs, family = family, bold = bold)
        nrow(info) == length(glyphs) &&
          all(!is.na(info$index) & info$index > 0)
      }, logical(1L))), error = function(e) FALSE)
      if (okay) return(family)
    }
    fail(paste(
      "No installed candidate font covers every displayed character in",
      "regular and bold. Install a suitable local font and supply",
      "font_family. Missing glyphs are not silently replaced.",
      "No output was created."
    ))
  }

  palette <- c(
    Functional = "#13B9C9", Omics = "#8550D1", Patient = "#FFBD3D",
    missing = "#E9EEF4", score = "#F46A57"
  )
  grounding_palette <- c(
    complete = "#16856B", partial = "#D97706", none = "#C2415D",
    no_findings = "#66758A", unavailable = "#66758A"
  )
  plausibility_palette <- c(
    strongly_plausible = "#00A6A6", plausible = "#63B700",
    mixed_context_dependent = "#F58A07", limited_evidence = "#E83E73",
    not_assessable = "#66758A"
  )
  plausibility_wrapped <- vapply(
    rows$Plausibility.Display, wrap, "", width = 137L
  )
  summary_wrapped <- vapply(rows$Summary.Display, wrap, "", width = 129L)
  context_wrapped <- vapply(rows$Context.Display, wrap, "", width = 139L)
  process_wrapped <- vapply(rows$Processes.Display, wrap, "", width = 137L)
  reference_wrapped <- vapply(references, wrap, "", width = 144L)
  grounding_wrapped <- vapply(
    grounding_counts$Display, wrap, "", width = 144L
  )
  pathway_wrapped <- vapply(rows$Pathway, wrap, "", width = 32L)
  line_count <- function(x) {
    ifelse(nzchar(x), lengths(strsplit(x, "\n", fixed = TRUE)), 0L)
  }
  row_mm <- pmax(
    28.5,
    15 + line_count(summary_wrapped) * 3.1 +
      line_count(context_wrapped) * 2.85 +
      line_count(plausibility_wrapped) * 3.1 +
      line_count(process_wrapped) * 2.9 +
      line_count(reference_wrapped) * 2.7 +
      line_count(grounding_wrapped) * 2.7
  )
  width_mm <- 190
  height_mm <- max(122, sum(row_mm) + 46)
  max_png_pixels <- 30000
  if (height_mm / 25.4 * 300 > max_png_pixels) {
    fail(paste(
      "The requested figure is too tall for safe 300-dpi rendering.",
      "Use compact_english or shorten the validated narrative."
    ))
  }
  notices <- if (height_mm > 210) {
    paste(
      "Full exploration text increases the figure height beyond 210 mm.",
      "Nothing is silently truncated; use a concise validated synthesis for",
      "a compact main figure."
    )
  } else {
    character()
  }
  grounding_caption <- if (quote_grounding_available) {
    paste0(
      "Audit labels separate quote-quarantined rows (all findings minus strictly verified quotes) from disease-excluded rows (verified quotes minus evidence-eligible findings). ",
      "Both are excluded from channel counts, Key PMIDs and source-linked public synthesis. ",
      "When any pathway finding is ineligible, model-level biological summary, context and hypothesis are retained only in Model.* audit fields; for partial eligibility, the public summary uses only eligible Source.Linked.Observations. ",
      "Disease-excluded counts use model-supplied Study.Disease and Disease.Match values plus lexical target-label consistency; they are not an independent semantic disease classification from the article. ",
      "Pathway/PMID membership and exact quote grounding do not semantically prove that a model Observation or claim is supported by the excerpt. "
    )
  } else {
    ""
  }
  text_mode_caption <- if (identical(text_mode, "compact_english")) {
    paste0(
      "Compact English figure mode retains structured pathway, process, source-count, scope, PMID and grounding fields; ",
      "the complete source-linked synthesis, context and hypotheses remain in the accompanying report. "
    )
  } else {
    ""
  }
  caption <- paste0(
    "Tumor-biology exploration of the Top ", nrow(rows), " ", direction,
    "-scoring metabolic pathways. Disease: ", disease,
    ". Original analysis context only: endpoint ", endpoint,
    "; time unit ", time_unit,
    ". Pathway order and signed Difference are unchanged. Difference is not pathway activity, metabolic flux, or a hazard ratio. ",
    "Functional, omics and patient channels are non-exclusive study-design groupings, not evidence grades, positive votes or prediction-validation results. ",
    "Counts use unique PMID records within each channel; repeated PMIDs contribute once to the overall unique non-review source count, which also includes other eligible non-review study types. ",
    "Zero means no records assigned to that channel in the assessed selection, not proof that no such literature exists; NA/not assessed remains distinct. ",
    "Reviews are background sources, not independent functional experiments. Component-, metabolite-, enzyme- or gene-level studies may inform biology without testing an entire pathway. ",
    "The patient channel can include tissue, observational or interventional studies; its presence is not evidence of a prospective clinical trial. ",
    "Scope markers are descriptive: P = pathway, C = component (including genes, enzymes or metabolites), ? = unclear. ",
    "The displayed biological-plausibility label is a model-assisted synthesis of the screened, quote-verified selection; it is not a causal grade, a pathway-activity measurement, independent replication, or validation that the InterMeta score direction is correct. ",
    "Biological processes and public synthesis are copied from the structurally validated report's eligible, source-linked fields; this does not independently verify the semantic accuracy of a model/assistant Observation or claim. Null and differing-context findings are not removed when eligible. ",
    grounding_caption, text_mode_caption,
    "Key references are selected illustrations, not the complete evidence inventory. Research hypotheses remain in the report and are not established conclusions. ",
    "Interpretation provider: ", provider_name, "; model: ", provider_model,
    ". Assessment source: ", paste(sources, collapse = "; "),
    ". The figure renderer makes no model call or new inference. Human review and independent experiments remain necessary."
  )
  if (synthetic_fixture) {
    caption <- paste(
      "SYNTHETIC SOFTWARE FIXTURE ONLY.",
      "The displayed pathway order and signed Difference are copied from the supplied InterMeta report.",
      "No PubMed retrieval or AI model call is represented by the evidence panel.",
      "All record identifiers, literature-like fields, classifications and summaries in that panel are generated test fixtures; they are not publications or biological evidence.",
      "Signed Difference is not pathway activity, metabolic flux, a hazard ratio or a biological truth score.",
      paste0("Fixture provider: ", provider_name, "; model label: ", provider_model, "."),
      "The fixture validates deterministic report structure and rendering only."
    )
  }

  draw <- function() {
    if (!font_checked) {
      fail("Font coverage has not been checked. Render through the helper with render = TRUE.")
    }
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      width = grid::unit(width_mm, "mm"),
      height = grid::unit(height_mm, "mm"),
      xscale = c(0, width_mm), yscale = c(0, height_mm)
    ))
    on.exit(grid::popViewport(), add = TRUE)
    ink <- "#172B43"
    muted <- "#52647B"
    txt <- function(
        label, x, y, size = 7.7, col = ink, just = "left", bold = FALSE) {
      grid::grid.text(
        label, x = grid::unit(x, "native"), y = grid::unit(y, "native"),
        just = just,
        gp = grid::gpar(
          fontfamily = drawing_font, fontsize = size, col = col,
          fontface = if (bold) "bold" else "plain", lineheight = 1.1
        )
      )
    }
    rect <- function(x, y, width, height, fill) {
      grid::grid.rect(
        x = grid::unit(x, "native"), y = grid::unit(y, "native"),
        width = grid::unit(width, "native"),
        height = grid::unit(height, "native"),
        gp = grid::gpar(fill = fill, col = NA)
      )
    }
    line <- function(x, y, col = "#D7E1EC") {
      grid::grid.lines(
        x = grid::unit(x, "native"), y = grid::unit(y, "native"),
        gp = grid::gpar(col = col, lwd = 0.6)
      )
    }
    txt(
      paste0(
        if (synthetic_fixture) "INTERMETA PATHWAY SCORES: " else "",
        "Top ", nrow(rows), " ", direction, " metabolic pathways"
      ),
      5, height_mm - 7, size = 12.7, bold = TRUE
    )
    txt(
      if (synthetic_fixture) synthetic_banner else paste0(
        disease, "  |  Literature-derived biological context  |  Original ",
        direction, " score selection"
      ),
      5, height_mm - 13, size = 7.5,
      col = if (synthetic_fixture) "#C2415D" else muted,
      bold = synthetic_fixture
    )
    hy <- height_mm - 22
    rect(95, hy, 180, 11, "#F0F5FB")
    txt("Metabolic pathway", 7, hy, size = 8.1, bold = TRUE)
    txt("Signed Difference", 78, hy, size = 8, just = "centre", bold = TRUE)
    txt(
      if (synthetic_fixture) {
        "Synthetic fixture channels - generated record IDs"
      } else {
        "Evidence channels - unique PMID records"
      },
      143, hy, size = 8, just = "centre", bold = TRUE
    )
    edges <- c(height_mm - 29, height_mm - 29 - cumsum(row_mm))
    score_range <- range(c(0, rows$Difference))
    if (diff(score_range) == 0) score_range <- c(-1, 1)
    score_range <- score_range + c(-1, 1) * diff(score_range) * 0.04
    score_x <- function(x) {
      63 + (x - score_range[1]) / diff(score_range) * 30
    }
    for (i in seq_len(nrow(rows))) {
      y <- edges[i] - 6
      rect(
        95, mean(edges[c(i, i + 1L)]), 180, row_mm[i] - 0.5,
        if (i %% 2L) "#FFFFFF" else "#F8FAFD"
      )
      txt(pathway_wrapped[i], 7, y, size = 8.1, bold = TRUE)
      line(c(score_x(0), score_x(0)), y + c(-3.5, 3.1), "#A9B8C8")
      rect(
        mean(c(score_x(0), score_x(rows$Difference[i]))), y - 2.0,
        max(abs(score_x(rows$Difference[i]) - score_x(0)), 0.35),
        2.5, palette["score"]
      )
      txt(
        sprintf("%+.2f", rows$Difference[i]), 78, y + 2.4,
        size = 8.3, just = "centre", bold = TRUE
      )
      channel_names <- c("Functional", "Omics", "Patient")
      for (j in seq_along(channel_names)) {
        channel <- channel_names[j]
        count <- channel_counts[[channel]][i]
        center <- c(113, 141.5, 170)[j]
        fill <- if (!is.na(count) && count > 0) {
          palette[channel]
        } else {
          palette["missing"]
        }
        rect(center, y, 25.5, 9.5, fill)
        label <- if (is.na(count)) {
          paste0(
            channel, "\n",
            if (rows$Exploration.Status[i] == "input_failed") {
              "Input failed"
            } else {
              "Not assessed"
            }
          )
        } else {
          paste0(
            channel, " ", count, "\n",
            if (nzchar(channel_scopes[i, channel])) {
              paste("Scope:", channel_scopes[i, channel])
            } else {
              "None assigned"
            }
          )
        }
        txt(
          label, center, y, size = 7.5, just = "centre", bold = TRUE,
          col = if (channel == "Omics" && !is.na(count) && count > 0) {
            "#FFFFFF"
          } else {
            ink
          }
        )
      }
      top <- edges[i] - 12.4
      txt(
        plausibility_wrapped[i], 7,
        top - line_count(plausibility_wrapped[i]) * 1.55,
        size = 7.5,
        col = unname(plausibility_palette[rows$Biological.Plausibility[i]]),
        bold = TRUE
      )
      top <- top - line_count(plausibility_wrapped[i]) * 3.1 - 0.4
      txt(
        process_wrapped[i], 7,
        top - line_count(process_wrapped[i]) * 1.45,
        size = 7.3, col = "#136B82", bold = TRUE
      )
      top <- top - line_count(process_wrapped[i]) * 2.9 - 0.6
      txt(
        summary_wrapped[i], 7,
        top - line_count(summary_wrapped[i]) * 1.55, size = 7.7
      )
      top <- top - line_count(summary_wrapped[i]) * 3.1
      if (nzchar(context_wrapped[i])) {
        txt(
          context_wrapped[i], 7,
          top - line_count(context_wrapped[i]) * 1.425,
          size = 7.15, col = muted
        )
        top <- top - line_count(context_wrapped[i]) * 2.85
      }
      txt(
        reference_wrapped[i], 7,
        top - line_count(reference_wrapped[i]) * 1.35,
        size = 7, col = "#435B79"
      )
      top <- top - line_count(reference_wrapped[i]) * 2.7
      if (nzchar(grounding_wrapped[i])) {
        txt(
          grounding_wrapped[i], 7,
          top - line_count(grounding_wrapped[i]) * 1.35,
          size = 7,
          col = unname(grounding_palette[
            grounding_counts$Eligibility.Status[i]
          ]),
          bold = TRUE
        )
      }
      line(c(5, 185), rep(edges[i + 1L], 2))
    }
    footer <- edges[length(edges)] - 4
    txt(
      if (synthetic_fixture) {
        paste(
          "AI narrative, channels and record IDs are synthetic fixtures,",
          "not literature evidence."
        )
      } else {
        paste(
          "Channels overlap and describe study designs, not evidence quality,",
          "positive findings or prognostic validation."
        )
      },
      5, footer, size = 7.1, bold = TRUE,
      col = if (synthetic_fixture) "#C2415D" else ink
    )
    txt(
      paste(
        "Scope: P = pathway; C = component / gene / metabolite; ? = unclear.",
        "Zero = none assigned. Reviews are background."
      ),
      5, footer - 4, size = 7, col = muted
    )
    txt(
      paste(
        "Difference is not pathway activity or flux. Summaries remain",
        "exploratory; source-level audit and hypotheses accompany the report."
      ),
      5, footer - 8, size = 7, col = muted
    )
  }

  output <- list(pdf = NULL, png = NULL)
  if (render) {
    if (is.null(pdf_file) && is.null(png_file)) {
      fail("Supply a new pdf_file or png_file, or use render = FALSE.")
    }
    if (!is.null(pdf_file)) {
      output$pdf <- output_path(pdf_file, "pdf_file", "pdf")
    }
    if (!is.null(png_file)) {
      output$png <- output_path(png_file, "png_file", "png")
    }
    if (!font_checked) {
      drawing_font <- verify_font()
      font_checked <- TRUE
    }

    staged <- list()
    committed <- character()
    on.exit({
      remaining <- unlist(staged, use.names = FALSE)
      if (length(remaining)) unlink(remaining[file.exists(remaining)], force = TRUE)
    }, add = TRUE)
    write_plot <- function(target, type) {
      stage <- tempfile(
        pattern = paste0(".", basename(target), "."),
        tmpdir = dirname(target), fileext = paste0(".", type)
      )
      keep_stage <- FALSE
      on.exit({
        if (!keep_stage && file.exists(stage)) unlink(stage, force = TRUE)
      }, add = TRUE)
      if (identical(type, "pdf")) {
        if (capabilities("cairo")) {
          grDevices::cairo_pdf(
            stage, width = width_mm / 25.4, height = height_mm / 25.4,
            family = drawing_font, onefile = TRUE
          )
        } else {
          grDevices::pdf(
            stage, width = width_mm / 25.4, height = height_mm / 25.4,
            family = "Helvetica", useDingbats = FALSE
          )
        }
      } else if (non_ascii || !is.null(font_family)) {
        grDevices::png(
          stage, width = width_mm, height = height_mm, units = "mm", res = 300,
          bg = "white", type = "cairo", family = drawing_font
        )
      } else {
        grDevices::png(
          stage, width = width_mm, height = height_mm, units = "mm", res = 300,
          bg = "white"
        )
      }
      device <- grDevices::dev.cur()
      on.exit({
        open_devices <- grDevices::dev.list()
        if (!is.null(open_devices) && device %in% open_devices) {
          grDevices::dev.off(device)
        }
      }, add = TRUE)
      draw()
      grDevices::dev.off(device)
      if (!file.exists(stage) || is.na(file.info(stage)$size) ||
          file.info(stage)$size <= 0) {
        fail("The ", toupper(type), " graphics device produced no usable output.")
      }
      keep_stage <- TRUE
      stage
    }
    if (!is.null(output$pdf)) staged$pdf <- write_plot(output$pdf, "pdf")
    if (!is.null(output$png)) staged$png <- write_plot(output$png, "png")

    for (type in names(staged)) {
      target <- output[[type]]
      if (file.exists(target) || dir.exists(target) ||
          !file.rename(staged[[type]], target)) {
        if (length(committed)) unlink(committed, force = TRUE)
        fail("Could not atomically publish the ", toupper(type), " output: ", target)
      }
      committed <- c(committed, target)
      staged[[type]] <- NULL
    }
  }

  invisible(list(
    pathways = rows,
    channel_counts = channel_counts,
    channel_sets = channel_sets,
    channel_scopes = channel_scopes,
    grounding_counts = grounding_counts,
    analysis = report$analysis,
    assessment_sources = sources,
    provider = provider,
    caption = caption,
    notices = notices,
    palette = palette,
    grounding_palette = grounding_palette,
    plausibility_palette = plausibility_palette,
    text_mode = text_mode,
    dimensions_mm = c(width = width_mm, height = height_mm),
    font = list(
      requested = font_family, selected = drawing_font,
      non_ascii = non_ascii,
      coverage_required = non_ascii || !is.null(font_family),
      coverage_checked = font_checked &&
        (non_ascii || !is.null(font_family))
    ),
    evidence_source = list(
      source_kind = source_kind,
      database = evidence_database,
      synthetic_fixture = synthetic_fixture
    ),
    output = output,
    draw = draw
  ))
}
