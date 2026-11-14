# Deterministic, privacy-preserving context for literature-backed reporting.

.InterMetaAIScalarCharacter <- function(x, name, allow_null = FALSE){
  if (allow_null && is.null(x)){
    return("unknown")
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))){
    stop('"', name, '" must be one non-empty character value.')
  }
  trimws(x)
}

.InterMetaAIHash <- function(x){
  if (is.null(x)){
    return("unknown")
  }
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

.InterMetaAINormalizeGroup <- function(x){
  value <- tolower(trimws(as.character(x)))
  value[value %in% c("f", "favorable")] <- "F"
  value[value %in% c("u", "unfavorable")] <- "U"
  value[value %in% c("n", "neutral")] <- "N"
  value
}

.InterMetaAINormalizeScore <- function(score){
  if (!is.data.frame(score) || nrow(score) == 0L){
    stop('"score" must be a non-empty data frame.')
  }
  required <- c("Pathway", "Difference", "Group")
  missing <- setdiff(required, names(score))
  if (length(missing) > 0L){
    stop("The score data are missing required column(s): ",
         paste(missing, collapse = ", "), ".")
  }
  pathway <- trimws(as.character(score$Pathway))
  if (anyNA(pathway) || any(!nzchar(pathway))){
    stop('"score$Pathway" must not contain missing or empty values.')
  }
  if (anyDuplicated(pathway)){
    stop('"score$Pathway" must contain unique pathway names.')
  }
  if (!is.numeric(score$Difference) || any(!is.finite(score$Difference))){
    stop('"score$Difference" must contain only finite numeric values.')
  }
  group <- .InterMetaAINormalizeGroup(score$Group)
  if (anyNA(group) || any(!group %in% c("F", "U", "N"))){
    stop('"score$Group" must contain only F, U, N, favorable, ',
         'unfavorable, or neutral.')
  }
  inconsistent <- (group == "F" & score$Difference <= 0) |
    (group == "U" & score$Difference >= 0) |
    (group == "N" & score$Difference != 0)
  if (any(inconsistent)){
    stop('"score$Group" is inconsistent with the sign of "score$Difference".')
  }

  safe_names <- intersect(
    c("Index", "Rank", "Signed.Rank", "Pathway", "OutDegree",
      "Favorable.Score", "InDegree", "Unfavorable.Score", "Difference",
      "Group"),
    names(score)
  )
  safe <- score[, safe_names, drop = FALSE]
  safe$Pathway <- pathway
  safe$Difference <- as.numeric(score$Difference)
  safe$Group <- group
  row.names(safe) <- NULL
  safe
}

.InterMetaAINormalizeEnrichment <- function(enrichment){
  if (is.null(enrichment)){
    return(NULL)
  }
  if (!is.data.frame(enrichment) || nrow(enrichment) == 0L){
    stop('"enrichment" must be NULL or a non-empty data frame.')
  }
  required <- c("Favorable.Pathway", "Unfavorable.Pathway")
  missing <- setdiff(required, names(enrichment))
  if (length(missing) > 0L){
    stop("The enrichment data are missing required column(s): ",
         paste(missing, collapse = ", "), ".")
  }
  for (column in required){
    value <- trimws(as.character(enrichment[[column]]))
    if (anyNA(value) || any(!nzchar(value))){
      stop('"enrichment$', column,
           '" must not contain missing or empty values.')
    }
    enrichment[[column]] <- value
  }
  if (!"Pathway.Interaction" %in% names(enrichment)){
    enrichment$Pathway.Interaction <- paste(
      enrichment$Favorable.Pathway,
      enrichment$Unfavorable.Pathway,
      sep = "__"
    )
  }
  enrichment$Pathway.Interaction <- trimws(
    as.character(enrichment$Pathway.Interaction)
  )
  if (anyNA(enrichment$Pathway.Interaction) ||
      any(!nzchar(enrichment$Pathway.Interaction))){
    stop('"enrichment$Pathway.Interaction" must not contain missing or ',
         "empty values.")
  }
  if ("Shared.MGP" %in% names(enrichment)){
    shared <- as.character(enrichment$Shared.MGP)
    shared[is.na(shared)] <- ""
    enrichment$Shared.MGP <- shared
  }
  numeric_columns <- intersect(
    c("Index", "Rank", "No.Total.MGP", "No.Shared.MGP", "Enrichment.Ratio",
      "P.Value", "Adjust.P.Value", "Permutation.FDR"),
    names(enrichment)
  )
  for (column in numeric_columns){
    if (!is.numeric(enrichment[[column]])){
      stop('"enrichment$', column, '" must be numeric.')
    }
  }
  safe_names <- intersect(
    c("Index", "Rank", "Pathway.Interaction", "Favorable.Pathway",
      "Unfavorable.Pathway", "No.Total.MGP", "No.Shared.MGP",
      "Enrichment.Ratio", "P.Value", "Adjust.P.Value", "Permutation.FDR",
      "Shared.MGP"),
    names(enrichment)
  )
  safe <- enrichment[, safe_names, drop = FALSE]
  row.names(safe) <- NULL
  safe
}

.InterMetaAINormalizeSynonyms <- function(pathway_synonyms){
  if (is.null(pathway_synonyms)){
    return(list())
  }
  if (is.data.frame(pathway_synonyms)){
    if (!all(c("Pathway", "Synonym") %in% names(pathway_synonyms))){
      stop('A data-frame "pathway_synonyms" must contain Pathway and ',
           "Synonym columns.")
    }
    pathway <- trimws(as.character(pathway_synonyms$Pathway))
    synonym <- trimws(as.character(pathway_synonyms$Synonym))
    if (anyNA(pathway) || anyNA(synonym) || any(!nzchar(pathway)) ||
        any(!nzchar(synonym))){
      stop('"pathway_synonyms" must not contain missing or empty values.')
    }
    result <- split(synonym, pathway)
  } else {
    if (!is.list(pathway_synonyms) && !is.character(pathway_synonyms)){
      stop('"pathway_synonyms" must be a named list, named character ',
           "vector, or data frame.")
    }
    if (is.null(names(pathway_synonyms)) ||
        anyNA(names(pathway_synonyms)) ||
        any(!nzchar(trimws(names(pathway_synonyms))))){
      stop('"pathway_synonyms" must have non-empty pathway names.')
    }
    result <- lapply(pathway_synonyms, function(x){
      if (!is.character(x) || anyNA(x) || any(!nzchar(trimws(x)))){
        stop('Each "pathway_synonyms" entry must contain non-empty ',
             "character values.")
      }
      trimws(x)
    })
    names(result) <- trimws(names(pathway_synonyms))
  }
  if (anyDuplicated(names(result))){
    combined <- split(unlist(result, use.names = FALSE),
                      rep(names(result), lengths(result)))
    result <- combined
  }
  result <- lapply(result, function(x) sort(unique(x)))
  result[order(names(result))]
}

.InterMetaAIParameter <- function(parameters, names){
  for (name in names){
    value <- parameters[[name]]
    if (!is.null(value) && length(value) == 1L && !is.na(value)){
      if (is.character(value)){
        value <- trimws(value)
        if (!nzchar(value)){
          next
        }
      }
      return(value)
    }
  }
  "unknown"
}

.InterMetaAIStabilitySummary <- function(stability){
  if (is.null(stability)){
    return(NULL)
  }
  summary <- if (is.data.frame(stability)) {
    stability
  } else if (is.list(stability) && is.data.frame(stability$summary)) {
    stability$summary
  } else {
    stop('"stability" must be an InterMetaStability object, a list with a ',
         'data-frame "summary", or a summary data frame.')
  }
  if (!"Pathway" %in% names(summary)){
    stop('The stability summary is missing required column: Pathway.')
  }
  pathway <- trimws(as.character(summary$Pathway))
  if (anyNA(pathway) || any(!nzchar(pathway)) || anyDuplicated(pathway)){
    stop('"stability$summary$Pathway" must contain unique, non-empty names.')
  }
  safe_names <- intersect(
    c("Pathway", "Selection.Count", "Selection.Frequency",
      "Favorable.Frequency", "Unfavorable.Frequency", "Neutral.Frequency",
      "Dominant.Direction", "Direction.Consistency", "Mean.Difference",
      "SD.Difference", "CV.Difference", "Mean.Rank", "SD.Rank",
      "Mean.Magnitude.Rank", "SD.Magnitude.Rank",
      grep(
        "^Top[0-9]+\\.(Favorable\\.|Unfavorable\\.)?Frequency$",
        names(summary), value = TRUE
      ),
      "Stability.Class"),
    names(summary)
  )
  safe <- summary[, safe_names, drop = FALSE]
  safe$Pathway <- pathway
  row.names(safe) <- NULL
  safe
}

.InterMetaAIStabilityRuns <- function(stability, summary){
  if (is.null(stability)){
    return(NA_integer_)
  }
  if (is.list(stability) && is.data.frame(stability$run_status) &&
      "Status" %in% names(stability$run_status)){
    return(sum(tolower(as.character(stability$run_status$Status)) == "success",
               na.rm = TRUE))
  }
  if (is.list(stability) && !is.null(stability$runs)){
    if (!is.list(stability$runs)){
      stop('"stability$runs" must be a list.')
    }
    return(sum(vapply(stability$runs, function(x){
      is.list(x) && (is.null(x$status) || identical(x$status, "success"))
    }, logical(1))))
  }
  if (is.list(stability) && is.list(stability$parameters)){
    count <- stability$parameters$n_resample
    if (is.numeric(count) && length(count) == 1L && is.finite(count) &&
        count >= 0 && count == floor(count)){
      return(as.integer(count))
    }
  }
  if (!is.null(summary) &&
      all(c("Selection.Count", "Selection.Frequency") %in% names(summary))){
    usable <- is.finite(summary$Selection.Count) &
      is.finite(summary$Selection.Frequency) & summary$Selection.Frequency > 0
    estimates <- summary$Selection.Count[usable] /
      summary$Selection.Frequency[usable]
    if (length(estimates) > 0L && all(is.finite(estimates))){
      return(as.integer(round(max(estimates))))
    }
  }
  NA_integer_
}

.InterMetaAISelectPathways <- function(score, direction, top_n){
  selected <- list()
  if (direction %in% c("both", "unfavorable")){
    unfavorable <- score[score$Group == "U", , drop = FALSE]
    unfavorable <- unfavorable[
      order(unfavorable$Difference, unfavorable$Pathway,
            method = "radix"),
      , drop = FALSE
    ]
    unfavorable <- head(unfavorable, top_n)
    unfavorable$Direction <- rep("unfavorable", nrow(unfavorable))
    unfavorable$Direction.Rank <- seq_len(nrow(unfavorable))
    selected[[length(selected) + 1L]] <- unfavorable
  }
  if (direction %in% c("both", "favorable")){
    favorable <- score[score$Group == "F", , drop = FALSE]
    favorable <- favorable[
      order(-favorable$Difference, favorable$Pathway,
            method = "radix"),
      , drop = FALSE
    ]
    favorable <- head(favorable, top_n)
    favorable$Direction <- rep("favorable", nrow(favorable))
    favorable$Direction.Rank <- seq_len(nrow(favorable))
    selected[[length(selected) + 1L]] <- favorable
  }
  result <- do.call(rbind, selected)
  row.names(result) <- NULL
  result
}

.InterMetaAIRelevantInteractions <- function(enrichment, pathways,
                                             score_condition, score_cutoff){
  if (is.null(enrichment) || nrow(pathways) == 0L){
    return(NULL)
  }
  favorable <- pathways$Pathway[pathways$Direction == "favorable"]
  unfavorable <- pathways$Pathway[pathways$Direction == "unfavorable"]
  keep <- enrichment$Favorable.Pathway %in% favorable |
    enrichment$Unfavorable.Pathway %in% unfavorable
  if (is.character(score_condition) && length(score_condition) == 1L &&
      !identical(score_condition, "unknown") &&
      score_condition %in% names(enrichment) &&
      is.numeric(score_cutoff) && length(score_cutoff) == 1L &&
      is.finite(score_cutoff)){
    keep <- keep & is.finite(enrichment[[score_condition]]) &
      enrichment[[score_condition]] < score_cutoff
  }
  result <- enrichment[keep, , drop = FALSE]
  if (nrow(result) > 0L){
    key <- paste(result$Pathway.Interaction,
                 if ("Index" %in% names(result)) result$Index else seq_len(nrow(result)),
                 sep = "\r")
    result <- result[order(key, method = "radix"), , drop = FALSE]
  }
  row.names(result) <- NULL
  result
}

.InterMetaAIDriverPairs <- function(interactions, pathways){
  empty <- data.frame(
    Pathway = character(),
    Direction = character(),
    Gene.Pair = character(),
    Source.Gene = character(),
    Target.Gene = character(),
    Pathway.Interaction = character(),
    stringsAsFactors = FALSE
  )
  if (is.null(interactions) || nrow(interactions) == 0L ||
      !"Shared.MGP" %in% names(interactions)){
    return(empty)
  }
  rows <- vector("list", nrow(interactions) * 2L)
  cursor <- 0L
  for (i in seq_len(nrow(interactions))){
    pairs <- trimws(strsplit(interactions$Shared.MGP[i], "/", fixed = TRUE)[[1]])
    pairs <- sort(unique(pairs[nzchar(pairs)]))
    if (length(pairs) == 0L){
      next
    }
    parsed <- .InterMetaParseGenePairs(pairs)
    mapping <- list(
      favorable = interactions$Favorable.Pathway[i],
      unfavorable = interactions$Unfavorable.Pathway[i]
    )
    for (entry_direction in names(mapping)){
      pathway <- mapping[[entry_direction]]
      if (!any(pathways$Pathway == pathway &
               pathways$Direction == entry_direction)){
        next
      }
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        Pathway = pathway,
        Direction = entry_direction,
        Gene.Pair = parsed$Pair,
        Source.Gene = parsed$Source.Gene,
        Target.Gene = parsed$Target.Gene,
        Pathway.Interaction = interactions$Pathway.Interaction[i],
        stringsAsFactors = FALSE
      )
    }
  }
  if (cursor == 0L){
    return(empty)
  }
  result <- unique(do.call(rbind, rows[seq_len(cursor)]))
  result <- result[
    order(result$Direction, result$Pathway, result$Gene.Pair,
          result$Pathway.Interaction, method = "radix"),
    , drop = FALSE
  ]
  row.names(result) <- NULL
  result
}

#' Create a deterministic, privacy-preserving AI analysis context
#'
#' Selects the strongest favorable and/or unfavorable pathways and constructs
#' a compact context suitable for literature retrieval and language-model
#' interpretation. Patient-level expression, survival records, sample IDs, and
#' resampling run payloads are never copied into the returned object.
#'
#' @param enrichment An optional InterMeta enrichment data frame.
#' @param score An optional InterMeta pathway-score data frame.
#' @param stability An optional `InterMetaStability` object or stability summary.
#' @param analysis An optional full-sample result list containing `enrich`,
#'   `score`, and `parameters`. Non-NULL direct arguments take precedence.
#' @param disease One non-empty disease name.
#' @param endpoint Survival or clinical endpoint; use `"unknown"` if not known.
#' @param time_unit Time unit used for the endpoint.
#' @param cohort Optional scalar cohort label, never a vector of sample IDs.
#' @param top_n Positive integer giving the number selected per direction.
#' @param direction One of `"both"`, `"unfavorable"`, or `"favorable"`.
#' @param pathway_synonyms Optional named character vector, named list, or data
#'   frame with `Pathway` and `Synonym` columns.
#'
#' @return An object of class `InterMetaAIContext` with deterministic analysis,
#'   pathway, interaction, driver-gene-pair, warning, and provenance fields.
#' @export
InterMetaAIContext <- function(enrichment = NULL,
                               score = NULL,
                               stability = NULL,
                               analysis = NULL,
                               disease,
                               endpoint = "unknown",
                               time_unit = "unknown",
                               cohort = NULL,
                               top_n = 5L,
                               direction = c("both", "unfavorable", "favorable"),
                               pathway_synonyms = NULL){
  disease <- .InterMetaAIScalarCharacter(disease, "disease")
  endpoint <- .InterMetaAIScalarCharacter(endpoint, "endpoint")
  time_unit <- .InterMetaAIScalarCharacter(time_unit, "time_unit")
  cohort <- .InterMetaAIScalarCharacter(cohort, "cohort", allow_null = TRUE)
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
      top_n < 1 || top_n != floor(top_n)){
    stop('"top_n" must be one positive integer.')
  }
  top_n <- as.integer(top_n)
  direction <- match.arg(direction)
  if (!is.null(analysis) && (!is.list(analysis) || is.data.frame(analysis))){
    stop('"analysis" must be NULL or a full-sample result list.')
  }
  if (!is.null(analysis) && is.null(names(analysis))){
    stop('"analysis" must be a named list.')
  }

  if (is.null(enrichment) && !is.null(analysis)){
    enrichment <- if (!is.null(analysis$enrich)) {
      analysis$enrich
    } else {
      analysis$enrichment
    }
  }
  if (is.null(score) && !is.null(analysis)){
    score <- analysis$score
  }
  if (is.null(stability) && !is.null(analysis)){
    stability <- analysis$stability
  }
  parameters <- if (!is.null(analysis) && is.list(analysis$parameters)) {
    analysis$parameters
  } else {
    list()
  }

  score_safe <- .InterMetaAINormalizeScore(score)
  enrichment_safe <- .InterMetaAINormalizeEnrichment(enrichment)
  synonyms <- .InterMetaAINormalizeSynonyms(pathway_synonyms)
  stability_summary <- .InterMetaAIStabilitySummary(stability)
  stability_runs <- .InterMetaAIStabilityRuns(stability, stability_summary)

  score_condition <- .InterMetaAIParameter(
    parameters, c("score_condition", "condition")
  )
  score_cutoff <- .InterMetaAIParameter(
    parameters, c("score_cutoff", "cutoff")
  )
  pair_condition <- .InterMetaAIParameter(parameters, "pair_condition")
  pair_cutoff <- .InterMetaAIParameter(parameters, "pair_cutoff")
  analysis_profile <- .InterMetaAIParameter(
    parameters, c("analysis_profile", "profile")
  )

  pathways <- .InterMetaAISelectPathways(score_safe, direction, top_n)
  pathways$Synonyms <- vapply(pathways$Pathway, function(pathway){
    value <- synonyms[[pathway]]
    if (is.null(value)) "" else paste(value, collapse = " | ")
  }, character(1))

  interactions <- .InterMetaAIRelevantInteractions(
    enrichment_safe, pathways, score_condition, score_cutoff
  )
  driver_gene_pairs <- .InterMetaAIDriverPairs(interactions, pathways)
  pathways$Driver.Genes <- vapply(seq_len(nrow(pathways)), function(i){
    subset <- driver_gene_pairs[
      driver_gene_pairs$Pathway == pathways$Pathway[i] &
        driver_gene_pairs$Direction == pathways$Direction[i],
      , drop = FALSE
    ]
    genes <- if (pathways$Direction[i] == "favorable") {
      subset$Source.Gene
    } else {
      subset$Target.Gene
    }
    paste(sort(unique(genes[nzchar(genes)])), collapse = ",")
  }, character(1))
  # Keep a syntactic-name alias for consumers that normalize punctuation.
  pathways$DriverGenes <- pathways$Driver.Genes

  if (!is.null(stability_summary)){
    selected_stability <- stability_summary[
      match(pathways$Pathway, stability_summary$Pathway),
      setdiff(names(stability_summary), "Pathway"),
      drop = FALSE
    ]
    pathways <- cbind(pathways, selected_stability)
    row.names(pathways) <- NULL
  }

  warnings <- character()
  if (tolower(endpoint) == "unknown"){
    warnings <- c(warnings,
      "The clinical endpoint is unknown; biological direction cannot be interpreted precisely.")
  }
  if (tolower(time_unit) == "unknown"){
    warnings <- c(warnings,
      "The endpoint time unit is unknown and must be confirmed before reporting.")
  }
  if (identical(toupper(as.character(pair_condition)), "PV")){
    warnings <- c(warnings,
      "Gene pairs were selected with unadjusted P-values (pair_condition = PV); treat findings as exploratory.")
  }
  if (identical(pair_condition, "unknown")){
    warnings <- c(warnings,
      "The gene-pair selection condition is unknown.")
  }
  if (identical(score_condition, "unknown") || identical(score_cutoff, "unknown")){
    warnings <- c(warnings,
      "The pathway-score condition or cutoff is unknown; related enrichment rows were not significance-filtered.")
  } else if (!is.null(enrichment_safe) &&
             !score_condition %in% names(enrichment_safe)){
    warnings <- c(warnings,
      paste0("The enrichment table lacks score condition column ",
             score_condition, "; related rows were not significance-filtered."))
  }
  if (is.null(stability)){
    warnings <- c(warnings,
      "Stability results were not supplied; pathway ranking robustness is unassessed.")
  } else if (is.na(stability_runs)){
    warnings <- c(warnings,
      "The number of stability resamples could not be determined.")
  } else {
    warnings <- c(warnings,
      paste0("Only ", stability_runs,
             " prespecified internal resamples were available; interpret this as an internal consistency check, not as external validation or a population-level uncertainty estimate."))
  }
  requested_groups <- switch(
    direction,
    both = c(U = "unfavorable", F = "favorable"),
    unfavorable = c(U = "unfavorable"),
    favorable = c(F = "favorable")
  )
  for (group in names(requested_groups)){
    available <- sum(score_safe$Group == group)
    if (available < top_n){
      warnings <- c(warnings,
        paste0("Only ", available, " ", requested_groups[[group]],
               " pathways were available; fewer than top_n = ", top_n,
               " were returned."))
    }
  }
  if (!is.null(interactions) && nrow(interactions) > 0L &&
      !"Shared.MGP" %in% names(interactions)){
    warnings <- c(warnings,
      "Shared.MGP is absent from enrichment results; driver gene pairs could not be extracted.")
  }
  warnings <- unique(warnings)

  analysis_context <- list(
    disease = disease,
    endpoint = endpoint,
    time_unit = time_unit,
    cohort = cohort,
    direction = direction,
    top_n = top_n,
    score_condition = score_condition,
    score_cutoff = score_cutoff,
    pair_condition = pair_condition,
    pair_cutoff = pair_cutoff,
    analysis_profile = analysis_profile,
    stability_resamples = if (is.na(stability_runs)) "unknown" else stability_runs
  )
  stability_hash_input <- if (is.null(stability_summary)) {
    NULL
  } else {
    list(summary = stability_summary, resamples = stability_runs)
  }
  provenance <- list(
    context_schema_version = "intermeta-ai-context-v1",
    enrichment_hash = .InterMetaAIHash(enrichment_safe),
    score_hash = .InterMetaAIHash(score_safe),
    stability_hash = .InterMetaAIHash(stability_hash_input),
    pathway_synonyms_hash = .InterMetaAIHash(
      if (length(synonyms) == 0L) NULL else synonyms
    )
  )
  provenance$input_hash <- .InterMetaAIHash(list(
    analysis = analysis_context,
    pathways = pathways,
    interactions = interactions,
    driver_gene_pairs = driver_gene_pairs,
    warnings = warnings
  ))

  result <- list(
    analysis = analysis_context,
    pathways = pathways,
    interactions = if (is.null(interactions)) {
      data.frame(
        Pathway.Interaction = character(),
        Favorable.Pathway = character(),
        Unfavorable.Pathway = character(),
        stringsAsFactors = FALSE
      )
    } else {
      interactions
    },
    driver_gene_pairs = driver_gene_pairs,
    warnings = warnings,
    provenance = provenance
  )
  class(result) <- c("InterMetaAIContext", "list")
  result
}
