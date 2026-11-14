# Portions of this file were adapted from TimiGP 1.3.0:
# https://github.com/CSkylarL/TimiGP
# The upstream source identifies Chenyang Skylar Li as author.
# Modified for InterMeta; changes through 2026-11-03.
# SPDX-License-Identifier: GPL-3.0-only
# See inst/NOTICE and COPYING for provenance and license details.

## InterMeta pathway-level wrappers and metabolic interaction utilities.

.InterMetaPAdjustMethods <- c("holm", "hochberg", "hommel", "bonferroni",
                              "BH", "BY", "fdr", "none")

.InterMetaValidatePAdjust <- function(method){
  if (length(method) != 1 || !method %in% .InterMetaPAdjustMethods){
    stop('"p.adj" must be one of stats::p.adjust.methods.')
  }
  method
}

.InterMetaAnalysisProfile <- function(analysis_profile){
  if (!is.character(analysis_profile) || length(analysis_profile) != 1L ||
      is.na(analysis_profile) || analysis_profile != "TimiGP-legacy") {
    stop('Only analysis_profile = "TimiGP-legacy" is supported. ',
         'Omit analysis_profile or explicitly use "TimiGP-legacy". ',
         'To reproduce an archived "InterMeta"-profile run, use its original ',
         'package version; results are not converted automatically.',
         call. = FALSE)
  }
  "TimiGP-legacy"
}

.InterMetaLegacyGenesetAttribute <- "InterMeta.legacy.permutation.rows"

.InterMetaNormalizeGenesetRows <- function(pathwayset,
                                            dataset = NULL,
                                            deduplicate = TRUE){
  if (is.null(pathwayset) || !is.data.frame(pathwayset)){
    stop('"pathwayset" must be a data frame.')
  }
  ## Data-frame subsetting preserves arbitrary attributes.  Strip the hidden
  ## legacy slots here so that a stale attribute can never influence ordinary
  ## validation or de-duplication.
  attr(pathwayset, .InterMetaLegacyGenesetAttribute) <- NULL
  names(pathwayset) <- trimws(names(pathwayset))
  if (!"Pathway" %in% names(pathwayset)){
    if ("CellType" %in% names(pathwayset)){
      names(pathwayset)[names(pathwayset) == "CellType"] <- "Pathway"
    } else {
      stop('The geneset must contain a "Pathway" or "CellType" column.')
    }
  }
  if (!"Gene" %in% names(pathwayset)){
    stop('The geneset must contain a "Gene" column.')
  }
  if (!is.null(dataset)){
    if (!"Dataset" %in% names(pathwayset)){
      stop('Filtering by "dataset" requires a "Dataset" column.')
    }
    pathwayset <- pathwayset[pathwayset$Dataset %in% dataset, , drop = FALSE]
    if (nrow(pathwayset) == 0){
      stop("No rows matched the requested dataset.")
    }
  }
  pathwayset$Pathway <- trimws(as.character(pathwayset$Pathway))
  pathwayset$Gene <- trimws(as.character(pathwayset$Gene))
  keep <- !is.na(pathwayset$Pathway) & !is.na(pathwayset$Gene) &
    nzchar(pathwayset$Pathway) & nzchar(pathwayset$Gene)
  pathwayset <- pathwayset[keep, , drop = FALSE]
  if (deduplicate){
    pathwayset <- unique(pathwayset)
  }
  if (nrow(pathwayset) == 0){
    stop("The geneset has no valid pathway-gene rows.")
  }
  if (any(grepl("_", pathwayset$Gene, fixed = TRUE))){
    stop("Gene identifiers containing underscores are not supported because ",
         "gene-pair identifiers use an underscore separator.")
  }
  pathwayset
}

.InterMetaNormalizeGeneset <- function(pathwayset,
                                       dataset = NULL){
  .InterMetaNormalizeGenesetRows(
    pathwayset = pathwayset,
    dataset = dataset,
    deduplicate = TRUE
  )
}

.InterMetaGenesetSlotKey <- function(pathwayset,
                                     columns){
  values <- lapply(columns, function(column){
    value <- as.character(pathwayset[[column]])
    paste0(nchar(value, type = "bytes"), ":", value)
  })
  do.call(paste, c(values, sep = "|"))
}

.InterMetaLegacyPermutationGeneset <- function(pathwayset,
                                                normalized = NULL){
  input_rows <- .InterMetaNormalizeGenesetRows(
    pathwayset,
    deduplicate = FALSE
  )
  if (is.null(normalized)){
    normalized <- .InterMetaNormalizeGeneset(pathwayset)
  }

  legacy_rows <- attr(
    pathwayset,
    .InterMetaLegacyGenesetAttribute,
    exact = TRUE
  )
  if (is.null(legacy_rows)){
    return(input_rows)
  }
  legacy_rows <- tryCatch(
    .InterMetaNormalizeGenesetRows(legacy_rows, deduplicate = FALSE),
    error = function(e) NULL
  )
  if (is.null(legacy_rows)){
    return(input_rows)
  }

  ## Prefer Dataset when it is present on both sides, but permit users to
  ## select only Pathway/Gene columns.  Length-prefixed keys avoid delimiter
  ## collisions in custom pathway names.
  identity_columns <- c("Pathway", "Gene")
  if ("Dataset" %in% names(input_rows) &&
      "Dataset" %in% names(legacy_rows)){
    identity_columns <- c(identity_columns, "Dataset")
  }
  input_key <- .InterMetaGenesetSlotKey(input_rows, identity_columns)
  normalized_key <- .InterMetaGenesetSlotKey(normalized, identity_columns)
  legacy_key <- .InterMetaGenesetSlotKey(legacy_rows, identity_columns)

  ## Explicit duplicate rows in a modified/custom data frame are authoritative.
  ## The loader itself exposes a unique table, so this does not suppress its
  ## hidden raw multiplicity.
  if (anyDuplicated(input_key)){
    return(input_rows)
  }
  requested_key <- unique(normalized_key)
  if (!all(requested_key %in% legacy_key)){
    return(input_rows)
  }
  selected <- legacy_rows[legacy_key %in% requested_key, , drop = FALSE]
  selected_key <- .InterMetaGenesetSlotKey(selected, identity_columns)
  if (!setequal(unique(selected_key), requested_key)){
    return(input_rows)
  }
  row.names(selected) <- NULL
  selected
}

.InterMetaEmptyEnrichment <- function(include_permutation = FALSE){
  out <- data.frame(
    Index = integer(),
    Rank = numeric(),
    Pathway.Interaction = character(),
    Favorable.Pathway = character(),
    Unfavorable.Pathway = character(),
    No.Total.MGP = integer(),
    No.Shared.MGP = integer(),
    Enrichment.Ratio = numeric(),
    P.Value = numeric(),
    Adjust.P.Value = numeric(),
    stringsAsFactors = FALSE
  )
  if (include_permutation){
    out$Permutation.FDR <- numeric()
  }
  out$Shared.MGP <- character()
  out$Total.MGP <- character()
  out
}

.InterMetaEmptyScore <- function(){
  data.frame(
    Index = integer(),
    Rank = numeric(),
    Signed.Rank = numeric(),
    Pathway = character(),
    OutDegree = integer(),
    Favorable.Score = numeric(),
    InDegree = integer(),
    Unfavorable.Score = numeric(),
    Difference = numeric(),
    Group = character(),
    stringsAsFactors = FALSE
  )
}

.InterMetaParseGenePairs <- function(pairs){
  pairs <- unique(as.character(pairs))
  pairs <- pairs[!is.na(pairs) & nzchar(pairs)]
  if (length(pairs) == 0){
    return(data.frame(
      Pair = character(),
      Source.Gene = character(),
      Target.Gene = character(),
      stringsAsFactors = FALSE
    ))
  }
  separator <- regexpr("_", pairs, fixed = TRUE)
  invalid <- separator <= 1L | separator >= nchar(pairs)
  remainder <- substring(pairs, separator + 1L)
  multiple <- regexpr("_", remainder, fixed = TRUE) > 0L
  if (any(invalid | multiple)){
    stop("Every gene-pair identifier must contain exactly one underscore. ",
         "Gene identifiers containing underscores are not supported.")
  }
  data.frame(
    Pair = pairs,
    Source.Gene = substring(pairs, 1L, separator - 1L),
    Target.Gene = remainder,
    stringsAsFactors = FALSE
  )
}

InterMetaAvailableGenesets <- function(){
  c("glucose", "lipid", "amino", "all", "kegg", "reactome")
}

InterMetaLoadGeneset <- function(dataset = "all",
                                 file = NULL){
  bundled <- is.null(file)
  if (bundled){
    dataset <- match.arg(dataset, InterMetaAvailableGenesets())
    file <- system.file("extdata", "metabolic-genesets",
                        paste0("geneset-", dataset, ".txt"),
                        package = "InterMeta")
    if (!nzchar(file)){
      stop("Cannot find the bundled metabolic geneset file.")
    }
  }
  if (!file.exists(file)){
    stop("Geneset file does not exist: ", file)
  }

  geneset <- utils::read.delim(file, stringsAsFactors = FALSE,
                               check.names = FALSE)
  keep_column <- vapply(geneset, function(x){
    any(!is.na(x) & nzchar(trimws(as.character(x))))
  }, logical(1))
  geneset <- geneset[, keep_column, drop = FALSE]
  geneset <- geneset[, nzchar(trimws(names(geneset))), drop = FALSE]
  if (!"Dataset" %in% names(geneset)){
    geneset$Dataset <- if (bundled) dataset else "custom"
  }
  raw_geneset <- .InterMetaNormalizeGenesetRows(
    geneset,
    deduplicate = FALSE
  )
  raw_geneset <- raw_geneset[, c("Pathway", "Gene", "Dataset"), drop = FALSE]
  row.names(raw_geneset) <- NULL
  geneset <- unique(raw_geneset)
  row.names(geneset) <- NULL
  attr(geneset, .InterMetaLegacyGenesetAttribute) <- raw_geneset
  geneset
}

InterMetaCheckEvent <- function(info = NULL,
                                event_col = 1,
                                time_col = 2,
                                min_events = 1,
                                min_censored = 1){
  if (is.null(info) || !is.data.frame(info) || ncol(info) < 2){
    stop('"info" must be a data frame with event and time columns.')
  }
  if (is.null(row.names(info)) || anyDuplicated(row.names(info))){
    stop("Clinical sample IDs must be unique row names.")
  }
  resolve_column <- function(column, label){
    if (is.character(column)){
      if (length(column) != 1 || !column %in% names(info)){
        stop("Unknown ", label, " column: ", paste(column, collapse = ", "))
      }
      match(column, names(info))
    } else {
      column <- as.integer(column)
      if (length(column) != 1 || is.na(column) ||
          column < 1 || column > ncol(info)){
        stop("Invalid ", label, " column index.")
      }
      column
    }
  }
  event_col <- resolve_column(event_col, "event")
  time_col <- resolve_column(time_col, "time")
  if (event_col == time_col){
    stop("Event and time must use different columns.")
  }

  event <- suppressWarnings(as.numeric(as.character(info[[event_col]])))
  time <- suppressWarnings(as.numeric(as.character(info[[time_col]])))
  invalid_event <- !is.na(event) & !event %in% c(0, 1)
  if (any(invalid_event)){
    stop("Event values must be coded as 0 (censored) and 1 (event).")
  }
  keep <- is.finite(event) & is.finite(time) & time > 0
  removed <- sum(!keep)
  info <- info[keep, , drop = FALSE]
  event <- event[keep]
  time <- time[keep]
  info[[event_col]] <- event
  info[[time_col]] <- time
  if (sum(event == 1) < min_events){
    stop("Clinical data contain fewer than ", min_events, " events.")
  }
  if (sum(event == 0) < min_censored){
    stop("Clinical data contain fewer than ", min_censored,
         " censored observations.")
  }
  message(removed, " individuals were filtered out due to missing/invalid ",
          "event data or time <= 0.")
  info
}

InterMetaPreprocess <- function(pathwayset = NULL,
                                marker = NULL,
                                rna = NULL,
                                cohort = NULL,
                                log = TRUE,
                                GMNorm = TRUE,
                                low_expression_filter = TRUE,
                                analysis_profile = "TimiGP-legacy"){
  analysis_profile <- .InterMetaAnalysisProfile(analysis_profile)
  if (is.null(rna) || !(is.matrix(rna) || is.data.frame(rna))){
    stop('"rna" must be a numeric matrix or data frame.')
  }
  if (is.null(row.names(rna)) || anyDuplicated(row.names(rna))){
    stop("RNA gene identifiers must be unique row names.")
  }
  if (is.null(colnames(rna)) || anyDuplicated(colnames(rna))){
    stop("RNA sample identifiers must be unique column names.")
  }
  if (!is.logical(log) || length(log) != 1 || is.na(log)){
    stop('"log" must be TRUE or FALSE.')
  }
  if (!is.logical(GMNorm) || length(GMNorm) != 1 || is.na(GMNorm)){
    stop('"GMNorm" must be TRUE or FALSE.')
  }
  if (!is.logical(low_expression_filter) ||
      length(low_expression_filter) != 1 ||
      is.na(low_expression_filter)){
    stop('"low_expression_filter" must be TRUE or FALSE.')
  }

  if (is.null(marker)){
    if (is.null(pathwayset)){
      stop('Either "pathwayset" or "marker" is required.')
    }
    pathwayset <- .InterMetaNormalizeGeneset(pathwayset)
    marker <- unique(pathwayset$Gene)
  }
  marker <- unique(trimws(as.character(marker)))
  marker <- marker[!is.na(marker) & nzchar(marker)]

  if (!is.null(cohort)){
    cohort <- unique(as.character(cohort))
    missing_samples <- setdiff(cohort, colnames(rna))
    if (length(missing_samples) > 0){
      warning(length(missing_samples),
              " requested samples were absent from the RNA matrix.")
    }
    cohort <- cohort[cohort %in% colnames(rna)]
    if (length(cohort) < 5){
      stop("Fewer than 5 requested samples were found in the RNA matrix.")
    }
    legacy_cohort <- colnames(rna)[colnames(rna) %in% cohort]
    rna <- rna[, legacy_cohort, drop = FALSE]
  }

  matched_marker <- intersect(row.names(rna), marker)
  if (length(matched_marker) == 0){
    stop("No geneset genes were found in the RNA matrix.")
  }
  rna <- rna[matched_marker, , drop = FALSE]
  original <- rna
  rna <- suppressWarnings(matrix(
    as.numeric(as.matrix(rna)),
    nrow = nrow(rna),
    dimnames = dimnames(rna)
  ))
  if (anyNA(rna) && !anyNA(as.matrix(original))){
    stop("RNA values must be numeric.")
  }
  if (any(!is.finite(rna))){
    stop("RNA values must be finite and cannot contain NA, NaN, or Inf.")
  }
  if (low_expression_filter){
    expressed <- rowSums(rna > 0) >= ncol(rna) * 0.5
    if (!any(expressed)){
      stop("All genes were removed by the TimiGP legacy low-expression filter.")
    }
    message(sum(!expressed),
            " genes with low expression were filtered out.")
    rna <- rna[expressed, , drop = FALSE]
  }
  if (log){
    log_safe <- apply(rna, 1, function(x) all(x > -1))
    if (!any(log_safe)){
      stop("All genes would produce non-finite values under the TimiGP ",
           "legacy log10(x + 1) transformation.")
    }
    if (any(!log_safe)){
      message(sum(!log_safe),
              " genes that would produce non-finite legacy log values ",
              "were filtered out safely.")
      rna <- rna[log_safe, , drop = FALSE]
    }
    rna <- log10(rna + 1)
  }
  if (GMNorm){
    rna <- sweep(rna, 1, apply(rna, 1, stats::median), FUN = "-")
  }
  rna[order(row.names(rna)), , drop = FALSE]
}

InterMetaGenePair <- function(rna = NULL,
                              both = FALSE,
                              cont = FALSE,
                              chunk_size = 10000,
                              analysis_profile = "TimiGP-legacy"){
  analysis_profile <- .InterMetaAnalysisProfile(analysis_profile)
  if (is.null(rna) || !(is.matrix(rna) || is.data.frame(rna))){
    stop('"rna" must be a numeric matrix or data frame.')
  }
  if (is.null(row.names(rna)) || is.null(colnames(rna))){
    stop('"rna" must have gene row names and sample column names.')
  }
  if (anyDuplicated(row.names(rna))){
    stop("RNA gene row names must be unique.")
  }
  if (!is.logical(both) || length(both) != 1 || is.na(both) ||
      !is.logical(cont) || length(cont) != 1 || is.na(cont)){
    stop('"both" and "cont" must be TRUE or FALSE.')
  }
  chunk_size <- as.integer(chunk_size)
  if (is.na(chunk_size) || chunk_size < 1){
    stop('"chunk_size" must be a positive integer.')
  }
  rna <- as.matrix(rna)
  storage.mode(rna) <- "double"
  if (any(!is.finite(rna))){
    stop("RNA values must be finite.")
  }
  if (nrow(rna) < 2){
    stop("At least two genes are required to generate gene pairs.")
  }
  if (cont && all(apply(rna, 1, stats::median) == 0)){
    stop("Continuous gene-pair scores require GMNorm = FALSE.")
  }

  order_gene <- order(row.names(rna))
  rna <- rna[order_gene, , drop = FALSE]
  pair_index <- which(
    lower.tri(matrix(FALSE, nrow = nrow(rna), ncol = nrow(rna))),
    arr.ind = TRUE
  )
  pair_count <- nrow(pair_index)
  keep <- logical(pair_count)
  starts <- seq.int(1L, pair_count, by = chunk_size)
  threshold <- as.integer(round(ncol(rna) * 0.1))
  for (start in starts){
    index <- start:min(start + chunk_size - 1L, pair_count)
    difference <- rna[pair_index[index, 1], , drop = FALSE] -
      rna[pair_index[index, 2], , drop = FALSE]
    if (cont){
      keep[index] <- apply(difference, 1, function(x){
        diff(range(x)) > 0.1
      })
    } else {
      positive <- rowSums(difference > 0)
      keep[index] <- positive >= threshold &
        (ncol(rna) - positive) >= threshold
    }
  }
  pair_index <- pair_index[keep, , drop = FALSE]
  output_count <- nrow(pair_index)
  message(pair_count - output_count, " gene pairs were filtered out.")
  if (output_count == 0){
    return(matrix(
      if (cont) numeric() else logical(),
      nrow = 0,
      ncol = ncol(rna),
      dimnames = list(character(), colnames(rna))
    ))
  }

  output <- matrix(
    if (cont) 0 else FALSE,
    nrow = output_count,
    ncol = ncol(rna),
    dimnames = list(NULL, colnames(rna))
  )
  starts <- seq.int(1L, output_count, by = chunk_size)
  for (start in starts){
    index <- start:min(start + chunk_size - 1L, output_count)
    difference <- rna[pair_index[index, 1], , drop = FALSE] -
      rna[pair_index[index, 2], , drop = FALSE]
    output[index, ] <- if (cont) difference else difference > 0
  }
  row.names(output) <- paste(
    row.names(rna)[pair_index[, 1]],
    row.names(rna)[pair_index[, 2]],
    sep = "_"
  )
  if (both){
    reverse <- if (cont) -output else !output
    reverse_names <- paste(
      row.names(rna)[pair_index[, 2]],
      row.names(rna)[pair_index[, 1]],
      sep = "_"
    )
    row.names(reverse) <- reverse_names
    output <- rbind(output, reverse)
  }
  message(nrow(output), " gene pairs were produced.")
  output
}

InterMetaCox <- function(mps = NULL,
                         info = NULL,
                         p.adj = "BH",
                         parallel = FALSE,
                         core = 1,
                         keep_mps = TRUE,
                         analysis_profile = "TimiGP-legacy"){
  analysis_profile <- .InterMetaAnalysisProfile(analysis_profile)
  info <- InterMetaCheckEvent(info)
  TimiCOX(mps = mps, info = info, p.adj = p.adj,
          parallel = parallel, core = core, keep_mps = keep_mps,
          analysis_profile = analysis_profile)
}

InterMetaBG <- function(marker.pair = NULL){
  if (is.null(marker.pair)){
    stop('"marker.pair" is required.')
  }
  parsed <- .InterMetaParseGenePairs(marker.pair)
  reverse_pair <- paste(parsed$Target.Gene, parsed$Source.Gene, sep = "_")
  background <- unique(c(parsed$Pair, reverse_pair))
  message(length(background), " background pairs were generated from ",
          nrow(parsed), " input pairs.")
  background
}

InterMetaPathwayPair <- function(pathwayset = NULL,
                                 dataset = NULL,
                                 core = 1){
  pathwayset <- .InterMetaNormalizeGeneset(pathwayset, dataset = dataset)
  pathwayset <- unique(pathwayset[, c("Pathway", "Gene"), drop = FALSE])
  pathways <- unique(pathwayset$Pathway)
  if (length(pathways) < 2){
    stop("At least two pathways are required.")
  }
  genes <- unique(pathwayset$Gene)
  membership <- data.frame(
    Pathway.ID = match(pathwayset$Pathway, pathways),
    Gene.ID = match(pathwayset$Gene, genes),
    stringsAsFactors = FALSE
  )
  pathway_genes <- split(pathwayset$Gene,
                         factor(pathwayset$Pathway, levels = pathways))
  out <- list(
    pathways = pathways,
    genes = genes,
    pathway_genes = pathway_genes,
    membership = membership,
    geneset_hash = digest::digest(pathwayset, algo = "xxhash64"),
    n_pathway_interactions = length(pathways) * (length(pathways) - 1L)
  )
  class(out) <- "InterMetaPathwayPair"
  out
}

print.InterMetaPathwayPair <- function(x, ...){
  cat("<InterMetaPathwayPair>\n")
  cat("  Pathways:", length(x$pathways), "\n")
  cat("  Genes:", length(x$genes), "\n")
  cat("  Directed pathway interactions:", x$n_pathway_interactions, "\n")
  invisible(x)
}

.InterMetaPreparePairContext <- function(gene,
                                         background){
  background <- .InterMetaParseGenePairs(background)
  query <- unique(as.character(gene))
  query <- query[!is.na(query) & nzchar(query)]
  background$Is.Query <- background$Pair %in% query
  background
}

.InterMetaPairMembers <- function(pair_context,
                                  source_genes,
                                  target_genes){
  keep <- pair_context$Source.Gene %in% source_genes &
    pair_context$Target.Gene %in% target_genes
  pair_context$Pair[keep]
}

.InterMetaEnrichCompact <- function(pathway.pair,
                                    pair_context,
                                    p.adj,
                                    keep_gene_pairs){
  pathways <- pathway.pair$pathways
  genes <- pathway.pair$genes
  n_pathway <- length(pathways)
  n_gene <- length(genes)
  membership <- matrix(0, nrow = n_pathway, ncol = n_gene)
  membership[cbind(pathway.pair$membership$Pathway.ID,
                   pathway.pair$membership$Gene.ID)] <- 1

  source_id <- match(pair_context$Source.Gene, genes)
  target_id <- match(pair_context$Target.Gene, genes)
  valid <- !is.na(source_id) & !is.na(target_id) &
    source_id != target_id
  path_count <- colSums(membership)
  sole_path <- rep(NA_integer_, n_gene)
  single_gene <- which(path_count == 1)
  if (length(single_gene) > 0){
    single_rows <- pathway.pair$membership[
      pathway.pair$membership$Gene.ID %in% single_gene, , drop = FALSE
    ]
    sole_path[single_rows$Gene.ID] <- single_rows$Pathway.ID
  }
  both_single_same <- rep(FALSE, length(valid))
  candidate <- which(valid & path_count[source_id] == 1 &
                       path_count[target_id] == 1)
  if (length(candidate) > 0){
    both_single_same[candidate] <-
      sole_path[source_id[candidate]] == sole_path[target_id[candidate]]
  }
  valid <- valid & !both_single_same
  pair_context <- pair_context[valid, , drop = FALSE]
  source_id <- source_id[valid]
  target_id <- target_id[valid]

  N <- nrow(pair_context)
  n <- sum(pair_context$Is.Query)
  adjacency <- matrix(0, nrow = n_gene, ncol = n_gene)
  if (N > 0){
    adjacency[cbind(source_id, target_id)] <- 1
  }
  total_count <- (membership %*% adjacency) %*% t(membership)
  adjacency[] <- 0
  if (n > 0){
    query_rows <- which(pair_context$Is.Query)
    adjacency[cbind(source_id[query_rows], target_id[query_rows])] <- 1
  }
  shared_count <- (membership %*% adjacency) %*% t(membership)
  rm(adjacency)

  grid <- expand.grid(
    Source.ID = seq_len(n_pathway),
    Target.ID = seq_len(n_pathway),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- grid[grid$Source.ID != grid$Target.ID, , drop = FALSE]
  matrix_index <- cbind(grid$Source.ID, grid$Target.ID)
  M <- as.integer(total_count[matrix_index])
  k <- as.integer(shared_count[matrix_index])
  ratio <- rep(NA_real_, length(M))
  testable <- M > 0 & n > 0 & N > 0
  ratio[testable] <- k[testable] * N / M[testable] / n
  p_value <- rep(1, length(M))
  p_value[testable] <- stats::phyper(
    k[testable] - 1,
    M[testable],
    N - M[testable],
    n,
    lower.tail = FALSE
  )

  source_pathway <- pathways[grid$Source.ID]
  target_pathway <- pathways[grid$Target.ID]
  shared_members <- rep("", nrow(grid))
  total_members <- rep(NA_character_, nrow(grid))
  if (keep_gene_pairs %in% c("shared", "all") && n > 0){
    query_context <- pair_context[pair_context$Is.Query, , drop = FALSE]
    for (i in which(k > 0)){
      shared_members[i] <- paste(.InterMetaPairMembers(
        query_context,
        pathway.pair$pathway_genes[[grid$Source.ID[i]]],
        pathway.pair$pathway_genes[[grid$Target.ID[i]]]
      ), collapse = "/")
    }
  }
  if (keep_gene_pairs == "all"){
    warning('keep_gene_pairs = "all" can create very large output strings.')
    for (i in which(M > 0)){
      total_members[i] <- paste(.InterMetaPairMembers(
        pair_context,
        pathway.pair$pathway_genes[[grid$Source.ID[i]]],
        pathway.pair$pathway_genes[[grid$Target.ID[i]]]
      ), collapse = "/")
    }
  }

  res <- data.frame(
    Pathway.Interaction = paste(source_pathway, target_pathway, sep = "__"),
    Favorable.Pathway = source_pathway,
    Unfavorable.Pathway = target_pathway,
    No.Total.MGP = M,
    No.Shared.MGP = k,
    Enrichment.Ratio = round(ratio, 3),
    P.Value = p_value,
    Adjust.P.Value = stats::p.adjust(p_value, method = p.adj),
    Shared.MGP = shared_members,
    Total.MGP = total_members,
    stringsAsFactors = FALSE
  )
  res$Rank <- rank(res$P.Value, ties.method = "min")
  res <- res[order(res$P.Value, res$Adjust.P.Value), , drop = FALSE]
  res$Index <- seq_len(nrow(res))
  res[, c("Index", "Rank", "Pathway.Interaction",
          "Favorable.Pathway", "Unfavorable.Pathway",
          "No.Total.MGP", "No.Shared.MGP", "Enrichment.Ratio",
          "P.Value", "Adjust.P.Value", "Shared.MGP", "Total.MGP")]
}

.InterMetaEnrichExpanded <- function(gene,
                                     background,
                                     pathway.pair,
                                     p.adj,
                                     keep_gene_pairs){
  required <- c("Pathway.Pair", "Source.Pathway", "Target.Pathway",
                "Metabolic.Gene.Pair")
  if (!all(required %in% names(pathway.pair))){
    stop("Expanded pathway-pair data are missing required columns.")
  }
  gene <- unique(as.character(gene))
  background <- unique(as.character(background))
  annotated <- intersect(
    background,
    unique(as.character(pathway.pair$Metabolic.Gene.Pair))
  )
  N <- length(annotated)
  gene <- intersect(gene, annotated)
  n <- length(gene)
  split_pair <- split(pathway.pair, pathway.pair$Pathway.Pair)
  rows <- lapply(split_pair, function(subset){
    members <- intersect(unique(subset$Metabolic.Gene.Pair), annotated)
    shared <- intersect(gene, members)
    M <- length(members)
    k <- length(shared)
    data.frame(
      Pathway.Interaction = subset$Pathway.Pair[1],
      Favorable.Pathway = subset$Source.Pathway[1],
      Unfavorable.Pathway = subset$Target.Pathway[1],
      No.Total.MGP = M,
      No.Shared.MGP = k,
      Enrichment.Ratio = if (M == 0 || n == 0) NA_real_ else k * N / M / n,
      P.Value = if (M == 0 || n == 0) 1 else
        stats::phyper(k - 1, M, N - M, n, lower.tail = FALSE),
      Shared.MGP = if (keep_gene_pairs == "none") "" else
        paste(shared, collapse = "/"),
      Total.MGP = if (keep_gene_pairs == "all") paste(members, collapse = "/")
      else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  if (length(rows) == 0){
    return(.InterMetaEmptyEnrichment())
  }
  res <- do.call(rbind, rows)
  res$Enrichment.Ratio <- round(res$Enrichment.Ratio, 3)
  res$Adjust.P.Value <- stats::p.adjust(res$P.Value, method = p.adj)
  res$Rank <- rank(res$P.Value, ties.method = "min")
  res <- res[order(res$P.Value, res$Adjust.P.Value), , drop = FALSE]
  res$Index <- seq_len(nrow(res))
  res[, c("Index", "Rank", "Pathway.Interaction",
          "Favorable.Pathway", "Unfavorable.Pathway",
          "No.Total.MGP", "No.Shared.MGP", "Enrichment.Ratio",
          "P.Value", "Adjust.P.Value", "Shared.MGP", "Total.MGP")]
}

InterMetaEnrich <- function(gene = NULL,
                            background = NULL,
                            pathway.pair = NULL,
                            p.adj = "BH",
                            core = 1,
                            keep_gene_pairs = c("shared", "none", "all")){
  if (is.null(gene)){
    stop('"gene" is required. Use character(0) when no pairs were selected.')
  }
  if (is.null(background)){
    stop('"background" is required.')
  }
  if (is.null(pathway.pair)){
    stop('"pathway.pair" is required.')
  }
  p.adj <- .InterMetaValidatePAdjust(p.adj)
  keep_gene_pairs <- match.arg(keep_gene_pairs)
  if (inherits(pathway.pair, "InterMetaPathwayPair")){
    context <- .InterMetaPreparePairContext(gene, background)
    return(.InterMetaEnrichCompact(
      pathway.pair = pathway.pair,
      pair_context = context,
      p.adj = p.adj,
      keep_gene_pairs = keep_gene_pairs
    ))
  }
  .InterMetaEnrichExpanded(
    gene = gene,
    background = background,
    pathway.pair = pathway.pair,
    p.adj = p.adj,
    keep_gene_pairs = keep_gene_pairs
  )
}

InterMetaPermFDR <- function(resdata = NULL,
                             pathwayset = NULL,
                             gene = NULL,
                             background = NULL,
                             niter = 100,
                             core = 1,
                             seed = 1,
                             analysis_profile = "TimiGP-legacy"){
  analysis_profile <- .InterMetaAnalysisProfile(analysis_profile)
  if (is.null(resdata) || !is.data.frame(resdata)){
    stop('"resdata" is required.')
  }
  normalized_pathwayset <- .InterMetaNormalizeGeneset(pathwayset)
  permutation_pathwayset <- .InterMetaLegacyPermutationGeneset(
    pathwayset,
    normalized = normalized_pathwayset
  )
  if (is.null(gene)){
    stop('"gene" is required.')
  }
  if (is.null(background)){
    stop('"background" is required.')
  }
  niter <- as.integer(niter)
  core <- as.integer(core)
  seed <- as.integer(seed)
  if (is.na(niter) || niter < 1){
    stop('"niter" must be a positive integer.')
  }
  if (is.na(core) || core < 1){
    stop('"core" must be a positive integer.')
  }
  if (is.na(seed)){
    stop('"seed" must be an integer.')
  }
  if (nrow(resdata) == 0){
    resdata$Permutation.FDR <- numeric()
    return(resdata)
  }
  if (!"Rank" %in% names(resdata) ||
      length(resdata$Rank) != nrow(resdata)){
    stop('TimiGP-legacy permutation FDR requires a "Rank" column.')
  }

  context <- .InterMetaPreparePairContext(gene, background)
  iteration_seeds <- seq_len(niter)
  run_iteration <- function(iteration_seed){
    set.seed(iteration_seed)
    permuted <- permutation_pathwayset
    permuted$Gene <- sample(permuted$Gene, replace = FALSE)
    pair <- InterMetaPathwayPair(permuted)
    result <- .InterMetaEnrichCompact(
      pathway.pair = pair,
      pair_context = context,
      p.adj = "none",
      keep_gene_pairs = "none"
    )
    result$P.Value
  }

  if (core > 1 && niter > 1){
    worker_count <- min(core, niter)
    if (.Platform$OS.type != "windows"){
      permutation_p <- unlist(parallel::mclapply(
        iteration_seeds,
        run_iteration,
        mc.cores = worker_count,
        mc.preschedule = TRUE
      ), use.names = FALSE)
    } else {
      cluster <- parallel::makeCluster(worker_count)
      on.exit({
        parallel::stopCluster(cluster)
        foreach::registerDoSEQ()
      }, add = TRUE)
      worker_libraries <- .libPaths()
      parallel::clusterCall(cluster, function(paths){
        .libPaths(paths)
        NULL
      }, worker_libraries)
      doParallel::registerDoParallel(cluster)
      permutation_p <- foreach::foreach(
        iteration_seed = iteration_seeds,
        .combine = c,
        .inorder = TRUE,
        .packages = "InterMeta"
      ) %dopar% {
        run_iteration(iteration_seed)
      }
      foreach::registerDoSEQ()
    }
  } else {
    permutation_p <- unlist(lapply(iteration_seeds, run_iteration),
                            use.names = FALSE)
  }

  observed_p <- resdata$P.Value
  observed_n <- suppressWarnings(as.numeric(resdata$Rank))
  raw_fdr <- vapply(seq_along(observed_p), function(i){
    if (!is.finite(observed_p[i]) || !is.finite(observed_n[i]) ||
        observed_n[i] <= 0){
      return(NA_real_)
    }
    min(sum(permutation_p <= observed_p[i], na.rm = TRUE) /
          niter / observed_n[i], 1)
  }, numeric(1))
  resdata$Permutation.FDR <- raw_fdr
  columns <- c("Index", "Rank", "Pathway.Interaction",
               "Favorable.Pathway", "Unfavorable.Pathway",
               "No.Total.MGP", "No.Shared.MGP", "Enrichment.Ratio",
               "P.Value", "Adjust.P.Value", "Permutation.FDR",
               "Shared.MGP", "Total.MGP")
  resdata[, intersect(columns, names(resdata)), drop = FALSE]
}

InterMetaPathwayNetwork <- function(resdata = NULL,
                                    select = NULL,
                                    pathwayset = NULL,
                                    condition = "Adjust.P.Value",
                                    cutoff = 0.05,
                                    export = TRUE,
                                    path = NULL){
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  if (!condition %in% c("P.Value", "Adjust.P.Value", "Permutation.FDR")){
    stop('Please choose one of c("P.Value", "Adjust.P.Value", "Permutation.FDR").')
  }
  if (!condition %in% names(resdata)){
    stop("The enrichment result does not contain column: ", condition)
  }

  if (is.null(select)){
    resdata <- resdata[is.finite(resdata[[condition]]) &
                         resdata[[condition]] < cutoff, , drop = FALSE]
  } else {
    resdata <- resdata[resdata$Index %in% select, , drop = FALSE]
  }
  if (nrow(resdata) == 0){
    warning("No pathway interactions were selected.")
    net <- data.frame(Source = character(), Interaction = character(),
                      Target = character(), stringsAsFactors = FALSE)
    node <- data.frame(Pathway = character(), stringsAsFactors = FALSE)
    edge <- resdata
    if (export){
      if (is.null(path)){
        path <- "./"
      }
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      utils::write.table(node, file = file.path(path, "pathway_node.txt"),
                         quote = FALSE, sep = "\t", row.names = FALSE)
      utils::write.table(net, file = file.path(path, "pathway_network.sif"),
                         quote = FALSE, sep = "\t", row.names = FALSE,
                         col.names = FALSE)
      utils::write.table(edge, file = file.path(path, "pathway_edge.txt"),
                         quote = FALSE, sep = "\t", row.names = FALSE)
    }
    return(list(network = net, node = node, edge = edge))
  }

  net <- data.frame(Source = resdata$Favorable.Pathway,
                    Interaction = "InterMeta",
                    Target = resdata$Unfavorable.Pathway,
                    stringsAsFactors = FALSE)
  edge <- resdata
  edge$Key <- paste0(edge$Favorable.Pathway, " (InterMeta) ", edge$Unfavorable.Pathway)
  edge <- edge[, c("Key", setdiff(names(edge), "Key"))]

  selected <- unique(c(net$Source, net$Target))
  node <- data.frame(Pathway = selected, stringsAsFactors = FALSE)
  if (!is.null(pathwayset)){
    counts <- table(pathwayset$Pathway)
    node$No.Genes <- as.integer(counts[node$Pathway])
    node$No.Genes[is.na(node$No.Genes)] <- 0
    if ("Dataset" %in% names(pathwayset)){
      node$Dataset <- vapply(node$Pathway, function(x){
        paste(unique(pathwayset$Dataset[pathwayset$Pathway == x]), collapse = "/")
      }, character(1))
    }
  }

  if (export == TRUE){
    if (is.null(path)){
      warning("No output directory was provided. Export to working directory.")
      path <- "./"
    }
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    utils::write.table(node, file = file.path(path, "pathway_node.txt"),
                       quote = FALSE, sep = "\t", row.names = FALSE)
    utils::write.table(net, file = file.path(path, "pathway_network.sif"),
                       quote = FALSE, sep = "\t", row.names = FALSE,
                       col.names = FALSE)
    utils::write.table(edge, file = file.path(path, "pathway_edge.txt"),
                       quote = FALSE, sep = "\t", row.names = FALSE)
  }

  list(network = net, node = node, edge = edge)
}

InterMetaPathwayScore <- function(resdata = NULL,
                                  condition = "Adjust.P.Value",
                                  cutoff = 0.05,
                                  analysis_profile = "TimiGP-legacy"){
  analysis_profile <- .InterMetaAnalysisProfile(analysis_profile)
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  if (!condition %in% c("P.Value", "Adjust.P.Value", "Permutation.FDR")){
    stop('Please choose one of c("P.Value", "Adjust.P.Value", "Permutation.FDR").')
  }
  if (!condition %in% names(resdata)){
    stop("The enrichment result does not contain column: ", condition)
  }

  resdata <- resdata[is.finite(resdata[[condition]]) &
                       resdata[[condition]] < cutoff, , drop = FALSE]
  if (nrow(resdata) == 0){
    warning("No pathway interactions were selected; returning an empty score.")
    return(.InterMetaEmptyScore())
  }
  total_pathway <- unique(c(resdata$Favorable.Pathway, resdata$Unfavorable.Pathway))
  fa <- as.data.frame(table(resdata$Favorable.Pathway), stringsAsFactors = FALSE)
  unfa <- as.data.frame(table(resdata$Unfavorable.Pathway), stringsAsFactors = FALSE)
  names(fa) <- c("Pathway", "OutDegree")
  names(unfa) <- c("Pathway", "InDegree")
  fa <- merge(data.frame(Pathway = total_pathway), fa, all.x = TRUE)
  unfa <- merge(data.frame(Pathway = total_pathway), unfa, all.x = TRUE)
  fa$OutDegree[is.na(fa$OutDegree)] <- 0
  unfa$InDegree[is.na(unfa$InDegree)] <- 0
  score <- merge(fa, unfa, by = "Pathway")
  score$Favorable.Score <- score$OutDegree / sum(score$OutDegree) * 100
  score$Unfavorable.Score <- score$InDegree / sum(score$InDegree) * 100
  score$Difference <- score$Favorable.Score - score$Unfavorable.Score
  score <- score[score$Difference != 0, , drop = FALSE]
  if (nrow(score) == 0){
    warning("Only neutral pathway scores were found; returning an empty score.")
    return(.InterMetaEmptyScore())
  }
  score$Group <- ifelse(score$Difference > 0, "F", "U")
  score$Rank <- rank(-score$Difference, ties.method = "min")
  score$Signed.Rank <- score$Rank
  score <- score[order(score$Rank), ]
  score$Index <- seq_len(nrow(score))
  score[, c("Index", "Rank", "Signed.Rank", "Pathway", "OutDegree",
            "Favorable.Score", "InDegree", "Unfavorable.Score",
            "Difference", "Group")]
}

InterMetaDotplot <- function(resdata = NULL,
                             select = NULL,
                             condition = "Adjust.P.Value",
                             cutoff = 0.05,
                             top = 10L){
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  if (!condition %in% c("P.Value", "Adjust.P.Value", "Permutation.FDR")){
    stop('Please choose one of c("P.Value", "Adjust.P.Value", "Permutation.FDR").')
  }
  if (!condition %in% names(resdata)){
    stop("The enrichment result does not contain column: ", condition)
  }

  required <- c("Index", "Favorable.Pathway", "Unfavorable.Pathway",
                "Enrichment.Ratio", "No.Shared.MGP")
  missing <- setdiff(required, names(resdata))
  if (length(missing) > 0L){
    stop("The enrichment result is missing required column(s): ",
         paste(missing, collapse = ", "))
  }
  if (is.null(select)){
    eligible <- which(
      is.finite(resdata[[condition]]) & resdata[[condition]] < cutoff
    )
    if (length(eligible) > 0L){
      eligible <- eligible[order(resdata$Index[eligible])]
    }
    if (!is.null(top)){
      if (length(top) != 1L || !is.numeric(top) || is.na(top) ||
          top < 1 || top != floor(top)){
        stop('The parameter "top" must be NULL or a positive integer.')
      }
      eligible <- utils::head(eligible, as.integer(top))
    }
    select <- resdata$Index[eligible]
  }
  select <- unique(select[!is.na(select)])
  se <- match(as.character(select), as.character(resdata$Index), nomatch = 0L)
  se <- se[se > 0L]
  plotdata <- resdata[se, , drop = FALSE]
  if (nrow(plotdata) == 0){
    return(ggplot2::ggplot() +
             ggplot2::theme_void() +
             ggplot2::labs(title = "No pathway interactions selected"))
  }
  plotdata$Display.Interaction <- paste0(
    plotdata$Favorable.Pathway, " ->\n", plotdata$Unfavorable.Pathway
  )
  plotdata$Display.Interaction <- factor(
    plotdata$Display.Interaction,
    levels = unique(plotdata$Display.Interaction)
  )
  plotdata$Condition.Value <- plotdata[[condition]]
  ymax <- max(1, max(ceiling(plotdata$Enrichment.Ratio), na.rm = TRUE) + 1)
  condition_label <- switch(
    condition,
    P.Value = "P-value",
    Adjust.P.Value = "Adjusted P-value",
    Permutation.FDR = "FDR",
    condition
  )
  ggplot2::ggplot(plotdata) +
    ggplot2::geom_point(
      ggplot2::aes(x = Display.Interaction, y = Enrichment.Ratio,
                   color = Condition.Value, size = No.Shared.MGP)) +
    ggplot2::scale_color_continuous(
      low = "#9970ab", high = "#5aae61",
      limits = c(0, cutoff), name = condition_label,
      guide = ggplot2::guide_colorbar(reverse = TRUE)
    ) +
    ggplot2::scale_size(range = c(4, 9), name = "No. Shared\nGene Pairs") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::annotate(
      "text", x = 0.55, y = ymax * 0.55,
      label = "Favorable Outcome", angle = 90,
      color = "#9970ab", size = 3.6
    ) +
    ggplot2::labs(x = "Pathway Interaction", y = "Enrichment Ratio",
                  title = "Pathway Interaction Enrichment") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(
        angle = 40, vjust = 1, hjust = 1, size = 9
      ),
      axis.line = ggplot2::element_line(
        arrow = grid::arrow(length = grid::unit(0.11, "cm"), type = "closed")
      ),
      legend.position = "right",
      plot.margin = ggplot2::margin(10, 18, 12, 105)
    ) +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(add = c(0.6, 0.35))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, ymax),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    )
}

InterMetaPathwayScoreBar <- function(score = NULL,
                                     select = NULL,
                                     top = 5L){
  if (is.null(score)){
    stop('The parameter "score" is required.')
  }
  if (!is.data.frame(score)){
    stop('The parameter "score" must be a data frame.')
  }
  required <- c("Index", "Pathway", "Favorable.Score",
                "Unfavorable.Score", "Difference")
  missing <- setdiff(required, names(score))
  if (length(missing) > 0L){
    stop("The score is missing required column(s): ",
         paste(missing, collapse = ", "))
  }
  if (!is.null(select)){
    select <- unique(select[!is.na(select)])
    position <- match(as.character(select), as.character(score$Index),
                      nomatch = 0L)
    score <- score[position[position > 0L], , drop = FALSE]
  } else if (!is.null(top)){
    if (length(top) != 1L || !is.numeric(top) || is.na(top) ||
        top < 1 || top != floor(top)){
      stop('The parameter "top" must be NULL or a positive integer.')
    }
    favorable <- which(is.finite(score$Difference) & score$Difference > 0)
    unfavorable <- which(is.finite(score$Difference) & score$Difference < 0)
    favorable <- favorable[order(-score$Difference[favorable])]
    unfavorable <- unfavorable[order(score$Difference[unfavorable])]
    favorable <- utils::head(favorable, as.integer(top))
    unfavorable <- utils::head(unfavorable, as.integer(top))
    score <- score[c(favorable, unfavorable), , drop = FALSE]
  }
  if (nrow(score) == 0){
    return(ggplot2::ggplot() +
             ggplot2::theme_void() +
             ggplot2::labs(title = "No pathway scores available"))
  }
  score <- score[order(-score$Difference), , drop = FALSE]
  score$Display.Pathway <- vapply(as.character(score$Pathway), function(x){
    paste(strwrap(x, width = 22), collapse = "\n")
  }, character(1))
  score$Display.Pathway <- factor(
    score$Display.Pathway,
    levels = unique(score$Display.Pathway)
  )
  limit <- max(c(score$Favorable.Score, score$Unfavorable.Score), na.rm = TRUE)
  limit <- max(2, ceiling(limit / 2) * 2)

  ggplot2::ggplot(score, ggplot2::aes(x = Display.Pathway)) +
    ggplot2::geom_bar(data = score,
                      ggplot2::aes(y = Favorable.Score),
                      width = 0.76, stat = "identity", alpha = 0.5,
                      fill = "#fc4e2a") +
    ggplot2::geom_bar(data = score,
                      ggplot2::aes(y = -Unfavorable.Score),
                      width = 0.76, stat = "identity", alpha = 0.5,
                      fill = "#253494") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.6) +
    ggplot2::annotate(
      "text", x = 0.55, y = limit * 0.55,
      label = "Favorable Score", angle = 90,
      color = "#fc4e2a", size = 3.5
    ) +
    ggplot2::annotate(
      "text", x = 0.55, y = -limit * 0.55,
      label = "Unfavorable Score", angle = 90,
      color = "#253494", size = 3.5
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(x = "Pathway Type", y = "Favorability Score",
                  title = "Favorability") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(
        angle = 55, vjust = 1, hjust = 1, size = 9
      ),
      axis.line = ggplot2::element_line(
        arrow = grid::arrow(length = grid::unit(0.11, "cm"), type = "closed")
      ),
      plot.margin = ggplot2::margin(10, 14, 12, 12)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-limit, limit),
      breaks = pretty(c(-limit, limit), n = 5),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    )
}

InterMetaStabilityScoreBar <- function(
    stability = NULL,
    select = NULL,
    top = 5L,
    min_selection_frequency = 1,
    min_direction_consistency = 1,
    error = c("sd", "se", "none")){
  if (is.null(stability)){
    stop('The parameter "stability" is required.')
  }
  summary <- if (inherits(stability, "InterMetaStability")) {
    stability$summary
  } else if (is.list(stability) && is.data.frame(stability$summary)) {
    stability$summary
  } else if (is.data.frame(stability)) {
    stability
  } else {
    stop('"stability" must be an InterMetaStability object or summary data frame.')
  }
  required <- c(
    "Pathway", "Selection.Count", "Selection.Frequency",
    "Dominant.Direction", "Direction.Consistency", "Mean.Difference",
    "SD.Difference"
  )
  missing <- setdiff(required, names(summary))
  if (length(missing) > 0L){
    stop("The stability summary is missing required column(s): ",
         paste(missing, collapse = ", "))
  }
  validate_frequency <- function(x, name){
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
        x < 0 || x > 1){
      stop('"', name, '" must be a finite number from 0 to 1.')
    }
  }
  validate_frequency(min_selection_frequency, "min_selection_frequency")
  validate_frequency(min_direction_consistency,
                     "min_direction_consistency")
  error <- match.arg(error)

  plotdata <- summary[
    is.finite(summary$Mean.Difference) &
      summary$Selection.Frequency >= min_selection_frequency &
      is.finite(summary$Direction.Consistency) &
      summary$Direction.Consistency >= min_direction_consistency &
      summary$Dominant.Direction %in% c("F", "U"),
    , drop = FALSE
  ]
  if (!is.null(select)){
    select <- unique(as.character(select[!is.na(select)]))
    position <- match(select, as.character(plotdata$Pathway), nomatch = 0L)
    plotdata <- plotdata[position[position > 0L], , drop = FALSE]
  } else if (!is.null(top)){
    if (length(top) != 1L || !is.numeric(top) || is.na(top) ||
        top < 1 || top != floor(top)){
      stop('The parameter "top" must be NULL or a positive integer.')
    }
    favorable <- which(plotdata$Dominant.Direction == "F")
    unfavorable <- which(plotdata$Dominant.Direction == "U")
    favorable <- favorable[order(-plotdata$Mean.Difference[favorable])]
    unfavorable <- unfavorable[order(plotdata$Mean.Difference[unfavorable])]
    favorable <- utils::head(favorable, as.integer(top))
    unfavorable <- utils::head(unfavorable, as.integer(top))
    plotdata <- plotdata[c(favorable, unfavorable), , drop = FALSE]
  }
  if (nrow(plotdata) == 0L){
    empty_plot <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "No stable pathway scores available")
    attr(empty_plot, "plot_data") <- plotdata
    return(empty_plot)
  }
  if (error != "none" && any(plotdata$Selection.Count < 2L)){
    stop(
      "At least two successful resample scores per displayed pathway are ",
      "required for SD or SE error bars. Use error = \"none\" for a ",
      "single-run summary."
    )
  }

  plotdata <- plotdata[order(-plotdata$Mean.Difference), , drop = FALSE]
  plotdata$Direction <- ifelse(
    plotdata$Dominant.Direction == "F", "Favorable", "Unfavorable"
  )
  plotdata$Error <- switch(
    error,
    sd = plotdata$SD.Difference,
    se = plotdata$SD.Difference / sqrt(plotdata$Selection.Count),
    none = rep(0, nrow(plotdata))
  )
  plotdata$Error[!is.finite(plotdata$Error)] <- 0
  plotdata$Lower <- plotdata$Mean.Difference - plotdata$Error
  plotdata$Upper <- plotdata$Mean.Difference + plotdata$Error
  plotdata$Display.Pathway <- vapply(
    as.character(plotdata$Pathway),
    function(x) paste(strwrap(x, width = 22), collapse = "\n"),
    character(1)
  )
  plotdata$Display.Pathway <- factor(
    plotdata$Display.Pathway,
    levels = unique(plotdata$Display.Pathway)
  )
  limit <- max(abs(c(plotdata$Lower, plotdata$Upper)), na.rm = TRUE)
  limit <- max(2, ceiling(limit / 2) * 2)

  plot <- ggplot2::ggplot(
    plotdata,
    ggplot2::aes(x = Display.Pathway, y = Mean.Difference, fill = Direction)
  ) +
    ggplot2::geom_col(
      width = 0.76, color = "black", linewidth = 0.25
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Lower, ymax = Upper),
      width = 0.24,
      linewidth = 0.55
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.6) +
    ggplot2::scale_fill_manual(
      values = c(Favorable = "#FDA694", Unfavorable = "#9299C9")
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-limit, limit),
      breaks = pretty(c(-limit, limit), n = 5),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      x = "Pathway Type",
      y = "Difference",
      fill = "Direction"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, vjust = 1, hjust = 1, size = 9
      ),
      axis.line = ggplot2::element_line(
        arrow = grid::arrow(length = grid::unit(0.11, "cm"), type = "closed")
      ),
      legend.position = "top",
      legend.direction = "horizontal",
      plot.margin = ggplot2::margin(8, 14, 12, 12)
    )
  attr(plot, "plot_data") <- plotdata
  plot
}

.InterMetaNormalizeScore <- function(score,
                                     run = NA_integer_){
  if (is.null(score)){
    score <- .InterMetaEmptyScore()
  }
  if ("Cell.Type" %in% names(score) && !"Pathway" %in% names(score)){
    names(score)[names(score) == "Cell.Type"] <- "Pathway"
  }
  required <- c("Pathway", "Rank", "Difference", "Group")
  missing <- setdiff(required, names(score))
  if (length(missing) > 0){
    stop("Score is missing required column(s): ",
         paste(missing, collapse = ", "))
  }
  score$Pathway <- as.character(score$Pathway)
  score$Rank <- as.numeric(score$Rank)
  score$Difference <- as.numeric(score$Difference)
  score$Group <- as.character(score$Group)
  if (anyNA(score$Pathway) || any(!nzchar(trimws(score$Pathway)))){
    stop("Score Pathway values must be non-missing and non-empty.")
  }
  score$Group[score$Group %in% c("Favorable", "favorable")] <- "F"
  score$Group[score$Group %in% c("Unfavorable", "unfavorable")] <- "U"
  invalid_group <- is.na(score$Group) | !score$Group %in% c("F", "U", "N")
  if (any(invalid_group)){
    stop(
      "Score Group values must be F, U, N, Favorable, or Unfavorable."
    )
  }
  if (anyDuplicated(score$Pathway)){
    stop("Each score table must contain at most one row per pathway.")
  }
  score$Run <- rep(run, nrow(score))
  score
}

.InterMetaJaccard <- function(x,
                              y){
  union_size <- length(union(x, y))
  if (union_size == 0){
    return(NA_real_)
  }
  length(intersect(x, y)) / union_size
}

.InterMetaMagnitudeOrder <- function(score,
                                     group = NULL){
  keep <- is.finite(score$Difference)
  if (!is.null(group)){
    keep <- keep & !is.na(score$Group) & score$Group == group
  }
  eligible <- which(keep)
  if (length(eligible) == 0L){
    return(integer())
  }
  eligible[order(
    -abs(score$Difference[eligible]),
    score$Pathway[eligible],
    method = "radix",
    na.last = TRUE
  )]
}

.InterMetaMagnitudeTopSet <- function(score,
                                      k,
                                      group = NULL){
  ordered <- .InterMetaMagnitudeOrder(score, group = group)
  if (length(ordered) == 0L){
    return(character())
  }
  as.character(score$Pathway[utils::head(ordered, min(k, length(ordered)))])
}

.InterMetaMagnitudeRank <- function(score){
  magnitude_rank <- rep(NA_real_, nrow(score))
  ordered <- .InterMetaMagnitudeOrder(score)
  if (length(ordered) > 0L){
    magnitude_rank[ordered] <- as.numeric(seq_along(ordered))
  }
  magnitude_rank
}

.InterMetaTopSetDefinition <- function(top){
  list(
    K = as.integer(top),
    eligibility = "finite Difference",
    ranking = "descending abs(Difference)",
    tie_break = "ascending Pathway using radix ordering",
    size = "exactly min(K, eligible pathway count)",
    scopes = c(
      generic = "all eligible pathways",
      favorable = "eligible pathways with Group == 'F'",
      unfavorable = "eligible pathways with Group == 'U'"
    ),
    empty_set_jaccard = "NA (undefined when neither run has an eligible pathway)",
    classification_top = as.integer(max(top)),
    classification_rule = paste0(
      "Stability.Class uses Top", max(top),
      ".Frequency plus Selection.Frequency and Direction.Consistency; ",
      "fewer than two successful runs are Not assessed"
    )
  )
}

.InterMetaEmptyPairwiseStability <- function(top){
  out <- data.frame(
    Run.A = character(),
    Run.B = character(),
    Common.Pathways = integer(),
    Difference.Pearson.Selected = numeric(),
    Difference.Spearman.Selected = numeric(),
    Difference.Pearson.All = numeric(),
    Difference.Spearman.All = numeric(),
    Direction.Agreement = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (k in top){
    out[[paste0("Top", k, ".Jaccard")]] <- numeric()
    out[[paste0("Top", k, ".Favorable.Jaccard")]] <- numeric()
    out[[paste0("Top", k, ".Unfavorable.Jaccard")]] <- numeric()
  }
  out
}

.InterMetaStabilitySummaryTemplate <- function(pathways,
                                                top){
  n_pathway <- length(pathways)
  out <- data.frame(
    Pathway = as.character(pathways),
    Selection.Count = rep(0L, n_pathway),
    Selection.Frequency = rep(0, n_pathway),
    Favorable.Frequency = rep(0, n_pathway),
    Unfavorable.Frequency = rep(0, n_pathway),
    Neutral.Frequency = rep(0, n_pathway),
    Dominant.Direction = rep("Not selected", n_pathway),
    Direction.Consistency = rep(NA_real_, n_pathway),
    Mean.Difference = rep(NA_real_, n_pathway),
    SD.Difference = rep(NA_real_, n_pathway),
    CV.Difference = rep(NA_real_, n_pathway),
    Mean.Rank = rep(NA_real_, n_pathway),
    SD.Rank = rep(NA_real_, n_pathway),
    Mean.Magnitude.Rank = rep(NA_real_, n_pathway),
    SD.Magnitude.Rank = rep(NA_real_, n_pathway),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (k in top){
    out[[paste0("Top", k, ".Frequency")]] <- rep(0, n_pathway)
    out[[paste0("Top", k, ".Favorable.Frequency")]] <- rep(0, n_pathway)
    out[[paste0("Top", k, ".Unfavorable.Frequency")]] <- rep(0, n_pathway)
  }
  out$Stability.Class <- rep("Not assessed", n_pathway)
  out
}

.InterMetaSafeCor <- function(x,
                              y,
                              method){
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2 || stats::sd(x[keep]) == 0 || stats::sd(y[keep]) == 0){
    return(NA_real_)
  }
  stats::cor(x[keep], y[keep], method = method)
}

.InterMetaPairwiseStability <- function(scores,
                                        top = c(5, 10, 20),
                                        pathways = NULL,
                                        run_ids = seq_along(scores)){
  if (length(scores) < 2){
    return(.InterMetaEmptyPairwiseStability(top))
  }
  pathways <- sort(unique(c(
    as.character(pathways),
    unlist(lapply(scores, function(x) x$Pathway), use.names = FALSE)
  )))
  pairs <- utils::combn(seq_along(scores), 2)
  pairwise <- lapply(seq_len(ncol(pairs)), function(i){
    a <- pairs[1, i]
    b <- pairs[2, i]
    score_a <- scores[[a]]
    score_b <- scores[[b]]
    common <- intersect(score_a$Pathway, score_b$Pathway)
    sub_a <- score_a[match(common, score_a$Pathway), ]
    sub_b <- score_b[match(common, score_b$Pathway), ]
    difference_a <- rep(0, length(pathways))
    difference_b <- rep(0, length(pathways))
    names(difference_a) <- pathways
    names(difference_b) <- pathways
    difference_a[score_a$Pathway] <- score_a$Difference
    difference_b[score_b$Pathway] <- score_b$Difference
    out <- data.frame(
      Run.A = as.character(run_ids[a]),
      Run.B = as.character(run_ids[b]),
      Common.Pathways = length(common),
      Difference.Pearson.Selected = .InterMetaSafeCor(
        sub_a$Difference, sub_b$Difference, "pearson"
      ),
      Difference.Spearman.Selected = .InterMetaSafeCor(
        sub_a$Difference, sub_b$Difference, "spearman"
      ),
      Difference.Pearson.All = .InterMetaSafeCor(
        difference_a, difference_b, "pearson"
      ),
      Difference.Spearman.All = .InterMetaSafeCor(
        difference_a, difference_b, "spearman"
      ),
      Direction.Agreement = if (length(common) > 0) {
        mean(sub_a$Group == sub_b$Group, na.rm = TRUE)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
    for (k in top){
      top_a <- .InterMetaMagnitudeTopSet(score_a, k)
      top_b <- .InterMetaMagnitudeTopSet(score_b, k)
      out[[paste0("Top", k, ".Jaccard")]] <- .InterMetaJaccard(top_a, top_b)
      favorable_a <- .InterMetaMagnitudeTopSet(score_a, k, group = "F")
      favorable_b <- .InterMetaMagnitudeTopSet(score_b, k, group = "F")
      out[[paste0("Top", k, ".Favorable.Jaccard")]] <-
        .InterMetaJaccard(favorable_a, favorable_b)
      unfavorable_a <- .InterMetaMagnitudeTopSet(score_a, k, group = "U")
      unfavorable_b <- .InterMetaMagnitudeTopSet(score_b, k, group = "U")
      out[[paste0("Top", k, ".Unfavorable.Jaccard")]] <-
        .InterMetaJaccard(unfavorable_a, unfavorable_b)
    }
    out
  })
  do.call(rbind, pairwise)
}

.InterMetaSafeMean <- function(x){
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

.InterMetaSafeSD <- function(x){
  x <- x[is.finite(x)]
  if (length(x) < 2) NA_real_ else stats::sd(x)
}

InterMetaSummarizeStability <- function(stability = NULL,
                                        scores = NULL,
                                        pathways = NULL,
                                        top = c(5, 10, 20),
                                        high.frequency = 0.8,
                                        moderate.frequency = 0.5,
                                        high.direction = 0.8,
                                        moderate.direction = 0.6){
  if (is.null(scores)){
    if (is.null(stability) || is.null(stability$runs)){
      stop('Either "stability" or "scores" is required.')
    }
    successful <- vapply(stability$runs, function(x){
      is.null(x$status) || identical(x$status, "success")
    }, logical(1))
    selected_runs <- stability$runs[successful]
    scores <- lapply(selected_runs, function(x) x$score)
    run_ids <- vapply(selected_runs, function(x){
      if (is.null(x$run_id)) as.character(x$run) else x$run_id
    }, character(1))
    if (is.null(pathways)){
      pathways <- stability$candidate_pathways
    }
  } else {
    run_ids <- if (is.null(names(scores))) {
      as.character(seq_along(scores))
    } else {
      names(scores)
    }
  }
  if (!is.numeric(top) || length(top) == 0L || anyNA(top) ||
      any(!is.finite(top)) || any(top < 1) || any(top != floor(top)) ||
      any(top > .Machine$integer.max)){
    stop('"top" must contain positive integers.')
  }
  top <- sort(unique(as.integer(top)))
  top_set_definition <- .InterMetaTopSetDefinition(top)
  scores <- lapply(seq_along(scores), function(i){
    .InterMetaNormalizeScore(scores[[i]], run = run_ids[i])
  })
  n_run <- length(scores)

  all_pathway <- sort(unique(c(
    as.character(pathways),
    unlist(lapply(scores, function(x) x$Pathway), use.names = FALSE)
  )), method = "radix")
  summary <- .InterMetaStabilitySummaryTemplate(all_pathway, top)
  if (n_run == 0){
    return(list(
      summary = summary,
      pairwise = .InterMetaEmptyPairwiseStability(top),
      scores = scores,
      successful_runs = 0L,
      top_set_definition = top_set_definition
    ))
  }
  if (length(all_pathway) == 0L){
    pairwise <- .InterMetaPairwiseStability(
      scores = scores,
      top = top,
      pathways = all_pathway,
      run_ids = run_ids
    )
    return(list(
      summary = summary,
      pairwise = pairwise,
      scores = scores,
      successful_runs = n_run,
      top_set_definition = top_set_definition
    ))
  }

  selection_count <- vapply(all_pathway, function(pathway){
    sum(vapply(scores, function(score) pathway %in% score$Pathway, logical(1)))
  }, integer(1))
  group_count <- function(pathway, group){
    sum(vapply(scores, function(score){
      idx <- match(pathway, score$Pathway)
      !is.na(idx) && identical(score$Group[idx], group)
    }, logical(1)))
  }
  favorable_count <- vapply(all_pathway, group_count, integer(1), group = "F")
  unfavorable_count <- vapply(all_pathway, group_count, integer(1), group = "U")
  neutral_count <- vapply(all_pathway, group_count, integer(1), group = "N")
  summary$Selection.Count <- selection_count
  summary$Selection.Frequency <- selection_count / n_run
  summary$Favorable.Frequency <- favorable_count / n_run
  summary$Unfavorable.Frequency <- unfavorable_count / n_run
  summary$Neutral.Frequency <- neutral_count / n_run
  direction_counts <- cbind(F = favorable_count,
                            U = unfavorable_count,
                            N = neutral_count)
  summary$Dominant.Direction <- vapply(seq_len(nrow(direction_counts)), function(i){
    if (selection_count[i] == 0){
      return("Not selected")
    }
    maximum <- max(direction_counts[i, ])
    winner <- names(which(direction_counts[i, ] == maximum))
    if (length(winner) == 1) winner else "Mixed"
  }, character(1))
  summary$Direction.Consistency <- ifelse(
    selection_count > 0,
    apply(direction_counts, 1, max) / selection_count,
    NA_real_
  )

  diff_list <- lapply(all_pathway, function(pathway){
    values <- vapply(scores, function(score){
      idx <- match(pathway, score$Pathway)
      if (is.na(idx)) NA_real_ else score$Difference[idx]
    }, numeric(1))
    values
  })
  rank_list <- lapply(all_pathway, function(pathway){
    values <- vapply(scores, function(score){
      idx <- match(pathway, score$Pathway)
      if (is.na(idx)) NA_real_ else score$Rank[idx]
    }, numeric(1))
    values
  })
  magnitude_rank_by_run <- lapply(scores, .InterMetaMagnitudeRank)
  magnitude_rank_list <- lapply(all_pathway, function(pathway){
    vapply(seq_along(scores), function(i){
      idx <- match(pathway, scores[[i]]$Pathway)
      if (is.na(idx)) NA_real_ else magnitude_rank_by_run[[i]][idx]
    }, numeric(1))
  })

  summary$Mean.Difference <- vapply(diff_list, .InterMetaSafeMean, numeric(1))
  summary$SD.Difference <- vapply(diff_list, .InterMetaSafeSD, numeric(1))
  summary$CV.Difference <- ifelse(abs(summary$Mean.Difference) > 0,
                                  summary$SD.Difference / abs(summary$Mean.Difference),
                                  NA_real_)
  summary$Mean.Rank <- vapply(rank_list, .InterMetaSafeMean, numeric(1))
  summary$SD.Rank <- vapply(rank_list, .InterMetaSafeSD, numeric(1))
  summary$Mean.Magnitude.Rank <- vapply(
    magnitude_rank_list, .InterMetaSafeMean, numeric(1)
  )
  summary$SD.Magnitude.Rank <- vapply(
    magnitude_rank_list, .InterMetaSafeSD, numeric(1)
  )

  for (k in top){
    generic_sets <- lapply(scores, .InterMetaMagnitudeTopSet, k = k)
    favorable_sets <- lapply(
      scores, .InterMetaMagnitudeTopSet, k = k, group = "F"
    )
    unfavorable_sets <- lapply(
      scores, .InterMetaMagnitudeTopSet, k = k, group = "U"
    )
    summary[[paste0("Top", k, ".Frequency")]] <- vapply(all_pathway, function(pathway){
      mean(vapply(generic_sets, function(top_set){
        pathway %in% top_set
      }, logical(1)))
    }, numeric(1))
    summary[[paste0("Top", k, ".Favorable.Frequency")]] <- vapply(
      all_pathway,
      function(pathway){
        mean(vapply(favorable_sets, function(top_set){
          pathway %in% top_set
        }, logical(1)))
      },
      numeric(1)
    )
    summary[[paste0("Top", k, ".Unfavorable.Frequency")]] <- vapply(
      all_pathway,
      function(pathway){
        mean(vapply(unfavorable_sets, function(top_set){
          pathway %in% top_set
        }, logical(1)))
      },
      numeric(1)
    )
  }

  if (n_run < 2L){
    summary$Stability.Class <- rep("Not assessed", nrow(summary))
  } else {
    top_col <- paste0("Top", max(top), ".Frequency")
    high <-
      summary$Selection.Frequency >= high.frequency &
        summary[[top_col]] >= high.frequency &
        summary$Direction.Consistency >= high.direction
    moderate <-
      summary$Selection.Frequency >= moderate.frequency &
        summary[[top_col]] >= moderate.frequency &
        summary$Direction.Consistency >= moderate.direction
    high[is.na(high)] <- FALSE
    moderate[is.na(moderate)] <- FALSE
    summary$Stability.Class <- rep("Low", nrow(summary))
    summary$Stability.Class[moderate] <- "Moderate"
    summary$Stability.Class[high] <- "High"
    summary$Stability.Class[
      summary$Dominant.Direction == "N"
    ] <- "Neutral"
    summary$Stability.Class[
      summary$Dominant.Direction == "Not selected"
    ] <- "Not selected"
    summary$Stability.Class[
      summary$Dominant.Direction == "Mixed"
    ] <- "Low"
  }
  summary <- summary[order(-summary$Selection.Frequency,
                           -summary$Direction.Consistency,
                           summary$Mean.Magnitude.Rank,
                           summary$Pathway,
                           method = "radix",
                           na.last = TRUE), , drop = FALSE]
  row.names(summary) <- NULL

  pairwise <- .InterMetaPairwiseStability(
    scores = scores,
    top = top,
    pathways = all_pathway,
    run_ids = run_ids
  )
  list(
    summary = summary,
    pairwise = pairwise,
    scores = scores,
    successful_runs = n_run,
    top_set_definition = top_set_definition
  )
}

.InterMetaSummarizeEdges <- function(runs,
                                     condition,
                                     cutoff){
  successful <- runs[vapply(runs, function(x){
    identical(x$status, "success")
  }, logical(1))]
  n_success <- length(successful)
  empty <- data.frame(
    Pathway.Interaction = character(),
    Favorable.Pathway = character(),
    Unfavorable.Pathway = character(),
    Selection.Count = integer(),
    Selection.Frequency = numeric(),
    Mean.Enrichment.Ratio = numeric(),
    Mean.P.Value = numeric(),
    Mean.Adjust.P.Value = numeric(),
    stringsAsFactors = FALSE
  )
  if (n_success == 0){
    return(empty)
  }
  selected <- lapply(successful, function(run){
    enrich <- run$enrich
    if (is.null(enrich) || !condition %in% names(enrich)){
      return(NULL)
    }
    keep <- is.finite(enrich[[condition]]) & enrich[[condition]] < cutoff
    enrich <- enrich[keep, , drop = FALSE]
    if (nrow(enrich) == 0){
      return(NULL)
    }
    enrich$Run.ID <- run$run_id
    enrich
  })
  selected <- selected[!vapply(selected, is.null, logical(1))]
  if (length(selected) == 0){
    return(empty)
  }
  selected <- do.call(rbind, selected)
  interactions <- unique(selected$Pathway.Interaction)
  out <- lapply(interactions, function(interaction){
    subset <- selected[selected$Pathway.Interaction == interaction, ,
                       drop = FALSE]
    data.frame(
      Pathway.Interaction = interaction,
      Favorable.Pathway = subset$Favorable.Pathway[1],
      Unfavorable.Pathway = subset$Unfavorable.Pathway[1],
      Selection.Count = length(unique(subset$Run.ID)),
      Selection.Frequency = length(unique(subset$Run.ID)) / n_success,
      Mean.Enrichment.Ratio = .InterMetaSafeMean(subset$Enrichment.Ratio),
      Mean.P.Value = .InterMetaSafeMean(subset$P.Value),
      Mean.Adjust.P.Value = .InterMetaSafeMean(subset$Adjust.P.Value),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- out[order(-out$Selection.Frequency, out$Mean.P.Value), ,
             drop = FALSE]
  row.names(out) <- NULL
  out
}

.InterMetaStratifiedSample <- function(samples,
                                       event,
                                       size,
                                       seed,
                                       stratify,
                                       min_events,
                                       min_censored){
  set.seed(seed)
  if (!stratify){
    return(sample(samples, size = size, replace = FALSE))
  }
  event_samples <- samples[event == 1]
  censored_samples <- samples[event == 0]
  lower <- min_events
  upper <- min(length(event_samples), size - min_censored)
  if (lower > upper || length(censored_samples) < min_censored){
    stop("The requested sample size cannot satisfy the minimum event and ",
         "censoring counts.")
  }
  target_event <- round(size * length(event_samples) / length(samples))
  target_event <- min(max(target_event, lower), upper)
  target_censored <- size - target_event
  if (target_censored > length(censored_samples)){
    target_censored <- length(censored_samples)
    target_event <- size - target_censored
  }
  selected <- c(
    sample(event_samples, target_event, replace = FALSE),
    sample(censored_samples, target_censored, replace = FALSE)
  )
  sample(selected, length(selected), replace = FALSE)
}

.InterMetaCacheKey <- function(rna,
                               info,
                               pathwayset,
                               pathway_pair,
                               parameters){
  version <- tryCatch(
    as.character(utils::packageVersion("InterMeta")),
    error = function(e) "development"
  )
  digest::digest(
    list(
      rna = rna,
      info = info,
      pathwayset = pathwayset,
      pathway_pair = pathway_pair,
      parameters = parameters,
      package_version = version
    ),
    algo = "xxhash64"
  )
}

.InterMetaPermutationSeed <- function(seed,
                                      run){
  hash <- digest::digest(
    list(seed = as.integer(seed), run = as.integer(run)),
    algo = "xxhash32"
  )
  strtoi(substr(hash, 1L, 7L), base = 16L)
}

.InterMetaWriteTableAtomic <- function(x,
                                       file){
  temporary_file <- tempfile(
    pattern = paste0(".", basename(file), "-"),
    tmpdir = dirname(file),
    fileext = ".tmp"
  )
  previous_file <- paste0(file, ".previous")
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)

  if (!file.exists(file) && file.exists(previous_file)){
    if (!file.rename(previous_file, file)){
      stop("Could not restore the previous manifest before writing: ", file)
    }
  }
  unlink(previous_file, force = TRUE)
  utils::write.table(
    x,
    file = temporary_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  temporary_info <- file.info(temporary_file)
  if (!file.exists(temporary_file) ||
      is.na(temporary_info$size) || temporary_info$size <= 0){
    stop("Could not create a complete temporary manifest: ", temporary_file)
  }

  if (.Platform$OS.type != "windows"){
    if (!file.rename(temporary_file, file)){
      stop("Could not replace the manifest atomically: ", file)
    }
    return(invisible(file))
  }

  had_previous <- file.exists(file)
  if (had_previous && !file.rename(file, previous_file)){
    stop("Could not preserve the previous manifest before replacement: ", file)
  }
  if (!file.rename(temporary_file, file)){
    if (had_previous && file.exists(previous_file) && !file.exists(file)){
      file.rename(previous_file, file)
    }
    stop("Could not replace the manifest atomically: ", file)
  }
  unlink(previous_file, force = TRUE)
  invisible(file)
}

.InterMetaRunStatus <- function(runs){
  rows <- lapply(runs, function(run){
    data.frame(
      Run = run$run,
      Run.ID = run$run_id,
      Seed = run$seed,
      Permutation.Seed = if (is.null(run$permutation_seed)) {
        NA_integer_
      } else {
        run$permutation_seed
      },
      Status = run$status,
      Samples = length(run$selected_samples),
      Events = run$event_count,
      Censored = run$censored_count,
      Query.Pairs = run$query_pair_count,
      Score.Pathways = nrow(run$score),
      Elapsed.Seconds = run$elapsed_seconds,
      Error = if (is.null(run$error_message)) {
        ""
      } else {
        gsub("[\t\r\n]+", " ", as.character(run$error_message))
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.InterMetaLightRun <- function(run){
  run$cox_res <- NULL
  run$enrich <- NULL
  run$query_pairs <- NULL
  run$network <- NULL
  run
}

InterMetaStability <- function(rna = NULL,
                               info = NULL,
                               pathwayset = NULL,
                               pathway.pair = NULL,
                               n_resample = 3,
                               sample_fraction = 0.8,
                               seeds = NULL,
                               p.adj = "BH",
                               pair_condition = "PV",
                               pair_cutoff = 0.05,
                               enrich_p.adj = "BH",
                               score_condition = "Adjust.P.Value",
                               score_cutoff = 0.05,
                               run_permutation = FALSE,
                               permutation_niter = 100,
                               run_network = FALSE,
                               log = TRUE,
                               GMNorm = TRUE,
                               parallel = FALSE,
                               core = 1,
                               outdir = NULL,
                               reuse = TRUE,
                               stratify_event = FALSE,
                               min_events = 5,
                               min_censored = 5,
                               continue_on_error = TRUE,
                               keep_details = NULL,
                               analysis_profile = "TimiGP-legacy"){
  analysis_profile <- .InterMetaAnalysisProfile(analysis_profile)
  if (is.null(rna)){
    stop('The parameter "rna" is required.')
  }
  if (is.null(info)){
    stop('The parameter "info" is required.')
  }
  if (is.null(pathwayset)){
    stop('The parameter "pathwayset" is required.')
  }
  normalized_pathwayset <- .InterMetaNormalizeGeneset(pathwayset)
  permutation_pathwayset <- .InterMetaLegacyPermutationGeneset(
    pathwayset,
    normalized = normalized_pathwayset
  )
  pathwayset <- normalized_pathwayset
  p.adj <- .InterMetaValidatePAdjust(p.adj)
  enrich_p.adj <- .InterMetaValidatePAdjust(enrich_p.adj)
  if (!pair_condition %in% c("PV", "QV")){
    stop('Please choose one of c("PV", "QV") for pair_condition.')
  }
  if (pair_condition == "PV"){
    warning('pair_condition = "PV" uses unadjusted gene-pair P-values. ',
            'Use "QV" for confirmatory analysis or report a sensitivity analysis.')
  }
  if (!score_condition %in% c("P.Value", "Adjust.P.Value", "Permutation.FDR")){
    stop('Please choose one of c("P.Value", "Adjust.P.Value", "Permutation.FDR") for score_condition.')
  }
  if (score_condition == "Permutation.FDR" && run_permutation == FALSE){
    stop('score_condition = "Permutation.FDR" requires run_permutation = TRUE.')
  }
  if (is.null(rownames(info))){
    stop("The clinical info must have sample IDs as row names.")
  }
  if (is.null(colnames(rna))){
    stop("The RNA matrix must have sample IDs as column names.")
  }
  if (length(sample_fraction) != 1L || !is.numeric(sample_fraction) ||
      !is.finite(sample_fraction) || sample_fraction <= 0 ||
      sample_fraction > 1){
    stop("sample_fraction must be one finite number > 0 and <= 1.")
  }
  if (length(n_resample) != 1L || !is.numeric(n_resample) ||
      !is.finite(n_resample) || n_resample < 1 ||
      n_resample != floor(n_resample) ||
      n_resample > .Machine$integer.max){
    stop("n_resample must be one positive integer.")
  }
  n_resample <- as.integer(n_resample)
  if (!is.logical(stratify_event) || length(stratify_event) != 1 ||
      is.na(stratify_event)){
    stop('"stratify_event" must be TRUE or FALSE.')
  }
  if (!is.logical(continue_on_error) || length(continue_on_error) != 1 ||
      is.na(continue_on_error)){
    stop('"continue_on_error" must be TRUE or FALSE.')
  }
  positive_integer <- function(x){
    length(x) == 1L && is.numeric(x) && is.finite(x) &&
      x >= 1 && x == floor(x) && x <= .Machine$integer.max
  }
  if (!positive_integer(min_events) ||
      !positive_integer(min_censored) || !positive_integer(core)){
    stop("min_events, min_censored, and core must be positive integers.")
  }
  if (!positive_integer(permutation_niter)){
    stop("permutation_niter must be a positive integer.")
  }
  min_events <- as.integer(min_events)
  min_censored <- as.integer(min_censored)
  core <- as.integer(core)
  permutation_niter <- as.integer(permutation_niter)
  if (is.null(keep_details)){
    keep_details <- is.null(outdir)
  }
  if (!is.logical(keep_details) || length(keep_details) != 1 ||
      is.na(keep_details)){
    stop('"keep_details" must be TRUE or FALSE.')
  }
  if (is.null(outdir) && !keep_details){
    stop('keep_details = FALSE requires a non-NULL "outdir".')
  }

  common_samples <- intersect(colnames(rna), rownames(info))
  if (length(common_samples) < min_events + min_censored){
    stop("Too few common samples to satisfy minimum event/censoring counts.")
  }
  info <- InterMetaCheckEvent(
    info[common_samples, , drop = FALSE],
    min_events = min_events,
    min_censored = min_censored
  )
  common_samples <- intersect(common_samples, row.names(info))
  sample_size <- max(5, round(length(common_samples) * sample_fraction))
  sample_size <- min(sample_size, length(common_samples))
  if (sample_size < min_events + min_censored){
    stop("The resample size is smaller than min_events + min_censored.")
  }
  if (is.null(seeds)){
    seeds <- 20260701L + seq_len(n_resample) - 1L
  }
  if (length(seeds) < n_resample){
    stop("The length of seeds must be at least n_resample.")
  }
  seeds <- seeds[seq_len(n_resample)]
  if (!is.numeric(seeds) || any(!is.finite(seeds)) ||
      any(seeds != floor(seeds)) ||
      any(abs(seeds) > .Machine$integer.max)){
    stop("All used seeds must be finite integers in the R integer range.")
  }
  seeds <- as.integer(seeds)
  if (anyDuplicated(seeds)){
    warning("Duplicate seeds will generate duplicate subsamples.")
  }

  if (is.null(pathway.pair)){
    message("Preparing compact pathway-pair annotation once for all resamples.")
    pathway.pair <- InterMetaPathwayPair(pathwayset = pathwayset, core = core)
  } else if (inherits(pathway.pair, "InterMetaPathwayPair")){
    expected_hash <- digest::digest(
      unique(pathwayset[, c("Pathway", "Gene"), drop = FALSE]),
      algo = "xxhash64"
    )
    if (!identical(pathway.pair$geneset_hash, expected_hash)){
      stop("The supplied pathway.pair was created from a different geneset.")
    }
  }
  if (!is.null(outdir)){
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  parameters <- list(
    analysis_profile = analysis_profile,
    n_resample = n_resample,
    sample_fraction = sample_fraction,
    seeds = seeds,
    p.adj = p.adj,
    pair_condition = pair_condition,
    pair_cutoff = pair_cutoff,
    enrich_p.adj = enrich_p.adj,
    score_condition = score_condition,
    score_cutoff = score_cutoff,
    run_permutation = run_permutation,
    permutation_niter = permutation_niter,
    permutation_seed_scheme = "TimiGP-legacy-1:niter",
    run_network = run_network,
    log = log,
    GMNorm = GMNorm,
    stratify_event = stratify_event,
    min_events = min_events,
    min_censored = min_censored
  )
  cache_parameters <- parameters
  cache_parameters$permutation_geneset_rows <- nrow(permutation_pathwayset)
  cache_parameters$permutation_geneset_hash <- digest::digest(
    permutation_pathwayset[, c("Pathway", "Gene"), drop = FALSE],
    algo = "xxhash64"
  )
  cache_parameters$n_resample <- NULL
  cache_parameters$seeds <- NULL
  message("Computing the cache/provenance fingerprint.")
  analysis_hash <- .InterMetaCacheKey(
    rna = rna[, common_samples, drop = FALSE],
    info = info[common_samples, , drop = FALSE],
    pathwayset = pathwayset,
    pathway_pair = pathway.pair,
    parameters = cache_parameters
  )

  runs <- vector("list", n_resample)
  active_run_manifest <- data.frame(
    Analysis.Hash = rep(analysis_hash, n_resample),
    Run = seq_len(n_resample),
    Run.ID = sprintf("resample_%03d", seq_len(n_resample)),
    Seed = seeds,
    Permutation.Seed = rep(1L, n_resample),
    Status = rep("pending", n_resample),
    Result.File = file.path(
      sprintf("resample_%03d", seq_len(n_resample)),
      "run_result.rds"
    ),
    stringsAsFactors = FALSE
  )
  write_active_manifest <- function(){
    if (!is.null(outdir)){
      .InterMetaWriteTableAtomic(
        active_run_manifest,
        file = file.path(outdir, "active_run_manifest.txt")
      )
    }
  }
  write_active_manifest()
  for (i in seq_len(n_resample)){
    run_id <- sprintf("resample_%03d", i)
    run_dir <- if (is.null(outdir)) NULL else file.path(outdir, run_id)
    run_rds <- if (is.null(run_dir)) NULL else file.path(run_dir, "run_result.rds")
    persisted_result_file <- if (is.null(run_rds)) {
      NULL
    } else {
      active_run_manifest$Result.File[i]
    }
    permutation_seed <- 1L
    run_cache_key <- digest::digest(
      list(
        analysis_hash = analysis_hash,
        run = i,
        seed = seeds[i],
        permutation_seed = permutation_seed
      ),
      algo = "xxhash64"
    )
    if (!is.null(run_rds) && reuse && file.exists(run_rds)){
      cached <- tryCatch(readRDS(run_rds), error = function(e) NULL)
      if (!is.null(cached) &&
          identical(cached$cache_key, run_cache_key) &&
          identical(cached$status, "success")){
        message("Reusing validated result for ", run_id)
        if (!identical(cached$result_file, persisted_result_file)){
          cached$result_file <- persisted_result_file
          saveRDS(cached, run_rds)
        }
        stale_error_file <- file.path(run_dir, "error.txt")
        if (file.exists(stale_error_file)){
          unlink(stale_error_file)
        }
        runs[[i]] <- cached
        active_run_manifest$Status[i] <- "success"
        write_active_manifest()
        next
      }
      message("Ignoring stale or incompatible cache for ", run_id)
    }
    if (!is.null(run_dir)){
      dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    }

    selected_samples <- .InterMetaStratifiedSample(
      samples = common_samples,
      event = info[common_samples, 1],
      size = sample_size,
      seed = seeds[i],
      stratify = stratify_event,
      min_events = min_events,
      min_censored = min_censored
    )
    selected_samples <- sort(selected_samples)
    message("Running ", run_id, " with ", length(selected_samples),
            " samples (seed = ", seeds[i], ").")
    started_at <- Sys.time()
    event_count <- sum(info[selected_samples, 1] == 1)
    censored_count <- sum(info[selected_samples, 1] == 0)
    run <- tryCatch({
      if (event_count < min_events || censored_count < min_censored){
        stop("The selected subsample does not satisfy event/censoring minima.")
      }
      sub_info <- info[selected_samples, , drop = FALSE]
      sub_rna <- rna[, selected_samples, drop = FALSE]
      preprocessed <- InterMetaPreprocess(
        pathwayset = pathwayset,
        rna = sub_rna,
        cohort = row.names(sub_info),
        log = log,
        GMNorm = GMNorm,
        analysis_profile = analysis_profile
      )
      mps <- InterMetaGenePair(
        preprocessed,
        analysis_profile = analysis_profile
      )
      cox <- InterMetaCox(
        mps = mps,
        info = sub_info,
        p.adj = p.adj,
        parallel = parallel,
        core = core,
        keep_mps = FALSE,
        analysis_profile = analysis_profile
      )
      cox_res <- cox$cox_res
      rm(mps, cox)
      invisible(gc())
      selected_pair <- is.finite(cox_res[[pair_condition]]) &
        cox_res[[pair_condition]] < pair_cutoff
      query_pairs <- row.names(cox_res)[selected_pair]
      background <- InterMetaBG(row.names(cox_res))
      enrich <- InterMetaEnrich(
        gene = query_pairs,
        background = background,
        pathway.pair = pathway.pair,
        p.adj = enrich_p.adj,
        core = core,
        keep_gene_pairs = "none"
      )
      if (run_permutation){
        enrich <- InterMetaPermFDR(
          resdata = enrich,
          pathwayset = permutation_pathwayset,
          gene = query_pairs,
          background = background,
          niter = permutation_niter,
          core = core,
          seed = permutation_seed,
          analysis_profile = analysis_profile
        )
      }
      score <- InterMetaPathwayScore(
        resdata = enrich,
        condition = score_condition,
        cutoff = score_cutoff,
        analysis_profile = analysis_profile
      )
      network <- NULL
      if (run_network){
        network <- InterMetaPathwayNetwork(
          resdata = enrich,
          pathwayset = pathwayset,
          condition = score_condition,
          cutoff = score_cutoff,
          export = !is.null(run_dir),
          path = run_dir
        )
      }
      finished_at <- Sys.time()
      list(
        run = i,
        run_id = run_id,
        seed = seeds[i],
        permutation_seed = permutation_seed,
        cache_key = run_cache_key,
        result_file = persisted_result_file,
        status = "success",
        error_message = NULL,
        selected_samples = selected_samples,
        event_count = event_count,
        censored_count = censored_count,
        score = score,
        enrich = enrich,
        cox_res = cox_res,
        query_pairs = query_pairs,
        query_pair_count = length(query_pairs),
        network = network,
        started_at = started_at,
        finished_at = finished_at,
        elapsed_seconds = as.numeric(difftime(
          finished_at, started_at, units = "secs"
        ))
      )
    }, error = function(error){
      finished_at <- Sys.time()
      list(
        run = i,
        run_id = run_id,
        seed = seeds[i],
        permutation_seed = permutation_seed,
        cache_key = run_cache_key,
        result_file = persisted_result_file,
        status = "error",
        error_message = conditionMessage(error),
        selected_samples = selected_samples,
        event_count = event_count,
        censored_count = censored_count,
        score = .InterMetaEmptyScore(),
        enrich = NULL,
        cox_res = NULL,
        query_pairs = character(),
        query_pair_count = 0L,
        network = NULL,
        started_at = started_at,
        finished_at = finished_at,
        elapsed_seconds = as.numeric(difftime(
          finished_at, started_at, units = "secs"
        ))
      )
    })
    runs[[i]] <- run

    if (!is.null(run_dir)){
      utils::write.table(data.frame(Sample = selected_samples),
                         file = file.path(run_dir, "selected_samples.txt"),
                         quote = FALSE, sep = "\t", row.names = FALSE)
      utils::write.table(run$score,
                         file = file.path(run_dir, "score.txt"),
                         quote = FALSE, sep = "\t", row.names = FALSE)
      saveRDS(run, run_rds)
      error_file <- file.path(run_dir, "error.txt")
      if (run$status == "error"){
        writeLines(run$error_message, error_file)
      } else if (file.exists(error_file)){
        unlink(error_file)
      }
    }
    active_run_manifest$Status[i] <- run$status
    write_active_manifest()
    if (run$status == "error"){
      warning(run_id, " failed: ", run$error_message)
      if (!continue_on_error){
        stop(run$error_message)
      }
    }
  }

  candidate_pathways <- sort(unique(pathwayset$Pathway))
  stability <- list(
    runs = runs,
    pathway_pair = pathway.pair,
    candidate_pathways = candidate_pathways,
    parameters = parameters,
    provenance = list(
      analysis_hash = analysis_hash,
      package_version = tryCatch(
        as.character(utils::packageVersion("InterMeta")),
        error = function(e) "development"
      ),
      created_at = Sys.time(),
      session_info = utils::sessionInfo()
    )
  )
  summary <- InterMetaSummarizeStability(
    stability = stability,
    pathways = candidate_pathways
  )
  stability$summary <- summary$summary
  stability$pairwise <- summary$pairwise
  stability$top_set_definition <- summary$top_set_definition
  stability$edge_summary <- .InterMetaSummarizeEdges(
    runs = runs,
    condition = score_condition,
    cutoff = score_cutoff
  )
  stability$run_status <- .InterMetaRunStatus(runs)
  stability$active_run_manifest <- active_run_manifest
  successful_runs <- runs[vapply(runs, function(run){
    identical(run$status, "success")
  }, logical(1))]
  score_rows <- lapply(successful_runs, function(run){
    score <- run$score
    if (nrow(score) == 0L){
      return(NULL)
    }
    data.frame(
      Run = run$run,
      Run.ID = run$run_id,
      Seed = run$seed,
      score,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  score_rows <- score_rows[!vapply(score_rows, is.null, logical(1))]
  stability$score_long <- if (length(score_rows) > 0L) {
    do.call(rbind, score_rows)
  } else {
    data.frame(
      Run = integer(), Run.ID = character(), Seed = integer(),
      .InterMetaEmptyScore(),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  row.names(stability$score_long) <- NULL
  if (!keep_details){
    stability$runs <- lapply(stability$runs, .InterMetaLightRun)
  }
  class(stability) <- c("InterMetaStability", class(stability))

  if (!is.null(outdir)){
    utils::write.table(stability$summary,
                       file = file.path(outdir, "stability_summary.txt"),
                       quote = FALSE, sep = "\t", row.names = FALSE)
    utils::write.table(stability$pairwise,
                       file = file.path(outdir, "stability_pairwise.txt"),
                       quote = FALSE, sep = "\t", row.names = FALSE)
    utils::write.table(stability$edge_summary,
                       file = file.path(outdir, "stability_edges.txt"),
                       quote = FALSE, sep = "\t", row.names = FALSE)
    utils::write.table(stability$run_status,
                       file = file.path(outdir, "run_status.txt"),
                       quote = FALSE, sep = "\t", row.names = FALSE)
    .InterMetaWriteTableAtomic(
      stability$active_run_manifest,
      file = file.path(outdir, "active_run_manifest.txt")
    )
    utils::write.table(stability$score_long,
                       file = file.path(outdir, "stability_scores_long.txt"),
                       quote = FALSE, sep = "\t", row.names = FALSE)
    saveRDS(stability, file.path(outdir, "stability_result.rds"))
  }
  stability
}

print.InterMetaStability <- function(x, ...){
  cat("<InterMetaStability>\n")
  cat("  Runs:", nrow(x$run_status), "\n")
  cat("  Successful:", sum(x$run_status$Status == "success"), "\n")
  cat("  Failed:", sum(x$run_status$Status != "success"), "\n")
  cat("  Candidate pathways:", length(x$candidate_pathways), "\n")
  cat("  Recurrent interactions:", nrow(x$edge_summary), "\n")
  invisible(x)
}
