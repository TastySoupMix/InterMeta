# Portions of this file were adapted from TimiGP 1.3.0:
# https://github.com/CSkylarL/TimiGP
# The upstream source identifies Chenyang Skylar Li as author.
# Modified for InterMeta; changes through 2026-09-13.
# SPDX-License-Identifier: GPL-3.0-only
# See inst/NOTICE and COPYING for provenance and license details.

.InterMetaPathwayChordLayout <- function(
    nodes,
    order,
    canvas_limit = NULL,
    sector_floor_quantile = NULL){
  dense <- length(order) >= 30L
  if (is.null(canvas_limit)){
    canvas_limit <- 1.45
  }

  sector_xmax <- NULL
  if (dense && !is.null(sector_floor_quantile)){
    node_weight <- stats::setNames(
      nodes$OutWeight + nodes$InWeight,
      nodes$Pathway
    )[order]
    weight_floor <- unname(stats::quantile(
      node_weight,
      probs = sector_floor_quantile,
      names = FALSE,
      type = 8
    ))
    sector_xmax <- stats::setNames(
      pmax(node_weight, weight_floor),
      order
    )
  }

  list(
    dense = dense,
    canvas_limit = canvas_limit,
    sector_xmax = sector_xmax
  )
}

.InterMetaPathwayChordLabelCex <- function(
    labels,
    label_cex,
    canvas_limit,
    adaptive_labels = FALSE){
  if (!adaptive_labels || length(labels) < 30L){
    return(label_cex)
  }

  density_cex <- 0.7 * sqrt(40 / length(labels))
  effective_cex <- min(label_cex, density_cex)
  device_size <- suppressWarnings(min(graphics::par("din")))
  label_width <- suppressWarnings(max(graphics::strwidth(
    labels,
    units = "inches",
    cex = 1
  )))
  if (is.finite(device_size) && device_size > 0 &&
      is.finite(label_width) && label_width > 0){
    radial_margin <- device_size / 2 *
      (canvas_limit - 1) / canvas_limit * 0.9
    effective_cex <- min(effective_cex, radial_margin / label_width)
  }
  max(effective_cex, .Machine$double.eps)
}

InterMetaPathwayChord <- function(
    resdata = NULL,
    select = NULL,
    condition = "Adjust.P.Value",
    cutoff = 0.05,
    max_edges = NULL,
    link_weight = c("No.Shared.MGP", "Enrichment.Ratio", "significance",
                    "uniform"),
    group = NULL,
    color = NULL,
    order = NULL,
    transparency = 0.4,
    label_cex = 0.7,
    start_degree = 90,
    adaptive_labels = FALSE,
    canvas_limit = NULL,
    sector_floor_quantile = NULL,
    draw = TRUE){
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  if (!is.data.frame(resdata)){
    stop('The parameter "resdata" must be a data frame.')
  }
  if (!is.character(condition) || length(condition) != 1L ||
      is.na(condition) || !nzchar(condition)){
    stop('The parameter "condition" must name one column.')
  }
  required <- c("Index", "Favorable.Pathway", "Unfavorable.Pathway",
                condition)
  missing <- setdiff(required, names(resdata))
  if (length(missing) > 0L){
    stop("The enrichment result is missing required column(s): ",
         paste(missing, collapse = ", "))
  }
  if (!is.numeric(resdata[[condition]])){
    stop("The selection column must be numeric: ", condition)
  }
  if (length(cutoff) != 1L || !is.numeric(cutoff) ||
      !is.finite(cutoff) || cutoff < 0 || cutoff > 1){
    stop('The parameter "cutoff" must be a finite number from 0 to 1.')
  }
  if (anyNA(resdata$Index) || anyDuplicated(resdata$Index)){
    stop('The enrichment column "Index" must contain unique, non-missing values.')
  }
  if (!is.null(max_edges) &&
      (length(max_edges) != 1L || !is.numeric(max_edges) ||
       is.na(max_edges) || max_edges < 1 ||
       (is.finite(max_edges) && max_edges != floor(max_edges)))){
    stop('The parameter "max_edges" must be NULL, Inf, or a positive integer.')
  }
  if (length(transparency) != 1L || !is.numeric(transparency) ||
      !is.finite(transparency) || transparency < 0 || transparency > 1){
    stop('The parameter "transparency" must be between 0 and 1.')
  }
  if (length(label_cex) != 1L || !is.numeric(label_cex) ||
      !is.finite(label_cex) || label_cex <= 0){
    stop('The parameter "label_cex" must be a positive finite number.')
  }
  if (length(start_degree) != 1L || !is.numeric(start_degree) ||
      !is.finite(start_degree)){
    stop('The parameter "start_degree" must be a finite number.')
  }
  if (length(adaptive_labels) != 1L || is.na(adaptive_labels) ||
      !is.logical(adaptive_labels)){
    stop('The parameter "adaptive_labels" must be TRUE or FALSE.')
  }
  if (!is.null(canvas_limit) &&
      (length(canvas_limit) != 1L || !is.numeric(canvas_limit) ||
       !is.finite(canvas_limit) || canvas_limit <= 1)){
    stop('The parameter "canvas_limit" must be NULL or a finite number greater than 1.')
  }
  if (!is.null(sector_floor_quantile) &&
      (length(sector_floor_quantile) != 1L ||
       !is.numeric(sector_floor_quantile) ||
       !is.finite(sector_floor_quantile) ||
       sector_floor_quantile < 0 || sector_floor_quantile > 1)){
    stop('The parameter "sector_floor_quantile" must be NULL or a finite number from 0 to 1.')
  }
  if (length(draw) != 1L || is.na(draw) || !is.logical(draw)){
    stop('The parameter "draw" must be TRUE or FALSE.')
  }
  link_weight <- match.arg(link_weight)

  selection_mode <- if (is.null(select)) "condition" else "index"
  eligible_count <- 0L
  significance_floor <- NA_real_
  if (is.null(select)){
    keep <- which(is.finite(resdata[[condition]]) &
                    resdata[[condition]] < cutoff)
    eligible_count <- length(keep)
    if (length(keep) > 0L){
      shared <- if ("No.Shared.MGP" %in% names(resdata)) {
        suppressWarnings(as.numeric(resdata$No.Shared.MGP[keep]))
      } else {
        rep(0, length(keep))
      }
      shared[!is.finite(shared)] <- 0
      keep <- keep[order(resdata[[condition]][keep], -shared, keep)]
    }
  } else {
    select <- unique(select[!is.na(select)])
    unknown <- setdiff(as.character(select), as.character(resdata$Index))
    if (length(unknown) > 0L){
      warning("Ignoring unknown enrichment Index value(s): ",
              paste(unknown, collapse = ", "))
    }
    keep <- match(as.character(select), as.character(resdata$Index),
                  nomatch = 0L)
    keep <- keep[keep > 0L]
    eligible_count <- length(keep)
  }
  selected <- resdata[keep, , drop = FALSE]

  empty_result <- function(selected_rows){
    edges <- data.frame(
      Favorable.Pathway = character(),
      Unfavorable.Pathway = character(),
      Weight = numeric(),
      No.Interactions = integer(),
      stringsAsFactors = FALSE
    )
    nodes <- data.frame(
      Pathway = character(), Group = character(), Color = character(),
      OutDegree = integer(), InDegree = integer(),
      OutWeight = numeric(), InWeight = numeric(),
      stringsAsFactors = FALSE
    )
    structure(
      list(
        selected = selected_rows,
        edges = edges,
        nodes = nodes,
        order = character(),
        group = structure(character(), names = character()),
        color = structure(character(), names = character()),
        condition = condition,
        cutoff = cutoff,
        selection = selection_mode,
        eligible_count = as.integer(eligible_count),
        selected_count = nrow(selected_rows),
        displayed_count = 0L,
        max_edges = max_edges,
        link_weight = link_weight,
        significance_floor = significance_floor,
        layout = NULL
      ),
      class = "InterMetaPathwayChord"
    )
  }

  if (nrow(selected) == 0L){
    warning("No pathway interactions were selected; no chord diagram was drawn.")
    return(invisible(empty_result(selected)))
  }

  source <- trimws(as.character(selected$Favorable.Pathway))
  target <- trimws(as.character(selected$Unfavorable.Pathway))
  valid_pathway <- !is.na(source) & nzchar(source) &
    !is.na(target) & nzchar(target)
  if (any(!valid_pathway)){
    warning(sum(!valid_pathway),
            " selected interaction(s) with a missing pathway name were ignored.")
  }
  selected <- selected[valid_pathway, , drop = FALSE]
  source <- source[valid_pathway]
  target <- target[valid_pathway]

  self_loop <- source == target
  if (any(self_loop)){
    warning(sum(self_loop), " self-interaction(s) were ignored.")
    selected <- selected[!self_loop, , drop = FALSE]
    source <- source[!self_loop]
    target <- target[!self_loop]
  }
  if (nrow(selected) == 0L){
    warning("No drawable pathway interactions remained; no chord diagram was drawn.")
    return(invisible(empty_result(selected)))
  }

  if (link_weight == "uniform"){
    weight <- rep(1, nrow(selected))
  } else if (link_weight == "significance"){
    probability <- selected[[condition]]
    invalid_probability <- !is.finite(probability) |
      probability < 0 | probability > 1
    if (any(invalid_probability)){
      stop('Selected rows must contain finite probabilities from 0 to 1 in "',
           condition, '" when link_weight = "significance".')
    }
    positive_probability <- probability[probability > 0]
    significance_floor <- if (length(positive_probability) > 0L) {
      max(min(positive_probability) / 10, .Machine$double.xmin)
    } else {
      .Machine$double.eps
    }
    weight <- -log10(pmax(probability, significance_floor))
    weight <- pmax(weight, .Machine$double.eps)
  } else {
    if (!link_weight %in% names(selected)){
      stop("The enrichment result does not contain link-weight column: ",
           link_weight)
    }
    weight <- suppressWarnings(as.numeric(selected[[link_weight]]))
  }
  valid_weight <- is.finite(weight) & weight > 0
  if (any(!valid_weight)){
    warning(sum(!valid_weight),
            " selected interaction(s) with a non-positive or non-finite ",
            "link weight were ignored.")
    selected <- selected[valid_weight, , drop = FALSE]
    source <- source[valid_weight]
    target <- target[valid_weight]
    weight <- weight[valid_weight]
  }
  if (nrow(selected) == 0L){
    warning("No drawable pathway interactions remained; no chord diagram was drawn.")
    return(invisible(empty_result(selected)))
  }

  if (selection_mode == "condition" && !is.null(max_edges) &&
      is.finite(max_edges) && nrow(selected) > max_edges){
    message(
      eligible_count, " interactions met ", condition, " < ", cutoff,
      "; drawing the top ", as.integer(max_edges),
      ". Set max_edges = NULL to draw all."
    )
    keep_display <- seq_len(as.integer(max_edges))
    selected <- selected[keep_display, , drop = FALSE]
    source <- source[keep_display]
    target <- target[keep_display]
    weight <- weight[keep_display]
  }

  key <- paste(source, target, sep = "\r")
  unique_key <- unique(key)
  first <- match(unique_key, key)
  edges <- data.frame(
    Favorable.Pathway = source[first],
    Unfavorable.Pathway = target[first],
    Weight = vapply(unique_key, function(x) sum(weight[key == x]), numeric(1)),
    No.Interactions = vapply(unique_key, function(x) sum(key == x), integer(1)),
    stringsAsFactors = FALSE
  )
  pathways <- unique(as.vector(t(as.matrix(edges[, c(
    "Favorable.Pathway", "Unfavorable.Pathway"
  )]))))

  group_names <- names(group)
  if (is.null(group)){
    group <- rep("Metabolic pathways", length(pathways))
    names(group) <- pathways
  } else {
    group <- as.character(group)
    names(group) <- group_names
    if (is.null(names(group))){
      if (length(group) != length(pathways)){
        stop('An unnamed "group" must have one value per selected pathway.')
      }
      names(group) <- pathways
    } else {
      if (any(!nzchar(names(group))) || anyDuplicated(names(group))){
        stop('A named "group" must have unique, non-empty pathway names.')
      }
      missing_group <- setdiff(pathways, names(group))
      if (length(missing_group) > 0L){
        stop('The "group" mapping is missing selected pathway(s): ',
             paste(missing_group, collapse = ", "))
      }
      group <- group[pathways]
    }
    if (anyNA(group) || any(!nzchar(trimws(group)))){
      stop('The "group" mapping contains a missing or empty group label.')
    }
  }
  group <- factor(group, levels = unique(as.character(group)))
  names(group) <- pathways

  base_order <- pathways
  if (is.null(order)){
    group_order <- unique(as.character(group[base_order]))
    order <- unlist(lapply(group_order, function(x){
      base_order[as.character(group[base_order]) == x]
    }), use.names = FALSE)
  } else {
    order <- as.character(order)
    if (anyNA(order) || any(!nzchar(order))){
      stop('The parameter "order" cannot contain missing or empty names.')
    }
    selected_order <- order[order %in% pathways]
    if (anyDuplicated(selected_order)){
      stop('The parameter "order" contains duplicated selected pathways.')
    }
    missing_order <- setdiff(pathways, selected_order)
    if (length(missing_order) > 0L){
      stop('The parameter "order" is missing selected pathway(s): ',
           paste(missing_order, collapse = ", "))
    }
    order <- selected_order
    group_runs <- rle(as.character(group[order]))$values
    if (anyDuplicated(group_runs)){
      stop('The parameter "order" must keep pathways from each group contiguous.')
    }
  }

  color_names <- names(color)
  if (is.null(color)){
    color <- vapply(pathways, function(pathway){
      hash <- digest::digest(enc2utf8(pathway), algo = "xxhash32",
                             serialize = FALSE)
      value <- strtoi(substr(hash, 1, 7), base = 16L)
      grDevices::hcl(
        h = value %% 360,
        c = 65 + (value %/% 360) %% 16,
        l = 60 + (value %/% (360 * 16)) %% 11,
        fixup = TRUE
      )
    }, character(1))
    names(color) <- pathways
  } else {
    color <- as.character(color)
    names(color) <- color_names
    group_levels <- unique(as.character(group[pathways]))
    if (is.null(names(color))){
      if (length(color) == length(pathways)){
        names(color) <- pathways
      } else if (length(color) == length(group_levels)){
        color <- stats::setNames(color, group_levels)[as.character(group)]
        names(color) <- pathways
      } else {
        stop('An unnamed "color" must have one value per selected pathway ',
             "or one value per selected group.")
      }
    } else {
      if (any(!nzchar(names(color))) || anyDuplicated(names(color))){
        stop('A named "color" must have unique, non-empty names.')
      }
      if (all(pathways %in% names(color))){
        color <- color[pathways]
      } else if (all(group_levels %in% names(color))){
        color <- color[as.character(group)]
        names(color) <- pathways
      } else {
        missing_color <- setdiff(pathways, names(color))
        stop('The "color" mapping must cover every selected pathway or group',
             if (length(missing_color) > 0L) paste0(": ",
               paste(missing_color, collapse = ", ")) else ".")
      }
    }
    if (anyNA(color) || any(!nzchar(color))){
      stop('The "color" mapping contains a missing or empty color.')
    }
    tryCatch(
      grDevices::col2rgb(unname(color)),
      error = function(e){
        stop('The "color" mapping contains an invalid R color: ',
             conditionMessage(e), call. = FALSE)
      }
    )
  }
  color <- color[pathways]

  out_degree <- table(factor(edges$Favorable.Pathway, levels = order))
  in_degree <- table(factor(edges$Unfavorable.Pathway, levels = order))
  out_weight <- vapply(order, function(x){
    sum(edges$Weight[edges$Favorable.Pathway == x])
  }, numeric(1))
  in_weight <- vapply(order, function(x){
    sum(edges$Weight[edges$Unfavorable.Pathway == x])
  }, numeric(1))
  nodes <- data.frame(
    Pathway = order,
    Group = as.character(group[order]),
    Color = unname(color[order]),
    OutDegree = as.integer(out_degree),
    InDegree = as.integer(in_degree),
    OutWeight = unname(out_weight),
    InWeight = unname(in_weight),
    stringsAsFactors = FALSE
  )
  adaptive_layout <- .InterMetaPathwayChordLayout(
    nodes = nodes,
    order = order,
    canvas_limit = canvas_limit,
    sector_floor_quantile = sector_floor_quantile
  )

  result <- structure(
    list(
      selected = selected,
      edges = edges,
      nodes = nodes,
      order = order,
      group = group[order],
      color = color[order],
      condition = condition,
      cutoff = cutoff,
      selection = selection_mode,
      eligible_count = as.integer(eligible_count),
      selected_count = nrow(selected),
      displayed_count = nrow(edges),
      max_edges = max_edges,
      link_weight = link_weight,
      significance_floor = significance_floor,
      layout = NULL
    ),
    class = "InterMetaPathwayChord"
  )

  if (draw){
    circlize::circos.clear()
    on.exit(circlize::circos.clear(), add = TRUE)
    circlize::circos.par(
      start.degree = start_degree,
      canvas.xlim = c(-adaptive_layout$canvas_limit,
                      adaptive_layout$canvas_limit),
      canvas.ylim = c(-adaptive_layout$canvas_limit,
                      adaptive_layout$canvas_limit),
      points.overflow.warning = FALSE
    )
    effective_label_cex <- .InterMetaPathwayChordLabelCex(
      labels = order,
      label_cex = label_cex,
      canvas_limit = adaptive_layout$canvas_limit,
      adaptive_labels = adaptive_labels
    )
    plot_data <- data.frame(
      from = edges$Favorable.Pathway,
      to = edges$Unfavorable.Pathway,
      value = edges$Weight,
      stringsAsFactors = FALSE
    )
    chord_arguments <- list(
      x = plot_data,
      grid.col = color[order],
      transparency = transparency,
      order = order,
      directional = 1,
      direction.type = c("diffHeight", "arrows"),
      diffHeight = 0.04,
      link.arr.type = "big.arrow",
      link.arr.length = 0.04,
      link.sort = TRUE,
      link.decreasing = TRUE,
      link.zindex = rank(edges$Weight, ties.method = "first"),
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = 0.16)
    )
    if (!is.null(adaptive_layout$sector_xmax)){
      chord_arguments$xmax <- adaptive_layout$sector_xmax
    }
    n_groups <- length(unique(as.character(group[order])))
    if (n_groups > 1L && n_groups < length(order)){
      plot_group <- factor(
        as.character(group[order]),
        levels = unique(as.character(group[order]))
      )
      names(plot_group) <- order
      chord_arguments$group <- plot_group
      small_gap <- min(1, 180 / length(order))
      available_group_gap <- 180 -
        small_gap * (length(order) - n_groups)
      chord_arguments$big.gap <- min(10, available_group_gap / n_groups)
      chord_arguments$small.gap <- small_gap
    }
    result$layout <- do.call(circlize::chordDiagram, chord_arguments)
    circlize::circos.trackPlotRegion(
      track.index = 1,
      bg.border = NA,
      panel.fun = function(x, y){
        sector <- circlize::get.cell.meta.data("sector.index")
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        circlize::circos.text(
          mean(xlim), mean(ylim), sector,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = effective_label_cex
        )
      }
    )
  }

  invisible(result)
}

InterMetaGeneChord <- function(
    resdata = NULL,
    select = 1L,
    condition = "Adjust.P.Value",
    color = NULL,
    transparency = 0.4,
    label_cex = 0.8,
    start_degree = 180,
    draw = TRUE){
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  if (!is.data.frame(resdata)){
    stop('The parameter "resdata" must be a data frame.')
  }
  if (!is.character(condition) || length(condition) != 1L ||
      is.na(condition) || !nzchar(condition)){
    stop('The parameter "condition" must name one column.')
  }
  required <- c("Index", "Favorable.Pathway", "Unfavorable.Pathway",
                "Shared.MGP", condition)
  missing <- setdiff(required, names(resdata))
  if (length(missing) > 0L){
    stop("The enrichment result is missing required column(s): ",
         paste(missing, collapse = ", "))
  }
  if (length(select) != 1L || is.na(select)){
    stop('The parameter "select" must be one non-missing enrichment Index.')
  }
  if (anyNA(resdata$Index) || anyDuplicated(resdata$Index)){
    stop('The enrichment column "Index" must contain unique, non-missing values.')
  }
  selected_position <- match(as.character(select), as.character(resdata$Index))
  if (is.na(selected_position)){
    stop("The selected enrichment Index was not found: ", select)
  }
  if (!is.numeric(resdata[[condition]])){
    stop("The selection column must be numeric: ", condition)
  }
  if (length(transparency) != 1L || !is.numeric(transparency) ||
      !is.finite(transparency) || transparency < 0 || transparency > 1){
    stop('The parameter "transparency" must be between 0 and 1.')
  }
  if (length(label_cex) != 1L || !is.numeric(label_cex) ||
      !is.finite(label_cex) || label_cex <= 0){
    stop('The parameter "label_cex" must be a positive finite number.')
  }
  if (length(start_degree) != 1L || !is.numeric(start_degree) ||
      !is.finite(start_degree)){
    stop('The parameter "start_degree" must be a finite number.')
  }
  if (length(draw) != 1L || is.na(draw) || !is.logical(draw)){
    stop('The parameter "draw" must be TRUE or FALSE.')
  }

  selected <- resdata[selected_position, , drop = FALSE]
  shared_text <- trimws(as.character(selected$Shared.MGP[1]))
  retained_count <- if ("No.Shared.MGP" %in% names(selected)) {
    suppressWarnings(as.numeric(selected$No.Shared.MGP[1]))
  } else {
    NA_real_
  }
  if (is.na(shared_text) || !nzchar(shared_text)){
    if (is.finite(retained_count) && retained_count > 0){
      stop(
        "The selected enrichment row does not retain Shared.MGP. ",
        "Re-run InterMetaEnrich(..., keep_gene_pairs = \"shared\")."
      )
    }
    stop("The selected pathway interaction contains no shared gene pairs.")
  }
  pair_tokens <- trimws(strsplit(shared_text, "/", fixed = TRUE)[[1]])
  if (length(pair_tokens) == 0L || any(!nzchar(pair_tokens))){
    stop('The selected "Shared.MGP" value contains an empty gene-pair token.')
  }
  parsed <- .InterMetaParseGenePairs(pair_tokens)
  if (nrow(parsed) == 0L){
    stop("The selected pathway interaction contains no shared gene pairs.")
  }
  if (is.finite(retained_count) && retained_count != nrow(parsed)){
    warning(
      'The selected row reports No.Shared.MGP = ', retained_count,
      ", but ", nrow(parsed), " unique pair(s) were parsed from Shared.MGP."
    )
  }

  parsed$Source.Sector <- paste0("F::", parsed$Source.Gene)
  parsed$Target.Sector <- paste0("U::", parsed$Target.Gene)
  parsed$Weight <- 1
  source_sector <- unique(parsed$Source.Sector)
  target_sector <- unique(parsed$Target.Sector)
  order <- c(source_sector, target_sector)
  node_label <- c(
    stats::setNames(unique(parsed$Source.Gene), source_sector),
    stats::setNames(unique(parsed$Target.Gene), target_sector)
  )
  group <- c(
    stats::setNames(rep("Favorable pathway", length(source_sector)),
                    source_sector),
    stats::setNames(rep("Unfavorable pathway", length(target_sector)),
                    target_sector)
  )
  group <- factor(group[order],
                  levels = c("Favorable pathway", "Unfavorable pathway"))
  names(group) <- order

  color_names <- names(color)
  if (is.null(color)){
    source_color <- vapply(source_sector, function(sector){
      hash <- digest::digest(enc2utf8(sector), algo = "xxhash32",
                             serialize = FALSE)
      value <- strtoi(substr(hash, 1, 7), base = 16L)
      grDevices::hcl(
        h = 5 + value %% 46,
        c = 60 + (value %/% 46) %% 16,
        l = 58 + (value %/% (46 * 16)) %% 13,
        fixup = TRUE
      )
    }, character(1))
    target_color <- vapply(target_sector, function(sector){
      hash <- digest::digest(enc2utf8(sector), algo = "xxhash32",
                             serialize = FALSE)
      value <- strtoi(substr(hash, 1, 7), base = 16L)
      grDevices::hcl(
        h = 205 + value %% 66,
        c = 55 + (value %/% 66) %% 16,
        l = 58 + (value %/% (66 * 16)) %% 13,
        fixup = TRUE
      )
    }, character(1))
    color <- c(stats::setNames(source_color, source_sector),
               stats::setNames(target_color, target_sector))
  } else {
    color <- as.character(color)
    names(color) <- color_names
    if (is.null(names(color))){
      if (length(color) != length(order)){
        stop('An unnamed "color" must have one value per gene sector.')
      }
      names(color) <- order
    } else {
      if (any(!nzchar(names(color))) || anyDuplicated(names(color))){
        stop('A named "color" must have unique, non-empty sector names.')
      }
      missing_color <- setdiff(order, names(color))
      if (length(missing_color) > 0L){
        stop(
          'The "color" mapping is missing gene sector(s): ',
          paste(missing_color, collapse = ", "),
          '. Use names such as "F::GENE" and "U::GENE".'
        )
      }
      color <- color[order]
    }
    if (anyNA(color) || any(!nzchar(color))){
      stop('The "color" mapping contains a missing or empty color.')
    }
    tryCatch(
      grDevices::col2rgb(unname(color)),
      error = function(e){
        stop('The "color" mapping contains an invalid R color: ',
             conditionMessage(e), call. = FALSE)
      }
    )
  }
  color <- color[order]

  nodes <- data.frame(
    Sector = order,
    Gene = unname(node_label[order]),
    Side = as.character(group[order]),
    Color = unname(color[order]),
    Degree = as.integer(vapply(order, function(sector){
      sum(parsed$Source.Sector == sector | parsed$Target.Sector == sector)
    }, integer(1))),
    stringsAsFactors = FALSE
  )
  result <- structure(
    list(
      selected = selected,
      edges = parsed[, c("Pair", "Source.Gene", "Target.Gene",
                         "Source.Sector", "Target.Sector", "Weight")],
      nodes = nodes,
      order = order,
      group = group,
      color = color,
      condition = condition,
      condition_value = selected[[condition]][1],
      pathway_interaction = paste(
        selected$Favorable.Pathway[1], selected$Unfavorable.Pathway[1],
        sep = " -> "
      ),
      reported_pair_count = retained_count,
      displayed_pair_count = nrow(parsed),
      layout = NULL
    ),
    class = "InterMetaGeneChord"
  )

  if (draw){
    circlize::circos.clear()
    on.exit(circlize::circos.clear(), add = TRUE)
    circlize::circos.par(
      start.degree = start_degree,
      canvas.xlim = c(-1.35, 1.35),
      canvas.ylim = c(-1.35, 1.35),
      points.overflow.warning = FALSE
    )
    plot_group <- factor(as.character(group[order]),
                         levels = levels(group))
    names(plot_group) <- order
    plot_data <- data.frame(
      from = parsed$Source.Sector,
      to = parsed$Target.Sector,
      value = parsed$Weight,
      stringsAsFactors = FALSE
    )
    result$layout <- circlize::chordDiagram(
      x = plot_data,
      grid.col = color[order],
      transparency = transparency,
      order = order,
      group = plot_group,
      big.gap = 10,
      small.gap = min(1, 90 / length(order)),
      directional = 1,
      direction.type = c("diffHeight", "arrows"),
      diffHeight = 0.04,
      link.arr.type = "big.arrow",
      link.arr.length = 0.04,
      link.sort = TRUE,
      link.decreasing = TRUE,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = 0.16)
    )
    circlize::circos.trackPlotRegion(
      track.index = 1,
      bg.border = NA,
      panel.fun = function(x, y){
        sector <- circlize::get.cell.meta.data("sector.index")
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        circlize::circos.text(
          mean(xlim), mean(ylim), node_label[[sector]],
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = label_cex
        )
      }
    )
    metric <- result$condition_value
    metric_label <- if (is.finite(metric)) {
      format(metric, digits = 3, scientific = metric < 0.001)
    } else {
      "NA"
    }
    graphics::title(
      main = paste0(
        result$pathway_interaction, "\n", condition, " = ", metric_label
      ),
      line = -1,
      cex.main = 0.9
    )
  }

  invisible(result)
}

print.InterMetaGeneChord <- function(x, ...){
  cat("<InterMetaGeneChord>\n")
  cat("  Pathway interaction:", x$pathway_interaction, "\n")
  cat("  Shared gene pairs:", nrow(x$edges), "\n")
  cat("  Favorable-side genes:",
      sum(x$nodes$Side == "Favorable pathway"), "\n")
  cat("  Unfavorable-side genes:",
      sum(x$nodes$Side == "Unfavorable pathway"), "\n")
  cat("  Selection metric:", x$condition, "=", x$condition_value, "\n")
  invisible(x)
}

print.InterMetaPathwayChord <- function(x, ...){
  cat("<InterMetaPathwayChord>\n")
  cat("  Eligible enrichment rows:", x$eligible_count, "\n")
  cat("  Selected enrichment rows:", nrow(x$selected), "\n")
  cat("  Drawable directional links:", nrow(x$edges), "\n")
  cat("  Pathways:", nrow(x$nodes), "\n")
  cat("  Selection:", x$selection, "\n")
  cat("  Link weight:", x$link_weight, "\n")
  invisible(x)
}
