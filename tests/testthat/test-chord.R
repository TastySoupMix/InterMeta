make_chord_enrichment <- function(){
  data.frame(
    Index = 1:6,
    Favorable.Pathway = c("A", "B", "A", "C", "D", "A"),
    Unfavorable.Pathway = c("B", "C", "C", "A", "A", "B"),
    No.Shared.MGP = c(4, 3, 2, 1, 1, 1),
    Enrichment.Ratio = c(4, 3, 2, 1.5, 1, 2),
    P.Value = c(0, 0.002, 0.003, 0.004, 0.005, 0.8),
    Adjust.P.Value = c(0, 0.01, 0.02, 0.03, 0.04, 0.8),
    Permutation.FDR = c(0, 0.02, 0.03, 0.04, 0.05, 0.9),
    stringsAsFactors = FALSE
  )
}

test_that("pathway chord selection and colors are deterministic", {
  enrichment <- make_chord_enrichment()
  first <- InterMetaPathwayChord(
    enrichment,
    cutoff = 0.05,
    max_edges = 3,
    draw = FALSE
  )
  second <- InterMetaPathwayChord(
    enrichment,
    cutoff = 0.05,
    max_edges = 3,
    draw = FALSE
  )
  subset_chord <- InterMetaPathwayChord(
    enrichment,
    select = 1,
    draw = FALSE
  )

  expect_s3_class(first, "InterMetaPathwayChord")
  expect_equal(first$selected$Index, 1:3)
  expect_equal(first$edges$Favorable.Pathway, c("A", "B", "A"))
  expect_equal(first$edges$Unfavorable.Pathway, c("B", "C", "C"))
  expect_equal(first$edges$Weight, c(4, 3, 2))
  expect_equal(first$eligible_count, 5L)
  expect_equal(first$selected_count, 3L)
  expect_equal(first$displayed_count, 3L)
  expect_equal(first$color, second$color)
  expect_equal(first$color["A"], subset_chord$color["A"])
  expect_identical(names(first$color), first$order)
  expect_identical(names(first$group), first$order)
  expect_setequal(first$nodes$Pathway, c("A", "B", "C"))
})

test_that("pathway chord preserves direction and supports finite significance weights", {
  enrichment <- make_chord_enrichment()
  enrichment$Adjust.P.Value[6] <- 0.01
  chord <- InterMetaPathwayChord(
    enrichment,
    select = c(1, 6),
    link_weight = "significance",
    draw = FALSE
  )

  expect_equal(nrow(chord$selected), 2L)
  expect_equal(nrow(chord$edges), 1L)
  expect_equal(chord$edges$Favorable.Pathway, "A")
  expect_equal(chord$edges$Unfavorable.Pathway, "B")
  expect_true(is.finite(chord$edges$Weight))
  expect_equal(chord$significance_floor, 0.001)
  expect_equal(chord$edges$Weight, 5)
})

test_that("extra group, color, and order mappings are safely cropped", {
  enrichment <- make_chord_enrichment()
  pathway_group <- factor(
    c("source", "target", "target", "unused"),
    levels = c("source", "target", "unused")
  )
  names(pathway_group) <- c("A", "B", "C", "unused pathway")
  group_color <- c(source = "#1B9E77", target = "#D95F02", unused = "grey80")

  chord <- InterMetaPathwayChord(
    enrichment,
    select = 1:3,
    group = pathway_group,
    color = group_color,
    order = c("unused pathway", "A", "B", "C"),
    draw = FALSE
  )

  expect_equal(chord$order, c("A", "B", "C"))
  expect_equal(as.character(chord$group), c("source", "target", "target"))
  expect_equal(nlevels(chord$group), 2L)
  expect_equal(
    unname(chord$color),
    unname(group_color[c("source", "target", "target")])
  )
  expect_false("unused pathway" %in% names(chord$group))
})

test_that("empty and invalid chord selections are handled clearly", {
  enrichment <- make_chord_enrichment()
  expect_warning(
    empty <- InterMetaPathwayChord(
      enrichment,
      cutoff = 0,
      draw = FALSE
    ),
    "No pathway interactions were selected"
  )
  expect_s3_class(empty, "InterMetaPathwayChord")
  expect_equal(nrow(empty$edges), 0L)
  expect_equal(nrow(empty$nodes), 0L)

  expect_error(
    InterMetaPathwayChord(
      enrichment,
      select = 1:3,
      group = c(A = "one", B = "two"),
      draw = FALSE
    ),
    "missing selected pathway"
  )
  expect_error(
    InterMetaPathwayChord(
      enrichment,
      select = 1:3,
      color = c(A = "red", B = "blue"),
      draw = FALSE
    ),
    "cover every selected pathway or group"
  )
  expect_error(
    InterMetaPathwayChord(
      enrichment,
      adaptive_labels = NA,
      draw = FALSE
    ),
    "adaptive_labels"
  )
  expect_error(
    InterMetaPathwayChord(
      enrichment,
      canvas_limit = 1,
      draw = FALSE
    ),
    "canvas_limit"
  )
  expect_error(
    InterMetaPathwayChord(
      enrichment,
      sector_floor_quantile = 1.1,
      draw = FALSE
    ),
    "sector_floor_quantile"
  )
  observed_warnings <- character()
  enrichment$No.Shared.MGP[5] <- 0
  zero_weight <- withCallingHandlers(
    InterMetaPathwayChord(
      enrichment,
      select = 5,
      draw = FALSE
    ),
    warning = function(w){
      observed_warnings <<- c(observed_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("non-positive or non-finite", observed_warnings)))
  expect_true(any(grepl("No drawable", observed_warnings)))
  expect_equal(nrow(zero_weight$edges), 0L)
})

test_that("automatic truncation is applied after invalid links are removed", {
  enrichment <- make_chord_enrichment()
  enrichment$No.Shared.MGP[1] <- 0
  observed_warnings <- character()
  chord <- withCallingHandlers(
    InterMetaPathwayChord(
      enrichment,
      cutoff = 0.05,
      max_edges = 2,
      draw = FALSE
    ),
    warning = function(w){
      observed_warnings <<- c(observed_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(chord$eligible_count, 5L)
  expect_equal(chord$selected$Index, c(2L, 3L))
  expect_equal(chord$selected_count, 2L)
  expect_equal(chord$displayed_count, 2L)
  expect_true(any(grepl("non-positive or non-finite", observed_warnings)))
})

test_that("PRAD top interactions tolerate a mapping with unused factor levels", {
  data("InterMetaPRADEnrich", package = "InterMeta", envir = environment())
  all_pathways <- sort(unique(c(
    InterMetaPRADEnrich$Favorable.Pathway,
    InterMetaPRADEnrich$Unfavorable.Pathway
  )))
  pathway_group <- factor(
    rep("Metabolic pathways", length(all_pathways)),
    levels = c("Metabolic pathways", "Unused level")
  )
  names(pathway_group) <- all_pathways
  pathway_color <- stats::setNames(
    grDevices::hcl.colors(length(all_pathways), "Dark 3"),
    all_pathways
  )

  chord <- InterMetaPathwayChord(
    InterMetaPRADEnrich,
    condition = "Permutation.FDR",
    cutoff = 0.2,
    max_edges = 10,
    group = pathway_group,
    color = pathway_color,
    draw = FALSE
  )

  expect_equal(nrow(chord$edges), 10L)
  expect_equal(length(chord$order), 11L)
  expect_equal(nlevels(chord$group), 1L)
  expect_setequal(names(chord$color), chord$order)
})

test_that("PRAD pathway chord draws every interaction meeting the cutoff by default", {
  data("InterMetaPRADEnrich", package = "InterMeta", envir = environment())
  chord <- InterMetaPathwayChord(
    InterMetaPRADEnrich,
    condition = "Permutation.FDR",
    cutoff = 0.2,
    draw = FALSE
  )

  expect_null(chord$max_edges)
  expect_equal(chord$eligible_count, 104L)
  expect_equal(chord$selected_count, 104L)
  expect_equal(chord$displayed_count, 104L)
  expect_equal(nrow(chord$nodes), 39L)
})

test_that("gene chord expands one retained pathway interaction", {
  data("InterMetaExampleEnrich", package = "InterMeta", envir = environment())
  set.seed(123)
  rng_before <- .Random.seed
  first <- InterMetaGeneChord(
    InterMetaExampleEnrich,
    select = InterMetaExampleEnrich$Index[1],
    condition = "Permutation.FDR",
    draw = FALSE
  )
  rng_after <- .Random.seed
  second <- InterMetaGeneChord(
    InterMetaExampleEnrich,
    select = InterMetaExampleEnrich$Index[1],
    condition = "Permutation.FDR",
    draw = FALSE
  )

  expect_s3_class(first, "InterMetaGeneChord")
  expect_equal(nrow(first$edges), 25L)
  expect_equal(first$displayed_pair_count, 25L)
  expect_equal(first$reported_pair_count, 25L)
  expect_equal(first$edges$Pair[1], "LDHA_GCLC")
  expect_equal(first$edges$Source.Gene[1], "LDHA")
  expect_equal(first$edges$Target.Gene[1], "GCLC")
  expect_identical(first$color, second$color)
  expect_identical(rng_before, rng_after)
})

test_that("gene chord keeps source and target copies of a shared gene distinct", {
  enrichment <- data.frame(
    Index = 7L,
    Favorable.Pathway = "Source",
    Unfavorable.Pathway = "Target",
    No.Shared.MGP = 2L,
    Adjust.P.Value = 0.01,
    Shared.MGP = "A_A/A_B",
    stringsAsFactors = FALSE
  )
  chord <- InterMetaGeneChord(enrichment, select = 7, draw = FALSE)

  expect_true(all(c("F::A", "U::A") %in% chord$nodes$Sector))
  expect_equal(chord$nodes$Gene[chord$nodes$Sector == "F::A"], "A")
  expect_equal(chord$nodes$Gene[chord$nodes$Sector == "U::A"], "A")
  expect_false(
    chord$nodes$Side[chord$nodes$Sector == "F::A"] ==
      chord$nodes$Side[chord$nodes$Sector == "U::A"]
  )
})

test_that("gene chord reports invalid selections and omitted shared pairs clearly", {
  data("InterMetaPRADEnrich", package = "InterMeta", envir = environment())
  retained_index <- InterMetaPRADEnrich$Index[
    which(InterMetaPRADEnrich$No.Shared.MGP > 0)[1]
  ]

  expect_error(
    InterMetaGeneChord(InterMetaPRADEnrich, select = c(1, 2), draw = FALSE),
    "must be one"
  )
  expect_error(
    InterMetaGeneChord(InterMetaPRADEnrich, select = -1, draw = FALSE),
    "was not found"
  )
  expect_error(
    InterMetaGeneChord(
      InterMetaPRADEnrich,
      select = retained_index,
      condition = c("P.Value", "Adjust.P.Value"),
      draw = FALSE
    ),
    "condition.*must name one column"
  )
  expect_error(
    InterMetaGeneChord(
      InterMetaPRADEnrich,
      select = retained_index,
      draw = FALSE
    ),
    "keep_gene_pairs.*shared"
  )
})

test_that("gene chord draws to a graphics device and clears circlize state", {
  data("InterMetaExampleEnrich", package = "InterMeta", envir = environment())
  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output, width = 8, height = 8)
  device <- grDevices::dev.cur()
  on.exit({
    if (grDevices::dev.cur() == device) {
      grDevices::dev.off()
    }
    unlink(output)
  }, add = TRUE)

  chord <- InterMetaGeneChord(
    InterMetaExampleEnrich,
    select = InterMetaExampleEnrich$Index[1]
  )
  grDevices::dev.off()

  expect_s3_class(chord, "InterMetaGeneChord")
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
  expect_length(circlize::get.all.sector.index(), 0L)
})

test_that("pathway chord draws to a graphics device and clears circlize state", {
  enrichment <- make_chord_enrichment()
  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output, width = 8, height = 8)
  device <- grDevices::dev.cur()
  on.exit({
    if (grDevices::dev.cur() == device) {
      grDevices::dev.off()
    }
    unlink(output)
  }, add = TRUE)

  chord <- InterMetaPathwayChord(
    enrichment,
    select = 1:3,
    draw = TRUE
  )
  grDevices::dev.off()

  expect_s3_class(chord, "InterMetaPathwayChord")
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
  expect_length(circlize::get.all.sector.index(), 0L)
})

test_that("classic defaults keep proportional sectors and dense layout is opt-in", {
  defaults <- formals(InterMetaPathwayChord)
  expect_identical(defaults$adaptive_labels, FALSE)
  expect_null(defaults$sector_floor_quantile)
  expect_null(defaults$canvas_limit)
  expect_null(defaults$max_edges)
  expect_equal(defaults$label_cex, 0.7)
  expect_equal(defaults$start_degree, 90)
  pathways <- paste0("P", seq_len(40))
  nodes <- data.frame(
    Pathway = pathways,
    OutWeight = seq_len(40),
    InWeight = rev(seq_len(40))^2
  )
  layout <- InterMeta:::.InterMetaPathwayChordLayout(nodes, pathways)
  expect_equal(layout$canvas_limit, 1.45)
  expect_null(layout$sector_xmax)
  expect_equal(
    InterMeta:::.InterMetaPathwayChordLabelCex(pathways, 0.7, 1.45),
    0.7
  )
})

test_that("opt-in dense pathway chord preserves all links and fits long labels", {
  pathways <- sprintf(
    "Pathway %02d with a deliberately long metabolic pathway label",
    seq_len(52)
  )
  candidate <- expand.grid(
    Favorable.Pathway = pathways,
    Unfavorable.Pathway = pathways,
    stringsAsFactors = FALSE
  )
  candidate <- candidate[
    candidate$Favorable.Pathway != candidate$Unfavorable.Pathway,
    ,
    drop = FALSE
  ]
  candidate <- candidate[seq_len(663), , drop = FALSE]
  enrichment <- data.frame(
    Index = seq_len(nrow(candidate)),
    candidate,
    No.Shared.MGP = 1 + (seq_len(nrow(candidate)) %% 37)^2,
    Adjust.P.Value = rep(0.001, nrow(candidate)),
    stringsAsFactors = FALSE
  )

  prepared <- InterMetaPathwayChord(enrichment, draw = FALSE)
  adaptive_layout <- InterMeta:::.InterMetaPathwayChordLayout(
    nodes = prepared$nodes,
    order = prepared$order,
    canvas_limit = 1.65,
    sector_floor_quantile = 0.75
  )
  node_weight <- stats::setNames(
    prepared$nodes$OutWeight + prepared$nodes$InWeight,
    prepared$nodes$Pathway
  )[prepared$order]

  expect_equal(nrow(prepared$selected), 663L)
  expect_equal(nrow(prepared$edges), 663L)
  expect_equal(nrow(prepared$nodes), 52L)
  expect_equal(prepared$displayed_count, 663L)
  expect_identical(
    names(prepared),
    c(
      "selected", "edges", "nodes", "order", "group", "color",
      "condition", "cutoff", "selection", "eligible_count",
      "selected_count", "displayed_count", "max_edges", "link_weight",
      "significance_floor", "layout"
    )
  )
  expect_equal(adaptive_layout$canvas_limit, 1.65)
  expect_true(all(adaptive_layout$sector_xmax >= node_weight))
  expect_equal(
    min(adaptive_layout$sector_xmax),
    unname(stats::quantile(node_weight, 0.75, names = FALSE, type = 8))
  )

  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output, width = 9, height = 9)
  device <- grDevices::dev.cur()
  on.exit({
    if (grDevices::dev.cur() == device) {
      grDevices::dev.off()
    }
    unlink(output)
  }, add = TRUE)
  fitted_cex <- InterMeta:::.InterMetaPathwayChordLabelCex(
    labels = prepared$order,
    label_cex = 0.7,
    canvas_limit = adaptive_layout$canvas_limit,
    adaptive_labels = TRUE
  )
  radial_margin <- min(graphics::par("din")) / 2 *
    (adaptive_layout$canvas_limit - 1) /
    adaptive_layout$canvas_limit * 0.9
  fitted_width <- max(graphics::strwidth(
    prepared$order,
    units = "inches",
    cex = fitted_cex
  ))
  expect_lt(fitted_cex, 0.7)
  expect_lte(fitted_width, radial_margin + 1e-8)

  classic <- InterMetaPathwayChord(enrichment)
  drawn <- InterMetaPathwayChord(
    enrichment,
    adaptive_labels = TRUE,
    canvas_limit = 1.65,
    sector_floor_quantile = 0.75
  )
  grDevices::dev.off()
  expect_equal(classic$edges, drawn$edges)
  expect_equal(classic$nodes, drawn$nodes)
  expect_equal(nrow(classic$layout), 663L)
  expect_equal(nrow(drawn$edges), 663L)
  expect_equal(nrow(drawn$nodes), 52L)
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
  expect_length(circlize::get.all.sector.index(), 0L)
})

test_that("many pathway groups cannot exhaust the circular gap", {
  pathways <- paste0("P", seq_len(40))
  enrichment <- data.frame(
    Index = seq_len(40),
    Favorable.Pathway = pathways,
    Unfavorable.Pathway = c(pathways[-1], pathways[1]),
    No.Shared.MGP = rep(1, 40),
    Enrichment.Ratio = rep(2, 40),
    Adjust.P.Value = seq(0.001, 0.04, length.out = 40),
    stringsAsFactors = FALSE
  )
  unique_group <- factor(pathways, levels = c(pathways, "unused"))
  names(unique_group) <- pathways
  almost_unique_group <- factor(
    c("G1", "G1", paste0("G", 2:39)),
    levels = c(paste0("G", 1:39), "unused")
  )
  names(almost_unique_group) <- pathways

  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output, width = 8, height = 8)
  device <- grDevices::dev.cur()
  on.exit({
    if (grDevices::dev.cur() == device) {
      grDevices::dev.off()
    }
    unlink(output)
  }, add = TRUE)

  unique_chord <- InterMetaPathwayChord(
    enrichment,
    select = enrichment$Index,
    group = unique_group,
    label_cex = 0.3
  )
  almost_unique_chord <- InterMetaPathwayChord(
    enrichment,
    select = enrichment$Index,
    group = almost_unique_group,
    label_cex = 0.3
  )
  grDevices::dev.off()

  expect_s3_class(unique_chord, "InterMetaPathwayChord")
  expect_s3_class(almost_unique_chord, "InterMetaPathwayChord")
  expect_equal(nlevels(unique_chord$group), 40L)
  expect_equal(nlevels(almost_unique_chord$group), 39L)
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
  expect_length(circlize::get.all.sector.index(), 0L)
})
