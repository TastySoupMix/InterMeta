make_expanded_pairs <- function(pathwayset) {
  pathways <- unique(pathwayset$Pathway)
  grid <- expand.grid(
    Source.Pathway = pathways,
    Target.Pathway = pathways,
    stringsAsFactors = FALSE
  )
  grid <- grid[grid$Source.Pathway != grid$Target.Pathway, , drop = FALSE]
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    source <- grid$Source.Pathway[i]
    target <- grid$Target.Pathway[i]
    source_genes <- pathwayset$Gene[pathwayset$Pathway == source]
    target_genes <- pathwayset$Gene[pathwayset$Pathway == target]
    pair <- as.vector(outer(source_genes, target_genes, paste, sep = "_"))
    pair <- pair[
      sub("_.*$", "", pair) != sub("^.*_", "", pair)
    ]
    data.frame(
      Pathway.Pair = paste(source, target, sep = "__"),
      Source.Pathway = source,
      Target.Pathway = target,
      Metabolic.Gene.Pair = unique(pair),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

test_that("compact enrichment matches expanded legacy annotation", {
  geneset <- data.frame(
    Pathway = c("A", "A", "B", "B", "C", "C"),
    Gene = c("G1", "G2", "G3", "G4", "G2", "G5"),
    Dataset = "toy",
    stringsAsFactors = FALSE
  )
  genes <- paste0("G", 1:5)
  background <- as.vector(outer(genes, genes, paste, sep = "_"))
  background <- background[
    sub("_.*$", "", background) != sub("^.*_", "", background)
  ]
  query <- c("G1_G3", "G2_G4", "G3_G1", "G5_G4")

  compact <- InterMetaEnrich(
    gene = query,
    background = background,
    pathway.pair = InterMetaPathwayPair(geneset),
    p.adj = "BH",
    keep_gene_pairs = "none"
  )
  expanded <- InterMetaEnrich(
    gene = query,
    background = background,
    pathway.pair = make_expanded_pairs(geneset),
    p.adj = "BH",
    keep_gene_pairs = "none"
  )
  compact <- compact[order(compact$Pathway.Interaction), ]
  expanded <- expanded[order(expanded$Pathway.Interaction), ]
  expect_equal(compact$Pathway.Interaction, expanded$Pathway.Interaction)
  expect_equal(compact$No.Total.MGP, expanded$No.Total.MGP)
  expect_equal(compact$No.Shared.MGP, expanded$No.Shared.MGP)
  expect_equal(compact$P.Value, expanded$P.Value, tolerance = 1e-12)
})

test_that("balanced and empty pathway scores are handled", {
  enrichment <- data.frame(
    Index = 1:2,
    Rank = 1:2,
    Pathway.Interaction = c("A__B", "B__A"),
    Favorable.Pathway = c("A", "B"),
    Unfavorable.Pathway = c("B", "A"),
    No.Total.MGP = c(2, 2),
    No.Shared.MGP = c(1, 1),
    Enrichment.Ratio = c(2, 2),
    P.Value = c(0.01, 0.01),
    Adjust.P.Value = c(0.01, 0.01),
    Shared.MGP = c("", ""),
    Total.MGP = c(NA, NA)
  )
  expect_warning(score <- InterMetaPathwayScore(enrichment), "neutral")
  expect_equal(nrow(score), 0L)

  directional <- InterMetaPathwayScore(enrichment[1, , drop = FALSE])
  expect_equal(directional$Rank, c(1, 2))
  expect_equal(directional$Group, c("F", "U"))

  enrichment$Adjust.P.Value <- 1
  empty <- suppressWarnings(InterMetaPathwayScore(enrichment))
  expect_equal(nrow(empty), 0)
  expect_named(
    empty,
    c("Index", "Rank", "Signed.Rank", "Pathway", "OutDegree",
      "Favorable.Score", "InDegree", "Unfavorable.Score", "Difference",
      "Group")
  )
})

test_that("permutation FDR is deterministic", {
  geneset <- data.frame(
    Pathway = c("A", "A", "B", "B", "C", "C"),
    Gene = paste0("G", 1:6),
    Dataset = "toy"
  )
  genes <- paste0("G", 1:6)
  background <- as.vector(outer(genes, genes, paste, sep = "_"))
  background <- background[
    sub("_.*$", "", background) != sub("^.*_", "", background)
  ]
  query <- c("G1_G3", "G2_G4", "G5_G1")
  enrichment <- InterMetaEnrich(
    query,
    background,
    InterMetaPathwayPair(geneset),
    keep_gene_pairs = "none"
  )
  first <- InterMetaPermFDR(
    enrichment, geneset, query, background,
    niter = 3, core = 1, seed = 11
  )
  second <- InterMetaPermFDR(
    enrichment, geneset, query, background,
    niter = 3, core = 1, seed = 11
  )
  expect_equal(first$Permutation.FDR, second$Permutation.FDR)
})
