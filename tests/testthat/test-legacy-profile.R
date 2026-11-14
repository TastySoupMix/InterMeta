test_that("TimiGP legacy preprocessing safely reproduces the old filters", {
  samples <- paste0("S", 1:6)
  geneset <- data.frame(
    Pathway = rep(c("A", "B"), each = 2),
    Gene = paste0("G", 1:4),
    Dataset = "toy"
  )
  rna <- rbind(
    G1 = c(-0.5, 0.1, 1, 2, 3, 4),
    G2 = c(-2, 1, 2, 3, 4, 5),
    G3 = c(-0.5, -0.2, 0, 1, 0, 0),
    G4 = c(-0.9, 0.2, 0.3, 0.4, 0.5, 0.6)
  )
  colnames(rna) <- samples

  current <- suppressMessages(InterMetaPreprocess(
    pathwayset = geneset,
    rna = rna
  ))
  legacy <- suppressMessages(InterMetaPreprocess(
    pathwayset = geneset,
    rna = rna,
    analysis_profile = "TimiGP-legacy"
  ))

  expect_equal(current, legacy)
  expect_equal(rownames(legacy), c("G1", "G4"))
  expect_true(all(is.finite(legacy)))
  expect_equal(as.numeric(apply(legacy, 1, median)), c(0, 0))
  expect_error(
    InterMetaPreprocess(
      pathwayset = geneset,
      rna = rna,
      analysis_profile = "unknown"
    ),
    "Only analysis_profile"
  )
})

test_that("TimiGP legacy preprocessing recovers 861 PRAD genes", {
  data("InterMetaPRADRNA", package = "InterMeta", envir = environment())
  legacy <- suppressMessages(InterMetaPreprocess(
    pathwayset = InterMetaLoadGeneset("all"),
    rna = InterMetaPRADRNA,
    analysis_profile = "TimiGP-legacy"
  ))

  expect_equal(dim(legacy), c(861L, 463L))
  expect_true(all(is.finite(legacy)))
  expect_equal(max(abs(apply(legacy, 1, median))), 0)
})

test_that("TimiGP legacy gene-pair threshold uses round", {
  samples <- paste0("S", 1:11)
  rna <- rbind(
    G1 = rep(0, 11),
    G2 = c(1, rep(-1, 10))
  )
  colnames(rna) <- samples

  current <- suppressMessages(InterMetaGenePair(rna))
  legacy <- suppressMessages(InterMetaGenePair(
    rna,
    analysis_profile = "TimiGP-legacy"
  ))

  expect_equal(current, legacy)
  expect_equal(nrow(legacy), 1L)
  expect_identical(rownames(legacy), "G2_G1")
})

test_that("TimiGP legacy logical Cox uses factor encoding", {
  samples <- paste0("S", 1:20)
  tag <- rep(c(FALSE, TRUE), each = 10)
  mps <- matrix(
    tag,
    nrow = 1,
    dimnames = list("G2_G1", samples)
  )
  info <- data.frame(
    event = rep(c(1, 0, 1, 1, 0), 4),
    time = c(2, 5, 7, 11, 14, 1, 4, 8, 10, 15,
             3, 6, 9, 13, 18, 2.5, 4.5, 9.5, 12, 20),
    row.names = samples
  )
  model_tag <- as.factor(tag)
  expected_fit <- survival::coxph(
    survival::Surv(info$time, info$event) ~ model_tag,
    data = data.frame(info, model_tag = model_tag)
  )
  expected_summary <- summary(expected_fit)
  expected_hr <- unname(expected_summary$conf.int[1, "exp(coef)"])
  expected_pv <- unname(
    expected_summary$coefficients[1, "Pr(>|z|)"]
  )
  expected_name <- "G2_G1"
  if (expected_hr > 1){
    expected_hr <- 1 / expected_hr
    expected_name <- "G1_G2"
  }

  result <- suppressMessages(InterMetaCox(
    mps,
    info,
    keep_mps = FALSE,
    analysis_profile = "TimiGP-legacy"
  ))$cox_res

  expect_identical(rownames(result), expected_name)
  expect_equal(result$HR, expected_hr, tolerance = 0)
  expect_equal(result$PV, expected_pv, tolerance = 0)
})

test_that("TimiGP legacy score removes neutral pathways and ranks signed", {
  enrichment <- data.frame(
    Index = 1:3,
    Rank = 1:3,
    Pathway.Interaction = c("A__B", "A__C", "C__B"),
    Favorable.Pathway = c("A", "A", "C"),
    Unfavorable.Pathway = c("B", "C", "B"),
    No.Total.MGP = 3L,
    No.Shared.MGP = 1L,
    Enrichment.Ratio = 1,
    P.Value = 0.01,
    Adjust.P.Value = 0.01,
    Shared.MGP = "",
    Total.MGP = NA_character_
  )

  current <- InterMetaPathwayScore(enrichment)
  legacy <- InterMetaPathwayScore(
    enrichment,
    analysis_profile = "TimiGP-legacy"
  )

  expect_equal(current, legacy)
  expect_equal(legacy$Pathway, c("A", "B"))
  expect_equal(legacy$Difference, c(200 / 3, -200 / 3))
  expect_equal(legacy$Rank, c(1, 2))
  expect_equal(legacy$Signed.Rank, legacy$Rank)
  expect_equal(legacy$Group, c("F", "U"))

  neutral <- enrichment[1:2, , drop = FALSE]
  neutral$Pathway.Interaction <- c("A__B", "B__A")
  neutral$Favorable.Pathway <- c("A", "B")
  neutral$Unfavorable.Pathway <- c("B", "A")
  expect_warning(
    empty <- InterMetaPathwayScore(
      neutral,
      analysis_profile = "TimiGP-legacy"
    ),
    "neutral"
  )
  expect_equal(nrow(empty), 0L)
})

test_that("TimiGP legacy permutation FDR uses seeds 1:n and legacy ranks", {
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
    niter = 3, core = 1, seed = 11,
    analysis_profile = "TimiGP-legacy"
  )
  second <- InterMetaPermFDR(
    enrichment, geneset, query, background,
    niter = 3, core = 1, seed = 999,
    analysis_profile = "TimiGP-legacy"
  )

  normalized <- InterMeta:::.InterMetaNormalizeGeneset(geneset)
  context <- InterMeta:::.InterMetaPreparePairContext(query, background)
  permutation_p <- unlist(lapply(1:3, function(iteration){
    set.seed(iteration)
    permuted <- normalized
    permuted$Gene <- sample(permuted$Gene, replace = FALSE)
    result <- InterMeta:::.InterMetaEnrichCompact(
      pathway.pair = InterMetaPathwayPair(permuted),
      pair_context = context,
      p.adj = "none",
      keep_gene_pairs = "none"
    )
    result$P.Value
  }), use.names = FALSE)
  expected <- vapply(seq_len(nrow(enrichment)), function(i){
    min(sum(permutation_p <= enrichment$P.Value[i]) /
          enrichment$Rank[i] / 3, 1)
  }, numeric(1))

  expect_equal(first$Permutation.FDR, expected)
  expect_equal(second$Permutation.FDR, expected)
  expect_error(
    InterMetaPermFDR(
      enrichment[, setdiff(names(enrichment), "Rank"), drop = FALSE],
      geneset, query, background,
      niter = 1, core = 1,
      analysis_profile = "TimiGP-legacy"
    ),
    "Rank"
  )
})

test_that("legacy permutation shuffles raw duplicate marker slots", {
  raw_geneset <- data.frame(
    Pathway = c("A", "A", "A", "B", "B", "C", "C"),
    Gene = c("G1", "G1", "G2", "G3", "G4", "G5", "G6"),
    Dataset = "toy"
  )
  geneset_file <- tempfile(fileext = ".txt")
  on.exit(unlink(geneset_file), add = TRUE)
  utils::write.table(
    raw_geneset,
    geneset_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  visible_geneset <- InterMetaLoadGeneset(file = geneset_file)
  expect_equal(nrow(visible_geneset), 6L)
  expect_equal(
    nrow(InterMeta:::.InterMetaLegacyPermutationGeneset(visible_geneset)),
    7L
  )

  genes <- paste0("G", 1:6)
  background <- as.vector(outer(genes, genes, paste, sep = "_"))
  background <- background[
    sub("_.*$", "", background) != sub("^.*_", "", background)
  ]
  query <- c("G2_G1", "G3_G1", "G4_G1", "G5_G1")
  enrichment <- InterMetaEnrich(
    query,
    background,
    InterMetaPathwayPair(visible_geneset),
    keep_gene_pairs = "none"
  )
  legacy <- InterMetaPermFDR(
    enrichment,
    visible_geneset,
    query,
    background,
    niter = 3,
    core = 1,
    seed = 999,
    analysis_profile = "TimiGP-legacy"
  )

  context <- InterMeta:::.InterMetaPreparePairContext(query, background)
  permutation_p <- unlist(lapply(1:3, function(iteration){
    set.seed(iteration)
    permuted <- raw_geneset
    permuted$Gene <- sample(permuted$Gene, replace = FALSE)
    result <- InterMeta:::.InterMetaEnrichCompact(
      pathway.pair = InterMetaPathwayPair(permuted),
      pair_context = context,
      p.adj = "none",
      keep_gene_pairs = "none"
    )
    result$P.Value
  }), use.names = FALSE)
  expected <- vapply(seq_len(nrow(enrichment)), function(i){
    min(sum(permutation_p <= enrichment$P.Value[i]) /
          enrichment$Rank[i] / 3, 1)
  }, numeric(1))
  expect_equal(legacy$Permutation.FDR, expected)

  unique_geneset <- unique(raw_geneset)
  legacy_unique <- InterMetaPermFDR(
    enrichment,
    unique_geneset,
    query,
    background,
    niter = 3,
    core = 1,
    analysis_profile = "TimiGP-legacy"
  )
  expect_false(isTRUE(all.equal(
    legacy$Permutation.FDR,
    legacy_unique$Permutation.FDR
  )))

  subset_geneset <- visible_geneset[
    visible_geneset$Pathway != "C",
    ,
    drop = FALSE
  ]
  expect_equal(
    nrow(InterMeta:::.InterMetaLegacyPermutationGeneset(subset_geneset)),
    sum(raw_geneset$Pathway != "C")
  )
  stale_geneset <- visible_geneset
  stale_geneset$Gene[1] <- "G7"
  expect_equal(
    nrow(InterMeta:::.InterMetaLegacyPermutationGeneset(stale_geneset)),
    nrow(stale_geneset)
  )
})

test_that("TimiGP legacy stability defaults reproduce old sampling controls", {
  set.seed(2026)
  samples <- paste0("S", 1:20)
  genes <- paste0("G", 1:4)
  rna <- matrix(
    stats::runif(length(samples) * length(genes), 0.2, 10),
    nrow = length(genes),
    dimnames = list(genes, samples)
  )
  info <- data.frame(
    event = rep(c(0, 1), each = 10),
    time = seq_len(20) + stats::runif(20),
    row.names = samples
  )
  geneset <- data.frame(
    Pathway = rep(c("P1", "P2"), each = 2),
    Gene = genes,
    Dataset = "toy"
  )

  stability <- suppressWarnings(suppressMessages(InterMetaStability(
    rna = rna,
    info = info,
    pathwayset = geneset,
    n_resample = 3,
    sample_fraction = 0.8,
    pair_cutoff = 1.01,
    score_cutoff = 1.01,
    parallel = FALSE,
    min_events = 3,
    min_censored = 3,
    continue_on_error = FALSE
  )))

  expect_identical(stability$parameters$analysis_profile, "TimiGP-legacy")
  expect_identical(stability$parameters$pair_condition, "PV")
  expect_identical(stability$parameters$seeds, 20260701:20260703)
  expect_true(stability$parameters$log)
  expect_false(stability$parameters$stratify_event)
  expect_identical(
    stability$parameters$permutation_seed_scheme,
    "TimiGP-legacy-1:niter"
  )
  expect_identical(
    stability$active_run_manifest$Permutation.Seed,
    rep(1L, 3)
  )
  expected_samples <- lapply(20260701:20260703, function(seed){
    set.seed(seed)
    sort(sample(samples, 16, replace = FALSE))
  })
  expect_true(all(vapply(seq_len(3), function(i){
    identical(stability$runs[[i]]$selected_samples, expected_samples[[i]])
  }, logical(1))))

  current_key <- InterMeta:::.InterMetaCacheKey(
    rna, info, geneset, stability$pathway_pair,
    list(analysis_profile = "InterMeta")
  )
  legacy_key <- InterMeta:::.InterMetaCacheKey(
    rna, info, geneset, stability$pathway_pair,
    list(analysis_profile = "TimiGP-legacy")
  )
  expect_false(identical(current_key, legacy_key))
})
