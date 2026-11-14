make_score <- function(pathway, difference, group, rank = 1) {
  data.frame(
    Pathway = pathway,
    Rank = rank,
    Difference = difference,
    Group = group,
    stringsAsFactors = FALSE
  )
}

test_that("stability defaults to exploratory PV", {
  expect_identical(formals(InterMetaStability)$pair_condition, "PV")
})

test_that("direction consistency is conditional on pathway selection", {
  scores <- list(
    make_score("A", 2, "F"),
    make_score("A", 1, "F"),
    make_score(character(), numeric(), character(), numeric())
  )
  result <- InterMetaSummarizeStability(
    scores = scores,
    pathways = c("A", "B"),
    top = 1
  )
  pathway_a <- result$summary[result$summary$Pathway == "A", ]
  pathway_b <- result$summary[result$summary$Pathway == "B", ]
  expect_equal(pathway_a$Selection.Frequency, 2 / 3)
  expect_equal(pathway_a$Direction.Consistency, 1)
  expect_equal(pathway_b$Selection.Frequency, 0)
  expect_true(is.na(pathway_b$Direction.Consistency))
})

test_that("stability top sets use exact magnitude order without changing legacy ranks", {
  legacy_score <- make_score(
    pathway = c("B_favorable", "A_unfavorable", "C_small", "D_nonfinite"),
    difference = c(10, -10, 1, Inf),
    group = c("F", "U", "F", "U"),
    rank = c(1, 4, 2, 3)
  )
  original_rank <- legacy_score$Rank
  result <- InterMetaSummarizeStability(
    scores = list(run_one = legacy_score, run_two = legacy_score),
    pathways = legacy_score$Pathway,
    top = 1
  )

  summary <- result$summary
  unfavorable <- summary[summary$Pathway == "A_unfavorable", ]
  favorable <- summary[summary$Pathway == "B_favorable", ]
  nonfinite <- summary[summary$Pathway == "D_nonfinite", ]

  expect_identical(result$scores[[1]]$Rank, original_rank)
  expect_identical(result$scores[[2]]$Rank, original_rank)
  expect_equal(sum(summary$Top1.Frequency), 1)
  expect_equal(unfavorable$Top1.Frequency, 1)
  expect_equal(favorable$Top1.Frequency, 0)
  expect_equal(favorable$Top1.Favorable.Frequency, 1)
  expect_equal(unfavorable$Top1.Unfavorable.Frequency, 1)
  expect_equal(nonfinite$Top1.Unfavorable.Frequency, 0)
  expect_equal(unfavorable$Mean.Magnitude.Rank, 1)
  expect_equal(favorable$Mean.Magnitude.Rank, 2)
  expect_true(is.na(nonfinite$Mean.Magnitude.Rank))
  expect_equal(unfavorable$SD.Magnitude.Rank, 0)
  expect_identical(summary$Pathway[1], "A_unfavorable")
  expect_identical(unfavorable$Stability.Class, "High")
  expect_identical(favorable$Stability.Class, "Low")
  expect_identical(result$top_set_definition$K, 1L)
  expect_identical(
    result$top_set_definition$ranking,
    "descending abs(Difference)"
  )
})

test_that("pairwise magnitude overlap is reported overall and by direction", {
  first <- make_score(
    c("A_unfavorable", "B_favorable"),
    c(-10, 8),
    c("U", "F"),
    c(2, 1)
  )
  second <- make_score(
    c("A_unfavorable", "C_favorable"),
    c(-9, 7),
    c("U", "F"),
    c(2, 1)
  )
  result <- InterMetaSummarizeStability(
    scores = list(first = first, second = second),
    top = 1
  )

  expect_equal(result$pairwise$Top1.Jaccard, 1)
  expect_equal(result$pairwise$Top1.Favorable.Jaccard, 0)
  expect_equal(result$pairwise$Top1.Unfavorable.Jaccard, 1)
  expect_equal(
    result$summary$Top1.Favorable.Frequency[
      result$summary$Pathway == "B_favorable"
    ],
    0.5
  )
})

test_that("zero-run, one-run, and empty-direction outputs keep stable types", {
  zero <- InterMetaSummarizeStability(
    scores = list(), pathways = c("A", "B"), top = 1
  )
  expected_pairwise <- c(
    "Run.A", "Run.B", "Common.Pathways",
    "Difference.Pearson.Selected", "Difference.Spearman.Selected",
    "Difference.Pearson.All", "Difference.Spearman.All",
    "Direction.Agreement", "Top1.Jaccard",
    "Top1.Favorable.Jaccard", "Top1.Unfavorable.Jaccard"
  )
  expect_identical(names(zero$pairwise), expected_pairwise)
  expect_identical(zero$successful_runs, 0L)
  expect_type(zero$summary$Selection.Count, "integer")
  expect_type(zero$summary$Selection.Frequency, "double")
  expect_type(zero$summary$Mean.Magnitude.Rank, "double")
  expect_type(zero$summary$Stability.Class, "character")
  expect_type(zero$pairwise$Run.A, "character")
  expect_type(zero$pairwise$Common.Pathways, "integer")
  expect_type(zero$pairwise$Top1.Favorable.Jaccard, "double")

  one <- InterMetaSummarizeStability(
    scores = list(make_score("A", 2, "F")), top = 1
  )
  expect_identical(names(one$pairwise), expected_pairwise)
  expect_equal(nrow(one$pairwise), 0L)
  expect_true(all(one$summary$Stability.Class == "Not assessed"))

  empty_score <- make_score(
    character(), numeric(), character(), numeric()
  )
  empty <- InterMetaSummarizeStability(
    scores = list(empty_score, empty_score), top = 1
  )
  expect_equal(nrow(empty$summary), 0L)
  expect_type(empty$summary$Stability.Class, "character")
  expect_equal(nrow(empty$pairwise), 1L)
  expect_true(is.na(empty$pairwise$Top1.Jaccard))
  expect_true(is.na(empty$pairwise$Top1.Favorable.Jaccard))
  expect_true(is.na(empty$pairwise$Top1.Unfavorable.Jaccard))

  favorable_only <- InterMetaSummarizeStability(
    scores = list(
      make_score("A", 2, "F"),
      make_score("A", 1, "F")
    ),
    top = 1
  )
  expect_type(
    favorable_only$summary$Top1.Unfavorable.Frequency,
    "double"
  )
  expect_equal(favorable_only$summary$Top1.Unfavorable.Frequency, 0)
  expect_type(
    favorable_only$pairwise$Top1.Unfavorable.Jaccard,
    "double"
  )
  expect_true(is.na(favorable_only$pairwise$Top1.Unfavorable.Jaccard))
})

test_that("magnitude Top-K is sign symmetric and row-order deterministic", {
  original <- make_score(
    c("B", "A", "C"), c(5, -5, 1), c("F", "U", "F"), c(1, 3, 2)
  )
  flipped <- original[c(3, 1, 2), , drop = FALSE]
  flipped$Difference <- -flipped$Difference
  flipped$Group <- ifelse(flipped$Group == "F", "U", "F")

  first <- InterMetaSummarizeStability(
    scores = list(original, original), top = 1
  )
  second <- InterMetaSummarizeStability(
    scores = list(flipped, flipped), top = 1
  )
  expect_identical(
    first$summary$Pathway[first$summary$Top1.Frequency == 1],
    second$summary$Pathway[second$summary$Top1.Frequency == 1]
  )
  expect_identical(
    first$summary$Top1.Favorable.Frequency,
    second$summary$Top1.Unfavorable.Frequency
  )
  expect_identical(first$top_set_definition$classification_top, 1L)
})

test_that("invalid public score labels fail instead of becoming neutral", {
  bad_group <- make_score("A", 1, "typo", 1)
  expect_error(
    InterMetaSummarizeStability(scores = list(bad_group), top = 1),
    "Score Group values"
  )
  bad_pathway <- make_score("", 1, "F", 1)
  expect_error(
    InterMetaSummarizeStability(scores = list(bad_pathway), top = 1),
    "Pathway values"
  )
})

test_that("active manifests are replaced through complete staging files", {
  manifest_file <- tempfile("intermeta-active-manifest-", fileext = ".txt")
  previous_file <- paste0(manifest_file, ".previous")
  temporary_pattern <- file.path(
    dirname(manifest_file),
    paste0(".", basename(manifest_file), "-*.tmp")
  )

  InterMeta:::.InterMetaWriteTableAtomic(
    data.frame(Status = "pending", stringsAsFactors = FALSE),
    manifest_file
  )
  InterMeta:::.InterMetaWriteTableAtomic(
    data.frame(Status = "success", stringsAsFactors = FALSE),
    manifest_file
  )

  expect_identical(
    utils::read.delim(manifest_file, stringsAsFactors = FALSE)$Status,
    "success"
  )
  expect_length(Sys.glob(temporary_pattern), 0L)
  expect_false(file.exists(previous_file))
  unlink(c(manifest_file, previous_file), force = TRUE)
})

test_that("three-run stability workflow handles neutral scores and checkpoints", {
  set.seed(42)
  samples <- paste0("S", seq_len(40))
  genes <- paste0("G", seq_len(8))
  rna <- matrix(
    stats::rnorm(length(genes) * length(samples)),
    nrow = length(genes),
    dimnames = list(genes, samples)
  )
  info <- data.frame(
    event = rep(c(0, 1), each = 20),
    time = seq_len(40) + stats::runif(40),
    row.names = samples
  )
  geneset <- data.frame(
    Pathway = rep(paste0("P", 1:4), each = 2),
    Gene = genes,
    Dataset = "toy"
  )
  outdir <- tempfile("intermeta-stability-")

  neutral_warnings <- character()
  stability <- withCallingHandlers(suppressMessages(InterMetaStability(
    rna = rna,
    info = info,
    pathwayset = geneset,
    n_resample = 3,
    sample_fraction = 0.8,
    seeds = 1:3,
    pair_condition = "QV",
    pair_cutoff = 1.01,
    score_condition = "Adjust.P.Value",
    score_cutoff = 1.01,
    log = FALSE,
    parallel = FALSE,
    outdir = outdir,
    min_events = 3,
    min_censored = 3,
    keep_details = FALSE
  )), warning = function(w) {
    neutral_warnings <<- c(neutral_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(neutral_warnings, 3L)
  expect_true(all(grepl("Only neutral pathway scores", neutral_warnings)))
  expect_true(all(stability$run_status$Status == "success"))
  expect_identical(stability$parameters$pair_condition, "QV")
  expect_equal(nrow(stability$summary), 4L)
  expect_true(all(stability$summary$Selection.Count == 0L))
  expect_true(nrow(stability$edge_summary) > 0)
  expect_equal(nrow(stability$score_long), 0L)
  expect_true(all(c(
    "Top5.Frequency", "Top5.Favorable.Frequency",
    "Top5.Unfavorable.Frequency", "Mean.Magnitude.Rank",
    "SD.Magnitude.Rank"
  ) %in% names(stability$summary)))
  expect_true(all(c(
    "Top5.Jaccard", "Top5.Favorable.Jaccard",
    "Top5.Unfavorable.Jaccard"
  ) %in% names(stability$pairwise)))
  expect_identical(stability$top_set_definition$K, c(5L, 10L, 20L))
  expect_true(file.exists(file.path(outdir, "run_status.txt")))
  expect_true(file.exists(file.path(outdir, "active_run_manifest.txt")))
  expect_true(file.exists(file.path(outdir, "stability_scores_long.txt")))
  expect_equal(nrow(stability$active_run_manifest), 3L)
  expect_equal(
    stability$active_run_manifest$Run.ID,
    sprintf("resample_%03d", 1:3)
  )
  expect_identical(
    stability$parameters$permutation_seed_scheme,
    "TimiGP-legacy-1:niter"
  )
  expect_true(all(!is.na(
    stability$active_run_manifest$Permutation.Seed
  )))

  failed_cache_file <- file.path(
    outdir, "resample_001", "run_result.rds"
  )
  failed_cache <- readRDS(failed_cache_file)
  relative_result_file <- file.path("resample_001", "run_result.rds")
  expect_identical(failed_cache$result_file, relative_result_file)
  expect_identical(stability$runs[[1]]$result_file, relative_result_file)
  failed_cache$status <- "error"
  failed_cache$error_message <- "simulated interrupted run"
  saveRDS(failed_cache, failed_cache_file)
  failed_error_file <- file.path(
    outdir, "resample_001", "error.txt"
  )
  writeLines(failed_cache$error_message, failed_error_file)
  recovered <- suppressWarnings(suppressMessages(InterMetaStability(
    rna = rna,
    info = info,
    pathwayset = geneset,
    n_resample = 3,
    sample_fraction = 0.8,
    seeds = 1:3,
    pair_condition = "QV",
    pair_cutoff = 1.01,
    score_condition = "Adjust.P.Value",
    score_cutoff = 1.01,
    log = FALSE,
    parallel = FALSE,
    outdir = outdir,
    min_events = 3,
    min_censored = 3,
    keep_details = FALSE
  )))
  expect_true(all(recovered$run_status$Status == "success"))
  expect_identical(readRDS(failed_cache_file)$status, "success")
  expect_false(file.exists(failed_error_file))

  first_key <- readRDS(
    file.path(outdir, "resample_001", "run_result.rds")
  )$cache_key
  portable_cache <- readRDS(failed_cache_file)
  portable_cache$result_file <- normalizePath(
    failed_cache_file, winslash = "/", mustWork = TRUE
  )
  saveRDS(portable_cache, failed_cache_file)
  writeLines("stale error after a completed checkpoint", failed_error_file)
  extended <- suppressWarnings(suppressMessages(InterMetaStability(
    rna = rna,
    info = info,
    pathwayset = geneset,
    n_resample = 4,
    sample_fraction = 0.8,
    seeds = 1:4,
    pair_condition = "QV",
    pair_cutoff = 1.01,
    score_condition = "Adjust.P.Value",
    score_cutoff = 1.01,
    log = FALSE,
    parallel = FALSE,
    outdir = outdir,
    min_events = 3,
    min_censored = 3,
    keep_details = FALSE
  )))
  extended_key <- readRDS(
    file.path(outdir, "resample_001", "run_result.rds")
  )$cache_key
  expect_identical(first_key, extended_key)
  expect_true(all(extended$run_status$Status == "success"))
  expect_false(file.exists(failed_error_file))
  expect_true(file.exists(file.path(outdir, "resample_004", "score.txt")))
  expect_identical(
    readRDS(failed_cache_file)$result_file,
    relative_result_file
  )

  changed <- suppressWarnings(suppressMessages(InterMetaStability(
    rna = rna,
    info = info,
    pathwayset = geneset,
    n_resample = 3,
    sample_fraction = 0.8,
    seeds = 1:3,
    pair_condition = "QV",
    pair_cutoff = 1,
    score_condition = "Adjust.P.Value",
    score_cutoff = 1.01,
    log = FALSE,
    parallel = FALSE,
    outdir = outdir,
    min_events = 3,
    min_censored = 3,
    keep_details = FALSE
  )))
  second_key <- readRDS(
    file.path(outdir, "resample_001", "run_result.rds")
  )$cache_key
  expect_false(identical(first_key, second_key))
  expect_true(all(changed$run_status$Status == "success"))
})

test_that("stability cache key includes the supplied pathway-pair object", {
  samples <- paste0("S", 1:6)
  genes <- paste0("G", 1:4)
  rna <- matrix(
    seq_len(length(samples) * length(genes)),
    nrow = length(genes),
    dimnames = list(genes, samples)
  )
  info <- data.frame(
    event = rep(c(0, 1), each = 3),
    time = seq_len(6),
    row.names = samples
  )
  geneset <- data.frame(
    Pathway = rep(c("P1", "P2"), each = 2),
    Gene = genes,
    Dataset = "toy"
  )
  pair_one <- InterMetaPathwayPair(geneset)
  pair_two <- pair_one
  pair_two$pathway_genes[[1]] <- rev(pair_two$pathway_genes[[1]])

  key_one <- InterMeta:::.InterMetaCacheKey(
    rna, info, geneset, pair_one, list(pair_cutoff = 0.05)
  )
  key_two <- InterMeta:::.InterMetaCacheKey(
    rna, info, geneset, pair_two, list(pair_cutoff = 0.05)
  )
  expect_false(identical(key_one, key_two))
})

test_that("stability validates scalar controls and derives safe seeds", {
  samples <- paste0("S", 1:10)
  genes <- paste0("G", 1:4)
  rna <- matrix(
    stats::rnorm(length(samples) * length(genes)),
    nrow = length(genes),
    dimnames = list(genes, samples)
  )
  info <- data.frame(
    event = rep(c(0, 1), each = 5),
    time = seq_len(10),
    row.names = samples
  )
  geneset <- data.frame(
    Pathway = rep(c("P1", "P2"), each = 2),
    Gene = genes,
    Dataset = "toy"
  )

  expect_error(
    InterMetaStability(
      rna, info, geneset, n_resample = 1.5,
      pair_condition = "QV", log = FALSE, parallel = FALSE
    ),
    "n_resample must be one positive integer"
  )
  expect_error(
    InterMetaStability(
      rna, info, geneset, sample_fraction = c(0.7, 0.8),
      pair_condition = "QV", log = FALSE, parallel = FALSE
    ),
    "sample_fraction must be one finite number"
  )
  expect_error(
    InterMetaStability(
      rna, info, geneset, seeds = c(1, 2.5, 3),
      pair_condition = "QV", log = FALSE, parallel = FALSE,
      min_events = 2, min_censored = 2
    ),
    "All used seeds must be finite integers"
  )

  derived <- vapply(
    c(-.Machine$integer.max, 0, .Machine$integer.max),
    function(seed) InterMeta:::.InterMetaPermutationSeed(seed, 1L),
    integer(1)
  )
  expect_true(all(!is.na(derived)))
  expect_true(all(derived >= 0L))
})

test_that("run-status export keeps one physical row per failed run", {
  runs <- list(list(
    run = 1L,
    run_id = "resample_001",
    seed = 20260701L,
    permutation_seed = 1L,
    status = "error",
    selected_samples = c("S1", "S2"),
    event_count = 1L,
    censored_count = 1L,
    query_pair_count = 0L,
    score = data.frame(),
    elapsed_seconds = 0.1,
    error_message = "first\tfield\r\nsecond line"
  ))
  status <- InterMeta:::.InterMetaRunStatus(runs)

  expect_identical(status$Error, "first field second line")
  expect_false(grepl("[\t\r\n]", status$Error))

  output <- tempfile(fileext = ".txt")
  on.exit(unlink(output), add = TRUE)
  utils::write.table(
    status, output, quote = FALSE, sep = "\t", row.names = FALSE
  )
  expect_length(readLines(output, warn = FALSE), 2L)
})
