test_that("bundled TCGA-PRAD data are aligned and documented", {
  data("InterMetaPRADInfo", package = "InterMeta", envir = environment())
  data("InterMetaPRADRNA", package = "InterMeta", envir = environment())
  data("InterMetaPRADEnrich", package = "InterMeta", envir = environment())
  data("InterMetaPRADScore", package = "InterMeta", envir = environment())
  data("InterMetaPRADProvenance", package = "InterMeta",
       envir = environment())

  expect_equal(dim(InterMetaPRADInfo), c(463L, 2L))
  expect_equal(dim(InterMetaPRADRNA), c(1260L, 463L))
  expect_identical(colnames(InterMetaPRADRNA), row.names(InterMetaPRADInfo))
  expect_equal(sum(InterMetaPRADInfo$event == 1), 99L)
  expect_equal(sum(InterMetaPRADInfo$event == 0), 364L)
  expect_true(all(is.finite(InterMetaPRADRNA)))
  expect_true(any(InterMetaPRADRNA < 0))

  expect_equal(nrow(InterMetaPRADEnrich), 2652L)
  expect_equal(sum(InterMetaPRADEnrich$Adjust.P.Value < 0.05), 617L)
  expect_equal(sum(InterMetaPRADEnrich$Permutation.FDR < 0.2), 104L)
  expect_equal(sum(InterMetaPRADEnrich$Permutation.FDR < 0.05), 50L)
  expect_equal(nrow(InterMetaPRADScore), 50L)

  expect_identical(InterMetaPRADProvenance$cohort, "TCGA-PRAD")
  expect_named(
    InterMetaPRADProvenance$source,
    c(
      "database", "cohort", "url", "github", "doi", "access_date",
      "database_version", "snapshot"
    )
  )
  expect_identical(InterMetaPRADProvenance$source$database, "PCaDB")
  expect_identical(InterMetaPRADProvenance$source$cohort, "TCGA-PRAD")
  expect_identical(
    InterMetaPRADProvenance$source$url,
    "http://bioinfo.jialab-ucr.org/PCaDB/"
  )
  expect_identical(
    InterMetaPRADProvenance$source$github,
    "https://github.com/rli012/PCaDB"
  )
  expect_identical(
    InterMetaPRADProvenance$source$doi,
    "10.1101/2021.06.29.449134"
  )
  expect_true(is.na(InterMetaPRADProvenance$source$access_date))
  expect_true(is.na(InterMetaPRADProvenance$source$database_version))
  expect_match(
    InterMetaPRADProvenance$source$snapshot,
    "were not recorded"
  )
  expect_identical(InterMetaPRADProvenance$time_unit, "days")
  expect_equal(InterMetaPRADProvenance$cox_gene_pairs, 370230L)
  expect_identical(InterMetaPRADProvenance$permutation_iterations, 100L)
  expect_match(
    InterMetaPRADProvenance$permutation_fdr$source,
    "fresh InterMeta 1[.]1[.]3 example09"
  )
  expect_match(
    InterMetaPRADProvenance$permutation_fdr$source,
    "100 TimiGP-legacy annotation permutations"
  )
  expect_identical(
    InterMetaPRADProvenance$permutation_fdr$iterations,
    100L
  )
  expect_true(InterMetaPRADProvenance$permutation_fdr$recomputed)
  expect_true(
    InterMetaPRADProvenance$permutation_fdr$
      historical_reference_used_for_validation_only
  )
  fdr_generation <- InterMetaPRADProvenance$permutation_fdr$generation
  expect_true(fdr_generation$precomputed_cox_input)
  expect_identical(fdr_generation$analysis_profile, "TimiGP-legacy")
  expect_identical(fdr_generation$niter, 100L)
  expect_identical(fdr_generation$seed, 1L)
  expect_lte(
    fdr_generation$historical_reference_max_abs_difference,
    1e-12
  )
  fdr_reproduction <-
    InterMetaPRADProvenance$permutation_fdr$fresh_reproduction
  expect_identical(
    fdr_reproduction$validation_id,
    paste0(
      "prad-complete-full-1.1.3/",
      "run_20260917_101835_fresh/",
      "InterMeta_TCGA_PRAD_complete_1_1_3"
    )
  )
  expect_identical(fdr_reproduction$package_version, "1.1.3")
  expect_false(fdr_reproduction$precomputed_statistics_used_for_analysis)
  expect_false(fdr_reproduction$reuse_checkpoints)
  expect_identical(fdr_reproduction$analysis_profile, "TimiGP-legacy")
  expect_identical(fdr_reproduction$pair_condition, "PV")
  expect_identical(fdr_reproduction$pair_cutoff, 0.05)
  expect_identical(
    fdr_reproduction$score_condition,
    "Adjust.P.Value"
  )
  expect_identical(fdr_reproduction$score_cutoff, 0.05)
  expect_identical(fdr_reproduction$niter, 100L)
  expect_identical(fdr_reproduction$permutation_seed, 1L)
  expect_identical(fdr_reproduction$cox_gene_pairs, 370230L)
  expect_identical(fdr_reproduction$selected_gene_pairs, 101323L)
  expect_identical(fdr_reproduction$enrichment_rows, 2652L)
  expect_identical(fdr_reproduction$score_rows, 50L)
  expect_identical(fdr_reproduction$parity_checks_passed, 12L)
  expect_identical(fdr_reproduction$resamples_requested, 3L)
  expect_identical(fdr_reproduction$resamples_successful, 3L)
  expect_identical(fdr_reproduction$sample_fraction, 0.8)
  expect_identical(
    fdr_reproduction$resample_seeds,
    20260701:20260703
  )
  expect_identical(fdr_reproduction$output_manifest_entries, 45L)
  expect_identical(
    fdr_reproduction$output_manifest_sha256,
    paste0(
      "c8eeab6e34859adcd982ac47a7c0ae5c32772fa5f06a36f9c",
      "c6909fafb9ce82b"
    )
  )
  expect_identical(
    fdr_reproduction$acceptance_candidate_archive_sha256,
    paste0(
      "9c2502a06645a0ca58e0879bd05eaad96f8700e10322f6494a",
      "8d8aaf5d661f3d"
    )
  )
  expect_match(
    fdr_reproduction$archive_binding,
    "did not independently record"
  )
  expect_lte(
    fdr_reproduction$historical_reference_max_abs_difference,
    1e-12
  )
  expect_lte(
    fdr_reproduction$
      prior_bundle_max_abs_difference_at_first_publication,
    1e-12
  )
  expect_named(
    InterMetaPRADProvenance$input_sha256,
    c(
      "conditions", "expression", "geneset", "prad_cox",
      "prad_enrichment_reference"
    )
  )
  expect_equal(length(InterMetaPRADProvenance$input_sha256), 5L)
  expect_true(all(nchar(InterMetaPRADProvenance$input_sha256) == 64L))
})

test_that("public PRAD generator recomputes rather than copies permutation FDR", {
  generator <- system.file(
    "scripts", "generate_prad_data.R", package = "InterMeta"
  )
  expect_true(nzchar(generator))
  script <- paste(readLines(generator, warn = FALSE), collapse = "\n")
  expect_match(script, "InterMetaPermFDR\\s*\\(", perl = TRUE)
  expect_false(grepl(
    paste0(
      "InterMetaPRADEnrich\\$Permutation[.]FDR\\s*<-\\s*",
      "prad_enrich_reference\\$Permutation[.]FDR"
    ),
    script,
    perl = TRUE
  ))
  expect_match(script, "permutation_cores\\s*<-\\s*1L", perl = TRUE)
  expect_match(
    script,
    "historical_reference_used_for_validation_only\\s*=\\s*TRUE",
    perl = TRUE
  )
})

test_that("public PRAD provenance updater is fail-closed", {
  updater <- system.file(
    "scripts", "update_prad_provenance_impl.R", package = "InterMeta"
  )
  expect_true(nzchar(updater))
  script <- paste(readLines(updater, warn = FALSE), collapse = "\n")
  expect_match(script, "Reuse checkpoints", fixed = TRUE)
  expect_match(script, "expected_checks <- c", fixed = TRUE)
  expect_match(script, "output_manifest_sha256.tsv", fixed = TRUE)
  expect_match(
    script,
    "prior_bundle_max_abs_difference > 1e-12",
    fixed = TRUE
  )
  expect_match(script, "current_fingerprint", fixed = TRUE)
  expect_match(script, "run_provenance_file", fixed = TRUE)
  expect_match(script, "resamples_successful", fixed = TRUE)
})

test_that("TCGA-PRAD expression passes InterMeta preprocessing", {
  data("InterMetaPRADInfo", package = "InterMeta", envir = environment())
  data("InterMetaPRADRNA", package = "InterMeta", envir = environment())

  pathwayset <- InterMetaLoadGeneset("all")
  info <- InterMetaCheckEvent(
    InterMetaPRADInfo,
    min_events = 5,
    min_censored = 5
  )
  rna <- InterMetaPreprocess(
    pathwayset = pathwayset,
    rna = InterMetaPRADRNA,
    cohort = row.names(info),
    log = FALSE,
    GMNorm = TRUE
  )

  # Explicit log=FALSE still retains the legacy half-sample positivity filter.
  expect_equal(dim(rna), c(1018L, 463L))
  expect_true(all(is.finite(rna)))
  expect_equal(max(abs(apply(rna, 1, stats::median))), 0)
})
