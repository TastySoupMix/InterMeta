test_that("complete TCGA-PRAD example has a safe static contract", {
  example_file <- system.file(
    "examples",
    "example09_tcga_prad_complete_workflow.R",
    package = "InterMeta"
  )
  expect_true(nzchar(example_file))
  if (!nzchar(example_file)) {
    return(invisible(NULL))
  }

  parsed <- NULL
  expect_silent(
    parsed <- parse(example_file, keep.source = TRUE)
  )
  expect_type(parsed, "expression")

  lines <- readLines(example_file, warn = FALSE, encoding = "UTF-8")
  script <- paste(lines, collapse = "\n")

  # Extract only syntactic function-call tokens. Object names such as the
  # bundled data objects are intentionally not treated as function calls.
  parse_data <- getParseData(parsed)
  calls <- unique(parse_data$text[
    parse_data$token == "SYMBOL_FUNCTION_CALL"
  ])
  intermeta_calls <- calls[grepl("^InterMeta", calls)]
  missing_exports <- setdiff(
    intermeta_calls,
    getNamespaceExports("InterMeta")
  )
  expect_length(missing_exports, 0L)

  required_calls <- c(
    "InterMetaCheckEvent",
    "InterMetaPreprocess",
    "InterMetaGenePair",
    "InterMetaCox",
    "InterMetaBG",
    "InterMetaEnrich",
    "InterMetaPermFDR",
    "InterMetaPathwayScore",
    "InterMetaStability",
    "InterMetaAIContext",
    "InterMetaBuildPubMedQueries",
    "InterMetaAIWorkflow",
    "InterMetaAIResume",
    "InterMetaAIStatus"
  )
  expect_setequal(intersect(calls, required_calls), required_calls)

  data_pattern <- "data\\s*\\(\\s*['\"]([^'\"]+)['\"]"
  data_hits <- regmatches(
    script,
    gregexpr(data_pattern, script, perl = TRUE)
  )[[1L]]
  data_objects <- sub(data_pattern, "\\1", data_hits, perl = TRUE)
  public_prad_objects <- c(
    "InterMetaPRADInfo",
    "InterMetaPRADRNA",
    "InterMetaPRADEnrich",
    "InterMetaPRADScore",
    "InterMetaPRADProvenance"
  )
  expect_true(all(data_objects %in% public_prad_objects))
  expect_true(all(c(
    "InterMetaPRADInfo", "InterMetaPRADRNA", "InterMetaPRADProvenance"
  ) %in%
                    data_objects))

  # Operational environment overrides must not alter the fixed scientific
  # parameters, and a clean acceptance run can explicitly disable reuse/AI.
  expect_match(script, "INTERMETA_OUTPUT_DIR", fixed = TRUE)
  expect_match(script, "INTERMETA_RUN_FULL_SAMPLE", fixed = TRUE)
  expect_match(script, "INTERMETA_RUN_STABILITY", fixed = TRUE)
  expect_match(script, "INTERMETA_RUN_AI_WORKFLOW", fixed = TRUE)
  expect_match(script, "INTERMETA_REUSE_CHECKPOINTS", fixed = TRUE)

  # Resume regression: validate the aggregate and Cox checkpoints up front.
  # Either missing object triggers recovery, while each expensive substage is
  # reconstructed only behind its own missing-object guard.
  full_file_assignment <- regexpr(
    "(?m)^full_result_file\\s*<-\\s*file\\.path\\(",
    script,
    perl = TRUE
  )[[1L]]
  cox_file_assignment <- regexpr(
    "(?m)^cox_file\\s*<-\\s*file\\.path\\(",
    script,
    perl = TRUE
  )[[1L]]
  enrich_file_assignment <- regexpr(
    "(?m)^enrich_file\\s*<-\\s*file\\.path\\(",
    script,
    perl = TRUE
  )[[1L]]
  full_checkpoint_load <- regexpr(
    paste0(
      "(?m)^full_result\\s*<-\\s*load_checkpoint\\s*\\(\\s*",
      "full_result_file\\s*,"
    ),
    script,
    perl = TRUE
  )[[1L]]
  cox_checkpoint_load <- regexpr(
    "(?m)^cox_res\\s*<-\\s*load_checkpoint\\s*\\(\\s*cox_file\\s*,",
    script,
    perl = TRUE
  )[[1L]]
  need_assignment <- regexpr(
    paste0(
      "(?m)^need_full_computation\\s*<-\\s*is.null\\(full_result\\)",
      "\\s*\\|\\|\\s*is.null\\(cox_res\\)"
    ),
    script,
    perl = TRUE
  )[[1L]]
  recovery_gate <- regexpr(
    "if\\s*\\(need_full_computation\\s*&&\\s*!run_full_sample\\)",
    script,
    perl = TRUE
  )[[1L]]
  recovery_block <- regexpr(
    "if\\s*\\(need_full_computation\\)\\s*\\{",
    script,
    perl = TRUE
  )[[1L]]
  cox_rebuild_guard <- regexpr(
    "(?m)^  if\\s*\\(is.null\\(cox_res\\)\\)\\s*\\{",
    script,
    perl = TRUE
  )[[1L]]
  full_rebuild_guard <- regexpr(
    "(?m)^if\\s*\\(is.null\\(full_result\\)\\)\\s*\\{",
    script,
    perl = TRUE
  )[[1L]]
  enrich_checkpoint_load <- regexpr(
    paste0(
      "(?m)^  enrich\\s*<-\\s*load_checkpoint\\s*\\(\\s*",
      "enrich_file\\s*,"
    ),
    script,
    perl = TRUE
  )[[1L]]
  expect_gt(full_file_assignment, 0L)
  expect_gt(cox_file_assignment, 0L)
  expect_gt(enrich_file_assignment, 0L)
  expect_gt(full_checkpoint_load, full_file_assignment)
  expect_gt(cox_checkpoint_load, cox_file_assignment)
  expect_gt(need_assignment, full_checkpoint_load)
  expect_gt(need_assignment, cox_checkpoint_load)
  expect_gt(recovery_gate, need_assignment)
  expect_gt(recovery_block, recovery_gate)
  expect_gt(cox_rebuild_guard, recovery_block)
  expect_gt(full_rebuild_guard, cox_rebuild_guard)
  expect_gt(enrich_checkpoint_load, full_rebuild_guard)

  # The packaged reference result may be used for parity checks, but the AI
  # context must be bound to the results recomputed by this workflow.
  expect_match(
    script,
    "analysis_for_ai\\s*<-\\s*list\\s*\\([\\s\\S]*?enrich\\s*=\\s*enrich",
    perl = TRUE
  )
  expect_match(
    script,
    "analysis_for_ai\\s*<-\\s*list\\s*\\([\\s\\S]*?score\\s*=\\s*score",
    perl = TRUE
  )
  expect_match(
    script,
    paste0(
      "analysis_for_ai\\s*<-\\s*list\\s*\\([\\s\\S]*?",
      "stability\\s*=\\s*stability"
    ),
    perl = TRUE
  )
  for (parameter in c(
    "analysis_profile", "pair_condition", "pair_cutoff",
    "score_condition", "score_cutoff"
  )) {
    expect_match(
      script,
      paste0(
        "analysis_for_ai\\s*<-\\s*list\\s*\\([\\s\\S]{0,1000}?",
        parameter, "\\s*=\\s*", parameter
      ),
      perl = TRUE
    )
  }
  expect_match(
    script,
    paste0(
      "InterMetaAIContext\\s*\\(\\s*",
      "analysis\\s*=\\s*analysis_for_ai"
    ),
    perl = TRUE
  )

  expect_match(script, "permutation_niter\\s*<-\\s*100L?", perl = TRUE)
  expect_match(script, "top_n\\s*=\\s*5L?", perl = TRUE)
  expect_match(script, "direction\\s*=\\s*['\"]unfavorable['\"]")
  expect_match(
    script,
    "endpoint\\s*=\\s*['\"]biochemical recurrence \\(BCR\\)['\"]"
  )
  expect_match(script, "time_unit\\s*=\\s*['\"]days['\"]")
  expect_match(script, "cohort\\s*=\\s*['\"]TCGA-PRAD['\"]")
  expect_match(script, "RNGkind\\s*\\(", perl = TRUE)
  expect_match(script, "output_manifest_sha256[.]tsv", perl = TRUE)
  expect_match(script, "prad_data_provenance[.]txt", perl = TRUE)

  # Paid requests and retries must both default to off, with credentials read
  # only through an environment-variable name.
  expect_match(script, "allow_model_calls\\s*<-\\s*FALSE")
  expect_match(script, "retry_failed_model_calls\\s*<-\\s*FALSE")
  expect_match(
    script,
    "allow_model_calls\\s*=\\s*allow_model_calls"
  )
  expect_match(
    script,
    "retry_failed\\s*=\\s*retry_failed_model_calls"
  )
  expect_match(script, "Sys.getenv\\(ai_key_env")
  expect_false(grepl("sk-[[:alnum:]_.-]{10,}", script, perl = TRUE))

  # A model or prompt-affecting setting must lead to a distinct AI run
  # directory; authorization, retry, and credential controls must not affect
  # the scientific/configuration identity of that directory.
  expect_match(
    script,
    paste0(
      "ai_configuration\\s*<-\\s*list\\s*\\([\\s\\S]{0,1500}?",
      "model\\s*=\\s*ai_model"
    ),
    perl = TRUE
  )
  expect_match(
    script,
    paste0(
      "ai_configuration_hash\\s*<-\\s*digest::digest\\s*\\(\\s*",
      "ai_configuration"
    ),
    perl = TRUE
  )
  expect_match(
    script,
    paste0(
      "paste0\\s*\\(\\s*['\"]run-['\"]\\s*,\\s*context_id\\s*,",
      "\\s*['\"]-['\"]\\s*,\\s*substr\\s*\\(\\s*",
      "ai_configuration_hash"
    ),
    perl = TRUE
  )
  config_start <- regexpr(
    "ai_configuration\\s*<-\\s*list\\s*\\(", script, perl = TRUE
  )[[1L]]
  config_end <- regexpr(
    "ai_configuration_hash\\s*<-", script, perl = TRUE
  )[[1L]]
  expect_gt(config_start, 0L)
  expect_gt(config_end, config_start)
  if (config_start > 0L && config_end > config_start) {
    config_text <- substr(script, config_start, config_end - 1L)
    expect_false(grepl("allow_model_calls", config_text, fixed = TRUE))
    expect_false(grepl("retry_failed_model_calls", config_text, fixed = TRUE))
    expect_false(grepl("ai_key_env", config_text, fixed = TRUE))
  }

  # The human-readable summary is portable: it records the AI run beneath
  # the workflow root without serializing the server's absolute output path.
  expect_match(
    script,
    paste0(
      "AI run directory \\(relative to workflow root\\): [\"']",
      "[\\s\\S]{0,200}?file\\.path\\(\\s*[\"']03_ai_workflow[\"']",
      "\\s*,\\s*basename\\(ai_run_dir\\)"
    ),
    perl = TRUE
  )
  expect_false(grepl(
    "paste0\\(\\s*[\"']AI run directory: [\"']\\s*,\\s*ai_run_dir",
    script,
    perl = TRUE
  ))
  # Omitting this argument selects the workflow's run-parent-relative cache
  # specification; passing file.path(ai_root_dir, ...) would freeze an
  # absolute server path in workflow_config.rds.
  expect_false(grepl("pubmed_cache_dir\\s*=", script, perl = TRUE))

  # Construct the private-cohort markers at test time so the source package
  # itself never contains those identifiers.
  private_ascii <- intToUtf8(c(115L, 117L, 110L))
  private_non_ascii <- intToUtf8(c(0x771fL, 0x5b9eL, 0x4e16L, 0x754cL))
  expect_false(grepl(
    paste0("\\b", private_ascii, "\\b"),
    script,
    ignore.case = TRUE,
    perl = TRUE
  ))
  expect_false(grepl(private_non_ascii, script, fixed = TRUE))
})
