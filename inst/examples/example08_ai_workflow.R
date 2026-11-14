# Complete InterMeta 1.1.3 AI workflow using only bundled public TCGA-PRAD
# results. The workflow performs deterministic PubMed retrieval, bounded
# article screening, pathway-level synthesis, local validation, reporting,
# and fixed figure rendering.
#
# Configure DEEPSEEK_API_KEY in the server or RStudio environment before
# starting R. Never paste a credential into this script, the package source,
# an R history file, or a saved workflow directory.

library(InterMeta)

if (utils::packageVersion("InterMeta") < "1.1.3") {
  stop("InterMeta 1.1.3 or newer is required.")
}

data("InterMetaPRADEnrich")
data("InterMetaPRADScore")

if (!nzchar(Sys.getenv("DEEPSEEK_API_KEY", unset = ""))) {
  stop(
    "DEEPSEEK_API_KEY is not configured. Set it outside this script, ",
    "restart R, and run the example again."
  )
}

context <- InterMetaAIContext(
  enrichment = InterMetaPRADEnrich,
  score = InterMetaPRADScore,
  disease = "prostate adenocarcinoma",
  endpoint = "biochemical recurrence (BCR)",
  time_unit = "days",
  cohort = "TCGA-PRAD",
  top_n = 5L,
  direction = "unfavorable"
)

# The first call creates this directory. A later invocation resumes the same
# hash-bound run instead of repeating completed PubMed or model stages.
run_dir <- "InterMeta_AI_PRAD_1_1_3"

if (!file.exists(run_dir) && !dir.exists(run_dir)) {
  ai_result <- InterMetaAIWorkflow(
    context = context,
    run_dir = run_dir,
    allow_model_calls = TRUE,
    language = "zh",
    purpose = "biology_exploration",
    model = "deepseek-v4-pro",
    api_key_env = "DEEPSEEK_API_KEY",
    screen_batch_size = 20L,
    max_findings_per_pathway = 5L,
    render = TRUE,
    figure_text_mode = "compact_english"
  )
} else if (dir.exists(run_dir)) {
  ai_result <- InterMetaAIResume(
    run_dir = run_dir,
    allow_model_calls = TRUE,
    retry_failed = FALSE
  )
} else {
  stop("run_dir exists but is not a directory: ", run_dir)
}

# Status inspection is local and makes no PubMed or model request.
workflow_status <- InterMetaAIStatus(run_dir)
print(workflow_status$status)
print(workflow_status$stages)

if (identical(workflow_status$status, "complete")) {
  print(ai_result$report$exploration_summary)

  # The workflow has already written its PDF and PNG. These are their paths.
  figure_files <- lapply(ai_result$figure, function(x) x$output)
  print(figure_files)

  # Validate and inspect the fixed figure data without creating a new file.
  # This is entirely local and makes no external request.
  figure_check <- InterMetaAIExplorationFigure(
    ai_result$report,
    render = FALSE,
    text_mode = "compact_english"
  )
  print(figure_check$caption)
}

# Safe continuation after an interruption:
#
# ai_result <- InterMetaAIResume(
#   run_dir,
#   allow_model_calls = TRUE,
#   retry_failed = FALSE
# )
#
# If a stage is failed or delivery-uncertain, inspect InterMetaAIStatus() and
# the saved provider metadata first. A new, potentially billable attempt must
# be authorized deliberately with retry_failed = TRUE.

# Safe offline replay after a completed run. Choose a NEW outdir each time.
# Replay revalidates the saved evidence and synthesis and redraws the figure;
# it makes zero PubMed calls and zero model calls.
#
# replay_result <- InterMetaAIReplay(
#   run_dir,
#   outdir = "InterMeta_AI_PRAD_1_1_3_replay",
#   render = TRUE,
#   text_mode = "compact_english"
# )
