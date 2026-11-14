# Optional post-analysis PubMed and AI evidence layer.
#
# The default path is fully offline: it creates a privacy-preserving context,
# deterministic PubMed queries, and a report skeleton. Set
# INTERMETA_RUN_LIVE_PUBMED=true to perform the explicit NCBI request block.

library(InterMeta)

data("InterMetaExampleEnrich")
data("InterMetaExampleScore")

context <- InterMetaAIContext(
  enrichment = InterMetaExampleEnrich,
  score = InterMetaExampleScore,
  disease = "example cancer",
  endpoint = "overall survival",
  time_unit = "days",
  cohort = "InterMeta synthetic example",
  top_n = 5L,
  direction = "unfavorable"
)

queries <- InterMetaBuildPubMedQueries(context)
report_skeleton <- InterMetaAIReport(context, language = "en")

print(context$pathways)
print(queries[, c("pathway", "layer", "query")])
print(report_skeleton$warnings)

run_live_pubmed <- identical(
  tolower(Sys.getenv("INTERMETA_RUN_LIVE_PUBMED", "false")),
  "true"
)

if (run_live_pubmed) {
  ncbi_email <- Sys.getenv("NCBI_EMAIL", unset = "")
  if (!nzchar(ncbi_email)) {
    ncbi_email <- NULL
  }
  evidence <- InterMetaPubMedSearch(
    context,
    queries = queries,
    retmax = 20L,
    email = ncbi_email,
    cache_dir = file.path(tempdir(), "InterMeta_pubmed_cache")
  )

  evidence_report <- InterMetaAIReport(
    context,
    evidence = evidence,
    language = "en"
  )
  print(evidence_report$warnings)
}

# To add model-written text, securely configure DEEPSEEK_API_KEY in the R
# process environment. Never paste a key into this file or an R history.
# This block is commented out to prevent implicit billable model calls.
#
# provider <- InterMetaDeepSeekProvider(model = "deepseek-v4-pro")
#
# ai_report <- InterMetaAIReport(
#   context,
#   evidence = evidence,
#   provider = provider,
#   language = "en",
#   outdir = "InterMeta_AI_report"
# )
