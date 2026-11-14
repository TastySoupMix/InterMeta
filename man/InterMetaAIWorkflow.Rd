\name{InterMetaAIWorkflow}
\alias{InterMetaAIWorkflow}
\alias{InterMetaAIResume}
\alias{InterMetaAIStatus}
\alias{InterMetaAIReplay}
\title{Run, Resume, Inspect, or Replay the InterMeta AI Workflow}
\usage{
InterMetaAIWorkflow(
  context,
  run_dir,
  evidence = NULL,
  queries = NULL,
  provider = NULL,
  allow_model_calls = FALSE,
  language = c("en", "zh"),
  purpose = c("biology_exploration", "prognostic_validation"),
  disease_mesh = NULL,
  pathway_synonyms = NULL,
  pathway_components = NULL,
  include_components = TRUE,
  retmax = 100L,
  page_size = 100L,
  sort = c("relevance", "pub_date"),
  pubmed_cache_dir = file.path(dirname(run_dir), "pubmed_cache"),
  ncbi_email_env = "NCBI_EMAIL",
  ncbi_api_key_env = "NCBI_API_KEY",
  pubmed_tool = "InterMeta",
  screen_batch_size = 20L,
  max_findings_per_pathway = 5L,
  model = "deepseek-v4-pro",
  api_key_env = "DEEPSEEK_API_KEY",
  model_max_tokens = 32768L,
  model_timeout = 900,
  thinking = c("enabled", "disabled"),
  reasoning_effort = c("low", "high", "max"),
  continue_on_error = TRUE,
  render = TRUE,
  figure_text_mode = c("compact_english", "full"),
  font_family = NULL
)

InterMetaAIResume(
  run_dir,
  provider = NULL,
  allow_model_calls = FALSE,
  retry_failed = FALSE,
  pubmed_cache_dir = NULL
)

InterMetaAIStatus(run_dir)

InterMetaAIReplay(
  run_dir,
  outdir = file.path(run_dir, "replay"),
  render = TRUE,
  text_mode = c("compact_english", "full"),
  font_family = NULL
)
}
\arguments{
  \item{context}{An \code{InterMetaAIContext} containing the aggregate Top
  pathways and analysis metadata. Patient-level input is not retained.}
  \item{run_dir}{A new workflow directory for
  \code{InterMetaAIWorkflow()}, or an existing workflow directory for the
  other functions. A new workflow never overwrites an existing path.}
  \item{evidence}{Optional frozen \code{InterMetaPubMedEvidence}. Supplying it
  avoids live PubMed retrieval, but its query plan must match.}
  \item{queries}{Optional query table from
  \code{InterMetaBuildPubMedQueries()}.}
  \item{provider}{Optional provider adapter. If omitted, the workflow creates
  the configured DeepSeek provider when a model stage is authorized. A
  provider function is never serialized into workflow configuration.}
  \item{allow_model_calls}{Logical safety gate. New model requests are disabled
  by default and require an explicit \code{TRUE}.}
  \item{language}{Language for screening annotations and synthesis text.}
  \item{purpose}{Workflow goal. The package-native two-stage workflow currently supports only
  \code{"biology_exploration"}; use \code{InterMetaAIInterpret()} directly
  for legacy endpoint-specific prognostic validation.}
  \item{disease_mesh, pathway_synonyms, pathway_components,
  include_components}{Query-building controls passed to
  \code{InterMetaBuildPubMedQueries()}. Component terms expand retrieval and
  are not treated as whole-pathway synonyms.}
  \item{retmax}{Positive integer maximum PubMed identifiers retrieved per
  query.}
  \item{page_size}{Positive integer ESearch page size.}
  \item{sort}{PubMed ordering, either \code{"relevance"} or
  \code{"pub_date"}.}
  \item{pubmed_cache_dir}{For a new workflow, the directory for validated
  query-level PubMed caches. When omitted, fixed configuration stores the
  relative path \code{"pubmed_cache"} anchored to the parent of
  \code{run_dir}. For \code{InterMetaAIResume()}, an optional operational
  cache override used only if the PubMed stage must run. This override is not
  persisted and does not change the fixed configuration or its hash.}
  \item{ncbi_email_env}{Name of an environment variable containing an
  optional NCBI contact email.}
  \item{ncbi_api_key_env}{Name of an environment variable containing an
  optional NCBI API key.}
  \item{pubmed_tool}{Tool label supplied to NCBI E-utilities.}
  \item{screen_batch_size}{Maximum candidates in one pathway-specific
  screening request.}
  \item{max_findings_per_pathway}{Integer from 1 through 5 giving the maximum
  number of findings selected in each pathway synthesis.}
  \item{model}{DeepSeek model identifier saved as non-secret configuration.}
  \item{api_key_env}{Name of the environment variable containing the DeepSeek
  credential. Pass a variable name, never an API-key value.}
  \item{model_max_tokens}{Positive integer output-token budget for each model
  request.}
  \item{model_timeout}{Positive model request timeout in seconds.}
  \item{thinking}{DeepSeek thinking-mode setting.}
  \item{reasoning_effort}{Requested DeepSeek reasoning effort.}
  \item{continue_on_error}{Logical; allow later batches or pathways to run
  after an isolated failure. An incomplete stage is not published as a
  complete report.}
  \item{render}{Logical; create the fixed local figure after report validation,
  or recreate it during replay.}
  \item{figure_text_mode}{Figure text mode used by the original workflow.}
  \item{font_family}{Optional installed font-family name for figure rendering.}
  \item{retry_failed}{Logical; together with
  \code{allow_model_calls = TRUE}, explicitly permits a new request after a
  failed or delivery-uncertain model attempt.}
  \item{outdir}{New directory to receive an offline replay. It is never
  overwritten.}
  \item{text_mode}{Figure text mode used by the replay.}
}
\value{
\code{InterMetaAIWorkflow()} and \code{InterMetaAIResume()} return an
\code{InterMetaAIWorkflowResult} with the run directory, overall status,
stage table, and any report and figure completed in that call. A result can be
\code{"paused"}, \code{"failed"}, or \code{"complete"}.

\code{InterMetaAIStatus()} returns an \code{InterMetaAIWorkflowStatus} with
per-stage status and hash-integrity checks, without network or model access.

\code{InterMetaAIReplay()} invisibly returns an
\code{InterMetaAIWorkflowReplay}. It revalidates saved inputs, publishes a new
report/figure directory, and makes zero PubMed and zero model calls.
}
\description{
Orchestrate the recoverable InterMeta PubMed, article-screening,
pathway-synthesis, report, and fixed-figure pipeline with hash-bound stage
artifacts and explicit authorization for model requests.
}
\details{
The eight stages are context, queries, PubMed retrieval, article screening,
verified evidence, pathway synthesis, report, and figure. Query construction
and PubMed retrieval are deterministic R stages. Screening and synthesis are
separate model stages. Screening admits or quarantines individual
pathway--PMID candidates; synthesis can select only verified candidate IDs.
R maps the frozen quotations back to their source records and validates the
final interpretation before publication. A structurally complete screening
with individually quarantined rows is recorded as
\code{"complete_with_quarantine"}; those rows remain in the audit and never
enter verified evidence or synthesis.

If one or more successful PubMed queries reach their configured retrieval
limit, the PubMed stage is \code{"complete_with_truncation"}. Downstream
stages may continue, but the result is a bounded evidence snapshot rather than
an exhaustive literature search.

When \code{allow_model_calls = FALSE}, the workflow pauses before the first
new model request. Calling \code{InterMetaAIResume()} with the same default
continues only through stages that need no model call. If a previous process
ended during a model stage, the stage is marked delivery-uncertain because the
provider may have received and billed the request. No retry occurs unless both
\code{allow_model_calls = TRUE} and \code{retry_failed = TRUE} are supplied.
Returned screening and synthesis responses are retained unchanged as
pre-validation artifacts, separately from their local validation results. When
a compatible package changes only a local validation revision, resume can write
a new immutable attempt by revalidating the saved response, with zero new model
calls and no paid retry authorization.

The saved configuration contains provider/model labels, settings, and the
names of credential environment variables, but never credential values. The
provider resolves credentials at invocation. Do not put keys in a script,
provider object, report, or package source.

Every completed stage is bound to an input and artifact SHA-256. Compatible
completed work is reused; an altered or corrupt artifact is not silently
accepted. \code{InterMetaAIReplay()} uses the saved context, verified evidence,
and parsed synthesis to rebuild a report and figure under the current local
validators, writing only to a new output directory.

New workflow configurations store the ordinary sibling PubMed cache as a
run-parent-relative path. If an older run is moved, its conventional sibling
\code{pubmed_cache} is resolved at the new location when possible; a moved
unfinished run with a nonstandard cache can instead supply
\code{pubmed_cache_dir} to \code{InterMetaAIResume()}. Neither compatibility
mapping nor an operational override rewrites the saved configuration. Replay
metadata records its source relative to the replay directory. These path-only
changes do not alter searches, model requests, validation, or published
scientific results, and completed InterMeta 1.1.1 runs remain usable.
Because the legacy configuration did not record whether a directory named
code{pubmed_cache} was the default or a same-named custom location, supply an
explicit resume override when that rare distinction matters.

The workflow-stage synthesis aggregate stores checkpoint references as
code{checkpoints/calls/...} relative to the directory containing
code{05_synthesis/synthesis.rds}, identified by
code{checkpoint_path_base = "synthesis_artifact_dir"}. The lower-level
synthesis aggregate inside code{05_synthesis/checkpoints} retains references
relative to its own code{outdir}. Replay validates the saved workflow
configuration SHA before using its publication settings.

Within an ordinary workflow run, each regenerated report and figure is written
to its own attempt directory under
\code{06_report/outputs/attempt-NNN} and
\code{07_figure/outputs/attempt-NNN}, respectively. This prevents a later
revalidation attempt from reusing stale publication files.

All conclusions remain exploratory and evidence-bounded. The workflow does not
establish causality, pathway activity, clinical utility, or independent
validation of the prognostic score direction.
}
\seealso{
\code{InterMetaAIContext()}, \code{InterMetaBuildPubMedQueries()},
\code{InterMetaAIScreenEvidence()}, \code{InterMetaAISynthesizeScreened()},
\code{InterMetaAIExplorationFigure()}
}
\examples{
\dontrun{
data("InterMetaPRADScore")
context <- InterMetaAIContext(
  score = InterMetaPRADScore,
  disease = "prostate adenocarcinoma",
  endpoint = "biochemical recurrence-free survival",
  time_unit = "days",
  cohort = "TCGA-PRAD",
  top_n = 5L,
  direction = "unfavorable"
)

# Set DEEPSEEK_API_KEY outside this script before starting R.
run <- InterMetaAIWorkflow(
  context,
  run_dir = "InterMeta_AI_PRAD_1_1_3",
  allow_model_calls = TRUE,
  language = "zh"
)

InterMetaAIStatus(run$run_dir)

# Ordinary continuation reuses completed work and makes no failed-call retry.
run <- InterMetaAIResume(
  run$run_dir,
  allow_model_calls = TRUE
)

# A completed run can be revalidated and redrawn with no external request.
replay <- InterMetaAIReplay(
  run$run_dir,
  outdir = "InterMeta_AI_PRAD_replay"
)
}
}
