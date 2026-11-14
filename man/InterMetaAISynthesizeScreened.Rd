\name{InterMetaAISynthesizeScreened}
\alias{InterMetaAISynthesizeScreened}
\title{Synthesize Locally Verified PubMed Evidence by Pathway}
\usage{
InterMetaAISynthesizeScreened(
  context,
  screening,
  provider = NULL,
  language = c("en", "zh"),
  outdir = NULL,
  max_findings_per_pathway = 5L,
  reuse = TRUE,
  retry_failed = FALSE,
  continue_on_error = TRUE
)
}
\arguments{
  \item{context}{An \code{InterMetaAIContext} object.}
  \item{screening}{A usable \code{InterMetaAIScreening} object with status
  \code{"complete"} or \code{"complete_with_quarantine"}, or an
  \code{InterMetaAIScreenedEvidence} object produced from it. Quarantined rows
  are excluded before this stage.}
  \item{provider}{An optional \code{InterMetaAIProvider}. When \code{NULL},
  a DeepSeek provider is constructed and its credential is read only at model
  invocation time.}
  \item{language}{Language used for the human-readable synthesis.}
  \item{outdir}{Optional directory for attempt-preserving pathway checkpoints,
  the merged interpretation, and pathway status table.}
  \item{max_findings_per_pathway}{Integer from 1 through 5 giving the maximum
  number of verified candidate findings selected for each pathway.}
  \item{reuse}{Logical; reuse a compatible successful pathway checkpoint.}
  \item{retry_failed}{Logical; explicitly permit one new attempt for a
  previously failed pathway. There is no automatic model retry.}
  \item{continue_on_error}{Logical; preserve completed pathway results and
  continue after one pathway fails.}
}
\value{
An \code{InterMetaAISynthesisRun} containing \code{interpretation}, the
validated merged \code{InterMetaAIInterpretation}; \code{pathway_status};
normalized \code{pathway_responses}; checkpoint references relative to
\code{outdir}, identified by \code{checkpoint_path_base = "outdir"};
provenance; and warnings. The \code{complete} field is true when every pathway
was either synthesized or had no locally verified evidence.
}
\description{
Run the second, independent model stage of the InterMeta two-stage workflow. It
uses only articles that passed first-stage screening and creates one bounded
biological synthesis per selected pathway.
}
\details{
The provider receives structured verified findings and may select only existing
candidate identifiers. It reports an evidence-bounded biological-plausibility
category, summary, context dependence, research hypothesis, limitations, and
up to five selected candidate identifiers.

The second model stage is not asked to reproduce abstracts. R maps every
selected identifier back to the frozen first-stage record and revalidates the
stored quotation when the merged interpretation is built. This removes a
second opportunity for citation transcription drift while keeping the article
and decision trail auditable.

Biological plausibility is a model-assisted synthesis of the retrieved,
screened subset. It is not a causal result, an evidence grade, proof that the
literature search is exhaustive, or independent validation of the sign of an
InterMeta prognostic score. The cited papers still require expert review.

Provider calls may incur charges even when transport or validation fails.
Every attempt is preserved when \code{outdir} is supplied, and another paid
attempt is made only when \code{retry_failed = TRUE} is explicit. The returned
parsed provider response is retained unchanged as the pre-validation artifact,
separately from its locally normalized result and including when local
validation fails. If a compatible later package build changes only the
synthesis validation revision, that saved response is revalidated into a new
immutable attempt checkpoint with zero new model calls.

When \code{outdir} is supplied, returned checkpoint references are relative to
that directory and remain valid after the complete output tree is moved.
Internal reads, writes, validation, and hashing continue to use resolved paths.
Legacy aggregate objects containing absolute checkpoint references remain
readable; compatible checkpoint discovery is based on the checkpoint files,
not on stale aggregate paths. In \code{pathway_status},
\code{Model.Call.Attempted} describes the current function invocation and is
therefore false when a checkpoint is reused;
\code{Paid.May.Have.Occurred} retains the historical billing-risk flag stored
with that checkpoint.
}
\seealso{
\code{InterMetaAIScreenEvidence()}, \code{InterMetaAIReport()},
\code{InterMetaAIWorkflow()}
}
\examples{
\dontrun{
synthesis <- InterMetaAISynthesizeScreened(
  context,
  screening,
  outdir = "ai_run/synthesis",
  max_findings_per_pathway = 5L
)

synthesis$pathway_status
synthesis$interpretation$exploration_summary
}
}
