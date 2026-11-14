\name{InterMetaAIScreening}
\alias{InterMetaAIScreenEvidence}
\alias{InterMetaAIScreenedEvidence}
\title{Screen PubMed Evidence Before Pathway Synthesis}
\usage{
InterMetaAIScreenEvidence(
  context,
  evidence,
  provider = NULL,
  language = c("en", "zh"),
  batch_size = 20L,
  checkpoint_dir = NULL,
  resume = TRUE,
  retry_failed = FALSE,
  continue_on_error = TRUE
)

InterMetaAIScreenedEvidence(screening, require_complete = TRUE)
}
\arguments{
  \item{context}{An \code{InterMetaAIContext} object containing only the
  aggregate pathway-analysis context.}
  \item{evidence}{An \code{InterMetaPubMedEvidence} object returned by
  \code{InterMetaPubMedSearch()}.}
  \item{provider}{An optional \code{InterMetaAIProvider}. When \code{NULL},
  a DeepSeek provider is constructed; its credential is resolved only when a
  non-empty model batch is invoked.}
  \item{language}{Language used for human-readable screening annotations.}
  \item{batch_size}{Positive integer maximum number of candidates in one
  single-pathway provider call.}
  \item{checkpoint_dir}{Optional directory in which attempt-preserving batch
  checkpoints are stored.}
  \item{resume}{Logical; reuse checkpoints whose frozen input and provider
  bindings match the current call.}
  \item{retry_failed}{Logical; explicitly permit a new attempt for a failed or
  partial checkpoint. There is no automatic retry within an attempt.}
  \item{continue_on_error}{Logical; continue to later batches after a provider
  or whole-response structural failure.}
  \item{screening}{An \code{InterMetaAIScreening} object.}
  \item{require_complete}{Logical; require that no whole provider batch failed
  or remained incomplete. The default accepts screening status
  \code{"complete"} or \code{"complete_with_quarantine"}, but refuses a
  failed or incomplete whole batch. Individually quarantined candidates never
  enter the returned evidence.}
}
\value{
\code{InterMetaAIScreenEvidence()} returns an \code{InterMetaAIScreening}
object. Its \code{included}, \code{excluded}, and \code{unresolved} tables are
mutually exclusive candidate-level audit sets, and \code{batch_status} records
each model attempt.

\code{InterMetaAIScreenedEvidence()} returns an
\code{InterMetaAIScreenedEvidence} object inheriting from
\code{InterMetaPubMedEvidence}. Its \code{records} contain only locally
validated included candidates. The original query hits remain available as
\code{candidate_records}, and the complete screening audit is retained in
\code{screening}.
}
\description{
Run the first model stage of the InterMeta two-stage evidence workflow. Candidate
articles are screened in deterministic, pathway-specific batches, then the
structured decisions are checked locally before any article can reach pathway
synthesis.
}
\details{
Candidate identity is the pathway--PMID pair. Local validation checks candidate
membership, allowed structured labels, lexical consistency of the supplied
target-disease label, and contiguous Abstract quotations under the documented
whitespace and Unicode-typography rules. These checks establish structural and
source grounding only. They do not prove relevance, study design, effect, or
biological truth.

An invalid article row is quarantined as \code{"unresolved"} instead of being
silently repaired or invalidating otherwise usable rows. If no whole provider
batch failed and one or more candidate rows remain quarantined, the aggregate
screening status is \code{"complete_with_quarantine"}; this status remains eligible for
\code{InterMetaAIScreenedEvidence(..., require_complete = TRUE)}. Quarantined
rows are retained in the audit but never enter synthesis. A failed or
incomplete whole provider batch instead leaves the screening \code{"partial"}
or \code{"failed"}, which the default conversion refuses. With
\code{continue_on_error = TRUE}, later batches can still finish and be
preserved for inspection.

Provider calls can incur charges. A timed-out or failed request may still have
reached the service. Failed attempts are preserved, and a further call requires
\code{retry_failed = TRUE}; the function never retries a paid request
automatically. Provider credentials are not stored in the returned object or
checkpoint configuration. The returned raw structured provider response is
retained unchanged in its own checkpoint field before any locally validated
decision is accepted. If only the local validation revision changes in a
compatible later build, the saved response can be revalidated into a new
immutable attempt checkpoint with zero new model calls.
}
\seealso{
\code{InterMetaPubMedSearch()}, \code{InterMetaAISynthesizeScreened()},
\code{InterMetaAIWorkflow()}
}
\examples{
\dontrun{
# context and evidence were created by InterMetaAIContext() and
# InterMetaPubMedSearch(), respectively.
screening <- InterMetaAIScreenEvidence(
  context,
  evidence,
  checkpoint_dir = "ai_run/screening",
  batch_size = 20L
)

screening$status
screening$included
screening$unresolved

# complete_with_quarantine is usable: unresolved rows remain audit-only.
# A failed or incomplete whole provider batch is deliberately refused.
screened_evidence <- InterMetaAIScreenedEvidence(screening)
}
}
