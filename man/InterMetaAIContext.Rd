\name{InterMetaAIContext}
\alias{InterMetaAIContext}
\title{Create a Privacy-Preserving Context for Evidence Reporting}
\usage{
InterMetaAIContext(
  enrichment = NULL,
  score = NULL,
  stability = NULL,
  analysis = NULL,
  disease,
  endpoint = "unknown",
  time_unit = "unknown",
  cohort = NULL,
  top_n = 5L,
  direction = c("both", "unfavorable", "favorable"),
  pathway_synonyms = NULL
)
}
\arguments{
  \item{enrichment}{Optional result from \code{InterMetaEnrich()}. A directly
  supplied value takes precedence over \code{analysis$enrich}.}
  \item{score}{Pathway score table from \code{InterMetaPathwayScore()}. A
  directly supplied value takes precedence over \code{analysis$score}.}
  \item{stability}{Optional \code{InterMetaStability} result, list containing a
  \code{summary} table, or stability-summary data frame.}
  \item{analysis}{Optional named full-analysis list containing \code{enrich},
  \code{score}, \code{stability}, and/or \code{parameters}. Patient-level
  fields in this list are ignored.}
  \item{disease}{One non-empty disease label used to construct literature
  queries.}
  \item{endpoint}{One endpoint label, or \code{"unknown"}.}
  \item{time_unit}{One endpoint time-unit label, or \code{"unknown"}.}
  \item{cohort}{Optional scalar cohort label. Vectors of sample identifiers are
  rejected.}
  \item{top_n}{Positive integer number of pathways selected per requested
  direction. The default is five.}
  \item{direction}{Select both directions, unfavorable pathways only, or
  favorable pathways only.}
  \item{pathway_synonyms}{Optional named list, named character vector, or data
  frame with \code{Pathway} and \code{Synonym} columns.}
}
\value{
An object of class \code{InterMetaAIContext}. It contains aggregate analysis
metadata, deterministically ranked pathways, relevant pathway interactions,
driver gene pairs, scientific warnings, and input hashes. It does not retain
expression matrices, survival rows, sample identifiers, or resampling payloads.
}
\description{
Selects the most extreme pathway-score results for an optional, post-analysis
literature and AI reporting layer. Unfavorable pathways are ordered by
increasing \code{Difference}; favorable pathways are ordered by decreasing
\code{Difference}; pathway name provides a deterministic tie break.
}
\details{
This function does not perform a network request or call a language model.
The returned object is a restricted hand-off object intended for
\code{InterMetaBuildPubMedQueries()}, \code{InterMetaPubMedSearch()}, and
\code{InterMetaAIReport()}.
}
\examples{
data("InterMetaExampleEnrich")
data("InterMetaExampleScore")

context <- InterMetaAIContext(
  enrichment = InterMetaExampleEnrich,
  score = InterMetaExampleScore,
  disease = "example cancer",
  endpoint = "overall survival",
  time_unit = "days",
  top_n = 5L,
  direction = "unfavorable"
)
context$pathways
}
