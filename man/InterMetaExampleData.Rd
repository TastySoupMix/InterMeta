\name{InterMetaExampleData}
\alias{InterMetaExampleData}
\alias{InterMetaExampleInfo}
\alias{InterMetaExampleRNA}
\alias{InterMetaExamplePathwayset}
\alias{InterMetaExampleCox}
\alias{InterMetaExampleEnrich}
\alias{InterMetaExampleScore}
\docType{data}
\title{Synthetic Example Data for the InterMeta Prognosis Workflow}
\description{
A compact, reproducible synthetic cohort and its reference results. The data
are intended for learning, testing an installation, and checking an analysis
workflow. They must not be used as biological or clinical evidence.
}
\usage{
data("InterMetaExampleInfo")
data("InterMetaExampleRNA")
data("InterMetaExamplePathwayset")
data("InterMetaExampleCox")
data("InterMetaExampleEnrich")
data("InterMetaExampleScore")
}
\format{
Six objects are provided:
\describe{
  \item{\code{InterMetaExampleInfo}}{A data frame with 120 samples as row
  names. The columns are binary event status and positive follow-up time.}
  \item{\code{InterMetaExampleRNA}}{A numeric matrix containing normalized
  expression for 30 genes (rows) across 120 samples (columns).}
  \item{\code{InterMetaExamplePathwayset}}{A data frame assigning five genes
  to each of six metabolic pathways.}
  \item{\code{InterMetaExampleCox}}{The 435 gene-pair Cox results generated
  from the example cohort.}
  \item{\code{InterMetaExampleEnrich}}{Directed pathway-interaction
  enrichment results, including retained shared gene pairs and 100-iteration
  permutation FDR values. This object can be used directly with
  \code{InterMetaGeneChord()}.}
  \item{\code{InterMetaExampleScore}}{The pathway-level favorable and
  unfavorable prognostic score table.}
}
}
\details{
The expression and survival values were generated with the fixed seed
\code{20260728}. A latent risk variable was used to create an intentionally
strong, ordered pathway signal so that a small example produces non-empty
results quickly. The complete generation procedure is available in
\file{scripts/generate_example_data.R} in the installed package and in
\file{inst/scripts/generate_example_data.R} in the source tarball. The
development repository also keeps a matching \file{data-raw} copy.
}
\source{
Synthetic data generated for InterMeta; no patient records are included.
}
\keyword{datasets}
