\name{InterMetaPubMed}
\alias{InterMetaBuildPubMedQueries}
\alias{InterMetaPubMedSearch}
\title{Build and Run Traceable PubMed Searches}
\usage{
InterMetaBuildPubMedQueries(
  context,
  disease_mesh = NULL,
  pathway_synonyms = NULL,
  max_genes = 5L,
  pathway_components = NULL,
  include_components = TRUE,
  purpose = c("biology_exploration", "prognostic_validation")
)

InterMetaPubMedSearch(
  context,
  queries = NULL,
  retmax = 20L,
  email = NULL,
  tool = "InterMeta",
  api_key = Sys.getenv("NCBI_API_KEY"),
  cache_dir = NULL,
  refresh = FALSE,
  request_fun = NULL,
  sleep_fun = Sys.sleep,
  page_size = NULL,
  sort = NULL,
  purpose = NULL
)
}
\arguments{
  \item{context}{An \code{InterMetaAIContext} object, or a compatible list
  containing \code{analysis$disease}, \code{analysis$endpoint}, and a
  \code{pathways} table with \code{Pathway} and \code{Direction}.}
  \item{disease_mesh}{Optional disease MeSH term or disease-to-term mapping.}
  \item{pathway_synonyms}{Optional pathway-to-synonym mapping supplied as a
  named list, named character vector, or data frame. Synonyms stored in
  \code{context$pathways$Synonyms} are also used automatically.}
  \item{max_genes}{Maximum number of pathway driver genes included in an
  optional gene-bridge query.}
  \item{pathway_components}{Optional pathway-to-component mapping. Component
  terms are retrieval expansions, not strict synonyms of the whole pathway.}
  \item{include_components}{Logical; include separate component retrieval
  layers. Defaults to \code{TRUE}.}
  \item{purpose}{The default query purpose is \code{biology_exploration}:
  broad disease/pathway, biological-process, and review queries, plus separate
  component queries. \code{prognostic_validation} retains the version 2.1 plan.
  Search infers purpose from supplied queries; an explicit purpose must match.}
  \item{queries}{Optional table produced by
  \code{InterMetaBuildPubMedQueries()}. If omitted, queries are built from
  \code{context}.}
  \item{retmax}{Positive integer maximum number of PubMed identifiers retrieved
  per query. The default is deliberately bounded at 20.}
  \item{email}{Optional contact email supplied to NCBI E-utilities.}
  \item{tool}{Tool label supplied to NCBI E-utilities.}
  \item{api_key}{Optional NCBI API key. By default it is read at call time from
  \code{NCBI_API_KEY}; it is never stored in the returned object or cache.}
  \item{cache_dir}{Optional directory for query-level response caches.}
  \item{refresh}{Logical; ignore valid cache entries when \code{TRUE}.}
  \item{request_fun}{Optional two-argument request function for controlled
  testing. Ordinary use should leave this as \code{NULL}.}
  \item{sleep_fun}{Rate-limiting function. Ordinary use should retain
  \code{Sys.sleep}.}
  \item{page_size}{Optional positive integer page size. Activates paged
  ESearch with a total per-query \code{retmax} cap of 10000.}
  \item{sort}{Optional \code{relevance} or \code{pub_date}. Either this or
  \code{page_size} activates the separate version 3.1 paged cache. Leaving both
  NULL retains the legacy single-page transport and cache behavior.}
}
\value{
\code{InterMetaBuildPubMedQueries()} returns a deterministic query data frame.
\code{InterMetaPubMedSearch()} returns an \code{InterMetaPubMedEvidence}
object containing queries, per-query retrieval status, article records,
response hashes, and provenance.
}
\description{
Builds field-tagged disease, pathway, prognosis, review, and optional driver
gene queries, then retrieves article metadata through NCBI E-utilities only
when \code{InterMetaPubMedSearch()} is called explicitly.
}
\details{
Query version 3.1 defaults to exploratory tumour biology without a prognosis
filter. The clinical endpoint remains context, not an admission criterion.
Version 2.1 is retained for explicit prognostic validation. Both use a conservative synonym/component
dictionary. Pathway identity and scores are never changed by query expansion.
Bare, nonspecific single-amino-acid component terms are avoided where a
pathway-specific intermediate can be used; article screening remains the
record-level relevance gate.
BCR expands to biochemical recurrence and biochemical relapse; unknown
endpoint labels are not included literally. Queries, searches, and records
retain \code{term_scope}, \code{search_terms}, \code{query_version}, and
\code{query_warning}. Legacy supplied queries remain labelled legacy.
Component-level retrieval cannot itself establish whole-pathway direction.

API failures remain \code{NA} with an error status and are never represented as
zero evidence. Successful zero-hit queries retain a genuine count of zero.
The search table records whether the configured \code{retmax} truncated a
query. Retracted publications are excluded from generated queries.

Requests have bounded timeout and retry behavior. Cache entries distinguish
the live NCBI transport from custom test transports and validate stored hashes
and parsed records before reuse. Response bodies are redacted for the API key
before storage; \code{body_redacted} and \code{response_hash_scope} identify
this explicitly. The absolute cache path is represented by a hash in
provenance, not included as plain text.

The optional packages \pkg{httr2}, \pkg{jsonlite}, and \pkg{xml2} are required
for live retrieval. Query construction itself does not require network access.
}
\examples{
data("InterMetaExampleScore")
context <- InterMetaAIContext(
  score = InterMetaExampleScore,
  disease = "example cancer",
  endpoint = "overall survival",
  direction = "unfavorable"
)
queries <- InterMetaBuildPubMedQueries(context)
head(queries)

\dontrun{
evidence <- InterMetaPubMedSearch(
  context,
  queries = queries,
  email = Sys.getenv("NCBI_EMAIL"),
  cache_dir = "pubmed_cache"
)
}
}
