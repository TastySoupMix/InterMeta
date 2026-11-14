\name{InterMetaAIExplorationFigure}
\alias{InterMetaAIExplorationFigure}
\title{Render the Fixed InterMeta AI Exploration Figure}
\usage{
InterMetaAIExplorationFigure(
  report,
  pdf_file = NULL,
  png_file = NULL,
  render = TRUE,
  font_family = NULL,
  text_mode = c("compact_english", "full")
)
}
\arguments{
  \item{report}{A validated \code{InterMetaAIReport} whose purpose is
  \code{"biology_exploration"}. Its \code{top_pathways} must contain one to
  five pathways from one prognostic direction.}
  \item{pdf_file}{Optional path to a new PDF file.}
  \item{png_file}{Optional path to a new PNG file.}
  \item{render}{Logical; render requested files. Set to \code{FALSE} to
  validate and inspect the prepared figure data without opening a graphics
  device.}
  \item{font_family}{\code{NULL}, or one installed font-family name. ASCII-only
  output uses the ordinary \code{sans} device font. Non-ASCII output and
  custom fonts require Cairo, the optional \pkg{systemfonts} package, and
  complete glyph coverage in regular and bold faces.}
  \item{text_mode}{\code{"compact_english"}, the publication default, or
  \code{"full"}. Compact mode leaves the complete narrative in the report and
  displays a short source-count summary.}
}
\value{
Invisibly, a list containing validated display tables, caption, calculated
dimensions, font and evidence-source audits, output paths, and a local drawing
function. When \code{render = TRUE}, the requested PDF and/or PNG is created
only after validation succeeds.
}
\description{
Render a fixed, publication-style summary of a validated tumour-biology
exploration report. Rendering is entirely local and performs no PubMed
retrieval or model request.
}
\details{
The figure combines the selected InterMeta score direction with evidence
channels, biological processes, the evidence-bounded biological-plausibility
label, grounding/eligibility state, and source-linked references. Functional,
omics, and patient-study counts are non-exclusive unique-PMID counts; they are
not evidence grades, quality scores, or independent replication counts.

The plausibility label summarizes only the retrieved and locally verified
literature subset. It is model-assisted, does not establish causal importance,
and does not independently validate the prognostic score direction. Read the
full report and original articles before using the figure in a publication.

Output paths must not already exist. This prevents accidental replacement of a
previous figure. For a report containing both prognostic directions, use
\code{InterMetaAIWorkflow()} or create one direction-specific report per
figure.
}
\seealso{
\code{InterMetaAIReport()}, \code{InterMetaAIWorkflow()},
\code{InterMetaAIReplay()}
}
\examples{
\dontrun{
# report is a completed biology_exploration InterMetaAIReport.
figure_result <- InterMetaAIExplorationFigure(
  report,
  pdf_file = "InterMeta_AI_exploration.pdf",
  png_file = "InterMeta_AI_exploration.png",
  text_mode = "compact_english"
)
}
}
