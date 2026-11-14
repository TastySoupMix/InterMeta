\name{InterMetaPathwayChord}
\alias{InterMetaPathwayChord}
\alias{print.InterMetaPathwayChord}
\title{Directional Metabolic Pathway Chord Diagram}
\description{
Draws a directional chord diagram from an InterMeta pathway-interaction
enrichment table. Chords point from favorable pathways to unfavorable
pathways, and their widths encode a selected interaction weight.
}
\usage{
InterMetaPathwayChord(
  resdata = NULL,
  select = NULL,
  condition = "Adjust.P.Value",
  cutoff = 0.05,
  max_edges = NULL,
  link_weight = c("No.Shared.MGP", "Enrichment.Ratio", "significance", "uniform"),
  group = NULL,
  color = NULL,
  order = NULL,
  transparency = 0.4,
  label_cex = 0.7,
  start_degree = 90,
  adaptive_labels = FALSE,
  canvas_limit = NULL,
  sector_floor_quantile = NULL,
  draw = TRUE
)

\method{print}{InterMetaPathwayChord}(x, ...)
}
\arguments{
  \item{resdata}{Data frame returned by \code{InterMetaEnrich()} or
  \code{InterMetaPermFDR()}. It must contain \code{Index},
  \code{Favorable.Pathway}, \code{Unfavorable.Pathway}, and the selected
  significance and link-weight columns.}
  \item{select}{Optional vector of enrichment \code{Index} values. When
  supplied, these interactions are used in the requested order and
  \code{max_edges} is ignored. Unknown values are reported and ignored.}
  \item{condition}{Numeric significance column used when \code{select} is
  \code{NULL} and when \code{link_weight = "significance"}.}
  \item{cutoff}{Strict significance cutoff from 0 to 1.}
  \item{max_edges}{Optional maximum number of the most significant
  interactions drawn when \code{select} is \code{NULL}. The default
  \code{NULL} draws every interaction meeting \code{condition < cutoff},
  matching the full-interaction behavior of the legacy chord plot. Set a
  positive integer only when an explicitly reduced diagram is wanted.}
  \item{link_weight}{Chord-width definition. \code{"No.Shared.MGP"} uses
  shared prognostic gene-pair counts, \code{"Enrichment.Ratio"} uses the
  enrichment ratio, \code{"significance"} uses
  \code{-log10(condition)}, and \code{"uniform"} gives all links equal
  width. Zero probabilities are assigned a finite floor derived from the
  smallest positive selected probability.}
  \item{group}{Optional pathway-to-group mapping. A named vector may include
  pathways not present in the selected interactions; extras are cropped. An
  unnamed vector must contain one value per selected pathway.}
  \item{color}{Optional colors mapped either to selected pathway names or to
  selected group names. An unnamed vector must contain one color per pathway
  or group. A deterministic, pathway-specific color is generated when this is
  \code{NULL}; the same pathway therefore keeps its default color across
  selections.}
  \item{order}{Optional character vector giving circular pathway order. It may
  contain unused pathway names, but it must cover every selected pathway. With
  multiple groups, pathways in each group must remain contiguous.}
  \item{transparency}{Chord transparency from 0 to 1.}
  \item{label_cex}{Positive text-size multiplier for pathway labels. When
  \code{adaptive_labels = TRUE}, this is the maximum permitted label size.}
  \item{start_degree}{Starting angle in degrees.}
  \item{adaptive_labels}{Whether dense diagrams should reduce label size to
  fit the current graphics device and longest pathway name. The default
  \code{FALSE} keeps the classic fixed label size. This never removes a pathway
  or interaction.}
  \item{canvas_limit}{Optional circular canvas extent greater than 1. The
  default \code{NULL} uses the classic extent of 1.45. Set to 1.65 to leave
  additional room around a dense diagram.}
  \item{sector_floor_quantile}{Optional quantile from 0 to 1 used as a minimum
  sector width for diagrams with at least 30 pathways. The default \code{NULL}
  preserves sectors proportional to their link weights without width
  expansion. An explicit value such as 0.75 gives low-weight pathways more
  sector space while retaining every chord.}
  \item{draw}{Whether to draw on the current base graphics device. Set to
  \code{FALSE} to prepare and inspect the chord data without drawing.}
  \item{x}{An object returned by \code{InterMetaPathwayChord()}.}
  \item{...}{Additional arguments ignored by the print method.}
}
\details{
Arrows always follow the InterMeta direction
\code{Favorable.Pathway -> Unfavorable.Pathway}. The function does not use the
legacy \code{dataset} or \code{geneset} arguments. User mappings are reduced to
the pathways that are actually drawn, and unused factor levels are removed
before calling \pkg{circlize}. This avoids failures caused by pathways that
are present in a full annotation but absent from a selected interaction set.

With the default \code{select = NULL} and \code{max_edges = NULL}, all
interactions satisfying the significance cutoff are drawn. Pathways without a
qualifying link are not shown as isolated sectors. To draw every enrichment
row regardless of significance, explicitly use \code{select = resdata$Index}.

The default drawing uses the classic proportional-sector layout, fixed label
size, and a canvas extent of 1.45. Dense-layout adjustments are opt-in through
\code{adaptive_labels = TRUE}, \code{canvas_limit = 1.65}, and
\code{sector_floor_quantile = 0.75}. These controls affect only the drawing;
they do not change selected interactions, weights, colors, arrows, or analysis
results. The returned \code{selected_count} and \code{displayed_count} remain
auditable and unchanged by the layout controls.

The example retains a 14 by 14 inch PDF page to leave room for long pathway
names; restoring the classic layout does not require a smaller export device.

Duplicate directional pathway pairs are combined by summing their weights.
Self-interactions and interactions without a positive finite weight are
reported and omitted. The global circular plotting state is cleared after
drawing, including when plotting fails.
}
\value{
Invisibly returns an object of class \code{InterMetaPathwayChord} containing:
\itemize{
  \item \code{selected}: selected enrichment rows retained for drawing;
  \item \code{edges}: drawable directional links and aggregated weights;
  \item \code{nodes}: pathway groups, colors, degrees, and weighted degrees;
  \item \code{order}, \code{group}, and \code{color}: the exact displayed
  pathway mappings;
  \item \code{condition}, \code{cutoff}, \code{selection}, and
  \code{link_weight}: selection metadata;
  \item \code{eligible_count}, \code{selected_count},
  \code{displayed_count}, \code{max_edges}, and
  \code{significance_floor}: auditable selection and weighting details;
  \item \code{layout}: the layout returned by \code{circlize::chordDiagram()}
  when \code{draw = TRUE}, otherwise \code{NULL}.
}
}
\examples{
\dontrun{
library(InterMeta)
data("InterMetaExampleEnrich")

pdf("pathway_chord.pdf", width = 14, height = 14, useDingbats = FALSE)
chord <- InterMetaPathwayChord(
  InterMetaExampleEnrich,
  condition = "Permutation.FDR",
  cutoff = 0.05,
  max_edges = NULL,
  label_cex = 0.7,
  start_degree = 90,
  adaptive_labels = FALSE,
  canvas_limit = 1.45,
  sector_floor_quantile = NULL
)
dev.off()

print(chord)
}
}
\seealso{
\code{\link{InterMetaEnrich}}, \code{\link{InterMetaPermFDR}},
\code{\link{InterMetaPathwayNetwork}}, \code{\link{InterMetaGeneChord}}
}
\keyword{hplot}
