# Validate and publish PRAD release provenance from an InterMeta source tree.
# Usage: Rscript data-raw/update_prad_provenance.R <output-dir> <tarball>
# or:    Rscript inst/scripts/update_prad_provenance.R <output-dir> <tarball>
# The public implementation is shipped in inst/scripts for source-tarball use.

implementation <- file.path("inst", "scripts", "update_prad_provenance_impl.R")
if (!file.exists(implementation)) {
  stop("Run this script from an InterMeta source-tree root.")
}
source(implementation, local = globalenv())
