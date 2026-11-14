# Inspect the bundled real TCGA-PRAD cohort and precomputed results.

library(InterMeta)

data("InterMetaPRADInfo")
data("InterMetaPRADRNA")
data("InterMetaPRADEnrich")
data("InterMetaPRADScore")
data("InterMetaPRADProvenance")

pathwayset <- InterMetaLoadGeneset("all")

## 1. Validate the real input data --------------------------------------------

info <- InterMetaCheckEvent(
  InterMetaPRADInfo,
  min_events = 5,
  min_censored = 5
)
rna <- InterMetaPreprocess(
  pathwayset = pathwayset,
  rna = InterMetaPRADRNA,
  cohort = row.names(info),
  log = FALSE,
  GMNorm = TRUE
)

dim(info)
table(info$event)
dim(rna)
InterMetaPRADProvenance

## 2. Inspect the migrated enrichment and score results -----------------------

head(InterMetaPRADEnrich)
InterMetaPRADScore

dotplot <- InterMetaDotplot(
  InterMetaPRADEnrich,
  condition = "Permutation.FDR",
  cutoff = 0.2,
  top = 10
)
score_bar <- InterMetaPathwayScoreBar(InterMetaPRADScore)
network <- InterMetaPathwayNetwork(
  InterMetaPRADEnrich,
  pathwayset = pathwayset,
  condition = "Permutation.FDR",
  cutoff = 0.2,
  export = FALSE
)
chord <- InterMetaPathwayChord(
  InterMetaPRADEnrich,
  condition = "Permutation.FDR",
  cutoff = 0.2,
  draw = interactive()
)

print(dotplot)
print(score_bar)
print(chord)
head(network$network)
