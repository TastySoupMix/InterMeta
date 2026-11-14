# Inspect InterMeta reference results without rerunning Cox or permutation.

library(InterMeta)

data("InterMetaExamplePathwayset")
data("InterMetaExampleCox")
data("InterMetaExampleEnrich")
data("InterMetaExampleScore")

head(InterMetaExampleCox)
head(InterMetaExampleEnrich)
InterMetaExampleScore

dotplot <- InterMetaDotplot(
  InterMetaExampleEnrich,
  condition = "Permutation.FDR",
  cutoff = 0.05,
  top = 10
)
score_bar <- InterMetaPathwayScoreBar(InterMetaExampleScore)
network <- InterMetaPathwayNetwork(
  InterMetaExampleEnrich,
  pathwayset = InterMetaExamplePathwayset,
  condition = "Permutation.FDR",
  cutoff = 0.05,
  export = FALSE
)
chord <- InterMetaPathwayChord(
  InterMetaExampleEnrich,
  condition = "Permutation.FDR",
  cutoff = 0.05,
  draw = interactive()
)
gene_chord <- InterMetaGeneChord(
  InterMetaExampleEnrich,
  select = InterMetaExampleEnrich$Index[1],
  condition = "Permutation.FDR",
  draw = interactive()
)

print(dotplot)
print(score_bar)
print(chord)
print(gene_chord)
head(network$network)
