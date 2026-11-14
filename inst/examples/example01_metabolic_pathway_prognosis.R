# Complete InterMeta workflow using the bundled synthetic prognosis cohort.
# The example is deliberately small and normally finishes within a minute.

library(InterMeta)

data("InterMetaExampleInfo")
data("InterMetaExampleRNA")
data("InterMetaExamplePathwayset")

clinical_info <- InterMetaExampleInfo
bulk_expression <- InterMetaExampleRNA
pathwayset <- InterMetaExamplePathwayset

info <- InterMetaCheckEvent(clinical_info)

rna <- InterMetaPreprocess(
  pathwayset = pathwayset,
  rna = bulk_expression,
  cohort = rownames(info),
  log = FALSE,
  GMNorm = TRUE
)

mps <- InterMetaGenePair(rna, chunk_size = 10000)

cox <- InterMetaCox(
  mps = mps,
  info = info,
  p.adj = "BH",
  parallel = FALSE,
  core = 1,
  keep_mps = FALSE
)

cox_res <- cox$cox_res
query_pairs <- rownames(cox_res)[
  is.finite(cox_res$PV) & cox_res$PV < 0.05
]
background <- InterMetaBG(rownames(cox_res))

pathway_pair <- InterMetaPathwayPair(
  pathwayset = pathwayset,
  core = 1
)

enrich <- InterMetaEnrich(
  gene = query_pairs,
  background = background,
  pathway.pair = pathway_pair,
  p.adj = "BH",
  core = 1,
  keep_gene_pairs = "shared"
)

# The synthetic example is small enough for a 100-iteration permutation test.
enrich <- InterMetaPermFDR(
  resdata = enrich,
  pathwayset = pathwayset,
  gene = query_pairs,
  background = background,
  niter = 100,
  core = 1,
  seed = 20260728
)

net <- InterMetaPathwayNetwork(
  resdata = enrich,
  pathwayset = pathwayset,
  condition = "Permutation.FDR",
  cutoff = 0.05,
  export = FALSE
)

score <- InterMetaPathwayScore(
  enrich,
  condition = "Permutation.FDR",
  cutoff = 0.05
)
p1 <- InterMetaDotplot(
  enrich,
  condition = "Permutation.FDR",
  cutoff = 0.05,
  top = 10
)
p2 <- InterMetaPathwayScoreBar(score)
pathway_chord <- InterMetaPathwayChord(
  enrich,
  condition = "Permutation.FDR",
  cutoff = 0.05,
  draw = interactive()
)
gene_chord <- InterMetaGeneChord(
  enrich,
  select = enrich$Index[1],
  condition = "Permutation.FDR",
  draw = interactive()
)

head(cox_res)
head(enrich)
score
head(net$network)
print(p1)
print(p2)
print(pathway_chord)
print(gene_chord)
