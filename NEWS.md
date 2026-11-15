# InterMeta 1.1.4

- Uses TimiGP-legacy as the default statistical implementation and removes
  alternative statistical branches.
- Makes preprocessing, ranking, permutation and three-resample defaults explicit.
  Preprocessing remains configurable for the input expression scale.
- Updates bundled TCGA-PRAD pathway scores from the existing enrichment results,
  retaining 50 nonzero scores. Synthetic reference results are regenerated using
  TimiGP-legacy.
- Updates examples and help, retaining earlier full-analysis and permutation
  provenance with a score-migration record.
- AI algorithms and historical AI archive metadata are unchanged.
- Keeps developer guidance and the manuscript template in `docs/`, outside the
  installable package.
