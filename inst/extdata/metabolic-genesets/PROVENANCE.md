# Metabolic gene-set sources

Source record updated: 2026-11-10. The data files and their hashes are unchanged.

| File | Source and preparation |
| --- | --- |
| geneset-kegg.txt | KEGG metabolism annotations converted to long pathway-gene format. |
| geneset-reactome.txt | Reactome metabolism annotations converted to long pathway-gene format. |
| geneset-all.txt | KEGG-derived energy-metabolism pathways selected, organized and merged by Zeheng Tan. |
| geneset-glucose.txt | Author-curated glucose-metabolism subset of the KEGG-derived collection. |
| geneset-lipid.txt | Author-curated lipid-metabolism subset of the KEGG-derived collection. |
| geneset-amino.txt | Author-curated amino-acid-metabolism subset of the KEGG-derived collection. |

`CellType` is a historical column name and denotes a metabolic pathway here.
The filename is `geneset-glucose.txt`, not `geneset-gluscose.txt`.

## Source references and reproducibility

- KEGG: https://www.kegg.jp/ ; https://www.kegg.jp/kegg/docs/intro.html
- Reactome: https://reactome.org/ ; https://reactome.org/cite-us
- scMetabolism reference implementation: https://github.com/wu-yc/scMetabolism
- MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/

Local comparison found 2,913 pathway-entry pairs shared with the scMetabolism
0.2.1 KEGG GMT and 5,639 with its Reactome GMT, after pathway-label and
entry-punctuation conversion. The matching Reactome GMT contains MSigDB card
URLs. This documents a matching resource; it does not establish which site was
the immediate download route. The original database releases, download dates
and a complete scripted reconstruction of the manual curation were not recorded.
The table below identifies the exact files used by this package.

The GPL-3 license applies to InterMeta code; source-database data retain their
applicable terms. Reactome data are CC0. Source terms are available at
https://www.kegg.jp/kegg/legal.html, https://reactome.org/license and
https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp.
This source record does not grant additional database rights.

## Current file fingerprints

SHA-256 fingerprints of the pre-release files:

| File | SHA-256 |
| --- | --- |
| geneset-all.txt | 5c11db0d40107e745bfe0322a418211bbaa98e819d24eed13f0ab84b3e15ffb3 |
| geneset-amino.txt | bc63316b856902ac36af18306d2ad40372524b5862a69c5c6f068ef4ec549b63 |
| geneset-glucose.txt | 19ebccdde496744901e73fdc7f2ae9b2fe9f715ccb6752ddcbad241626f733f1 |
| geneset-kegg.txt | 7d149a24e762c563785f8d1361f0ba068b5071e47a0dcaa4104a6ea6c6faf592 |
| geneset-lipid.txt | be458f9fd7fc867a3682b20a3974005dd8a84e2648daed65697ae326f33d8d14 |
| geneset-reactome.txt | 1c4ce1707ff15fd6258c2186704467c642ad838f769072172125d97b29491af0 |
