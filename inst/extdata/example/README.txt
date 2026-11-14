InterMeta synthetic example files

These files contain generated data only. They do not contain patient records
and must not be interpreted as biological or clinical evidence.

conditions.txt
  Tab-delimited clinical data. Column 1 is the sample ID, column 2 is event
  status (0 = censored, 1 = event), and column 3 is positive follow-up time.

exp_tumor.txt
  Tab-delimited normalized expression. Column 1 is the gene symbol, and each
  remaining column is one sample. Use log = FALSE with these values.

geneset-example.txt
  Tab-delimited Pathway, Gene, and Dataset columns for six metabolic pathways.

The same content is available after library(InterMeta) through:
  data("InterMetaExampleInfo")
  data("InterMetaExampleRNA")
  data("InterMetaExamplePathwayset")

See system.file("examples", package = "InterMeta") for runnable scripts.
