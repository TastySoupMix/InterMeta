# Load the InterMeta synthetic example as package data or text files.

library(InterMeta)

## Option 1: data objects ------------------------------------------------------

data("InterMetaExampleInfo")
data("InterMetaExampleRNA")
data("InterMetaExamplePathwayset")

dim(InterMetaExampleInfo)
dim(InterMetaExampleRNA)
table(InterMetaExampleInfo$event)
table(InterMetaExamplePathwayset$Pathway)

## Option 2: conditions.txt and exp_tumor.txt --------------------------------

example_dir <- system.file("extdata", "example", package = "InterMeta")

clinical_table <- read.delim(
  file.path(example_dir, "conditions.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
clinical_info <- clinical_table[-1]
row.names(clinical_info) <- clinical_table[[1]]

expression_table <- read.delim(
  file.path(example_dir, "exp_tumor.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
bulk_expression <- as.matrix(expression_table[-1])
storage.mode(bulk_expression) <- "double"
row.names(bulk_expression) <- expression_table[[1]]

pathwayset <- InterMetaLoadGeneset(
  file = file.path(example_dir, "geneset-example.txt")
)

stopifnot(
  isTRUE(all.equal(clinical_info, InterMetaExampleInfo)),
  isTRUE(all.equal(bulk_expression, InterMetaExampleRNA)),
  identical(
    pathwayset[, c("Pathway", "Gene", "Dataset")],
    InterMetaExamplePathwayset
  )
)

message("Bundled data objects and text files are consistent.")
