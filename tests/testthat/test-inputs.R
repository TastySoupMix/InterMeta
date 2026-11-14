test_that("bundled genesets load and use pathway terminology", {
  available <- InterMetaAvailableGenesets()
  expect_setequal(
    available,
    c("glucose", "lipid", "amino", "all", "kegg", "reactome")
  )
  geneset <- InterMetaLoadGeneset("all")
  expect_named(geneset, c("Pathway", "Gene", "Dataset"))
  expect_equal(length(unique(geneset$Pathway)), 52)
  expect_equal(length(unique(geneset$Gene)), 1292)
  expect_equal(nrow(geneset), 1999L)
  expect_equal(
    nrow(InterMeta:::.InterMetaLegacyPermutationGeneset(geneset)),
    2022L
  )
  expect_false(anyDuplicated(geneset) > 0)
})

test_that("legacy preprocessing allows explicit no-log input and rejects an empty safe set", {
  geneset <- data.frame(
    Pathway = c("A", "A", "B"),
    Gene = c("G1", "G2", "G3"),
    Dataset = "toy"
  )
  rna <- matrix(
    c(
      -2, -1, 0, 1, 2, 3,
      3, 2, 1, 0, -1, -2,
      -1, 0, 1, 2, 3, 4
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("G1", "G2", "G3"), paste0("S", 1:6))
  )
  processed <- suppressMessages(InterMetaPreprocess(
    pathwayset = geneset,
    rna = rna,
    log = FALSE,
    GMNorm = TRUE
  ))
  expect_true(all(is.finite(processed)))
  expect_equal(as.numeric(apply(processed, 1, median)), rep(0, 3))
  expect_error(
    InterMetaPreprocess(
      pathwayset = geneset,
      rna = rna,
      log = TRUE
    ),
    "All genes would produce non-finite"
  )
})

test_that("clinical input is validated explicitly", {
  info <- data.frame(
    event = c(0, 1, 0, 1, 0, 1),
    time = 1:6,
    row.names = paste0("S", 1:6)
  )
  expect_equal(
    suppressMessages(InterMetaCheckEvent(info, min_events = 2)),
    info
  )
  invalid <- info
  invalid$event[1] <- 2
  expect_error(InterMetaCheckEvent(invalid), "0.*1")
})
