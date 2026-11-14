test_that("bundled example inputs are coherent and runnable", {
  data("InterMetaExampleInfo", package = "InterMeta", envir = environment())
  data("InterMetaExampleRNA", package = "InterMeta", envir = environment())
  data("InterMetaExamplePathwayset", package = "InterMeta",
       envir = environment())
  data("InterMetaExampleEnrich", package = "InterMeta",
       envir = environment())
  data("InterMetaExampleScore", package = "InterMeta",
       envir = environment())

  expect_equal(dim(InterMetaExampleInfo), c(120L, 2L))
  expect_equal(dim(InterMetaExampleRNA), c(30L, 120L))
  expect_setequal(colnames(InterMetaExampleRNA),
                  row.names(InterMetaExampleInfo))
  expect_equal(length(unique(InterMetaExamplePathwayset$Pathway)), 6L)
  expect_true(any(InterMetaExampleEnrich$Adjust.P.Value < 0.05))
  expect_true(all(is.finite(InterMetaExampleEnrich$Permutation.FDR)))
  expect_equal(min(InterMetaExampleEnrich$Permutation.FDR), 0.09)
  expect_gte(nrow(InterMetaExampleScore), 2L)

  info <- InterMetaCheckEvent(InterMetaExampleInfo)
  rna <- InterMetaPreprocess(
    pathwayset = InterMetaExamplePathwayset,
    rna = InterMetaExampleRNA,
    cohort = row.names(info),
    log = FALSE,
    GMNorm = TRUE
  )
  mps <- InterMetaGenePair(rna)
  expect_equal(nrow(mps), 435L)

  cox <- InterMetaCox(
    mps = mps,
    info = info,
    p.adj = "BH",
    parallel = FALSE,
    keep_mps = FALSE
  )
  # The bundled enrichment object is the explicit adjusted-QV fixture; the
  # public InterMetaStability() default is tested separately as PV.
  query <- row.names(cox$cox_res)[
    is.finite(cox$cox_res$QV) & cox$cox_res$QV < 0.05
  ]
  expect_gte(length(query), 10L)

  enrich <- InterMetaEnrich(
    gene = query,
    background = InterMetaBG(row.names(cox$cox_res)),
    pathway.pair = InterMetaPathwayPair(InterMetaExamplePathwayset),
    p.adj = "BH"
  )
  expect_true(any(enrich$Adjust.P.Value < 0.05))
})

test_that("bundled example text files are installed", {
  example_dir <- system.file("extdata", "example", package = "InterMeta")
  expect_true(nzchar(example_dir))
  expect_true(all(file.exists(file.path(
    example_dir,
    c("conditions.txt", "exp_tumor.txt", "geneset-example.txt", "README.txt")
  ))))
})

test_that("dataset generation and provenance scripts are installed", {
  script_dir <- system.file("scripts", package = "InterMeta")
  expect_true(nzchar(script_dir))
  expect_true(all(file.exists(file.path(
    script_dir,
    c(
      "generate_example_data.R",
      "generate_prad_data.R",
      "update_prad_provenance.R",
      "update_prad_provenance_impl.R"
    )
  ))))
})
