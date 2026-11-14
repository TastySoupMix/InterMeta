ai_figure_fun <- function(...) {
  getFromNamespace("InterMetaAIExplorationFigure", "InterMeta")(...)
}

ai_figure_fixture <- function(n = 2L) {
  pathways <- c(
    "Taurine and hypotaurine metabolism",
    "Oxidative phosphorylation",
    "Inositol phosphate metabolism",
    "Glycerophospholipid metabolism",
    "Lysine degradation"
  )[seq_len(n)]
  top <- data.frame(
    Pathway = pathways,
    Direction = rep("unfavorable", n),
    Direction.Rank = seq_len(n),
    Difference = -seq.int(7, length.out = n),
    stringsAsFactors = FALSE
  )
  finding_rows <- lapply(seq_len(n), function(i) {
    primary <- as.character(11000000L + i)
    review <- as.character(21000000L + i)
    data.frame(
      Pathway = rep(pathways[i], 4L),
      Finding.ID = c("component", "unclear", "pathway", "review"),
      PMID = c(primary, primary, primary, review),
      Scope = c("component", "unclear", "pathway", "unclear"),
      Eligible.For.Summary = rep(TRUE, 4L),
      Quote.Verified = rep(TRUE, 4L),
      Functional = c(TRUE, TRUE, TRUE, FALSE),
      Omics = c(TRUE, FALSE, FALSE, FALSE),
      Patient = c(FALSE, FALSE, TRUE, FALSE),
      Review = c(FALSE, FALSE, FALSE, TRUE),
      stringsAsFactors = FALSE
    )
  })
  findings <- do.call(rbind, c(finding_rows, list(make.row.names = FALSE)))
  summary <- top
  summary$Biological.Summary <- paste(
    "Eligible source-linked findings provide exploratory biological context."
  )
  summary$Biological.Processes <- "proliferation;metabolic_adaptation"
  summary$Context.Dependence <-
    "The observed relationship can depend on model and measurement context."
  summary$Research.Hypothesis <-
    "Test the pathway using an independent perturbation and metabolite readout."
  primary_ids <- as.character(11000000L + seq_len(n))
  review_ids <- as.character(21000000L + seq_len(n))
  summary$Functional.PMIDs <- primary_ids
  summary$Omics.PMIDs <- primary_ids
  summary$Patient.PMIDs <- primary_ids
  summary$Review.PMIDs <- review_ids
  summary$Key.PMIDs <- paste(primary_ids, review_ids, sep = ";")
  summary$Unique.PMID.Count <- rep(1L, n)
  summary$Assessment.Source <- "model_assisted"
  summary$Exploration.Status <- "provided"
  summary$Quote.Grounding <- "complete"
  summary$Evidence.Eligibility <- "complete"
  records <- unique(data.frame(
    pathway = findings$Pathway,
    PMID = findings$PMID,
    stringsAsFactors = FALSE
  ))
  structure(
    list(
      schema_version = "4.2.0",
      purpose = "biology_exploration",
      analysis = list(
        disease = "prostate adenocarcinoma", endpoint = "BCR",
        time_unit = "months"
      ),
      top_pathways = top,
      exploration_summary = summary,
      exploration_findings = findings,
      evidence = list(
        records = records,
        provenance = list(source_kind = "pubmed", database = "PubMed")
      ),
      interpretation = list(
        provider = list(name = "offline-test-provider", model = "fixture-v1")
      )
    ),
    class = c("InterMetaAIReport", "list")
  )
}

test_that("AI exploration figure preparation is offline and non-mutating", {
  report <- ai_figure_fixture()
  original <- report
  before_device <- grDevices::dev.cur()

  figure <- ai_figure_fun(report, render = FALSE)

  expect_identical(report, original)
  expect_identical(grDevices::dev.cur(), before_device)
  expect_null(figure$output$pdf)
  expect_null(figure$output$png)
  expect_true(is.function(figure$draw))
  expect_identical(figure$text_mode, "compact_english")
  expect_identical(figure$pathways$Pathway, report$top_pathways$Pathway)
  expect_identical(
    figure$pathways$Difference, report$top_pathways$Difference
  )
  expect_false(figure$evidence_source$synthetic_fixture)
  expect_true(all(
    figure$pathways$Biological.Plausibility == "not_assessable"
  ))
})

test_that("biological plausibility is validated and displayed", {
  report <- ai_figure_fixture(2L)
  report$exploration_summary$Biological.Plausibility <- c(
    "strongly_plausible", "mixed_context_dependent"
  )

  figure <- ai_figure_fun(report, render = FALSE)

  expect_identical(
    figure$pathways$Biological.Plausibility,
    c("strongly_plausible", "mixed_context_dependent")
  )
  expect_match(
    figure$pathways$Plausibility.Display[1L],
    "Strongly plausible", fixed = TRUE
  )
  expect_match(
    figure$caption,
    "not a causal grade", fixed = TRUE
  )

  report$exploration_summary$Biological.Plausibility[1L] <- "certain"
  expect_error(
    ai_figure_fun(report, render = FALSE),
    "Biological.Plausibility"
  )
})

test_that("channel counts are non-exclusive and scope order is fixed", {
  figure <- ai_figure_fun(ai_figure_fixture(), render = FALSE)

  expect_true(all(figure$channel_counts$Functional == 1L))
  expect_true(all(figure$channel_counts$Omics == 1L))
  expect_true(all(figure$channel_counts$Patient == 1L))
  expect_true(all(figure$channel_counts$Review == 1L))
  expect_true(all(figure$channel_counts$Unique == 1L))
  # Findings deliberately arrive as component, unclear, pathway.  Display order
  # is semantic and deterministic rather than dependent on provider row order.
  expect_true(all(figure$channel_scopes[, "Functional"] == "P/C/?"))
  expect_true(all(figure$channel_scopes[, "Omics"] == "C"))
  expect_true(all(figure$channel_scopes[, "Patient"] == "P"))

  shuffled <- ai_figure_fixture()
  shuffled$exploration_findings <- shuffled$exploration_findings[
    rev(seq_len(nrow(shuffled$exploration_findings))), , drop = FALSE
  ]
  again <- ai_figure_fun(shuffled, render = FALSE)
  expect_identical(again$channel_scopes, figure$channel_scopes)
  expect_identical(again$channel_counts, figure$channel_counts)
})

test_that("synthetic evidence receives a generic, dataset-independent warning", {
  report <- ai_figure_fixture(1L)
  report$evidence$provenance <- list(
    source_kind = "synthetic_fixture", database = "synthetic_not_pubmed"
  )
  report$interpretation$provider$name <- "arbitrary local adapter"

  figure <- ai_figure_fun(report, render = FALSE, text_mode = "full")

  expect_true(figure$evidence_source$synthetic_fixture)
  expect_match(figure$caption, "SYNTHETIC SOFTWARE FIXTURE ONLY", fixed = TRUE)
  expect_match(
    figure$pathways$Summary.Display,
    "Synthetic AI-evidence fixture", fixed = TRUE
  )
  expect_match(
    figure$pathways$Key.Evidence.Display,
    "Generated record IDs", fixed = TRUE
  )
  expect_false(grepl(
    "TCGA-PRAD|463 cases|99 BCR events|1,260-gene",
    figure$caption, perl = TRUE
  ))
})

test_that("quote quarantine and disease exclusion remain distinct", {
  quarantined <- ai_figure_fixture(1L)
  row <- quarantined$exploration_findings$Finding.ID == "review"
  quarantined$exploration_findings$Quote.Verified[row] <- FALSE
  quarantined$exploration_findings$Eligible.For.Summary[row] <- FALSE
  quarantined$exploration_summary$Review.PMIDs <- ""
  quarantined$exploration_summary$Key.PMIDs <-
    quarantined$exploration_summary$Functional.PMIDs
  quarantined$exploration_summary$Quote.Grounding <- "partial"
  quarantined$exploration_summary$Evidence.Eligibility <- "partial"
  figure <- ai_figure_fun(quarantined, render = FALSE)
  expect_identical(figure$grounding_counts$Quote.Quarantined, 1L)
  expect_identical(figure$grounding_counts$Disease.Excluded, 0L)

  excluded <- ai_figure_fixture(1L)
  row <- excluded$exploration_findings$Finding.ID == "review"
  excluded$exploration_findings$Eligible.For.Summary[row] <- FALSE
  excluded$exploration_summary$Review.PMIDs <- ""
  excluded$exploration_summary$Key.PMIDs <-
    excluded$exploration_summary$Functional.PMIDs
  excluded$exploration_summary$Evidence.Eligibility <- "partial"
  figure <- ai_figure_fun(excluded, render = FALSE)
  expect_identical(figure$grounding_counts$Quote.Quarantined, 0L)
  expect_identical(figure$grounding_counts$Disease.Excluded, 1L)
})

test_that("compact ASCII mode does not require a CJK font", {
  report <- ai_figure_fixture(1L)
  report$exploration_summary$Biological.Summary <-
    "\u4e2d\u6587\u751f\u7269\u5b66\u6458\u8981\u4ec5\u4fdd\u7559\u5728\u62a5\u544a\u4e2d\u3002"
  report$exploration_summary$Context.Dependence <- "\u4e2d\u6587\u60c5\u5883\u8bf4\u660e\u3002"

  compact <- ai_figure_fun(report, render = FALSE)
  full <- ai_figure_fun(report, render = FALSE, text_mode = "full")

  expect_false(compact$font$non_ascii)
  expect_identical(compact$font$selected, "sans")
  expect_true(all(compact$pathways$Context.Display == ""))
  expect_true(full$font$non_ascii)
  expect_false(full$font$coverage_checked)
  expect_error(full$draw(), "Font coverage")
})

test_that("invalid report contracts fail before rendering", {
  report <- ai_figure_fixture()
  report$purpose <- "prognostic_validation"
  expect_error(ai_figure_fun(report, render = FALSE), "purpose")

  report <- ai_figure_fixture()
  report$schema_version <- "3.9.0"
  expect_error(ai_figure_fun(report, render = FALSE), "schema")

  report <- ai_figure_fixture()
  report$top_pathways$Difference[1L] <- 1
  report$exploration_summary$Difference[1L] <- 1
  expect_error(ai_figure_fun(report, render = FALSE), "sign of Difference")

  report <- ai_figure_fixture()
  report$top_pathways$Direction.Rank <- rev(report$top_pathways$Direction.Rank)
  report$exploration_summary$Direction.Rank <-
    report$top_pathways$Direction.Rank
  expect_error(ai_figure_fun(report, render = FALSE), "unique increasing")

  report <- ai_figure_fixture()
  report$exploration_summary$Unique.PMID.Count[1L] <- 2L
  expect_error(ai_figure_fun(report, render = FALSE), "unique eligible")

  report <- ai_figure_fixture()
  report$evidence$records <- report$evidence$records[-1L, , drop = FALSE]
  expect_error(ai_figure_fun(report, render = FALSE), "absent from retrieved")

  report <- ai_figure_fixture()
  report$exploration_findings$Quote.Verified[1L] <- FALSE
  expect_error(ai_figure_fun(report, render = FALSE), "not verified")

  expect_error(ai_figure_fun(ai_figure_fixture()), "Supply a new")
})

test_that("PDF and PNG outputs are staged and published without leftovers", {
  skip_if_not(capabilities("png"))
  output_dir <- tempfile("intermeta-ai-figure-")
  expect_true(dir.create(output_dir))
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  pdf <- file.path(output_dir, "main.pdf")
  png <- file.path(output_dir, "main.png")

  figure <- ai_figure_fun(
    ai_figure_fixture(1L), pdf_file = pdf, png_file = png
  )

  expect_identical(figure$output$pdf, normalizePath(pdf, winslash = "/"))
  expect_identical(figure$output$png, normalizePath(png, winslash = "/"))
  expect_gt(file.info(pdf)$size, 1000)
  expect_gt(file.info(png)$size, 1000)
  expect_identical(readBin(pdf, what = "raw", n = 4L), charToRaw("%PDF"))
  expect_identical(
    readBin(png, what = "raw", n = 8L),
    as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
  )
  leftovers <- list.files(
    output_dir, all.files = TRUE, no.. = TRUE,
    pattern = "^\\.(main[.](pdf|png))[.]"
  )
  expect_length(leftovers, 0L)
  expect_error(
    ai_figure_fun(ai_figure_fixture(1L), pdf_file = pdf),
    "Refusing to overwrite"
  )
})

test_that("all output targets are checked before a graphics device is opened", {
  output_dir <- tempfile("intermeta-ai-preflight-")
  expect_true(dir.create(output_dir))
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  valid_pdf <- file.path(output_dir, "must-not-appear.pdf")
  missing_png <- file.path(output_dir, "missing", "main.png")
  before_device <- grDevices::dev.cur()

  expect_error(
    ai_figure_fun(
      ai_figure_fixture(1L), pdf_file = valid_pdf, png_file = missing_png
    ),
    "Output directory does not exist"
  )
  expect_false(file.exists(valid_pdf))
  expect_identical(grDevices::dev.cur(), before_device)

  missing_font_pdf <- file.path(output_dir, "missing-font.pdf")
  expect_error(
    ai_figure_fun(
      ai_figure_fixture(1L), pdf_file = missing_font_pdf,
      font_family = "__InterMeta_font_that_does_not_exist__"
    ),
    "Cairo|systemfonts|not installed"
  )
  expect_false(file.exists(missing_font_pdf))
})
