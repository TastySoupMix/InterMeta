test_that("dot plot defaults to first ten qualifying interactions in Index order", {
  data("InterMetaPRADEnrich", package = "InterMeta", envir = environment())
  plot <- InterMetaDotplot(
    InterMetaPRADEnrich,
    condition = "Permutation.FDR",
    cutoff = 0.2
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$data), 10L)
  expect_true(all(grepl("->", as.character(plot$data$Display.Interaction),
                        fixed = TRUE)))
  eligible <- which(
    is.finite(InterMetaPRADEnrich$Permutation.FDR) &
      InterMetaPRADEnrich$Permutation.FDR < 0.2
  )
  eligible <- eligible[order(InterMetaPRADEnrich$Index[eligible])]
  expect_equal(
    as.character(plot$data$Index),
    as.character(InterMetaPRADEnrich$Index[utils::head(eligible, 10)])
  )
})

test_that("score bar selects five pathways in each prognostic direction", {
  data("InterMetaPRADScore", package = "InterMeta", envir = environment())
  plot <- InterMetaPathwayScoreBar(InterMetaPRADScore)

  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$data), 10L)
  expect_equal(sum(plot$data$Difference > 0), 5L)
  expect_equal(sum(plot$data$Difference < 0), 5L)
  expect_true(all(diff(plot$data$Difference) <= 0))
})

test_that("stability score bar uses recurrent direction-consistent mean and SD", {
  summary <- data.frame(
    Pathway = c("F1", "F2", "F3", "U1", "U2", "U3", "Mixed", "Missing"),
    Selection.Count = c(3, 3, 3, 3, 3, 3, 3, 2),
    Selection.Frequency = c(rep(1, 7), 2 / 3),
    Dominant.Direction = c("F", "F", "F", "U", "U", "U", "Mixed", "F"),
    Direction.Consistency = c(rep(1, 6), 2 / 3, 1),
    Mean.Difference = c(5, 3, 1, -2, -4, -6, 0.5, 7),
    SD.Difference = c(0.5, 0.3, 0.1, 0.2, 0.4, 0.6, 1, 0.7),
    stringsAsFactors = FALSE
  )
  plot <- InterMetaStabilityScoreBar(summary, top = 2)
  plotdata <- attr(plot, "plot_data")

  expect_s3_class(plot, "ggplot")
  expect_equal(plotdata$Pathway, c("F1", "F2", "U2", "U3"))
  expect_equal(plotdata$Error, c(0.5, 0.3, 0.4, 0.6))
  expect_equal(plotdata$Lower, plotdata$Mean.Difference - plotdata$Error)
  expect_equal(plotdata$Upper, plotdata$Mean.Difference + plotdata$Error)
  expect_equal(plotdata$Direction, c("Favorable", "Favorable",
                                     "Unfavorable", "Unfavorable"))

  one_run <- summary[1, , drop = FALSE]
  one_run$Selection.Count <- 1L
  one_run$Selection.Frequency <- 1
  one_run$SD.Difference <- NA_real_
  expect_error(
    InterMetaStabilityScoreBar(one_run, error = "sd"),
    "At least two successful resample scores"
  )
  expect_s3_class(
    InterMetaStabilityScoreBar(one_run, error = "none"),
    "ggplot"
  )
})
