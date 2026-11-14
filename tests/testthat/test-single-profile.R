test_that("every statistical entry point defaults to the sole legacy mode", {
  entry_points <- c("InterMetaPreprocess", "InterMetaGenePair", "InterMetaCox",
                    "InterMetaPermFDR", "InterMetaPathwayScore", "InterMetaStability",
                    "TimiCOX")
  for (name in entry_points) {
    fn <- get(name, envir = asNamespace("InterMeta"))
    expect_identical(formals(fn)$analysis_profile, "TimiGP-legacy", info = name)
    for (bad in list("InterMeta", "TimiGP", "unknown", NA_character_, NULL,
                     c("TimiGP-legacy", "InterMeta"))) {
      expect_error(do.call(fn, list(analysis_profile = bad)),
                   "Only analysis_profile", info = name)
    }
  }
})

test_that("the bundled PRAD score is exactly the default legacy view", {
  data("InterMetaPRADEnrich", package = "InterMeta", envir = environment())
  data("InterMetaPRADScore", package = "InterMeta", envir = environment())
  score <- InterMetaPathwayScore(InterMetaPRADEnrich)
  expect_equal(score, InterMetaPRADScore)
  expect_equal(nrow(score), 50L)
  expect_false(any(score$Difference == 0))
  expect_equal(score$Rank, score$Signed.Rank)
})
