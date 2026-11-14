# Portions of this file were adapted from TimiGP 1.3.0:
# https://github.com/CSkylarL/TimiGP
# The upstream source identifies Chenyang Skylar Li as author.
# Modified for InterMeta; changes through 2026-11-03.
# SPDX-License-Identifier: GPL-3.0-only
# See inst/NOTICE and COPYING for provenance and license details.

## Internal Cox implementation inherited from TimiGP and hardened for
## long-running InterMeta analyses.

TimiCOX <- function(mps = NULL,
                    info = NULL,
                    p.adj = "BH",
                    parallel = FALSE,
                    core = 1,
                    keep_mps = TRUE,
                    analysis_profile = "TimiGP-legacy"){
  analysis_profile <- .InterMetaAnalysisProfile(analysis_profile)
  if (is.null(mps) || !(is.matrix(mps) || is.data.frame(mps))){
    stop('"mps" must be a matrix or data frame.')
  }
  if (is.null(info) || !is.data.frame(info) || ncol(info) < 2){
    stop('"info" must contain event and time columns.')
  }
  if (is.null(row.names(mps)) || is.null(colnames(mps))){
    stop('"mps" must have gene-pair row names and sample column names.')
  }
  if (!p.adj %in% stats::p.adjust.methods){
    stop('"p.adj" must be one of stats::p.adjust.methods.')
  }
  if (!is.logical(parallel) || length(parallel) != 1 || is.na(parallel)){
    stop('"parallel" must be TRUE or FALSE.')
  }
  if (!is.logical(keep_mps) || length(keep_mps) != 1 || is.na(keep_mps)){
    stop('"keep_mps" must be TRUE or FALSE.')
  }
  core <- as.integer(core)
  if (is.na(core) || core < 1){
    stop('"core" must be a positive integer.')
  }

  common_samples <- intersect(colnames(mps), row.names(info))
  if (length(common_samples) < 5){
    stop("Fewer than 5 common samples were found between mps and info.")
  }
  if (!identical(colnames(mps), common_samples)){
    mps <- mps[, common_samples, drop = FALSE]
  }
  if (!identical(row.names(info), common_samples)){
    info <- info[common_samples, , drop = FALSE]
  }
  event <- as.numeric(info[[1]])
  time <- as.numeric(info[[2]])
  if (any(!is.finite(event)) || any(!is.finite(time)) ||
      any(!event %in% c(0, 1)) || any(time <= 0)){
    stop("Clinical event/time values are invalid.")
  }
  if (!(is.logical(mps) || is.numeric(mps))){
    stop('"mps" must contain logical or numeric values.')
  }

  logical_mps <- is.logical(mps)
  fit_one <- function(tag){
    numeric_tag <- as.numeric(tag)
    if (any(!is.finite(numeric_tag)) || length(unique(numeric_tag)) < 2){
      return(c(HR = NA_real_, PV = NA_real_))
    }
    model_tag <- if (logical_mps) {
      as.factor(tag)
    } else {
      numeric_tag
    }
    tryCatch({
      fit <- survival::coxph(
        survival::Surv(time, event) ~ model_tag,
        data = data.frame(time = time, event = event,
                          model_tag = model_tag)
      )
      fit_summary <- summary(fit)
      c(
        HR = unname(fit_summary$conf.int[1, "exp(coef)"]),
        PV = unname(fit_summary$coefficients[1, "Pr(>|z|)"])
      )
    }, error = function(error){
      c(HR = NA_real_, PV = NA_real_)
    })
  }

  message("Calculating univariate Cox models for ", nrow(mps),
          " gene pairs.")
  if (nrow(mps) == 0){
    stop("No gene pairs are available for Cox analysis.")
  }
  if (parallel && core > 1 && nrow(mps) > 1){
    worker_count <- min(core, nrow(mps))
    row_chunks <- split(
      seq_len(nrow(mps)),
      cut(seq_len(nrow(mps)), breaks = worker_count, labels = FALSE)
    )
    if (.Platform$OS.type != "windows"){
      result_blocks <- parallel::mclapply(
        row_chunks,
        function(index){
          t(vapply(index, function(i) fit_one(mps[i, ]), numeric(2)))
        },
        mc.cores = worker_count,
        mc.preschedule = TRUE
      )
      result_rows <- do.call(rbind, result_blocks)
    } else if (as.numeric(utils::object.size(mps)) * 2 > 2 * 1024^3){
      warning("Parallel Cox was disabled on Windows because chunk transfer ",
              "would require more than approximately 2 GB of extra memory.")
      result_rows <- t(vapply(
        seq_len(nrow(mps)),
        function(i) fit_one(mps[i, ]),
        numeric(2)
      ))
    } else {
      mps_chunks <- lapply(row_chunks, function(index){
        mps[index, , drop = FALSE]
      })
      cluster <- parallel::makeCluster(worker_count)
      on.exit({
        parallel::stopCluster(cluster)
        foreach::registerDoSEQ()
      }, add = TRUE)
      worker_libraries <- .libPaths()
      parallel::clusterCall(cluster, function(paths){
        .libPaths(paths)
        NULL
      }, worker_libraries)
      doParallel::registerDoParallel(cluster)
      result_rows <- foreach::foreach(
        chunk = mps_chunks,
        .combine = rbind,
        .inorder = TRUE,
        .packages = "survival"
      ) %dopar% {
        t(vapply(
          seq_len(nrow(chunk)),
          function(i) fit_one(chunk[i, ]),
          numeric(2)
        ))
      }
      foreach::registerDoSEQ()
      rm(mps_chunks)
    }
  } else {
    result_rows <- t(vapply(
      seq_len(nrow(mps)),
      function(i) fit_one(mps[i, ]),
      numeric(2)
    ))
  }

  cox_res <- data.frame(
    HR = result_rows[, "HR"],
    PV = result_rows[, "PV"],
    stringsAsFactors = FALSE,
    row.names = row.names(mps)
  )
  cox_res$QV <- stats::p.adjust(cox_res$PV, method = p.adj)
  failed <- sum(!is.finite(cox_res$HR) | !is.finite(cox_res$PV))
  if (failed > 0){
    warning(failed, " gene-pair Cox models failed and were recorded as NA.")
  }

  reverse <- which(is.finite(cox_res$HR) & cox_res$HR > 1)
  message("The direction of ", length(reverse), " gene pairs was reversed.")
  if (length(reverse) > 0){
    parsed <- .InterMetaParseGenePairs(row.names(cox_res)[reverse])
    reversed_name <- paste(parsed$Target.Gene, parsed$Source.Gene, sep = "_")
    cox_res$HR[reverse] <- 1 / cox_res$HR[reverse]
    row.names(cox_res)[reverse] <- reversed_name
    if (keep_mps){
      row.names(mps)[reverse] <- reversed_name
      if (is.logical(mps)){
        mps[reverse, ] <- !mps[reverse, , drop = FALSE]
      } else {
        mps[reverse, ] <- -mps[reverse, , drop = FALSE]
      }
    }
  }

  order_index <- order(cox_res$PV, na.last = TRUE)
  list(
    mps = if (keep_mps) mps else NULL,
    cox_res = cox_res[order_index, , drop = FALSE]
  )
}
