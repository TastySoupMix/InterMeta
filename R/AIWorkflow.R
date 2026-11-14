# Recoverable orchestration for the PubMed -> screening -> synthesis workflow.
#
# The workflow deliberately stores only non-secret configuration.  Provider
# credentials continue to be resolved by the provider at call time.

.InterMetaAIWorkflowSchemaVersion <- "1.1.0"
.InterMetaAIWorkflowStages <- c(
  "context", "queries", "pubmed", "screening", "verified_evidence",
  "synthesis", "report", "figure"
)
.InterMetaAIWorkflowSuccessfulStatuses <- c(
  "complete", "complete_with_quarantine", "complete_with_truncation"
)
.InterMetaAIWorkflowReusableStatuses <- c(
  .InterMetaAIWorkflowSuccessfulStatuses, "skipped"
)

.InterMetaAIWorkflowNow <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

.InterMetaAIWorkflowScalar <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !nzchar(trimws(x))) {
    stop(name, " must be one non-empty character value.", call. = FALSE)
  }
  trimws(x)
}

.InterMetaAIWorkflowFlag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.InterMetaAIWorkflowInteger <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < minimum || x != floor(x) || x > .Machine$integer.max) {
    stop(name, " must be one integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.InterMetaAIWorkflowHashObject <- function(x) {
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

.InterMetaAIWorkflowHashFile <- function(path) {
  if (!file.exists(path) || isTRUE(file.info(path)$isdir)) return(NA_character_)
  digest::digest(path, algo = "sha256", file = TRUE)
}

.InterMetaAIWorkflowAtomicRDS <- function(x, path) {
  parent <- dirname(path)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create workflow directory: ", parent, call. = FALSE)
  }
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = parent)
  on.exit(if (file.exists(temporary)) unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(x, temporary, version = 3L)
  if (!file.copy(temporary, path, overwrite = TRUE, copy.mode = TRUE,
                 copy.date = FALSE)) {
    stop("Could not publish workflow file: ", path, call. = FALSE)
  }
  invisible(path)
}

.InterMetaAIWorkflowAtomicTSV <- function(x, path) {
  if (!is.data.frame(x)) return(invisible(NULL))
  parent <- dirname(path)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create workflow directory: ", parent, call. = FALSE)
  }
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = parent)
  on.exit(if (file.exists(temporary)) unlink(temporary, force = TRUE), add = TRUE)
  utils::write.table(
    x, file = temporary, sep = "\t", quote = TRUE, row.names = FALSE,
    col.names = TRUE, na = "NA", fileEncoding = "UTF-8"
  )
  if (!file.copy(temporary, path, overwrite = TRUE, copy.mode = TRUE,
                 copy.date = FALSE)) {
    stop("Could not publish workflow table: ", path, call. = FALSE)
  }
  invisible(path)
}

.InterMetaAIWorkflowPlainData <- function(x, depth = 0L) {
  if (depth > 64L) return(FALSE)
  if (is.null(x)) return(TRUE)
  if (is.function(x) || is.environment(x) || is.language(x) ||
      is.expression(x) || isS4(x) ||
      typeof(x) %in% c("externalptr", "weakref")) return(FALSE)
  if (!is.atomic(x) && !is.list(x)) return(FALSE)
  if (is.numeric(x) && any(!is.finite(x), na.rm = TRUE)) return(FALSE)
  if (is.list(x) && !all(vapply(
    x, .InterMetaAIWorkflowPlainData, logical(1L), depth = depth + 1L
  ))) return(FALSE)
  TRUE
}

.InterMetaAIWorkflowContainsSecret <- function(x, secrets) {
  secrets <- unique(as.character(secrets))
  secrets <- secrets[!is.na(secrets) & nzchar(secrets)]
  if (!length(secrets)) return(FALSE)
  raw <- serialize(x, NULL, ascii = TRUE, version = 3L)
  text <- rawToChar(raw)
  any(vapply(secrets, function(secret) grepl(secret, text, fixed = TRUE),
             logical(1L)))
}

.InterMetaAIWorkflowPaths <- function(run_dir) {
  list(
    config = file.path(run_dir, "workflow_config.rds"),
    state = file.path(run_dir, "workflow_state.rds"),
    state_tsv = file.path(run_dir, "workflow_state.tsv"),
    supplied_queries = file.path(run_dir, "inputs", "supplied_queries.rds"),
    supplied_evidence = file.path(run_dir, "inputs", "supplied_evidence.rds"),
    context = file.path(run_dir, "00_context", "context.rds"),
    queries = file.path(run_dir, "01_queries", "queries.rds"),
    evidence = file.path(run_dir, "02_pubmed", "pubmed_evidence.rds"),
    screening = file.path(run_dir, "03_screening", "screening.rds"),
    screened = file.path(run_dir, "04_verified", "evidence_bundle.rds"),
    synthesis = file.path(run_dir, "05_synthesis", "synthesis.rds"),
    report = file.path(run_dir, "06_report", "report.rds"),
    figure = file.path(run_dir, "07_figure", "figure.rds")
  )
}

.InterMetaAIWorkflowNewState <- function() {
  stages <- data.frame(
    Stage = .InterMetaAIWorkflowStages,
    Status = rep("pending", length(.InterMetaAIWorkflowStages)),
    Attempt = rep(0L, length(.InterMetaAIWorkflowStages)),
    Input.SHA256 = rep("", length(.InterMetaAIWorkflowStages)),
    Output.SHA256 = rep("", length(.InterMetaAIWorkflowStages)),
    Artifact = rep("", length(.InterMetaAIWorkflowStages)),
    Started.UTC = rep("", length(.InterMetaAIWorkflowStages)),
    Completed.UTC = rep("", length(.InterMetaAIWorkflowStages)),
    Message = rep("", length(.InterMetaAIWorkflowStages)),
    stringsAsFactors = FALSE
  )
  list(
    schema_version = .InterMetaAIWorkflowSchemaVersion,
    status = "initialized", created_at = .InterMetaAIWorkflowNow(),
    updated_at = .InterMetaAIWorkflowNow(), stages = stages
  )
}

.InterMetaAIWorkflowSaveState <- function(state, run_dir) {
  state$updated_at <- .InterMetaAIWorkflowNow()
  paths <- .InterMetaAIWorkflowPaths(run_dir)
  .InterMetaAIWorkflowAtomicRDS(state, paths$state)
  .InterMetaAIWorkflowAtomicTSV(state$stages, paths$state_tsv)
  state
}

.InterMetaAIWorkflowSetStage <- function(state, run_dir, stage, status,
                                         input_sha256 = NULL,
                                         output_sha256 = NULL,
                                         artifact = NULL, message = NULL,
                                         increment = FALSE) {
  index <- match(stage, state$stages$Stage)
  if (is.na(index)) stop("Unknown workflow stage: ", stage, call. = FALSE)
  if (increment) {
    state$stages$Attempt[index] <- state$stages$Attempt[index] + 1L
    state$stages$Started.UTC[index] <- .InterMetaAIWorkflowNow()
    state$stages$Completed.UTC[index] <- ""
  }
  state$stages$Status[index] <- status
  if (!is.null(input_sha256)) state$stages$Input.SHA256[index] <- input_sha256
  if (!is.null(output_sha256)) state$stages$Output.SHA256[index] <- output_sha256
  if (!is.null(artifact)) state$stages$Artifact[index] <- artifact
  if (!is.null(message)) state$stages$Message[index] <- message
  if (status %in% .InterMetaAIWorkflowReusableStatuses) {
    state$stages$Completed.UTC[index] <- .InterMetaAIWorkflowNow()
  }
  .InterMetaAIWorkflowSaveState(state, run_dir)
}

.InterMetaAIWorkflowRelative <- function(path, run_dir) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root, "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
}

.InterMetaAIWorkflowAbsolutePath <- function(path) {
  is.character(path) && length(path) == 1L && !is.na(path) &&
    grepl("^(?:[A-Za-z]:[/\\\\]|[/\\\\])", path, perl = TRUE)
}

.InterMetaAIWorkflowRelativeCachePath <- function(path, name) {
  path <- .InterMetaAIWorkflowScalar(path, name)
  path <- gsub("\\\\", "/", path)
  if (.InterMetaAIWorkflowAbsolutePath(path) ||
      grepl("^[A-Za-z]:", path, perl = TRUE) || startsWith(path, "~")) {
    stop(name, " must be a relative path below the run directory parent.",
         call. = FALSE)
  }
  pieces <- strsplit(path, "/", fixed = TRUE)[[1L]]
  pieces <- pieces[nzchar(pieces) & pieces != "."]
  if (!length(pieces) || any(pieces == "..")) {
    stop(name, " must be a safe relative path without '..'.", call. = FALSE)
  }
  paste(pieces, collapse = "/")
}

.InterMetaAIWorkflowCacheSpec <- function(path, portable_default) {
  if (portable_default) {
    return(list(cache_dir = "pubmed_cache", cache_anchor = "run_parent"))
  }
  path <- .InterMetaAIWorkflowScalar(path, "pubmed_cache_dir")
  if (.InterMetaAIWorkflowAbsolutePath(path)) {
    return(list(
      cache_dir = normalizePath(path, winslash = "/", mustWork = FALSE),
      cache_anchor = "absolute"
    ))
  }
  list(
    cache_dir = .InterMetaAIWorkflowRelativeCachePath(
      path, "pubmed_cache_dir"
    ),
    cache_anchor = "run_parent"
  )
}

.InterMetaAIWorkflowOperationalCachePath <- function(path, run_dir, name) {
  path <- .InterMetaAIWorkflowScalar(path, name)
  if (.InterMetaAIWorkflowAbsolutePath(path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  relative <- .InterMetaAIWorkflowRelativeCachePath(path, name)
  normalizePath(file.path(dirname(run_dir), relative), winslash = "/",
                mustWork = FALSE)
}

.InterMetaAIWorkflowResolveCacheDir <- function(config, run_dir,
                                                 override = NULL) {
  if (!is.null(override)) {
    return(.InterMetaAIWorkflowOperationalCachePath(
      override, run_dir, "pubmed_cache_dir"
    ))
  }
  settings <- config$pubmed_settings
  path <- .InterMetaAIWorkflowScalar(
    settings$cache_dir, "saved pubmed_settings$cache_dir"
  )
  anchor <- settings$cache_anchor
  if (identical(anchor, "run_parent")) {
    relative <- .InterMetaAIWorkflowRelativeCachePath(
      path, "saved pubmed_settings$cache_dir"
    )
    return(normalizePath(
      file.path(dirname(run_dir), relative), winslash = "/",
      mustWork = FALSE
    ))
  }
  if (identical(anchor, "absolute")) {
    if (!.InterMetaAIWorkflowAbsolutePath(path)) {
      stop("The saved absolute PubMed cache path is malformed.",
           call. = FALSE)
    }
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  if (!is.null(anchor)) {
    stop("The saved PubMed cache anchor is unsupported; pass pubmed_cache_dir to InterMetaAIResume().",
         call. = FALSE)
  }

  # InterMeta 1.1.1 stored only an absolute cache path.  Preserve that path
  # when it is still available.  A missing historical default is safely
  # relocatable because its basename was fixed; custom paths require an
  # explicit operational override rather than a guess.
  if (.InterMetaAIWorkflowAbsolutePath(path) && dir.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  legacy_name <- basename(gsub("\\\\", "/", path))
  if (identical(tolower(legacy_name), "pubmed_cache")) {
    return(normalizePath(
      file.path(dirname(run_dir), "pubmed_cache"), winslash = "/",
      mustWork = FALSE
    ))
  }
  stop(paste(
    "The saved legacy PubMed cache path is unavailable and cannot be",
    "mapped reliably. Pass pubmed_cache_dir to InterMetaAIResume()."
  ), call. = FALSE)
}

.InterMetaAIWorkflowPathIdentity <- function(path) {
  path <- gsub("\\\\", "/", path)
  if (grepl("^[A-Za-z]:/", path)) {
    root <- tolower(substr(path, 1L, 2L))
    remainder <- substring(path, 4L)
    insensitive <- TRUE
  } else if (startsWith(path, "//")) {
    pieces <- strsplit(sub("^//", "", path), "/", fixed = FALSE)[[1L]]
    if (length(pieces) < 2L) return(NULL)
    root <- paste0("//", tolower(pieces[[1L]]), "/",
                   tolower(pieces[[2L]]))
    remainder <- paste(pieces[-c(1L, 2L)], collapse = "/")
    insensitive <- TRUE
  } else if (startsWith(path, "/")) {
    root <- "/"
    remainder <- substring(path, 2L)
    insensitive <- FALSE
  } else {
    return(NULL)
  }
  components <- strsplit(remainder, "/", fixed = TRUE)[[1L]]
  components <- components[nzchar(components)]
  list(root = root, components = components,
       compare = if (insensitive) tolower(components) else components)
}

.InterMetaAIWorkflowRelativeBetween <- function(target, base) {
  target <- normalizePath(target, winslash = "/", mustWork = FALSE)
  base <- normalizePath(base, winslash = "/", mustWork = FALSE)
  target_id <- .InterMetaAIWorkflowPathIdentity(target)
  base_id <- .InterMetaAIWorkflowPathIdentity(base)
  if (is.null(target_id) || is.null(base_id) ||
      !identical(target_id$root, base_id$root)) return(NULL)
  limit <- min(length(target_id$compare), length(base_id$compare))
  common <- 0L
  if (limit > 0L) {
    for (index in seq_len(limit)) {
      if (!identical(target_id$compare[[index]],
                     base_id$compare[[index]])) break
      common <- index
    }
  }
  upward <- rep("..", length(base_id$components) - common)
  downward <- if (common < length(target_id$components)) {
    target_id$components[seq.int(common + 1L, length(target_id$components))]
  } else character()
  relative <- c(upward, downward)
  if (!length(relative)) "." else paste(relative, collapse = "/")
}

.InterMetaAIWorkflowReplaySourceForStorage <- function(run_dir, replay_dir) {
  relative <- .InterMetaAIWorkflowRelativeBetween(run_dir, replay_dir)
  if (is.null(relative)) {
    return(list(
      schema_version = "1", kind = "content_identity", anchor = "none",
      path = NULL, name = basename(run_dir)
    ))
  }
  list(
    schema_version = "1", kind = "relative_path", anchor = "replay_dir",
    path = relative, name = basename(run_dir)
  )
}

.InterMetaAIWorkflowReplaySourceForSession <- function(source_run,
                                                        replay_dir) {
  # Compatibility with replay.rds files written through InterMeta 1.1.1,
  # whose source_run field was a character path (normally absolute).
  if (is.character(source_run) && length(source_run) == 1L &&
      !is.na(source_run) && nzchar(source_run)) {
    candidate <- if (.InterMetaAIWorkflowAbsolutePath(source_run)) {
      source_run
    } else file.path(replay_dir, source_run)
    return(normalizePath(candidate, winslash = "/",
                         mustWork = dir.exists(candidate)))
  }
  if (!is.list(source_run) ||
      !identical(source_run$schema_version, "1") ||
      !identical(source_run$kind, "relative_path") ||
      !identical(source_run$anchor, "replay_dir") ||
      !is.character(source_run$path) || length(source_run$path) != 1L ||
      is.na(source_run$path) || !nzchar(source_run$path) ||
      .InterMetaAIWorkflowAbsolutePath(source_run$path) ||
      grepl("^[A-Za-z]:", source_run$path, perl = TRUE)) {
    return(NA_character_)
  }
  candidate <- file.path(replay_dir, source_run$path)
  normalizePath(candidate, winslash = "/", mustWork = dir.exists(candidate))
}

.InterMetaAIWorkflowSynthesisForStageArtifact <- function(synthesis) {
  if (!is.list(synthesis) ||
      !identical(synthesis$checkpoint_path_base, "outdir")) {
    return(synthesis)
  }
  checkpoints <- synthesis$checkpoints
  if (is.null(checkpoints)) checkpoints <- character()
  if (!is.character(checkpoints) || anyNA(checkpoints)) {
    stop("Synthesis checkpoint references are malformed.", call. = FALSE)
  }
  normalize_reference <- function(path) {
    path <- gsub("\\\\", "/", path)
    pieces <- strsplit(path, "/", fixed = TRUE)[[1L]]
    pieces <- pieces[nzchar(pieces) & pieces != "."]
    if (!length(pieces) || .InterMetaAIWorkflowAbsolutePath(path) ||
        startsWith(path, "~") || any(pieces == "..")) {
      stop("Synthesis checkpoint references must be safe relative paths.",
           call. = FALSE)
    }
    paste(c("checkpoints", pieces), collapse = "/")
  }
  synthesis$checkpoints <- vapply(
    checkpoints, normalize_reference, character(1L), USE.NAMES = FALSE
  )
  # The workflow stores this aggregate one directory above the checkpoint
  # engine's own aggregate.  Make the exposed base explicit so a reader can
  # resolve each reference from dirname(05_synthesis/synthesis.rds).
  synthesis$checkpoint_path_base <- "synthesis_artifact_dir"
  synthesis
}

.InterMetaAIWorkflowAnchoredRelative <- function(path, anchors) {
  normalized <- gsub("\\\\", "/", path)
  pieces <- strsplit(normalized, "/", fixed = TRUE)[[1L]]
  positions <- which(pieces %in% anchors)
  if (!length(positions)) return(NULL)
  paste(pieces[positions[[length(positions)]]:length(pieces)], collapse = "/")
}

.InterMetaAIWorkflowMapFigurePaths <- function(figure, mapper) {
  visit <- function(value) {
    if (!is.list(value)) return(value)
    if ("output" %in% names(value) && is.list(value$output)) {
      for (format in intersect(c("pdf", "png"), names(value$output))) {
        path <- value$output[[format]]
        if (is.character(path) && length(path) == 1L && !is.na(path) &&
            nzchar(path)) {
          value$output[[format]] <- mapper(path)
        }
      }
    }
    nested <- which(vapply(value, is.list, logical(1L)))
    value_names <- names(value)
    if (!is.null(value_names)) {
      nested_names <- value_names[nested]
      nested <- nested[is.na(nested_names) | nested_names != "output"]
    }
    for (index in nested) value[[index]] <- visit(value[[index]])
    value
  }
  visit(figure)
}

.InterMetaAIWorkflowFigureForStorage <- function(
    figure, root, anchors, strict = TRUE) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  .InterMetaAIWorkflowMapFigurePaths(figure, function(path) {
    if (!.InterMetaAIWorkflowAbsolutePath(path)) {
      return(gsub("\\\\", "/", path))
    }
    relative <- .InterMetaAIWorkflowRelative(path, root)
    if (!.InterMetaAIWorkflowAbsolutePath(relative)) return(relative)
    anchored <- .InterMetaAIWorkflowAnchoredRelative(path, anchors)
    if (!is.null(anchored)) return(anchored)
    if (strict) {
      stop("A rendered figure output path is outside its artifact root.",
           call. = FALSE)
    }
    normalizePath(path, winslash = "/", mustWork = FALSE)
  })
}

.InterMetaAIWorkflowFigureForSession <- function(figure, root, anchors) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  .InterMetaAIWorkflowMapFigurePaths(figure, function(path) {
    if (!.InterMetaAIWorkflowAbsolutePath(path)) {
      return(normalizePath(file.path(root, path), winslash = "/",
                           mustWork = FALSE))
    }
    relative <- .InterMetaAIWorkflowRelative(path, root)
    if (!.InterMetaAIWorkflowAbsolutePath(relative)) {
      return(normalizePath(file.path(root, relative), winslash = "/",
                           mustWork = FALSE))
    }
    anchored <- .InterMetaAIWorkflowAnchoredRelative(path, anchors)
    if (!is.null(anchored)) {
      return(normalizePath(file.path(root, anchored), winslash = "/",
                           mustWork = FALSE))
    }
    normalizePath(path, winslash = "/", mustWork = FALSE)
  })
}

.InterMetaAIWorkflowComplete <- function(state, run_dir, stage, input_sha256,
                                         artifact,
                                         status = "complete",
                                         message = "") {
  output_sha256 <- .InterMetaAIWorkflowHashFile(artifact)
  if (is.na(output_sha256)) {
    stop("Workflow stage did not create its artifact: ", artifact,
         call. = FALSE)
  }
  marker <- list(
    schema_version = .InterMetaAIWorkflowSchemaVersion,
    stage = stage, status = status, input_sha256 = input_sha256,
    output_sha256 = output_sha256,
    artifact = .InterMetaAIWorkflowRelative(artifact, run_dir),
    completed_at = .InterMetaAIWorkflowNow()
  )
  .InterMetaAIWorkflowAtomicRDS(
    marker, file.path(dirname(artifact), "_SUCCESS.rds")
  )
  .InterMetaAIWorkflowSetStage(
    state, run_dir, stage, status, input_sha256, output_sha256,
    marker$artifact, message
  )
}

.InterMetaAIWorkflowStageValid <- function(state, run_dir, stage,
                                           input_sha256 = NULL) {
  index <- match(stage, state$stages$Stage)
  if (is.na(index) || !state$stages$Status[index] %in%
      .InterMetaAIWorkflowReusableStatuses) return(FALSE)
  artifact <- state$stages$Artifact[index]
  if (!nzchar(artifact)) return(FALSE)
  path <- if (grepl("^(?:[A-Za-z]:[/\\\\]|[/\\\\])", artifact, perl = TRUE)) {
    artifact
  } else file.path(run_dir, artifact)
  observed <- .InterMetaAIWorkflowHashFile(path)
  if (is.na(observed) || !identical(observed, state$stages$Output.SHA256[index])) {
    return(FALSE)
  }
  if (!is.null(input_sha256) &&
      !identical(input_sha256, state$stages$Input.SHA256[index])) return(FALSE)
  marker_path <- file.path(dirname(path), "_SUCCESS.rds")
  marker <- tryCatch(readRDS(marker_path), error = function(e) NULL)
  is.list(marker) && identical(marker$stage, stage) &&
    identical(marker$output_sha256, observed) &&
    (is.null(input_sha256) || identical(marker$input_sha256, input_sha256))
}

.InterMetaAIWorkflowDependency <- function(name) {
  hooks <- getOption("InterMeta.ai_workflow_dependencies", NULL)
  if (is.list(hooks) && is.function(hooks[[name]])) return(hooks[[name]])
  environment <- environment(InterMetaAIWorkflow)
  if (!exists(name, envir = environment, mode = "function", inherits = TRUE)) {
    stop("Required InterMeta workflow dependency is unavailable: ", name,
         call. = FALSE)
  }
  get(name, envir = environment, mode = "function", inherits = TRUE)
}

.InterMetaAIWorkflowCall <- function(name, arguments) {
  fun <- .InterMetaAIWorkflowDependency(name)
  formal_names <- names(formals(fun))
  if (!"..." %in% formal_names) {
    arguments <- arguments[names(arguments) %in% formal_names]
  }
  do.call(fun, arguments)
}

.InterMetaAIWorkflowProviderDescriptor <- function(provider, model) {
  if (is.null(provider)) return(list(kind = "deepseek", name = "DeepSeek", model = model))
  if (!inherits(provider, "InterMetaAIProvider") || !is.list(provider) ||
      !is.function(provider$generate) || !is.character(provider$name) ||
      length(provider$name) != 1L || !is.character(provider$model) ||
      length(provider$model) != 1L) {
    stop("provider must be NULL or an InterMetaAIProvider.", call. = FALSE)
  }
  list(kind = "supplied", name = provider$name, model = provider$model)
}

.InterMetaAIWorkflowResolveProvider <- function(config, provider) {
  expected <- config$provider
  if (is.null(provider)) {
    if (!identical(expected$kind, "deepseek")) {
      stop("This run used a supplied provider; pass the same provider when resuming model stages.",
           call. = FALSE)
    }
    settings <- config$model_settings
    return(.InterMetaAIWorkflowCall("InterMetaDeepSeekProvider", list(
      model = expected$model, api_key_env = settings$api_key_env,
      max_tokens = settings$max_tokens, timeout = settings$timeout,
      thinking = settings$thinking,
      reasoning_effort = settings$reasoning_effort
    )))
  }
  observed <- .InterMetaAIWorkflowProviderDescriptor(provider, expected$model)
  if (!identical(observed$name, expected$name) ||
      !identical(observed$model, expected$model)) {
    stop("The provider name/model does not match the frozen workflow configuration.",
         call. = FALSE)
  }
  provider
}

.InterMetaAIWorkflowCandidates <- function(evidence) {
  records <- evidence$records
  if (!is.data.frame(records) || !nrow(records) ||
      !all(c("pathway", "PMID") %in% names(records))) {
    return(data.frame(
      Candidate.ID = character(), pathway = character(), PMID = character(),
      stringsAsFactors = FALSE
    ))
  }
  keys <- paste(as.character(records$pathway), as.character(records$PMID),
                sep = "\r")
  first <- !duplicated(keys)
  core <- intersect(
    c("pathway", "direction", "PMID", "Title", "Abstract", "Journal", "Year",
      "DOI", "PublicationType", "MeSH"), names(records)
  )
  out <- records[first, core, drop = FALSE]
  unique_keys <- keys[first]
  collapse <- function(column) {
    if (!column %in% names(records)) return(rep("", length(unique_keys)))
    vapply(unique_keys, function(key) {
      values <- unique(as.character(records[[column]][keys == key]))
      values <- values[!is.na(values) & nzchar(values)]
      paste(values, collapse = ";")
    }, character(1L))
  }
  out$Query.IDs <- collapse("query_id")
  out$Layers <- collapse("layer")
  out$Term.Scopes <- collapse("term_scope")
  out$Candidate.ID <- sprintf("c%06d", seq_len(nrow(out)))
  out <- out[, c("Candidate.ID", setdiff(names(out), "Candidate.ID")), drop = FALSE]
  rownames(out) <- NULL
  out
}

.InterMetaAIWorkflowPubMedReady <- function(evidence) {
  searches <- evidence$searches
  is.data.frame(searches) && nrow(searches) > 0L &&
    all(c("query_id", "api_status") %in% names(searches)) &&
    !anyDuplicated(searches$query_id) && !anyNA(searches$api_status) &&
    all(searches$api_status == "ok")
}

.InterMetaAIWorkflowPubMedTruncated <- function(evidence) {
  searches <- evidence$searches
  is.data.frame(searches) && "truncated" %in% names(searches) &&
    any(searches$truncated %in% TRUE)
}

.InterMetaAIWorkflowQueriesMatch <- function(queries, evidence_queries) {
  if (!is.data.frame(queries) || !is.data.frame(evidence_queries) ||
      !all(c("query_id", "query") %in% names(queries)) ||
      !all(c("query_id", "query") %in% names(evidence_queries)) ||
      anyNA(queries$query_id) || anyNA(evidence_queries$query_id) ||
      anyDuplicated(queries$query_id) || anyDuplicated(evidence_queries$query_id) ||
      !setequal(as.character(queries$query_id),
                as.character(evidence_queries$query_id))) return(FALSE)
  columns <- intersect(
    c("query_id", "pathway", "direction", "layer", "term_scope",
      "search_terms", "query_version", "query", "purpose"),
    intersect(names(queries), names(evidence_queries))
  )
  normalize <- function(x) {
    x <- x[order(as.character(x$query_id), method = "radix"), columns,
           drop = FALSE]
    x[] <- lapply(x, as.character)
    rownames(x) <- NULL
    x
  }
  identical(normalize(queries), normalize(evidence_queries))
}

.InterMetaAIWorkflowCompleteFlag <- function(x) {
  if (!is.list(x)) return(FALSE)
  if (is.logical(x$complete) && length(x$complete) == 1L &&
      !is.na(x$complete)) return(x$complete)
  status <- x$status
  is.character(status) && length(status) == 1L && !is.na(status) &&
    status %in% c(.InterMetaAIWorkflowSuccessfulStatuses,
                  "completed", "success", "ok")
}

.InterMetaAIWorkflowScreeningNoCallReady <- function(
    run_dir, context, candidates, provider_descriptor, batch_size) {
  input <- tryCatch(
    .InterMetaAIScreeningInputValidation(
      candidates, as.character(context$pathways$Pathway)
    ),
    error = function(e) NULL
  )
  if (is.null(input)) return(FALSE)
  eligible <- candidates[input$eligible, , drop = FALSE]
  if (!nrow(eligible)) return(TRUE)
  pathway_order <- as.character(context$pathways$Pathway)
  expected <- character()
  for (pathway_index in seq_along(pathway_order)) {
    count <- sum(eligible$Pathway == pathway_order[[pathway_index]])
    if (!count) next
    expected <- c(expected, sprintf(
      "screen_p%03d_b%03d", pathway_index,
      seq_len(ceiling(count / batch_size))
    ))
  }
  checkpoint_dir <- file.path(run_dir, "03_screening", "checkpoints")
  all(vapply(expected, function(batch_id) {
    files <- .InterMetaAIScreeningCheckpointFiles(checkpoint_dir, batch_id)
    if (!length(files)) return(FALSE)
    checkpoint <- tryCatch(readRDS(utils::tail(files, 1L)), error = function(e) NULL)
    intact <- is.list(checkpoint) &&
      identical(checkpoint$checkpoint_schema_version,
                .InterMetaAIScreeningSchemaVersion) &&
      identical(checkpoint$batch_id, batch_id) &&
      identical(checkpoint$provider$name, provider_descriptor$name) &&
      identical(checkpoint$provider$model, provider_descriptor$model) &&
      is.character(checkpoint$content_sha256) &&
      identical(
        checkpoint$content_sha256,
        .InterMetaAIScreeningCheckpointContentSHA256(checkpoint)
      )
    intact && (
      checkpoint$status %in% c("complete", "partial") ||
      (!is.null(checkpoint$raw_response) &&
       !identical(checkpoint$validation_revision,
                  .InterMetaAIScreeningValidationRevision))
    )
  }, logical(1L)))
}

.InterMetaAIWorkflowSynthesisNoCallReady <- function(
    run_dir, context, screened, provider_descriptor) {
  records <- screened$records
  if (!is.data.frame(records) || !nrow(records)) return(TRUE)
  pathways <- as.character(context$pathways$Pathway)
  all(vapply(seq_along(pathways), function(i) {
    candidates <- tryCatch(
      .InterMetaAIScreenedSynthesisCandidateTable(records, pathways[[i]]),
      error = function(e) NULL
    )
    if (is.null(candidates)) return(FALSE)
    if (!nrow(candidates)) return(TRUE)
    files <- .InterMetaAIScreenedSynthesisCheckpointFiles(
      file.path(run_dir, "05_synthesis", "checkpoints"), i
    )
    if (!length(files)) return(FALSE)
    checkpoint <- tryCatch(readRDS(utils::tail(files, 1L)), error = function(e) NULL)
    intact <- is.list(checkpoint) &&
      identical(checkpoint$schema_version,
                .InterMetaAIScreenedSynthesisSchemaVersion) &&
      identical(checkpoint$pathway, pathways[[i]]) &&
      identical(checkpoint$provider$name, provider_descriptor$name) &&
      identical(checkpoint$provider$model, provider_descriptor$model) &&
      is.character(checkpoint$content_sha256) &&
      identical(
        checkpoint$content_sha256,
        .InterMetaAIScreenedSynthesisCheckpointContentSHA256(checkpoint)
      )
    intact && (
      identical(checkpoint$status, "complete") ||
      (!is.null(checkpoint$raw_parsed_response) &&
       !identical(checkpoint$validation_revision,
                  .InterMetaAIScreenedSynthesisValidationRevision))
    )
  }, logical(1L)))
}

.InterMetaAIWorkflowBlockProviderCalls <- function(provider) {
  InterMetaAIProvider(
    generate = function(messages, schema) {
      stop(
        "A model call was required but allow_model_calls is FALSE; no request was sent.",
        call. = FALSE
      )
    },
    name = provider$name, model = provider$model
  )
}

.InterMetaAIWorkflowWriteScreening <- function(screening, directory) {
  fields <- c(
    decisions = "article_screening.tsv", validation = "screening_validation.tsv",
    included = "included.tsv", excluded = "excluded_articles.tsv",
    unresolved = "unresolved.tsv", batches = "batches.tsv",
    candidates = "candidate_articles.tsv",
    batch_status = "batches.tsv"
  )
  for (field in names(fields)) {
    if (is.data.frame(screening[[field]])) {
      .InterMetaAIWorkflowAtomicTSV(
        screening[[field]], file.path(directory, fields[[field]])
      )
    }
  }
  invisible(NULL)
}

.InterMetaAIWorkflowWriteSynthesis <- function(synthesis, directory) {
  candidates <- list(
    synthesis$pathway_evidence_summary,
    synthesis$pathway_summary,
    synthesis$pathway_status,
    if (is.list(synthesis$interpretation)) synthesis$interpretation$exploration_summary else NULL
  )
  table <- NULL
  for (candidate in candidates) {
    if (is.data.frame(candidate)) {
      table <- candidate
      break
    }
  }
  if (!is.null(table)) {
    .InterMetaAIWorkflowAtomicTSV(
      table, file.path(directory, "pathway_evidence_summary.tsv")
    )
  }
  invisible(NULL)
}

.InterMetaAIWorkflowSubsetReport <- function(report, direction) {
  pathways <- as.character(report$top_pathways$Pathway[
    as.character(report$top_pathways$Direction) == direction
  ])
  if (!length(pathways)) stop("No pathways were available for figure direction: ", direction,
                              call. = FALSE)
  pathways <- head(pathways, 5L)
  out <- report
  out$top_pathways <- report$top_pathways[
    match(pathways, as.character(report$top_pathways$Pathway)), , drop = FALSE
  ]
  if (is.data.frame(report$exploration_summary) &&
      "Pathway" %in% names(report$exploration_summary)) {
    out$exploration_summary <- report$exploration_summary[
      match(pathways, as.character(report$exploration_summary$Pathway)), , drop = FALSE
    ]
  }
  if (is.data.frame(report$exploration_findings) &&
      "Pathway" %in% names(report$exploration_findings)) {
    out$exploration_findings <- report$exploration_findings[
      as.character(report$exploration_findings$Pathway) %in% pathways, , drop = FALSE
    ]
  }
  for (field in c("queries", "searches", "records")) {
    table <- out$evidence[[field]]
    if (is.data.frame(table) && "pathway" %in% names(table)) {
      out$evidence[[field]] <- table[as.character(table$pathway) %in% pathways,
                                    , drop = FALSE]
    }
  }
  out
}

.InterMetaAIWorkflowRenderFigures <- function(report, output_dir,
                                               font_family = NULL,
                                               text_mode = "compact_english") {
  directions <- unique(as.character(report$top_pathways$Direction))
  directions <- directions[!is.na(directions) & nzchar(directions)]
  if (!length(directions)) stop("The report contains no figure direction.", call. = FALSE)
  multiple <- length(directions) > 1L
  figures <- list()
  for (direction in directions) {
    suffix <- if (multiple) paste0("_", gsub("[^A-Za-z0-9_-]", "_", direction)) else ""
    subset <- .InterMetaAIWorkflowSubsetReport(report, direction)
    value <- .InterMetaAIWorkflowCall("InterMetaAIExplorationFigure", list(
      report = subset,
      pdf_file = file.path(output_dir, paste0("Top5", suffix, "_biology_exploration.pdf")),
      png_file = file.path(output_dir, paste0("Top5", suffix, "_biology_exploration.png")),
      render = TRUE, font_family = font_family, text_mode = text_mode
    ))
    if (is.list(value)) value$draw <- NULL
    figures[[direction]] <- value
  }
  figures
}

.InterMetaAIWorkflowReturn <- function(run_dir, state, report = NULL,
                                       figure = NULL) {
  structure(
    list(
      run_dir = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
      status = state$status, stages = state$stages,
      report = report, figure = figure
    ),
    class = c("InterMetaAIWorkflowResult", "list")
  )
}

.InterMetaAIWorkflowRun <- function(run_dir, config, state, provider = NULL,
                                    context = NULL, queries = NULL,
                                    evidence = NULL,
                                    allow_model_calls = FALSE,
                                    retry_failed = FALSE,
                                    pubmed_cache_dir_override = NULL) {
  paths <- .InterMetaAIWorkflowPaths(run_dir)
  run_stage <- function(stage, input_hash, expression, artifact) {
    if (.InterMetaAIWorkflowStageValid(state, run_dir, stage, input_hash)) {
      return(readRDS(artifact))
    }
    index <- match(stage, state$stages$Stage)
    if (state$stages$Status[index] %in%
        .InterMetaAIWorkflowReusableStatuses) {
      stop("Workflow checkpoint integrity failed for stage: ", stage,
           call. = FALSE)
    }
    state <<- .InterMetaAIWorkflowSetStage(
      state, run_dir, stage, "running", input_sha256 = input_hash,
      message = "stage started", increment = TRUE
    )
    value <- tryCatch(
      force(expression),
      error = function(e) {
        state <<- .InterMetaAIWorkflowSetStage(
          state, run_dir, stage, "failed", input_sha256 = input_hash,
          message = conditionMessage(e)
        )
        stop(e)
      }
    )
    .InterMetaAIWorkflowAtomicRDS(value, artifact)
    state <<- .InterMetaAIWorkflowComplete(
      state, run_dir, stage, input_hash, artifact
    )
    value
  }
  pause_before <- function(stage, reason) {
    index <- match(stage, state$stages$Stage)
    state$stages$Status[index] <- "paused"
    state$stages$Message[index] <- reason
    state$status <- "paused"
    state <<- .InterMetaAIWorkflowSaveState(state, run_dir)
    .InterMetaAIWorkflowReturn(run_dir, state)
  }

  context_input_hash <- if (!is.null(context)) {
    .InterMetaAIWorkflowHashObject(context)
  } else state$stages$Input.SHA256[match("context", state$stages$Stage)]
  if (is.null(context)) {
    if (!.InterMetaAIWorkflowStageValid(state, run_dir, "context",
                                        context_input_hash)) {
      stop("The saved context checkpoint is missing or was modified.", call. = FALSE)
    }
    context <- readRDS(paths$context)
  } else {
    context <- run_stage("context", context_input_hash, {
      if (!inherits(context, "InterMetaAIContext")) {
        stop("context must be created by InterMetaAIContext().", call. = FALSE)
      }
      context
    }, paths$context)
    .InterMetaAIWorkflowAtomicTSV(
      context$pathways, file.path(run_dir, "00_context", "pathways.tsv")
    )
  }

  query_input_hash <- .InterMetaAIWorkflowHashObject(list(
    context = .InterMetaAIWorkflowHashFile(paths$context),
    settings = config$query_settings,
    supplied = config$inputs$supplied_queries_sha256
  ))
  if (is.null(queries) && file.exists(paths$supplied_queries)) {
    queries <- readRDS(paths$supplied_queries)
  }
  queries <- run_stage("queries", query_input_hash, {
    if (is.null(queries)) {
      queries <- .InterMetaAIWorkflowCall("InterMetaBuildPubMedQueries", list(
        context = context,
        disease_mesh = config$query_settings$disease_mesh,
        pathway_synonyms = config$query_settings$pathway_synonyms,
        pathway_components = config$query_settings$pathway_components,
        include_components = config$query_settings$include_components,
        purpose = config$purpose
      ))
    }
    if (!is.data.frame(queries) || !nrow(queries) ||
        !all(c("query_id", "pathway", "layer", "query") %in% names(queries))) {
      stop("The PubMed query plan is empty or invalid.", call. = FALSE)
    }
    queries
  }, paths$queries)
  .InterMetaAIWorkflowAtomicTSV(
    queries, file.path(run_dir, "01_queries", "planned_queries.tsv")
  )

  pubmed_input_hash <- .InterMetaAIWorkflowHashObject(list(
    queries = .InterMetaAIWorkflowHashFile(paths$queries),
    settings = config$pubmed_settings,
    supplied = config$inputs$supplied_evidence_sha256
  ))
  if (is.null(evidence) && file.exists(paths$supplied_evidence)) {
    evidence <- readRDS(paths$supplied_evidence)
  }
  evidence <- run_stage("pubmed", pubmed_input_hash, {
    if (is.null(evidence)) {
      email <- trimws(Sys.getenv(config$pubmed_settings$email_env, unset = ""))
      ncbi_key <- trimws(Sys.getenv(config$pubmed_settings$api_key_env,
                                    unset = ""))
      cache_dir <- .InterMetaAIWorkflowResolveCacheDir(
        config, run_dir, override = pubmed_cache_dir_override
      )
      evidence <- .InterMetaAIWorkflowCall("InterMetaPubMedSearch", list(
        context = context, queries = queries,
        retmax = config$pubmed_settings$retmax,
        page_size = config$pubmed_settings$page_size,
        sort = config$pubmed_settings$sort,
        email = if (nzchar(email)) email else NULL,
        tool = config$pubmed_settings$tool,
        api_key = ncbi_key,
        cache_dir = cache_dir,
        refresh = FALSE, purpose = config$purpose
      ))
    }
    if (!inherits(evidence, "InterMetaPubMedEvidence") || !is.list(evidence)) {
      stop("evidence must be an InterMetaPubMedEvidence object.", call. = FALSE)
    }
    if (!.InterMetaAIWorkflowQueriesMatch(queries, evidence$queries)) {
      stop("Supplied PubMed evidence does not match the frozen query plan.",
           call. = FALSE)
    }
    evidence
  }, paths$evidence)
  pubmed_dir <- dirname(paths$evidence)
  for (entry in c("queries", "searches", "records", "response_hashes")) {
    if (is.data.frame(evidence[[entry]])) {
      .InterMetaAIWorkflowAtomicTSV(
        evidence[[entry]], file.path(pubmed_dir, paste0(entry, ".tsv"))
      )
    }
  }
  candidates <- .InterMetaAIWorkflowCandidates(evidence)
  .InterMetaAIWorkflowAtomicTSV(
    candidates, file.path(pubmed_dir, "candidate_articles.tsv")
  )
  if (is.data.frame(evidence$records)) {
    link_columns <- intersect(
      c("query_id", "pathway", "direction", "layer", "term_scope", "PMID"),
      names(evidence$records)
    )
    .InterMetaAIWorkflowAtomicTSV(
      evidence$records[, link_columns, drop = FALSE],
      file.path(pubmed_dir, "retrieval_links.tsv")
    )
  }
  if (!.InterMetaAIWorkflowPubMedReady(evidence)) {
    state <- .InterMetaAIWorkflowSetStage(
      state, run_dir, "pubmed", "failed", input_sha256 = pubmed_input_hash,
      output_sha256 = .InterMetaAIWorkflowHashFile(paths$evidence),
      artifact = .InterMetaAIWorkflowRelative(paths$evidence, run_dir),
      message = paste(
        "PubMed retrieval is incomplete. No model call was made; resume reuses",
        "valid query caches."
      )
    )
    state$status <- "failed"
    state <- .InterMetaAIWorkflowSaveState(state, run_dir)
    stop("PubMed retrieval is incomplete; no model call was made.", call. = FALSE)
  }
  if (.InterMetaAIWorkflowPubMedTruncated(evidence)) {
    truncated_count <- sum(evidence$searches$truncated %in% TRUE)
    truncation_message <- sprintf(
      paste(
        "PubMed retrieval completed, but %d of %d searches reached the",
        "configured retrieval limit. This is a bounded evidence snapshot,",
        "not an exhaustive literature search."
      ),
      truncated_count, nrow(evidence$searches)
    )
    pubmed_index <- match("pubmed", state$stages$Stage)
    if (!identical(state$stages$Status[pubmed_index],
                   "complete_with_truncation") ||
        !identical(state$stages$Message[pubmed_index], truncation_message)) {
      state <- .InterMetaAIWorkflowComplete(
        state, run_dir, "pubmed", pubmed_input_hash, paths$evidence,
        status = "complete_with_truncation", message = truncation_message
      )
    }
  }

  screening_input_hash <- .InterMetaAIWorkflowHashObject(list(
    context = .InterMetaAIWorkflowHashFile(paths$context),
    evidence = .InterMetaAIWorkflowHashFile(paths$evidence),
    provider = config$provider, settings = config$screening_settings,
    contract = config$contract_revision
  ))
  screening_index <- match("screening", state$stages$Stage)
  screening_no_call <- .InterMetaAIWorkflowScreeningNoCallReady(
    run_dir, context, candidates, config$provider,
    config$screening_settings$batch_size
  )
  if (!.InterMetaAIWorkflowStageValid(
        state, run_dir, "screening", screening_input_hash)) {
    prior_status <- state$stages$Status[screening_index]
    if (!allow_model_calls && !screening_no_call) {
      return(pause_before(
        "screening",
        "Set allow_model_calls = TRUE to authorize new DeepSeek screening calls."
      ))
    }
    if (prior_status %in% c("failed", "delivery_uncertain") &&
        !retry_failed && !screening_no_call) {
      return(pause_before(
        "screening",
        paste("Screening has a failed or uncertain attempt; explicitly set",
              "retry_failed = TRUE to authorize another paid attempt.")
      ))
    }
  }
  selected_provider <- NULL
  if (!.InterMetaAIWorkflowStageValid(
        state, run_dir, "screening", screening_input_hash)) {
    selected_provider <- .InterMetaAIWorkflowResolveProvider(config, provider)
    if (!allow_model_calls) {
      selected_provider <- .InterMetaAIWorkflowBlockProviderCalls(
        selected_provider
      )
    }
  }
  screening <- run_stage("screening", screening_input_hash, {
    .InterMetaAIWorkflowCall("InterMetaAIScreenEvidence", list(
      context = context, evidence = evidence, provider = selected_provider,
      language = config$language,
      batch_size = config$screening_settings$batch_size,
      checkpoint_dir = file.path(run_dir, "03_screening", "checkpoints"),
      resume = TRUE, retry_failed = retry_failed,
      continue_on_error = config$screening_settings$continue_on_error
    ))
  }, paths$screening)
  .InterMetaAIWorkflowWriteScreening(screening, dirname(paths$screening))
  if (identical(screening$status, "complete_with_quarantine")) {
    state <- .InterMetaAIWorkflowComplete(
      state, run_dir, "screening", screening_input_hash, paths$screening,
      status = "complete_with_quarantine",
      message = paste(
        "Screening completed; invalid article rows were quarantined and",
        "excluded from synthesis."
      )
    )
  }
  if (!.InterMetaAIWorkflowCompleteFlag(screening)) {
    state <- .InterMetaAIWorkflowSetStage(
      state, run_dir, "screening", "failed",
      input_sha256 = screening_input_hash,
      output_sha256 = .InterMetaAIWorkflowHashFile(paths$screening),
      artifact = .InterMetaAIWorkflowRelative(paths$screening, run_dir),
      message = "Screening returned an incomplete result; no synthesis call was made."
    )
    state$status <- "failed"
    state <- .InterMetaAIWorkflowSaveState(state, run_dir)
    return(.InterMetaAIWorkflowReturn(run_dir, state))
  }

  verified_input_hash <- .InterMetaAIWorkflowHashObject(list(
    evidence = .InterMetaAIWorkflowHashFile(paths$evidence),
    screening = .InterMetaAIWorkflowHashFile(paths$screening)
  ))
  screened <- run_stage("verified_evidence", verified_input_hash, {
    .InterMetaAIWorkflowCall("InterMetaAIScreenedEvidence", list(
      screening = screening, require_complete = TRUE
    ))
  }, paths$screened)
  verified_dir <- dirname(paths$screened)
  verified_table <- if (is.data.frame(screening$included)) {
    screening$included
  } else if (is.data.frame(screened$records)) screened$records else NULL
  if (is.data.frame(verified_table)) {
    .InterMetaAIWorkflowAtomicTSV(
      verified_table, file.path(verified_dir, "verified_findings.tsv")
    )
  }
  if (is.data.frame(screening$excluded)) {
    .InterMetaAIWorkflowAtomicTSV(
      screening$excluded, file.path(verified_dir, "excluded_articles.tsv")
    )
  }

  synthesis_input_hash <- .InterMetaAIWorkflowHashObject(list(
    context = .InterMetaAIWorkflowHashFile(paths$context),
    screened = .InterMetaAIWorkflowHashFile(paths$screened),
    screening = .InterMetaAIWorkflowHashFile(paths$screening),
    provider = config$provider, settings = config$synthesis_settings,
    contract = config$contract_revision
  ))
  synthesis_index <- match("synthesis", state$stages$Stage)
  synthesis_no_call <- .InterMetaAIWorkflowSynthesisNoCallReady(
    run_dir, context, screened, config$provider
  )
  if (!.InterMetaAIWorkflowStageValid(
        state, run_dir, "synthesis", synthesis_input_hash)) {
    prior_status <- state$stages$Status[synthesis_index]
    if (!allow_model_calls && !synthesis_no_call) {
      return(pause_before(
        "synthesis",
        "Set allow_model_calls = TRUE to authorize new DeepSeek synthesis calls."
      ))
    }
    if (prior_status %in% c("failed", "delivery_uncertain") &&
        !retry_failed && !synthesis_no_call) {
      return(pause_before(
        "synthesis",
        paste("Synthesis has a failed or uncertain attempt; explicitly set",
              "retry_failed = TRUE to authorize another paid attempt.")
      ))
    }
    if (is.null(selected_provider)) {
      selected_provider <- .InterMetaAIWorkflowResolveProvider(config, provider)
    }
    if (!allow_model_calls) {
      selected_provider <- .InterMetaAIWorkflowBlockProviderCalls(
        selected_provider
      )
    }
  }
  synthesis <- run_stage("synthesis", synthesis_input_hash, {
    synthesis_result <- .InterMetaAIWorkflowCall("InterMetaAISynthesizeScreened", list(
      context = context, screening = screening, provider = selected_provider,
      language = config$language,
      outdir = file.path(run_dir, "05_synthesis", "checkpoints"),
      max_findings_per_pathway = config$synthesis_settings$max_findings_per_pathway,
      reuse = TRUE, retry_failed = retry_failed,
      continue_on_error = config$synthesis_settings$continue_on_error
    ))
    .InterMetaAIWorkflowSynthesisForStageArtifact(synthesis_result)
  }, paths$synthesis)
  .InterMetaAIWorkflowWriteSynthesis(synthesis, dirname(paths$synthesis))
  if (!.InterMetaAIWorkflowCompleteFlag(synthesis)) {
    state <- .InterMetaAIWorkflowSetStage(
      state, run_dir, "synthesis", "failed",
      input_sha256 = synthesis_input_hash,
      output_sha256 = .InterMetaAIWorkflowHashFile(paths$synthesis),
      artifact = .InterMetaAIWorkflowRelative(paths$synthesis, run_dir),
      message = "Synthesis returned an incomplete result; no final report was published."
    )
    state$status <- "failed"
    state <- .InterMetaAIWorkflowSaveState(state, run_dir)
    return(.InterMetaAIWorkflowReturn(run_dir, state))
  }

  interpretation <- if (is.list(synthesis) &&
                        !is.null(synthesis$interpretation)) {
    synthesis$interpretation
  } else synthesis
  report_input_hash <- .InterMetaAIWorkflowHashObject(list(
    context = .InterMetaAIWorkflowHashFile(paths$context),
    screened = .InterMetaAIWorkflowHashFile(paths$screened),
    synthesis = .InterMetaAIWorkflowHashFile(paths$synthesis),
    language = config$language, purpose = config$purpose,
    report_schema_version = .InterMetaAIReportSchemaVersion
  ))
  report <- run_stage("report", report_input_hash, {
    report_index <- match("report", state$stages$Stage)
    report_bundle <- file.path(
      run_dir, "06_report", "outputs",
      paste0("attempt-", sprintf("%03d", state$stages$Attempt[report_index]))
    )
    .InterMetaAIWorkflowCall("InterMetaAIReport", list(
      context = context, evidence = screened, interpretation = interpretation,
      language = config$language, purpose = config$purpose,
      outdir = report_bundle
    ))
  }, paths$report)

  figure_input_hash <- .InterMetaAIWorkflowHashObject(list(
    report = .InterMetaAIWorkflowHashFile(paths$report),
    render = config$figure_settings$render,
    text_mode = config$figure_settings$text_mode,
    font_family = config$figure_settings$font_family
  ))
  if (!isTRUE(config$figure_settings$render)) {
    figure <- list(rendered = FALSE, reason = "render = FALSE")
    if (!.InterMetaAIWorkflowStageValid(state, run_dir, "figure",
                                        figure_input_hash)) {
      .InterMetaAIWorkflowAtomicRDS(figure, paths$figure)
      state <- .InterMetaAIWorkflowComplete(
        state, run_dir, "figure", figure_input_hash, paths$figure,
        status = "skipped", message = "Figure rendering was disabled."
      )
    }
  } else {
    figure <- run_stage("figure", figure_input_hash, {
      figure_index <- match("figure", state$stages$Stage)
      output_dir <- file.path(
        run_dir, "07_figure", "outputs",
        paste0("attempt-", sprintf("%03d", state$stages$Attempt[figure_index]))
      )
      if (!dir.exists(output_dir) &&
          !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
        stop("Could not create the figure output directory.", call. = FALSE)
      }
      if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
        stop("Figure attempt directory is unexpectedly non-empty.", call. = FALSE)
      }
      rendered <- .InterMetaAIWorkflowRenderFigures(
        report, output_dir,
        font_family = config$figure_settings$font_family,
        text_mode = config$figure_settings$text_mode
      )
      .InterMetaAIWorkflowFigureForStorage(
        rendered, run_dir, anchors = "07_figure"
      )
    }, paths$figure)
  }

  portable_figure <- .InterMetaAIWorkflowFigureForStorage(
    figure, run_dir, anchors = "07_figure", strict = FALSE
  )
  if (!identical(portable_figure, figure)) {
    figure_index <- match("figure", state$stages$Stage)
    figure_status <- state$stages$Status[figure_index]
    figure_message <- state$stages$Message[figure_index]
    .InterMetaAIWorkflowAtomicRDS(portable_figure, paths$figure)
    state <- .InterMetaAIWorkflowComplete(
      state, run_dir, "figure", figure_input_hash, paths$figure,
      status = figure_status, message = figure_message
    )
  }
  figure <- portable_figure

  state$status <- "complete"
  state <- .InterMetaAIWorkflowSaveState(state, run_dir)
  .InterMetaAIWorkflowReturn(
    run_dir, state, report,
    .InterMetaAIWorkflowFigureForSession(
      figure, run_dir, anchors = "07_figure"
    )
  )
}

#' Run the recoverable InterMeta PubMed and DeepSeek workflow
#'
#' Model calls are disabled unless `allow_model_calls = TRUE` is supplied.
#' The run directory contains a hash-bound context, query plan, PubMed snapshot,
#' screening, verified evidence, synthesis, report, and figure checkpoints.
#' API-key values are never written to the run directory.
#'
#' @export
InterMetaAIWorkflow <- function(
    context, run_dir, evidence = NULL, queries = NULL, provider = NULL,
    allow_model_calls = FALSE, language = c("en", "zh"),
    purpose = c("biology_exploration", "prognostic_validation"),
    disease_mesh = NULL, pathway_synonyms = NULL, pathway_components = NULL,
    include_components = TRUE, retmax = 100L, page_size = 100L,
    sort = c("relevance", "pub_date"),
    pubmed_cache_dir = file.path(dirname(run_dir), "pubmed_cache"),
    ncbi_email_env = "NCBI_EMAIL", ncbi_api_key_env = "NCBI_API_KEY",
    pubmed_tool = "InterMeta", screen_batch_size = 20L,
    max_findings_per_pathway = 5L, model = "deepseek-v4-pro",
    api_key_env = "DEEPSEEK_API_KEY", model_max_tokens = 32768L,
    model_timeout = 900, thinking = c("enabled", "disabled"),
    reasoning_effort = c("low", "high", "max"),
    continue_on_error = TRUE, render = TRUE,
    figure_text_mode = c("compact_english", "full"), font_family = NULL) {
  portable_cache_default <- missing(pubmed_cache_dir)
  run_dir <- .InterMetaAIWorkflowScalar(run_dir, "run_dir")
  allow_model_calls <- .InterMetaAIWorkflowFlag(
    allow_model_calls, "allow_model_calls"
  )
  include_components <- .InterMetaAIWorkflowFlag(
    include_components, "include_components"
  )
  continue_on_error <- .InterMetaAIWorkflowFlag(
    continue_on_error, "continue_on_error"
  )
  render <- .InterMetaAIWorkflowFlag(render, "render")
  language <- match.arg(language)
  purpose <- match.arg(purpose)
  if (!identical(purpose, "biology_exploration")) {
    stop(paste(
      "The two-stage screening workflow currently supports",
      "purpose = 'biology_exploration' only; use InterMetaAIInterpret()",
      "for the legacy prognostic_validation workflow."
    ), call. = FALSE)
  }
  sort <- match.arg(sort)
  thinking <- match.arg(thinking)
  reasoning_effort <- match.arg(reasoning_effort)
  figure_text_mode <- match.arg(figure_text_mode)
  retmax <- .InterMetaAIWorkflowInteger(retmax, "retmax")
  page_size <- .InterMetaAIWorkflowInteger(page_size, "page_size")
  screen_batch_size <- .InterMetaAIWorkflowInteger(
    screen_batch_size, "screen_batch_size"
  )
  max_findings_per_pathway <- .InterMetaAIWorkflowInteger(
    max_findings_per_pathway, "max_findings_per_pathway"
  )
  if (max_findings_per_pathway > 5L) {
    stop("max_findings_per_pathway must be between 1 and 5.", call. = FALSE)
  }
  model_max_tokens <- .InterMetaAIWorkflowInteger(
    model_max_tokens, "model_max_tokens"
  )
  if (!is.numeric(model_timeout) || length(model_timeout) != 1L ||
      is.na(model_timeout) || !is.finite(model_timeout) || model_timeout <= 0) {
    stop("model_timeout must be one positive number.", call. = FALSE)
  }
  model <- .InterMetaAIWorkflowScalar(model, "model")
  api_key_env <- .InterMetaAIWorkflowScalar(api_key_env, "api_key_env")
  ncbi_email_env <- .InterMetaAIWorkflowScalar(ncbi_email_env,
                                                "ncbi_email_env")
  ncbi_api_key_env <- .InterMetaAIWorkflowScalar(
    ncbi_api_key_env, "ncbi_api_key_env"
  )
  env_names <- c(
    api_key_env = api_key_env,
    ncbi_email_env = ncbi_email_env,
    ncbi_api_key_env = ncbi_api_key_env
  )
  invalid_env <- !grepl("^[A-Za-z_][A-Za-z0-9_]*$", env_names)
  if (any(invalid_env)) {
    stop(
      names(env_names)[which(invalid_env)[1L]],
      " must be an environment-variable name, not a credential value.",
      call. = FALSE
    )
  }
  pubmed_tool <- .InterMetaAIWorkflowScalar(pubmed_tool, "pubmed_tool")
  if (!is.null(font_family)) {
    font_family <- .InterMetaAIWorkflowScalar(font_family, "font_family")
  }
  if (!inherits(context, "InterMetaAIContext")) {
    stop("context must be created by InterMetaAIContext().", call. = FALSE)
  }
  descriptor <- .InterMetaAIWorkflowProviderDescriptor(provider, model)
  if (file.exists(run_dir) || dir.exists(run_dir)) {
    stop("run_dir already exists; use InterMetaAIResume(): ", run_dir,
         call. = FALSE)
  }
  if (!dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create run_dir: ", run_dir, call. = FALSE)
  }
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  paths <- .InterMetaAIWorkflowPaths(run_dir)
  if (!is.null(queries)) .InterMetaAIWorkflowAtomicRDS(queries, paths$supplied_queries)
  if (!is.null(evidence)) .InterMetaAIWorkflowAtomicRDS(evidence, paths$supplied_evidence)
  cache_spec <- .InterMetaAIWorkflowCacheSpec(
    if (portable_cache_default) "pubmed_cache" else pubmed_cache_dir,
    portable_default = portable_cache_default
  )
  config <- list(
    schema_version = .InterMetaAIWorkflowSchemaVersion,
    contract_revision = "intermeta-ai-workflow-v1",
    purpose = purpose, language = language, provider = descriptor,
    model_settings = list(
      api_key_env = api_key_env, max_tokens = model_max_tokens,
      timeout = as.numeric(model_timeout), thinking = thinking,
      reasoning_effort = reasoning_effort
    ),
    query_settings = list(
      disease_mesh = disease_mesh, pathway_synonyms = pathway_synonyms,
      pathway_components = pathway_components,
      include_components = include_components
    ),
    pubmed_settings = list(
      retmax = retmax, page_size = page_size, sort = sort,
      cache_dir = cache_spec$cache_dir,
      cache_anchor = cache_spec$cache_anchor,
      email_env = ncbi_email_env, api_key_env = ncbi_api_key_env,
      tool = pubmed_tool
    ),
    screening_settings = list(
      batch_size = screen_batch_size,
      continue_on_error = continue_on_error
    ),
    synthesis_settings = list(
      max_findings_per_pathway = max_findings_per_pathway,
      continue_on_error = continue_on_error
    ),
    figure_settings = list(
      render = render, text_mode = figure_text_mode,
      font_family = font_family
    ),
    inputs = list(
      supplied_queries_sha256 = if (file.exists(paths$supplied_queries))
        .InterMetaAIWorkflowHashFile(paths$supplied_queries) else NULL,
      supplied_evidence_sha256 = if (file.exists(paths$supplied_evidence))
        .InterMetaAIWorkflowHashFile(paths$supplied_evidence) else NULL
    )
  )
  if (!.InterMetaAIWorkflowPlainData(config)) {
    stop("Workflow configuration must contain plain serializable data only.",
         call. = FALSE)
  }
  secrets <- c(
    Sys.getenv(api_key_env, unset = ""),
    Sys.getenv(ncbi_api_key_env, unset = "")
  )
  if (.InterMetaAIWorkflowContainsSecret(config, secrets)) {
    stop("Refusing to persist workflow configuration containing credential material.",
         call. = FALSE)
  }
  .InterMetaAIWorkflowAtomicRDS(config, paths$config)
  state <- .InterMetaAIWorkflowNewState()
  state$config_sha256 <- .InterMetaAIWorkflowHashFile(paths$config)
  state <- .InterMetaAIWorkflowSaveState(state, run_dir)
  .InterMetaAIWorkflowRun(
    run_dir, config, state, provider = provider, context = context,
    queries = queries, evidence = evidence,
    allow_model_calls = allow_model_calls, retry_failed = FALSE
  )
}

#' Resume an existing InterMeta AI workflow
#'
#' Failed or delivery-uncertain model calls are not retried unless both
#' `allow_model_calls` and `retry_failed` are explicitly `TRUE`.
#'
#' @export
InterMetaAIResume <- function(run_dir, provider = NULL,
                              allow_model_calls = FALSE,
                              retry_failed = FALSE,
                              pubmed_cache_dir = NULL) {
  run_dir <- .InterMetaAIWorkflowScalar(run_dir, "run_dir")
  allow_model_calls <- .InterMetaAIWorkflowFlag(
    allow_model_calls, "allow_model_calls"
  )
  retry_failed <- .InterMetaAIWorkflowFlag(retry_failed, "retry_failed")
  if (!is.null(pubmed_cache_dir)) {
    pubmed_cache_dir <- .InterMetaAIWorkflowScalar(
      pubmed_cache_dir, "pubmed_cache_dir"
    )
  }
  if (retry_failed && !allow_model_calls) {
    stop("retry_failed = TRUE also requires allow_model_calls = TRUE.",
         call. = FALSE)
  }
  if (!dir.exists(run_dir)) stop("run_dir does not exist.", call. = FALSE)
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  paths <- .InterMetaAIWorkflowPaths(run_dir)
  if (!file.exists(paths$config) || !file.exists(paths$state)) {
    stop("run_dir is not an InterMeta AI workflow directory.", call. = FALSE)
  }
  config <- readRDS(paths$config)
  state <- readRDS(paths$state)
  if (!identical(config$schema_version, .InterMetaAIWorkflowSchemaVersion) ||
      !identical(state$schema_version, .InterMetaAIWorkflowSchemaVersion)) {
    stop("Unsupported workflow schema version.", call. = FALSE)
  }
  if (!identical(.InterMetaAIWorkflowHashFile(paths$config),
                 state$config_sha256)) {
    stop("Workflow configuration failed its SHA256 check.", call. = FALSE)
  }
  if (!is.data.frame(state$stages) ||
      !identical(as.character(state$stages$Stage),
                 .InterMetaAIWorkflowStages)) {
    stop("Workflow state is malformed.", call. = FALSE)
  }
  running <- state$stages$Status == "running" &
    state$stages$Stage %in% c("screening", "synthesis")
  if (any(running)) {
    state$stages$Status[running] <- "delivery_uncertain"
    state$stages$Message[running] <- paste(
      "The previous process ended during a model stage; delivery/billing is",
      "uncertain and no automatic retry will be made."
    )
    state <- .InterMetaAIWorkflowSaveState(state, run_dir)
  }
  .InterMetaAIWorkflowRun(
    run_dir, config, state, provider = provider,
    allow_model_calls = allow_model_calls, retry_failed = retry_failed,
    pubmed_cache_dir_override = pubmed_cache_dir
  )
}

#' Inspect an InterMeta AI workflow without contacting external services
#'
#' @export
InterMetaAIStatus <- function(run_dir) {
  run_dir <- .InterMetaAIWorkflowScalar(run_dir, "run_dir")
  if (!dir.exists(run_dir)) stop("run_dir does not exist.", call. = FALSE)
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  paths <- .InterMetaAIWorkflowPaths(run_dir)
  if (!file.exists(paths$config) || !file.exists(paths$state)) {
    stop("run_dir is not an InterMeta AI workflow directory.", call. = FALSE)
  }
  config <- readRDS(paths$config)
  state <- readRDS(paths$state)
  integrity <- vapply(.InterMetaAIWorkflowStages, function(stage) {
    index <- match(stage, state$stages$Stage)
    if (is.na(index) || !state$stages$Status[index] %in%
        .InterMetaAIWorkflowReusableStatuses) return(NA)
    .InterMetaAIWorkflowStageValid(
      state, run_dir, stage, state$stages$Input.SHA256[index]
    )
  }, logical(1L))
  stages <- state$stages
  stages$Integrity.OK <- unname(integrity)
  structure(
    list(
      run_dir = run_dir, status = state$status,
      config_sha256_ok = identical(
        .InterMetaAIWorkflowHashFile(paths$config), state$config_sha256
      ),
      provider = config$provider, purpose = config$purpose,
      language = config$language, stages = stages,
      updated_at = state$updated_at
    ),
    class = c("InterMetaAIWorkflowStatus", "list")
  )
}

#' Replay a completed InterMeta AI workflow without PubMed or model calls
#'
#' The saved context, verified evidence, and parsed synthesis are revalidated by
#' the current local report code.  Outputs are written to a new directory.
#'
#' @export
InterMetaAIReplay <- function(run_dir,
                              outdir = file.path(run_dir, "replay"),
                              render = TRUE,
                              text_mode = c("compact_english", "full"),
                              font_family = NULL) {
  run_dir <- .InterMetaAIWorkflowScalar(run_dir, "run_dir")
  outdir <- .InterMetaAIWorkflowScalar(outdir, "outdir")
  render <- .InterMetaAIWorkflowFlag(render, "render")
  text_mode <- match.arg(text_mode)
  if (!is.null(font_family)) {
    font_family <- .InterMetaAIWorkflowScalar(font_family, "font_family")
  }
  status <- InterMetaAIStatus(run_dir)
  if (!isTRUE(status$config_sha256_ok)) {
    stop("Workflow configuration failed its SHA256 check.", call. = FALSE)
  }
  required <- c("context", "verified_evidence", "synthesis")
  indices <- match(required, status$stages$Stage)
  if (anyNA(indices) ||
      any(!status$stages$Status[indices] %in%
          .InterMetaAIWorkflowSuccessfulStatuses) ||
      any(!status$stages$Integrity.OK[indices])) {
    stop("Replay requires intact completed context, verified evidence, and synthesis checkpoints.",
         call. = FALSE)
  }
  if (file.exists(outdir) || dir.exists(outdir)) {
    stop("outdir already exists; replay never overwrites outputs.", call. = FALSE)
  }
  parent <- dirname(outdir)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create replay parent directory.", call. = FALSE)
  }
  parent <- normalizePath(parent, winslash = "/", mustWork = TRUE)
  outdir <- file.path(parent, basename(outdir))
  stage <- tempfile(paste0(".", basename(outdir), "-staging-"), tmpdir = parent)
  if (!dir.create(stage, showWarnings = FALSE)) {
    stop("Could not create replay staging directory.", call. = FALSE)
  }
  published <- FALSE
  on.exit({
    if (!published && dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  paths <- .InterMetaAIWorkflowPaths(normalizePath(run_dir, winslash = "/",
                                                   mustWork = TRUE))
  config <- readRDS(paths$config)
  context <- readRDS(paths$context)
  screened <- readRDS(paths$screened)
  synthesis <- readRDS(paths$synthesis)
  interpretation <- if (is.list(synthesis) &&
                        !is.null(synthesis$interpretation)) {
    synthesis$interpretation
  } else synthesis
  report <- .InterMetaAIWorkflowCall("InterMetaAIReport", list(
    context = context, evidence = screened, interpretation = interpretation,
    language = config$language, purpose = config$purpose,
    outdir = file.path(stage, "report")
  ))
  figure <- list(rendered = FALSE)
  if (render) {
    figure_dir <- file.path(stage, "figures")
    if (!dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create replay figure directory.", call. = FALSE)
    }
    figure <- .InterMetaAIWorkflowRenderFigures(
      report, figure_dir, font_family = font_family, text_mode = text_mode
    )
  }
  stored_figure <- .InterMetaAIWorkflowFigureForStorage(
    figure, stage, anchors = "figures"
  )
  replay <- list(
    schema_version = .InterMetaAIWorkflowSchemaVersion,
    source_run = .InterMetaAIWorkflowReplaySourceForStorage(
      run_dir, outdir
    ),
    source_hashes = list(
      context = .InterMetaAIWorkflowHashFile(paths$context),
      verified_evidence = .InterMetaAIWorkflowHashFile(paths$screened),
      synthesis = .InterMetaAIWorkflowHashFile(paths$synthesis)
    ),
    replayed_at = .InterMetaAIWorkflowNow(),
    model_calls = 0L, pubmed_calls = 0L, report = report,
    figure = stored_figure
  )
  .InterMetaAIWorkflowAtomicRDS(replay, file.path(stage, "replay.rds"))
  files <- list.files(stage, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[!file.info(files)$isdir]
  manifest <- data.frame(
    File = substring(
      normalizePath(files, winslash = "/", mustWork = TRUE),
      nchar(normalizePath(stage, winslash = "/", mustWork = TRUE)) + 2L
    ),
    Bytes = as.numeric(file.info(files)$size),
    SHA256 = unname(vapply(files, .InterMetaAIWorkflowHashFile, character(1L))),
    stringsAsFactors = FALSE
  )
  .InterMetaAIWorkflowAtomicTSV(manifest, file.path(stage, "manifest.tsv"))
  if (!file.rename(stage, outdir)) {
    stop("Could not atomically publish replay directory.", call. = FALSE)
  }
  published <- TRUE
  source_run_session <- .InterMetaAIWorkflowReplaySourceForSession(
    replay$source_run, outdir
  )
  if (is.na(source_run_session)) {
    source_run_session <- normalizePath(
      run_dir, winslash = "/", mustWork = TRUE
    )
  }
  invisible(structure(list(
    run_dir = source_run_session,
    source_run = source_run_session,
    replay_dir = normalizePath(outdir, winslash = "/", mustWork = TRUE),
    report = report,
    figure = .InterMetaAIWorkflowFigureForSession(
      stored_figure, outdir, anchors = "figures"
    )
  ), class = c("InterMetaAIWorkflowReplay", "list")))
}
