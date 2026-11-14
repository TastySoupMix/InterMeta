# Optional DeepSeek adapter. Only aggregate messages supplied by AIReport.R
# are sent; the constructor never reads, retains, or serializes credentials.

.InterMetaDeepSeekEndpoint <- "https://api.deepseek.com/chat/completions"

.InterMetaDeepSeekRequire <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The optional package '", package,
         "' is required to generate a DeepSeek interpretation.", call. = FALSE)
  }
  invisible(NULL)
}

.InterMetaDeepSeekJSON <- function(x) {
  tryCatch(
    as.character(jsonlite::toJSON(
      x, auto_unbox = TRUE, dataframe = "rows", null = "null", na = "null",
      digits = NA, force = TRUE
    )),
    error = function(e) stop("DeepSeek request data could not be encoded as JSON.",
                             call. = FALSE)
  )
}

.InterMetaDeepSeekSchemaArrays <- function(x) {
  if (!is.list(x)) return(x)
  for (i in seq_along(x)) {
    field <- if (is.null(names(x))) "" else names(x)[i]
    if (!is.na(field) && field %in% c("required", "enum") &&
        is.atomic(x[[i]])) {
      # jsonlite auto_unbox must not turn a singleton JSON Schema array into
      # a string (the top-level required = "pathways" is such a singleton).
      x[[i]] <- as.list(x[[i]])
    } else {
      x[[i]] <- .InterMetaDeepSeekSchemaArrays(x[[i]])
    }
  }
  x
}

.InterMetaDeepSeekMessages <- function(messages, schema) {
  if (!is.list(messages) || !length(messages) ||
      !is.list(schema) || !length(schema)) {
    stop("DeepSeek requires a non-empty messages list and JSON schema list.",
         call. = FALSE)
  }
  out <- lapply(messages, function(message) {
    if (!is.list(message) || !is.character(message$role) ||
        length(message$role) != 1L || is.na(message$role) ||
        !message$role %in% c("system", "user", "assistant")) {
      stop("DeepSeek messages require a valid system, user, or assistant role.",
           call. = FALSE)
    }
    content <- message$content
    if (is.list(content)) {
      content <- .InterMetaDeepSeekJSON(content)
    }
    if (!is.character(content) || length(content) != 1L || is.na(content) ||
        !nzchar(trimws(content))) {
      stop("DeepSeek message content must be a non-empty string or list.",
           call. = FALSE)
    }
    list(role = message$role, content = content)
  })
  biology_example <- if (!is.null(schema$properties$pathways$items$properties$biology_synthesis)) {
    paste0(',"biology_synthesis":{"mechanism_summary":"No mechanism evidence is available.",',
           '"limitations":"No biological inference can be made.","findings":[]}')
  } else ""
  example <- paste0(
    '{"pathways":[{"pathway":"EXACT_PATHWAY_FROM_CONTEXT",',
    '"literature_consistency":"insufficient_evidence",',
    '"summary":"The supplied evidence is insufficient.","claims":[],"assessments":[]',
    biology_example, '}]}'
  )
  if (!is.null(schema$properties$pathways$items$properties$exploration)) {
    example <- paste0('{"pathways":[{"pathway":"EXACT_PATHWAY_FROM_CONTEXT",',
      '"exploration":{"biological_summary":"No source-grounded findings are available.",',
      '"context_dependence":"Not assessed.","research_hypothesis":"No hypothesis is supported yet.",',
      '"findings":[]}}]}')
  }
  if (identical(schema$name, "intermeta_screened_pathway_synthesis")) {
    example <- paste0('{"pathways":[{"pathway":"EXACT_PATHWAY_FROM_CONTEXT",',
      '"biological_plausibility":"limited_evidence",',
      '"biological_summary":"The verified findings provide limited pathway-relevant evidence.",',
      '"context_dependence":"Evidence is context dependent.",',
      '"research_hypothesis":"Test the pathway in a disease-relevant perturbation model.",',
      '"limitations":"Only the supplied verified findings were assessed.",',
      '"selected_candidate_ids":["EXACT_CANDIDATE_ID"]}]}')
  }
  if (!is.null(schema$name) &&
      schema$name %in% c("intermeta_article_screening", "intermeta_literature_screening")) {
    example <- paste0('{"pathways":[{"pathway":"EXACT_PATHWAY_FROM_BATCH",',
      '"articles":[{"candidate_id":"EXACT_CANDIDATE_ID",',
      '"pathway":"EXACT_PATHWAY_FROM_BATCH","pmid":"EXACT_PMID",',
      '"include":false,"exclusion_reason":"insufficient_source_text",',
      '"study_disease":"unknown","disease_match":"unknown","scope":"unclear",',
      '"evidence_types":["other"],"biological_processes":["other"],',
      '"effect":"unclear","observation":"No usable observation was available.",',
      '"study_context":"unknown","abstract_quote":"",',
      '"relevance_reason":"The supplied abstract was insufficient."}]}]}')
  }
  json_instruction <- paste(
    "Return one valid json object only, without Markdown fences or commentary.",
    "Match this JSON schema:",
    .InterMetaDeepSeekJSON(.InterMetaDeepSeekSchemaArrays(schema)),
    "JSON output example (replace the placeholder with a supplied pathway;",
    "never invent evidence or citations):", example,
    sep = "\n"
  )
  system_index <- which(vapply(out, function(x) identical(x$role, "system"),
                               logical(1)))
  if (length(system_index)) {
    index <- system_index[1L]
    out[[index]]$content <- paste(out[[index]]$content, json_instruction, sep = "\n\n")
  } else {
    out <- c(list(list(role = "system", content = json_instruction)), out)
  }
  unname(out)
}

.InterMetaDeepSeekContainsSecret <- function(x, api_key) {
  if (is.character(x)) {
    return(any(grepl(api_key, x[!is.na(x)], fixed = TRUE)))
  }
  if (is.list(x)) {
    if (!is.null(names(x)) &&
        .InterMetaDeepSeekContainsSecret(names(x), api_key)) return(TRUE)
    return(any(vapply(x, .InterMetaDeepSeekContainsSecret, logical(1),
                      api_key = api_key)))
  }
  FALSE
}

# This small transport boundary is mockable for entirely offline testing.
# No redirects, retries, request logging, or fallback model calls are enabled.
.InterMetaDeepSeekRequest <- function(endpoint, body, api_key, timeout) {
  if (!identical(endpoint, .InterMetaDeepSeekEndpoint)) {
    stop("DeepSeek requests must use the fixed official endpoint.", call. = FALSE)
  }
  .InterMetaDeepSeekRequire("httr2")
  tryCatch(suppressWarnings({
    request <- httr2::request(endpoint)
    request <- httr2::req_method(request, "POST")
    request <- httr2::req_headers(
      request, Authorization = paste("Bearer", api_key), Accept = "application/json"
    )
    request <- httr2::req_body_raw(
      request, .InterMetaDeepSeekJSON(body), type = "application/json"
    )
    request <- httr2::req_timeout(request, seconds = timeout)
    request <- httr2::req_options(
      request, followlocation = FALSE, maxredirs = 0L, verbose = FALSE
    )
    request <- httr2::req_retry(
      request, max_tries = 1L, is_transient = function(response) FALSE
    )
    request <- httr2::req_error(request, is_error = function(response) FALSE)
    response <- httr2::req_perform(request)
    list(status_code = httr2::resp_status(response),
         body = httr2::resp_body_string(response))
  }), error = function(e) {
    stop(paste(
      "DeepSeek request failed (network, timeout, or transport error);",
      "no automatic retry was attempted."
    ), call. = FALSE)
  })
}

.InterMetaDeepSeekHTTPError <- function(status) {
  detail <- switch(as.character(status),
    `400` = "invalid request", `401` = "authentication failed",
    `402` = "insufficient account balance", `403` = "access denied",
    `404` = "requested model or endpoint unavailable", `408` = "request timeout",
    `422` = "invalid request parameters", `429` = "rate limit reached",
    if (status >= 300 && status < 400) "redirect refused" else
      if (status >= 500) "service unavailable" else "unsuccessful response"
  )
  stop("DeepSeek HTTP ", status, " (", detail,
       "); no automatic retry was attempted.", call. = FALSE)
}

.InterMetaDeepSeekDecode <- function(text, message) {
  tryCatch(
    # parse_json accepts text only: fromJSON also treats strings as URLs or
    # local paths, which an untrusted completion must never be able to open.
    suppressWarnings(jsonlite::parse_json(text, simplifyVector = FALSE)),
    error = function(e) stop(message, call. = FALSE)
  )
}

.InterMetaDeepSeekSafeFinishReason <- function(x) {
  allowed <- c("stop", "length", "content_filter", "tool_calls",
               "insufficient_system_resource")
  if (is.character(x) && length(x) == 1L && !is.na(x) && x %in% allowed) {
    return(x)
  }
  "unknown"
}

.InterMetaDeepSeekFinishReasonMessage <- function(reason) {
  guidance <- switch(reason,
    length = paste(
      "The output token budget was reached; increase max_tokens or request a",
      "more concise report before making a manual cached-evidence retry."
    ),
    content_filter = "The provider filtered the completion; do not retry it as a length problem.",
    tool_calls = "An unexpected tool call was returned even though this report does not use tools.",
    insufficient_system_resource = paste(
      "The provider interrupted generation because inference resources were insufficient;",
      "a later manual cached-evidence retry may succeed."
    ),
    "The provider returned an unrecognized or missing finish reason."
  )
  paste0(
    "DeepSeek completion was not finished normally (finish_reason=", reason,
    "); truncated or refused output is rejected. ", guidance
  )
}

.InterMetaDeepSeekMetadata <- function(response, body, request_json,
                                      response_json, finish_reason) {
  out <- list(
    endpoint = .InterMetaDeepSeekEndpoint, requested_model = body$model,
    finish_reason = finish_reason, thinking = body$thinking$type,
    max_tokens = body$max_tokens,
    request_sha256 = digest::digest(request_json, algo = "sha256", serialize = FALSE),
    response_sha256 = digest::digest(response_json, algo = "sha256", serialize = FALSE)
  )
  if (!is.null(body$reasoning_effort)) out$reasoning_effort <- body$reasoning_effort
  for (field in c("model", "id")) {
    value <- response[[field]]
    if (is.character(value) && length(value) == 1L && !is.na(value) &&
        nzchar(value) && nchar(value) <= 512L &&
        grepl("^[A-Za-z0-9][A-Za-z0-9._:/-]*$", value) &&
        !grepl("sk-", value, fixed = TRUE)) {
      out[[if (field == "model") "response_model" else "response_id"]] <- value
    }
  }
  scalar_number <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
  }
  if (scalar_number(response$created)) out$created <- response$created
  if (is.list(response$usage)) {
    fields <- c("prompt_tokens", "completion_tokens", "total_tokens",
                "prompt_cache_hit_tokens", "prompt_cache_miss_tokens")
    usage <- response$usage[intersect(fields, names(response$usage))]
    usage <- Filter(scalar_number, usage)
    if (length(usage)) out$usage <- usage
  }
  out
}

.InterMetaDeepSeekGenerate <- function(messages, schema, model, api_key_env,
                                       max_tokens, timeout, thinking,
                                       reasoning_effort) {
  # Resolve credentials in this invocation frame, never in the provider closure.
  api_key <- Sys.getenv(api_key_env, unset = "")
  if (!nzchar(trimws(api_key))) {
    stop("DeepSeek API key is missing; set the configured API-key environment variable.",
         call. = FALSE)
  }
  if (grepl("[[:space:][:cntrl:]]", api_key)) {
    stop("DeepSeek API key must be a single token without whitespace.", call. = FALSE)
  }
  .InterMetaDeepSeekRequire("jsonlite")
  .InterMetaDeepSeekRequire("httr2")
  body <- list(
    model = model, messages = .InterMetaDeepSeekMessages(messages, schema),
    response_format = list(type = "json_object"), stream = FALSE,
    max_tokens = max_tokens, thinking = list(type = thinking)
  )
  if (identical(thinking, "enabled")) {
    body$reasoning_effort <- reasoning_effort
  } else {
    body$temperature <- 0
  }
  request_json <- .InterMetaDeepSeekJSON(body)
  if (.InterMetaDeepSeekContainsSecret(body, api_key) ||
      .InterMetaDeepSeekContainsSecret(request_json, api_key)) {
    stop("DeepSeek request was rejected because it contained credential material.",
         call. = FALSE)
  }
  result <- tryCatch(
    suppressWarnings(.InterMetaDeepSeekRequest(
      .InterMetaDeepSeekEndpoint, body, api_key, timeout
    )),
    error = function(e) stop(paste(
      "DeepSeek request failed (network, timeout, or transport error);",
      "no automatic retry was attempted."
    ), call. = FALSE)
  )
  if (!is.list(result) || !is.numeric(result$status_code) ||
      length(result$status_code) != 1L || !is.finite(result$status_code) ||
      result$status_code != floor(result$status_code) ||
      result$status_code < 100 || result$status_code > 599) {
    stop("DeepSeek returned an invalid HTTP response.", call. = FALSE)
  }
  if (result$status_code < 200 || result$status_code >= 300) {
    .InterMetaDeepSeekHTTPError(result$status_code)
  }
  response_json <- result$body
  if (!is.character(response_json) || length(response_json) != 1L ||
      is.na(response_json) || !nzchar(trimws(response_json))) {
    stop("DeepSeek returned an empty response body.", call. = FALSE)
  }
  if (.InterMetaDeepSeekContainsSecret(response_json, api_key)) {
    stop("DeepSeek response was rejected because it contained credential material.",
         call. = FALSE)
  }
  response <- .InterMetaDeepSeekDecode(
    response_json, "DeepSeek returned an invalid JSON response envelope."
  )
  if (.InterMetaDeepSeekContainsSecret(response, api_key)) {
    stop("DeepSeek response was rejected because it contained credential material.",
         call. = FALSE)
  }
  if (!is.list(response) || !is.list(response$choices) ||
      !length(response$choices) || !is.list(response$choices[[1L]])) {
    stop("DeepSeek response is missing a completion choice.", call. = FALSE)
  }
  choice <- response$choices[[1L]]
  if (!identical(choice$finish_reason, "stop")) {
    reason <- .InterMetaDeepSeekSafeFinishReason(choice$finish_reason)
    stop(.InterMetaDeepSeekFinishReasonMessage(reason), call. = FALSE)
  }
  content <- if (is.list(choice$message)) choice$message$content else NULL
  if (!is.character(content) || length(content) != 1L || is.na(content) ||
      !nzchar(trimws(content))) {
    stop("DeepSeek returned empty completion content.", call. = FALSE)
  }
  parsed <- .InterMetaDeepSeekDecode(
    content, "DeepSeek completion content is not valid JSON."
  )
  if (.InterMetaDeepSeekContainsSecret(parsed, api_key)) {
    stop("DeepSeek response was rejected because it contained credential material.",
         call. = FALSE)
  }
  if (!is.list(parsed) || !"pathways" %in% names(parsed) ||
      !is.list(parsed$pathways)) {
    stop("DeepSeek completion must contain a pathways collection.", call. = FALSE)
  }
  # Semantic checks (including pathway-specific PMID membership and enums)
  # remain the responsibility of InterMetaAIInterpret(), not the HTTP adapter.
  attr(parsed, "intermeta_provider_metadata") <- .InterMetaDeepSeekMetadata(
    response, body, request_json, response_json, choice$finish_reason
  )
  parsed
}

#' Create the official DeepSeek provider for InterMeta AI interpretation
#'
#' Construction is offline. Credentials are resolved from an environment
#' variable only when generate() is called. Requests use JSON-object output,
#' the fixed official endpoint, and no automatic retries or model fallback.
InterMetaDeepSeekProvider <- function(
    model = "deepseek-v4-pro", api_key_env = "DEEPSEEK_API_KEY",
    max_tokens = 16384L, timeout = 180,
    thinking = c("enabled", "disabled"),
    reasoning_effort = c("high", "low", "max")) {
  if (!is.character(model) || length(model) != 1L || is.na(model) ||
      nchar(model) > 128L || !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]*$", model)) {
    stop("model must be one non-empty model identifier of at most 128 characters.",
         call. = FALSE)
  }
  if (!is.character(api_key_env) || length(api_key_env) != 1L ||
      is.na(api_key_env) || !grepl("^[A-Za-z_][A-Za-z0-9_]*$", api_key_env)) {
    stop("api_key_env must be one valid environment-variable name.", call. = FALSE)
  }
  if (!is.numeric(max_tokens) || length(max_tokens) != 1L ||
      !is.finite(max_tokens) || max_tokens != floor(max_tokens) ||
      max_tokens < 1 || max_tokens > .Machine$integer.max) {
    stop("max_tokens must be one positive integer within R's integer range.",
         call. = FALSE)
  }
  if (!is.numeric(timeout) || length(timeout) != 1L ||
      !is.finite(timeout) || timeout <= 0) {
    stop("timeout must be one finite positive number of seconds.", call. = FALSE)
  }
  if (missing(thinking)) thinking <- "enabled"
  if (!is.character(thinking) || length(thinking) != 1L || is.na(thinking) ||
      !thinking %in% c("enabled", "disabled")) {
    stop("thinking must be exactly 'enabled' or 'disabled'.", call. = FALSE)
  }
  if (missing(reasoning_effort)) reasoning_effort <- "high"
  if (!is.character(reasoning_effort) || length(reasoning_effort) != 1L ||
      is.na(reasoning_effort) || !reasoning_effort %in% c("high", "low", "max")) {
    stop("reasoning_effort must be exactly 'high', 'low', or 'max'.", call. = FALSE)
  }
  max_tokens <- as.integer(max_tokens)
  # Every constructor argument has been forced and validated before this
  # closure is created. Caller environments and credentials are not captured.
  generate <- function(messages, schema) {
    .InterMetaDeepSeekGenerate(
      messages, schema, model, api_key_env, max_tokens, timeout, thinking,
      reasoning_effort
    )
  }
  InterMetaAIProvider(generate = generate, name = "DeepSeek", model = model)
}
