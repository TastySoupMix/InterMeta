# Optional fixed-endpoint adapters. Construction/catalog access are offline;
# credentials are resolved only in a generate() invocation, never its closure.

#' List example provider and model presets without contacting a service
InterMetaAIModels <- function() {
  data.frame(
    provider = c("deepseek", "deepseek", "openai", "openai", "anthropic",
                 "gemini", "qwen", "moonshot"),
    model = c("deepseek-v4-pro", "deepseek-v4-flash", "gpt-4.1-mini", "gpt-4.1",
              "claude-sonnet-4-5", "gemini-2.5-flash", "qwen-plus", "kimi-k2.6"),
    api_key_env = c("DEEPSEEK_API_KEY", "DEEPSEEK_API_KEY", "OPENAI_API_KEY",
                    "OPENAI_API_KEY", "ANTHROPIC_API_KEY", "GEMINI_API_KEY",
                    "DASHSCOPE_API_KEY", "MOONSHOT_API_KEY"),
    endpoint = c(rep("https://api.deepseek.com/chat/completions", 2L),
                 rep("https://api.openai.com/v1/responses", 2L),
                 "https://api.anthropic.com/v1/messages",
                 "https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent",
                 "https://dashscope.aliyuncs.com/compatible-mode/v1/chat/completions",
                 "https://api.moonshot.ai/v1/chat/completions"),
    protocol = c("chat_completions", "chat_completions", "responses", "responses",
                 "messages", "generate_content", "chat_completions", "chat_completions"),
    default = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
    documentation_url = c(rep("https://api-docs.deepseek.com/quick_start/pricing", 2L),
      "https://developers.openai.com/api/docs/models/gpt-4.1-mini",
      "https://developers.openai.com/api/docs/models/gpt-4.1",
      "https://platform.claude.com/docs/en/build-with-claude/structured-outputs",
      "https://ai.google.dev/gemini-api/docs/models/gemini-2.5-flash",
      "https://www.alibabacloud.com/help/en/model-studio/compatibility-of-openai-with-dashscope",
      "https://platform.kimi.ai/docs/api/models-overview"),
    checked_on = rep("2026-09-04", 8L),
    note = c(rep("Example model, not a latest-model or availability guarantee.", 6L),
      paste("Beijing legacy compatibility endpoint; use a Beijing-region key.",
            "Workspace hosts require a custom InterMetaAIProvider adapter.",
            "Example model; availability and endpoint policy can change."),
      "Example model on the international Moonshot endpoint; availability can change."),
    stringsAsFactors = FALSE
  )
}

.InterMetaLLMProviders <- c("deepseek", "openai", "anthropic", "gemini", "qwen", "moonshot")

.InterMetaLLMValidateModel <- function(model) {
  if (!is.character(model) || length(model) != 1L || is.na(model) ||
      nchar(model) > 128L || !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]*$", model)) {
    stop("model must be one non-empty model identifier of at most 128 characters.",
         call. = FALSE)
  }
  invisible(model)
}

.InterMetaLLMLabel <- function(provider) {
  labels <- c(deepseek = "DeepSeek", openai = "OpenAI", anthropic = "Anthropic",
              gemini = "Gemini", qwen = "Qwen", moonshot = "Moonshot")
  unname(labels[[provider]])
}

.InterMetaLLMEndpoint <- function(provider, model) {
  if (!is.character(provider) || length(provider) != 1L || is.na(provider) ||
      !provider %in% .InterMetaLLMProviders) {
    stop("provider must be one supported canonical provider name.", call. = FALSE)
  }
  .InterMetaLLMValidateModel(model)
  catalog <- InterMetaAIModels()
  endpoint <- catalog$endpoint[match(provider, catalog$provider)]
  sub("{model}", model, endpoint, fixed = TRUE)
}

.InterMetaLLMRequire <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The optional package '", package,
         "' is required to generate an LLM interpretation.", call. = FALSE)
  }
  invisible(NULL)
}

.InterMetaLLMRequest <- function(provider, model, body, api_key, timeout) {
  endpoint <- .InterMetaLLMEndpoint(provider, model)
  if (provider == "deepseek") {
    stop("Use the dedicated DeepSeek adapter for this provider.", call. = FALSE)
  }
  .InterMetaLLMRequire("httr2")
  .InterMetaLLMRequire("jsonlite")
  tryCatch(suppressWarnings({
    request <- httr2::request(endpoint)
    request <- httr2::req_method(request, "POST")
    headers <- switch(provider,
      anthropic = list(`x-api-key` = api_key, `anthropic-version` = "2023-06-01"),
      gemini = list(`x-goog-api-key` = api_key),
      list(Authorization = paste("Bearer", api_key)))
    request <- do.call(httr2::req_headers,
      c(list(request), headers, list(Accept = "application/json")))
    request <- httr2::req_body_raw(request, .InterMetaDeepSeekJSON(body),
                                  type = "application/json")
    request <- httr2::req_timeout(request, seconds = timeout)
    request <- httr2::req_options(request, followlocation = FALSE,
                                  maxredirs = 0L, verbose = FALSE)
    request <- httr2::req_retry(request, max_tries = 1L,
                                is_transient = function(response) FALSE)
    request <- httr2::req_error(request, is_error = function(response) FALSE)
    response <- httr2::req_perform(request)
    list(status_code = httr2::resp_status(response),
         body = httr2::resp_body_string(response))
  }), error = function(e) {
    stop("LLM request failed (network, timeout, or transport error); no automatic retry was attempted.",
         call. = FALSE)
  })
}

.InterMetaLLMSchema <- function(schema, provider) {
  out <- .InterMetaDeepSeekSchemaArrays(schema)
  out$name <- NULL
  if (provider %in% c("anthropic", "gemini")) {
    # Native grammar subsets do not support every numeric/string constraint.
    # Full constraints remain in the system prompt and InterMeta validation.
    strip <- function(x) {
      if (!is.list(x)) return(x)
      if (!is.null(names(x))) {
        x <- x[!names(x) %in% c("minLength", "maxLength", "minimum", "maximum",
                               "minItems", "maxItems", "pattern", "format")]
      }
      lapply(x, strip)
    }
    out <- strip(out)
  }
  out
}

.InterMetaLLMBody <- function(provider, model, messages, schema, max_tokens) {
  messages <- .InterMetaDeepSeekMessages(messages, schema)
  if (provider == "openai") {
    return(list(model = model, input = messages, stream = FALSE, store = FALSE,
      max_output_tokens = max_tokens, truncation = "disabled",
      text = list(format = list(type = "json_schema", name = schema$name,
        strict = TRUE, schema = .InterMetaLLMSchema(schema, provider)))))
  }
  system_rows <- vapply(messages, function(x) identical(x$role, "system"), logical(1))
  system_text <- paste(vapply(messages[system_rows], `[[`, character(1), "content"),
                       collapse = "\n\n")
  if (provider == "anthropic") {
    return(list(model = model, system = system_text, messages = messages[!system_rows],
      max_tokens = max_tokens, stream = FALSE,
      output_config = list(format = list(type = "json_schema",
        schema = .InterMetaLLMSchema(schema, provider)))))
  }
  if (provider == "gemini") {
    contents <- lapply(messages[!system_rows], function(x) {
      list(role = if (x$role == "assistant") "model" else "user",
           parts = list(list(text = x$content)))
    })
    return(list(systemInstruction = list(parts = list(list(text = system_text))),
      contents = contents, generationConfig = list(
        maxOutputTokens = max_tokens, candidateCount = 1L,
        responseMimeType = "application/json",
        responseJsonSchema = .InterMetaLLMSchema(schema, provider))))
  }
  body <- list(model = model, messages = messages, stream = FALSE,
    max_tokens = max_tokens, response_format = list(type = "json_object"))
  if (provider == "qwen") body$enable_thinking <- FALSE
  # This documented preset supports non-thinking JSON. Unknown overrides keep
  # provider defaults; no extra paid request or model fallback is attempted.
  if (provider == "moonshot" && identical(model, "kimi-k2.6")) {
    body$thinking <- list(type = "disabled")
  }
  body
}

.InterMetaLLMContent <- function(response, provider) {
  fail <- function() stop(
    "LLM completion was not finished normally; malformed, truncated, refused, or nonterminal output is rejected.",
    call. = FALSE)
  scalar_text <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  if (!is.list(response) || is.null(names(response)) || !is.null(response$error)) fail()
  if (provider == "openai") {
    if (!identical(response$status, "completed") || !is.null(response$incomplete_details) ||
        !is.list(response$output) || !length(response$output)) fail()
    texts <- character()
    for (item in response$output) {
      if (!is.list(item)) fail()
      if (identical(item$type, "reasoning")) next
      if (!identical(item$type, "message") || !identical(item$role, "assistant") ||
          !identical(item$status, "completed") || !is.list(item$content) ||
          !length(item$content)) fail()
      for (part in item$content) {
        if (!is.list(part) || !identical(part$type, "output_text") ||
            !scalar_text(part$text)) fail()
        texts <- c(texts, part$text)
      }
    }
    if (!length(texts)) fail()
    return(list(text = paste(texts, collapse = ""), finish_reason = "completed"))
  }
  if (provider == "anthropic") {
    if (!identical(response$stop_reason, "end_turn") ||
        !identical(response$type, "message") || !identical(response$role, "assistant") ||
        !is.null(response$stop_details) || !is.list(response$content) ||
        !length(response$content)) fail()
    texts <- character()
    for (part in response$content) {
      if (!is.list(part)) fail()
      if (part$type %in% c("thinking", "redacted_thinking")) next
      if (!identical(part$type, "text") || !scalar_text(part$text)) fail()
      texts <- c(texts, part$text)
    }
    if (!length(texts)) fail()
    return(list(text = paste(texts, collapse = ""), finish_reason = "end_turn"))
  }
  if (provider == "gemini") {
    if (!is.null(response$promptFeedback$blockReason) || !is.list(response$candidates) ||
        length(response$candidates) != 1L) fail()
    candidate <- response$candidates[[1L]]
    if (!is.list(candidate) || !identical(candidate$finishReason, "STOP") ||
        !is.list(candidate$content) || !identical(candidate$content$role, "model") ||
        !is.list(candidate$content$parts) || !length(candidate$content$parts)) fail()
    texts <- character()
    for (part in candidate$content$parts) {
      if (!is.list(part)) fail()
      if (isTRUE(part$thought)) next
      if (!scalar_text(part$text) || any(c("functionCall", "functionResponse",
          "inlineData", "fileData", "executableCode", "codeExecutionResult") %in% names(part))) fail()
      texts <- c(texts, part$text)
    }
    if (!length(texts)) fail()
    return(list(text = paste(texts, collapse = ""), finish_reason = "STOP"))
  }
  if (!is.list(response$choices) || length(response$choices) != 1L) fail()
  choice <- response$choices[[1L]]
  if (!is.list(choice) || !identical(choice$finish_reason, "stop") ||
      !is.list(choice$message) || !identical(choice$message$role, "assistant") ||
      !is.null(choice$message$refusal) || length(choice$message$tool_calls) ||
      !is.null(choice$message$function_call) || !scalar_text(choice$message$content)) fail()
  list(text = choice$message$content, finish_reason = "stop")
}

.InterMetaLLMMetadata <- function(response, provider, model, body, max_tokens,
                                 request_json, response_json, finish_reason) {
  out <- list(endpoint = .InterMetaLLMEndpoint(provider, model), requested_model = model,
    finish_reason = finish_reason, max_tokens = max_tokens,
    request_sha256 = digest::digest(request_json, algo = "sha256", serialize = FALSE),
    response_sha256 = digest::digest(response_json, algo = "sha256", serialize = FALSE))
  model_field <- if (provider == "gemini") "modelVersion" else "model"
  id_field <- if (provider == "gemini") "responseId" else "id"
  for (field in c(model_field, id_field)) {
    value <- response[[field]]
    if (is.character(value) && length(value) == 1L && !is.na(value) &&
        nchar(value) <= 512L && grepl("^[A-Za-z0-9][A-Za-z0-9._:/-]*$", value) &&
        !grepl("sk-", value, fixed = TRUE)) {
      out[[if (field == model_field) "response_model" else "response_id"]] <- value
    }
  }
  scalar_number <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
  created <- if (provider == "openai") response$created_at else response$created
  if (scalar_number(created)) out$created <- created
  usage <- if (provider == "gemini") response$usageMetadata else response$usage
  if (is.list(usage)) {
    mapping <- switch(provider,
      openai = c(prompt_tokens = "input_tokens", completion_tokens = "output_tokens",
                 total_tokens = "total_tokens"),
      anthropic = c(prompt_tokens = "input_tokens", completion_tokens = "output_tokens"),
      gemini = c(prompt_tokens = "promptTokenCount", completion_tokens = "candidatesTokenCount",
                 total_tokens = "totalTokenCount"),
      c(prompt_tokens = "prompt_tokens", completion_tokens = "completion_tokens",
        total_tokens = "total_tokens", prompt_cache_hit_tokens = "prompt_cache_hit_tokens",
        prompt_cache_miss_tokens = "prompt_cache_miss_tokens"))
    normalized <- list()
    for (name in names(mapping)) {
      value <- usage[[mapping[[name]]]]
      if (scalar_number(value)) normalized[[name]] <- value
    }
    if (length(normalized)) out$usage <- normalized
  }
  if (provider == "qwen") out$thinking <- "disabled"
  if (provider == "moonshot" && !is.null(body$thinking)) out$thinking <- body$thinking$type
  .InterMetaAIProviderCallMetadata(out)
}

.InterMetaLLMGenerate <- function(messages, schema, provider, model, api_key_env,
                                 max_tokens, timeout) {
  api_key <- Sys.getenv(api_key_env, unset = "")
  if (!nzchar(trimws(api_key))) {
    stop(.InterMetaLLMLabel(provider), " API key is missing; set the configured API-key environment variable.",
         call. = FALSE)
  }
  if (grepl("[[:space:][:cntrl:]]", api_key)) {
    stop("LLM API key must be a single token without whitespace.", call. = FALSE)
  }
  .InterMetaLLMRequire("jsonlite")
  .InterMetaLLMRequire("httr2")
  # Shared serializer/message helpers never see a credential argument.
  body <- tryCatch(.InterMetaLLMBody(provider, model, messages, schema, max_tokens),
    error = function(e) stop("LLM messages or schema are invalid and could not be encoded.", call. = FALSE))
  request_json <- .InterMetaDeepSeekJSON(body)
  if (.InterMetaDeepSeekContainsSecret(body, api_key) ||
      .InterMetaDeepSeekContainsSecret(request_json, api_key)) {
    stop("LLM request was rejected because it contained credential material.", call. = FALSE)
  }
  result <- tryCatch(suppressWarnings(.InterMetaLLMRequest(provider, model, body, api_key, timeout)),
    error = function(e) stop(
      "LLM request failed (network, timeout, or transport error); no automatic retry was attempted.",
      call. = FALSE))
  if (!is.list(result) || !is.numeric(result$status_code) || length(result$status_code) != 1L ||
      !is.finite(result$status_code) || result$status_code != floor(result$status_code) ||
      result$status_code < 100 || result$status_code > 599) {
    stop("LLM provider returned an invalid HTTP response.", call. = FALSE)
  }
  if (result$status_code < 200 || result$status_code >= 300) {
    stop(.InterMetaLLMLabel(provider), " HTTP ", result$status_code,
         "; no automatic retry or redirect was attempted.", call. = FALSE)
  }
  response_json <- result$body
  if (!is.character(response_json) || length(response_json) != 1L || is.na(response_json) ||
      !nzchar(trimws(response_json))) stop("LLM provider returned an empty response body.", call. = FALSE)
  reject_secret <- function(x) {
    if (.InterMetaDeepSeekContainsSecret(x, api_key)) {
      stop("LLM response was rejected because it contained credential material.", call. = FALSE)
    }
  }
  reject_secret(response_json)
  response <- .InterMetaDeepSeekDecode(response_json, "LLM provider returned an invalid JSON response envelope.")
  reject_secret(response)
  content <- tryCatch(.InterMetaLLMContent(response, provider),
    error = function(e) stop(
      "LLM completion was not finished normally; malformed, truncated, refused, or nonterminal output is rejected.",
      call. = FALSE))
  parsed <- .InterMetaDeepSeekDecode(content$text, "LLM completion content is not valid JSON.")
  reject_secret(parsed)
  if (!is.list(parsed) || !"pathways" %in% names(parsed) || !is.list(parsed$pathways)) {
    stop("LLM completion must contain a pathways collection.", call. = FALSE)
  }
  attr(parsed, "intermeta_provider_metadata") <- .InterMetaLLMMetadata(
    response, provider, model, body, max_tokens, request_json, response_json, content$finish_reason)
  parsed
}

#' Create an optional fixed-endpoint LLM provider without reading credentials
InterMetaLLMProvider <- function(
    provider = c("deepseek", "openai", "anthropic", "gemini", "qwen", "moonshot"),
    model = NULL, api_key_env = NULL, max_tokens = 16384L, timeout = 180) {
  if (missing(provider)) provider <- "deepseek"
  if (!is.character(provider) || length(provider) != 1L || is.na(provider) ||
      !provider %in% .InterMetaLLMProviders) {
    stop("provider must be one supported canonical provider name.", call. = FALSE)
  }
  catalog <- InterMetaAIModels()
  preset <- catalog[catalog$provider == provider & catalog$default, , drop = FALSE]
  if (is.null(model)) model <- preset$model[1L]
  .InterMetaLLMValidateModel(model)
  if (is.null(api_key_env)) api_key_env <- preset$api_key_env[1L]
  if (!is.character(api_key_env) || length(api_key_env) != 1L || is.na(api_key_env) ||
      !grepl("^[A-Za-z_][A-Za-z0-9_]*$", api_key_env)) {
    stop("api_key_env must be one valid environment-variable name.", call. = FALSE)
  }
  if (!is.numeric(max_tokens) || length(max_tokens) != 1L || !is.finite(max_tokens) ||
      max_tokens != floor(max_tokens) || max_tokens < 1 || max_tokens > .Machine$integer.max) {
    stop("max_tokens must be one positive integer within R's integer range.", call. = FALSE)
  }
  if (!is.numeric(timeout) || length(timeout) != 1L || !is.finite(timeout) || timeout <= 0) {
    stop("timeout must be one finite positive number of seconds.", call. = FALSE)
  }
  max_tokens <- as.integer(max_tokens)
  if (provider == "deepseek") {
    return(InterMetaDeepSeekProvider(model = model, api_key_env = api_key_env,
                                    max_tokens = max_tokens, timeout = timeout))
  }
  generate <- function(messages, schema) {
    .InterMetaLLMGenerate(messages, schema, provider, model, api_key_env, max_tokens, timeout)
  }
  InterMetaAIProvider(generate, name = .InterMetaLLMLabel(provider), model = model)
}
