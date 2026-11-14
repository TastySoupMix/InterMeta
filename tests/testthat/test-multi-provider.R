multi_fn <- function(name) getFromNamespace(name, "InterMeta")
multi_providers <- c("openai", "anthropic", "gemini", "qwen", "moonshot")
multi_env <- "INTERMETA_OFFLINE_MULTI_KEY"
multi_key <- "offline-multi-key-NOT-A-REAL-CREDENTIAL"

multi_setup <- function(.env = parent.frame()) {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("httr2")
  withr::local_envvar(stats::setNames(multi_key, multi_env), .local_envir = .env)
}

multi_messages <- function() list(
  list(role = "system", content = "Interpret only the supplied aggregate evidence."),
  list(role = "user", content = list(context = list(pathways = data.frame(
    Pathway = "P1", Direction = "Unfavorable", stringsAsFactors = FALSE))))
)

multi_content <- function() '{"pathways":[{"pathway":"P1","literature_consistency":"insufficient_evidence","summary":"Insufficient evidence.","claims":[],"assessments":[]}]}'

multi_envelope <- function(provider, content = multi_content()) {
  switch(provider,
    openai = list(id = "resp_offline", model = "gpt-4.1-mini", status = "completed",
      created_at = 123456, output = list(list(type = "reasoning", summary = list()),
        list(type = "message", role = "assistant", status = "completed",
          content = list(list(type = "output_text", text = content)))),
      usage = list(input_tokens = 12L, output_tokens = 34L, total_tokens = 46L)),
    anthropic = list(id = "msg_offline", model = "claude-sonnet-4-5", type = "message",
      role = "assistant", stop_reason = "end_turn",
      content = list(list(type = "thinking", thinking = "PRIVATE_REASONING"),
        list(type = "text", text = content)),
      usage = list(input_tokens = 12L, output_tokens = 34L)),
    gemini = list(responseId = "gemini-offline", modelVersion = "gemini-2.5-flash",
      candidates = list(list(finishReason = "STOP", content = list(role = "model",
        parts = list(list(thought = TRUE, text = "PRIVATE_REASONING"), list(text = content))))),
      usageMetadata = list(promptTokenCount = 12L, candidatesTokenCount = 34L, totalTokenCount = 46L)),
    list(id = "chat-offline", model = if (provider == "qwen") "qwen-plus" else "kimi-k2.6",
      choices = list(list(finish_reason = "stop", message = list(role = "assistant",
        content = content, reasoning_content = "PRIVATE_REASONING"))),
      usage = list(prompt_tokens = 12L, completion_tokens = 34L, total_tokens = 46L)))
}

multi_response <- function(provider, envelope = multi_envelope(provider), status = 200L) {
  list(status_code = status,
       body = as.character(jsonlite::toJSON(envelope, auto_unbox = TRUE, null = "null")))
}

multi_provider <- function(provider, ...) multi_fn("InterMetaLLMProvider")(
  provider = provider, api_key_env = multi_env, ...)

multi_generate <- function(provider, object = multi_provider(provider)) {
  object$generate(multi_messages(), multi_fn(".InterMetaAIInterpretationSchema")())
}

test_that("offline model catalog has documented examples and one default per provider", {
  catalog <- multi_fn("InterMetaAIModels")()
  expect_true(all(c("provider", "model", "api_key_env", "endpoint", "protocol",
                    "default", "documentation_url", "checked_on", "note") %in% names(catalog)))
  expect_setequal(catalog$provider, c("deepseek", multi_providers))
  expect_true(all(tapply(catalog$default, catalog$provider, sum) == 1L))
  expect_true(all(catalog$checked_on == "2026-09-04"))
  expect_true(all(grepl("^https://", catalog$documentation_url)))
  expect_match(catalog$note[catalog$provider == "qwen"], "Beijing legacy")
  local_mocked_bindings(.InterMetaLLMRequest = function(...) stop("must not call"),
    .InterMetaDeepSeekRequest = function(...) stop("must not call"),
    .InterMetaLLMRequire = function(...) stop("must not load dependencies"), .package = "InterMeta")
  for (provider in c("deepseek", multi_providers)) {
    object <- multi_fn("InterMetaLLMProvider")(provider)
    expect_s3_class(object, "InterMetaAIProvider")
    expect_identical(names(object), c("generate", "name", "model"))
    expect_identical(object$model, catalog$model[catalog$provider == provider & catalog$default])
  }
  expect_identical(multi_fn("InterMetaLLMProvider")()$model, "deepseek-v4-pro")
  expect_identical(names(formals(multi_fn("InterMetaLLMProvider"))),
                   c("provider", "model", "api_key_env", "max_tokens", "timeout"))
})

test_that("generic constructor validates exact providers and safe parameter overrides", {
  constructor <- multi_fn("InterMetaLLMProvider")
  for (value in list("", "open", "OpenAI", "https://other.example", NA_character_,
                     c("openai", "qwen"), 1)) {
    expect_error(constructor(provider = value), "provider must")
  }
  for (value in list("", NA_character_, "../other", "model?key=value", "model#fragment",
                     "https://other.example", "model\nheader", c("a", "b"), 1,
                     paste(rep("m", 129L), collapse = ""))) {
    expect_error(constructor("openai", model = value), "model must")
  }
  for (value in list("", "1KEY", "A=B", "KEY\nHEADER", NA_character_, c("A", "B"))) {
    expect_error(constructor("gemini", api_key_env = value), "api_key_env must")
  }
  for (value in list(0, -1, Inf, NA_real_, TRUE, "100", c(1, 2), 1.5, 2^31)) {
    expect_error(constructor("anthropic", max_tokens = value), "max_tokens must")
  }
  for (value in list(0, -1, Inf, NA_real_, TRUE, "100", c(1, 2))) {
    expect_error(constructor("qwen", timeout = value), "timeout must")
  }
  expect_identical(constructor("openai", model = "gpt-4.1")$model, "gpt-4.1")
  expect_identical(constructor("gemini", model = "custom-model-identifier")$model,
                   "custom-model-identifier")
})

test_that("DeepSeek factory preserves the dedicated adapter defaults and overrides", {
  observed <- new.env(parent = emptyenv())
  local_mocked_bindings(InterMetaDeepSeekProvider = function(...) {
    observed$args <- list(...)
    structure(list(name = "DeepSeek", model = observed$args$model), class = "InterMetaAIProvider")
  }, .package = "InterMeta")
  multi_fn("InterMetaLLMProvider")()
  expect_identical(observed$args, list(model = "deepseek-v4-pro",
    api_key_env = "DEEPSEEK_API_KEY", max_tokens = 16384L, timeout = 180))
  multi_fn("InterMetaLLMProvider")("deepseek", model = "deepseek-v4-flash",
    api_key_env = multi_env, max_tokens = 42, timeout = 2.5)
  expect_identical(observed$args, list(model = "deepseek-v4-flash",
    api_key_env = multi_env, max_tokens = 42L, timeout = 2.5))
})

test_that("every native protocol serializes aggregate inputs and correct output controls", {
  multi_setup()
  observed <- new.env(parent = emptyenv())
  local_mocked_bindings(.InterMetaLLMRequest = function(provider, model, body, api_key, timeout) {
    observed$body <- body
    observed$model <- model
    expect_identical(api_key, multi_key)
    expect_identical(timeout, 180)
    multi_response(provider)
  }, .package = "InterMeta")
  for (provider in multi_providers) {
    result <- multi_generate(provider)
    expect_identical(result$pathways[[1L]]$pathway, "P1")
    body <- observed$body
    expect_false(grepl(multi_key, multi_fn(".InterMetaDeepSeekJSON")(body), fixed = TRUE))
    if (provider == "openai") {
      expect_identical(body$text$format$type, "json_schema")
      expect_true(body$text$format$strict)
      expect_false(body$store)
      expect_false(body$stream)
      expect_identical(body$truncation, "disabled")
      expect_identical(body$max_output_tokens, 16384L)
      expect_null(body$messages)
      expect_identical(body$input[[1L]]$role, "system")
      expect_identical(body$text$format$schema$required, list("pathways"))
      expect_false(body$text$format$schema$additionalProperties)
      expect_false(body$text$format$schema$properties$pathways$items$additionalProperties)
      expect_null(body$text$format$schema$name)
      expect_match(body$input[[2L]]$content, '"Pathway":"P1"', fixed = TRUE)
    } else if (provider == "anthropic") {
      expect_match(body$system, "valid json object")
      expect_identical(body$output_config$format$type, "json_schema")
      expect_identical(body$max_tokens, 16384L)
      expect_false(any(vapply(body$messages, function(x) x$role == "system", logical(1))))
      expect_false(body$stream)
      expect_identical(body$output_config$format$schema$required, list("pathways"))
      expect_false(body$output_config$format$schema$additionalProperties)
      expect_false(grepl('"maxLength":', multi_fn(".InterMetaDeepSeekJSON")(
        body$output_config$format$schema), fixed = TRUE))
      expect_match(body$system, '"maxLength":1000', fixed = TRUE)
    } else if (provider == "gemini") {
      expect_match(body$systemInstruction$parts[[1L]]$text, "valid json object")
      expect_identical(body$generationConfig$responseMimeType, "application/json")
      expect_identical(body$generationConfig$maxOutputTokens, 16384L)
      expect_identical(body$generationConfig$candidateCount, 1L)
      expect_identical(body$generationConfig$responseJsonSchema$required, list("pathways"))
      expect_false(body$generationConfig$responseJsonSchema$additionalProperties)
      expect_identical(body$contents[[1L]]$role, "user")
      expect_null(body$model)
      expect_null(body$messages)
    } else {
      expect_identical(body$response_format$type, "json_object")
      expect_false(body$stream)
      expect_identical(body$max_tokens, 16384L)
      expect_match(body$messages[[1L]]$content, "valid json object")
      if (provider == "qwen") expect_false(body$enable_thinking)
      if (provider == "moonshot") expect_identical(body$thinking$type, "disabled")
    }
    expect_null(body$temperature)
    meta <- attr(result, "intermeta_provider_metadata")
    expect_identical(meta$requested_model, observed$model)
    expect_match(meta$request_sha256, "^[a-f0-9]{64}$")
    expect_match(meta$response_sha256, "^[a-f0-9]{64}$")
    expect_equal(meta$usage$prompt_tokens, 12)
    expect_equal(meta$usage$completion_tokens, 34)
    if (provider == "anthropic") expect_null(meta$usage$total_tokens)
      else expect_equal(meta$usage$total_tokens, 46)
    expect_identical(multi_fn(".InterMetaAIProviderCallMetadata")(meta), meta)
    expect_false(grepl("PRIVATE_", rawToChar(serialize(result, NULL, ascii = TRUE)), fixed = TRUE))
  }
})

test_that("validated model override is sent unchanged without fallback or temperature assumptions", {
  multi_setup()
  observed <- new.env(parent = emptyenv())
  observed$calls <- 0L
  local_mocked_bindings(.InterMetaLLMRequest = function(provider, model, body, api_key, timeout) {
    observed$calls <- observed$calls + 1L
    expect_identical(model, "custom-safe-id-1")
    expect_identical(timeout, 12.5)
    expect_null(body$temperature)
    if (provider != "gemini") expect_identical(body$model, model)
    if (provider == "moonshot") expect_null(body$thinking)
    multi_response(provider)
  }, .package = "InterMeta")
  for (provider in multi_providers) {
    result <- multi_generate(provider, multi_provider(provider, model = "custom-safe-id-1",
                                                      max_tokens = 1024, timeout = 12.5))
    expect_identical(attr(result, "intermeta_provider_metadata")$requested_model, "custom-safe-id-1")
  }
  expect_identical(observed$calls, 5L)
})

test_that("provider credentials are lazy, rotatable, and absent from serialized closures", {
  multi_setup()
  observed <- new.env(parent = emptyenv())
  observed$keys <- character()
  local_mocked_bindings(.InterMetaLLMRequest = function(provider, model, body, api_key, timeout) {
    observed$keys <- c(observed$keys, api_key)
    multi_response(provider)
  }, .package = "InterMeta")
  for (provider in multi_providers) {
    Sys.unsetenv(multi_env)
    object <- multi_provider(provider)
    expect_error(multi_generate(provider, object), "API key is missing")
    for (key in c(multi_key, paste0(multi_key, "-ROTATED"))) {
      do.call(Sys.setenv, stats::setNames(list(key), multi_env))
      multi_generate(provider, object)
      expect_false(grepl(key, rawToChar(serialize(object, NULL, ascii = TRUE)), fixed = TRUE))
      expect_false("api_key" %in% ls(environment(object$generate), all.names = TRUE))
    }
    do.call(Sys.setenv, stats::setNames(list("bad token"), multi_env))
    expect_error(multi_generate(provider, object), "single token")
  }
  expect_identical(observed$keys, rep(c(multi_key, paste0(multi_key, "-ROTATED")), 5L))
})

test_that("HTTP transport uses fixed HTTPS hosts, proper auth headers and exactly one attempt", {
  multi_setup()
  observed <- new.env(parent = emptyenv())
  observed$calls <- 0L
  local_mocked_bindings(req_perform = function(req, ...) {
    observed$request <- req
    observed$calls <- observed$calls + 1L
    list(status = 200L, body = "{}")
  }, resp_status = function(x) x$status, resp_body_string = function(x) x$body,
  .package = "httr2")
  request <- multi_fn(".InterMetaLLMRequest")
  catalog <- multi_fn("InterMetaAIModels")()
  for (provider in multi_providers) {
    model <- catalog$model[catalog$provider == provider & catalog$default]
    result <- request(provider, model, list(model = model), multi_key, 21)
    expect_identical(result$status_code, 200L)
    req <- observed$request
    expect_identical(req$url, multi_fn(".InterMetaLLMEndpoint")(provider, model))
    expect_false(grepl(multi_key, req$url, fixed = TRUE))
    expect_identical(req$method, "POST")
    expect_identical(req$headers$Accept, "application/json")
    if (provider == "anthropic") {
      expect_identical(req$headers[["x-api-key"]], multi_key)
      expect_identical(req$headers[["anthropic-version"]], "2023-06-01")
      expect_null(req$headers$Authorization)
    } else if (provider == "gemini") {
      expect_identical(req$headers[["x-goog-api-key"]], multi_key)
      expect_null(req$headers$Authorization)
    } else expect_identical(req$headers$Authorization, paste("Bearer", multi_key))
    expect_false(req$options$followlocation)
    expect_equal(req$options$maxredirs, 0)
    expect_false(req$options$verbose)
    expect_equal(req$options$timeout_ms, 21000)
    expect_equal(req$policies$retry_max_tries, 1)
    expect_false(req$policies$retry_is_transient(list()))
  }
  expect_identical(observed$calls, 5L)
  expect_error(request("https://evil.example", "safe-model", list(), multi_key, 1), "provider must")
  expect_error(request("gemini", "safe-model?key=secret", list(), multi_key, 1), "model must")
  expect_identical(observed$calls, 5L)
})

test_that("HTTP and network errors never expose credential, body, or unsafe warnings", {
  multi_setup()
  observed <- new.env(parent = emptyenv())
  observed$status <- 401L
  observed$calls <- 0L
  local_mocked_bindings(.InterMetaLLMRequest = function(...) {
    observed$calls <- observed$calls + 1L
    list(status_code = observed$status, body = paste("PRIVATE_ERROR", multi_key))
  }, .package = "InterMeta")
  for (provider in multi_providers) for (status in c(301L, 400L, 401L, 403L, 429L, 500L)) {
    observed$status <- status
    error <- tryCatch(multi_generate(provider), error = identity)
    expect_match(conditionMessage(error), paste("HTTP", status), fixed = TRUE)
    expect_false(grepl("PRIVATE_|offline-multi-key", conditionMessage(error)))
    expect_null(conditionCall(error))
  }
  expect_identical(observed$calls, 30L)
  local_mocked_bindings(.InterMetaLLMRequest = function(...) {
    warning(paste("PRIVATE_WARNING", multi_key))
    stop(paste("PRIVATE_ERROR", multi_key))
  }, .package = "InterMeta")
  for (provider in multi_providers) {
    error <- tryCatch(multi_generate(provider), error = identity)
    expect_match(conditionMessage(error), "no automatic retry")
    expect_false(grepl("PRIVATE_|offline-multi-key", conditionMessage(error)))
    expect_null(conditionCall(error))
  }
})

test_that("nonterminal, truncated, refused and tool outputs fail closed for every protocol", {
  multi_setup()
  observed <- new.env(parent = emptyenv())
  local_mocked_bindings(.InterMetaLLMRequest = function(...) observed$response, .package = "InterMeta")
  for (provider in multi_providers) {
    mutations <- if (provider == "openai") list(
      function(x) { x$status <- "incomplete"; x },
      function(x) { x$status <- "queued"; x },
      function(x) { x$output[[2L]]$status <- "in_progress"; x },
      function(x) { x$output[[2L]]$content[[1L]]$type <- "refusal"; x },
      function(x) { x$output[[2L]]$type <- "function_call"; x },
      function(x) { x$incomplete_details <- list(reason = "max_output_tokens"); x }
    ) else if (provider == "anthropic") list(
      function(x) { x$stop_reason <- "max_tokens"; x },
      function(x) { x$stop_reason <- "refusal"; x },
      function(x) { x$stop_reason <- "pause_turn"; x },
      function(x) { x$content[[2L]]$type <- "tool_use"; x },
      function(x) { x$content[[2L]]$type <- NULL; x },
      function(x) { x$stop_details <- list(type = "refusal"); x }
    ) else if (provider == "gemini") list(
      function(x) { x$candidates[[1L]]$finishReason <- "MAX_TOKENS"; x },
      function(x) { x$candidates[[1L]]$finishReason <- "SAFETY"; x },
      function(x) { x$promptFeedback <- list(blockReason = "SAFETY"); x },
      function(x) { x$candidates[[1L]]$content$parts[[2L]]$functionCall <- list(name = "tool"); x }
    ) else list(
      function(x) { x$choices[[1L]]$finish_reason <- "length"; x },
      function(x) { x$choices[[1L]]$finish_reason <- "content_filter"; x },
      function(x) { x$choices[[1L]]$message$refusal <- "Refused"; x },
      function(x) { x$choices[[1L]]$message$tool_calls <- list(list(id = "tool")); x }
    )
    for (mutate in mutations) {
      observed$response <- multi_response(provider, mutate(multi_envelope(provider)))
      expect_error(multi_generate(provider), "not finished normally")
    }
    for (envelope in list(list(), list(error = list(message = "PRIVATE_ERROR")), NULL)) {
      observed$response <- multi_response(provider, envelope)
      expect_error(multi_generate(provider), "not finished normally")
    }
  }
})

test_that("malformed response bodies and completion JSON never trigger a fallback parser", {
  multi_setup()
  observed <- new.env(parent = emptyenv())
  local_mocked_bindings(.InterMetaLLMRequest = function(...) observed$response, .package = "InterMeta")
  for (provider in multi_providers) {
    for (content in c('{"pathways":', '```json\n{"pathways":[]}\n```',
                      'https://invalid.example/no-network', 'C:/not-a-real-file.json')) {
      observed$response <- multi_response(provider, multi_envelope(provider, content))
      expect_error(multi_generate(provider), "not valid JSON")
    }
    for (content in c('{}', '[]', '{"pathways":null}', '{"pathways":"wrong"}')) {
      observed$response <- multi_response(provider, multi_envelope(provider, content))
      expect_error(multi_generate(provider), "pathways collection")
    }
    observed$response <- list(status_code = "200", body = "PRIVATE_BODY")
    expect_error(multi_generate(provider), "invalid HTTP response")
    observed$response <- list(status_code = 200L, body = "")
    expect_error(multi_generate(provider), "empty response body")
    observed$response <- list(status_code = 200L, body = "PRIVATE_INVALID_JSON")
    expect_error(multi_generate(provider), "invalid JSON response envelope")
  }
})

test_that("credential echoes are rejected in raw, Unicode-decoded and JSON-key forms", {
  multi_setup()
  observed <- new.env(parent = emptyenv())
  observed$calls <- 0L
  local_mocked_bindings(.InterMetaLLMRequest = function(...) {
    observed$calls <- observed$calls + 1L
    observed$response
  }, .package = "InterMeta")
  encoded <- paste(sprintf("\\u%04x", utf8ToInt(multi_key)), collapse = "")
  contents <- c(paste0('{"pathways":[],"note":"', multi_key, '"}'),
    paste0('{"pathways":[],"note":"', encoded, '"}'),
    paste0('{"pathways":[],"', encoded, '":"secret key name"}'))
  for (provider in multi_providers) {
    for (content in contents) {
      observed$response <- multi_response(provider, multi_envelope(provider, content))
      error <- tryCatch(multi_generate(provider), error = identity)
      expect_match(conditionMessage(error), "contained credential material")
      expect_false(grepl(multi_key, conditionMessage(error), fixed = TRUE))
    }
    observed$response <- list(status_code = 200L,
      body = paste0('{"', encoded, '":"secret envelope field"}'))
    expect_error(multi_generate(provider), "contained credential material")
    messages <- multi_messages()
    messages[[2L]]$content$key <- multi_key
    before <- observed$calls
    expect_error(multi_provider(provider)$generate(messages,
      multi_fn(".InterMetaAIInterpretationSchema")()), "request was rejected")
    expect_identical(observed$calls, before)
  }
})

test_that("native provider metadata normalizes only safe scalar identifiers and token counts", {
  multi_setup()
  local_mocked_bindings(.InterMetaLLMRequest = function(provider, ...) {
    envelope <- multi_envelope(provider)
    envelope$id <- "PRIVATE_UNSAFE\nID"
    envelope$responseId <- "PRIVATE_UNSAFE\nID"
    envelope$model <- list("invalid")
    envelope$modelVersion <- list("invalid")
    envelope$usage <- list(input_tokens = -1, output_tokens = list(3),
      total_tokens = "not numeric", api_key = "PRIVATE_NOT_A_REAL_KEY")
    envelope$usageMetadata <- envelope$usage
    envelope$private_extra <- "PRIVATE_ENVELOPE"
    multi_response(provider, envelope)
  }, .package = "InterMeta")
  for (provider in multi_providers) {
    result <- multi_generate(provider)
    metadata <- attr(result, "intermeta_provider_metadata")
    expect_null(metadata$response_id)
    expect_null(metadata$response_model)
    expect_null(metadata$usage)
    expect_false(grepl("PRIVATE_", rawToChar(serialize(result, NULL, ascii = TRUE))))
  }
})
