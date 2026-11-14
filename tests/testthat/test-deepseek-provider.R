deepseek_test_function <- function(name) getFromNamespace(name, "InterMeta")

deepseek_test_key <- "dummy-intermeta-deepseek-key-NOT-REAL"
deepseek_test_env <- "INTERMETA_DEEPSEEK_TEST_KEY"

deepseek_test_setup <- function(.env = parent.frame()) {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("httr2")
  withr::local_envvar(
    stats::setNames(deepseek_test_key, deepseek_test_env), .local_envir = .env
  )
}

deepseek_test_provider <- function(...) {
  deepseek_test_function("InterMetaDeepSeekProvider")(
    api_key_env = deepseek_test_env, ...
  )
}

deepseek_test_messages <- function() {
  list(
    list(role = "system", content = "Use only the supplied aggregate evidence."),
    list(role = "user", content = list(
      language = "en", context = list(pathways = data.frame(
        Pathway = "P1", Direction = "Favorable", stringsAsFactors = FALSE
      )), evidence = list(records = data.frame(
        pathway = "P1", PMID = "12345", Title = "A supplied publication",
        stringsAsFactors = FALSE
      ))
    ))
  )
}

deepseek_test_schema <- function() {
  deepseek_test_function(".InterMetaAIInterpretationSchema")()
}

deepseek_test_content <- function(pmid = "12345", consistency = "supported") {
  as.character(jsonlite::toJSON(list(pathways = list(list(
    pathway = "P1", literature_consistency = consistency,
    summary = "The supplied evidence is consistent with this pathway.",
    claims = list(list(claim = "The publication discusses this pathway.",
                       pmids = list(pmid)))
  ))), auto_unbox = TRUE))
}

deepseek_test_response <- function(content = deepseek_test_content(),
                                   finish_reason = "stop", status_code = 200L) {
  envelope <- list(
    id = "offline-response-1", model = "deepseek-v4-pro", created = 1788393600,
    choices = list(list(
      finish_reason = finish_reason,
      message = list(role = "assistant", content = content,
                     reasoning_content = "PRIVATE_REASONING_NOT_TO_RETAIN")
    )),
    usage = list(prompt_tokens = 12L, completion_tokens = 34L, total_tokens = 46L,
                 prompt_cache_hit_tokens = 2L, prompt_cache_miss_tokens = 10L,
                 arbitrary = "PRIVATE_USAGE_NOT_TO_RETAIN"),
    arbitrary = "PRIVATE_ENVELOPE_NOT_TO_RETAIN"
  )
  list(status_code = status_code, body = as.character(jsonlite::toJSON(
    envelope, auto_unbox = TRUE, null = "null"
  )))
}

deepseek_test_mock <- function(response = deepseek_test_response(),
                              .env = parent.frame()) {
  testthat::local_mocked_bindings(
    .InterMetaDeepSeekRequest = function(endpoint, body, api_key, timeout) response,
    .package = "InterMeta", .env = .env
  )
}

deepseek_test_generate <- function(provider = deepseek_test_provider()) {
  provider$generate(deepseek_test_messages(), deepseek_test_schema())
}

deepseek_test_context <- function() {
  structure(list(
    analysis = list(disease = "prostate cancer", endpoint = "overall survival"),
    pathways = data.frame(Pathway = "P1", Direction = "Favorable",
                           stringsAsFactors = FALSE),
    interactions = data.frame(), driver_gene_pairs = data.frame(),
    warnings = character(), provenance = list(input_hash = "test-only")
  ), class = c("InterMetaAIContext", "list"))
}

deepseek_test_evidence <- function() {
  evidence <- deepseek_test_function(".InterMetaAIEmptyEvidence")()
  evidence$records <- data.frame(
    query_id = "q1", pathway = "P1", direction = "Favorable", layer = "biology",
    PMID = "12345", Title = "A supplied publication", Abstract = "Evidence.",
    stringsAsFactors = FALSE
  )
  evidence
}

test_that("DeepSeek constructor is offline and validates all public parameters", {
  constructor <- deepseek_test_function("InterMetaDeepSeekProvider")
  testthat::local_mocked_bindings(
    .InterMetaDeepSeekRequire = function(package) stop("must not load dependencies"),
    .InterMetaDeepSeekRequest = function(...) stop("must not request"),
    .package = "InterMeta"
  )
  provider <- constructor()
  expect_s3_class(provider, "InterMetaAIProvider")
  expect_identical(provider$name, "DeepSeek")
  expect_identical(provider$model, "deepseek-v4-pro")
  expect_identical(names(provider), c("generate", "name", "model"))
  expect_identical(names(formals(constructor)), c(
    "model", "api_key_env", "max_tokens", "timeout", "thinking", "reasoning_effort"
  ))
  for (value in list("", NA_character_, c("a", "b"), 1, "https://other.example",
                     "deepseek\nheader", " deepseek-v4-pro")) {
    expect_error(constructor(model = value), "model must")
  }
  for (value in list("", NA_character_, c("A", "B"), "A=B", "1KEY", "A\nB")) {
    expect_error(constructor(api_key_env = value), "api_key_env must")
  }
  for (value in list(0, -1, 1.5, NA_real_, Inf, "100", TRUE, c(1, 2), 2^31)) {
    expect_error(constructor(max_tokens = value), "max_tokens must")
  }
  for (value in list(0, -1, NA_real_, Inf, "10", TRUE, c(1, 2))) {
    expect_error(constructor(timeout = value), "timeout must")
  }
  for (value in list("en", "ENABLED", c("enabled", "disabled"), NA, "")) {
    expect_error(constructor(thinking = value), "thinking must")
  }
  for (value in list("h", "HIGH", c("high", "low"), NA, "")) {
    expect_error(constructor(reasoning_effort = value), "reasoning_effort must")
  }
})

test_that("credentials are read lazily and never retained by the provider", {
  deepseek_test_setup()
  withr::local_envvar(stats::setNames(NA_character_, deepseek_test_env))
  provider <- deepseek_test_provider()
  expect_error(deepseek_test_generate(provider), "API key is missing")
  observed <- new.env(parent = emptyenv())
  observed$keys <- character()
  testthat::local_mocked_bindings(
    .InterMetaDeepSeekRequest = function(endpoint, body, api_key, timeout) {
      observed$keys <- c(observed$keys, api_key)
      deepseek_test_response()
    }, .package = "InterMeta"
  )
  keys <- c(deepseek_test_key, paste0(deepseek_test_key, "-ROTATED"))
  for (key in keys) {
    do.call(Sys.setenv, stats::setNames(list(key), deepseek_test_env))
    deepseek_test_generate(provider)
    serialized <- rawToChar(serialize(provider, NULL, ascii = TRUE))
    expect_false(grepl(key, serialized, fixed = TRUE))
    expect_false("api_key" %in% ls(environment(provider$generate), all.names = TRUE))
  }
  expect_identical(observed$keys, keys)
  do.call(Sys.setenv, stats::setNames(list("  \n "), deepseek_test_env))
  expect_error(deepseek_test_generate(provider), "API key is missing")
  do.call(Sys.setenv, stats::setNames(list("bad token"), deepseek_test_env))
  expect_error(deepseek_test_generate(provider), "single token")
})

test_that("DeepSeek serializes messages, row tables, schema arrays, and JSON mode", {
  deepseek_test_setup()
  observed <- new.env(parent = emptyenv())
  response <- deepseek_test_response()
  testthat::local_mocked_bindings(
    .InterMetaDeepSeekRequest = function(endpoint, body, api_key, timeout) {
      observed$endpoint <- endpoint
      observed$body <- body
      observed$timeout <- timeout
      response
    }, .package = "InterMeta"
  )
  result <- deepseek_test_generate()
  body <- observed$body
  expect_identical(observed$endpoint, "https://api.deepseek.com/chat/completions")
  expect_equal(observed$timeout, 180)
  expect_identical(body$model, "deepseek-v4-pro")
  expect_identical(body$response_format, list(type = "json_object"))
  expect_identical(body$thinking, list(type = "enabled"))
  expect_identical(body$reasoning_effort, "high")
  expect_identical(body$max_tokens, 16384L)
  expect_false(body$stream)
  expect_null(body$temperature)
  expect_true(all(vapply(body$messages, function(x) {
    is.character(x$content) && length(x$content) == 1L
  }, logical(1))))
  expect_match(body$messages[[1L]]$content, "valid json object", fixed = TRUE)
  expect_match(body$messages[[1L]]$content, '"required":["pathways"]', fixed = TRUE)
  expect_match(body$messages[[1L]]$content, '"claims":[]', fixed = TRUE)
  user <- jsonlite::parse_json(body$messages[[2L]]$content)
  expect_identical(user$context$pathways[[1L]]$Pathway, "P1")
  expect_identical(user$evidence$records[[1L]]$PMID, "12345")
  expect_false(grepl(deepseek_test_key,
                     deepseek_test_function(".InterMetaDeepSeekJSON")(body), fixed = TRUE))
  metadata <- attr(result, "intermeta_provider_metadata")
  expect_identical(metadata$response_model, "deepseek-v4-pro")
  expect_identical(metadata$response_id, "offline-response-1")
  expect_identical(metadata$finish_reason, "stop")
  expect_equal(metadata$usage$total_tokens, 46)
  expect_setequal(names(metadata$usage), c(
    "prompt_tokens", "completion_tokens", "total_tokens",
    "prompt_cache_hit_tokens", "prompt_cache_miss_tokens"
  ))
  expect_identical(metadata$request_sha256, digest::digest(
    deepseek_test_function(".InterMetaDeepSeekJSON")(body),
    algo = "sha256", serialize = FALSE
  ))
  expect_identical(metadata$response_sha256,
                   digest::digest(response$body, algo = "sha256", serialize = FALSE))
  serialized <- rawToChar(serialize(result, NULL, ascii = TRUE))
  expect_false(grepl("PRIVATE_", serialized, fixed = TRUE))
  expect_false(grepl(deepseek_test_key, serialized, fixed = TRUE))
})

test_that("thinking modes and explicit model parameters are honored without fallback", {
  deepseek_test_setup()
  bodies <- list()
  testthat::local_mocked_bindings(
    .InterMetaDeepSeekRequest = function(endpoint, body, api_key, timeout) {
      bodies[[length(bodies) + 1L]] <<- body
      expect_equal(timeout, 15.5)
      deepseek_test_response()
    }, .package = "InterMeta"
  )
  for (effort in c("high", "low", "max")) {
    deepseek_test_generate(deepseek_test_provider(
      model = "deepseek-v4-flash", timeout = 15.5, max_tokens = 2048,
      reasoning_effort = effort
    ))
    body <- bodies[[length(bodies)]]
    expect_identical(body$model, "deepseek-v4-flash")
    expect_identical(body$max_tokens, 2048L)
    expect_identical(body$reasoning_effort, effort)
    expect_null(body$temperature)
  }
  result <- deepseek_test_generate(deepseek_test_provider(
    thinking = "disabled", reasoning_effort = "max", timeout = 15.5
  ))
  expect_identical(bodies[[4L]]$thinking$type, "disabled")
  expect_null(bodies[[4L]]$reasoning_effort)
  expect_identical(bodies[[4L]]$temperature, 0)
  expect_null(attr(result, "intermeta_provider_metadata")$reasoning_effort)
  expect_length(bodies, 4L)
})

test_that("HTTP transport is one fixed POST with headers, timeout, and no redirect", {
  deepseek_test_setup()
  observed <- new.env(parent = emptyenv())
  observed$calls <- 0L
  testthat::local_mocked_bindings(
    req_perform = function(req, ...) {
      observed$request <- req
      observed$calls <- observed$calls + 1L
      list(status = 200L, body = '{"choices":[]}')
    },
    resp_status = function(resp) resp$status,
    resp_body_string = function(resp) resp$body,
    .package = "httr2"
  )
  request <- deepseek_test_function(".InterMetaDeepSeekRequest")
  result <- request("https://api.deepseek.com/chat/completions",
                     list(model = "deepseek-v4-pro"), deepseek_test_key, 21)
  expect_equal(result$status_code, 200)
  req <- observed$request
  expect_identical(req$url, "https://api.deepseek.com/chat/completions")
  expect_identical(req$method, "POST")
  expect_identical(req$headers$Authorization, paste("Bearer", deepseek_test_key))
  expect_false(req$options$followlocation)
  expect_equal(req$options$maxredirs, 0)
  expect_false(req$options$verbose)
  expect_equal(req$options$timeout_ms, 21000)
  expect_equal(req$policies$retry_max_tries, 1)
  expect_false(req$policies$retry_is_transient(list()))
  expect_identical(observed$calls, 1L)
  expect_error(request("https://other.example", list(), deepseek_test_key, 1),
               "fixed official endpoint")
  expect_identical(observed$calls, 1L)
})

test_that("HTTP and network failures are safely redacted and never retried", {
  deepseek_test_setup()
  observed <- new.env(parent = emptyenv())
  observed$status <- 401L
  observed$calls <- 0L
  testthat::local_mocked_bindings(
    .InterMetaDeepSeekRequest = function(endpoint, body, api_key, timeout) {
      observed$calls <- observed$calls + 1L
      list(status_code = observed$status, body = paste("PRIVATE_ERROR", api_key))
    }, .package = "InterMeta"
  )
  statuses <- c(301L, 400L, 401L, 402L, 403L, 404L, 408L, 422L, 429L, 500L, 503L)
  for (status in statuses) {
    observed$status <- status
    error <- tryCatch(deepseek_test_generate(), error = identity)
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error), paste("DeepSeek HTTP", status), fixed = TRUE)
    expect_false(grepl(deepseek_test_key, conditionMessage(error), fixed = TRUE))
    expect_false(grepl("PRIVATE_ERROR", conditionMessage(error), fixed = TRUE))
    expect_null(conditionCall(error))
  }
  expect_identical(observed$calls, length(statuses))
  testthat::local_mocked_bindings(
    .InterMetaDeepSeekRequest = function(...) {
      warning(paste("PRIVATE_WARNING", deepseek_test_key))
      stop(paste("Timeout PRIVATE_ERROR", deepseek_test_key))
    }, .package = "InterMeta"
  )
  error <- tryCatch(deepseek_test_generate(), error = identity)
  expect_match(conditionMessage(error), "network, timeout, or transport error", fixed = TRUE)
  expect_false(grepl(deepseek_test_key, conditionMessage(error), fixed = TRUE))
  expect_false(grepl("PRIVATE_", conditionMessage(error), fixed = TRUE))
  expect_null(conditionCall(error))
})

test_that("truncated, empty, malformed, and non-JSON completions are rejected", {
  deepseek_test_setup()
  observed <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    .InterMetaDeepSeekRequest = function(...) observed$response,
    .package = "InterMeta"
  )
  for (reason in list("length", "content_filter", "tool_calls",
                      "insufficient_system_resource")) {
    observed$response <- deepseek_test_response(finish_reason = reason)
    error <- tryCatch(deepseek_test_generate(), error = identity)
    expect_match(conditionMessage(error), "not finished normally", fixed = TRUE)
    expect_match(conditionMessage(error), paste0("finish_reason=", reason), fixed = TRUE)
    expect_null(conditionCall(error))
  }
  for (reason in list("PRIVATE_REASON", NULL)) {
    observed$response <- deepseek_test_response(finish_reason = reason)
    error <- tryCatch(deepseek_test_generate(), error = identity)
    expect_match(conditionMessage(error), "finish_reason=unknown", fixed = TRUE)
    expect_false(grepl("PRIVATE_REASON", conditionMessage(error), fixed = TRUE))
    expect_null(conditionCall(error))
  }
  for (content in list("", " \n ", NULL, list(text = "wrong type"))) {
    observed$response <- deepseek_test_response(content = content)
    expect_error(deepseek_test_generate(), "empty completion content")
  }
  for (content in c('{"pathways":', '```json\n{"pathways":[]}\n```',
                    "PRIVATE_INVALID_JSON", "https://invalid.example/no-network")) {
    observed$response <- deepseek_test_response(content = content)
    error <- tryCatch(deepseek_test_generate(), error = identity)
    expect_match(conditionMessage(error), "not valid JSON", fixed = TRUE)
    expect_false(grepl(content, conditionMessage(error), fixed = TRUE))
  }
  for (content in c('{}', '{"pathways":null}', '{"pathways":"bad"}', '[]')) {
    observed$response <- deepseek_test_response(content = content)
    expect_error(deepseek_test_generate(), "pathways collection")
  }
  for (body in list("", "  ", NULL, 1)) {
    observed$response <- list(status_code = 200L, body = body)
    expect_error(deepseek_test_generate(), "empty response body")
  }
  observed$response <- list(status_code = 200L, body = "PRIVATE_INVALID_ENVELOPE")
  expect_error(deepseek_test_generate(), "invalid JSON response envelope")
  for (body in c('{}', '{"choices":[]}', '{"choices":[null]}', 'null')) {
    observed$response <- list(status_code = 200L, body = body)
    expect_error(deepseek_test_generate(), "missing a completion choice")
  }
  observed$response <- list(status_code = "200", body = "PRIVATE_STATUS")
  expect_error(deepseek_test_generate(), "invalid HTTP response")
})

test_that("credential echoes are rejected even after two JSON decoding layers", {
  deepseek_test_setup()
  observed <- new.env(parent = emptyenv())
  observed$calls <- 0L
  testthat::local_mocked_bindings(
    .InterMetaDeepSeekRequest = function(...) {
      observed$calls <- observed$calls + 1L
      observed$response
    }, .package = "InterMeta"
  )
  unicode_key <- paste(sprintf("\\u%04x", utf8ToInt(deepseek_test_key)), collapse = "")
  literal_content <- paste0('{"pathways":[],"note":"', deepseek_test_key, '"}')
  encoded_content <- paste0('{"pathways":[],"note":"', unicode_key, '"}')
  bodies <- list(
    deepseek_test_response(content = literal_content),
    deepseek_test_response(content = encoded_content),
    list(status_code = 200L, body = paste0('{"id":"', unicode_key, '"}'))
  )
  for (response in bodies) {
    observed$response <- response
    error <- tryCatch(deepseek_test_generate(), error = identity)
    expect_match(conditionMessage(error), "contained credential material", fixed = TRUE)
    expect_false(grepl(deepseek_test_key, conditionMessage(error), fixed = TRUE))
  }
  messages <- deepseek_test_messages()
  messages[[2L]]$content$note <- deepseek_test_key
  expect_error(deepseek_test_provider()$generate(messages, deepseek_test_schema()),
               "request was rejected")
  expect_identical(observed$calls, length(bodies))
})

test_that("only safe scalar provider metadata survives", {
  deepseek_test_setup()
  response <- deepseek_test_response()
  envelope <- jsonlite::parse_json(response$body)
  envelope$id <- "unsafe\nmetadata"
  envelope$model <- list("not scalar")
  envelope$created <- "not numeric"
  envelope$usage$prompt_tokens <- "not numeric"
  envelope$usage$completion_tokens <- -1
  envelope$usage$total_tokens <- list(1, 2)
  response$body <- as.character(jsonlite::toJSON(envelope, auto_unbox = TRUE))
  deepseek_test_mock(response)
  metadata <- attr(deepseek_test_generate(), "intermeta_provider_metadata")
  expect_null(metadata$response_id)
  expect_null(metadata$response_model)
  expect_null(metadata$created)
  expect_null(metadata$usage$prompt_tokens)
  expect_null(metadata$usage$completion_tokens)
  expect_null(metadata$usage$total_tokens)
  expect_setequal(names(metadata$usage), c("prompt_cache_hit_tokens", "prompt_cache_miss_tokens"))
  expect_true(all(names(metadata) %in% c(
    "endpoint", "requested_model", "response_model", "response_id", "created",
    "finish_reason", "thinking", "reasoning_effort", "max_tokens",
    "request_sha256", "response_sha256", "usage"
  )))
})

test_that("DeepSeek output goes through real InterMeta citation and enum validation", {
  deepseek_test_setup()
  observed <- new.env(parent = emptyenv())
  observed$response <- deepseek_test_response()
  testthat::local_mocked_bindings(
    .InterMetaDeepSeekRequest = function(...) observed$response,
    .package = "InterMeta"
  )
  interpret <- deepseek_test_function("InterMetaAIInterpret")
  run <- function() interpret(
    deepseek_test_context(), deepseek_test_evidence(), deepseek_test_provider(),
    purpose = "prognostic_validation"
  )
  result <- run()
  expect_s3_class(result, "InterMetaAIInterpretation")
  expect_identical(result$claims$PMID, "12345")
  expect_identical(result$provider$name, "DeepSeek")
  expect_identical(result$provider_call$response_id, "offline-response-1")
  observed$response <- deepseek_test_response(content = deepseek_test_content("99999"))
  expect_error(run(), "absent from pathway-specific evidence")
  observed$response <- deepseek_test_response(content = deepseek_test_content(consistency = "proven"))
  expect_error(run(), "Invalid literature_consistency")
})
