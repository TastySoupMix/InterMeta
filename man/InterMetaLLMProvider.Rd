\name{InterMetaLLMProvider}
\alias{InterMetaLLMProvider}
\alias{InterMetaAIModels}
\title{Select a Model Service for Citation-Constrained Reports}
\usage{
InterMetaLLMProvider(
  provider = c("deepseek", "openai", "anthropic", "gemini", "qwen", "moonshot"),
  model = NULL,
  api_key_env = NULL,
  max_tokens = 16384L,
  timeout = 180
)

InterMetaAIModels()
}
\arguments{
  \item{provider}{Model service identifier. The default is \code{"deepseek"}.
  Use \code{"anthropic"} for Claude, \code{"gemini"} for Google,
  \code{"qwen"} for Alibaba, and \code{"moonshot"} for Kimi.}
  \item{model}{Optional explicit model ID supported by the selected service
  and account. \code{NULL} chooses the catalog's default for that service.
  The model is not restricted to the catalog examples.}
  \item{api_key_env}{Optional name of an environment variable containing
  the user's service credential. \code{NULL} selects the service default.
  Pass an environment variable name, never the actual API key.}
  \item{max_tokens}{Positive integer generation token budget, mapped to the
  service's output-token field. Semantics, reasoning-token accounting, and
  model limits differ between services. Truncated outputs are rejected.}
  \item{timeout}{Positive finite request timeout in seconds.}
}
\value{
\code{InterMetaLLMProvider()} returns an \code{InterMetaAIProvider} object.
\code{InterMetaAIModels()} returns a data frame of offline model examples,
including service, model, credential variable, official endpoint, protocol,
default flag, documentation URL, and date checked. The catalog is not a live
list of account-accessible models or a ranking of model quality.
}
\description{
Construct a credential-free provider using a built-in service preset. Catalog
inspection and provider construction do not read credentials or use the network.
Generation starts only when a provider is explicitly invoked by a report call.
}
\details{
The default credential variables are \code{DEEPSEEK_API_KEY},
\code{OPENAI_API_KEY}, \code{ANTHROPIC_API_KEY}, \code{GEMINI_API_KEY},
\code{DASHSCOPE_API_KEY}, and \code{MOONSHOT_API_KEY}, respectively.
Credentials are read at invocation time and are not bundled or stored in the
provider, catalog, saved report, or call metadata. Existing credentials for
another service are not reused as a fallback. A provider cannot grant model
access, account credit, or regional service availability.

DeepSeek delegates to \code{InterMetaDeepSeekProvider()}, preserving its
existing defaults. OpenAI uses Responses; Anthropic uses Messages; Gemini
uses generateContent. Qwen and Moonshot use their official compatible chat
completion endpoints. Qwen's preset targets the Beijing legacy compatibility
host; workspace-specific or other regional hosts need a custom adapter.
Moonshot uses its international host. Check regional account compatibility.
Custom endpoints and other services can be implemented explicitly with
\code{InterMetaAIProvider()} rather than redirecting these built-in clients.

The optional packages \pkg{httr2} and \pkg{jsonlite} are required at invocation.
JSON or structured output is requested using the service's protocol. Transport
schemas may omit unsupported constraints; the local report and citation
validators remain authoritative. Model availability and parameter support can
change. No model/provider fallback, automatic retry, or redirect is performed.
Non-successful, refused, incomplete, malformed, or credential-containing
responses stop with sanitized errors. Failed validation may still incur a
service charge. Never retry blindly after an ambiguous timeout.

Only aggregate report context and selected public literature are sent. Local
schema, citation membership, endpoint, and direction checks apply identically
to every provider. A different model is not independent scientific validation;
human review of the cited sources remains necessary. Tests using mocked
transports do not verify live account access or actual model output quality.
}
\seealso{\code{InterMetaAIReport()}, \code{InterMetaAIProvider()},
\code{InterMetaDeepSeekProvider()}}
\examples{
# Offline only: no credential lookup, network request, or charge.
catalog <- InterMetaAIModels()
catalog[, c("provider", "model", "api_key_env", "default")]
provider <- InterMetaLLMProvider("deepseek")
provider$name
other <- InterMetaLLMProvider("openai", model = "gpt-4.1")
other$model
}
