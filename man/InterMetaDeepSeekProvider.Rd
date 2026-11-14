\name{InterMetaDeepSeekProvider}
\alias{InterMetaDeepSeekProvider}
\title{Create a DeepSeek JSON Report Provider}
\usage{
InterMetaDeepSeekProvider(
  model = "deepseek-v4-pro",
  api_key_env = "DEEPSEEK_API_KEY",
  max_tokens = 16384L,
  timeout = 180,
  thinking = c("enabled", "disabled"),
  reasoning_effort = c("high", "low", "max")
)
}
\arguments{
  \item{model}{Non-empty DeepSeek model identifier. The default is
  \code{deepseek-v4-pro}; \code{deepseek-v4-flash} can be selected explicitly.}
  \item{api_key_env}{Name of the environment variable containing the API key.
  The credential is read only when the provider is invoked. Do not pass a
  credential as this argument.}
  \item{max_tokens}{Positive integer output token limit. An insufficient
  budget may result in a truncated response, which is rejected.}
  \item{timeout}{Positive finite request timeout in seconds.}
  \item{thinking}{One of \code{"enabled"} or \code{"disabled"}.}
  \item{reasoning_effort}{One of \code{"high"}, \code{"low"}, or
  \code{"max"}. Only sent when thinking is enabled.}
}
\value{An \code{InterMetaAIProvider} object for use with
\code{InterMetaAIInterpret()} or \code{InterMetaAIReport()}.}
\description{
Creates a provider for the official DeepSeek chat-completions endpoint.
Construction does not read credentials or perform a network request.
}
\details{
The provider uses \url{https://api.deepseek.com/chat/completions}, not the
Anthropic compatibility endpoint. It requires the optional packages
\pkg{httr2} and \pkg{jsonlite} when invoked. Aggregate messages and the output
contract are sent using JSON mode; local validation is still required because
JSON mode alone does not enforce the schema or scientific correctness.

There are no automatic retries, model fallbacks, or redirects. A request may
incur API charges, including when later validation fails. Empty, malformed,
  incomplete, non-successful, or credential-containing responses are rejected
  without exposing response bodies in error messages. Only final content is
  parsed; reasoning content is not included in the report. The returned model
  may be a service-side version alias and is recorded separately.

  When a completion does not end with code{stop}, the error reports only the
  provider's allowlisted code{finish_reason}. Unknown values are replaced by
  code{unknown}; raw or partial completion text is never included.

The provider attaches allowlisted call metadata to its parsed result. Report
functions retain safe model identifiers, token usage, settings, and SHA-256
request/response hashes. They do not retain credentials, raw response bodies,
or reasoning text. Keep credentials out of source files, notebooks, command
history, and saved objects. Environment variables must be configured securely
in the process running R.
}
\seealso{\code{InterMetaAIReport()}, \code{InterMetaAIProvider()},
\url{https://api-docs.deepseek.com/guides/json_mode/}}
\examples{
# No network call or credential lookup occurs here.
provider <- InterMetaDeepSeekProvider()
provider$name
provider$model
}
