# Optional PubMed and AI workflow

InterMeta's AI layer is downstream of the statistical analysis. It constructs
bounded PubMed searches, screens retrieved abstracts, synthesizes selected
evidence, validates citation and excerpt structure, and renders a local report
and figure. It does not establish causality, clinical utility, or independent
validation, and all manuscript statements require human review.

## Offline first step

The offline example builds a context, deterministic queries, and a report
skeleton from synthetic data. It performs no network or model request unless
the documented environment switch is deliberately enabled.

```r
source(system.file(
  "examples", "example07_ai_evidence_report.R", package = "InterMeta"
))
```

Inspect the generated queries before enabling PubMed access. A bounded search
that reaches a configured result limit is not exhaustive and should be
reported as such.

## Full recoverable workflow

The PRAD workflow example uses bundled pathway-level results rather than
patient-level expression or survival records:

```r
script <- system.file(
  "examples", "example08_ai_workflow.R", package = "InterMeta"
)
file.show(script)
```

Before sourcing it:

1. configure the provider credential in the process environment, never in the
   script, R history, or run directory;
2. inspect the context and confirm that it contains only the intended
   pathway-level information;
3. use a new `run_dir` for the first call;
4. review provider data-handling, retention, cost, and model-version terms;
5. explicitly authorize model calls only when ready.

The bundled DeepSeek example reads `DEEPSEEK_API_KEY`. Other providers use
their own documented environment variables and should not be presented as
operationally validated unless a live smoke test has been completed for the
released adapter and model identifier.

## Resume and replay

Use `InterMetaAIStatus(run_dir)` to inspect local state. Continue an interrupted
run with `InterMetaAIResume()`; keep `retry_failed = FALSE` until a failed or
delivery-uncertain attempt has been inspected because a retry can create a new
billable request.

After completion, `InterMetaAIReplay()` revalidates saved evidence and redraws
the report in a new output directory without PubMed or model calls. Do not
overwrite the original run directory.

## Publication checklist

- retain the frozen query table, retrieved identifiers, screened evidence,
  synthesis checkpoints, validated report, and figure data;
- state the retrieval date, limits, provider, model identifier, prompts or
  contract revision, and whether any query was truncated;
- have a domain expert inspect every retained claim against the cited abstract
  and, for manuscript claims, the full article;
- report quarantined or unresolved rows instead of silently repairing them;
- do not send patient-level expression, event, or follow-up data to an external
  model service.


## Repeated model calls and pathway scope

Offline replay of a frozen run reuses its stored evidence and model responses.
A new model run is different: article ordering, screening choices and synthesis
can vary even with the same model name and candidate literature. Preserve the
exact run, prompts, provider/model metadata and outputs used in a publication.
Model-assisted classifications are not additional statistical tests.

Quote matching checks whether words occur in the stored abstract; it does not
prove that the study investigates the exact pathway. In particular, `Steroid
biosynthesis` and `Steroid hormone biosynthesis` are distinct annotations.
An androgen-synthesis paper must not automatically be counted as direct evidence
for the former. Review pathway membership and distinguish direct pathway
evidence from related biological background before using a claim or figure.
This is a known interpretation limitation of 1.1.4, not corrected by replay or
by changing a model sampling parameter.
