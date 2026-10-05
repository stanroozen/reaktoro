# Thermoinvert New-Chat Handoff

## Goal

Build a general thermodynamic inverse-modeling workflow, not an enthalpy-only
tool. The eventual demonstration will fit selected phase enthalpies against
phase-equilibrium data and write a corrected thermodynamic database.

## Current State

- Implementation: `thermoinvert/`; tests: `Testing/test_thermoinvert_*.py`.
- Usage and scientific boundaries: `docs/analysis/THERMOINVERT_USAGE.md`.
- Native-unit Parameter API with backward-compatible thermodynamic aliases.
- MAP, random-walk MCMC, multi-chain diagnostics, and parallel tempering.
- Thermodynamic integration with Monte Carlo bootstrap uncertainty and ladder
  quadrature sensitivity; this is an objective normalizer, not Bayesian evidence.
- Heuristic mode clustering with separation-threshold sensitivity diagnostics.
- Correlated priors, theoretical covariance scales, observation biases, joint
  thermodynamic/P-T inversion, and grouped calibration hierarchy.
- ALR logistic-normal and Dirichlet composition likelihoods; fixed or inferred
  phase concentration; explicit aggregate below-detection censoring for zeros.
- Grouped hierarchy combines database corrections, optional concentration,
  P-T offsets, continuous observation biases, and theoretical-error scales.
- JSON/configuration/result archives and experiment/group/assemblage/domain/
  database-family holdouts, including repeated grouped calibration validation.
- Latest tranche adds normalized observation log likelihoods and validation
  log scores; inspect `experiment_log_likelihood` and validation score fields.
  Distinguish plug-in MAP holdout scores from posterior predictive mixture scores.

## Roadmap

1. Evidence diagnostics: implemented baselines; actual normalized evidence/Bayes
   factors and rigorous mode identification remain scientific limitations.
2. Simplex uncertainty: concentration inference and censored zeros implemented;
   arbitrary covariance uses logistic-normal, not Dirichlet. Structural zeros
   require an explicit model choice rather than an automatic pseudocount.
3. Unified hierarchical nuisance inference: implemented for experiment-based
   grouped calibration. Joint bracket-based inference remains a separate gap.
4. Predictive validation scores: latest implementation passes tests; audit
   normalization, censoring, archive provenance, uncertainty/coverage, and the
   distinction between MAP and posterior predictive scoring before declaring
   the entire validation roadmap complete.
5. End-to-end forward-model/database provenance: next major unfinished item.
6. Stronger convergence diagnostics/proposals: rank-normalized split R-hat,
   bulk/tail ESS, Monte Carlo errors, and improved proposal handling remain.
7. Native-unit API cleanup: preserve legacy constructors and archives while
   removing inappropriate internal/public J/mol naming.

## Verification

Current snapshot: 214 tests passed on 2026-10-05. Run PowerShell:

```powershell
$tests = Get-ChildItem Testing/test_thermoinvert_*.py | ForEach-Object { $_.FullName }
& C:\Users\stanroozen\anaconda3\envs\reaktoro\Scripts\pytest.exe -o addopts='' -q $tests
```

The direct environment pytest executable works; recent `conda run` invocations
sometimes failed without useful output. PowerShell does not expand pytest path
wildcards automatically. The JS heap error was a worker/editor error, not a
failure of the Python regression suite.

## Working Rules

- Read current files before editing: concurrent user/formatter edits occurred.
- Make narrow patches and run a focused check immediately afterward.
- Avoid broad patches in archive.py, hierarchical.py, and validation.py.
- Do not equate passing smoke tests with validated statistical calibration.
- Preserve unrelated dirty-worktree changes. This snapshot publishes only
  thermoinvert code, its tests, and its related Markdown documentation.
- No unrelated C++ changes, tutorial databases, generated plots, logs, temporary
  outputs, or reference PDF are included in this commit.

## Suggested Prompt

Read docs/analysis/THERMOINVERT_HANDOFF.md and THERMOINVERT_USAGE.md. Continue
the seven-point general inversion roadmap, verify the point-4 predictive-score
boundaries, then implement point 5 (reproducible model provenance). Preserve
unrelated changes. Our final demonstration fits phase enthalpies using
phase-equilibrium observations, but the library must stay parameter-generic.