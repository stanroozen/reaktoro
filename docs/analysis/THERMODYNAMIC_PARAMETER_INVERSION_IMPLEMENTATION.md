# Thermodynamic Parameter Inversion Implementation

## 1. Purpose

Build a Python framework for estimating corrections to selected thermodynamic parameters from sparse experimental phase-equilibrium data, while preserving an existing internally consistent thermodynamic database as the prior model.

The initial parameterization is:

$$
H_i^{new} = H_i^{prior} + \delta H_i
$$

where the enthalpy correction $\delta H_i$ is inferred from data. The framework must avoid silently replacing the source database or interpreting numerical optimization as proof of physical identifiability.

Perple_X and Reaktoro are both supported as possible forward-model backends. The current repository implementation uses Reaktoro; the interface must remain backend-independent.

## 2. Scientific Scope

### Initial observations

Experiments may provide:

- pressure and temperature, with uncertainties
- bulk composition, with uncertainties
- stable and absent phases
- mineral or endmember compositions
- phase proportions
- experimental equilibrium brackets
- oxygen fugacity or other intensive variables
- study and laboratory metadata

Natural equilibrium samples are initially independent validation data. They should not receive calibration weight unless explicitly enabled.

### Initial unknowns

Start with a small set of selected standard/reference enthalpy corrections:

```text
[delta_H_phase_1, delta_H_phase_2, ...]
```

Later parameter types may include entropy, heat capacity, volume/EOS, and solution interaction parameters, but they must not be enabled by default.

## 3. Design Principles

1. Preserve the source database.
2. Use temporary model overlays for every parameter proposal.
3. Keep the forward thermodynamic engine behind an abstract interface.
4. Run identifiability diagnostics before expensive optimization.
5. Use priors and bounded optimization for sparse data.
6. Treat phase presence, compositions, modes, and brackets as separate likelihood terms.
7. Distinguish experimental uncertainty, thermodynamic prior uncertainty, inverse non-uniqueness, and model structural error.
8. Validate on withheld experiments, studies, and out-of-domain conditions.
9. Report parameter combinations and correlations when individual parameters are not identifiable.
10. Do not begin with MCMC; first establish a deterministic, tested forward loop.

## 4. Target Architecture

```text
thermoinvert/
├── data/
│   ├── experiments.py       # Experiment and observation models
│   ├── compositions.py      # Normalization and covariance utilities
│   └── validation.py        # Calibration/validation splits
├── thermo/
│   ├── database.py          # Safe temporary parameter overlays
│   ├── parameters.py        # Parameter and prior definitions
│   ├── forward_model.py     # Backend-independent contract
│   └── perplex.py           # Perple_X adapter, when enabled
├── likelihood/
│   ├── composition.py       # Covariance-aware composition likelihood
│   ├── assemblage.py        # Stable/absent phase likelihoods
│   ├── bracket.py           # Inequality and bracket terms
│   ├── phase_fraction.py    # Modal/fraction likelihoods
│   └── combined.py          # Weighted experiment likelihood
├── inversion/
│   ├── priors.py            # Gaussian and correlated priors
│   ├── identifiability.py   # Structural and numerical rank analysis
│   ├── sensitivity.py       # Finite-difference forward sensitivities
│   ├── map.py               # Prior-regularized MAP objective
│   ├── multistart.py        # Derivative-free multi-start optimization
│   └── posterior.py         # Later posterior sampling
├── diagnostics/
│   ├── correlations.py      # Covariance, correlation, null spaces
│   ├── posterior_predictive.py
│   ├── cross_validation.py
│   └── information.py       # Information gain and prior/posterior comparison
├── io/
│   ├── experimental_reader.py
│   ├── perplex_parser.py
│   └── results.py
└── tests/
```

## 5. Core Data Contracts

### Experiment

```python
@dataclass
class Experiment:
    id: str
    study: str
    pressure_bar: float
    pressure_sigma_bar: float | None
    temperature_k: float
    temperature_sigma_k: float | None
    bulk_composition: dict[str, float]
    bulk_sigma: dict[str, float] | None
    observed_phases: list[str]
    absent_phases: list[str]
    phase_compositions: dict[str, dict[str, float]] | None
    phase_composition_sigma: dict[str, dict[str, float]] | None
    phase_fractions: dict[str, float] | None
    phase_fraction_sigma: dict[str, float] | None
    metadata: dict[str, object]
```

### Thermodynamic parameter

```python
@dataclass
class ThermodynamicParameter:
    name: str
    phase: str
    parameter_type: str  # initially "dHf"
    prior_mean: float  # J/mol correction, normally zero
    prior_sigma: float  # J/mol
    lower_bound: float | None
    upper_bound: float | None
```

### Equilibrium result

```python
@dataclass
class EquilibriumResult:
    valid: bool
    stable_phases: tuple[str, ...]
    phase_compositions: dict[str, dict[str, float]]
    phase_fractions: dict[str, float]
    chemical_potentials: dict[str, float] | None
    phase_affinities: dict[str, float] | None
    diagnostics: dict[str, object]
```

### Forward-model contract

```python
class ThermodynamicForwardModel(Protocol):
    def predict(
        self,
        experiment: Experiment,
        parameters: Mapping[str, float],
    ) -> EquilibriumResult: ...
```

## 6. Safe Thermodynamic Parameter Overlays

For each proposal $\theta$:

1. Load the original database read-only.
2. Copy it into a temporary model representation.
3. Apply only explicitly selected corrections.
4. Run the forward calculation against the temporary model.
5. Delete the temporary model after the evaluation.
6. Never edit the user’s source database in place.

The implementation must record exactly which database fields are modified. A correction must not be applied simultaneously to multiple thermodynamic fields unless the database convention proves that this is required. For example, changing `Hf`, `H0`, `G0`, and Perple_X metadata independently can double-count a correction.

For a standard enthalpy correction, document whether the backend interprets the correction as:

$$
G_i(P,T;\delta H_i) = G_i^{prior}(P,T) + \delta H_i
$$

or as a change to a reference enthalpy term propagated through the equation of state.

## 7. Likelihood Model

The posterior target is:

$$
p(\theta, \eta \mid D) \propto p(D \mid \theta, \eta)\,p(\theta)\,p(\eta)
$$

where $\eta$ contains optional nuisance variables such as true pressure, true temperature, bulk uncertainty, or study offsets.

The negative log posterior is:

$$
\Phi(\theta,\eta) = -\log p(D\mid\theta,\eta) - \log p(\theta) - \log p(\eta)
$$

### Composition term

For observed composition $x_{obs}$ and prediction $x_{calc}$:

$$
\chi_X^2 = (x_{obs}-x_{calc})^T\Sigma_X^{-1}(x_{obs}-x_{calc})
$$

Use covariance when available. Do not assume composition components are independent if closure creates covariance.

### Phase-fraction term

Use uncertainty-weighted residuals initially, with the closure limitation recorded explicitly. A later implementation may use a Dirichlet, logistic-normal, or transformed-composition likelihood.

### Stable and absent phases

Provide modular terms for:

- observed phase present but not predicted
- predicted phase absent from observation
- optional smooth stability/affinity information

Avoid making an undocumented hard integer penalty the only phase likelihood.

### Experimental brackets

Represent bracket observations as inequalities or censored likelihoods. For example:

```text
At T_low: phase A stable
At T_high: phase B stable
```

can constrain a transition to:

$$
T_{low} < T_{boundary}(\theta) < T_{high}
$$

### P-T-X uncertainty

Treat reported conditions as observations when uncertainty is available:

$$
T_{true}\sim N(T_{reported},\sigma_T^2),\qquad
P_{true}\sim N(P_{reported},\sigma_P^2)
$$

Do not allow enthalpy corrections to absorb calibration uncertainty automatically.

## 8. Thermodynamic Priors

For selected enthalpy corrections:

$$
\delta H_i \sim N(0,\sigma_{H_i}^2)
$$

With correlated prior information:

$$
-\log p(\theta) =
\frac{1}{2}(\theta-\theta_0)^T\Sigma_H^{-1}(\theta-\theta_0)
$$

The prior must record:

- original database value
- correction prior mean
- prior standard deviation
- covariance, if known
- source and provenance

## 9. Identifiability Before Optimization

### Structural analysis

Construct reaction or stoichiometric sensitivity rows where possible. For example:

```text
          H_A  H_B  H_C  H_D
reaction1  -1   -1   +1   +1
reaction2   0   +1   -1    0
```

Use SVD:

$$
J = U\,S\,V^T
$$

Report rank, singular values, weak directions, and null-space vectors.

### Numerical analysis

Evaluate the actual forward model at perturbed parameters:

$$
J_{ij}\approx
\frac{y_i(\theta+\epsilon e_j)-y_i(\theta-\epsilon e_j)}{2\epsilon}
$$

Use one-sided differences near bounds and mark sensitivities near phase boundaries as potentially nonsmooth.

Do not report individual parameter precision when only a combination is constrained. Report combinations such as:

$$
\delta H_A - 0.82\,\delta H_B
$$

### Prior-vs-posterior shrinkage (uncertainty width is not identifiability)

A narrow posterior does not by itself mean a parameter is data-constrained --
it can simply reflect a narrow prior. `identifiability.identifiability_summary`
compares a `laplace_covariance`/`mcmc_sample` result run with
`uncertainty_source="all"` against the same diagnostic run with
`uncertainty_source="prior"`, and reports, per parameter:

$$
\text{shrinkage} = 1 - \frac{\sigma^2_{\text{posterior}}}{\sigma^2_{\text{prior}}}
$$

classified as `data_constrained` (shrinkage above a threshold, default 0.5),
`prior_dominated` (shrinkage below a threshold, default 0.1), or `ambiguous`
between them. It optionally cross-references `identify_sensitivity`'s
`weak_directions`: a parameter with a large component in any weak direction
is flagged `structurally_weak=True` regardless of its shrinkage, since that
reflects a forward-model/experiment-design limitation no amount of data can
fix. Validated against analytic quadratic cases (data-dominated and
prior-dominated) and a deliberately collinear two-parameter case.

## 10. MAP and Multi-Start Optimization

The first optimizer should minimize the prior-regularized objective:

```text
for theta_start in prior_or_bound_starts:
    optimize Phi(theta) with a bounded derivative-free method
retain and compare converged candidates
```

Recommended initial algorithms:

- bounded Nelder-Mead
- Powell
- differential evolution for diagnostic global searches
- direct two-dimensional grids for verification

Record every start, convergence status, objective, and final parameter vector. If substantially different solutions have similar objectives, report non-uniqueness.

## 11. Posterior and Uncertainty Stages

Do not rely only on an inverse-Hessian covariance estimate.

Later stages should support:

- posterior grids for one or two parameters (done: `grid_posterior`)
- `emcee`, `dynesty`, or another appropriate sampler (done, hand-rolled: `mcmc_sample`/`mcmc_sample_multi_chain`)
- asymmetric credible intervals (done: MCMC percentile-based `credible_interval_95`)
- covariance and correlation matrices (done: `laplace_covariance`)
- prior versus posterior width (done: `identifiability.identifiability_summary` shrinkage)
- information gain (done: `identifiability.information_gain_summary`, KL(posterior || prior) in nats, using each parameter's exact Gaussian prior directly)
- posterior predictive distributions (done: `posterior_predictive_check`)

Separate:

1. experimental uncertainty
2. thermodynamic prior uncertainty
3. inverse non-uniqueness
4. structural model uncertainty

## 12. Natural-Sample Validation

Natural samples should initially be excluded from calibration and used for validation:

```text
fit experimental calibration set
predict natural samples with posterior/MAP parameters
report phase and composition agreement
```

Later validation should include:

- withheld experiments
- leave-one-bracket-out tests
- leave-one-study-out tests
- out-of-calibration P-T sections
- original versus corrected database comparisons

A lower calibration residual alone is not evidence that a correction is physically valid.

## 13. Reproducibility and Outputs

Every run must record:

- source database path and hash
- backend and version
- experimental dataset and version
- parameter definitions and bounds
- prior distributions and covariance
- likelihood definitions and weights
- optimizer and settings
- random seed
- code version/commit
- calibration/validation split
- generated temporary or corrected database provenance

Target outputs:

```text
results/
├── map_result.json
├── parameter_summary.csv
├── covariance.csv
├── correlation.csv
├── sensitivity_matrix.csv
├── singular_values.csv
├── prior_scaled_information.csv
├── identifiability.json
├── predictions_calibration.csv
├── predictions_validation.csv
├── posterior_predictive.csv
└── corrected_database_<label>.json
```

The scientifically meaningful result is the posterior parameter report, not merely a modified database file.

## 14. Development Phases

### Phase 1: deterministic end-to-end loop

- experiment data model
- thermodynamic parameter model
- safe Reaktoro database overlay
- equilibrium prediction
- stable phase extraction
- simple phase/fraction likelihood
- Gaussian parameter prior
- one-parameter bounded optimization
- synthetic tests

### Phase 2: identifiability and multi-parameter MAP

- multiple enthalpy corrections
- numerical sensitivities
- structural sensitivity matrix
- SVD/rank diagnostics
- multi-start optimization
- parameter correlation reporting

### Phase 3: richer likelihoods and uncertainty

- composition covariance
- phase-fraction closure-aware likelihood
- phase presence/absence likelihoods
- bracket likelihoods
- uncertain P-T and bulk nuisance variables
- posterior sampling

### Phase 4: validation

- posterior predictive checks
- withheld experiments
- leave-one-study-out validation
- natural-sample validation
- out-of-domain phase-diagram checks

### Phase 5: advanced models

- hierarchical shrinkage priors
- study-specific calibration offsets
- additional thermodynamic parameter types
- checkpointing and caching
- process-isolated parallel forward evaluations

## 15. Current Repository Status

The repository currently contains an early Phase 1 scaffold:

- `thermoinvert/data.py`: data structures
- `thermoinvert/thermo.py`: temporary JSON overlay
- `thermoinvert/forward_model.py`: Reaktoro prediction wrapper
- `thermoinvert/likelihood.py`: phase/fraction misfit with optional explicitly keyed positive-definite analytical observation covariance, evaluated as a Mahalanobis term
- `logit_lsq` provides an opt-in transformed likelihood for independent compositions/fractions in $(0,1)$; covariance likelihoods remain on `weighted_lsq`
- `Experiment.observation_bias` supports known finite measurement/laboratory offsets for continuous observations and is serialized with the dataset
- `profile_observation_bias_by_group` estimates one fixed additive bias per study/laboratory metadata group by nested MAP profiling
- observation-bias profile results have JSON write/read helpers preserving group/observation keys, estimates, penalties, corrections, and profile records
- observation-bias profile results and archives retain bias bounds and Gaussian prior mean/sigma
- group observation-bias profiles re-evaluate the final coordinated bias vector before returning, so reported corrections and objective correspond to the reported group biases
- group observation-bias profiles validate that every experiment contains the requested continuous observation before nested refitting
- `Experiment.theoretical_covariance`: externally estimated forward-model uncertainty combined with analytical covariance as $C_D+C_T$ in the weighted likelihood; automatic estimation of $C_T$ remains intentionally out of scope
- `thermoinvert/inversion.py`: deterministic one-parameter MAP-like search
- `thermoinvert/sensitivity.py`: finite-difference sensitivities and SVD rank diagnostics
- `thermoinvert/posterior.py`: Laplace-approximation covariance/correlation diagnostic, brute-force grid posterior, random-walk Metropolis-Hastings sampler (with data/prior uncertainty-source isolation), multi-chain R-hat/effective-sample-size convergence diagnostics, and outlier/acceptance-threshold sample filtering
- `mcmc_sample` and `mcmc_sample_multi_chain` accept an optional positive-definite full proposal covariance for efficient exploration of correlated directions; it tunes sampling only and does not change the posterior prior covariance
- `MCMCResult` and MCMC archives retain the effective proposal covariance, random seed, and burn-in, preventing sampler choices from being lost in a reproducibility record
- `MultiChainMCMCResult` retains the starting points and seed schedule used for R-hat/ESS diagnostics, so apparent agreement can be assessed against chain dispersion
- `MultiChainMCMCResult` and its archive retain consistent aggregate parameter units across chains
- MCMC archive readers validate non-empty parameter units and permit legacy chain files that predate unit metadata
- LOOCV preserves the parent problem's `objective_criterion` and `misfit_form` during refitting and held-out scoring
- `leave_one_group_out_cross_validation` supports study/laboratory metadata holdout with grouped mean held-out misfits and preserved parent likelihood settings
- `LOOCVResult` and its archive retain the grouped metadata key, distinguishing study, laboratory, and other holdout analyses
- grouped LOOCV folds retain individual held-out experiment IDs; older three-column LOOCV archives remain readable
- grouped LOOCV folds retain per-experiment held-out misfits in addition to group means; older archives remain readable
- LOOCV archive readers validate grouped detail-array lengths, IDs, and finite misfits before reconstructing folds
- `plot_loocv_misfits` displays within-group held-out misfit spread when grouped fold details are available
- continuous posterior-predictive checks optionally add declared analytical/model observation noise; the default remains parameter-only prediction
- posterior-predictive checks support reproducible even or random-without-replacement posterior sample selection and retain the strategy in results/archives
- continuous and bracket posterior-predictive results have JSON archive/read helpers for reproducible reporting without rerunning the forward model
- continuous posterior-predictive results retain whether observation noise was included and the seed used to generate predictive replicates
- predictive observation-noise sampling validates covariance key coverage and matrix shape before drawing correlated residuals
- predictive observation-noise sampling uses eigen-factorization, so valid rank-deficient positive-semidefinite covariance is supported
- thermodynamic covariance diagnostics and perturbation sampling reject non-finite matrices and support valid rank-deficient positive-semidefinite matrices
- prior, proposal, and local-resolution covariance inputs consistently reject non-finite matrices before numerical factorization
- joint Gaussian KL/information diagnostics reject non-finite means and covariance matrices before determinant/precision evaluation
- MCMC rejects non-finite initial points and diagonal proposal sigmas before chain generation
- `PTInversionProblem` rejects non-finite and out-of-bound P-T objective inputs before invoking the forward model, matching database-parameter inversion support semantics
- `PTInversionProblem` exposes the shared prior-covariance contract for correlated pressure-temperature priors
- `ScalarParameter` provides native-unit scalar priors for P-T inversion while retaining compatibility aliases for shared posterior diagnostics
- `ScalarParameter` validates names, units, finite priors, positive sigma, and increasing finite explicit bounds at construction
- P-T configuration archives serialize `ScalarParameter` records with native units and retain backward compatibility with legacy thermodynamic parameter records
- Laplace and MCMC posterior results retain parameter units, deriving native units from `ScalarParameter` and defaulting thermodynamic parameters to `J/mol`
- correlated P-T prior covariance is used directly in the PT objective's Mahalanobis prior penalty and is therefore inherited by Laplace, MCMC, grid, and joint information diagnostics
- `write_pt_inversion_configuration` archives P-T-specific parameter order, correlated prior covariance, experiment/brackets, corrections, and likelihood settings
- P-T configuration archives identify `problem_mode` as `composition` or `bracket` while retaining full bracket endpoint geometry
- `read_pt_inversion_configuration` validates and reloads archived P-T provenance without requiring the runtime forward model
- both configuration readers enforce the live constructors' finite, positive-definite prior covariance contract
- archive readers reject duplicate parameter keys and inconsistent P-T problem-mode/bracket metadata
- archive readers validate top-level JSON object/array types before inspecting configuration fields
- archive readers validate experiment/bracket records as JSON objects with unique non-empty IDs
- archive readers validate required experiment P-T/composition/phase fields and bracket lower/upper structure
- archive readers validate parameter keys as non-empty strings and preserve their declared order
- database inversion archive readers validate parameter names/fields, finite prior means, positive sigmas, and increasing optional bounds
- database inversion archive readers enforce correction-key semantics derived from each parameter's phase and field
- database inversion archive readers restrict fields to the supported overlay set `H0`, `G0`, `V0`, `Hf`, and `GH`
- P-T archive readers validate pressure/temperature parameter names and finite positive prior sigmas
- P-T archive readers verify declared bracket IDs match serialized bracket records
- `theoretical_covariance_scale` provides fixed global model-error scaling with $C_D+s_T^2C_T$ in database and P-T objectives; the setting propagates through refits and archives
- `profile_theoretical_covariance_scale` estimates one global theoretical-error scale by nested MAP refits with an explicit Gaussian scale prior
- `InversionProblem` supports validated fixed theoretical-error scales by experiment metadata group, allowing separate study/laboratory $s_{T,g}$ values
- `profile_theoretical_covariance_group_scales` estimates metadata-group scales by coordinate-wise nested MAP refits with shared Gaussian scale priors
- archive readers validate penalties, objective criteria, misfit forms, and quadrature order against live constructor contracts
- P-T archive readers validate correction maps as finite numeric values keyed by non-empty strings
- sensitivity reports now persist baseline observable values, have a strict round-trip reader for Jacobian/SVD diagnostics, and accept the actual correlated prior covariance for prior-scaled information
- `prior_scaled_information` accepts the full correlated prior covariance and whitens the Jacobian before information/rank summaries; the independent-sigma behavior remains the default
- `local_resolution_matrix` accepts optional theoretical covariance and combines it with analytical covariance as $C_D+C_T$, matching the likelihood model
- `write_multi_chain_mcmc_result`/`read_multi_chain_mcmc_result` persist all chains plus aggregate R-hat, ESS, start, and seed provenance
- `thermoinvert/validation.py`: leave-one-out cross-validation and posterior predictive checks
- `thermoinvert/plotting.py`: optional matplotlib visualizations (MCMC trace/scatter, grid posterior density/contour, posterior predictive coverage, LOOCV misfit bars) plus a PDF/PNG export helper
- `thermoinvert/archive.py`: CSV/JSON round-trip persistence for MCMC, LOOCV, and grid-posterior results
- `thermoinvert/input_perturbation.py`: repeated analytical/thermodynamic/all input perturbation with MAP refitting, analogous to mc_fit `invunc` modes
- analytical perturbation honors an experiment's correlated analytical $C_D$ when supplied; theoretical $C_T$ remains a likelihood model-discrepancy term rather than sampled measurement noise
- `thermoinvert/input_perturbation.py` also loads packed Holland-Powell covariance JSON, reports PSD/rank/conditioning diagnostics, and propagates correlated thermodynamic shifts with jittered Cholesky sampling
- `InputPerturbationResult` retains covariance key order and PSD/rank/conditioning diagnostics for reproducible uncertainty provenance
- `archive.py` persists and reloads `InputPerturbationResult` samples, objectives, summary statistics, and covariance provenance
- `write_inversion_configuration`: archives parameter order/bounds, prior covariance, likelihood settings, and experiment/bracket uncertainty records for reproducible result interpretation
- `read_inversion_configuration` validates and reloads database-parameter archive provenance without requiring the runtime forward model
- `docs/analysis/THERMOINVERT_USAGE.md`: operational guide for parameter fields, independent/correlated $C_M$, $C_D$, $C_T$, thermodynamic perturbations, diagnostics, and archival
- `InversionProblem(prior_covariance_j_mol2=...)` now uses a formal positive-definite correlated Gaussian prior in the MAP objective; independent parameter sigmas remain the backward-compatible default
- explicit parameter bounds are enforced by the shared objective and MCMC initialization, giving MAP, grid, and MCMC one consistent truncated-Gaussian prior support
- `grid_posterior` rejects axes outside declared hard bounds and grids with no finite posterior objective, avoiding invalid normalization of off-support evaluations
- `laplace_covariance` is boundary-aware: it adapts interior finite-difference steps and explicitly withholds a symmetric covariance for a MAP at a hard bound
- `diagnose_map_globality`: clusters cold-start MAP endpoints and flags separated near-optimal basins, providing evidence against local uniqueness without claiming a global resolution operator
- `local_resolution_matrix`: explicit covariance-aware Tarantola resolution matrix for a supplied differentiable sensitivity Jacobian; it is deliberately limited to the local linear-Gaussian regime
- `diagnose_posterior_shape`: sample-based screening for skewed or non-normal MCMC marginals, so Laplace summaries are not treated as automatically adequate
- `diagnose_grid_modes`: direct local-maximum detection for 1-D/2-D bounded grid posteriors, including symmetric multimodality missed by moment-based screening
- `grid_information_gain`: numerical non-Gaussian KL information gain on an explicitly evaluated 1-D/2-D grid, retaining the configured correlated prior covariance
- `docs/analysis/THERMOINVERT_ALGORITHM_FOR_SCIENCE_PAPER.md`: paper-ready algorithm description, equations, pseudocode, diagnostics, assumptions, and limitations
- `thermoinvert/io.py`: JSON serialization for `Experiment`/`EquilibriumBracket` datasets, so a literature/experimental dataset can be curated once and reloaded independent of forward-model or optimizer code
- `EquilibriumBracket`: lower/upper endpoint phase-inequality data contract
- `Testing/test_thermoinvert_phase1.py`: Phase 1 tests
- `Testing/test_thermoinvert_integration.py`: end-to-end test combining brackets, condition-uncertainty marginalization, multi-parameter MAP, leave-one-out CV, Laplace covariance, MCMC sampling, and posterior predictive checks in one synthetic scenario
- `Testing/test_thermoinvert_real_forward_model.py`: computational smoke test running the real `ReaktoroForwardModel` and `numerical_sensitivity` against the actual DEW17HP622_Zn_2025 database on a Vazante-like bulk composition (not a scientific calibration claim)
- `Testing/test_thermoinvert_plotting.py`, `Testing/test_thermoinvert_archive.py`, `Testing/test_thermoinvert_io.py`: coverage for the plotting, archival, and dataset-serialization modules
- `Testing/test_thermoinvert_input_perturbation.py`: analytical/thermodynamic input perturbation propagation tests

The current scaffold has now been validated for the following concrete paths:

- temporary JSON overlays close their Windows file descriptors and leave the source database unchanged
- `dHf` proposals apply once to the active Reaktoro `StandardThermoModel.Constant.H0` field by default
- alternative `Hf`, `G0`, and Perple_X `GH` fields require explicit overlay-field selection
- the forward model discovers the repository-local `reaktoro4py` build
- a real Reaktoro equilibrium accepts a temporary zincite enthalpy correction
- the deterministic one-parameter synthetic MAP test passes
- numerical sensitivity diagnostics report singular values, rank, condition number, and weak parameter directions
- sensitivity diagnostics export the Jacobian, singular values, prior-scaled information matrix, and identifiability JSON report
- bounded multi-parameter MAP optimization uses deterministic coordinate-pattern searches from reproducible starts
- bracket likelihood tests score required and forbidden endpoint phases
- optional P-T nuisance uncertainty is marginalized with deterministic Gauss-Hermite quadrature
- `laplace_covariance` computes a numerical-Hessian local covariance/correlation approximation around a candidate MAP point, validated against an analytic quadratic objective
- `grid_posterior` evaluates a normalized brute-force posterior density and marginal moments on a regular 1-2 parameter grid, validated against analytic Gaussian moments
- `mcmc_sample` draws posterior samples with a random-walk Metropolis-Hastings sampler (no scipy/emcee dependency), reporting acceptance rate, posterior mean/std, and 95% credible intervals, validated against analytic Gaussian moments
- `leave_one_out_cross_validation` refits MAP corrections with each experiment held out in turn and scores it out-of-sample, validated on a synthetic multi-experiment recovery case
- `posterior_predictive_check` compares an experiment's observations to the predictive distribution implied by MCMC posterior samples, reporting predictive mean/std, 95% credible intervals, and per-observable coverage
- the full mc_fit-parity plan (Section 17) is complete: generalized multi-field parameter inversion, matplotlib visualizations, outlier/acceptance filtering, data-vs-prior uncertainty decomposition, hot-start multi-parameter MAP, multi-chain R-hat/ESS convergence diagnostics, chi-square/raw-LSQ/Bayes-penalty objective variants, and CSV/JSON archival for MCMC/LOOCV/grid-posterior results
- `Experiment`/`EquilibriumBracket` datasets can be saved to and loaded from JSON independent of any forward model or optimizer code
- `identifiability.identifiability_summary` combines prior-vs-posterior shrinkage with structural SVD weak-direction cross-referencing, per parameter (fixed a related bug: `identify_sensitivity` used the economy SVD, which silently drops null-space directions whenever there are fewer observables than parameters -- exactly the underdetermined case identifiability diagnostics exist to catch; now uses the full SVD)
- `identifiability.information_gain_summary` reports per-parameter KL(posterior || prior) in nats using each parameter's exact Gaussian prior, requiring only one posterior diagnostic result (no separate prior-only run needed)
- `identifiability.multivariate_information_gain` and `posterior_correlation_matrix` report joint covariance-aware information gain and parameter trade-offs, addressing the Tarantola requirement to analyze the posterior in joint model space rather than only by marginal standard deviations
- `run_input_perturbation` performs explicit input-space analytical, thermodynamic, or combined perturbation/refit replicates; unlike objective-term decomposition, this propagates perturbed inputs through repeated inverse solves and is the closer analogue of mc_fit `invunc`

Before extending it, verify:

1. the test import path and test collection configuration
2. that the overlay modifies the correct thermodynamic field exactly once
3. that all selected enthalpy parameters exist and are editable
4. that predictions expose compositions and fractions needed by the likelihood
5. that invalid or failed equilibria are reported separately from valid phase absence

## 16. Immediate Next Actions

1. Make the Phase 1 tests collect and pass.
2. Add an explicit overlay-field contract for `dHf` corrections.
3. Add a synthetic forward model test for one-parameter recovery.
4. Run one real Reaktoro prediction with a temporary correction.
5. Add richer likelihood terms and nuisance-condition handling.
6. Add a Laplace covariance/correlation diagnostic after MAP (done).
7. Add a brute-force grid posterior for 1-2 parameters as a verification tool (done).
8. Add a random-walk Metropolis-Hastings posterior sampler (done).
9. Add leave-one-out cross-validation and posterior predictive checks (done).
10. Add an end-to-end integration test combining all pipeline stages (done).
11. Add a computational smoke test against the real Reaktoro forward model and database (done).
12. Complete the mc_fit-parity plan in Section 17 (done).
13. Add JSON dataset serialization for `Experiment`/`EquilibriumBracket` so literature data can be curated independent of code (done).
14. Curate a real, literature-referenced experimental/natural bracket dataset for genuine scientific calibration.

The present Phase 1 forward result still needs phase-composition extraction,
explicit failed-equilibrium status in the public prediction contract, and a
validated real experimental dataset before it should be used for calibration.

`InversionProblem` can set `marginalize_condition_uncertainty=True` to evaluate
each experiment over deterministic pressure-temperature nuisance nodes. The
quadrature order is controlled by `condition_quadrature_order`; order $q$
uses $q^2$ forward evaluations per experiment. This is a practical first
approximation for sparse data, not a replacement for full nuisance-parameter
posterior sampling.

`laplace_covariance` differentiates the actual `InversionProblem.objective`
(priors, brackets, and penalties included) via a numerical Hessian at a
candidate solution. It reports whether the Hessian is positive definite before
exposing covariance/correlation, and it is explicitly documented as a local
Gaussian approximation rather than a posterior sample, per the guardrails in
Section 11.

`mcmc_sample` is a plain diagonal-proposal random-walk Metropolis-Hastings
sampler over the same `InversionProblem.objective`. It intentionally avoids a
scipy/emcee dependency (none is installed in the `reaktoro` conda
environment) and requires the caller to supply a reasonable starting point
and proposal scale; it does not perform automatic tuning, multiple chains, or
convergence diagnostics (e.g. R-hat), so acceptance rate and trace behaviour
should be checked before trusting reported summaries on real problems.

`leave_one_out_cross_validation` refits `optimize_one_parameter` or
`optimize_multi_parameter` on all experiments except one, then scores the
held-out experiment with `likelihood.experiment_misfit`. This is a generalization
check, not a replacement for domain review of individual brackets.
`posterior_predictive_check` reuses the same observable contract as
`sensitivity.result_observables` (phase presence, fraction, composition) but is
self-contained in `validation.py` so it can compare against an experiment's own
observed values rather than a prediction's.

`Testing/test_thermoinvert_real_forward_model.py` is deliberately restricted to a
small mineral set (`cc`, `mag`, `Znc`) and a single baseline solve plus a
one-parameter numerical Jacobian (3 real equilibrium solves total) so it stays
fast. All other tests use synthetic fake forward models; this is the only test
that exercises genuine Reaktoro equilibrium chemistry end-to-end, and it
intentionally does not assert specific enthalpy values or literature agreement.

The current multi-parameter optimizer is deliberately modest: it is a bounded,
derivative-free coordinate-pattern search intended for small nonsmooth problems.
It records reproducible starts and evaluations, but it is not a posterior sampler
and should be followed by sensitivity, correlation, and predictive-validation
diagnostics.

## 17. Plan of Action: mc_fit Feature Parity and Generalization

A re-read of Perple_X `mc_fit.f`/`mc_fit_ka.f` and its 8 MATLAB plotting
scripts identified functional and visualization gaps relative to this
framework. This section tracks the prioritized plan to close them, and the
generalized parameter architecture needed so the same inversion machinery
can invert more than one standard-state property.

### 17.1 Generalized parameter/field architecture (Phase 0 — done)

Previously every `ThermodynamicParameter` implicitly targeted a species'
`H0` field only, and `DatabaseOverlay.create()` applied one `fields` tuple
uniformly to every correction in a batch — there was no way to invert, say,
`V0` of one species and `H0` of another (or both fields of the same species)
in the same problem.

Implemented:

- `ThermodynamicParameter.field: str = "H0"` — selects which additive
  database record a parameter perturbs (`H0`, `G0`, `V0` are consumed by
  Reaktoro's `StandardThermoModel.Constant`; `Hf`/`GH` are informational
  copies retained for schema completeness only).
- `ThermodynamicParameter.key` property — `phase` for the default `H0`
  field (fully backward compatible), or `"phase:field"` otherwise, so
  multiple fields on the same species coexist in one corrections dict.
- `InversionProblem` uniqueness check now keys on `(phase, field)` pairs
  instead of `phase` alone.
- `DatabaseOverlay.create()` parses compound `"phase:field"` keys directly
  (no more uniform `fields=` override), grouping multiple corrections per
  species before applying them.
- All correction-dict builders (`InversionProblem.objective`,
  `optimize_multi_parameter`, `optimize_one_parameter`,
  `numerical_sensitivity`, `posterior_predictive_check`) now use
  `parameter.key` instead of `parameter.phase`.

Known model limitation (do not invert `S0`): this database's
`StandardThermoModel.Constant` record has no `S0` field, so an entropy
correction has no mechanism to reach the equilibrium solver here — it would
be a silently inert parameter. Margules/solution-model coefficients (w0, wT,
wP) are similarly out of reach until the forward model builds solid-solution
phases instead of pure `MineralPhases`; see 17.4.

### 17.2 Visualization module (done)

`thermoinvert/plotting.py` (matplotlib, optional import — the rest of the
package has no hard plotting dependency; `thermoinvert/__init__.py` falls
back to `None` names if matplotlib is unavailable):

1. `plot_mcmc_trace(mcmc_result)` — per-parameter trace + running mean;
   partially substitutes for missing R-hat/ESS diagnostics.
2. `plot_posterior_scatter(mcmc_result, laplace=None)` — 2-parameter
   scatter plot (or 1D histogram) with an optional Laplace covariance
   ellipse overlay (mc_fit's core "2-panel" plot).
3. `plot_grid_posterior(grid_result)` — 1D density-with-credible-band, or 2D
   filled-contour, for `GridPosteriorResult` (mc_fit's imprecision-band plot).
4. `plot_posterior_predictive_check(ppc_result)` and `plot_loocv_misfits(loocv_result)`
   — observed-vs-predictive-interval coverage plot and a per-experiment
   held-out misfit bar chart (split into two functions instead of one
   combined plot, since LOOCV and PPC operate on different result shapes).
5. `save_figure(fig, path)` — PDF (vector) or PNG (300 DPI) export by extension.

Each function accepts an already-computed result object (no recomputation)
and returns the created `Figure`. Covered by `Testing/test_thermoinvert_plotting.py`
(6 tests, matplotlib `Agg` backend).

### 17.3 Functional gaps, prioritized

1. **Outlier rejection / acceptance filtering (done)** — `filter_outliers_mcmc`
   iteratively drops samples whose objective is more than `z_threshold`
   standard deviations above the current mean (one-sided, recomputed each
   iteration), mirroring mc_fit's robust z-score filter.
   `filter_by_acceptance_threshold` keeps only samples at or below an
   explicit objective cutoff, mirroring mc_fit's `oktol`. Both return a new
   `MCMCResult` with recomputed mean/std/95% credible intervals (shared via
   a `_summarize_samples` helper also used by `mcmc_sample` itself).
2. **Uncertainty source decomposition (done)** — `InversionProblem.objective_components`
   now returns `(data_misfit, prior_penalty)` separately (`objective` is their
   sum, unchanged). `laplace_covariance`/`mcmc_sample` accept an
   `uncertainty_source` switch (`UNCERTAINTY_SOURCES = ("all", "data", "prior")`)
   that restricts exploration to one component, isolating each source's
   contribution to the inferred uncertainty. This decomposes by *objective
   term*, not by perturbing raw analytical/thermodynamic inputs like mc_fit
   does, so results are conceptually related but not numerically comparable
   to mc_fit's component ellipses. `uncertainty_source="prior"` is a useful
   sanity check on its own: it must reproduce the Gaussian prior regardless
   of the forward model (validated against analytic prior mean/std).
   `run_input_perturbation` complements this with explicit input-space
   perturbation/refit replicates for `analytical`, `thermodynamic`, and `all`
   sources, using declared observation sigmas and thermodynamic nuisance
   sigmas. This is the closer analogue of mc_fit `invunc`; correlated input
   covariance and automatic database-error estimation remain future work.
3. **Hot-start mode (done)** — `optimize_multi_parameter` accepts an optional
   `central_point` (e.g. a previous MAP result) plus `hot_start_spread_fraction`
   (default 0.1). When given, every search start is drawn near `central_point`
   (perturbed by that fraction of each parameter's bound width) instead of
   spanning the full prior range; omitting it keeps the previous cold-start
   behavior unchanged. `mcmc_sample`/`mcmc_sample_multi_chain` already take an
   explicit starting point per chain, so no separate hot-start flag was
   needed there -- callers pick the starting point(s) directly.
   `PTInversionProblem.optimize` now accepts the same `central_point` and
   `hot_start_spread_fraction` controls, with backward-compatible cold-start
   behavior when omitted; this is covered by the PT regression suite.
4. **MCMC convergence diagnostics (done)** — `mcmc_sample_multi_chain` runs
   several independent chains from different starting points and returns a
   `MultiChainMCMCResult` with per-parameter Gelman-Rubin `r_hat` (via
   `potential_scale_reduction`) and total `effective_sample_size` (via
   Geyer's initial-positive-sequence ESS estimator, summed across chains).
   At least two chains are required; both diagnostics are validated against
   synthetic well-mixed vs. poorly-mixed/autocorrelated cases.
5. **Alternative objective forms (done)** — `likelihood.experiment_misfit` takes
   a `misfit_form` (`MISFIT_FORMS = ("weighted_lsq", "chi_square", "raw_lsq")`):
   the default `weighted_lsq` is unchanged prior behavior; `chi_square` is
   mc_fit's classic (calculated - observed)² / calculated statistic
   (documented there as "dubious", included for completeness); `raw_lsq` is
   unweighted squared residuals. `InversionProblem(..., misfit_form=...)`
   threads the choice through. Separately, `InversionProblem(...,
   objective_criterion="bayes")` adds a fixed BIC-style
   `bayes_penalty()` (`0.5 * n_parameters * ln(n_experiments + n_brackets)`)
   to `objective()`, as an alternative best-model criterion to plain MAP;
   `posterior._component_objective`'s `"all"` branch now calls
   `problem.objective()` directly so MCMC/Laplace consistently include this
   penalty when configured.
6. **Persistent result archival (done)** — new `thermoinvert/archive.py`:
   `write_mcmc_result`/`read_mcmc_result` (samples.csv + log_posterior.csv +
   summary.json), `write_loocv_result`/`read_loocv_result` (folds.csv with
   corrections embedded as JSON + summary.json), and
   `write_grid_posterior`/`read_grid_posterior` (flattened density/log-posterior
   arrays with shape recorded in summary.json, supporting both 1D and 2D
   grids). All four result types round-trip exactly through disk, validated
   in `Testing/test_thermoinvert_archive.py` (4 tests), so results can be
   reloaded and re-plotted (via `plotting.py`) without rerunning the forward
   model. This completes every item planned in Section 17.3.

### 17.4 Explicitly out of scope for now

- **Margules/solution-model parameter inversion** requires
  `ReaktoroForwardModel._build_system` to construct solid-solution phases
  (not plain `MineralPhases`) — a forward-model architecture change, not
  just a parameter/overlay change. Track separately if a solid-solution
  database becomes available.
- **P-T-as-unknown (`INVPTX`) mode (done)** — `thermoinvert/pt_inversion.py`
   provides `PTInversionProblem`, the mc_fit `INVPTX` counterpart: it perturbs
   an experiment's own P-T (instead of a database parameter) at fixed
   database, scoring candidates via the same `likelihood.experiment_misfit`.
   Classic geothermobarometry (P-T from mineral compositions) uses this class
   directly. It duck-types the same surface as `InversionProblem`
   (`.parameters`, `.objective`, `.objective_components`), so
   `posterior.laplace_covariance`/`mcmc_sample`/`grid_posterior` and
   `identifiability.identifiability_summary`/`information_gain_summary` all
   work on it unchanged -- no code duplication needed for those diagnostics.
   `sensitivity.numerical_pt_sensitivity(problem, point)` is the direct P-T
   counterpart to `numerical_sensitivity`: it central-differences pressure
   and temperature, flattens the same composition/fraction observables, and
   reports the SVD rank, condition number, and weak P-T directions. It is
   intentionally rejected for bracket-only problems because binary
   phase-inequality observations are discontinuous; use posterior/Laplace
   geometry for those. `validation.posterior_predictive_check_pt` replays
   MCMC samples as candidate `[pressure_bar, temperature_c]` conditions and
   returns the same predictive mean/std/95% coverage result used for
   thermodynamic-parameter inversion. It applies to composition/fraction
   P-T observations and intentionally rejects bracket-only problems; the
   new `validation.posterior_predictive_check_pt_brackets` handles that
   bracket-only case separately by reporting the posterior fraction of
   samples satisfying every phase inequality, plus mean and 95% interval of
   the bracket misfit. `plotting.plot_bracket_posterior_predictive` visualizes
   that satisfaction fraction and misfit summary, while
   `plotting.plot_bracket_misfit_distribution` renders the full posterior
   misfit histogram retained by the result object. The existing
   `posterior_predictive_check` remains for database-field
   corrections. A single composition
   constraint generally cannot separate P from T (only their combination is
   identifiable, exactly as in real thermobarometry, which pairs an
   independent barometer + thermometer reaction) -- this was caught by the
   test suite itself (an initial single-observable synthetic model produced a
   near-degenerate MCMC posterior) and is documented in the module's docstring.
   `PTInversionProblem` also accepts `bracket=` or `brackets=`;
   it translates each endpoint geometry around every candidate P-T center and
   scores the endpoint phase inequalities with `bracket_misfit`. This is the
   direct inversion path for pure-phase reaction boundaries. One bracket
   normally constrains a boundary curve rather than a unique P-T point, so
   multiple independent brackets or an additional composition/thermometer
   constraint are required for separate pressure and temperature estimates.
   The PT tests cover both cases: one bracket validates boundary satisfaction,
   while two independent synthetic brackets recover a unique P-T intersection.

## 18. Scientific Guardrails

Always ask:

- Is the parameter identifiable?
- Is the posterior narrower than the prior?
- Is the correction correlated with another parameter?
- Is it physically plausible?
- Does it reproduce independent experiments?
- Does it behave reasonably outside the calibration domain?
- Could the residual come from P-T uncertainty?
- Could it come from an inadequate solution model?
- Could it reflect disequilibrium or kinetics?
- Does the added parameter provide information-supported improvement?

A result such as:

```text
The individual enthalpy is not independently identifiable.
The data constrain delta_H_A - 0.82 * delta_H_B.
```

is a valid scientific outcome, not an optimization failure.
