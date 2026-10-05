For strongly correlated parameters, use a symmetric positive-definite
`proposal_covariance_j_mol2` in `mcmc_sample` or `mcmc_sample_multi_chain`
to align the random walk with known posterior geometry. It is a proposal-tuning
The effective proposal covariance is retained in `MCMCResult` and
`write_mcmc_result` archives it together with the random seed and burn-in for
reproducibility. Posterior results also retain `parameter_units`; thermodynamic
`diagnose_grid_modes(grid_result)` reports separated one- or two-dimensional
grid maxima above a chosen relative-density threshold, including symmetric
bimodality that moment-based MCMC shape screening can miss. Its conclusion is
`grid_information_gain(problem, grid_result)` computes numerical
$D_{KL}(p(\theta|d)\,||\,p(\theta))$ from the actual gridded posterior and
the configured correlated prior, without a Gaussian posterior approximation.
Use it only for 1-D/2-D uniform grids. Check

# Thermoinvert Usage Guide

## Units and Parameter Fields

Use J/mol for thermodynamic corrections and (J/mol)^2 for covariance values.
`ThermodynamicParameter` identifies the fitted record by `phase` and `field`.
`H0` is the default active enthalpy field in the current Reaktoro JSON format;
`G0` and `V0` are also active constant fields. `Hf` and `GH` are retained for
record completeness but do not change Reaktoro equilibrium calculations in this
database format.

For non-thermodynamic scalar variables, use `ScalarParameter` with native
`prior_mean`, `prior_sigma`, bounds, and `unit`. It is currently supported by
P-T inversion, where pressure and temperature are represented in bar and degC
rather than J/mol. Its name and unit must be non-empty, its prior must be
finite with positive sigma, and explicit bounds must be finite and increasing.

```python
parameters = [
    ThermodynamicParameter(
        name="dH_zincite",
        phase="Znc",
        field="H0",
        prior_mean_j_mol=0.0,
        prior_sigma_j_mol=2_000.0,
        lower_bound_j_mol=-10_000.0,
        upper_bound_j_mol=10_000.0,
    ),
]
```

`parameter.key` controls correction-map and covariance order: `Znc` for `H0`,
or `Znc:G0` for another field. A phase/field pair can occur only once.

## Prior Covariance: $C_M$

### Independent

The default uses a diagonal $C_M$ from individual `prior_sigma_j_mol` values.
Use it when no defensible covariance exists between fitted corrections.

```python
problem = InversionProblem(forward_model, experiments, parameters)
```

### Correlated

Pass a positive-definite `prior_covariance` matrix in the exact `parameters` order. It changes the
MAP objective, Laplace approximation, MCMC posterior, and joint KL result.
The same option is available for `PTInversionProblem`, ordered as
`(pressure_prior, temperature_prior)`.

```python
prior_covariance = np.array([[4.0e6, 1.2e6], [1.2e6, 9.0e6]])
problem = InversionProblem(
    forward_model,
    experiments,
    parameters,
    prior_covariance_j_mol2=prior_covariance,
)
```

The covariance is unit-neutral: each entry `(i, j)` has units equal to the
product of parameter `i` and parameter `j` units. This matters for mixed
thermodynamic/P-T vectors, for example J/mol with bar. The legacy
`prior_covariance_j_mol2` keyword remains supported for thermodynamic-only
inversions; use `prior_covariance` for generic or mixed-unit parameters.

Do not silently pseudoinvert a singular Holland--Powell covariance. Subset or
reparameterize it to the fitted parameter set, or make and report a justified
regularization choice. Bounds are independent hard constraints and do not
alter $C_M$ inside the allowed range. `objective`, grid evaluation, MAP, and
MCMC all enforce them: the posterior is the Gaussian prior truncated to the
declared bounds, and an MCMC `initial_point` must lie within them.
All supplied prior and proposal covariance matrices must contain finite values;
non-finite covariance is rejected before eigendecomposition or factorization.
MCMC initial points and diagonal proposal sigmas must also be finite.
P-T inversion uses the same rule: non-finite and out-of-bound pressure or
temperature proposals receive zero posterior support before the forward solve.
Joint Gaussian KL diagnostics likewise reject non-finite means or covariance
matrices before evaluating determinants and precision terms.

## Observation Covariance: $C_D$

### Independent analytical uncertainty

Supply standard deviations for phase fractions or compositions. The likelihood
then sums independent weighted squared residuals.

```python
experiment = Experiment(
    id="run-1",
    pressure_bar=2000.0,
    temperature_c=500.0,
    bulk_composition={"Zn": 1.0},
    observed_phases=("Znc",),
    phase_compositions={"Znc": {"Zn": 0.80}},
    phase_composition_sigma={"Znc": {"Zn": 0.02}},
)
```

### Correlated analytical uncertainty

Use `observation_covariance` with keys for every continuous residual. Key order
must match covariance rows/columns.

```python
experiment = Experiment(
    ...,
    phase_fractions={"Znc": 0.60},
    phase_compositions={"Znc": {"Zn": 0.80}},
    observation_covariance=[[4e-4, 1e-4], [1e-4, 9e-4]],
    observation_covariance_keys=("fraction:Znc", "composition:Znc:Zn"),
)
```

The continuous likelihood is $r^T C_D^{-1}r$. Covariance likelihoods require
`misfit_form="weighted_lsq"`; they cannot be used with `raw_lsq` or
`chi_square`. Phase and bracket constraints are discrete penalties, not
Gaussian covariance observations.

For a known measurement or laboratory offset, use `observation_bias` with the
same continuous residual keys. A positive bias means the reported observation
is high by that amount; the likelihood compares against `observed - bias`.
Bias values must be finite, and unknown keys are rejected.
To estimate one bias per metadata group by nested MAP profiling, use
`profile_observation_bias_by_group(problem, "laboratory", observation_key, ...)`.
It returns group biases, refitted corrections, and profile records with a
Gaussian prior on each bias. This is a MAP profile, not a full hierarchical
MCMC bias model. Archive it with
`write_observation_bias_profile_result(result, "results/bias_profile")`.
The archive also records bias bounds and prior mean/sigma.
Every experiment must contain the profiled observation key; missing keys are
rejected before any nested refits begin.

## Theoretical Covariance: $C_T$

`theoretical_covariance` represents an externally quantified forward-model or
structural error. It must name the same complete residual set as $C_D$, though
it may declare a different key order. The continuous likelihood becomes
$r^T(C_D+C_T)^{-1}r$.
The optional `theoretical_covariance_scale` on `InversionProblem` and
`PTInversionProblem` applies the fixed-scale form $C_D+s_T^2C_T$. It must be
non-negative; `0` disables the supplied $C_T$ and `1` is the default. This is
a fixed scientific choice, not an estimated nuisance posterior. The chosen
scale is archived with the configuration.
For study/laboratory-specific fixed scales, provide
`theoretical_covariance_scale_metadata_key="study"` and a string-to-scale
`theoretical_covariance_scales` mapping. Every experiment carrying $C_T$ must
have a matching metadata group; missing groups are rejected.

To estimate one global scale by nested MAP refits, use
`profile_theoretical_covariance_scale`. It adds a Gaussian prior on the
non-negative scale and returns the selected scale, refitted corrections, and
profile records. This is a MAP profile, not an MCMC sample of model error.
For a joint posterior over corrections and one global model-error scale, use
`mcmc_sample_theoretical_covariance_scale`. The returned
`HierarchicalMCMCResult` includes a bounded, dimensionless `scale_C_T` column,
its Gaussian prior metadata, parameter units, and the usual MCMC summaries.
This sampler is intentionally global-only; metadata-specific scale maps must
use the existing group profiler until a full grouped hierarchy is added. Archive
the result with `write_hierarchical_mcmc_result`.
For one scale per study/laboratory group, use
`profile_theoretical_covariance_group_scales(problem, "study", ...)`. It
returns a scale map, refitted corrections, and group profile records. This is a
coordinate-wise MAP profile. For a joint posterior over corrections and those
group scales, use `mcmc_sample_theoretical_covariance_group_scales` with one
initial scale per metadata group. It returns dimensionless
`scale_C_T[group]` parameters and archives with
`write_grouped_hierarchical_mcmc_result`. The current grouped sampler uses
common bounds and Gaussian prior hyperparameters across groups.
For study or laboratory bias, use `mcmc_sample_group_observation_biases` with
one initial bias per metadata group and a continuous observation key. It uses
the same positive-bias convention as the MAP profiler, samples corrections and
all group biases jointly, and archives with
`write_grouped_observation_bias_mcmc_result`.
For several continuous keys, use
`mcmc_sample_group_observation_biases_multi_key` with initial values keyed by
`(observation_key, metadata_group)`. The resulting bias parameters are labeled
`bias[key|group]` and the archive retains the complete observation-key list.

If `observation_covariance` is absent, independent analytical sigmas construct
diagonal $C_D$. `thermoinvert` does not estimate $C_T$; document its source and
scientific justification in experiment metadata.

## Thermodynamic Input Perturbations

This is separate from $C_M$. `InputPerturbationConfig` perturbs the input
database and/or reported observations, refits repeatedly, and summarizes the
resulting fitted corrections. It does not change the prior density per fit.
When an experiment supplies `observation_covariance`, analytical replicates
draw its phase-fraction/composition observations jointly from scaled $C_D$.
Otherwise they use independent declared sigmas. $C_T$ is deliberately not
sampled here because it describes forward-model discrepancy, not a new
measurement realization.

```python
config = InputPerturbationConfig(
    n_iterations=200,
    source="thermodynamic",  # "analytical", "thermodynamic", or "all"
    thermodynamic_covariance_j_mol2=database_covariance,
    thermodynamic_covariance_keys=("Znc", "Cc"),
)
result = run_input_perturbation(problem, config)
```

These keys identify database entities to perturb and need not be only the
fitted parameters, but must map to valid correction keys. Use
`load_packed_covariance_json(path, scale=1_000)` for Holland--Powell exports
stored in (kJ/mol)^2, then inspect `covariance_diagnostics` before sampling.
Thermodynamic covariance sampling accepts rank-deficient positive-semidefinite
matrices and rejects non-finite or materially negative matrices.

## Fit, Diagnose, Archive

```python
map_result = problem.optimize_multi_parameter(starts=16, iterations=150, seed=7)
globality = diagnose_map_globality(problem, map_result)
mcmc = mcmc_sample(problem, list(map_result.corrections_j_mol.values()))
shape = diagnose_posterior_shape(mcmc)

write_inversion_configuration(problem, "results/configuration")
configuration = read_inversion_configuration("results/configuration")
write_mcmc_result(mcmc, "results/mcmc")
# For multi-chain runs:
write_multi_chain_mcmc_result(multi_chain_result, "results/multi_chain_mcmc")
# For predictive checks:
write_posterior_predictive_result(ppc, "results/ppc")
write_bracket_posterior_predictive_result(bracket_ppc, "results/bracket_ppc")
```

For grouped validation, such as study or laboratory holdout, use:

```python
grouped = leave_one_group_out_cross_validation(problem, "study")
```

For phase-assemblage holdout, use
`leave_one_phase_assemblage_out_cross_validation(problem)`. It derives groups
from each experiment's observed and absent phases, then preserves the held-out
experiment IDs and individual misfits. For an explicit out-of-domain test,
use `leave_one_out_of_domain_pt_cross_validation` with inclusive training
pressure and temperature bounds; experiments outside the rectangle are held
out together and scored against one fit from the in-domain training set.
For database-family holdout, assign labels such as `database_family` in
experiment metadata and call `leave_one_database_family_out_cross_validation`.
Alternatively, provide `database_path`; path-only records are grouped by the
database file's SHA-256. `database_file_provenance(path)` returns the resolved
path, size, and SHA-256, while `validate_database_provenance` checks recorded
hashes against the files. For repeated study/laboratory validation, use
`repeated_group_holdout_cross_validation(problem, "study", repeats=20,
test_fraction=0.25, seed=7)`. It holds out complete metadata groups. Bracket
endpoints must be assignable to groups or the fold is rejected to avoid leakage.

Every experiment must contain the requested metadata key, and at least two
groups must exist. The returned folds use group IDs and preserve the parent
objective and likelihood settings. The resulting `LOOCVResult` records the
metadata key used, and its archive preserves that provenance.
Each grouped fold also records all held-out experiment IDs.
Grouped folds also retain individual held-out misfits, not only their mean, so
study or laboratory-specific failures remain visible.
`plot_loocv_misfits` displays their within-group sample spread as error bars
when detailed grouped misfits are available.
`repeated_grouped_pt_calibration_validation` also reports per-held-out-
experiment `held_out_log_predictive_densities` and
`mean_log_predictive_density`. These are log-mean predictive densities over
posterior draws, using normalized Gaussian, logit-normal, ALR logistic-normal,
or Dirichlet densities as selected by `misfit_form`. Observed and absent phases
use a Bernoulli error model with `phase_error_probability` (default 0.01); choose
this value to reflect phase-identification reliability. The score is not
available for raw least squares, chi-square, unnormalized ALR, or bracket-only
constraints. In those cases validation retains its existing misfit results and
leaves the predictive-density fields empty rather than treating an unnormalized
objective as a probability.

Use cold starts to look for competing basins. Multiple separated near-optimal
basins invalidate a global interpretation of local Laplace covariance or
resolution. For 1-2 parameters, use `grid_posterior`; otherwise use multiple
MCMC chains with R-hat, ESS, traces, correlation, and shape checks.
`MultiChainMCMCResult` retains each initial point and its seed schedule. Use
over-dispersed, recorded starts: an R-hat value is not evidence of global
mixing when all chains were initialized in the same local basin.
For strongly correlated parameters, use a symmetric positive-definite
`proposal_covariance_j_mol2` in `mcmc_sample` or `mcmc_sample_multi_chain` to
align the random walk with known posterior geometry. It is a proposal-tuning
input only: it does not change $C_M$ or the target posterior. Do not also pass
`proposal_sigma_j_mol`; the latter remains the independent-proposal default.
The effective proposal covariance is retained in `MCMCResult` and
`write_mcmc_result` archives it together with the random seed and burn-in for
reproducibility. Posterior results also retain `parameter_units`; thermodynamic
parameters default to `J/mol`, while native-unit scalar parameters retain their
declared unit. Archive readers validate unit count and non-empty string values.
`diagnose_grid_modes(grid_result)` reports separated one- or two-dimensional
grid maxima above a chosen relative-density threshold, including symmetric
bimodality that moment-based MCMC shape screening can miss. Its conclusion is
limited to the supplied bounds and grid resolution.
`grid_information_gain(problem, grid_result)` computes numerical
$D_{KL}(p(\theta|d)\,||\,p(\theta))$ from the actual gridded posterior and
the configured correlated prior, without a Gaussian posterior approximation.
Use it only for 1-D/2-D uniform grids. Check
`grid_covers_parameter_bounds`: a false value means the result is conditional
on the selected grid window rather than the complete bounded prior support.
Grid axes must still lie inside the declared hard parameter bounds; use a
larger declared bound when a wider exploration domain is scientifically
justified.
For disconnected or phase-boundary posteriors, use
`parallel_tempering_sample` with over-dispersed initial points and an increasing
temperature ladder beginning at `1.0`. It returns cold-chain samples, adjacent
swap acceptance rates, mode occupancy, mode-specific means and best objectives,
and a `crossed_modes` flag. Mode labels are a distance-based diagnostic in the
declared bounded parameter domain, not a proof of complete global exploration;
inspect ladder swaps and repeat with different seeds. Use
`analyze_parallel_tempering_modes(result)` for shared basin identities across
temperatures. Its mode probabilities are cold-chain occupancy estimates, not
marginal likelihood or Bayes-factor estimates. The archive retains per-
temperature samples and shared labels. Archive with
`write_parallel_tempering_result`.
The sampler also evaluates the objective at `uniform_reference_samples` draws
from the declared bounded parameter box (default 128). Call
`thermodynamic_integration(result)` to estimate
`log E_uniform[exp(-objective)]` using the inverse-temperature ladder. This is
an objective-target normalizer relative to uniform measure on those bounds,
not a Bayes factor: the inversion objective does not include all normalized
likelihood and prior constants. Accuracy depends on chain mixing, the number
of uniform reference draws, and ladder resolution; check stability by
increasing all three and repeating with independent seeds. The result includes
a moving-block bootstrap standard error and 95% Monte Carlo interval that
resample the retained ladder trajectory and uniform-reference scores. This
interval reflects Monte Carlo sampling uncertainty only; it does not capture
temperature-ladder quadrature error or heuristic mode-clustering uncertainty.
The report also compares trapezoidal integration with quadratic panels over
the irregular inverse-temperature ladder and gives the maximum ladder step;
their difference is a quadrature sensitivity diagnostic, not a rigorous error
bound. `analyze_parallel_tempering_modes` reports mode counts and cold-chain
occupancy over nearby clustering thresholds, plus per-mode occupancy ranges.
Those ranges show sensitivity to the heuristic separation threshold, not
posterior probability intervals.
Continuous PPCs default to parameter-only predictions. Set
`include_observation_uncertainty=True` in `posterior_predictive_check` or
`posterior_predictive_check_pt` to add independent analytical and theoretical
observation noise using the experiment's sigmas and $C_D+C_T$; pass `seed` for
reproducible predictive replicates. The result records both choices. Phase-
presence indicators remain discrete.
PPCs default to evenly spaced posterior samples. Set
`sample_selection="random"` with a recorded `seed` when separated modes or
short-lived posterior regions need unbiased sample selection; the result
records the strategy.
When covariance is supplied, its keys must identify exactly every continuous
fraction/composition observable in the predictive check; mismatches are
rejected rather than silently dropping an observable.
Positive semidefinite, rank-deficient covariance is supported in predictive
sampling; only materially negative eigenvalues are rejected.
`laplace_covariance` adapts finite-difference steps near a bound, but returns
`boundary_limited=True` and no covariance for an exact-bound MAP. In that
case, use bounded-grid or MCMC quantiles rather than a symmetric Gaussian
uncertainty ellipse.
`local_resolution_matrix` is valid only for a differentiable local Jacobian.
Pass `observation_covariance` and optional `theoretical_covariance` to use the
same combined $C_D+C_T$ assumed by the likelihood; never apply it across
phase-boundary discontinuities.
To summarize resolution over an MCMC posterior without building Jacobians
yourself, use `posterior_resolution_from_mcmc` with either an `MCMCResult` or a
`MultiChainMCMCResult`. It evaluates experiment observables at evenly spaced
posterior draws (at most 100 by default), using central differences in the
interior and one-sided differences at parameter bounds. For multiple chains,
the sample budget is balanced across chains and must allow at least one draw
per chain. Any declared units on the posterior and individual chains must match
the inversion parameter order; legacy results with unspecified units remain
accepted. It requires a stable observable layout across the selected draws and
finite-difference evaluations, and excludes bracket-only problems. Supply
`observation_covariance` and `theoretical_covariance` only when overriding the
problem's experiment-level uncertainty model. Omitted matrices are assembled
as block-diagonal continuous-observation covariances from each experiment's
keyed $C_D$, analytical sigmas, and keyed $C_T$, including configured global or
metadata-group theoretical scales. Observations without a declared positive
analytical uncertainty are excluded when no keyed covariance covers them,
matching their zero contribution to the weighted likelihood. Explicit matrices
must use this flattened continuous-observation order: experiments in problem
order, fractions by sorted phase, then compositions by sorted phase and
component. Discrete phase-presence indicators are excluded from the Jacobian.
Invalid equilibrium solves and any stable-phase assemblage change across
selected draws or finite-difference perturbations are rejected, because the
resulting local derivative would not describe one smooth equilibrium branch.

```python
resolution = posterior_resolution_from_mcmc(
    problem,
    mcmc,
    observation_covariance=observation_covariance,
    theoretical_covariance=theoretical_covariance,
    max_samples=100,
)
write_posterior_resolution_result(resolution, "results/posterior_resolution")
plot_posterior_resolution_summary(resolution)
```

The summary contains the mean local resolution matrix and elementwise
credible bounds; these are descriptive intervals, not simultaneous matrix
credible regions. `plot_posterior_resolution_summary` displays the mean and
elementwise interval width. MCMC-derived results retain the selected chain and
draw indices, observable row names, finite-difference step fraction, requested
credible interval, and selection strategy; these are preserved by
`write_posterior_resolution_result`. The report also records whether analytical
and theoretical covariance came from explicit matrices or problem-assembled
experiment metadata; when present, the effective $C_D$ and $C_T$ matrices are
also archived, together with the effective prior covariance $C_M$. Never
interpret a posterior average as a global
resolution result when phase-boundary discontinuities are present.
Sensitivity diagnostics can be archived and reloaded without a forward-model
rerun:

```python
write_sensitivity_report(
    sensitivity,
    parameters,
    "results/sensitivity",
    prior_covariance=problem.prior_covariance_j_mol2,
)
restored = read_sensitivity_report("results/sensitivity")
```

For P-T inversions, archive the P-T-specific configuration separately:

```python
write_pt_inversion_configuration(problem, "results/pt_configuration")
configuration = read_pt_inversion_configuration("results/pt_configuration")
```

To infer thermodynamic corrections and sample pressure/temperature jointly,
use `JointInversionProblem`. Its parameter order is thermodynamic corrections,
pressure, then temperature; the latter use native `bar` and `degC` units. The
problem supports correlated priors, bounded MAP optimization, and the existing
MCMC/Laplace routines through the shared objective contract. Archive its
parameter order, units, covariance, experiment, and likelihood settings with
`write_joint_inversion_configuration`.

To include pressure and temperature calibration offsets, wrap a
`PTInversionProblem` in `PTCalibrationBiasProblem` with pressure- and
temperature-bias `Parameter` priors. The vector is reported pressure,
reported temperature, pressure offset, temperature offset; the forward model
uses reported conditions plus the corresponding inferred offsets. The wrapper
preserves correlated P-T priors and supports MAP, MCMC, and Laplace workflows.
Use `write_pt_inversion_configuration` for either the base P-T problem or the
calibration-bias wrapper; the archive records the four-parameter order, units,
and full covariance and remains distinguishable from legacy two-parameter
archives.

For study/laboratory-level calibration uncertainty across multiple experiments,
use `GroupedPTCalibrationProblem(forward_model, experiments, "study",
pressure_scale_hyperprior=..., temperature_scale_hyperprior=...)`. It samples
one pressure and temperature offset per group plus shared pressure/temperature
scale parameters. Group offsets are conditionally Normal, truncated at five
times their configured scale hyperprior, and each scale has a bounded
half-Normal hyperprior. The scale-dependent truncated-Normal normalization is
included in the posterior and held-out predictive offsets use the same bounds.
Optionally pass `thermodynamic_parameters` and
`thermodynamic_prior_covariance` to jointly infer database corrections with
calibration offsets; inferred correction keys must not also appear in the
fixed `corrections` mapping. Their parameter block comes first, followed by
all pressure offsets, all temperature offsets, then the two scales. Sample it
with `mcmc_sample` and archive its
configuration with `write_grouped_pt_calibration_configuration`; MCMC samples
can use the ordinary `write_mcmc_result` archive.

The same posterior can jointly include study-level observation biases and
theoretical-error scales. Pass `observation_bias_keys=("fraction:phase",)`
with `observation_bias_prior_sigma` and `observation_bias_bounds` to fit an
additive bias for that continuous observation in every group. Pass
`model_error_scale_prior_sigma` and `model_error_scale_bounds` to fit one
positive $s_T$ per group with theoretical covariance. The likelihood uses
$C_D+s_T^2C_T$ for each such experiment. The full parameter order is database
corrections, P offsets, T offsets, observation biases, C_T scales, then the two
P-T offset scales. For composition-only Dirichlet fits, add
`DirichletConcentrationParameter` values after database corrections; arbitrary
covariance is not a Dirichlet likelihood and should use
`logistic_normal_nll`. The configuration archive stores all enabled nuisance
blocks and their order.
For a bounded deterministic starting fit, call `problem.optimize(starts=8,
iterations=120, seed=7)`. It returns the best `MAPResult` and per-start records;
use this point to initialize chains only after checking that the optimization
is stable to different starts.

For convergence diagnostics across independent chains, use
`mcmc_sample_grouped_pt_calibration(problem, n_chains=4, n_samples=5000,
burn_in=500, seed=7)`. It generates reproducible in-bounds starts when none are
provided, draws correlated thermodynamic corrections jointly from their
configured prior covariance, and returns the common `MultiChainMCMCResult`
with per-parameter R-hat and total ESS. Inspect those diagnostics before
interpreting scale or thermodynamic-correction marginals; repeated grouped
validation uses this multi-chain sampler for its training hyperposterior as
well.

Use `repeated_grouped_pt_calibration_validation` to refit the hyperposterior
over repeated complete-group splits. Held-out group offsets are drawn from the
posterior scale distribution, producing predictive misfits rather than a
zero-offset plug-in score. Unseen observation biases and C_T scales are drawn
from their specified priors because those groups have no training posterior
draws. Record the seed and increase chain length and predictive draws to assess
stability.

It records the pressure/temperature parameter order, correlated $C_M$, bounds,
misfit settings, corrections, experiment, brackets, and whether the run is
composition-based or bracket-based. Complete bracket endpoint geometry is
stored in `brackets.json`.
When native-unit `ScalarParameter` objects are used, their `unit`, native prior,
and native bounds are preserved in the parameter records.
Both configuration readers validate that archived $C_M$ is finite and positive
definite before returning provenance. They also reject duplicate database
parameter keys and inconsistent P-T mode/bracket metadata.
P-T bracket IDs must also match the serialized bracket records exactly.
Archived experiments require `id`, finite pressure/temperature, bulk
composition, and observed-phase fields; archived brackets require lower and
upper experiment objects.
Missing required configuration fields are rejected with an explicit validation
error, as are wrong JSON container types for configuration, experiment, or
bracket files. Record objects must also have unique, non-empty string IDs.
Parameter keys must likewise be non-empty strings in the archived order.
Database-parameter records additionally require non-empty names/phases/fields,
finite positive prior sigmas, finite prior means, and increasing finite bounds
when bounds are specified. Their archived key must match the phase/field rule:
bare phase for `H0`, or `phase:field` otherwise.
Supported archived fields are `H0`, `G0`, `V0`, `Hf`, and `GH`.
P-T archives additionally require the expected pressure/temperature parameter
names and finite positive prior sigmas.
Objective settings are validated against the live constructors, including
non-negative penalties, valid criterion/misfit choices, and positive quadrature
order.
`logit_lsq` is available for independent fraction/composition observations in
(0, 1). It applies a logit transform with delta-method sigma and rejects
covariance likelihoods; use `weighted_lsq` for correlated $C_D/C_T$ data.
For closed mineral compositions, `alr_lsq` applies an additive-log-ratio
transform using the final sorted component as the reference. It supports
covariance-aware analytical and theoretical errors through the delta-method
Jacobian, requires strictly positive composition components, and rejects phase
fractions or single-component compositions. It is a composition-only
likelihood; use `weighted_lsq` for mixed fraction/composition observations.
For a normalized covariance-aware density in ALR coordinates, use
`logistic_normal_nll`. It includes the transformed covariance log determinant
and Gaussian normalization and supports analytical and theoretical covariance.
For a simplex likelihood using declared marginal composition sigmas, use
`dirichlet_nll`. It normalizes each phase composition, infers a concentration
from the predicted composition and component sigmas, and evaluates the
Dirichlet negative log likelihood. It requires multi-component compositions
and rejects covariance matrices; use `alr_lsq` or `logistic_normal_nll` when
transformed covariance is required.
To specify the Dirichlet precision directly, set
`phase_composition_concentration={"phase": kappa}` on `Experiment`. The value
is the total concentration (effective count) for that phase; component alpha
parameters are `kappa` times the normalized predicted proportions. Explicit
concentrations do not require marginal composition sigmas and are preserved in
experiment JSON and configuration archives. Phases without an explicit value
retain the legacy sigma-based concentration inference.
To infer concentration jointly with thermodynamic corrections, add
`DirichletConcentrationParameter("phase", prior_mean=..., prior_sigma=...,
lower_bound=..., upper_bound=...)` to the `InversionProblem` parameter list.
It participates in the same prior covariance and MAP/MCMC parameter vector but
is kept out of the Reaktoro database correction map. A fitted concentration
overrides a fixed experiment value for that phase.
For below-detection components, include those components with observed value
zero and set `phase_composition_detection_limit={"phase": limit}`. The limit
is for the total undetected fraction in that phase; the likelihood combines a
Beta-CDF probability for censored mass with a Dirichlet likelihood for the
conditional detected-component proportions. It does not add pseudocounts.
Zeros without an explicit detection limit are rejected by `dirichlet_nll`;
this API treats zero-marked values as censored, not as structural zeros.
P-T correction maps are also required to use non-empty string keys and finite
numeric values.

The report preserves observable values, the Jacobian, SVD rank and condition
diagnostics, weak directions, and prior-scaled information.
The JSON metadata also records the exact $C_M$ used for the prior-scaled
matrix, or `null` when independent marginal sigmas were used.
Pass `prior_covariance=problem.prior_covariance_j_mol2` to
`prior_scaled_information` when using a correlated prior; the Jacobian is then
whitened with the Cholesky factor of $C_M$ rather than scaled by independent
marginal sigmas.

The configuration archive stores parameter order, bounds, $C_M$, likelihood
settings, experiments, brackets, $C_D$, and $C_T$, and records the forward-model
implementation automatically. For `ReaktoroForwardModel`, this includes the
resolved database path, file size and SHA-256, engine module/version when
exposed by the binding, species lists, bulk-species mapping, activity model,
and phase-stability threshold. The database itself is not copied; preserve the
file and verify its hash with `database_file_provenance` when reproducing a run.
Custom forward models are identified by their qualified class name, but their
internal settings must still be archived by the application that defines them.
MCMC archives separately retain samples, objective values, summary statistics,
uncertainty-source choice, proposal covariance, random seed, and burn-in; save
the corresponding inversion configuration beside a result archive.
`write_multi_chain_mcmc_result` archives every chain in its own subdirectory
and records aggregate R-hat/ESS, initial points, and seed schedule in
`multi_chain_summary.json`, including aggregate `parameter_units`.
`write_posterior_predictive_result` and
`write_bracket_posterior_predictive_result` archive predictive intervals,
coverage flags, bracket satisfaction, and individual bracket misfits without
requiring a later forward-model rerun.

Current implementation boundary: global theoretical-error scale sampling,
parallel tempering, joint thermodynamic/P-T inversion, grouped and repeated
validation, database-file provenance checks, logistic-normal and Dirichlet
composition likelihoods, and P-T calibration-bias inversion are implemented.
Mode probabilities are cold-chain occupancy estimates, not evidence
calculations. Group-specific nuisance uncertainty beyond the implemented
samplers and profiled alternatives remains a modeling choice to document.