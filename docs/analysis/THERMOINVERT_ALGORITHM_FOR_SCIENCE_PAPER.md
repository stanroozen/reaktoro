# Thermoinvert: Bayesian Thermodynamic and P-T Inversion

## Abstract

`thermoinvert` is a Reaktoro-based framework for estimating thermodynamic
parameter corrections or pressure-temperature conditions from sparse
phase-equilibrium observations. The method combines equilibrium calculations,
explicit observational likelihoods, bounded priors, derivative-free MAP
optimization, local Laplace approximations, grid verification, and
Metropolis-Hastings posterior sampling. The same statistical machinery supports
classic geothermobarometry, in which pressure and temperature are inferred
from mineral compositions, and pure-phase reaction-boundary inversion, in
which phase inequalities constrain P-T conditions.

Practical configuration examples, including independent and correlated prior,
analytical, and theoretical covariance choices, are in
`THERMOINVERT_USAGE.md`.

The method is deliberately designed to distinguish three different claims:

1. whether a model fits the observations;
2. whether a parameter is constrained relative to its prior; and
3. whether the parameter is structurally identifiable from the experiment
   design and forward model.

## 1. Problem Definition

Let $\theta$ denote the unknown parameter vector. Depending on the problem,
$\theta$ can represent:

- additive corrections to thermodynamic database fields, e.g. $\delta H_i$;
- pressure and temperature, $\theta = (P,T)$;
- multiple P-T boundary coordinates constrained by reaction brackets.

For a database correction, the perturbed thermodynamic record is

$$
H_i^{new} = H_i^{prior} + \delta H_i.
$$

The forward model computes an equilibrium prediction

$$
\hat{y} = f(\theta; d),
$$

where $d$ contains the experiment pressure, temperature, bulk composition,
phase assemblage, and database. `ReaktoroForwardModel` creates an isolated
thermodynamic database overlay for each proposal, solves equilibrium, and
returns validity, stable phases, phase fractions, and phase compositions.

For P-T inversion, the database is fixed and the candidate conditions replace
the experiment's reported pressure and temperature before the equilibrium
solve.

## 2. Observations and Likelihood

An `Experiment` can contain:

- reported pressure and temperature, optionally with uncertainties;
- bulk composition;
- stable phases and absent phases;
- phase fractions and standard deviations;
- mineral or endmember compositions and standard deviations.

The default experiment misfit is a weighted least-squares statistic:

$$
\chi^2_{data} =
\sum_i \left(\frac{\hat{y}_i-y_i}{\sigma_i}\right)^2
+ \Pi_{phase},
$$

where $\Pi_{phase}$ penalizes missing, extra, or forbidden phases.

When correlated analytical uncertainty is available, an experiment can supply
an explicitly keyed positive-definite `observation_covariance`. The continuous
residual vector is then evaluated as $r^T C_D^{-1}r$, while phase-presence and
bracket inequality penalties remain separate non-Gaussian terms.

An explicitly keyed `theoretical_covariance` represents externally estimated
forward-model or structural uncertainty $C_T$. For continuous observations the
Gaussian covariance is $C_D+C_T$. The implementation does not infer $C_T$;
its provenance and adequacy remain a scientific-modeling responsibility.

Alternative objective forms are available:

- `weighted_lsq`: uncertainty-weighted squared residuals;
- `logit_lsq`: independent bounded composition/fraction residuals in logit
  space, with delta-method uncertainty;
- `chi_square`: $(\hat{y}-y)^2/\hat{y}$, included for compatibility with
  classic `mc_fit` workflows and documented as numerically fragile;
- `raw_lsq`: unweighted squared residuals.

  `logit_lsq` requires calculated and observed continuous values strictly inside
  $(0,1)$ and is not combined with covariance matrices. Correlated observations
  should use `weighted_lsq` with explicit $C_D$ and, when appropriate, $C_T$.

An `EquilibriumBracket` represents a phase transition bounded by two endpoint
experiments. Its likelihood is an inequality score:

$$
\Pi_{bracket} = w_{phase}\,N_{violations},
$$

where required phases must be present and forbidden phases must be absent at
the appropriate endpoint.

## 3. Prior and MAP Objective

Each thermodynamic or P-T parameter has a Gaussian prior and explicit bounds:

$$
\theta_j \sim \mathcal{N}(\mu_j,\sigma_j^2),
\qquad
\ell_j \leq \theta_j \leq u_j.
$$

The prior is evaluated on this bounded support: proposals outside the declared
bounds have zero posterior density. Thus the practical prior is a truncated
Gaussian even though the normalizing constant is omitted because it is fixed
for a given model.

The default negative log-posterior objective is

$$
\Phi(\theta) = \chi^2_{data}(\theta)
+ \frac{1}{2}\sum_j
\left(\frac{\theta_j-\mu_j}{\sigma_j}\right)^2.
$$

An optional BIC-style model-complexity term is available:

$$
\Phi_{Bayes}(\theta) = \Phi(\theta)
+ \frac{k}{2}\ln(n),
$$

where $k$ is the number of free parameters and $n$ is the number of
experiments plus brackets. This term is constant with respect to $\theta$ for
a fixed problem and therefore changes model comparison rather than the MAP
location within one fixed model.

When a thermodynamic covariance matrix is available, `InversionProblem` can
use it directly in the prior quadratic form:

$$
\Phi_{prior}=\frac{1}{2}(\theta-\mu)^T C_M^{-1}(\theta-\mu).
$$

Without an explicit covariance matrix, $C_M$ is diagonal and is constructed
from the individual parameter standard deviations. Correlated Holland--Powell
covariance in input perturbation propagates database uncertainty, whereas
`prior_covariance_j_mol2` changes the actual prior density used by MAP, Laplace,
and MCMC.

`InversionProblem.objective_components` returns the data misfit and prior
penalty separately. This supports both source diagnostics and explicit
objective decomposition.

## 4. Optimization

The primary optimizer is a bounded coordinate-pattern search:

1. construct the prior-mean start and additional random starts;
2. evaluate positive and negative coordinate steps;
3. accept an improving coordinate move;
4. halve all steps when no coordinate improves the objective;
5. stop after the requested iterations or when steps become negligible;
6. retain the best result and every-start diagnostic record.

The optimizer is derivative-free because phase stability can make the objective
nonsmooth. A hot-start mode restricts randomized starts to a neighborhood of a
user-supplied central estimate.

For one parameter, a bounded golden-section search is available.

## 5. P-T Inversion Modes

### 5.1 Mineral-composition geothermobarometry

`PTInversionProblem` replaces the candidate experiment's $P$ and $T$ with the
proposed values and scores predicted mineral compositions or fractions. The
same MAP, Laplace, grid, and MCMC functions operate on this problem through a
small common interface:

- `.parameters`;
- `.objective(values)`;
- `.objective_components(values)`.

A single mineral composition may constrain only a combination of P and T. An
independent barometer and thermometer, or multiple independent compositions,
are generally required for separate P and T estimates.

### 5.2 Pure-phase reaction boundaries

`PTInversionProblem` can also accept one or more `EquilibriumBracket` objects.
For each candidate P-T center, each bracket's endpoint geometry is translated
around that center and the endpoint phase inequalities are evaluated.

One bracket generally defines a curve in P-T space. Multiple independent
brackets can produce a unique intersection. Binary phase inequalities are
discontinuous, so bracket problems use posterior/Laplace geometry and bracket
posterior-predictive satisfaction rather than ordinary finite-difference
composition sensitivity.

## 6. Uncertainty Quantification

### 6.1 Laplace approximation

A central finite-difference Hessian is evaluated around a candidate MAP point:

$$
H_{ii} \approx \frac{\Phi(\theta+h_i e_i)-2\Phi(\theta)
+\Phi(\theta-h_i e_i)}{h_i^2},
$$

with analogous mixed derivatives for off-diagonal terms. If $H$ is positive
definite, the local covariance is

$$
\Sigma \approx H^{-1}.
$$

The correlation matrix is derived from $\Sigma$.

Finite differences are reduced near a hard bound to remain in posterior
support. At an exact-bound MAP, no symmetric local Gaussian approximation is
reported; bounded-grid or MCMC summaries are required instead.

### 6.2 Grid posterior

For one or two parameters, the objective is evaluated on a regular bounded grid.
The grid must be contained within the declared hard parameter support; grids
with no finite posterior objective are rejected rather than normalized.
The normalized density is

$$
q(\theta) =
\frac{\exp[-\Phi(\theta)+\Phi_{min}]}
{\int \exp[-\Phi(\theta)+\Phi_{min}]\,d\theta}.
$$

The grid result provides a direct verification of MAP and Laplace summaries,
especially for asymmetric or multimodal posteriors.

### 6.3 Metropolis-Hastings posterior

The baseline sampler uses a random-walk Gaussian proposal:

$$
\theta' = \theta + \epsilon,
\qquad
\epsilon_j \sim \mathcal{N}(0,s_j^2).
$$

An explicitly supplied positive-definite proposal covariance may replace the
diagonal proposal to improve exploration of correlated parameter directions;
it changes sampler efficiency but not the target posterior.

The proposal is accepted with probability

$$
\alpha = \min\{1,\exp[-\Phi(\theta')+\Phi(\theta)]\}.
$$

After burn-in, the sampler reports samples, acceptance rate, posterior means,
standard deviations, and 95% credible intervals. Multiple chains report
Gelman-Rubin $\hat{R}$ and autocorrelation-based effective sample size. Their
initial points and seed schedule are retained so convergence claims can be
audited against the required over-dispersed start design.

`diagnose_posterior_shape` screens MCMC marginals for skewness, excess
kurtosis, and standardized mean--median separation. A flagged marginal is
evidence against a Laplace/Gaussian summary and should be reported with sample
quantiles and posterior plots. This screening does not prove multimodality or
detect all joint non-Gaussian structure.

For explicitly evaluated one- or two-parameter grids, `diagnose_grid_modes`
reports local posterior maxima above a relative-density threshold. Multiple
reported modes give direct bounded-grid evidence of non-Gaussian global
structure, including symmetric bimodality that moment diagnostics can miss.

`grid_information_gain` directly integrates $D_{KL}(p(m|d)\,||\,p(m))$ over
such a grid, retaining non-Gaussian posterior shape and correlated $C_M$.
When the grid does not span the declared bounds, it is reported as conditional
on the selected grid domain rather than as full bounded-prior information.

## 7. Input-Space Uncertainty Propagation

The `run_input_perturbation` workflow is the closest implementation analogue
to Perple_X `mc_fit`'s `invunc` modes. For each replicate:

1. perturb analytical inputs using declared P-T, phase-fraction, and
   composition standard deviations;
2. optionally perturb thermodynamic nuisance corrections using user-supplied
   standard deviations;
3. construct an isolated perturbed inversion problem;
4. rerun MAP optimization;
5. store the fitted parameter vector and objective.

The source selector is:

- `analytical`: perturb observations only;
- `thermodynamic`: perturb thermodynamic nuisance corrections only;
- `all`: perturb both sources.

This is distinct from `uncertainty_source="data"`/`"prior"` in posterior
sampling. The latter separates objective terms; `run_input_perturbation`
propagates perturbations through the input space and reruns the inverse solve.

### 7.1 Correlated thermodynamic uncertainty

Thermodynamic database uncertainties may be correlated. The perturbation
configuration therefore accepts either independent standard deviations or a
full covariance matrix $\Sigma_{thermo}$ with an explicit ordered key list.
Correlated shifts are generated as

$$
\delta g = Lz, \qquad LL^T = \Sigma_{thermo}, \qquad z\sim N(0,I),
$$

using a numerically stabilized Cholesky factorization. The implementation also
loads the packed `Entities`/`PackedUpperTriangle` JSON format used by the
Holland--Powell covariance exports. The loader supports unit scaling, and
covariance diagnostics report positive semidefiniteness, eigenvalue range,
effective rank, and condition number.

This preserves shared database-reference errors and cross-entity correlations
through every perturbed equilibrium calculation. Replacing $\Sigma_{thermo}$
with only its diagonal is a different, less informative uncertainty model.
`InputPerturbationResult` retains the covariance key order and PSD/rank/
conditioning diagnostics so the uncertainty model used for a reported result
can be archived and audited.

## 8. Identifiability Diagnostics

### 8.1 Sensitivity SVD

For a numerical Jacobian $J$, compute

$$
J=USV^T.
$$

The singular values describe constrained directions. The rank and condition
number quantify structural identifiability. Rows of $V^T$ associated with the
null space are reported as weak directions. Full SVD is required when the
number of observables is smaller than the number of parameters.

`numerical_pt_sensitivity` provides the same calculation for P-T composition
observations.

### 8.2 Prior-to-posterior shrinkage

For each parameter,

$$
S_j = 1 - \frac{\operatorname{Var}(\theta_j|data)}
{\operatorname{Var}(\theta_j|prior)}.
$$

Near zero indicates prior domination; near one indicates strong data-driven
constraining. This is combined with SVD weak-direction flags in
`identifiability_summary`.

### 8.3 Information gain

For Gaussian marginal approximations, information gain is the KL divergence
from posterior to prior:

$$
D_{KL}(N_p\,||\,N_0)=\frac{1}{2}\left[
\ln\frac{\sigma_0^2}{\sigma_p^2}
+\frac{\sigma_p^2+(\mu_p-\mu_0)^2}{\sigma_0^2}-1\right].
$$

It is reported in nats by `information_gain_summary`.

### 8.4 Joint information and parameter trade-offs

Tarantola's formulation treats the inverse problem as a probability
distribution over the joint model space, not as a collection of unrelated
one-dimensional estimates. Accordingly, `gaussian_multivariate_kl_divergence`
and `multivariate_information_gain` compute the full
$D_{KL}(p(m|d)\,||\,p(m))$ using covariance determinants, precision-weighted
mean shifts, and the full posterior covariance. This can differ from the sum
of marginal KL values when parameters are correlated. When supplied,
`prior_covariance_j_mol2` is used as the full $C_M$ in this calculation rather
than reducing the prior to independent marginal variances.

`posterior_correlation_matrix` exposes the corresponding joint trade-offs.
The SVD weak directions provide the local resolution/null-space complement:
a small posterior variance caused only by a strong prior must not be reported
as data resolution. The framework therefore reports both prior-to-posterior
information gain and structural sensitivity directions.

### 8.5 Tarantola-derived methodological audit

The supplied Tarantola reference emphasizes five requirements that are now
explicit in this implementation:

1. define a forward operator and distinguish model space from data space;
2. represent experimental, prior, and theoretical/model uncertainties
  separately, preferably with covariance operators;
3. formulate the posterior as the conjunction of prior information and data
  information rather than as an optimizer-only result;
4. analyze resolution and null spaces in the joint model space; and
5. treat nonlinear uncertainty with posterior sampling or explicit
  transformation/Monte Carlo checks rather than relying only on a linearized
  covariance.

`thermoinvert` satisfies these through Reaktoro forward evaluation,
`objective_components`, correlated Holland--Powell covariance perturbations,
Laplace/grid/MCMC diagnostics, full-SVD weak directions, multivariate KL,
posterior-predictive validation, and explicitly supplied $C_D$/$C_T$
covariances. Remaining limitations are automatic model-error covariance
estimation and an explicit Tarantola-style resolution matrix for nonlinear
bounded problems; the current posterior covariance and SVD directions are the
appropriate local diagnostics, but should not be described as a complete
global resolution operator.

For one global model-error amplitude, `profile_theoretical_covariance_scale`
performs nested bounded MAP refits with a Gaussian scale prior. This provides a
practical estimated-$C_T$ amplitude while remaining distinct from a full
hierarchical posterior over experiment-, study-, or laboratory-specific
discrepancy.

Fixed metadata-group scales $s_{T,g}$ are supported when group assignments and
scales are supplied. Their values are not estimated jointly in the current
posterior; the remaining hierarchical gap is joint inference of group scales,
laboratory biases, and experiment-level discrepancy.

For differentiable continuous observables, `local_resolution_matrix` provides
the explicit local linear-Gaussian resolution matrix
$R=(G^T C_D^{-1}G+C_M^{-1})^{-1}G^T C_D^{-1}G$ and its diagonal. It must be
reported as a local approximation at the sensitivity point, not applied across
phase-boundary discontinuities or interpreted as global uniqueness.

Cold-start MAP records are additionally clustered in normalized parameter
bound space by `diagnose_map_globality`. Multiple separated basins within an
objective tolerance are reported as competing near-optimal solutions. One
near-optimal basin is supporting evidence that the attempted starts agree; it
does not establish global uniqueness or replace posterior exploration.

## 9. Validation and Reproducibility

The framework supports:

- leave-one-out cross-validation;
- composition/fraction posterior-predictive checks;
- P-T composition/fraction posterior-predictive checks;
- bracket posterior-predictive satisfaction fractions and misfit distributions;
- synthetic analytic recovery tests;
- real Reaktoro equilibrium smoke tests;
- JSON dataset serialization;
- CSV/JSON result archival;
- input-perturbation sample archival, including correlated covariance key order
  and PSD/rank/conditioning provenance;
- PDF/PNG diagnostic figures.

Every stochastic operation accepts a seed. Input datasets, inversion results,
MCMC samples, LOOCV results, and grid posteriors can be archived and reloaded.
For multi-chain inference, the archive retains every chain as well as R-hat,
ESS, starting points, and seed schedule, so convergence diagnostics remain
auditable rather than being reduced to aggregate summaries.

## 10. Scientific Guardrails

The implementation does not make a parameter scientifically identifiable merely
because an optimizer converged. Before interpreting a result:

1. inspect SVD rank, condition number, and weak directions;
2. compare posterior and prior widths;
3. report parameter correlations and KL information gain;
4. inspect multiple MCMC chains, $\hat{R}$, ESS, and trace plots;
5. inspect posterior-predictive coverage;
6. compare alternative databases and objective forms;
7. distinguish analytical uncertainty, thermodynamic uncertainty, inverse
   non-uniqueness, and structural model uncertainty;
8. avoid claiming database validation when the forward-model database or
   solution models differ from the reference method.

The framework supports a formal correlated Gaussian prior through
`InversionProblem(prior_covariance_j_mol2=...)`. The matrix must be positive
definite and ordered consistently with `parameters`; an external Holland--
Powell matrix may need subsetting or regularization before it can be used as
the fitted prior covariance.

## 11. Algorithm Pseudocode

```text
read experiment dataset and thermodynamic database
select problem type:
    database-parameter inversion, mineral-composition P-T inversion,
    or bracket-based P-T inversion
validate parameter bounds, priors, phase names, and observation uncertainties

for each proposed parameter vector theta:
    if database-parameter inversion:
        create isolated thermodynamic overlay
        solve equilibrium at reported P-T conditions
    if composition-based P-T inversion:
        replace experiment P-T with theta
        solve equilibrium
    if bracket-based P-T inversion:
        translate each bracket around theta
        solve both endpoint equilibria
    calculate phase, fraction, and composition likelihood terms
    calculate bracket inequality terms when applicable
    calculate Gaussian prior penalty
    return total objective

run bounded multi-start MAP search
run sensitivity/SVD diagnostics
run Laplace, grid, or MCMC posterior diagnostics
run posterior-predictive and cross-validation checks
calculate shrinkage, correlations, weak directions, and information gain
archive machine-readable results and create diagnostic figures
report estimates with uncertainty, identifiability, and model limitations
```

Continuous posterior-predictive checks distinguish parameter uncertainty from
full observation prediction. The default propagates posterior parameter samples
only; the opt-in observation-noise mode adds declared $C_D+C_T$ uncertainty.

## 12. Current Limitations

- Solid-solution/Margules parameter inversion requires a forward model that
  constructs and evaluates Reaktoro solid-solution phases.
- Input perturbation accepts declared analytical sigmas and either independent
  thermodynamic nuisance sigmas or an externally supplied correlated covariance
  matrix. It does not infer database-error covariance automatically.
- Bracket PPC reports phase-inequality satisfaction and misfit, not continuous
  residual intervals.
- The MCMC sampler is a transparent baseline, not a replacement for tuned HMC
  or ensemble samplers on high-dimensional problems.
- Real scientific calibration still requires curated literature-referenced
  experimental or natural datasets.
