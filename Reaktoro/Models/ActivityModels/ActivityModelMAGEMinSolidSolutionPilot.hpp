// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/PhaseList.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelGlobalizedSolidSolution.hpp>

namespace Reaktoro {

struct MAGEMinImportedConstrainedTernarySolutionOptions;

/// Local reduced constrained ternary model passed to optional pilot-local minimizers.
struct MAGEMinConstrainedTernaryLocalModel
{
    /// Stable MAGEMin model identifier associated with this local model instance.
    String modelId;

    /// Temperature used for this local model evaluation.
    real T = 0.0;

    /// Externally visible composition used to define the internal objective.
    ArrayXr visiblex;

    /// Objective callback over internal ternary composition coordinates.
    Fn<real(ArrayXrConstRef)> objective;

    /// Optional gradient callback over internal ternary composition coordinates (∇f, n-vector).
    Fn<ArrayXr(ArrayXrConstRef)> gradient;

    /// Lower bounds applied by bounded minimizers.
    ArrayXr lowerBounds;

    /// Upper bounds applied by bounded minimizers.
    ArrayXr upperBounds;

    /// Whether the local minimizer should enforce the unity simplex constraint.
    bool enforceUnityConstraint = true;

    /// Tolerance target for local minimization.
    real tolerance = 1.0e-10;

    /// Maximum number of local minimizer iterations.
    Index maxIterations = 256;

    ///@name Nonlinear Constraint Hooks (optional for advanced minimizers)
    ///@{

    /// Optional callback to evaluate nonlinear inequality constraints c(x) <= 0.
    /// Returns m-vector of constraint values at x. Absence indicates no nonlinear constraints.
    Fn<ArrayXr(ArrayXrConstRef)> constraints;

    /// Optional callback to evaluate the Jacobian of nonlinear constraints (m × n matrix).
    /// Required if constraints callback is provided and used by constraint-aware minimizers.
    /// Returns m × n matrix where row i = ∇c_i(x).
    Fn<MatrixXr(ArrayXrConstRef)> constraintJacobian;

    /// Lower bounds for nonlinear constraints c_i(x) >= constraintLowerBounds[i].
    ArrayXr constraintLowerBounds;

    /// Upper bounds for nonlinear constraints c_i(x) <= constraintUpperBounds[i].
    ArrayXr constraintUpperBounds;

    /// Quadratic penalty weight applied to nonlinear constraint violations in merit-based searches.
    real constraintPenaltyWeight = 1.0e3;

    /// Optional logarithmic barrier weight applied to feasible nonlinear constraints.
    /// A zero value disables the barrier term.
    real constraintBarrierWeight = 0.0;

    /// Whether line-search style minimizers should reject nonlinear-infeasible trial steps outright.
    bool requireFeasibleTrialPoints = false;

    ///@}

    ///@name Second-Order Information (optional for quasi-Newton and Newton-type methods)
    ///@{

    /// Optional callback to evaluate the Hessian of the objective function (n × n symmetric matrix).
    /// Arguments: (x, multipliers) where multipliers are m-vector of Lagrange multipliers for constraints.
    /// Returns ∇²f(x) + ∑_i lambda_i ∇²c_i(x) (Hessian of Lagrangian).
    Fn<MatrixXr(ArrayXrConstRef, ArrayXrConstRef)> objectiveHessian;

    /// Flag indicating whether second-order information (Hessian) should be used by the minimizer.
    /// Minimizers may ignore this if Hessian callback is not provided.
    bool useSecondOrderInfo = false;

    /// Trust-region radius applied to raw search displacements before projection.
    /// A non-positive value disables trust-region clipping.
    real trustRegionRadius = 0.0;

    ///@}
};

/// Function type for overriding the local constrained ternary minimizer used by imported MAGEMin pilots.
using MAGEMinConstrainedTernaryMinimizer = Fn<GlobalizedSolidSolutionInternalResult(
    MAGEMinImportedConstrainedTernarySolutionOptions const&,
    real,
    ArrayXrConstRef,
    Optional<ArrayXr>)>;

/// Function type for overriding constrained ternary minimization using a local reduced model contract.
using MAGEMinConstrainedTernaryLocalModelMinimizer = Fn<GlobalizedSolidSolutionInternalResult(
    MAGEMinConstrainedTernaryLocalModel const&,
    Optional<ArrayXr>)>;

/// Function type compatible with MAGEMin TC/NLopt multi-constraint callbacks.
///
/// Signature mirrors NLopt's `nlopt_add_inequality_mconstraint` callback shape:
/// `(m, result, n, x, grad, data)` where `grad` is an optional flattened row-major m×n Jacobian.
using MAGEMinTCMConstraintCallback = Fn<void(
    unsigned,
    double*,
    unsigned,
    const double*,
    double*,
    void*)>;

/// Bridge data for mapping MAGEMin TC flattened mconstraint callbacks to Reaktoro's dense constraint contract.
struct MAGEMinTCMConstraintBridge
{
    /// Number of constraints (`m`).
    unsigned numConstraints = 0;

    /// Number of local variables (`n`).
    unsigned numVariables = 0;

    /// Lower bounds for each mapped constraint value.
    ArrayXr constraintLowerBounds;

    /// Upper bounds for each mapped constraint value.
    ArrayXr constraintUpperBounds;

    /// MAGEMin/TC style mconstraint callback.
    MAGEMinTCMConstraintCallback callback;

    /// Optional map from native TC variables to visible pilot composition coordinates.
    ///
    /// When provided, the local-model objective/gradient are evaluated in visible composition
    /// space while minimization and nonlinear constraints remain in native TC variable space.
    Fn<ArrayXr(ArrayXrConstRef)> nativeToVisible;

    /// Optional Jacobian of `nativeToVisible` (visible_size x numVariables).
    ///
    /// When provided together with a visible-space gradient callback, the adapter uses
    /// chain-rule projection for native-space gradients. Otherwise finite-difference
    /// gradients are used in native space.
    Fn<MatrixXr(ArrayXrConstRef)> nativeToVisibleJacobian;

    /// Optional map from visible pilot composition coordinates to native TC variables.
    ///
    /// Used to convert visible warmstarts to native-variable warmstarts when dimensions differ.
    Fn<ArrayXr(ArrayXrConstRef)> visibleToNative;

    /// Optional lower bounds for native TC variables.
    ///
    /// If populated, size must equal `numVariables` and overrides local-model lower bounds.
    ArrayXr variableLowerBounds;

    /// Optional upper bounds for native TC variables.
    ///
    /// If populated, size must equal `numVariables` and overrides local-model upper bounds.
    ArrayXr variableUpperBounds;

    /// Optional override for enforcing the unity constraint in native TC variable space.
    ///
    /// For native TC coordinates, this is commonly disabled because constraints are already
    /// imposed through `callback`.
    Optional<bool> enforceUnityConstraint;

    /// Optional opaque pointer forwarded to `callback`.
    void* userData = nullptr;
};

/// Function type for attaching family-specific diagnostics payloads to local-model minimization outcomes.
using MAGEMinConstrainedTernaryLocalModelDiagnostics = Fn<Map<String, Any>(
    MAGEMinConstrainedTernaryLocalModel const&,
    GlobalizedSolidSolutionInternalResult const&)>;

/// Thin pilot wrapper for future MAGEMin-backed solid-solution models.
struct MAGEMinSolidSolutionPilotOptions
{
    /// Admissible branch set exposed to the outer problem.
    Vec<GlobalizedSolidSolutionBranch> branches;

    /// Separator used when duplicating branch-specific phases.
    String suffixSeparator = "#";
};

/// Reaktoro-side branch policy shared by imported MAGEMin pilot models.
struct MAGEMinSolidSolutionPilotBranchPolicyOptions
{
    /// Admissible branch set used for diagnostics and optional split requests.
    Vec<GlobalizedSolidSolutionBranch> branches;

    /// Tolerance used when screening externally visible compositions against `branches`.
    real branchTolerance = 1.0e-8;

    /// Half-width of the score hysteresis used to prefer the persisted branch when branch scores are similar.
    real branchScoreHysteresis = 1.0e-8;

    /// Optional candidate generator used to screen branches before branch selection.
    GlobalizedSolidSolutionCandidateGenerator candidateGenerator;

    /// Optional branch-local stability screen evaluated after candidate screening and before branch selection.
    GlobalizedSolidSolutionCandidateStabilityCriterion stabilityCriterion;
};

/// Imported thermodynamic parameters for a small MAGEMin binary phase model.
struct MAGEMinImportedBinarySolutionThermoModel
{
    /// Stable MAGEMin model identifier.
    String modelId;

    /// Stable MAGEMin identifier for endmember 0.
    String endmember0;

    /// Stable MAGEMin identifier for endmember 1.
    String endmember1;

    /// Binary interaction parameter imported from MAGEMin (in J/mol).
    real W = 0.0;

    /// Site multiplicity used by simple `m*(x0*ln x0 + x1*ln x1)` ideal mixing terms.
    real idealSiteMultiplicity = 1.0;
};

/// Combined options for imported MAGEMin binary pilot models.
struct MAGEMinImportedBinarySolutionOptions
{
    /// Imported thermodynamic parameters.
    MAGEMinImportedBinarySolutionThermoModel thermo;

    /// Reaktoro-side branch policy and split screening.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
};

/// Imported thermodynamic callbacks for a constrained ternary MAGEMin phase model.
struct MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    /// Stable MAGEMin model identifier.
    String modelId;

    /// Stable MAGEMin identifiers for the imported endmembers.
    Strings endmembers;

    /// Pairwise excess interaction between endmembers 0 and 1 (in J/mol).
    real W01 = 0.0;

    /// Pairwise excess interaction between endmembers 0 and 2 (in J/mol).
    real W02 = 0.0;

    /// Pairwise excess interaction between endmembers 1 and 2 (in J/mol).
    real W12 = 0.0;

    /// Optional callback returning excess chemical potentials for the current ternary composition.
    Fn<ArrayXr(ArrayXrConstRef)> excessChemicalPotentials;

    /// Return the ideal configurational contribution to Gibbs energy (in J/mol).
    Fn<real(real, ArrayXrConstRef)> idealGibbs;

    /// Return the ideal logarithmic activity terms for each endmember.
    Fn<ArrayXr(ArrayXrConstRef)> idealLnActivities;
};

/// Structured candidate-proposal options for imported constrained ternary models.
struct MAGEMinStructuredTernaryProposalOptions
{
    /// Whether to always include the externally visible composition as a proposal seed.
    bool includeVisibleCompositionSeed = true;

    /// Whether to add dominant-endmember proposal seeds.
    bool includeDominantEndmemberSeeds = true;

    /// Whether to add binary-edge midpoint proposal seeds.
    bool includeBinaryEdgeMidpointSeeds = true;

    /// Dominant-endmember proposal order. When empty, the natural endmember order is used.
    Indices dominantEndmemberOrder;

    /// Priority bias applied to the visible-composition seed.
    real visibleCompositionPriority = -1.0e-6;

    /// Base priority bias applied to dominant-endmember seeds.
    real dominantEndmemberPriority = -2.0e-6;

    /// Increment applied between successive dominant-endmember seeds.
    real dominantEndmemberPriorityStep = -1.0e-6;

    /// Priority bias applied to binary-edge midpoint seeds.
    real binaryEdgePriority = -5.0e-7;
};

/// Combined options for imported constrained ternary MAGEMin pilot models.
struct MAGEMinImportedConstrainedTernarySolutionOptions
{
    /// Imported thermodynamic parameters and callbacks.
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;

    /// Reaktoro-side branch policy and split screening.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;

    /// Structured candidate-proposal policy used before ternary branch-local minimization.
    MAGEMinStructuredTernaryProposalOptions proposals;

    /// Optional local minimizer overriding the default bounded coordinate search for ternary pilot models.
    MAGEMinConstrainedTernaryMinimizer minimizer;

    /// Optional local-model minimizer receiving a richer objective/gradient contract.
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;

    /// Optional NLopt-backed local-model adapter minimizer.
    ///
    /// This callback is intended as a bridge to MAGEMin's native NLopt local solvers
    /// (`SB_NLopt_opt_function.c`, `TC_database/NLopt_opt_function.c`) while keeping
    /// Reaktoro's outer constrained equilibrium formulation unchanged.
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;

    /// Prefer `nloptLocalModelMinimizer` over `localModelMinimizer` when both are provided.
    bool preferNLoptLocalModelMinimizer = false;

    /// Optional callback producing additional diagnostics payload for local-model minimization outcomes.
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;

    /// Built-in minimizer strategy used when `minimizer` is not supplied.
    String defaultMinimizerStrategy = "legacy";

    /// Whether the built-in projected-gradient rollout should compare itself against the legacy bounded search.
    bool compareProjectedGradientAgainstLegacy = false;

    /// Whether the built-in projected-gradient rollout should fall back to the legacy bounded search when it disagrees.
    bool fallbackToLegacyOnProjectedGradientDisagreement = false;

    /// Quadratic penalty that couples the internally minimized ternary composition to the visible composition.
    real externalCompositionPenalty = 25.0;

    /// Tolerance of the constrained ternary internal minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the constrained ternary internal minimizer.
    Index minimizerMaxIterations = 256;

    /// Whether to apply a tangent-plane distance stability check after all branch-local minimizations.
    ///
    /// When enabled, the pilot computes the tangent-plane distance (TPD) from the current
    /// branch-local minimum toward every other branch's local minimum. A negative TPD indicates
    /// that the single-phase state is thermodynamically unstable, and a split request is emitted.
    bool enableTangentPlaneStabilityCheck = false;

    /// Tolerance for the tangent-plane distance criterion (dimensionless, relative to R*T).
    ///
    /// A split request is emitted when the TPD at any competing branch minimum is more negative
    /// than `-tpdTolerance * R * T`. Reducing this value makes the criterion more sensitive.
    real tpdTolerance = 1.0e-4;
};

/// Options for the imported MAGEMin SB11 olivine (`fo`-`fa`) pilot model.
struct MAGEMinSB11OlivineOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
};

/// Options for the imported MAGEMin SB11 wadsleyite (`mgwa`-`fewa`) pilot model.
struct MAGEMinSB11WadsleyiteOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
};

/// Options for the imported MAGEMin SB11 akimotoite (`mgak`-`feak`-`co`) pilot model.
struct MAGEMinSB11AkimotoiteOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;

    /// Optional local minimizer overriding the default bounded coordinate search for this ternary pilot.
    MAGEMinConstrainedTernaryMinimizer minimizer;

    /// Optional local-model minimizer receiving a richer objective/gradient contract.
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;

    /// Optional NLopt-backed local-model adapter minimizer.
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;

    /// Prefer `nloptLocalModelMinimizer` over `localModelMinimizer` when both are provided.
    bool preferNLoptLocalModelMinimizer = false;

    /// Optional callback producing additional diagnostics payload for local-model minimization outcomes.
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;

    /// Quadratic penalty that couples the minimized internal composition to the visible ternary composition.
    real externalCompositionPenalty = 25.0;

    /// Tolerance of the constrained ternary internal minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the constrained ternary internal minimizer.
    Index minimizerMaxIterations = 256;
};

/// Options for the imported MAGEMin SB11 perovskite (`mgpv`-`fepv`-`alpv`) pilot model.
struct MAGEMinSB11PerovskiteOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;

    /// Optional local minimizer overriding the default bounded coordinate search for this ternary pilot.
    MAGEMinConstrainedTernaryMinimizer minimizer;

    /// Optional local-model minimizer receiving a richer objective/gradient contract.
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;

    /// Optional NLopt-backed local-model adapter minimizer.
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;

    /// Prefer `nloptLocalModelMinimizer` over `localModelMinimizer` when both are provided.
    bool preferNLoptLocalModelMinimizer = false;

    /// Optional callback producing additional diagnostics payload for local-model minimization outcomes.
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;

    /// Quadratic penalty that couples the minimized internal composition to the visible ternary composition.
    real externalCompositionPenalty = 25.0;

    /// Tolerance of the constrained ternary internal minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the constrained ternary internal minimizer.
    Index minimizerMaxIterations = 256;
};

/// Options for the imported MAGEMin SB11 calcioferrite (`mgcf`-`fecf`-`nacf`) pilot model.
struct MAGEMinSB11CalcioferriteOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;

    /// Optional local minimizer overriding the default bounded coordinate search for this ternary pilot.
    MAGEMinConstrainedTernaryMinimizer minimizer;

    /// Optional local-model minimizer receiving a richer objective/gradient contract.
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;

    /// Optional NLopt-backed local-model adapter minimizer.
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;

    /// Prefer `nloptLocalModelMinimizer` over `localModelMinimizer` when both are provided.
    bool preferNLoptLocalModelMinimizer = false;

    /// Optional callback producing additional diagnostics payload for local-model minimization outcomes.
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;

    /// Quadratic penalty that couples the minimized internal composition to the visible ternary composition.
    real externalCompositionPenalty = 25.0;

    /// Tolerance of the constrained ternary internal minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the constrained ternary internal minimizer.
    Index minimizerMaxIterations = 256;
};

/// Options for the imported MAGEMin SB21 spinel (`sp`-`hc`) pilot model.
struct MAGEMinSB21SpinelOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
};

/// Options for the imported MAGEMin SB21 NAL (`mnal`-`fnal`-`nnal`) pilot model.
struct MAGEMinSB21NALOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;

    /// Optional local minimizer overriding the default bounded coordinate search for this ternary pilot.
    MAGEMinConstrainedTernaryMinimizer minimizer;

    /// Optional local-model minimizer receiving a richer objective/gradient contract.
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;

    /// Optional NLopt-backed local-model adapter minimizer.
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;

    /// Prefer `nloptLocalModelMinimizer` over `localModelMinimizer` when both are provided.
    bool preferNLoptLocalModelMinimizer = false;

    /// Optional callback producing additional diagnostics payload for local-model minimization outcomes.
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;

    /// Quadratic penalty that couples the minimized internal composition to the visible ternary composition.
    real externalCompositionPenalty = 25.0;

    /// Tolerance of the constrained ternary internal minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the constrained ternary internal minimizer.
    Index minimizerMaxIterations = 256;
};

/// Options for the imported MAGEMin SB21 calcioferrite (`mgcf`-`fecf`-`nacf`) pilot model.
struct MAGEMinSB21CalcioferriteOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;

    /// Optional local minimizer overriding the default bounded coordinate search for this ternary pilot.
    MAGEMinConstrainedTernaryMinimizer minimizer;

    /// Optional local-model minimizer receiving a richer objective/gradient contract.
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;

    /// Optional NLopt-backed local-model adapter minimizer.
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;

    /// Prefer `nloptLocalModelMinimizer` over `localModelMinimizer` when both are provided.
    bool preferNLoptLocalModelMinimizer = false;

    /// Optional callback producing additional diagnostics payload for local-model minimization outcomes.
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;

    /// Quadratic penalty that couples the minimized internal composition to the visible ternary composition.
    real externalCompositionPenalty = 25.0;

    /// Tolerance of the constrained ternary internal minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the constrained ternary internal minimizer.
    Index minimizerMaxIterations = 256;

    /// Whether to apply a tangent-plane distance stability check after all branch-local minimizations.
    bool enableTangentPlaneStabilityCheck = false;

    /// Tolerance for the tangent-plane distance criterion (dimensionless, relative to R*T).
    real tpdTolerance = 1.0e-4;
};

/// Options for the imported MAGEMin SB21 OPX (`en`-`fs`-`mgts`-`odi`) pilot model.
struct MAGEMinSB21OPXOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;

    /// Optional local minimizer overriding the default projected-gradient minimizer for this pilot.
    MAGEMinConstrainedTernaryMinimizer minimizer;

    /// Optional local-model minimizer receiving a richer objective/gradient contract.
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;

    /// Optional NLopt-backed local-model adapter minimizer.
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;

    /// Prefer `nloptLocalModelMinimizer` over `localModelMinimizer` when both are provided.
    bool preferNLoptLocalModelMinimizer = false;

    /// Optional callback producing additional diagnostics payload for local-model minimization outcomes.
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;

    /// Optional TC-style flattened mconstraint bridge used to build `nloptLocalModelMinimizer`.
    Optional<MAGEMinTCMConstraintBridge> tcMConstraintBridge;

    /// Quadratic penalty that couples the minimized internal composition to the visible composition.
    real externalCompositionPenalty = 25.0;

    /// Tolerance of the internal minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the internal minimizer.
    Index minimizerMaxIterations = 256;
};

/// Options for the imported MAGEMin SB21 CPX (`di`-`he`-`cen`-`cats`-`jd`) pilot model.
struct MAGEMinSB21CPXOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;

    /// Optional local minimizer overriding the default projected-gradient minimizer for this pilot.
    MAGEMinConstrainedTernaryMinimizer minimizer;

    /// Optional local-model minimizer receiving a richer objective/gradient contract.
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;

    /// Optional NLopt-backed local-model adapter minimizer.
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;

    /// Prefer `nloptLocalModelMinimizer` over `localModelMinimizer` when both are provided.
    bool preferNLoptLocalModelMinimizer = false;

    /// Optional callback producing additional diagnostics payload for local-model minimization outcomes.
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;

    /// Optional TC-style flattened mconstraint bridge used to build `nloptLocalModelMinimizer`.
    Optional<MAGEMinTCMConstraintBridge> tcMConstraintBridge;

    /// Quadratic penalty that couples the minimized internal composition to the visible composition.
    real externalCompositionPenalty = 25.0;

    /// Tolerance of the internal minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the internal minimizer.
    Index minimizerMaxIterations = 256;
};

/// Options for the imported MAGEMin SB21 garnet-majorite (`py`-`alm`-`gr`-`mgmj`-`jdmj`) pilot model.
struct MAGEMinSB21GTMJOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;

    /// Optional local minimizer overriding the default projected-gradient minimizer for this pilot.
    MAGEMinConstrainedTernaryMinimizer minimizer;

    /// Optional local-model minimizer receiving a richer objective/gradient contract.
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;

    /// Optional NLopt-backed local-model adapter minimizer.
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;

    /// Prefer `nloptLocalModelMinimizer` over `localModelMinimizer` when both are provided.
    bool preferNLoptLocalModelMinimizer = false;

    /// Optional callback producing additional diagnostics payload for local-model minimization outcomes.
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;

    /// Optional TC-style flattened mconstraint bridge used to build `nloptLocalModelMinimizer`.
    Optional<MAGEMinTCMConstraintBridge> tcMConstraintBridge;

    /// Quadratic penalty that couples the minimized internal composition to the visible composition.
    real externalCompositionPenalty = 25.0;

    /// Tolerance of the internal minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the internal minimizer.
    Index minimizerMaxIterations = 256;
};

/// Options for the imported MAGEMin SB21 plagioclase (`an`-`ab`) binary pilot model.
struct MAGEMinSB21PLGOptions
{
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
};

/// Options for the imported MAGEMin SB21 olivine (`fo`-`fa`) binary pilot model.
struct MAGEMinSB21OlivineOptions
{
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
};

/// Options for the imported MAGEMin SB21 wadsleyite (`mgwa`-`fewa`) binary pilot model.
struct MAGEMinSB21WadsleyiteOptions
{
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
};

/// Options for the imported MAGEMin SB21 ringwoodite (`mgri`-`feri`) binary pilot model.
struct MAGEMinSB21RingwooditeOptions
{
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
};

/// Options for the imported MAGEMin SB21 high-P clinopyroxene (`hpcen`-`hpcfs`) binary pilot model.
struct MAGEMinSB21HPCPXOptions
{
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
};

/// Options for the imported MAGEMin SB21 akimotoite (`mgak`-`feak`-`co`) pilot model.
struct MAGEMinSB21AkimotoiteOptions
{
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
    MAGEMinConstrainedTernaryMinimizer minimizer;
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;
    bool preferNLoptLocalModelMinimizer = false;
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;
    Optional<MAGEMinTCMConstraintBridge> tcMConstraintBridge;
    real externalCompositionPenalty = 25.0;
    real minimizerTolerance = 1.0e-10;
    Index minimizerMaxIterations = 256;
};

/// Options for the imported MAGEMin SB21 perovskite (`mgpv`-`fepv`-`alpv`) pilot model.
struct MAGEMinSB21PerovskiteOptions
{
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
    MAGEMinConstrainedTernaryMinimizer minimizer;
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;
    bool preferNLoptLocalModelMinimizer = false;
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;
    Optional<MAGEMinTCMConstraintBridge> tcMConstraintBridge;
    real externalCompositionPenalty = 25.0;
    real minimizerTolerance = 1.0e-10;
    Index minimizerMaxIterations = 256;
};

/// Options for the imported MAGEMin SB21 post-perovskite (`mppv`-`fppv`-`appv`) pilot model.
struct MAGEMinSB21PostPerovskiteOptions
{
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
    MAGEMinConstrainedTernaryMinimizer minimizer;
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;
    bool preferNLoptLocalModelMinimizer = false;
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;
    Optional<MAGEMinTCMConstraintBridge> tcMConstraintBridge;
    real externalCompositionPenalty = 25.0;
    real minimizerTolerance = 1.0e-10;
    Index minimizerMaxIterations = 256;
};

/// Options for the imported MAGEMin SB21 magnesiowustite (`pe`-`wu`-`anao`) pilot model.
struct MAGEMinSB21MagnesiowustitesOptions
{
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;
    MAGEMinConstrainedTernaryMinimizer minimizer;
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;
    bool preferNLoptLocalModelMinimizer = false;
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;
    Optional<MAGEMinTCMConstraintBridge> tcMConstraintBridge;
    real externalCompositionPenalty = 25.0;
    real minimizerTolerance = 1.0e-10;
    Index minimizerMaxIterations = 256;
};

/// Options for the first Holland-Powell ig_opx xeos-native pilot model.
struct MAGEMinHPIGOPXOptions
{
    /// Reaktoro-side branch policy used for diagnostics and optional split requests.
    MAGEMinSolidSolutionPilotBranchPolicyOptions branchPolicy;

    /// Optional local-model minimizer receiving the xeos-native objective/gradient contract.
    MAGEMinConstrainedTernaryLocalModelMinimizer localModelMinimizer;

    /// Optional NLopt-backed local-model minimizer.
    MAGEMinConstrainedTernaryLocalModelMinimizer nloptLocalModelMinimizer;

    /// Prefer `nloptLocalModelMinimizer` over `localModelMinimizer` when both are provided.
    bool preferNLoptLocalModelMinimizer = false;

    /// Optional callback producing additional diagnostics payload for local-model minimization outcomes.
    MAGEMinConstrainedTernaryLocalModelDiagnostics localModelDiagnostics;

    /// Optional TC-style flattened mconstraint bridge used to inject opx_ig_c constraints.
    Optional<MAGEMinTCMConstraintBridge> tcMConstraintBridge;

    /// Optional callback that returns the per-endmember HP ig_opx reference-state Gibbs energies.
    Fn<ArrayXr(real, real)> referenceState;

    /// Quadratic penalty that couples minimized xeos coordinates to visible composition.
    real externalCompositionPenalty = 10.0;

    /// Tolerance of the xeos-native local minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the xeos-native local minimizer.
    Index minimizerMaxIterations = 256;
};

/// Return a reusable imported MAGEMin binary pilot model with simple binary ideal mixing.
auto MAGEMinSolidSolutionPilotModelImportedBinary(
    MAGEMinImportedBinarySolutionOptions options) -> GlobalizedSolidSolutionModel;

/// Return a reusable imported constrained ternary MAGEMin pilot model.
auto MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(
    MAGEMinImportedConstrainedTernarySolutionOptions options) -> GlobalizedSolidSolutionModel;

/// Return the first real MAGEMin-backed pilot model imported into Reaktoro.
///
/// This model reproduces the small binary `sb11_ol` olivine surface from MAGEMin's SB11
/// database using the real `fo`-`fa` interaction parameter and site-mixing multiplicity.
auto MAGEMinSolidSolutionPilotModelSB11Olivine(
    MAGEMinSB11OlivineOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb11_wa` wadsleyite binary pilot model.
auto MAGEMinSolidSolutionPilotModelSB11Wadsleyite(
    MAGEMinSB11WadsleyiteOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb11_ak` akimotoite constrained ternary pilot model.
auto MAGEMinSolidSolutionPilotModelSB11Akimotoite(
    MAGEMinSB11AkimotoiteOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb11_pv` perovskite constrained ternary pilot model.
auto MAGEMinSolidSolutionPilotModelSB11Perovskite(
    MAGEMinSB11PerovskiteOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb11_cf` calcioferrite constrained ternary pilot model.
auto MAGEMinSolidSolutionPilotModelSB11Calcioferrite(
    MAGEMinSB11CalcioferriteOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_sp` spinel binary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21Spinel(
    MAGEMinSB21SpinelOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_nal` constrained ternary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21NAL(
    MAGEMinSB21NALOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_cf` constrained ternary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21Calcioferrite(
    MAGEMinSB21CalcioferriteOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_opx` orthopyroxene 4-endmember pilot model.
auto MAGEMinSolidSolutionPilotModelSB21OPX(
    MAGEMinSB21OPXOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_cpx` clinopyroxene 5-endmember pilot model.
auto MAGEMinSolidSolutionPilotModelSB21CPX(
    MAGEMinSB21CPXOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_gtmj` garnet-majorite 5-endmember pilot model.
auto MAGEMinSolidSolutionPilotModelSB21GTMJ(
    MAGEMinSB21GTMJOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_plg` plagioclase binary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21PLG(
    MAGEMinSB21PLGOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_ol` olivine binary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21Olivine(
    MAGEMinSB21OlivineOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_wa` wadsleyite binary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21Wadsleyite(
    MAGEMinSB21WadsleyiteOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_ri` ringwoodite binary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21Ringwoodite(
    MAGEMinSB21RingwooditeOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_hpcpx` high-P clinopyroxene binary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21HPCPX(
    MAGEMinSB21HPCPXOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_ak` akimotoite ternary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21Akimotoite(
    MAGEMinSB21AkimotoiteOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_pv` perovskite ternary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21Perovskite(
    MAGEMinSB21PerovskiteOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_ppv` post-perovskite ternary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21PostPerovskite(
    MAGEMinSB21PostPerovskiteOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the imported MAGEMin `sb21_mw` magnesiowustite ternary pilot model.
auto MAGEMinSolidSolutionPilotModelSB21Magnesiowustites(
    MAGEMinSB21MagnesiowustitesOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Return the first Holland-Powell ig_opx xeos-native pilot model.
auto MAGEMinSolidSolutionPilotModelHPIGOPX(
    MAGEMinHPIGOPXOptions options = {}) -> GlobalizedSolidSolutionModel;

/// Bind a reduced MAGEMin-style solid-solution model to a single unsplit pilot phase.
auto MAGEMinSolidSolutionPilotPhase(
    Phase const& phase,
    GlobalizedSolidSolutionModel model) -> Phase;

/// Create a reusable definition for a MAGEMin-style pilot phase that can duplicate branches on demand.
auto MAGEMinSolidSolutionPilotDefinition(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    MAGEMinSolidSolutionPilotOptions options = {}) -> GlobalizedSolidSolutionPhaseDefinition;

/// Duplicate a MAGEMin-style pilot phase into one phase per admissible branch.
auto MAGEMinSolidSolutionPilotPhases(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    MAGEMinSolidSolutionPilotOptions options = {}) -> PhaseList;

/// Run the built-in projected-gradient solver on a reduced local constrained ternary model.
///
/// This utility exposes Reaktoro's projected-gradient step algorithm through the
/// `MAGEMinConstrainedTernaryLocalModelMinimizer` contract so external implementations
/// can use it with a caller-supplied gradient callback.  The model's `objective` and
/// `gradient` callbacks are called at every iteration.
///
/// Precondition: `model.gradient` must be populated.  A `std::runtime_error` is thrown
/// if the callback is empty.
auto MAGEMinProjectedGradientLocalModelMinimizer(
    MAGEMinConstrainedTernaryLocalModel const& model,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult;

/// Build a local-model adapter that maps TC-style flattened mconstraints to Reaktoro dense constraints.
///
/// The returned minimizer applies the bridge constraints to the provided local model, then delegates to
/// `fallback` (or to `MAGEMinProjectedGradientLocalModelMinimizer` if `fallback` is empty).
auto MAGEMinTCMConstraintBridgeLocalModelAdapter(
    MAGEMinTCMConstraintBridge bridge,
    MAGEMinConstrainedTernaryLocalModelMinimizer fallback = {}) -> MAGEMinConstrainedTernaryLocalModelMinimizer;

} // namespace Reaktoro