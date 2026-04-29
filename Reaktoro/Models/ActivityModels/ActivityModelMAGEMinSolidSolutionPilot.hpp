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

    /// Optional gradient callback over internal ternary composition coordinates.
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

    /// Quadratic penalty that couples the minimized internal composition to the visible ternary composition.
    real externalCompositionPenalty = 25.0;

    /// Tolerance of the constrained ternary internal minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the constrained ternary internal minimizer.
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

} // namespace Reaktoro