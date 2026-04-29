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

// C++ includes
#include <cstdint>

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Core/ActivityProps.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/PhaseList.hpp>
#include <Reaktoro/Models/ActivityModels/InternallyMinimizedSolidSolution.hpp>

namespace Reaktoro {

class ChemicalSystem;

/// Sentinel used when no branch is explicitly requested or selected.
constexpr auto GlobalizedSolidSolutionNoBranch = static_cast<Index>(-1);

/// Metadata describing one admissible branch of a reduced solid-solution model.
struct GlobalizedSolidSolutionBranch
{
    /// A stable identifier for the branch.
    String id;

    /// A display label for the branch.
    String label;

    /// Lower bounds on the externally visible composition coordinates for this branch.
    ArrayXr lowerBounds;

    /// Upper bounds on the externally visible composition coordinates for this branch.
    ArrayXr upperBounds;
};

/// Split request emitted by a reduced solid-solution model when the outer problem should duplicate branches.
struct GlobalizedSolidSolutionSplitRequest
{
    /// Whether the reduced model requests branch duplication in the outer problem.
    bool requested = false;

    /// The base phase name that should react to the split request, if known.
    String phaseName;

    /// The branch that triggered the split request.
    Index triggeringBranch = GlobalizedSolidSolutionNoBranch;

    /// The branch indices that should be duplicated in the outer problem.
    Indices branches;

    /// Optional branch identifiers matching `branches`.
    Strings branchIds;

    /// Optional reason attached to the split request.
    String reason;

    /// Extra diagnostics attached to the split request.
    Map<String, Any> extra;
};

/// Cache retained between successive reduced-model evaluations.
struct GlobalizedSolidSolutionState
{
    /// The most recent ChemicalProps state identifier seen by the reduced model.
    std::uint64_t chemicalPropsStateId = 0;

    /// The branch cached for the current ChemicalProps state.
    Index cachedBranchForState = GlobalizedSolidSolutionNoBranch;

    /// The minimized internal composition cached for the current ChemicalProps state.
    ArrayXr cachedInternalx;

    /// The branch selected during the previous evaluation.
    Index selectedBranch = GlobalizedSolidSolutionNoBranch;

    /// The number of times the reduced model has been evaluated.
    Index numEvaluations = 0;

    /// The last temperature used to evaluate the reduced model.
    real lastT = 0.0;

    /// The last pressure used to evaluate the reduced model.
    real lastP = 0.0;

    /// The last externally visible composition seen by the reduced model.
    ArrayXr lastx;

    /// The internal composition selected during the last successful branch minimization.
    ArrayXr lastInternalx;

    /// The most recent split request published by the reduced model.
    GlobalizedSolidSolutionSplitRequest lastSplitRequest;

    /// Extra internal state carried by the reduced model implementation.
    Map<String, Any> data;
};

/// Input passed to an internally minimized solid-solution model.
struct GlobalizedSolidSolutionInput
{
    /// The temperature of the phase (in K).
    real const& T;

    /// The pressure of the phase (in Pa).
    real const& P;

    /// The externally visible composition coordinates of the phase.
    ArrayXrConstRef x;

    /// Extra data carried over from the previous evaluation of the phase activity model.
    Map<String, Any> const& extra;

    /// Cache state that can be reused for warm-started internal solves.
    SharedPtr<GlobalizedSolidSolutionState> state;

    /// The branch that the caller wants to force, if any.
    Index requestedBranch = GlobalizedSolidSolutionNoBranch;
};

/// Output returned by an internally minimized solid-solution model.
struct GlobalizedSolidSolutionOutput
{
    /// The corrective molar Gibbs energy of the phase after internal minimization (in J/mol).
    real Gx = 0.0;

    /// The corrective molar enthalpy of the phase after internal minimization (in J/mol).
    real Hx = 0.0;

    /// The corrective molar isobaric heat capacity of the phase after internal minimization (in J/(mol K)).
    real Cpx = 0.0;

    /// The corrective molar volume of the phase after internal minimization (in m^3/mol).
    real Vx = 0.0;

    /// The temperature derivative of Vx at fixed pressure and external composition.
    real VxT = 0.0;

    /// The pressure derivative of Vx at fixed temperature and external composition.
    real VxP = 0.0;

    /// The species activity coefficients in natural log.
    ArrayXr ln_g;

    /// The species activities in natural log.
    ArrayXr ln_a;

    /// The composition derivatives of Vx at fixed temperature and pressure.
    ArrayXr Vxi;

    /// The diagnosed state of matter of the phase.
    StateOfMatter som = StateOfMatter::Solid;

    /// The branch selected during the evaluation.
    Index selectedBranch = GlobalizedSolidSolutionNoBranch;

    /// Metadata for the selected branch.
    GlobalizedSolidSolutionBranch branch;

    /// Metadata for all branches considered by the reduced model.
    Vec<GlobalizedSolidSolutionBranch> branches;

    /// Updated warm-start cache returned by the reduced model.
    SharedPtr<GlobalizedSolidSolutionState> state;

    /// Request to duplicate branches in the outer problem.
    GlobalizedSolidSolutionSplitRequest splitRequest;

    /// Extra diagnostics produced by the internal branch/minimization logic.
    Map<String, Any> extra;
};

/// One candidate internal state proposed before branch-local refinement.
struct GlobalizedSolidSolutionCandidate
{
    /// The branch this candidate belongs to.
    Index branch = GlobalizedSolidSolutionNoBranch;

    /// Optional warm-start for the internal coordinates of this candidate.
    ArrayXr initialInternalx;

    /// Bias applied during candidate screening. Smaller values are preferred.
    real priority = 0.0;

    /// Optional split request attached during candidate screening.
    GlobalizedSolidSolutionSplitRequest splitRequest;

    /// Extra diagnostics attached to the candidate.
    Map<String, Any> extra;
};

/// Result of screening one branch-local candidate for local stability.
struct GlobalizedSolidSolutionCandidateStability
{
    /// Whether the screened candidate is admissible for outer refinement.
    bool stable = true;

    /// Penalty added to the candidate score when ranking admissible candidates.
    real penalty = 0.0;

    /// Optional reason attached when the candidate is rejected or penalized.
    String reason;

    /// Optional split request emitted instead of rejecting the candidate outright.
    GlobalizedSolidSolutionSplitRequest splitRequest;

    /// Extra diagnostics produced by the stability screen.
    Map<String, Any> extra;
};

/// Function type for a reduced solid-solution thermodynamics model.
using GlobalizedSolidSolutionModel = Fn<GlobalizedSolidSolutionOutput(GlobalizedSolidSolutionInput)>;

/// Function type for a candidate generator that proposes branch-local internal states before refinement.
using GlobalizedSolidSolutionCandidateGenerator = Fn<Vec<GlobalizedSolidSolutionCandidate>(GlobalizedSolidSolutionInput const&, Vec<GlobalizedSolidSolutionBranch> const&)>;

/// Function type for screening a branch-local candidate after local refinement and before outer selection.
using GlobalizedSolidSolutionCandidateStabilityCriterion = Fn<GlobalizedSolidSolutionCandidateStability(GlobalizedSolidSolutionInput const&, GlobalizedSolidSolutionBranch const&, ArrayXrConstRef, real)>;

/// Result of composing and ranking one branch-local candidate.
struct GlobalizedSolidSolutionBranchSelection
{
    /// The selected branch index.
    Index branch = GlobalizedSolidSolutionNoBranch;

    /// The raw branch-local score before hysteresis and stability penalties.
    real score = 0.0;

    /// The internal composition retained for the selected branch.
    ArrayXr internalx;

    /// Extra diagnostics attached to the selected branch.
    Map<String, Any> extra;

    /// Whether the branch-local result reused a cached state.
    bool reusedState = false;

    /// Whether the branch-local result reused a warm start.
    bool usedWarmstart = false;

    /// Whether the selected candidate passed the optional stability screen.
    bool stable = true;

    /// Penalty applied by the optional stability screen.
    real stabilityPenalty = 0.0;

    /// Split request attached to the selected branch.
    GlobalizedSolidSolutionSplitRequest splitRequest;
};

/// Function type for evaluating one branch-local candidate.
using GlobalizedSolidSolutionBranchEvaluator = Fn<GlobalizedSolidSolutionBranchSelection(GlobalizedSolidSolutionCandidate const&, GlobalizedSolidSolutionInput const&, Vec<GlobalizedSolidSolutionBranch> const&)>;

/// Function type for generating a default split request after branch-local screening.
using GlobalizedSolidSolutionSplitRequestGenerator = Fn<GlobalizedSolidSolutionSplitRequest(GlobalizedSolidSolutionInput const&, Vec<GlobalizedSolidSolutionBranch> const&, Index)>;

/// Options controlling the shared default branch-candidate seeding policy.
struct GlobalizedSolidSolutionDefaultCandidateOptions
{
    /// Tolerance used when checking whether a composition lies in a branch.
    real branchTolerance = 1.0e-8;

    /// Priority assigned to a reusable cached-state candidate.
    real cachedStatePriority = -1.0;

    /// Priority assigned to the previously selected branch during branch screening.
    real preferredBranchPriority = 0.0;

    /// Whether a cached-state candidate requires cached internal coordinates to be available.
    bool requireCachedStateWarmstart = false;

    /// Diagnostic key used to store the candidate source string.
    String sourceKey;

    /// Diagnostic source string used for reusable cached-state candidates.
    String stateCacheSource = "state-cache";

    /// Diagnostic source string used for forced requested-branch candidates.
    String requestedBranchSource = "requested-branch";

    /// Diagnostic source string used for ordinary branch-screen candidates.
    String branchScreenSource = "branch-screen";

    /// Error raised when a forced branch lies outside the admissible branch set.
    String invalidRequestedBranchMessage = "Requested globalized solid-solution branch is out of range.";
};

/// Options controlling a shared branch-ambiguity stability screen.
struct GlobalizedSolidSolutionBranchAmbiguityStabilityOptions
{
    /// Tolerance used when checking branch membership and nearest-branch ambiguity.
    real branchTolerance = 1.0e-8;

    /// Key used to publish the current branch violation into diagnostics and split requests.
    String violationKey;

    /// Key used to publish the nearest-branch violation into diagnostics and split requests.
    String nearestViolationKey;

    /// Key used to publish the count of equally nearest branches into diagnostics.
    String nearestBranchCountKey;

    /// Minimum score mismatch that still counts as branch ambiguity.
    real ambiguousScoreTolerance = 1.0e-8;

    /// Extra penalty applied to admissible but ambiguous candidates.
    real ambiguityPenalty = 0.0;

    /// Whether ambiguous off-branch candidates should emit split requests instead of only a penalty.
    bool requestSplitOnAmbiguity = true;

    /// Reason used when a split request is emitted for an ambiguous candidate.
    String splitReason = "branch-ambiguity-instability";

    /// Reason used when an off-branch candidate is rejected without a split request.
    String unstableReason = "internal-composition-outside-branch";
};

/// Named seam-owned stability policies for reusable globalized solid-solution call sites.
enum class NamedGlobalizedSolidSolutionStabilityPolicy
{
    BranchAmbiguity,
    MAGEMinPilotBranchAmbiguity,
};

/// Definition of a globalized solid-solution phase family that can be duplicated on demand.
struct GlobalizedSolidSolutionPhaseDefinition
{
    /// Base phase name used to match the unsplit phase and its duplicates.
    String phaseName;

    /// Prototype phase bound to the unsplit reduced model.
    Phase prototype;

    /// Reduced model used by both the prototype and duplicates.
    GlobalizedSolidSolutionModel model;

    /// Branches admissible for the phase family.
    Vec<GlobalizedSolidSolutionBranch> branches;

    /// Separator used when naming duplicated branches.
    String suffixSeparator = "#";
};

/// A generic constrained internal minimization problem for reduced solid-solution models.
struct GlobalizedSolidSolutionInternalProblem
{
    /// Objective to be minimized over the internal coordinates.
    Fn<real(ArrayXrConstRef)> objective;

    /// Initial guess for the internal coordinates.
    ArrayXr initialx;

    /// Lower bounds for the internal coordinates.
    ArrayXr lowerBounds;

    /// Upper bounds for the internal coordinates.
    ArrayXr upperBounds;

    /// Initial step size for the coordinate-transfer search.
    real initialStep = 0.25;

    /// Termination tolerance on the step size.
    real tolerance = 1.0e-10;

    /// Maximum number of iterations.
    Index maxIterations = 256;

    /// Whether the internal coordinates must satisfy a unity-sum constraint.
    bool enforceUnityConstraint = true;
};

/// Result of the generic constrained internal minimization problem.
struct GlobalizedSolidSolutionInternalResult
{
    /// Best internal coordinates found.
    ArrayXr x;

    /// Objective value at the best internal coordinates.
    real objective = 0.0;

    /// Number of search iterations performed.
    Index iterations = 0;

    /// Whether the search terminated by reaching the tolerance.
    bool converged = false;
};

/// Minimize a constrained internal solid-solution objective over bounded coordinates.
auto MinimizeGlobalizedSolidSolutionInternalProblem(
    GlobalizedSolidSolutionInternalProblem const& problem) -> GlobalizedSolidSolutionInternalResult;

/// Return one default unconstrained branch covering the full visible composition simplex.
auto DefaultGlobalizedSolidSolutionBranches(Index numCoords) -> Vec<GlobalizedSolidSolutionBranch>;

/// Normalize branch metadata and fill missing bounds for the given visible composition size.
auto NormalizeGlobalizedSolidSolutionBranches(
    Vec<GlobalizedSolidSolutionBranch> branches,
    Index numCoords,
    String const& errorContext) -> Vec<GlobalizedSolidSolutionBranch>;

/// Return the squared violation of the given composition against branch bounds.
auto GlobalizedSolidSolutionBranchViolation(
    ArrayXrConstRef x,
    GlobalizedSolidSolutionBranch const& branch,
    real tolerance) -> real;

/// Return whether the given composition lies within branch bounds up to tolerance.
auto GlobalizedSolidSolutionBranchContains(
    ArrayXrConstRef x,
    GlobalizedSolidSolutionBranch const& branch,
    real tolerance) -> bool;

/// Merge two split requests while keeping branch ids and diagnostics unique.
auto MergeGlobalizedSolidSolutionSplitRequests(
    GlobalizedSolidSolutionSplitRequest lhs,
    GlobalizedSolidSolutionSplitRequest const& rhs) -> GlobalizedSolidSolutionSplitRequest;

/// Return the ChemicalProps state id stored in activity-model diagnostics, if any.
auto CurrentGlobalizedSolidSolutionChemicalPropsStateId(Map<String, Any> const& extra) -> std::uint64_t;

/// Return the shared default set of branch-local candidates for cache reuse, branch forcing, and branch screening.
auto DefaultGlobalizedSolidSolutionCandidates(
    GlobalizedSolidSolutionInput const& input,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    ArrayXrConstRef cachedStateComposition,
    Optional<ArrayXr> cachedStateWarmstart,
    Optional<ArrayXr> branchWarmstart,
    GlobalizedSolidSolutionDefaultCandidateOptions const& options) -> Vec<GlobalizedSolidSolutionCandidate>;

/// Return the default split request emitted when the external composition lies between branches.
auto DefaultGlobalizedSolidSolutionSplitRequest(
    GlobalizedSolidSolutionInput const& input,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    Index selectedBranch,
    real branchTolerance,
    String violationKey = "") -> GlobalizedSolidSolutionSplitRequest;

/// Return a shared stability screen that detects off-branch and ambiguous nearest-branch internal states.
auto BranchAmbiguityGlobalizedSolidSolutionStabilityCriterion(
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionBranchAmbiguityStabilityOptions const& options = {}) -> GlobalizedSolidSolutionCandidateStabilityCriterion;

/// Return a seam-owned named stability screen backed by shared globalized solid-solution policies.
auto NamedGlobalizedSolidSolutionStabilityCriterion(
    String const& name,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionBranchAmbiguityStabilityOptions const& options = {}) -> GlobalizedSolidSolutionCandidateStabilityCriterion;

/// Return a seam-owned named stability screen backed by a typed shared globalized solid-solution policy.
auto NamedGlobalizedSolidSolutionStabilityCriterion(
    NamedGlobalizedSolidSolutionStabilityPolicy policy,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionBranchAmbiguityStabilityOptions const& options = {}) -> GlobalizedSolidSolutionCandidateStabilityCriterion;

/// Compose branch-local candidates, apply optional stability screening, and select the best branch.
auto ComposeGlobalizedSolidSolutionBranch(
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionInput const& input,
    GlobalizedSolidSolutionCandidateGenerator const& candidateGenerator,
    GlobalizedSolidSolutionCandidateGenerator const& defaultCandidateGenerator,
    GlobalizedSolidSolutionCandidateStabilityCriterion const& stabilityCriterion,
    GlobalizedSolidSolutionBranchEvaluator const& evaluator,
    GlobalizedSolidSolutionSplitRequestGenerator const& splitRequestGenerator,
    String const& emptyCandidatesMessage,
    String const& invalidBranchMessage,
    String const& rejectedCandidatesMessage) -> GlobalizedSolidSolutionBranchSelection;

/// Construct an activity-model generator from a reduced internally minimized solid-solution model.
///
/// This adapter is the seam for future MAGEMin-backed solid-solution implementations. The reduced
/// model is responsible for handling internal coordinates, branch selection, and local refinement,
/// while Reaktoro continues to consume the resulting phase properties through its usual activity-model API.
auto ActivityModelGlobalizedSolidSolution(GlobalizedSolidSolutionModel model) -> ActivityModelGenerator;

/// Construct an activity-model generator from a reduced internally minimized solid-solution model
/// with phase-scoped diagnostics stored under the given phase name.
auto ActivityModelGlobalizedSolidSolution(GlobalizedSolidSolutionModel model, String phaseScope) -> ActivityModelGenerator;

/// Construct an activity-model generator that binds the reduced model to a specific branch.
auto ActivityModelGlobalizedSolidSolution(GlobalizedSolidSolutionModel model, Index requestedBranch) -> ActivityModelGenerator;

/// Construct an activity-model generator that binds the reduced model to a specific branch and phase scope.
auto ActivityModelGlobalizedSolidSolution(GlobalizedSolidSolutionModel model, Index requestedBranch, String phaseScope) -> ActivityModelGenerator;

/// Duplicate a phase into one phase per admissible branch of the reduced model.
///
/// Each duplicate keeps the same species and ideal activity model as the input phase, but binds the
/// globalized solid-solution activity model to a single branch so immiscible solutions can be represented
/// as separate phases in the outer Reaktoro minimization.
auto DuplicateGlobalizedSolidSolutionPhaseBranches(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    String suffixSeparator = "#") -> PhaseList;

/// Create a reusable phase definition that binds phase-scoped diagnostics to the reduced model.
auto MakeGlobalizedSolidSolutionPhaseDefinition(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    String suffixSeparator = "#") -> GlobalizedSolidSolutionPhaseDefinition;

/// Apply split requests published by globalized solid-solution phases and rebuild the phase list.
auto ApplyGlobalizedSolidSolutionSplitRequests(
    PhaseList const& phases,
    Map<String, Any> const& extra,
    Vec<GlobalizedSolidSolutionPhaseDefinition> const& definitions) -> PhaseList;

/// Apply split requests published by globalized solid-solution phases and rebuild the chemical system.
auto ApplyGlobalizedSolidSolutionSplitRequests(
    ChemicalSystem const& system,
    Map<String, Any> const& extra,
    Vec<GlobalizedSolidSolutionPhaseDefinition> const& definitions) -> ChemicalSystem;

/// Assemble multiple candidate phase instances from a single logical solid-solution model.
///
/// Creates one Reaktoro phase per candidate in `candidates`, each bound to the candidate's
/// branch so the outer equilibrium solver can simultaneously explore competing internal states
/// (e.g. different immiscible branches or endmember-dominated regions). This is the minimum
/// extension needed to support immiscibility and exsolution without rewriting the outer solver.
///
/// Each candidate phase is named `phase.name() + suffixSeparator + candidate.label`
/// (falling back to branch label/id or branch index when the candidate label is empty).
///
/// The `branches` vector must be the same set used by the `model`.
auto AssembleGlobalizedSolidSolutionCandidatePhases(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    Vec<SolidSolutionCandidateState> const& candidates,
    String suffixSeparator = "#") -> PhaseList;

/// Assemble candidate phases by evaluating a `SolidSolutionCandidateGenerator` at a reference condition.
///
/// Calls `generator(referenceT, referenceP, referencex)` to obtain the candidate list, then
/// delegates to the vector overload. Useful when the set of interesting candidates depends on
/// the pressure-temperature-composition window of interest.
auto AssembleGlobalizedSolidSolutionCandidatePhases(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    SolidSolutionCandidateGenerator const& generator,
    real referenceT,
    real referenceP,
    ArrayXrConstRef referencex,
    String suffixSeparator = "#") -> PhaseList;

} // namespace Reaktoro