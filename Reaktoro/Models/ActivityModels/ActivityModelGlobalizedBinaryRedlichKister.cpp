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

#include "ActivityModelGlobalizedBinaryRedlichKister.hpp"

// C++ includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIdealSolution.hpp>

namespace Reaktoro {

namespace {

using std::log;
using std::pow;

struct BranchCandidate
{
    Index branch = GlobalizedSolidSolutionNoBranch;
    real score = 0.0;
    ArrayXr internalx;
    Map<String, Any> extra;
    bool reusedState = false;
    bool usedWarmstart = false;
    bool stable = true;
    real stabilityPenalty = 0.0;
    GlobalizedSolidSolutionSplitRequest splitRequest;
};

constexpr auto CompositionFloor = 1.0e-12;

auto normalizedBranches(Vec<GlobalizedSolidSolutionBranch> branches) -> Vec<GlobalizedSolidSolutionBranch>
{
    return NormalizeGlobalizedSolidSolutionBranches(
        std::move(branches),
        2,
        "Globalized binary Redlich-Kister branches must define exactly two bounds.");
}

auto branchPenalty(real violation, real scale) -> real
{
    return scale * violation;
}

auto redlichKisterExcessGibbs(GlobalizedBinaryRedlichKisterOptions const& options, real T, real y1) -> real
{
    const auto y2 = 1.0 - y1;
    const auto RT = universalGasConstant * T;
    return (y1*y2*(options.a0 + options.a1*(y1 - y2) + options.a2*pow((y1 - y2), 2))) * RT;
}

auto redlichKisterExcessGibbs(GlobalizedBinaryRedlichKisterOptions const& options, real T, ArrayXrConstRef y) -> real
{
    return redlichKisterExcessGibbs(options, T, y[0]);
}

auto branchObjective(
    GlobalizedBinaryRedlichKisterOptions const& options,
    GlobalizedSolidSolutionBranch const& branch,
    real T,
    ArrayXrConstRef x,
    ArrayXrConstRef y) -> real
{
    const auto RT = universalGasConstant * T;
    const auto violation = GlobalizedSolidSolutionBranchViolation(y, branch, options.branchTolerance);
    const auto penalty = branchPenalty(violation, options.inactiveBranchPenalty) * RT;
    const auto mismatch = options.externalCompositionPenalty * RT * (y - x).matrix().squaredNorm();
    return redlichKisterExcessGibbs(options, T, y) + mismatch + penalty;
}

auto minimizeBranchObjective(
    GlobalizedBinaryRedlichKisterOptions const& options,
    GlobalizedSolidSolutionBranch const& branch,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> Pair<ArrayXr, real>
{
    GlobalizedSolidSolutionInternalProblem problem;
    problem.objective = [=](ArrayXrConstRef y)
    {
        return branchObjective(options, branch, T, x, y);
    };
    problem.initialx = warmstart.value_or(x);
    problem.lowerBounds = branch.lowerBounds;
    problem.upperBounds = branch.upperBounds;
    problem.tolerance = options.minimizerTolerance;
    problem.maxIterations = options.minimizerMaxIterations;
    problem.enforceUnityConstraint = true;

    const auto result = MinimizeGlobalizedSolidSolutionInternalProblem(problem);
    return {result.x, result.objective};
}

auto defaultCandidates(
    GlobalizedBinaryRedlichKisterOptions const& options,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionInput const& input) -> Vec<GlobalizedSolidSolutionCandidate>
{
    Optional<ArrayXr> cachedWarmstart = std::nullopt;
    Optional<ArrayXr> branchWarmstart = std::nullopt;
    if(input.state && input.state->cachedInternalx.size() == 2)
        cachedWarmstart = input.state->cachedInternalx;
    if(input.state && input.state->lastInternalx.size() == 2)
        branchWarmstart = input.state->lastInternalx;

    GlobalizedSolidSolutionDefaultCandidateOptions candidateOptions;
    candidateOptions.branchTolerance = options.branchTolerance;
    candidateOptions.cachedStatePriority = -1.0;
    candidateOptions.preferredBranchPriority = -options.branchScoreHysteresis;
    candidateOptions.requireCachedStateWarmstart = true;
    candidateOptions.sourceKey = "GlobalizedBinaryRedlichKister::CandidateSource";

    return DefaultGlobalizedSolidSolutionCandidates(
        input,
        branches,
        input.x,
        cachedWarmstart,
        branchWarmstart,
        candidateOptions);
}

auto selectBranch(
    GlobalizedBinaryRedlichKisterOptions const& options,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionInput const& input) -> BranchCandidate
{
    const auto defaultGenerator = [=](GlobalizedSolidSolutionInput const& screeningInput, Vec<GlobalizedSolidSolutionBranch> const& screeningBranches)
    {
        return defaultCandidates(options, screeningBranches, screeningInput);
    };

    const auto selection = ComposeGlobalizedSolidSolutionBranch(
        branches,
        input,
        options.candidateGenerator,
        defaultGenerator,
        options.stabilityCriterion,
        [=](GlobalizedSolidSolutionCandidate const& candidate, GlobalizedSolidSolutionInput const& evaluationInput, Vec<GlobalizedSolidSolutionBranch> const& evaluationBranches)
        {
            if(candidate.branch == GlobalizedSolidSolutionNoBranch || candidate.branch >= evaluationBranches.size())
                throw std::runtime_error("Globalized binary Redlich-Kister candidate generator returned an invalid branch.");

            const Optional<ArrayXr> branchWarmstart = candidate.initialInternalx.size() == 2
                ? Optional<ArrayXr>(candidate.initialInternalx)
                : std::nullopt;

            const auto [internalx, score] = minimizeBranchObjective(options, evaluationBranches[candidate.branch], evaluationInput.T, evaluationInput.x, branchWarmstart);

            GlobalizedSolidSolutionBranchSelection selection;
            selection.branch = candidate.branch;
            selection.score = score;
            selection.internalx = internalx;
            selection.usedWarmstart = branchWarmstart.has_value();
            if(const auto source = std::any_cast<String>(&candidate.extra.at("GlobalizedBinaryRedlichKister::CandidateSource")))
                selection.reusedState = (*source == "state-cache");
            return selection;
        },
        {},
        "Globalized binary Redlich-Kister candidate generator returned no candidates.",
        "Globalized binary Redlich-Kister candidate generator returned an invalid branch.",
        "Globalized binary Redlich-Kister stability screen rejected all branch candidates.");

    BranchCandidate best;
    best.branch = selection.branch;
    best.score = selection.score;
    best.internalx = selection.internalx;
    best.extra = selection.extra;
    best.reusedState = selection.reusedState;
    best.usedWarmstart = selection.usedWarmstart;
    best.stable = selection.stable;
    best.stabilityPenalty = selection.stabilityPenalty;
    best.splitRequest = selection.splitRequest;
    return best;
}

} // namespace

auto GlobalizedBinaryRedlichKisterModel(GlobalizedBinaryRedlichKisterOptions options) -> GlobalizedSolidSolutionModel
{
    options.branches = normalizedBranches(options.branches);

    return [=](GlobalizedSolidSolutionInput input)
    {
        if(input.x.size() != 2)
            throw std::runtime_error("Globalized binary Redlich-Kister model requires exactly two species.");

        auto state = input.state ? input.state : std::make_shared<GlobalizedSolidSolutionState>();
        const auto candidate = selectBranch(options, options.branches, input);
        const auto selectedBranch = candidate.branch;
        const auto& branch = options.branches[selectedBranch];
        ArrayXr internalx;
        if(candidate.internalx.size() == 2)
            internalx = candidate.internalx;
        else
            internalx = input.x;
        internalx[0] = std::clamp(static_cast<double>(internalx[0]), CompositionFloor, 1.0 - CompositionFloor);
        internalx[1] = 1.0 - internalx[0];
        const auto stateid = CurrentGlobalizedSolidSolutionChemicalPropsStateId(input.extra);

        const auto RT = universalGasConstant * input.T;
        const auto x1 = static_cast<double>(internalx[0]);
        const auto x2 = static_cast<double>(internalx[1]);

        GlobalizedSolidSolutionOutput output;
        output.branches = options.branches;
        output.selectedBranch = selectedBranch;
        output.branch = branch;
        output.state = state;
        output.som = StateOfMatter::Solid;
        output.splitRequest = candidate.splitRequest;
        output.Vxi = ArrayXr::Zero(2);
        output.ln_g = ArrayXr::Zero(2);
        output.ln_a = ArrayXr::Zero(2);

        output.ln_g[0] = x2*x2*(options.a0 + options.a1*(3*x1 - x2) + options.a2*(x1 - x2)*(5*x1 - x2));
        output.ln_g[1] = x1*x1*(options.a0 - options.a1*(3*x2 - x1) + options.a2*(x2 - x1)*(5*x2 - x1));
        output.ln_a = output.ln_g + log(internalx);

        output.Gx = (x1*x2*(options.a0 + options.a1*(x1 - x2) + options.a2*pow((x1 - x2), 2))) * RT;
        output.Hx = output.Gx;

        const auto mismatchPenalty = options.externalCompositionPenalty * (internalx - input.x).matrix().squaredNorm();
        output.Gx += mismatchPenalty * RT;
        output.Hx = output.Gx;

        state->chemicalPropsStateId = stateid;
        state->cachedBranchForState = selectedBranch;
        state->cachedInternalx = internalx;
        state->selectedBranch = selectedBranch;
        state->numEvaluations += 1;
        state->lastT = input.T;
        state->lastP = input.P;
        state->lastx = input.x;
        state->lastInternalx = internalx;
        state->data["GlobalizedBinaryRedlichKister::UsedWarmstart"] = candidate.usedWarmstart;
        state->data["GlobalizedBinaryRedlichKister::ReusedStateCache"] = candidate.reusedState;
        state->data["GlobalizedBinaryRedlichKister::BranchId"] = branch.id;
        state->data["GlobalizedBinaryRedlichKister::BranchLabel"] = branch.label;

        output.extra["GlobalizedBinaryRedlichKister::UsedWarmstart"] = candidate.usedWarmstart;
        output.extra["GlobalizedBinaryRedlichKister::ReusedStateCache"] = candidate.reusedState;
        output.extra["GlobalizedBinaryRedlichKister::BranchId"] = branch.id;
        output.extra["GlobalizedBinaryRedlichKister::BranchLabel"] = branch.label;
        output.extra["GlobalizedBinaryRedlichKister::BranchScore"] = candidate.score;
        output.extra["GlobalizedBinaryRedlichKister::InternalComposition"] = internalx;
        output.extra["GlobalizedBinaryRedlichKister::CompositionMismatchPenalty"] = mismatchPenalty;
        output.extra["GlobalizedSolidSolution::CandidateStable"] = candidate.stable;
        output.extra["GlobalizedSolidSolution::CandidateStabilityPenalty"] = candidate.stabilityPenalty;
        output.extra["GlobalizedSolidSolution::SplitRequested"] = candidate.splitRequest.requested;
        output.extra["GlobalizedSolidSolution::SplitRequest"] = candidate.splitRequest;
        if(candidate.splitRequest.requested)
            output.extra["GlobalizedSolidSolution::SplitReason"] = candidate.splitRequest.reason;
        for(const auto& [key, value] : candidate.extra)
            output.extra[key] = value;

        return output;
    };
}

auto ActivityModelGlobalizedBinaryRedlichKister(GlobalizedBinaryRedlichKisterOptions options) -> ActivityModelGenerator
{
    return ActivityModelGlobalizedSolidSolution(GlobalizedBinaryRedlichKisterModel(options));
}

auto DuplicateGlobalizedBinaryRedlichKisterPhaseBranches(
    Phase const& phase,
    GlobalizedBinaryRedlichKisterOptions options,
    String suffixSeparator) -> PhaseList
{
    const auto branches = normalizedBranches(options.branches);
    options.branches = branches;
    return DuplicateGlobalizedSolidSolutionPhaseBranches(
        phase,
        GlobalizedBinaryRedlichKisterModel(options),
        branches,
        suffixSeparator);
}

auto GlobalizedBinaryRedlichKisterSolidPhases(
    Database const& db,
    String name,
    Strings const& species,
    GlobalizedBinaryRedlichKisterOptions options,
    String suffixSeparator) -> PhaseList
{
    Vec<Species> phaseSpecies;
    phaseSpecies.reserve(species.size());
    for(const auto& speciesName : species)
        phaseSpecies.push_back(db.species().get(speciesName));

    Phase phase;
    phase = phase.withName(name);
    phase = phase.withSpecies(phaseSpecies);
    phase = phase.withStateOfMatter(StateOfMatter::Solid);
    phase = phase.withActivityModel(ActivityModelGlobalizedBinaryRedlichKister(options)(phase.species()));
    phase = phase.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(phase.species()));

    return DuplicateGlobalizedBinaryRedlichKisterPhaseBranches(phase, options, suffixSeparator);
}

} // namespace Reaktoro