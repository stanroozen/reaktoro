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

// C++ includes
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Species.hpp>

#include "ActivityModelGlobalizedSolidSolution.hpp"

namespace Reaktoro {

namespace {

auto globalizedPhaseScopePrefix(String const& phaseScope) -> String
{
    return "GlobalizedSolidSolution::Phase::" + phaseScope + "::";
}

auto describeGlobalizedSolidSolutionDiagnosticValue(std::any const& value) -> String
{
    if(const auto rendered = std::any_cast<String>(&value))
        return *rendered;
    if(const auto rendered = std::any_cast<const char*>(&value))
        return *rendered ? String(*rendered) : String();
    if(const auto rendered = std::any_cast<bool>(&value))
        return *rendered ? "true" : "false";
    if(const auto rendered = std::any_cast<double>(&value))
        return std::to_string(*rendered);
    if(const auto rendered = std::any_cast<std::uint64_t>(&value))
        return std::to_string(*rendered);
    if(const auto rendered = std::any_cast<std::int64_t>(&value))
        return std::to_string(*rendered);
    if(const auto rendered = std::any_cast<int>(&value))
        return std::to_string(*rendered);
    return {};
}

auto withPhaseName(GlobalizedSolidSolutionSplitRequest splitRequest, String const& phaseScope) -> GlobalizedSolidSolutionSplitRequest
{
    if(splitRequest.phaseName.empty())
        splitRequest.phaseName = phaseScope;
    return splitRequest;
}

auto splitRequestKey(String const& phaseScope) -> String
{
    return globalizedPhaseScopePrefix(phaseScope) + "SplitRequest";
}

auto splitRequestedKey(String const& phaseScope) -> String
{
    return globalizedPhaseScopePrefix(phaseScope) + "SplitRequested";
}

auto annotateCandidateSource(
    GlobalizedSolidSolutionCandidate& candidate,
    GlobalizedSolidSolutionDefaultCandidateOptions const& options,
    String const& source) -> void
{
    if(!options.sourceKey.empty())
        candidate.extra[options.sourceKey] = source;
}

auto duplicatePhaseNameMatches(String const& phaseName, String const& baseName, String const& suffixSeparator) -> bool
{
    if(phaseName == baseName)
        return true;

    const auto prefix = baseName + suffixSeparator;
    return phaseName.size() > prefix.size() && phaseName.compare(0, prefix.size(), prefix) == 0;
}

auto findSplitRequest(
    Map<String, Any> const& extra,
    GlobalizedSolidSolutionPhaseDefinition const& definition) -> std::optional<GlobalizedSolidSolutionSplitRequest>
{
    const auto scopedKey = splitRequestKey(definition.phaseName);
    const auto it = extra.find(scopedKey);
    if(it != extra.end())
    {
        if(const auto split = std::any_cast<GlobalizedSolidSolutionSplitRequest>(&it->second))
            return *split;
    }

    const auto generic = extra.find("GlobalizedSolidSolution::SplitRequest");
    if(generic != extra.end())
    {
        if(const auto split = std::any_cast<GlobalizedSolidSolutionSplitRequest>(&generic->second))
        {
            if(split->phaseName.empty() || split->phaseName == definition.phaseName)
                return *split;
        }
    }

    const auto stateIt = extra.find("GlobalizedSolidSolution::State");
    if(stateIt != extra.end())
    {
        if(const auto state = std::any_cast<SharedPtr<GlobalizedSolidSolutionState>>(&stateIt->second))
        {
            const auto& split = (*state)->lastSplitRequest;
            if(split.requested && (split.phaseName.empty() || split.phaseName == definition.phaseName))
                return split;
        }
    }

    return std::nullopt;
}

auto projectGlobalizedInternalCoordinates(
    ArrayXr x,
    ArrayXrConstRef lower,
    ArrayXrConstRef upper,
    bool enforceUnityConstraint,
    real tolerance) -> ArrayXr
{
    for(Index i = 0; i < x.size(); ++i)
    {
        auto value = static_cast<double>(x[i]);
        const auto loweri = static_cast<double>(lower[i]);
        const auto upperi = static_cast<double>(upper[i]);
        if(value < loweri)
            value = loweri;
        if(value > upperi)
            value = upperi;
        x[i] = value;
    }

    if(!enforceUnityConstraint)
        return x;

    for(Index iter = 0; iter < x.size() * 8; ++iter)
    {
        const auto residual = 1.0 - static_cast<double>(x.sum());
        if(std::abs(residual) <= tolerance)
            break;

        real capacity = 0.0;
        for(Index i = 0; i < x.size(); ++i)
        {
            const auto slack = residual > 0.0 ? upper[i] - x[i] : x[i] - lower[i];
            if(slack > tolerance)
                capacity += slack;
        }

        if(capacity <= tolerance)
            break;

        for(Index i = 0; i < x.size(); ++i)
        {
            const auto slack = residual > 0.0 ? upper[i] - x[i] : x[i] - lower[i];
            if(slack <= tolerance)
                continue;

            const auto delta = residual * (slack / capacity);
            auto candidate = static_cast<double>(x[i]) + delta;
            const auto loweri = static_cast<double>(lower[i]);
            const auto upperi = static_cast<double>(upper[i]);
            if(candidate < loweri)
                candidate = loweri;
            if(candidate > upperi)
                candidate = upperi;
            x[i] = candidate;
        }
    }

    return x;
}

auto fetchGlobalizedSolidSolutionState(Map<String, Any> const& extra) -> SharedPtr<GlobalizedSolidSolutionState>
{
    const auto it = extra.find("GlobalizedSolidSolution::State");
    if(it == extra.end() || !it->second.has_value())
        return {};

    try
    {
        return std::any_cast<SharedPtr<GlobalizedSolidSolutionState>>(it->second);
    }
    catch(const std::bad_any_cast&)
    {
        return {};
    }
}

auto mergeGlobalizedSolidSolutionExtra(
    Map<String, Any> extra,
    GlobalizedSolidSolutionOutput const& output,
    SharedPtr<GlobalizedSolidSolutionState> const& state,
    String const& phaseScope) -> Map<String, Any>
{
    for(const auto& [key, value] : output.extra)
        extra[key] = value;

    extra["GlobalizedSolidSolution::State"] = state;

    if(output.selectedBranch != GlobalizedSolidSolutionNoBranch)
        extra["GlobalizedSolidSolution::SelectedBranchIndex"] = static_cast<std::uint64_t>(output.selectedBranch);

    if(!output.branch.id.empty())
        extra["GlobalizedSolidSolution::SelectedBranchId"] = output.branch.id;

    if(!output.branch.label.empty())
        extra["GlobalizedSolidSolution::SelectedBranchLabel"] = output.branch.label;

    if(output.branch.lowerBounds.size())
        extra["GlobalizedSolidSolution::SelectedBranchLowerBounds"] = output.branch.lowerBounds;

    if(output.branch.upperBounds.size())
        extra["GlobalizedSolidSolution::SelectedBranchUpperBounds"] = output.branch.upperBounds;

    extra["GlobalizedSolidSolution::BranchCount"] = static_cast<std::uint64_t>(output.branches.size());

    auto splitRequest = withPhaseName(output.splitRequest, phaseScope);
    state->lastSplitRequest = splitRequest;
    extra["GlobalizedSolidSolution::SplitRequested"] = splitRequest.requested;
    extra["GlobalizedSolidSolution::SplitRequest"] = splitRequest;

    if(!phaseScope.empty())
    {
        extra[splitRequestedKey(phaseScope)] = splitRequest.requested;
        extra[splitRequestKey(phaseScope)] = splitRequest;
    }

    return extra;
}

auto ActivityModelGlobalizedSolidSolutionImpl(GlobalizedSolidSolutionModel model, Index requestedBranch, String phaseScope) -> ActivityModelGenerator
{
    return [=](SpeciesList const& species)
    {
        return [=](ActivityPropsRef props, ActivityModelArgs args)
        {
            auto state = fetchGlobalizedSolidSolutionState(props.extra);
            if(!state)
                state = std::make_shared<GlobalizedSolidSolutionState>();

            const auto output = model({args.T, args.P, args.x, props.extra, state, requestedBranch});
            const auto outputState = output.state ? output.state : state;

            if(output.selectedBranch != GlobalizedSolidSolutionNoBranch)
            {
                outputState->selectedBranch = output.selectedBranch;
                outputState->cachedBranchForState = output.selectedBranch;
            }

            if(output.ln_g.size() != species.size())
                throw std::runtime_error("ActivityModelGlobalizedSolidSolution received ln_g with invalid size.");

            if(output.ln_a.size() != species.size())
                throw std::runtime_error("ActivityModelGlobalizedSolidSolution received ln_a with invalid size.");

            if(output.Vxi.size() != species.size())
                throw std::runtime_error("ActivityModelGlobalizedSolidSolution received Vxi with invalid size.");

            props.Vx = output.Vx;
            props.VxT = output.VxT;
            props.VxP = output.VxP;
            props.Vxi = output.Vxi;
            props.Gx = output.Gx;
            props.Hx = output.Hx;
            props.Cpx = output.Cpx;
            props.ln_g = output.ln_g;
            props.ln_a = output.ln_a;
            props.som = output.som;
            props.extra = mergeGlobalizedSolidSolutionExtra(props.extra, output, outputState, phaseScope);
        };
    };
}

} // namespace

auto MinimizeGlobalizedSolidSolutionInternalProblem(
    GlobalizedSolidSolutionInternalProblem const& problem) -> GlobalizedSolidSolutionInternalResult
{
    if(!problem.objective)
        throw std::runtime_error("Globalized solid-solution internal minimization requires an objective.");

    const auto size = problem.initialx.size();
    if(size == 0)
        throw std::runtime_error("Globalized solid-solution internal minimization requires at least one coordinate.");
    if(problem.lowerBounds.size() != size || problem.upperBounds.size() != size)
        throw std::runtime_error("Globalized solid-solution internal minimization bounds must match the coordinate size.");

    auto bestx = projectGlobalizedInternalCoordinates(
        problem.initialx,
        problem.lowerBounds,
        problem.upperBounds,
        problem.enforceUnityConstraint,
        problem.tolerance);

    auto bestObjective = problem.objective(bestx);
    auto step = static_cast<double>(problem.initialStep);
    Index iterations = 0;

    while(iterations < problem.maxIterations && step > static_cast<double>(problem.tolerance))
    {
        bool improved = false;
        auto candidatex = bestx;
        auto candidateObjective = bestObjective;

        if(problem.enforceUnityConstraint)
        {
            for(Index i = 0; i < size; ++i)
            {
                for(Index j = 0; j < size; ++j)
                {
                    if(i == j)
                        continue;

                    const auto delta = std::min({step,
                        static_cast<double>(problem.upperBounds[i] - bestx[i]),
                        static_cast<double>(bestx[j] - problem.lowerBounds[j])});

                    if(delta <= static_cast<double>(problem.tolerance))
                        continue;

                    auto trialx = bestx;
                    trialx[i] += delta;
                    trialx[j] -= delta;
                    trialx = projectGlobalizedInternalCoordinates(
                        trialx,
                        problem.lowerBounds,
                        problem.upperBounds,
                        true,
                        problem.tolerance);

                    const auto trialObjective = problem.objective(trialx);
                    if(trialObjective + static_cast<double>(problem.tolerance) < candidateObjective)
                    {
                        candidateObjective = trialObjective;
                        candidatex = trialx;
                        improved = true;
                    }
                }
            }
        }
        else
        {
            for(Index i = 0; i < size; ++i)
            {
                for(const auto direction : {-1.0, 1.0})
                {
                    auto trialx = bestx;
                    trialx[i] = std::clamp(
                        static_cast<double>(trialx[i]) + direction * step,
                        static_cast<double>(problem.lowerBounds[i]),
                        static_cast<double>(problem.upperBounds[i]));
                    const auto trialObjective = problem.objective(trialx);
                    if(trialObjective + static_cast<double>(problem.tolerance) < candidateObjective)
                    {
                        candidateObjective = trialObjective;
                        candidatex = trialx;
                        improved = true;
                    }
                }
            }
        }

        ++iterations;
        if(improved)
        {
            bestx = candidatex;
            bestObjective = candidateObjective;
        }
        else
        {
            step *= 0.5;
        }
    }

    return {
        bestx,
        bestObjective,
        iterations,
        step <= static_cast<double>(problem.tolerance),
    };
}

auto DefaultGlobalizedSolidSolutionBranches(Index numCoords) -> Vec<GlobalizedSolidSolutionBranch>
{
    GlobalizedSolidSolutionBranch branch;
    branch.id = "global";
    branch.label = "global";
    branch.lowerBounds = ArrayXr::Zero(numCoords);
    branch.upperBounds = ArrayXr::Ones(numCoords);
    return {branch};
}

auto NormalizeGlobalizedSolidSolutionBranches(
    Vec<GlobalizedSolidSolutionBranch> branches,
    Index numCoords,
    String const& errorContext) -> Vec<GlobalizedSolidSolutionBranch>
{
    if(branches.empty())
        return DefaultGlobalizedSolidSolutionBranches(numCoords);

    for(Index i = 0; i < branches.size(); ++i)
    {
        auto& branch = branches[i];
        if(branch.id.empty())
            branch.id = "branch" + std::to_string(i);
        if(branch.label.empty())
            branch.label = branch.id;
        if(branch.lowerBounds.size() == 0)
            branch.lowerBounds = ArrayXr::Zero(numCoords);
        if(branch.upperBounds.size() == 0)
            branch.upperBounds = ArrayXr::Ones(numCoords);
        if(branch.lowerBounds.size() != numCoords || branch.upperBounds.size() != numCoords)
            throw std::runtime_error(errorContext);
    }

    return branches;
}

auto GlobalizedSolidSolutionBranchViolation(
    ArrayXrConstRef x,
    GlobalizedSolidSolutionBranch const& branch,
    real tolerance) -> real
{
    real violation = 0.0;
    for(Index i = 0; i < x.size(); ++i)
    {
        const auto below = std::max(0.0, static_cast<double>(branch.lowerBounds[i] - x[i] - tolerance));
        const auto above = std::max(0.0, static_cast<double>(x[i] - branch.upperBounds[i] - tolerance));
        violation += below * below + above * above;
    }
    return violation;
}

auto GlobalizedSolidSolutionBranchContains(
    ArrayXrConstRef x,
    GlobalizedSolidSolutionBranch const& branch,
    real tolerance) -> bool
{
    return GlobalizedSolidSolutionBranchViolation(x, branch, tolerance) == 0.0;
}

auto MergeGlobalizedSolidSolutionSplitRequests(
    GlobalizedSolidSolutionSplitRequest lhs,
    GlobalizedSolidSolutionSplitRequest const& rhs) -> GlobalizedSolidSolutionSplitRequest
{
    lhs.requested = lhs.requested || rhs.requested;
    if(lhs.phaseName.empty())
        lhs.phaseName = rhs.phaseName;
    if(lhs.triggeringBranch == GlobalizedSolidSolutionNoBranch)
        lhs.triggeringBranch = rhs.triggeringBranch;
    for(auto branch : rhs.branches)
    {
        if(std::find(lhs.branches.begin(), lhs.branches.end(), branch) == lhs.branches.end())
            lhs.branches.push_back(branch);
    }
    for(auto const& branchId : rhs.branchIds)
    {
        if(std::find(lhs.branchIds.begin(), lhs.branchIds.end(), branchId) == lhs.branchIds.end())
            lhs.branchIds.push_back(branchId);
    }
    if(lhs.reason.empty())
        lhs.reason = rhs.reason;
    for(const auto& [key, value] : rhs.extra)
        lhs.extra[key] = value;
    return lhs;
}

auto CurrentGlobalizedSolidSolutionChemicalPropsStateId(Map<String, Any> const& extra) -> std::uint64_t
{
    const auto it = extra.find("Reaktoro::ChemicalProps::StateId");
    if(it == extra.end())
        return 0;

    if(const auto value = std::any_cast<std::uint64_t>(&it->second))
        return *value;

    return 0;
}

auto DefaultGlobalizedSolidSolutionCandidates(
    GlobalizedSolidSolutionInput const& input,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    ArrayXrConstRef cachedStateComposition,
    Optional<ArrayXr> cachedStateWarmstart,
    Optional<ArrayXr> branchWarmstart,
    GlobalizedSolidSolutionDefaultCandidateOptions const& options) -> Vec<GlobalizedSolidSolutionCandidate>
{
    Vec<GlobalizedSolidSolutionCandidate> candidates;

    const auto stateid = CurrentGlobalizedSolidSolutionChemicalPropsStateId(input.extra);
    const auto hasReusableCachedState = input.state
        && input.state->cachedBranchForState < branches.size()
        && (!options.requireCachedStateWarmstart || cachedStateWarmstart.has_value())
        && (stateid != 0
            ? input.state->chemicalPropsStateId == stateid
            : GlobalizedSolidSolutionBranchContains(cachedStateComposition, branches[input.state->cachedBranchForState], options.branchTolerance));

    if(hasReusableCachedState)
    {
        GlobalizedSolidSolutionCandidate candidate;
        candidate.branch = input.state->cachedBranchForState;
        if(cachedStateWarmstart.has_value())
            candidate.initialInternalx = *cachedStateWarmstart;
        candidate.priority = options.cachedStatePriority;
        annotateCandidateSource(candidate, options, options.stateCacheSource);
        candidates.push_back(std::move(candidate));
        return candidates;
    }

    if(input.requestedBranch != GlobalizedSolidSolutionNoBranch)
    {
        if(input.requestedBranch >= branches.size())
            throw std::runtime_error(options.invalidRequestedBranchMessage);

        GlobalizedSolidSolutionCandidate candidate;
        candidate.branch = input.requestedBranch;
        if(branchWarmstart.has_value())
            candidate.initialInternalx = *branchWarmstart;
        annotateCandidateSource(candidate, options, options.requestedBranchSource);
        candidates.push_back(std::move(candidate));
        return candidates;
    }

    Index preferredBranch = GlobalizedSolidSolutionNoBranch;
    if(input.state && input.state->selectedBranch < branches.size())
        preferredBranch = input.state->selectedBranch;

    candidates.reserve(branches.size());
    for(Index i = 0; i < branches.size(); ++i)
    {
        GlobalizedSolidSolutionCandidate candidate;
        candidate.branch = i;
        if(branchWarmstart.has_value() && (preferredBranch == i || GlobalizedSolidSolutionBranchContains(*branchWarmstart, branches[i], options.branchTolerance)))
            candidate.initialInternalx = *branchWarmstart;
        if(preferredBranch == i)
            candidate.priority = options.preferredBranchPriority;
        annotateCandidateSource(candidate, options, options.branchScreenSource);
        candidates.push_back(std::move(candidate));
    }

    return candidates;
}

auto DefaultGlobalizedSolidSolutionSplitRequest(
    GlobalizedSolidSolutionInput const& input,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    Index selectedBranch,
    real branchTolerance,
    String violationKey) -> GlobalizedSolidSolutionSplitRequest
{
    if(input.requestedBranch != GlobalizedSolidSolutionNoBranch || branches.size() < 2)
        return {};

    Vec<real> violations(branches.size(), 0.0);
    auto minViolation = std::numeric_limits<double>::infinity();
    for(Index i = 0; i < branches.size(); ++i)
    {
        violations[i] = GlobalizedSolidSolutionBranchViolation(input.x, branches[i], branchTolerance);
        minViolation = std::min(minViolation, static_cast<double>(violations[i]));
    }

    if(minViolation <= static_cast<double>(branchTolerance))
        return {};

    GlobalizedSolidSolutionSplitRequest splitRequest;
    for(Index i = 0; i < branches.size(); ++i)
    {
        if(std::abs(static_cast<double>(violations[i] - minViolation)) <= static_cast<double>(branchTolerance))
        {
            splitRequest.branches.push_back(i);
            splitRequest.branchIds.push_back(branches[i].id);
        }
    }

    if(splitRequest.branches.size() < 2)
        return {};

    splitRequest.requested = true;
    splitRequest.triggeringBranch = selectedBranch;
    splitRequest.reason = "external-composition-between-branches";
    if(!violationKey.empty())
        splitRequest.extra[violationKey] = minViolation;
    return splitRequest;
}

auto BranchAmbiguityGlobalizedSolidSolutionStabilityCriterion(
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionBranchAmbiguityStabilityOptions const& options) -> GlobalizedSolidSolutionCandidateStabilityCriterion
{
    return [=](GlobalizedSolidSolutionInput const&, GlobalizedSolidSolutionBranch const& branch, ArrayXrConstRef internalx, real score)
    {
        GlobalizedSolidSolutionCandidateStability stability;

        const auto currentViolation = static_cast<double>(GlobalizedSolidSolutionBranchViolation(internalx, branch, options.branchTolerance));
        if(!options.violationKey.empty())
            stability.extra[options.violationKey] = currentViolation;

        if(currentViolation <= static_cast<double>(options.branchTolerance))
            return stability;

        auto nearestViolation = std::numeric_limits<double>::infinity();
        Indices nearestBranches;
        Strings nearestBranchIds;
        auto triggeringBranch = GlobalizedSolidSolutionNoBranch;

        for(Index i = 0; i < branches.size(); ++i)
        {
            if(branch.id == branches[i].id)
                triggeringBranch = i;

            const auto violation = static_cast<double>(GlobalizedSolidSolutionBranchViolation(internalx, branches[i], options.branchTolerance));
            if(violation + static_cast<double>(options.ambiguousScoreTolerance) < nearestViolation)
            {
                nearestViolation = violation;
                nearestBranches = {i};
                nearestBranchIds = {branches[i].id};
                continue;
            }

            if(std::abs(violation - nearestViolation) <= static_cast<double>(options.ambiguousScoreTolerance))
            {
                nearestBranches.push_back(i);
                nearestBranchIds.push_back(branches[i].id);
            }
        }

        if(!options.nearestViolationKey.empty() && std::isfinite(nearestViolation))
            stability.extra[options.nearestViolationKey] = nearestViolation;
        if(!options.nearestBranchCountKey.empty())
            stability.extra[options.nearestBranchCountKey] = static_cast<std::uint64_t>(nearestBranches.size());

        const auto ambiguous = nearestBranches.size() >= 2
            && nearestViolation <= currentViolation + static_cast<double>(options.ambiguousScoreTolerance)
            && score <= currentViolation + static_cast<double>(options.ambiguousScoreTolerance);

        if(ambiguous)
        {
            stability.penalty = options.ambiguityPenalty;
            stability.reason = options.splitReason;
            if(options.requestSplitOnAmbiguity)
            {
                stability.splitRequest.requested = true;
                stability.splitRequest.triggeringBranch = triggeringBranch;
                stability.splitRequest.branches = nearestBranches;
                stability.splitRequest.branchIds = nearestBranchIds;
                stability.splitRequest.reason = options.splitReason;
                for(const auto& [key, value] : stability.extra)
                    stability.splitRequest.extra[key] = value;
            }
            else
            {
                stability.stable = false;
            }

            return stability;
        }

        stability.stable = false;
        stability.reason = options.unstableReason;
        return stability;
    };
}

auto NamedGlobalizedSolidSolutionStabilityCriterion(
    String const& name,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionBranchAmbiguityStabilityOptions const& options) -> GlobalizedSolidSolutionCandidateStabilityCriterion
{
    if(name.empty())
        return {};

    if(name == "branch-ambiguity")
        return BranchAmbiguityGlobalizedSolidSolutionStabilityCriterion(branches, options);

    if(name == "magemin-pilot-branch-ambiguity")
    {
        auto namedOptions = options;
        const GlobalizedSolidSolutionBranchAmbiguityStabilityOptions defaults;
        if(namedOptions.nearestViolationKey.empty())
            namedOptions.nearestViolationKey = "MAGEMinSolidSolutionPilot::NearestBranchViolation";
        if(namedOptions.nearestBranchCountKey.empty())
            namedOptions.nearestBranchCountKey = "MAGEMinSolidSolutionPilot::NearestBranchCount";
        if(namedOptions.ambiguousScoreTolerance == defaults.ambiguousScoreTolerance)
            namedOptions.ambiguousScoreTolerance = namedOptions.branchTolerance;
        if(namedOptions.splitReason == defaults.splitReason)
            namedOptions.splitReason = "stability-screen-between-branches";
        if(namedOptions.unstableReason == defaults.unstableReason)
            namedOptions.unstableReason = "internal-composition-outside-branch";

        auto shared = BranchAmbiguityGlobalizedSolidSolutionStabilityCriterion(branches, namedOptions);
        return [=](GlobalizedSolidSolutionInput const& input, GlobalizedSolidSolutionBranch const& branch, ArrayXrConstRef internalx, real score)
        {
            auto stability = shared(input, branch, internalx, score);

            if(!stability.splitRequest.requested && !stability.stable)
            {
                auto triggeringBranch = GlobalizedSolidSolutionNoBranch;
                for(Index i = 0; i < branches.size(); ++i)
                {
                    if(branches[i].id == branch.id)
                    {
                        triggeringBranch = i;
                        break;
                    }
                }

                auto splitRequest = DefaultGlobalizedSolidSolutionSplitRequest(
                    input,
                    branches,
                    triggeringBranch,
                    namedOptions.branchTolerance,
                    namedOptions.violationKey);

                if(splitRequest.requested)
                {
                    stability.reason = namedOptions.splitReason;
                    stability.splitRequest = std::move(splitRequest);
                    stability.splitRequest.reason = namedOptions.splitReason;
                    for(const auto& [key, value] : stability.extra)
                        stability.splitRequest.extra[key] = value;
                }

                if(!stability.splitRequest.requested && triggeringBranch != GlobalizedSolidSolutionNoBranch)
                {
                    const auto currentViolation = static_cast<double>(GlobalizedSolidSolutionBranchViolation(
                        internalx,
                        branches[triggeringBranch],
                        namedOptions.branchTolerance));

                    if(currentViolation > static_cast<double>(namedOptions.branchTolerance))
                    {
                        auto nearestAlternativeViolation = std::numeric_limits<double>::infinity();
                        Indices nearestAlternativeBranches;
                        Strings nearestAlternativeBranchIds;

                        for(Index i = 0; i < branches.size(); ++i)
                        {
                            if(i == triggeringBranch)
                                continue;

                            const auto violation = static_cast<double>(GlobalizedSolidSolutionBranchViolation(
                                internalx,
                                branches[i],
                                namedOptions.branchTolerance));

                            if(violation + static_cast<double>(namedOptions.ambiguousScoreTolerance) < nearestAlternativeViolation)
                            {
                                nearestAlternativeViolation = violation;
                                nearestAlternativeBranches = {i};
                                nearestAlternativeBranchIds = {branches[i].id};
                                continue;
                            }

                            if(std::abs(violation - nearestAlternativeViolation) <= static_cast<double>(namedOptions.ambiguousScoreTolerance))
                            {
                                nearestAlternativeBranches.push_back(i);
                                nearestAlternativeBranchIds.push_back(branches[i].id);
                            }
                        }

                        if(!nearestAlternativeBranches.empty())
                        {
                            stability.reason = namedOptions.splitReason;
                            stability.splitRequest.requested = true;
                            stability.splitRequest.triggeringBranch = triggeringBranch;
                            stability.splitRequest.branches = {triggeringBranch};
                            stability.splitRequest.branchIds = {branch.id};
                            stability.splitRequest.branches.insert(
                                stability.splitRequest.branches.end(),
                                nearestAlternativeBranches.begin(),
                                nearestAlternativeBranches.end());
                            stability.splitRequest.branchIds.insert(
                                stability.splitRequest.branchIds.end(),
                                nearestAlternativeBranchIds.begin(),
                                nearestAlternativeBranchIds.end());
                            stability.splitRequest.reason = namedOptions.splitReason;
                            for(const auto& [key, value] : stability.extra)
                                stability.splitRequest.extra[key] = value;
                        }
                    }
                }
            }

            if(!stability.extra.count("MAGEMinSolidSolutionPilot::BranchViolation")
                && !namedOptions.violationKey.empty()
                && stability.extra.count(namedOptions.violationKey))
                stability.extra["MAGEMinSolidSolutionPilot::BranchViolation"] = stability.extra.at(namedOptions.violationKey);
            return stability;
        };
    }

    throw std::runtime_error("Unsupported globalized solid-solution stability policy: " + name);
}

auto NamedGlobalizedSolidSolutionStabilityCriterion(
    NamedGlobalizedSolidSolutionStabilityPolicy policy,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionBranchAmbiguityStabilityOptions const& options) -> GlobalizedSolidSolutionCandidateStabilityCriterion
{
    switch(policy)
    {
    case NamedGlobalizedSolidSolutionStabilityPolicy::BranchAmbiguity:
        return NamedGlobalizedSolidSolutionStabilityCriterion("branch-ambiguity", branches, options);
    case NamedGlobalizedSolidSolutionStabilityPolicy::MAGEMinPilotBranchAmbiguity:
        return NamedGlobalizedSolidSolutionStabilityCriterion("magemin-pilot-branch-ambiguity", branches, options);
    }

    throw std::runtime_error("Unsupported typed globalized solid-solution stability policy.");
}

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
    String const& rejectedCandidatesMessage) -> GlobalizedSolidSolutionBranchSelection
{
    const auto generator = candidateGenerator ? candidateGenerator : defaultCandidateGenerator;
    const auto candidates = generator(input, branches);
    if(candidates.empty())
        throw std::runtime_error(emptyCandidatesMessage);

    GlobalizedSolidSolutionBranchSelection best;
    auto bestEffectiveScore = std::numeric_limits<double>::infinity();
    Vec<String> rejectedCandidateDiagnostics;

    for(auto candidate : candidates)
    {
        if(candidate.branch == GlobalizedSolidSolutionNoBranch || candidate.branch >= branches.size())
            throw std::runtime_error(invalidBranchMessage);

        auto selection = evaluator(candidate, input, branches);
        if(selection.branch == GlobalizedSolidSolutionNoBranch)
            selection.branch = candidate.branch;

        for(const auto& [key, value] : candidate.extra)
            selection.extra[key] = value;

        selection.splitRequest = MergeGlobalizedSolidSolutionSplitRequests(candidate.splitRequest, selection.splitRequest);

        auto effectiveScore = selection.score + candidate.priority;

        const auto internalx = selection.internalx.size() ? selection.internalx : ArrayXr(input.x);
        if(stabilityCriterion)
        {
            const auto stability = stabilityCriterion(input, branches[candidate.branch], internalx, selection.score);
            for(const auto& [key, value] : stability.extra)
                selection.extra[key] = value;

            selection.stable = stability.stable;
            selection.stabilityPenalty = stability.penalty;
            selection.extra["GlobalizedSolidSolution::CandidateStable"] = stability.stable;
            selection.extra["GlobalizedSolidSolution::CandidateStabilityPenalty"] = stability.penalty;
            if(!stability.reason.empty())
                selection.extra["GlobalizedSolidSolution::CandidateStabilityReason"] = stability.reason;
            selection.splitRequest = MergeGlobalizedSolidSolutionSplitRequests(selection.splitRequest, stability.splitRequest);

            if(!stability.stable && !selection.splitRequest.requested)
            {
                std::ostringstream diagnostics;
                diagnostics << "branch=" << (!branches[candidate.branch].id.empty() ? branches[candidate.branch].id : std::to_string(candidate.branch));
                diagnostics << " score=" << static_cast<double>(selection.score);
                if(!stability.reason.empty())
                    diagnostics << " reason=" << stability.reason;
                for(const auto& [key, value] : stability.extra)
                {
                    const auto rendered = describeGlobalizedSolidSolutionDiagnosticValue(value);
                    if(rendered.empty())
                        continue;
                    diagnostics << ' ' << key << '=' << rendered;
                }
                rejectedCandidateDiagnostics.push_back(diagnostics.str());
                continue;
            }

            effectiveScore += stability.penalty;
        }

        if(splitRequestGenerator)
            selection.splitRequest = MergeGlobalizedSolidSolutionSplitRequests(
                selection.splitRequest,
                splitRequestGenerator(input, branches, candidate.branch));

        selection.extra["GlobalizedSolidSolution::SplitRequested"] = selection.splitRequest.requested;
        if(selection.splitRequest.requested)
        {
            selection.extra["GlobalizedSolidSolution::SplitReason"] = selection.splitRequest.reason;
            selection.extra["GlobalizedSolidSolution::SplitRequest"] = selection.splitRequest;
        }

        if(effectiveScore < bestEffectiveScore)
        {
            best = selection;
            bestEffectiveScore = effectiveScore;
        }
    }

    if(!std::isfinite(bestEffectiveScore))
    {
        auto message = rejectedCandidatesMessage;
        if(!rejectedCandidateDiagnostics.empty())
        {
            message += " Diagnostics: ";
            for(Index i = 0; i < rejectedCandidateDiagnostics.size(); ++i)
            {
                if(i != 0)
                    message += " | ";
                message += rejectedCandidateDiagnostics[i];
            }
        }
        throw std::runtime_error(message);
    }

    return best;
}

auto ActivityModelGlobalizedSolidSolution(GlobalizedSolidSolutionModel model) -> ActivityModelGenerator
{
    return ActivityModelGlobalizedSolidSolutionImpl(model, GlobalizedSolidSolutionNoBranch, "");
}

auto ActivityModelGlobalizedSolidSolution(GlobalizedSolidSolutionModel model, String phaseScope) -> ActivityModelGenerator
{
    return ActivityModelGlobalizedSolidSolutionImpl(model, GlobalizedSolidSolutionNoBranch, phaseScope);
}

auto ActivityModelGlobalizedSolidSolution(GlobalizedSolidSolutionModel model, Index requestedBranch) -> ActivityModelGenerator
{
    return ActivityModelGlobalizedSolidSolutionImpl(model, requestedBranch, "");
}

auto ActivityModelGlobalizedSolidSolution(GlobalizedSolidSolutionModel model, Index requestedBranch, String phaseScope) -> ActivityModelGenerator
{
    return ActivityModelGlobalizedSolidSolutionImpl(model, requestedBranch, phaseScope);
}

auto DuplicateGlobalizedSolidSolutionPhaseBranches(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    String suffixSeparator) -> PhaseList
{
    if(branches.empty())
        throw std::runtime_error("DuplicateGlobalizedSolidSolutionPhaseBranches requires at least one branch.");

    PhaseList duplicated;
    for(Index i = 0; i < branches.size(); ++i)
    {
        const auto suffix = !branches[i].label.empty() ? branches[i].label : branches[i].id;
        const auto name = suffix.empty()
            ? phase.name() + suffixSeparator + std::to_string(i)
            : phase.name() + suffixSeparator + suffix;
        auto branchPhase = phase.clone();

        duplicated.append(
            branchPhase
                .withName(name)
                .withActivityModel(ActivityModelGlobalizedSolidSolution(model, i, phase.name())(phase.species())));
    }

    return duplicated;
}

auto MakeGlobalizedSolidSolutionPhaseDefinition(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    String suffixSeparator) -> GlobalizedSolidSolutionPhaseDefinition
{
    GlobalizedSolidSolutionPhaseDefinition definition;
    definition.phaseName = phase.name();
    definition.model = model;
    definition.branches = branches;
    definition.suffixSeparator = suffixSeparator;

    auto prototype = phase.clone();
    prototype = prototype.withActivityModel(ActivityModelGlobalizedSolidSolution(model, phase.name())(phase.species()));
    definition.prototype = prototype;
    return definition;
}

auto ApplyGlobalizedSolidSolutionSplitRequests(
    PhaseList const& phases,
    Map<String, Any> const& extra,
    Vec<GlobalizedSolidSolutionPhaseDefinition> const& definitions) -> PhaseList
{
    PhaseList rebuilt;
    Set<String> handled;

    for(auto const& phase : phases)
    {
        auto matched = false;
        for(auto const& definition : definitions)
        {
            if(!duplicatePhaseNameMatches(phase.name(), definition.phaseName, definition.suffixSeparator))
                continue;

            matched = true;
            if(!handled.insert(definition.phaseName).second)
                break;

            const auto splitRequest = findSplitRequest(extra, definition).value_or(GlobalizedSolidSolutionSplitRequest{});
            if(splitRequest.requested)
            {
                const auto duplicated = DuplicateGlobalizedSolidSolutionPhaseBranches(
                    definition.prototype,
                    definition.model,
                    definition.branches,
                    definition.suffixSeparator);
                for(auto const& duplicatedPhase : duplicated)
                    rebuilt.append(duplicatedPhase);
            }
            else
            {
                rebuilt.append(definition.prototype);
            }
            break;
        }

        if(!matched)
            rebuilt.append(phase);
    }

    return rebuilt;
}

auto ApplyGlobalizedSolidSolutionSplitRequests(
    ChemicalSystem const& system,
    Map<String, Any> const& extra,
    Vec<GlobalizedSolidSolutionPhaseDefinition> const& definitions) -> ChemicalSystem
{
    const auto phases = ApplyGlobalizedSolidSolutionSplitRequests(system.phases(), extra, definitions);
    return ChemicalSystem(system.database(), phases, system.reactions(), system.surfaces());
}

auto AssembleGlobalizedSolidSolutionCandidatePhases(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    Vec<SolidSolutionCandidateState> const& candidates,
    String suffixSeparator) -> PhaseList
{
    if(candidates.empty())
        throw std::runtime_error("AssembleGlobalizedSolidSolutionCandidatePhases: candidate list is empty.");

    PhaseList result;
    for(Index i = 0; i < candidates.size(); ++i)
    {
        const auto& candidate = candidates[i];
        const auto branchIndex = candidate.branch;

        if(branchIndex != GlobalizedSolidSolutionNoBranch && branchIndex >= branches.size())
            throw std::runtime_error("AssembleGlobalizedSolidSolutionCandidatePhases: "
                "candidate branch index " + std::to_string(branchIndex) + " is out of range.");

        // Determine phase name suffix: prefer candidate label, then branch label/id, then index.
        String suffix = candidate.label;
        if(suffix.empty() && branchIndex != GlobalizedSolidSolutionNoBranch)
        {
            const auto& branch = branches[branchIndex];
            suffix = !branch.label.empty() ? branch.label : branch.id;
        }
        if(suffix.empty())
            suffix = std::to_string(i);

        const auto name = phase.name() + suffixSeparator + suffix;

        auto candidatePhase = phase.clone();
        const auto requestedBranch = branchIndex != GlobalizedSolidSolutionNoBranch
            ? branchIndex
            : Index(0);

        result.append(
            candidatePhase
                .withName(name)
                .withActivityModel(
                    ActivityModelGlobalizedSolidSolution(model, requestedBranch, phase.name())(phase.species())));
    }

    return result;
}

auto AssembleGlobalizedSolidSolutionCandidatePhases(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    SolidSolutionCandidateGenerator const& generator,
    real referenceT,
    real referenceP,
    ArrayXrConstRef referencex,
    String suffixSeparator) -> PhaseList
{
    const auto candidates = generator(referenceT, referenceP, referencex);
    return AssembleGlobalizedSolidSolutionCandidatePhases(
        phase, model, branches, candidates, suffixSeparator);
}

} // namespace Reaktoro