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

#include "ActivityModelMAGEMinSolidSolutionPilot.hpp"

// C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

namespace {

using std::log;

constexpr auto CompositionFloor = 1.0e-12;

constexpr auto CandidateSeedTolerance = 1.0e-8;

constexpr auto ProjectedGradientArmijo = 1.0e-4;

constexpr auto ProjectedGradientAgreementTolerance = 1.0e-6;

constexpr auto BranchStabilityObjectiveTolerance = 1.0e-3;

constexpr auto BranchStabilitySeedGapTolerance = 5.0e-2;

constexpr auto SplitCandidateStatesKey = "MAGEMinSolidSolutionPilot::SplitCandidates";

constexpr auto SplitCandidateObjectiveGapKey = "MAGEMinSolidSolutionPilot::CompetingStableBranchObjectiveGap";

constexpr auto SplitCandidateCountKey = "MAGEMinSolidSolutionPilot::CompetingStableBranchCount";

constexpr auto InternalObjectiveKey = "MAGEMinSolidSolutionPilot::InternalObjective";

constexpr auto PrecomputedCandidateIndexKey = "MAGEMinSolidSolutionPilot::PrecomputedCandidateIndex";

constexpr auto BuiltinLegacyMinimizerStrategy = "legacy";

constexpr auto BuiltinProjectedGradientMinimizerStrategy = "projected-gradient";

struct ConstrainedTernaryMinimizationOutcome
{
    GlobalizedSolidSolutionInternalResult result;
    Map<String, Any> extra;
};

auto sb11OlivineThermo() -> MAGEMinImportedBinarySolutionThermoModel
{
    MAGEMinImportedBinarySolutionThermoModel thermo;
    thermo.modelId = "sb11_ol";
    thermo.endmember0 = "fo";
    thermo.endmember1 = "fa";
    thermo.W = 7813.22;
    thermo.idealSiteMultiplicity = 2.0;
    return thermo;
}

auto sb11WadsleyiteThermo() -> MAGEMinImportedBinarySolutionThermoModel
{
    MAGEMinImportedBinarySolutionThermoModel thermo;
    thermo.modelId = "sb11_wa";
    thermo.endmember0 = "mgwa";
    thermo.endmember1 = "fewa";
    thermo.W = 16747.18;
    thermo.idealSiteMultiplicity = 2.0;
    return thermo;
}

auto sb21SpinelThermo() -> MAGEMinImportedBinarySolutionThermoModel
{
    MAGEMinImportedBinarySolutionThermoModel thermo;
    thermo.modelId = "sb21_sp";
    thermo.endmember0 = "sp";
    thermo.endmember1 = "hc";
    thermo.W = -533.21;
    thermo.idealSiteMultiplicity = 2.0;
    return thermo;
}

auto sb21NALThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb21_nal";
    thermo.endmembers = {"mnal", "fnal", "nnal"};
    thermo.W01 = 0.0;
    thermo.W02 = -60781.47;
    thermo.W12 = -60781.47;
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        const auto siteA = (5.0/6.0)*y[0] + (5.0/6.0)*y[1] + 0.5*y[2];
        const auto siteB = (1.0/6.0)*y[0] + (1.0/6.0)*y[1] + 0.5*y[2];
        const auto sum = y[0] + y[1] + y[2];

        return universalGasConstant * T * (
            2.0*y[0]*log(y[0])
            + sum*log(sum)
            + 6.0*siteA*log(siteA)
            + 6.0*siteB*log(siteB)
            + 2.0*y[1]*log(y[1])
            + 2.0*y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        const auto siteA = (5.0/6.0)*y[0] + (5.0/6.0)*y[1] + 0.5*y[2];
        const auto siteB = (1.0/6.0)*y[0] + (1.0/6.0)*y[1] + 0.5*y[2];
        const auto sum = y[0] + y[1] + y[2];

        ArrayXr ln_a(3);
        ln_a[0] = 2.0*log(y[0]) + log(sum) + 5.0*log(siteA) + log(siteB);
        ln_a[1] = log(sum) + 5.0*log(siteA) + log(siteB) + 2.0*log(y[1]);
        ln_a[2] = log(sum) + 3.0*log(siteA) + 3.0*log(siteB) + 2.0*log(y[2]);
        return ln_a;
    };
    return thermo;
}

auto sb21CalcioferriteThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb21_cf";
    thermo.endmembers = {"mgcf", "fecf", "nacf"};
    thermo.W01 = 0.0;
    thermo.W02 = 60825.08;
    thermo.W12 = 60825.08;
    thermo.excessChemicalPotentials = [](ArrayXrConstRef y)
    {
        const ArrayXr volumes = (ArrayXr(3) << 1.0, 1.0, 4.4532).finished();
        const auto sumv = y.matrix().dot(volumes.matrix());

        ArrayXr phi(3);
        phi[0] = y[0]*volumes[0]/sumv;
        phi[1] = y[1]*volumes[1]/sumv;
        phi[2] = y[2]*volumes[2]/sumv;

        const ArrayXr interactions = (ArrayXr(3) << 0.0, 60825.08, 60825.08).finished();

        ArrayXr mu(3);
        for(Index i = 0; i < 3; ++i)
        {
            real Gex = 0.0;
            Index interaction = 0;
            for(Index j = 0; j < 3; ++j)
            {
                const auto tmp = ((i == j) ? 1.0 : 0.0) - static_cast<double>(phi[j]);
                for(Index k = j + 1; k < 3; ++k)
                {
                    const auto delta = ((i == k) ? 1.0 : 0.0) - static_cast<double>(phi[k]);
                    Gex -= tmp*delta*(interactions[interaction]*2.0*volumes[i]/(volumes[j] + volumes[k]));
                    ++interaction;
                }
            }
            mu[i] = Gex;
        }

        return mu;
    };
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        return universalGasConstant * T * (
            y[0]*log(y[0])
            + (y[0] + y[1])*log(y[0] + y[1])
            + y[1]*log(y[1])
            + 2.0*y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        ArrayXr ln_a(3);
        ln_a[0] = log(y[0]) + log(y[0] + y[1]);
        ln_a[1] = log(y[0] + y[1]) + log(y[1]);
        ln_a[2] = 2.0*log(y[2]);
        return ln_a;
    };
    return thermo;
}

auto sb11AkimotoiteThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb11_ak";
    thermo.endmembers = {"mgak", "feak", "co"};
    thermo.W01 = 0.0;
    thermo.W02 = 66000.0;
    thermo.W12 = 0.0;
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        return universalGasConstant * T * (
            2.0*y[0]*log(y[0])
            + (y[1] + y[2])*log(y[1] + y[2])
            + y[1]*log(y[1])
            + y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        ArrayXr ln_a(3);
        ln_a[0] = 2.0*log(y[0]);
        ln_a[1] = log(y[1] + y[2]) + log(y[1]);
        ln_a[2] = log(y[1] + y[2]) + log(y[2]);
        return ln_a;
    };
    return thermo;
}

auto sb11PerovskiteThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb11_pv";
    thermo.endmembers = {"mgpv", "fepv", "alpv"};
    thermo.W01 = 0.0;
    thermo.W02 = 116000.0;
    thermo.W12 = 0.0;
    thermo.excessChemicalPotentials = [](ArrayXrConstRef y)
    {
        const auto y0 = static_cast<double>(y[0]);
        const auto y1 = static_cast<double>(y[1]);
        const auto y2 = static_cast<double>(y[2]);
        const auto v0 = 1.0;
        const auto v1 = 1.0;
        const auto v2 = 0.39;
        const auto sumv = y0*v0 + y1*v1 + y2*v2;

        ArrayXr phi(3);
        phi[0] = y0*v0/sumv;
        phi[1] = y1*v1/sumv;
        phi[2] = y2*v2/sumv;

        ArrayXr mu(3);
        mu[0] = 2.0*116000.0*v0*phi[2]*(1.0 - phi[0])/(v0 + v2);
        mu[1] = -2.0*116000.0*v1*phi[0]*phi[2]/(v0 + v2);
        mu[2] = 2.0*116000.0*v2*phi[0]*(1.0 - phi[2])/(v0 + v2);
        return mu;
    };
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        return universalGasConstant * T * (
            2.0*y[0]*log(y[0])
            + (y[1] + y[2])*log(y[1] + y[2])
            + y[1]*log(y[1])
            + y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        ArrayXr ln_a(3);
        ln_a[0] = 2.0*log(y[0]);
        ln_a[1] = log(y[1] + y[2]) + log(y[1]);
        ln_a[2] = log(y[1] + y[2]) + log(y[2]);
        return ln_a;
    };
    return thermo;
}

auto sb11CalcioferriteThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb11_cf";
    thermo.endmembers = {"mgcf", "fecf", "nacf"};
    thermo.W01 = 0.0;
    thermo.W02 = 0.0;
    thermo.W12 = 0.0;
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        return universalGasConstant * T * (
            y[0]*log(y[0])
            + (y[0] + y[2])*log(y[0] + y[2])
            + 2.0*y[1]*log(y[1])
            + y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        ArrayXr ln_a(3);
        ln_a[0] = log(y[0]) + log(y[0] + y[2]);
        ln_a[1] = 2.0*log(y[1]);
        ln_a[2] = log(y[2]) + log(y[0] + y[2]);
        return ln_a;
    };
    return thermo;
}

auto normalizedMAGEMinPilotBranches(Vec<GlobalizedSolidSolutionBranch> branches, Index numCoords) -> Vec<GlobalizedSolidSolutionBranch>
{
    return NormalizeGlobalizedSolidSolutionBranches(
        std::move(branches),
        numCoords,
        "MAGEMin pilot branches must define bounds matching the visible composition size.");
}

auto projectSeedToBranch(ArrayXr seed, GlobalizedSolidSolutionBranch const& branch) -> ArrayXr
{
    if(seed.size() != branch.lowerBounds.size() || seed.size() != branch.upperBounds.size())
        return seed;

    for(Index i = 0; i < seed.size(); ++i)
        seed[i] = std::clamp(static_cast<double>(seed[i]), static_cast<double>(branch.lowerBounds[i]), static_cast<double>(branch.upperBounds[i]));

    for(Index iter = 0; iter < seed.size() * 8; ++iter)
    {
        const auto residual = 1.0 - static_cast<double>(seed.sum());
        if(std::abs(residual) <= CandidateSeedTolerance)
            break;

        real capacity = 0.0;
        for(Index i = 0; i < seed.size(); ++i)
        {
            const auto slack = residual > 0.0 ? branch.upperBounds[i] - seed[i] : seed[i] - branch.lowerBounds[i];
            if(slack > CandidateSeedTolerance)
                capacity += slack;
        }

        if(capacity <= CandidateSeedTolerance)
            break;

        for(Index i = 0; i < seed.size(); ++i)
        {
            const auto slack = residual > 0.0 ? branch.upperBounds[i] - seed[i] : seed[i] - branch.lowerBounds[i];
            if(slack <= CandidateSeedTolerance)
                continue;

            const auto delta = residual * (slack / capacity);
            seed[i] = std::clamp(
                static_cast<double>(seed[i] + delta),
                static_cast<double>(branch.lowerBounds[i]),
                static_cast<double>(branch.upperBounds[i]));
        }
    }

    seed = seed.max(CompositionFloor);
    seed /= seed.sum();
    return seed;
}

auto seedsEquivalent(ArrayXrConstRef lhs, ArrayXrConstRef rhs, real tolerance) -> bool
{
    return lhs.size() == rhs.size() && (lhs - rhs).cwiseAbs().maxCoeff() <= tolerance;
}

auto appendDistinctCandidate(
    Vec<GlobalizedSolidSolutionCandidate>& candidates,
    GlobalizedSolidSolutionCandidate candidate,
    String const& seedLabel,
    real extraPriority = 0.0) -> void
{
    if(candidate.initialInternalx.size() == 0)
        return;

    for(const auto& existing : candidates)
    {
        if(existing.branch != candidate.branch)
            continue;
        if(existing.initialInternalx.size() == candidate.initialInternalx.size()
            && seedsEquivalent(existing.initialInternalx, candidate.initialInternalx, CandidateSeedTolerance))
            return;
    }

    candidate.priority += extraPriority;
    candidate.extra["MAGEMinSolidSolutionPilot::CandidateSeedLabel"] = seedLabel;
    candidates.push_back(std::move(candidate));
}

struct MAGEMinTernaryCandidateSeedSpec
{
    ArrayXr seed;
    String label;
    real priority = 0.0;
};

auto dominantEndmemberSeed(Index size, Index dominant) -> ArrayXr
{
    ArrayXr seed = ArrayXr::Constant(size, CompositionFloor);
    seed[dominant] = 1.0 - (size - 1)*CompositionFloor;
    return seed;
}

auto binaryEdgeMidpointSeed(Index size, Index lhs, Index rhs) -> ArrayXr
{
    ArrayXr seed = ArrayXr::Constant(size, CompositionFloor);
    const auto retained = 1.0 - (size - 2)*CompositionFloor;
    seed[lhs] = 0.5*retained;
    seed[rhs] = 0.5*retained;
    return seed;
}

auto normalizedDominantEndmemberOrder(
    MAGEMinStructuredTernaryProposalOptions const& options,
    Index size) -> Indices
{
    Indices order;
    order.reserve(size);

    for(const auto index : options.dominantEndmemberOrder)
    {
        if(index >= size)
            continue;
        if(std::find(order.begin(), order.end(), index) == order.end())
            order.push_back(index);
    }

    for(Index i = 0; i < size; ++i)
    {
        if(std::find(order.begin(), order.end(), i) == order.end())
            order.push_back(i);
    }

    return order;
}

auto defaultTernaryCandidateSeedSpecs(
    MAGEMinImportedConstrainedTernarySolutionThermoModel const& thermo,
    MAGEMinStructuredTernaryProposalOptions const& options,
    ArrayXrConstRef visiblex) -> Vec<MAGEMinTernaryCandidateSeedSpec>
{
    Vec<MAGEMinTernaryCandidateSeedSpec> specs;

    if(options.includeVisibleCompositionSeed)
        specs.push_back({ArrayXr(visiblex), "visible-composition", options.visibleCompositionPriority});

    const auto dominantOrder = normalizedDominantEndmemberOrder(options, visiblex.size());

    auto appendSeed = [&](ArrayXr seed, String label, real priority)
    {
        specs.push_back({std::move(seed), std::move(label), priority});
    };

    if(options.includeDominantEndmemberSeeds)
    {
        for(Index i = 0; i < dominantOrder.size(); ++i)
        {
            const auto dominant = dominantOrder[i];
            appendSeed(
                dominantEndmemberSeed(visiblex.size(), dominant),
                "dominant::" + thermo.endmembers[dominant],
                options.dominantEndmemberPriority + i*options.dominantEndmemberPriorityStep);
        }
    }

    if(options.includeBinaryEdgeMidpointSeeds)
    {
        Index pairRank = 0;
        for(Index i = 0; i < dominantOrder.size(); ++i)
        {
            for(Index j = i + 1; j < dominantOrder.size(); ++j)
            {
                const auto lhs = dominantOrder[i];
                const auto rhs = dominantOrder[j];
                appendSeed(
                    binaryEdgeMidpointSeed(visiblex.size(), lhs, rhs),
                    "edge::" + thermo.endmembers[lhs] + "-" + thermo.endmembers[rhs],
                    options.binaryEdgePriority - pairRank*1.0e-7);
                pairRank += 1;
            }
        }
    }

    return specs;
}

auto augmentDefaultConstrainedTernaryCandidates(
    MAGEMinImportedConstrainedTernarySolutionThermoModel const& thermo,
    MAGEMinStructuredTernaryProposalOptions const& proposalOptions,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    ArrayXrConstRef visiblex,
    Vec<GlobalizedSolidSolutionCandidate> candidates) -> Vec<GlobalizedSolidSolutionCandidate>
{
    if(candidates.empty())
        return candidates;

    const auto sourceIt = candidates[0].extra.find("MAGEMinSolidSolutionPilot::CandidateSource");
    const auto source = sourceIt != candidates[0].extra.end()
        ? std::any_cast<String>(&sourceIt->second)
        : nullptr;

    if(!(source && (*source == "branch-screen" || *source == "requested-branch")))
        return candidates;

    const auto seedSpecs = defaultTernaryCandidateSeedSpecs(thermo, proposalOptions, visiblex);
    Strings generatedSeedLabels;
    generatedSeedLabels.reserve(seedSpecs.size());
    for(const auto& spec : seedSpecs)
        generatedSeedLabels.push_back(spec.label);

    Vec<GlobalizedSolidSolutionCandidate> augmented = candidates;
    for(const auto& candidate : candidates)
    {
        const auto& branch = branches[candidate.branch];
        for(const auto& spec : seedSpecs)
        {
            auto seeded = candidate;
            seeded.initialInternalx = projectSeedToBranch(spec.seed, branch);
            appendDistinctCandidate(augmented, std::move(seeded), spec.label, spec.priority);
        }
    }

    for(auto& candidate : augmented)
    {
        candidate.extra["MAGEMinSolidSolutionPilot::GeneratedCandidateCount"] = static_cast<std::uint64_t>(augmented.size());
        candidate.extra["MAGEMinSolidSolutionPilot::GeneratedCandidateSeedLabels"] = generatedSeedLabels;
    }

    return augmented;
}

auto defaultCandidates(
    MAGEMinSolidSolutionPilotBranchPolicyOptions const& options,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionInput const& input,
    ArrayXrConstRef seedx) -> Vec<GlobalizedSolidSolutionCandidate>
{
    Optional<ArrayXr> cachedWarmstart = std::nullopt;
    Optional<ArrayXr> branchWarmstart = std::nullopt;
    if(input.state && input.state->cachedInternalx.size() == seedx.size())
        cachedWarmstart = input.state->cachedInternalx;
    if(input.state && input.state->lastInternalx.size() == seedx.size())
        branchWarmstart = input.state->lastInternalx;

    GlobalizedSolidSolutionDefaultCandidateOptions candidateOptions;
    candidateOptions.branchTolerance = options.branchTolerance;
    candidateOptions.cachedStatePriority = -1.0;
    candidateOptions.preferredBranchPriority = -options.branchScoreHysteresis;
    candidateOptions.requireCachedStateWarmstart = true;
    candidateOptions.sourceKey = "MAGEMinSolidSolutionPilot::CandidateSource";
    candidateOptions.invalidRequestedBranchMessage = "Requested MAGEMin imported solid-solution pilot branch is out of range.";

    // During duplicated-phase requested-branch evaluations, disable cached-state reuse so
    // requested branch selection is not short-circuited by stale cached branch state.
    // Keep branch warmstart enabled to preserve solver stability.
    if(input.requestedBranch != GlobalizedSolidSolutionNoBranch)
    {
        cachedWarmstart = std::nullopt;
    }

    ArrayXr initialSeed = cachedWarmstart.has_value()
        ? ArrayXr(*cachedWarmstart)
        : ArrayXr(seedx);

    auto candidates = DefaultGlobalizedSolidSolutionCandidates(
        input,
        branches,
        initialSeed,
        cachedWarmstart,
        branchWarmstart,
        candidateOptions);

    for(auto& candidate : candidates)
        candidate.extra["MAGEMinSolidSolutionPilot::GeneratedCandidateCount"] = static_cast<std::uint64_t>(candidates.size());

    return candidates;
}

auto branchBoundsEqual(ArrayXrConstRef lhs, ArrayXrConstRef rhs, real tolerance) -> bool
{
    if(lhs.size() != rhs.size())
        return false;
    if(lhs.size() == 0)
        return true;
    return (lhs - rhs).cwiseAbs().maxCoeff() <= tolerance;
}

auto findBranchIndex(
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionBranch const& branch,
    real tolerance) -> Index
{
    for(Index i = 0; i < branches.size(); ++i)
    {
        if(!branch.id.empty() && branches[i].id == branch.id)
            return i;

        if(branchBoundsEqual(branches[i].lowerBounds, branch.lowerBounds, tolerance)
            && branchBoundsEqual(branches[i].upperBounds, branch.upperBounds, tolerance))
            return i;
    }

    return GlobalizedSolidSolutionNoBranch;
}

auto defaultPilotStabilityCriterion(
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    real branchTolerance,
    String const& splitViolationKey) -> GlobalizedSolidSolutionCandidateStabilityCriterion
{
    GlobalizedSolidSolutionBranchAmbiguityStabilityOptions options;
    options.branchTolerance = branchTolerance;
    options.violationKey = splitViolationKey;
    return NamedGlobalizedSolidSolutionStabilityCriterion(
        NamedGlobalizedSolidSolutionStabilityPolicy::MAGEMinPilotBranchAmbiguity,
        branches,
        options);
}

struct PrecomputedConstrainedTernaryCandidateEvaluation
{
    GlobalizedSolidSolutionBranchSelection selection;
};

auto branchCandidateLabel(GlobalizedSolidSolutionBranch const& branch, Index branchIndex) -> String
{
    if(!branch.label.empty())
        return branch.label;
    if(!branch.id.empty())
        return branch.id;
    return std::to_string(branchIndex);
}

auto candidateReusedState(GlobalizedSolidSolutionCandidate const& candidate) -> bool
{
    const auto it = candidate.extra.find("MAGEMinSolidSolutionPilot::CandidateSource");
    if(it == candidate.extra.end())
        return false;
    if(const auto source = std::any_cast<String>(&it->second))
        return (*source == "state-cache");
    return false;
}

auto solveConstrainedTernaryInternalProblem(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult;

auto solveConstrainedTernaryInternalProblemWithDiagnostics(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> ConstrainedTernaryMinimizationOutcome;

auto candidateInternalObjective(GlobalizedSolidSolutionBranchSelection const& selection) -> real
{
    const auto it = selection.extra.find(InternalObjectiveKey);
    if(it == selection.extra.end())
        return std::numeric_limits<double>::infinity();
    if(const auto objective = std::any_cast<real>(&it->second))
        return *objective;
    return std::numeric_limits<double>::infinity();
}

auto formatArray(ArrayXrConstRef values) -> String
{
    std::ostringstream out;
    out << "[";
    for(Index i = 0; i < values.size(); ++i)
    {
        if(i)
            out << ", ";
        out << static_cast<double>(values[i]);
    }
    out << "]";
    return out.str();
}

auto makeConstrainedTernarySelection(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    GlobalizedSolidSolutionCandidate const& candidate,
    GlobalizedSolidSolutionInput const& evaluationInput,
    Vec<GlobalizedSolidSolutionBranch> const& evaluationBranches,
    ArrayXrConstRef visiblex,
    real branchTolerance) -> GlobalizedSolidSolutionBranchSelection
{
    const Optional<ArrayXr> warmstart = candidate.initialInternalx.size() == visiblex.size()
        ? Optional<ArrayXr>(candidate.initialInternalx)
        : std::nullopt;

    const auto minimized = solveConstrainedTernaryInternalProblemWithDiagnostics(options, evaluationInput.T, visiblex, warmstart);
    ArrayXr internalx = minimized.result.x;
    internalx = internalx.max(CompositionFloor);
    internalx /= internalx.sum();
    auto requestedSeedProjectionApplied = false;

    if(evaluationInput.requestedBranch != GlobalizedSolidSolutionNoBranch
        && candidate.branch == evaluationInput.requestedBranch
        && candidate.branch < evaluationBranches.size())
    {
        auto requestedSeed = warmstart.has_value()
            ? projectSeedToBranch(ArrayXr(*warmstart), evaluationBranches[candidate.branch])
            : projectSeedToBranch(ArrayXr(visiblex), evaluationBranches[candidate.branch]);

        if(requestedSeed.size() == internalx.size())
        {
            const auto minimizedViolation = GlobalizedSolidSolutionBranchViolation(internalx, evaluationBranches[candidate.branch], branchTolerance);
            const auto seedViolation = GlobalizedSolidSolutionBranchViolation(requestedSeed, evaluationBranches[candidate.branch], branchTolerance);
            if(seedViolation + branchTolerance < minimizedViolation)
            {
                internalx = requestedSeed;
                requestedSeedProjectionApplied = true;
            }
        }
    }

    GlobalizedSolidSolutionBranchSelection selection;
    selection.branch = candidate.branch;
    selection.score = GlobalizedSolidSolutionBranchViolation(internalx, evaluationBranches[candidate.branch], branchTolerance);
    selection.internalx = internalx;
    selection.usedWarmstart = warmstart.has_value();
    selection.reusedState = candidateReusedState(candidate);
    selection.extra["MAGEMinSolidSolutionPilot::InternalMinimizerIterations"] = static_cast<std::uint64_t>(minimized.result.iterations);
    selection.extra["MAGEMinSolidSolutionPilot::InternalMinimizerConverged"] = minimized.result.converged;
    selection.extra[InternalObjectiveKey] = minimized.result.objective;
    selection.extra["MAGEMinSolidSolutionPilot::RequestedBranchProjectedSeedApplied"] = requestedSeedProjectionApplied;

    for(const auto& [key, value] : minimized.extra)
        selection.extra[key] = value;

    for(const auto& [key, value] : candidate.extra)
        selection.extra[key] = value;

    if(evaluationInput.requestedBranch != GlobalizedSolidSolutionNoBranch)
    {
        const auto seedUsed = warmstart.has_value() ? *warmstart : ArrayXr(visiblex);
        std::cerr
            << "[MAGEMinPilotDiag] model=" << options.thermo.modelId
            << " requestedBranch=" << evaluationInput.requestedBranch
            << " candidateBranch=" << candidate.branch
            << " selectedBranch=" << selection.branch
            << " seed=" << formatArray(seedUsed)
            << " internalx=" << formatArray(selection.internalx)
            << " objective=" << static_cast<double>(minimized.result.objective)
            << " projectedSeedApplied=" << (requestedSeedProjectionApplied ? 1 : 0)
            << "\n";
    }

    return selection;
}

auto constrainedTernarySplitCandidates(
    Vec<PrecomputedConstrainedTernaryCandidateEvaluation> const& evaluations,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionInput const& input,
    real branchTolerance) -> Vec<SolidSolutionCandidateState>
{
    struct BranchBestCandidate
    {
        bool valid = false;
        real objective = std::numeric_limits<double>::infinity();
        ArrayXr internalx;
    };

    Vec<BranchBestCandidate> best(branches.size());
    for(const auto& evaluation : evaluations)
    {
        const auto branchIndex = evaluation.selection.branch;
        if(branchIndex == GlobalizedSolidSolutionNoBranch || branchIndex >= branches.size())
            continue;

        const auto branchViolation = GlobalizedSolidSolutionBranchViolation(
            evaluation.selection.internalx,
            branches[branchIndex],
            branchTolerance);
        if(branchViolation > branchTolerance)
            continue;

        const auto objective = candidateInternalObjective(evaluation.selection);
        if(!best[branchIndex].valid || objective < best[branchIndex].objective)
        {
            best[branchIndex].valid = true;
            best[branchIndex].objective = objective;
            best[branchIndex].internalx = evaluation.selection.internalx;
        }
    }

    auto bestObjective = std::numeric_limits<double>::infinity();
    for(const auto& candidate : best)
        if(candidate.valid)
            bestObjective = std::min(bestObjective, static_cast<double>(candidate.objective));

    if(!std::isfinite(bestObjective))
        return {};

    Vec<SolidSolutionCandidateState> candidates;
    ArrayXr referenceSeed;
    for(Index i = 0; i < branches.size(); ++i)
    {
        if(!best[i].valid)
            continue;

        const auto objectiveGap = static_cast<double>(best[i].objective - bestObjective);
        if(objectiveGap > BranchStabilityObjectiveTolerance)
            continue;

        if(referenceSeed.size() == 0)
        {
            referenceSeed = best[i].internalx;
        }
        else if((best[i].internalx - referenceSeed).matrix().norm() <= BranchStabilitySeedGapTolerance)
        {
            continue;
        }

        candidates.push_back({i, best[i].internalx, objectiveGap, branchCandidateLabel(branches[i], i)});
    }

    if(candidates.size() < 2)
        return {};

    auto visibleBetweenBranches = false;
    for(const auto& candidate : candidates)
    {
        const auto branchViolation = GlobalizedSolidSolutionBranchViolation(input.x, branches[candidate.branch], branchTolerance);
        if(branchViolation > branchTolerance)
        {
            visibleBetweenBranches = true;
            break;
        }
    }

    return visibleBetweenBranches ? candidates : Vec<SolidSolutionCandidateState>{};
}

auto constrainedTernarySplitRequest(
    GlobalizedSolidSolutionInput const& input,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    Indices const& candidateBranches,
    Vec<SolidSolutionCandidateState> const& candidates,
    Index triggeringBranch,
    real branchTolerance,
    String const& splitViolationKey) -> GlobalizedSolidSolutionSplitRequest
{
    // Build the working seed list. Prefer minimizer-found seeds; if none are
    // within branch bounds (bulk in spinodal / two-phase region), synthesize seeds
    // by projecting the bulk composition onto each branch's feasible simplex.
    Vec<SolidSolutionCandidateState> workingCandidates = candidates;

    if(workingCandidates.size() < 2)
    {
        // Determine which branch indices to generate seeds for.
        // Use explicitly nominated candidateBranches if available; otherwise all branches.
        workingCandidates.clear();
        const auto n = input.x.size();
        const auto numBranches = static_cast<Index>(branches.size());
        const bool useCandidateBranches = static_cast<Index>(candidateBranches.size()) >= 2;

        for(Index bi = 0; bi < (useCandidateBranches ? static_cast<Index>(candidateBranches.size()) : numBranches); ++bi)
        {
            const auto branchIndex = useCandidateBranches ? candidateBranches[bi] : bi;
            if(branchIndex >= static_cast<Index>(branches.size()))
                continue;
            const auto& branch = branches[branchIndex];
            auto seed = projectSeedToBranch(ArrayXr(input.x), branch);

            const auto seedSum = static_cast<double>(seed.sum());
            if(seedSum <= 0.0)
                continue;
            seed /= seedSum;
            SolidSolutionCandidateState cs;
            cs.branch = branchIndex;
            cs.seedx = seed;
            cs.priority = 0.0;
            cs.label = branchCandidateLabel(branch, branchIndex) + "-projected";
            workingCandidates.push_back(std::move(cs));
        }
    }

    if(workingCandidates.size() < 2)
        return {};

    auto splitRequest = DefaultGlobalizedSolidSolutionSplitRequest(
        input,
        branches,
        triggeringBranch,
        branchTolerance,
        splitViolationKey);

    // Collect branch indices from workingCandidates (may differ from candidateBranches
    // when seeds were synthesized from all branches).
    Indices workingBranchIndices;
    workingBranchIndices.reserve(workingCandidates.size());
    for(const auto& wc : workingCandidates)
        workingBranchIndices.push_back(wc.branch);

    splitRequest.requested = true;
    splitRequest.triggeringBranch = triggeringBranch;
    splitRequest.branches = workingBranchIndices;
    splitRequest.branchIds.clear();
    for(const auto branchIndex : workingBranchIndices)
        if(branchIndex < static_cast<Index>(branches.size()))
            splitRequest.branchIds.push_back(branches[branchIndex].id);
    splitRequest.reason = "branch-stability-between-branches";
    splitRequest.extra[SplitCandidateStatesKey] = workingCandidates;
    splitRequest.extra[SplitCandidateCountKey] = static_cast<std::uint64_t>(workingCandidates.size());

    auto objectiveGap = 0.0;
    for(const auto& candidate : workingCandidates)
        objectiveGap = std::max(objectiveGap, static_cast<double>(candidate.priority));
    splitRequest.extra[SplitCandidateObjectiveGapKey] = objectiveGap;

    return splitRequest;
}


auto selectBranch(
    MAGEMinSolidSolutionPilotBranchPolicyOptions const& options,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionInput const& input,
    ArrayXrConstRef internalx,
    String const& emptyCandidatesMessage,
    String const& invalidBranchMessage,
    String const& rejectedCandidatesMessage,
    String const& splitViolationKey) -> GlobalizedSolidSolutionBranchSelection
{
    const auto defaultGenerator = [=](GlobalizedSolidSolutionInput const& screeningInput, Vec<GlobalizedSolidSolutionBranch> const& screeningBranches)
    {
        return defaultCandidates(options, screeningBranches, screeningInput, internalx);
    };

    const auto stabilityCriterion = options.stabilityCriterion
        ? options.stabilityCriterion
        : defaultPilotStabilityCriterion(branches, options.branchTolerance, splitViolationKey);

    return ComposeGlobalizedSolidSolutionBranch(
        branches,
        input,
        options.candidateGenerator,
        defaultGenerator,
        stabilityCriterion,
        [=](GlobalizedSolidSolutionCandidate const& candidate, GlobalizedSolidSolutionInput const&, Vec<GlobalizedSolidSolutionBranch> const& evaluationBranches)
        {
            GlobalizedSolidSolutionBranchSelection selection;
            selection.branch = candidate.branch;
            selection.score = GlobalizedSolidSolutionBranchViolation(internalx, evaluationBranches[candidate.branch], options.branchTolerance);
            selection.reusedState = candidateReusedState(candidate);
            return selection;
        },
        {},
        emptyCandidatesMessage,
        invalidBranchMessage,
        rejectedCandidatesMessage);
}

auto selectConstrainedTernaryBranch(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    MAGEMinSolidSolutionPilotBranchPolicyOptions const& branchPolicy,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionInput const& input,
    ArrayXrConstRef visiblex,
    String const& emptyCandidatesMessage,
    String const& invalidBranchMessage,
    String const& rejectedCandidatesMessage,
    String const& splitViolationKey) -> GlobalizedSolidSolutionBranchSelection
{
    const auto defaultGenerator = [=](GlobalizedSolidSolutionInput const& screeningInput, Vec<GlobalizedSolidSolutionBranch> const& screeningBranches)
    {
        return augmentDefaultConstrainedTernaryCandidates(
            options.thermo,
            options.proposals,
            screeningBranches,
            visiblex,
            defaultCandidates(branchPolicy, screeningBranches, screeningInput, visiblex));
    };

    const auto generator = branchPolicy.candidateGenerator
        ? branchPolicy.candidateGenerator
        : defaultGenerator;
    auto candidates = generator(input, branches);
    if(candidates.empty())
        throw std::runtime_error(emptyCandidatesMessage);

    Vec<PrecomputedConstrainedTernaryCandidateEvaluation> evaluations(candidates.size());
    for(Index i = 0; i < candidates.size(); ++i)
    {
        candidates[i].extra[PrecomputedCandidateIndexKey] = static_cast<std::uint64_t>(i);
        evaluations[i].selection = makeConstrainedTernarySelection(
            options,
            candidates[i],
            input,
            branches,
            visiblex,
            branchPolicy.branchTolerance);
    }

    const auto splitCandidates = constrainedTernarySplitCandidates(evaluations, branches, input, branchPolicy.branchTolerance);
    Indices splitCandidateBranches;
    splitCandidateBranches.reserve(splitCandidates.size());
    for(const auto& candidate : splitCandidates)
        splitCandidateBranches.push_back(candidate.branch);

    const auto stabilityCriterion = branchPolicy.stabilityCriterion
        ? branchPolicy.stabilityCriterion
        : defaultPilotStabilityCriterion(branches, branchPolicy.branchTolerance, splitViolationKey);

    return ComposeGlobalizedSolidSolutionBranch(
        branches,
        input,
        [candidates](GlobalizedSolidSolutionInput const&, Vec<GlobalizedSolidSolutionBranch> const&) { return candidates; },
        defaultGenerator,
        stabilityCriterion,
        [evaluations](GlobalizedSolidSolutionCandidate const& candidate, GlobalizedSolidSolutionInput const&, Vec<GlobalizedSolidSolutionBranch> const&)
        {
            const auto it = candidate.extra.find(PrecomputedCandidateIndexKey);
            if(it == candidate.extra.end())
                throw std::runtime_error("Precomputed constrained ternary candidate is missing its evaluation index.");
            const auto* index = std::any_cast<std::uint64_t>(&it->second);
            if(!index || *index >= evaluations.size())
                throw std::runtime_error("Precomputed constrained ternary candidate has an invalid evaluation index.");
            return evaluations[static_cast<Index>(*index)].selection;
        },
        [=](GlobalizedSolidSolutionInput const& splitInput, Vec<GlobalizedSolidSolutionBranch> const& splitBranches, Index triggeringBranch)
        {
            return constrainedTernarySplitRequest(
                splitInput,
                splitBranches,
                splitCandidateBranches,
                splitCandidates,
                triggeringBranch,
                branchPolicy.branchTolerance,
                splitViolationKey);
        },
        emptyCandidatesMessage,
        invalidBranchMessage,
        rejectedCandidatesMessage);
}

auto regularTernaryExcessChemicalPotentials(MAGEMinImportedConstrainedTernarySolutionThermoModel const& thermo, ArrayXrConstRef y) -> ArrayXr
{
    if(thermo.excessChemicalPotentials)
        return thermo.excessChemicalPotentials(y);

    ArrayXr mu(3);
    const auto y0 = static_cast<double>(y[0]);
    const auto y1 = static_cast<double>(y[1]);
    const auto y2 = static_cast<double>(y[2]);

    mu[0] = y1*(1.0 - y0)*thermo.W01 + y2*(1.0 - y0)*thermo.W02 - y1*y2*thermo.W12;
    mu[1] = y0*(1.0 - y1)*thermo.W01 - y0*y2*thermo.W02 + y2*(1.0 - y1)*thermo.W12;
    mu[2] = -y0*y1*thermo.W01 + y0*(1.0 - y2)*thermo.W02 + y1*(1.0 - y2)*thermo.W12;
    return mu;
}

auto normalizedInternalComposition(ArrayXr y) -> ArrayXr
{
    if(y.size() == 0)
        return y;

    for(Index iter = 0; iter < 8; ++iter)
    {
        y = y.max(CompositionFloor);
        const auto sum = static_cast<double>(y.sum());
        if(sum <= 0.0)
            y = ArrayXr::Constant(y.size(), 1.0/static_cast<double>(y.size()));
        else y /= sum;
    }

    return y;
}

auto constrainedTernaryObjective(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    ArrayXrConstRef y) -> real
{
    const auto muEx = regularTernaryExcessChemicalPotentials(options.thermo, y);
    const auto Gex = y.matrix().dot(muEx.matrix());
    const auto Gid = options.thermo.idealGibbs ? options.thermo.idealGibbs(T, y) : real(0.0);
    return Gex + Gid + options.externalCompositionPenalty*universalGasConstant*T*(y - x).matrix().squaredNorm();
}

auto constrainedTernaryGradient(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    ArrayXrConstRef y) -> ArrayXr
{
    const auto RT = universalGasConstant*T;
    ArrayXr gradient = regularTernaryExcessChemicalPotentials(options.thermo, y);
    if(options.thermo.idealLnActivities)
        gradient += RT*options.thermo.idealLnActivities(y);
    gradient += 2.0*options.externalCompositionPenalty*RT*(y - x);
    return gradient;
}

auto makeConstrainedTernaryLocalModel(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x) -> MAGEMinConstrainedTernaryLocalModel
{
    MAGEMinConstrainedTernaryLocalModel model;
    model.modelId = options.thermo.modelId;
    model.T = T;
    model.visiblex = ArrayXr(x);
    model.objective = [=](ArrayXrConstRef y) -> real { return constrainedTernaryObjective(options, T, x, y); };
    model.gradient = [=](ArrayXrConstRef y) -> ArrayXr { return constrainedTernaryGradient(options, T, x, y); };
    model.lowerBounds = ArrayXr::Constant(3, CompositionFloor);
    model.upperBounds = ArrayXr::Constant(3, 1.0 - CompositionFloor);
    model.enforceUnityConstraint = true;
    model.tolerance = options.minimizerTolerance;
    model.maxIterations = options.minimizerMaxIterations;
    return model;
}

auto projectedGradientConstrainedTernaryMinimizer(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
{
    auto current = normalizedInternalComposition(warmstart.value_or(ArrayXr(x)));
    auto currentObjective = constrainedTernaryObjective(options, T, x, current);

    Index iterations = 0;
    auto converged = false;

    for(; iterations < options.minimizerMaxIterations; ++iterations)
    {
        const auto gradient = constrainedTernaryGradient(options, T, x, current);
        auto projectedGradient = gradient.array() - gradient.mean();

        if(projectedGradient.cwiseAbs().maxCoeff() <= options.minimizerTolerance)
        {
            converged = true;
            break;
        }

        auto step = 0.25;
        auto accepted = false;
        const auto referenceSlope = projectedGradient.matrix().squaredNorm();

        for(Index backtrack = 0; backtrack < 20; ++backtrack)
        {
            auto trial = normalizedInternalComposition(current - step*projectedGradient);
            const auto trialObjective = constrainedTernaryObjective(options, T, x, trial);

            if(trialObjective <= currentObjective - ProjectedGradientArmijo*step*referenceSlope)
            {
                current = std::move(trial);
                currentObjective = trialObjective;
                accepted = true;
                break;
            }

            step *= 0.5;
        }

        if(!accepted)
        {
            converged = true; // Armijo failure: stuck at local minimum
            break;
        }
    }

    return {
        current,
        currentObjective,
        iterations,
        converged,
    };
}

auto defaultConstrainedTernaryMinimizer(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult;

auto constrainedTernaryBuiltinMinimizerOutcome(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> ConstrainedTernaryMinimizationOutcome
{
    ConstrainedTernaryMinimizationOutcome outcome;

    if(options.defaultMinimizerStrategy != BuiltinProjectedGradientMinimizerStrategy)
    {
        outcome.result = defaultConstrainedTernaryMinimizer(options, T, x, warmstart);
        outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String(BuiltinLegacyMinimizerStrategy);
        outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount"] = static_cast<std::uint64_t>(0);
        outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(outcome.result.iterations);
        return outcome;
    }

    const auto projected = projectedGradientConstrainedTernaryMinimizer(options, T, x, warmstart);
    outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount"] = static_cast<std::uint64_t>(projected.iterations);
    outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientConverged"] = projected.converged;

    if(!options.compareProjectedGradientAgainstLegacy)
    {
        outcome.result = projected;
        outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String(BuiltinProjectedGradientMinimizerStrategy);
        outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = projected.converged;
        outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(0);
        return outcome;
    }

    const auto legacy = defaultConstrainedTernaryMinimizer(options, T, x, warmstart);
    const auto sameShape = projected.x.size() == legacy.x.size();
    const auto compositionDelta = sameShape
        ? static_cast<double>((projected.x - legacy.x).cwiseAbs().maxCoeff())
        : std::numeric_limits<double>::infinity();
    const auto objectiveDelta = std::abs(static_cast<double>(projected.objective - legacy.objective));
    const auto agreement = projected.converged
        && sameShape
        && compositionDelta <= ProjectedGradientAgreementTolerance
        && objectiveDelta <= ProjectedGradientAgreementTolerance;
    const auto fallbackToLegacy = options.fallbackToLegacyOnProjectedGradientDisagreement && !agreement;

    outcome.result = fallbackToLegacy ? legacy : projected;
    outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String(fallbackToLegacy ? BuiltinLegacyMinimizerStrategy : BuiltinProjectedGradientMinimizerStrategy);
    outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = true;
    outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = fallbackToLegacy;
    outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = agreement;
    outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(legacy.iterations);
    outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientLegacyCompositionDelta"] = compositionDelta;
    outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientLegacyObjectiveDelta"] = objectiveDelta;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientHasLowerObjective"] = bool(projected.objective < legacy.objective);
    return outcome;
}

auto defaultConstrainedTernaryMinimizer(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
{
    GlobalizedSolidSolutionInternalProblem problem;
    problem.objective = [=](ArrayXrConstRef y) -> real { return constrainedTernaryObjective(options, T, x, y); };
    problem.initialx = warmstart.value_or(x);
    problem.lowerBounds = ArrayXr::Constant(3, CompositionFloor);
    problem.upperBounds = ArrayXr::Constant(3, 1.0 - CompositionFloor);
    problem.tolerance = options.minimizerTolerance;
    problem.maxIterations = options.minimizerMaxIterations;
    problem.enforceUnityConstraint = true;
    return MinimizeGlobalizedSolidSolutionInternalProblem(problem);
}

auto solveConstrainedTernaryInternalProblem(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
{
    if(options.localModelMinimizer)
        return options.localModelMinimizer(makeConstrainedTernaryLocalModel(options, T, x), warmstart);

    if(options.minimizer)
        return options.minimizer(options, T, x, warmstart);

    return defaultConstrainedTernaryMinimizer(options, T, x, warmstart);
}

auto solveConstrainedTernaryInternalProblemWithDiagnostics(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> ConstrainedTernaryMinimizationOutcome
{
    if(options.localModelMinimizer)
    {
        ConstrainedTernaryMinimizationOutcome outcome;
        const auto model = makeConstrainedTernaryLocalModel(options, T, x);
        outcome.result = options.localModelMinimizer(model, warmstart);
        outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String("custom-local-model");
        outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount"] = static_cast<std::uint64_t>(0);
        outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(0);

        if(options.localModelDiagnostics)
        {
            const auto payload = options.localModelDiagnostics(model, outcome.result);
            for(const auto& [key, value] : payload)
                outcome.extra[key] = value;
        }

        return outcome;
    }

    if(options.minimizer)
    {
        ConstrainedTernaryMinimizationOutcome outcome;
        outcome.result = options.minimizer(options, T, x, warmstart);
        outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String("custom");
        outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount"] = static_cast<std::uint64_t>(0);
        outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(0);
        return outcome;
    }

    return constrainedTernaryBuiltinMinimizerOutcome(options, T, x, warmstart);
}

} // namespace

auto MAGEMinProjectedGradientLocalModelMinimizer(
    MAGEMinConstrainedTernaryLocalModel const& model,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
{
    if(!model.gradient)
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: model.gradient callback must be populated.");

    auto current = normalizedInternalComposition(warmstart.value_or(ArrayXr(model.visiblex)));
    auto currentObjective = model.objective(current);

    Index iterations = 0;
    auto converged = false;

    for(; iterations < model.maxIterations; ++iterations)
    {
        const auto gradient = model.gradient(current);
        auto projectedGradient = gradient.array() - gradient.mean();

        if(projectedGradient.cwiseAbs().maxCoeff() <= model.tolerance)
        {
            converged = true;
            break;
        }

        auto step = 0.25;
        auto accepted = false;
        const auto referenceSlope = projectedGradient.matrix().squaredNorm();

        for(Index backtrack = 0; backtrack < 20; ++backtrack)
        {
            auto trial = normalizedInternalComposition(current - step*projectedGradient);
            const auto trialObjective = model.objective(trial);

            if(trialObjective <= currentObjective - ProjectedGradientArmijo*step*referenceSlope)
            {
                current = std::move(trial);
                currentObjective = trialObjective;
                accepted = true;
                break;
            }

            step *= 0.5;
        }

        if(!accepted)
        {
            converged = true; // Armijo failure: stuck at local minimum
            break;
        }
    }

    return {
        current,
        currentObjective,
        iterations,
        converged,
    };
}

auto MAGEMinSolidSolutionPilotModelImportedBinary(
    MAGEMinImportedBinarySolutionOptions options) -> GlobalizedSolidSolutionModel
{
    const auto thermo = options.thermo;
    auto branchPolicy = options.branchPolicy;
    branchPolicy.branches = normalizedMAGEMinPilotBranches(branchPolicy.branches, 2);

    return [=](GlobalizedSolidSolutionInput input)
    {
        if(input.x.size() != 2)
            throw std::runtime_error("MAGEMin imported binary pilot model requires exactly two species.");

        auto state = input.state ? input.state : std::make_shared<GlobalizedSolidSolutionState>();

        ArrayXr internalx = input.x;
        internalx[0] = std::clamp(static_cast<double>(internalx[0]), CompositionFloor, 1.0 - CompositionFloor);
        internalx[1] = 1.0 - internalx[0];

        const auto selected = selectBranch(
            branchPolicy,
            branchPolicy.branches,
            input,
            internalx,
            "MAGEMin imported binary pilot candidate generator returned no candidates.",
            "MAGEMin imported binary pilot candidate generator returned an invalid branch.",
            "MAGEMin imported binary pilot stability screen rejected all branch candidates.",
            "MAGEMinSolidSolutionPilot::SplitViolation");
        const auto selectedBranch = selected.branch;

        const auto RT = universalGasConstant * input.T;
        const auto x0 = static_cast<double>(internalx[0]);
        const auto x1 = static_cast<double>(internalx[1]);
        const auto muEx0 = x1*x1*thermo.W;
        const auto muEx1 = x0*x0*thermo.W;

        GlobalizedSolidSolutionOutput output;
        output.branches = branchPolicy.branches;
        output.selectedBranch = selectedBranch;
        output.branch = branchPolicy.branches[selectedBranch];
        output.state = state;
        output.som = StateOfMatter::Solid;
        output.Vxi = ArrayXr::Zero(2);
        output.ln_g = ArrayXr::Zero(2);
        output.ln_a = ArrayXr::Zero(2);

        output.ln_g[0] = muEx0/RT;
        output.ln_g[1] = muEx1/RT;
        output.ln_a = output.ln_g + thermo.idealSiteMultiplicity*log(internalx);

        output.Gx = x0*x1*thermo.W;
        output.Hx = output.Gx;
        output.splitRequest = selected.splitRequest;

        const auto stateid = CurrentGlobalizedSolidSolutionChemicalPropsStateId(input.extra);
        state->chemicalPropsStateId = stateid;
        state->selectedBranch = selectedBranch;
        state->cachedBranchForState = selectedBranch;
        state->cachedInternalx = internalx;
        state->numEvaluations += 1;
        state->lastT = input.T;
        state->lastP = input.P;
        state->lastx = input.x;
        state->lastInternalx = internalx;
        state->lastSplitRequest = output.splitRequest;
        state->data["MAGEMinSolidSolutionPilot::ModelId"] = thermo.modelId;
        state->data["MAGEMinSolidSolutionPilot::Endmember0"] = thermo.endmember0;
        state->data["MAGEMinSolidSolutionPilot::Endmember1"] = thermo.endmember1;

        output.extra["MAGEMinSolidSolutionPilot::ModelId"] = thermo.modelId;
        output.extra["MAGEMinSolidSolutionPilot::Endmember0"] = thermo.endmember0;
        output.extra["MAGEMinSolidSolutionPilot::Endmember1"] = thermo.endmember1;
        output.extra["MAGEMinSolidSolutionPilot::InternalComposition"] = internalx;
        output.extra["MAGEMinSolidSolutionPilot::BinaryInteractionParameterW"] = thermo.W;
        output.extra["MAGEMinSolidSolutionPilot::IdealSiteMultiplicity"] = thermo.idealSiteMultiplicity;
        output.extra["MAGEMinSolidSolutionPilot::UsedStateCache"] = selected.reusedState;
        for(const auto& [key, value] : selected.extra)
            output.extra[key] = value;
        output.extra["GlobalizedSolidSolution::SplitRequested"] = output.splitRequest.requested;
        output.extra["GlobalizedSolidSolution::SplitRequest"] = output.splitRequest;
        if(output.splitRequest.requested)
            output.extra["GlobalizedSolidSolution::SplitReason"] = output.splitRequest.reason;

        if(input.requestedBranch != GlobalizedSolidSolutionNoBranch)
        {
            std::cerr
                << "[MAGEMinPilotDiag] model=" << thermo.modelId
                << " requestedBranch=" << input.requestedBranch
                << " selectedBranch=" << output.selectedBranch
                << " finalInternalx=" << formatArray(internalx)
                << " splitRequested=" << (output.splitRequest.requested ? 1 : 0)
                << "\n";
        }

        return output;
    };
}

auto MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(
    MAGEMinImportedConstrainedTernarySolutionOptions options) -> GlobalizedSolidSolutionModel
{
    const auto thermo = options.thermo;
    auto branchPolicy = options.branchPolicy;
    branchPolicy.branches = normalizedMAGEMinPilotBranches(branchPolicy.branches, 3);

    return [=](GlobalizedSolidSolutionInput input)
    {
        if(input.x.size() != 3)
            throw std::runtime_error("MAGEMin imported constrained ternary pilot model requires exactly three species.");

        auto state = input.state ? input.state : std::make_shared<GlobalizedSolidSolutionState>();

        ArrayXr visiblex = input.x;
        visiblex = visiblex.max(CompositionFloor);
        visiblex /= visiblex.sum();

        const auto selected = selectConstrainedTernaryBranch(
            options,
            branchPolicy,
            branchPolicy.branches,
            input,
            visiblex,
            "MAGEMin imported constrained ternary pilot candidate generator returned no candidates.",
            "MAGEMin imported constrained ternary pilot candidate generator returned an invalid branch.",
            "MAGEMin imported constrained ternary pilot stability screen rejected all branch candidates.",
            "MAGEMinSolidSolutionPilot::SplitViolation");

        ArrayXr internalx = selected.internalx;
        internalx = internalx.max(CompositionFloor);
        internalx /= internalx.sum();

        const auto muEx = regularTernaryExcessChemicalPotentials(thermo, internalx);
        const auto RT = universalGasConstant * input.T;
        const auto Gex = internalx.matrix().dot(muEx.matrix());
        const auto Gid = thermo.idealGibbs ? thermo.idealGibbs(input.T, internalx) : real(0.0);
        const auto idealLnA = thermo.idealLnActivities ? thermo.idealLnActivities(internalx) : ArrayXr::Zero(3);

        GlobalizedSolidSolutionOutput output;
        output.branches = branchPolicy.branches;
        output.selectedBranch = selected.branch;
        output.branch = branchPolicy.branches[selected.branch];
        output.state = state;
        output.som = StateOfMatter::Solid;
        output.Vxi = ArrayXr::Zero(3);
        output.ln_g = muEx/RT;
        output.ln_a = output.ln_g + idealLnA;
        output.Gx = Gex + Gid;
        output.Hx = Gex;
        output.splitRequest = selected.splitRequest;

        const auto stateid = CurrentGlobalizedSolidSolutionChemicalPropsStateId(input.extra);
        state->chemicalPropsStateId = stateid;
        state->selectedBranch = selected.branch;
        state->cachedBranchForState = selected.branch;
        state->cachedInternalx = internalx;
        state->numEvaluations += 1;
        state->lastT = input.T;
        state->lastP = input.P;
        state->lastx = input.x;
        state->lastInternalx = internalx;
        state->lastSplitRequest = output.splitRequest;
        state->data["MAGEMinSolidSolutionPilot::ModelId"] = thermo.modelId;
        state->data["MAGEMinSolidSolutionPilot::Endmembers"] = thermo.endmembers;

        output.extra["MAGEMinSolidSolutionPilot::ModelId"] = thermo.modelId;
        output.extra["MAGEMinSolidSolutionPilot::Endmembers"] = thermo.endmembers;
        output.extra["MAGEMinSolidSolutionPilot::InternalComposition"] = internalx;
        output.extra["MAGEMinSolidSolutionPilot::UsedStateCache"] = selected.reusedState;
        output.extra["MAGEMinSolidSolutionPilot::UsedWarmstart"] = selected.usedWarmstart;
        output.extra["MAGEMinSolidSolutionPilot::BinaryInteraction01"] = thermo.W01;
        output.extra["MAGEMinSolidSolutionPilot::BinaryInteraction02"] = thermo.W02;
        output.extra["MAGEMinSolidSolutionPilot::BinaryInteraction12"] = thermo.W12;
        for(const auto& [key, value] : selected.extra)
            output.extra[key] = value;
        output.extra["GlobalizedSolidSolution::SplitRequested"] = output.splitRequest.requested;
        output.extra["GlobalizedSolidSolution::SplitRequest"] = output.splitRequest;
        if(output.splitRequest.requested)
            output.extra["GlobalizedSolidSolution::SplitReason"] = output.splitRequest.reason;

        return output;
    };
}

auto MAGEMinSolidSolutionPilotModelSB11Olivine(
    MAGEMinSB11OlivineOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedBinarySolutionOptions imported;
    imported.thermo = sb11OlivineThermo();
    imported.branchPolicy = options.branchPolicy;
    return MAGEMinSolidSolutionPilotModelImportedBinary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB11Wadsleyite(
    MAGEMinSB11WadsleyiteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedBinarySolutionOptions imported;
    imported.thermo = sb11WadsleyiteThermo();
    imported.branchPolicy = options.branchPolicy;
    return MAGEMinSolidSolutionPilotModelImportedBinary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB11Akimotoite(
    MAGEMinSB11AkimotoiteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb11AkimotoiteThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1};
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.compareProjectedGradientAgainstLegacy = true;
    imported.fallbackToLegacyOnProjectedGradientDisagreement = true;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB11Perovskite(
    MAGEMinSB11PerovskiteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb11PerovskiteThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1};
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.compareProjectedGradientAgainstLegacy = true;
    imported.fallbackToLegacyOnProjectedGradientDisagreement = true;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB11Calcioferrite(
    MAGEMinSB11CalcioferriteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb11CalcioferriteThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1};
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.compareProjectedGradientAgainstLegacy = true;
    imported.fallbackToLegacyOnProjectedGradientDisagreement = true;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21Spinel(
    MAGEMinSB21SpinelOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedBinarySolutionOptions imported;
    imported.thermo = sb21SpinelThermo();
    imported.branchPolicy = options.branchPolicy;
    return MAGEMinSolidSolutionPilotModelImportedBinary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21NAL(
    MAGEMinSB21NALOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb21NALThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1}; // nnal-first: largest W interactions
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.compareProjectedGradientAgainstLegacy = true;
    imported.fallbackToLegacyOnProjectedGradientDisagreement = true;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21Calcioferrite(
    MAGEMinSB21CalcioferriteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb21CalcioferriteThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1};
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotPhase(
    Phase const& phase,
    GlobalizedSolidSolutionModel model) -> Phase
{
    auto pilot = phase.clone();
    pilot = pilot.withActivityModel(ActivityModelGlobalizedSolidSolution(model, phase.name())(phase.species()));
    return pilot;
}

auto MAGEMinSolidSolutionPilotDefinition(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    MAGEMinSolidSolutionPilotOptions options) -> GlobalizedSolidSolutionPhaseDefinition
{
    return MakeGlobalizedSolidSolutionPhaseDefinition(
        MAGEMinSolidSolutionPilotPhase(phase, model),
        model,
        options.branches,
        options.suffixSeparator);
}

auto MAGEMinSolidSolutionPilotPhases(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    MAGEMinSolidSolutionPilotOptions options) -> PhaseList
{
    return DuplicateGlobalizedSolidSolutionPhaseBranches(
        MAGEMinSolidSolutionPilotPhase(phase, model),
        model,
        options.branches,
        options.suffixSeparator);
}

} // namespace Reaktoro