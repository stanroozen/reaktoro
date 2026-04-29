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

#include "EquilibriumUtils.hpp"

// C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {

namespace {

auto solveEquilibrium(
    ChemicalState& state,
    EquilibriumSpecs const& specs,
    EquilibriumConditions const& conditions,
    EquilibriumRestrictions const& restrictions,
    EquilibriumOptions const& options) -> EquilibriumResult;

auto solveClosedSystemEquilibrium(
    ChemicalState& state,
    EquilibriumRestrictions const& restrictions,
    EquilibriumOptions const& options,
    ArrayXdConstRef b0) -> EquilibriumResult
{
    EquilibriumSpecs specs(state.system());
    specs.temperature();
    specs.pressure();

    EquilibriumConditions conditions(specs);
    conditions.temperature(state.temperature());
    conditions.pressure(state.pressure());
    conditions.setInitialComponentAmounts(b0);

    return solveEquilibrium(state, specs, conditions, restrictions, options);
}

auto solveEquilibrium(
    ChemicalState& state,
    EquilibriumSpecs const& specs,
    EquilibriumConditions const& conditions,
    EquilibriumRestrictions const& restrictions,
    EquilibriumOptions const& options) -> EquilibriumResult
{
    EquilibriumOptions opts(options);

    EquilibriumSolver solver(specs);

    opts.use_ideal_activity_models = true; // force ideal activity models for the first computation
    solver.setOptions(opts);

    auto result = solver.solve(state, conditions, restrictions);

    // Skip the second computation if the first one using ideal activity models has already failed.
    if(result.failed())
        return result;

    opts.use_ideal_activity_models = options.use_ideal_activity_models; // for the second computation, use what user wants (maybe ideal model again, in which case the calculation below will converge immediately).
    solver.setOptions(opts);

    result += solver.solve(state, conditions, restrictions);

    return result;
}

auto extractBracketTokens(String const& input) -> Strings
{
    Strings tokens;
    std::size_t pos = 0;
    while(true)
    {
        const auto open = input.find('[', pos);
        if(open == String::npos)
            break;
        const auto close = input.find(']', open + 1);
        if(close == String::npos)
            break;
        tokens.push_back(input.substr(open + 1, close - open - 1));
        pos = close + 1;
    }
    return tokens;
}

auto rebuildEquilibriumSpecsForSystem(
    EquilibriumSpecs const& original,
    ChemicalSystem const& rebuiltSystem) -> EquilibriumSpecs
{
    EquilibriumSpecs rebuilt(rebuiltSystem);

    for(auto const& input : original.inputs())
    {
        if(input == "T") rebuilt.temperature();
        else if(input == "P") rebuilt.pressure();
        else if(input == "V") rebuilt.volume();
        else if(input == "U") rebuilt.internalEnergy();
        else if(input == "H") rebuilt.enthalpy();
        else if(input == "G") rebuilt.gibbsEnergy();
        else if(input == "A") rebuilt.helmholtzEnergy();
        else if(input == "S") rebuilt.entropy();
        else if(input == "charge") rebuilt.charge();
        else if(input == "pH") rebuilt.pH();
        else if(input == "pMg") rebuilt.pMg();
        else if(input == "pE") rebuilt.pE();
        else if(input == "Eh") rebuilt.Eh();
        else if(input.rfind("elementAmountInPhase[", 0) == 0)
        {
            const auto tokens = extractBracketTokens(input);
            if(tokens.size() != 2)
                throw std::runtime_error("Unsupported equilibrium input while rebuilding split-retry specs: " + input);
            rebuilt.elementAmountInPhase(tokens[0], tokens[1]);
        }
        else if(input.rfind("elementMassInPhase[", 0) == 0)
        {
            const auto tokens = extractBracketTokens(input);
            if(tokens.size() != 2)
                throw std::runtime_error("Unsupported equilibrium input while rebuilding split-retry specs: " + input);
            rebuilt.elementMassInPhase(tokens[0], tokens[1]);
        }
        else if(input.rfind("elementAmount[", 0) == 0)
        {
            const auto tokens = extractBracketTokens(input);
            if(tokens.size() != 1)
                throw std::runtime_error("Unsupported equilibrium input while rebuilding split-retry specs: " + input);
            rebuilt.elementAmount(tokens[0]);
        }
        else if(input.rfind("elementMass[", 0) == 0)
        {
            const auto tokens = extractBracketTokens(input);
            if(tokens.size() != 1)
                throw std::runtime_error("Unsupported equilibrium input while rebuilding split-retry specs: " + input);
            rebuilt.elementMass(tokens[0]);
        }
        else if(input.rfind("phaseAmount[", 0) == 0)
        {
            const auto tokens = extractBracketTokens(input);
            if(tokens.size() != 1)
                throw std::runtime_error("Unsupported equilibrium input while rebuilding split-retry specs: " + input);
            rebuilt.phaseAmount(tokens[0]);
        }
        else if(input.rfind("phaseMass[", 0) == 0)
        {
            const auto tokens = extractBracketTokens(input);
            if(tokens.size() != 1)
                throw std::runtime_error("Unsupported equilibrium input while rebuilding split-retry specs: " + input);
            rebuilt.phaseMass(tokens[0]);
        }
        else if(input.rfind("phaseVolume[", 0) == 0)
        {
            const auto tokens = extractBracketTokens(input);
            if(tokens.size() != 1)
                throw std::runtime_error("Unsupported equilibrium input while rebuilding split-retry specs: " + input);
            rebuilt.phaseVolume(tokens[0]);
        }
        else if(input.rfind("u[", 0) == 0)
        {
            const auto tokens = extractBracketTokens(input);
            if(tokens.size() != 1)
                throw std::runtime_error("Unsupported equilibrium input while rebuilding split-retry specs: " + input);
            rebuilt.chemicalPotential(tokens[0]);
        }
        else if(input.rfind("ln(a[", 0) == 0)
        {
            const auto tokens = extractBracketTokens(input);
            if(tokens.size() != 1)
                throw std::runtime_error("Unsupported equilibrium input while rebuilding split-retry specs: " + input);
            rebuilt.lnActivity(tokens[0]);
        }
        else if(input.rfind("f[", 0) == 0)
        {
            const auto tokens = extractBracketTokens(input);
            if(tokens.size() != 1)
                throw std::runtime_error("Unsupported equilibrium input while rebuilding split-retry specs: " + input);
            rebuilt.fugacity(tokens[0]);
        }
        else
        {
            throw std::runtime_error("Unsupported equilibrium input while rebuilding split-retry specs: " + input);
        }
    }

    for(auto const& pvar : original.namesControlVariablesP())
    {
        if(pvar == "T" || pvar == "P")
            continue;

        const auto currentPVars = rebuilt.namesControlVariablesP();
        if(std::find(currentPVars.begin(), currentPVars.end(), pvar) == currentPVars.end())
        {
            if(pvar.size() >= 2 && pvar.front() == '[' && pvar.back() == ']')
                rebuilt.openTo(pvar.substr(1, pvar.size() - 2));
            else
                rebuilt.openTo(pvar);
        }
    }

    return rebuilt;
}

auto rebuildEquilibriumConditionsForSystem(
    EquilibriumConditions const& original,
    EquilibriumSpecs const& rebuiltSpecs) -> EquilibriumConditions
{
    EquilibriumConditions rebuilt(rebuiltSpecs);

    const auto& names = original.inputNames();
    const auto& values = original.inputValues();
    for(Index i = 0; i < names.size(); ++i)
    {
        if(std::isnan(values[i].val()))
            continue;
        rebuilt.setInputVariable(names[i], values[i]);
    }

    if(original.initialComponentAmounts().size())
        rebuilt.setInitialComponentAmounts(original.initialComponentAmounts());

    if(original.lowerBoundsControlVariablesP().size())
        rebuilt.setLowerBoundsControlVariablesP(original.lowerBoundsControlVariablesP());
    if(original.upperBoundsControlVariablesP().size())
        rebuilt.setUpperBoundsControlVariablesP(original.upperBoundsControlVariablesP());

    return rebuilt;
}

auto rebuildEquilibriumRestrictionsForSystem(
    EquilibriumRestrictions const& original,
    ChemicalSystem const& rebuiltSystem) -> EquilibriumRestrictions
{
    EquilibriumRestrictions rebuilt(rebuiltSystem);

    const auto applyToMatchingSpecies = [&](Index sourceIndex, auto apply)
    {
        const auto name = original.system().species(sourceIndex).name();
        for(Index ispecies = 0; ispecies < rebuiltSystem.species().size(); ++ispecies)
            if(rebuiltSystem.species(ispecies).name() == name)
                apply(ispecies);
    };

    for(auto const& ispecies : original.speciesCannotIncrease())
        applyToMatchingSpecies(ispecies, [&](Index j){ rebuilt.cannotIncrease(j); });
    for(auto const& ispecies : original.speciesCannotDecrease())
        applyToMatchingSpecies(ispecies, [&](Index j){ rebuilt.cannotDecrease(j); });
    for(auto const& [ispecies, upper] : original.speciesCannotIncreaseAbove())
        applyToMatchingSpecies(ispecies, [&](Index j){ rebuilt.cannotIncreaseAbove(j, upper); });
    for(auto const& [ispecies, lower] : original.speciesCannotDecreaseBelow())
        applyToMatchingSpecies(ispecies, [&](Index j){ rebuilt.cannotDecreaseBelow(j, lower); });

    return rebuilt;
}

auto duplicatePhaseNameMatches(String const& phaseName, String const& baseName, String const& suffixSeparator) -> bool
{
    if(phaseName == baseName)
        return true;

    const auto prefix = baseName + suffixSeparator;
    return phaseName.size() > prefix.size() && phaseName.compare(0, prefix.size(), prefix) == 0;
}

auto findDefinitionForPhase(
    String const& phaseName,
    Vec<GlobalizedSolidSolutionPhaseDefinition> const& definitions) -> GlobalizedSolidSolutionPhaseDefinition const*
{
    for(auto const& definition : definitions)
        if(duplicatePhaseNameMatches(phaseName, definition.phaseName, definition.suffixSeparator))
            return &definition;
    return nullptr;
}

auto splitRequestKey(String const& phaseScope) -> String
{
    return "GlobalizedSolidSolution::Phase::" + phaseScope + "::SplitRequest";
}

auto findSplitRequestForDefinition(
    Map<String, Any> const& extra,
    GlobalizedSolidSolutionPhaseDefinition const& definition) -> Optional<GlobalizedSolidSolutionSplitRequest>
{
    const auto scoped = extra.find(splitRequestKey(definition.phaseName));
    if(scoped != extra.end())
        if(const auto* split = std::any_cast<GlobalizedSolidSolutionSplitRequest>(&scoped->second))
            return *split;

    const auto generic = extra.find("GlobalizedSolidSolution::SplitRequest");
    if(generic != extra.end())
        if(const auto* split = std::any_cast<GlobalizedSolidSolutionSplitRequest>(&generic->second))
            if(split->phaseName.empty() || split->phaseName == definition.phaseName)
                return *split;

    return {};
}

auto extractSplitCandidates(GlobalizedSolidSolutionSplitRequest const& splitRequest) -> Vec<SolidSolutionCandidateState>
{
    const auto it = splitRequest.extra.find("MAGEMinSolidSolutionPilot::SplitCandidates");
    if(it == splitRequest.extra.end())
        return {};
    if(const auto* candidates = std::any_cast<Vec<SolidSolutionCandidateState>>(&it->second))
        return *candidates;
    return {};
}

auto anyToReal(Any const& value) -> Optional<real>
{
    if(const auto* v = std::any_cast<real>(&value)) return *v;
    if(const auto* v = std::any_cast<double>(&value)) return static_cast<real>(*v);
    if(const auto* v = std::any_cast<float>(&value)) return static_cast<real>(*v);
    if(const auto* v = std::any_cast<int>(&value)) return static_cast<real>(*v);
    if(const auto* v = std::any_cast<long>(&value)) return static_cast<real>(*v);
    if(const auto* v = std::any_cast<long long>(&value)) return static_cast<real>(*v);
    if(const auto* v = std::any_cast<unsigned>(&value)) return static_cast<real>(*v);
    if(const auto* v = std::any_cast<unsigned long>(&value)) return static_cast<real>(*v);
    if(const auto* v = std::any_cast<unsigned long long>(&value)) return static_cast<real>(*v);
    return {};
}

auto splitRequestObjectiveGap(GlobalizedSolidSolutionSplitRequest const& splitRequest) -> Optional<real>
{
    const auto it = splitRequest.extra.find("MAGEMinSolidSolutionPilot::CompetingStableBranchObjectiveGap");
    if(it == splitRequest.extra.end())
        return {};
    return anyToReal(it->second);
}

auto shouldAcceptSplitRetry(
    EquilibriumResult const& result,
    Map<String, Any> const& extra,
    Vec<GlobalizedSolidSolutionPhaseDefinition> const& definitions,
    GlobalizedSolidSolutionEquilibriumRetryOptions const& options) -> bool
{
    if(result.iterations() < options.minIterationsForSplitRetry)
        return false;

    bool hasRequestedSplit = false;
    for(const auto& definition : definitions)
    {
        const auto splitRequest = findSplitRequestForDefinition(extra, definition);
        if(!splitRequest || !splitRequest->requested)
            continue;

        hasRequestedSplit = true;

        if(options.splitImprovementTolerance > 0.0)
        {
            const auto objectiveGap = splitRequestObjectiveGap(*splitRequest);
            if(objectiveGap && std::isfinite(static_cast<double>(*objectiveGap))
                && static_cast<double>(*objectiveGap) > static_cast<double>(options.splitImprovementTolerance))
                continue;
        }

        return true;
    }

    return hasRequestedSplit ? false : true;
}

auto duplicatePhaseBranchIndex(
    String const& phaseName,
    GlobalizedSolidSolutionPhaseDefinition const& definition) -> Index
{
    for(Index i = 0; i < definition.branches.size(); ++i)
    {
        const auto suffix = !definition.branches[i].label.empty()
            ? definition.branches[i].label
            : (!definition.branches[i].id.empty() ? definition.branches[i].id : std::to_string(i));
        if(phaseName == definition.phaseName + definition.suffixSeparator + suffix)
            return i;
    }

    return GlobalizedSolidSolutionNoBranch;
}

auto normalizedSplitWeights(ArrayXr weights) -> ArrayXr
{
    if(weights.size() == 0)
        return weights;

    for(Index i = 0; i < weights.size(); ++i)
        if(weights[i] < 0.0)
            weights[i] = 0.0;

    const auto sum = static_cast<double>(weights.sum());
    if(sum <= 0.0)
        return ArrayXr::Constant(weights.size(), 1.0/static_cast<double>(weights.size()));
    return weights/sum;
}

auto nudgeSeedInsideBranch(
    ArrayXr seed,
    GlobalizedSolidSolutionBranch const& branch,
    real margin = 1.0e-3) -> ArrayXr
{
    if(seed.size() != branch.lowerBounds.size() || seed.size() != branch.upperBounds.size())
        return seed;

    auto lower = ArrayXr(branch.lowerBounds);
    auto upper = ArrayXr(branch.upperBounds);
    for(Index i = 0; i < seed.size(); ++i)
    {
        const auto span = static_cast<double>(upper[i] - lower[i]);
        if(span <= 2.0*margin)
            continue;
        lower[i] += margin;
        upper[i] -= margin;
    }

    for(Index i = 0; i < seed.size(); ++i)
        seed[i] = std::clamp(static_cast<double>(seed[i]), static_cast<double>(lower[i]), static_cast<double>(upper[i]));

    for(Index iter = 0; iter < seed.size() * 8; ++iter)
    {
        const auto residual = 1.0 - static_cast<double>(seed.sum());
        if(std::abs(residual) <= 1.0e-12)
            break;

        real capacity = 0.0;
        for(Index i = 0; i < seed.size(); ++i)
        {
            const auto slack = residual > 0.0 ? upper[i] - seed[i] : seed[i] - lower[i];
            if(slack > 1.0e-12)
                capacity += slack;
        }

        if(capacity <= 1.0e-12)
            break;

        for(Index i = 0; i < seed.size(); ++i)
        {
            const auto slack = residual > 0.0 ? upper[i] - seed[i] : seed[i] - lower[i];
            if(slack <= 1.0e-12)
                continue;

            const auto delta = residual * (slack / capacity);
            seed[i] = std::clamp(
                static_cast<double>(seed[i] + delta),
                static_cast<double>(lower[i]),
                static_cast<double>(upper[i]));
        }
    }

    return seed/seed.sum();
}

auto seededMixingObjective(
    Vec<ArrayXr> const& seeds,
    ArrayXrConstRef bulkx,
    ArrayXrConstRef weights) -> real
{
    ArrayXr mixed = ArrayXr::Zero(bulkx.size());
    for(Index i = 0; i < weights.size(); ++i)
        mixed += weights[i]*seeds[i];
    return (mixed - bulkx).matrix().squaredNorm();
}

auto solveSplitWeights(
    Vec<ArrayXr> const& seeds,
    ArrayXrConstRef bulkx) -> ArrayXr
{
    if(seeds.empty())
        return {};

    if(seeds.size() == 1)
        return ArrayXr::Ones(1);


    if(seeds.size() == 2)
    {
        Index bestCoord = 0;
        auto bestDelta = 0.0;
        for(Index i = 0; i < bulkx.size(); ++i)
        {
            const auto delta = std::abs(static_cast<double>(seeds[0][i] - seeds[1][i]));
            if(delta > bestDelta)
            {
                bestDelta = delta;
                bestCoord = i;
            }
        }

        if(bestDelta > 1.0e-12)
        {
            const auto w0 = std::clamp(
                static_cast<double>((bulkx[bestCoord] - seeds[1][bestCoord])/(seeds[0][bestCoord] - seeds[1][bestCoord])),
                0.0,
                1.0);
            ArrayXr weights(2);
            weights << w0, 1.0 - w0;
            return weights;
        }
    }
    ArrayXr weights = ArrayXr::Constant(seeds.size(), 1.0/static_cast<double>(seeds.size()));
    real currentObjective = seededMixingObjective(seeds, bulkx, weights);

    for(Index iter = 0; iter < 64; ++iter)
    {
        ArrayXr mixed = ArrayXr::Zero(bulkx.size());
        for(Index i = 0; i < weights.size(); ++i)
            mixed += weights[i]*seeds[i];
        const auto residual = mixed - bulkx;

        ArrayXr gradient(weights.size());
        for(Index i = 0; i < weights.size(); ++i)
            gradient[i] = 2.0*residual.matrix().dot(seeds[i].matrix());
        auto projected = gradient.array() - gradient.mean();
        if(projected.cwiseAbs().maxCoeff() <= 1.0e-12)
            break;

        auto step = 0.5;
        auto accepted = false;
        for(Index backtrack = 0; backtrack < 16; ++backtrack)
        {
            ArrayXr trial = normalizedSplitWeights(weights - step*projected);
            const auto trialObjective = seededMixingObjective(seeds, bulkx, trial);
            if(trialObjective < currentObjective)
            {
                weights = trial;
                currentObjective = trialObjective;
                accepted = true;
                break;
            }
            step *= 0.5;
        }

        if(!accepted)
            break;
    }

    return weights;
}

auto allocateSeededFamilyAmounts(
    ArrayXrConstRef familyAmounts,
    Vec<ArrayXr> const& seeds) -> Vec<ArrayXr>
{
    if(seeds.empty())
        return {};

    const auto totalAmount = static_cast<double>(familyAmounts.sum());
    if(totalAmount <= 0.0)
        return {};

    const auto bulkx = familyAmounts/totalAmount;
    const auto weights = solveSplitWeights(seeds, bulkx);

    Vec<ArrayXr> allocations;
    allocations.reserve(seeds.size());
    for(Index i = 0; i < weights.size(); ++i)
        allocations.push_back(totalAmount*weights[i]*seeds[i]);
    return allocations;
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

auto formatIndices(Indices const& values) -> String
{
    std::ostringstream out;
    out << "[";
    for(Index i = 0; i < static_cast<Index>(values.size()); ++i)
    {
        if(i)
            out << ", ";
        out << values[i];
    }
    out << "]";
    return out.str();
}

auto rebuildStateForGlobalizedSolidSolutionSystem(
    ChemicalState const& source,
    ChemicalSystem const& rebuiltSystem,
    Vec<GlobalizedSolidSolutionPhaseDefinition> const& definitions) -> ChemicalState
{
    ChemicalState rebuilt(rebuiltSystem);
    rebuilt.setTemperature(source.temperature());
    rebuilt.setPressure(source.pressure());

    ArrayXr rebuiltAmounts = ArrayXr::Zero(rebuiltSystem.species().size());
    Map<String, ArrayXr> aggregatedFamilyAmounts;
    Map<String, Index> rebuiltFamilyCounts;
    Map<String, Indices> rebuiltFamilyPhaseIndices;
    Map<Index, ArrayXr> seededPhaseAmounts;

    for(Index iphase = 0; iphase < source.system().phases().size(); ++iphase)
    {
        const auto& phase = source.system().phase(iphase);
        const auto* definition = findDefinitionForPhase(phase.name(), definitions);
        if(!definition)
            continue;

        auto& amounts = aggregatedFamilyAmounts[definition->phaseName];
        if(amounts.size() == 0)
            amounts = ArrayXr::Zero(phase.species().size());
        amounts += source.speciesAmountsInPhase(iphase);
    }

    for(Index iphase = 0; iphase < rebuiltSystem.phases().size(); ++iphase)
    {
        const auto& phase = rebuiltSystem.phase(iphase);
        const auto* definition = findDefinitionForPhase(phase.name(), definitions);
        if(definition)
        {
            rebuiltFamilyCounts[definition->phaseName] += 1;
            rebuiltFamilyPhaseIndices[definition->phaseName].push_back(iphase);
        }
    }

    const auto sourceExtra = source.props().extra();
    for(const auto& [phaseName, phaseIndices] : rebuiltFamilyPhaseIndices)
    {
        if(phaseIndices.size() < 2)
            continue;

        const auto* definition = findDefinitionForPhase(rebuiltSystem.phase(phaseIndices[0]).name(), definitions);
        if(!definition)
            continue;

        const auto oldAmounts = aggregatedFamilyAmounts.find(phaseName);
        if(oldAmounts == aggregatedFamilyAmounts.end() || oldAmounts->second.size() == 0)
            continue;

        const auto splitRequest = findSplitRequestForDefinition(sourceExtra, *definition);
        if(!splitRequest || !splitRequest->requested)
            continue;

        std::cerr
            << "[EquilibriumSplitDiag] phaseFamily=" << phaseName
            << " triggerBranch=" << splitRequest->triggeringBranch
            << " requestedBranches=" << formatIndices(splitRequest->branches)
            << " reason=" << splitRequest->reason
            << "\n";

        const auto splitCandidates = extractSplitCandidates(*splitRequest);
        if(splitCandidates.size() < 2)
            continue;

        Vec<ArrayXr> seeds;
        Indices seededPhaseIndices;
        seeds.reserve(phaseIndices.size());
        seededPhaseIndices.reserve(phaseIndices.size());

        for(const auto iphase : phaseIndices)
        {
            const auto branchIndex = duplicatePhaseBranchIndex(rebuiltSystem.phase(iphase).name(), *definition);
            if(branchIndex == GlobalizedSolidSolutionNoBranch)
                continue;

            const auto& branch = definition->branches[branchIndex];

            for(const auto& candidate : splitCandidates)
            {
                if(candidate.branch != branchIndex || candidate.seedx.size() != oldAmounts->second.size())
                    continue;

                auto seed = candidate.seedx;
                for(Index i = 0; i < seed.size(); ++i)
                    if(seed[i] < 0.0)
                        seed[i] = 0.0;
                const auto sum = static_cast<double>(seed.sum());
                if(sum <= 0.0)
                    continue;
                seed /= sum;
                // Nudge seed towards center of composition space to avoid numerical instability at branch boundaries
                seed *= 0.99;
                seed /= seed.sum();
                seed = nudgeSeedInsideBranch(std::move(seed), branch);

                seeds.push_back(seed);
                seededPhaseIndices.push_back(iphase);
                std::cerr
                    << "[EquilibriumSplitDiag] phase=" << rebuiltSystem.phase(iphase).name()
                    << " requestedBranch=" << branchIndex
                    << " seed=" << formatArray(seed)
                    << "\n";
                break;
            }
        }

        if(seeds.size() < 2 || seededPhaseIndices.size() != seeds.size())
            continue;

        const auto allocations = allocateSeededFamilyAmounts(oldAmounts->second, seeds);
        if(allocations.size() != seededPhaseIndices.size())
            continue;

        for(Index i = 0; i < seededPhaseIndices.size(); ++i)
        {
            seededPhaseAmounts[seededPhaseIndices[i]] = allocations[i];
            std::cerr
                << "[EquilibriumSplitDiag] allocatedPhase=" << rebuiltSystem.phase(seededPhaseIndices[i]).name()
                << " requestedBranch=" << duplicatePhaseBranchIndex(rebuiltSystem.phase(seededPhaseIndices[i]).name(), *definition)
                << " allocatedSpeciesAmounts=" << formatArray(allocations[i])
                << "\n";
        }
    }

    for(Index iphase = 0; iphase < rebuiltSystem.phases().size(); ++iphase)
    {
        const auto& phase = rebuiltSystem.phase(iphase);
        const auto offset = rebuiltSystem.phases().numSpeciesUntilPhase(iphase);
        const auto count = phase.species().size();
        const auto* definition = findDefinitionForPhase(phase.name(), definitions);

        if(definition)
        {
            const auto seeded = seededPhaseAmounts.find(iphase);
            if(seeded != seededPhaseAmounts.end())
            {
                rebuiltAmounts.segment(offset, count) = seeded->second;
                continue;
            }

            const auto oldAmounts = aggregatedFamilyAmounts.find(definition->phaseName);
            const auto newCount = rebuiltFamilyCounts.find(definition->phaseName);
            if(oldAmounts != aggregatedFamilyAmounts.end() && newCount != rebuiltFamilyCounts.end() && newCount->second > 0)
                rebuiltAmounts.segment(offset, count) = oldAmounts->second/static_cast<real>(newCount->second);
            continue;
        }

        const auto oldPhaseIndex = source.system().phases().find(phase.name());
        if(oldPhaseIndex != source.system().phases().size())
            rebuiltAmounts.segment(offset, count) = source.speciesAmountsInPhase(oldPhaseIndex);
    }

    rebuilt.setSpeciesAmounts(rebuiltAmounts);
    return rebuilt;
}

auto canonicalizeDuplicatedPhaseBranchAssignments(
    ChemicalState& state,
    Vec<GlobalizedSolidSolutionPhaseDefinition> const& definitions) -> void
{
    auto amounts = ArrayXr(state.speciesAmounts());
    auto changed = false;

    for(const auto& definition : definitions)
    {
        Indices phaseIndices;
        for(Index iphase = 0; iphase < state.system().phases().size(); ++iphase)
        {
            if(duplicatePhaseNameMatches(state.system().phase(iphase).name(), definition.phaseName, definition.suffixSeparator))
                phaseIndices.push_back(iphase);
        }

        if(phaseIndices.size() != 2)
            continue;

        const auto leftIndex = phaseIndices[0];
        const auto rightIndex = phaseIndices[1];
        const auto leftBranch = duplicatePhaseBranchIndex(state.system().phase(leftIndex).name(), definition);
        const auto rightBranch = duplicatePhaseBranchIndex(state.system().phase(rightIndex).name(), definition);
        if(leftBranch == GlobalizedSolidSolutionNoBranch || rightBranch == GlobalizedSolidSolutionNoBranch)
            continue;

        const auto leftAmounts = ArrayXr(state.speciesAmountsInPhase(leftIndex));
        const auto rightAmounts = ArrayXr(state.speciesAmountsInPhase(rightIndex));
        const auto leftTotal = static_cast<double>(leftAmounts.sum());
        const auto rightTotal = static_cast<double>(rightAmounts.sum());
        if(leftTotal <= 0.0 || rightTotal <= 0.0)
            continue;

        const auto leftx = leftAmounts/leftTotal;
        const auto rightx = rightAmounts/rightTotal;
        const auto currentCost = static_cast<double>(
            GlobalizedSolidSolutionBranchViolation(leftx, definition.branches[leftBranch], 1.0e-12)
            + GlobalizedSolidSolutionBranchViolation(rightx, definition.branches[rightBranch], 1.0e-12));
        const auto swappedCost = static_cast<double>(
            GlobalizedSolidSolutionBranchViolation(rightx, definition.branches[leftBranch], 1.0e-12)
            + GlobalizedSolidSolutionBranchViolation(leftx, definition.branches[rightBranch], 1.0e-12));

        if(swappedCost + 1.0e-12 >= currentCost)
            continue;

        const auto leftOffset = state.system().phases().numSpeciesUntilPhase(leftIndex);
        const auto rightOffset = state.system().phases().numSpeciesUntilPhase(rightIndex);
        const auto count = state.system().phase(leftIndex).species().size();
        if(state.system().phase(rightIndex).species().size() != count)
            continue;

        amounts.segment(leftOffset, count) = rightAmounts;
        amounts.segment(rightOffset, count) = leftAmounts;
        changed = true;
    }

    if(changed)
        state.setSpeciesAmounts(amounts);
}

} // namespace

auto equilibrate(ChemicalState& state) -> EquilibriumResult
{
    EquilibriumOptions options;
    EquilibriumRestrictions restrictions(state.system());
    return equilibrate(state, restrictions, options);
}

auto equilibrate(ChemicalState& state, const EquilibriumOptions& options) -> EquilibriumResult
{
    EquilibriumRestrictions restrictions(state.system());
    return equilibrate(state, restrictions, options);
}

auto equilibrate(ChemicalState& state, const EquilibriumRestrictions& restrictions) -> EquilibriumResult
{
    EquilibriumOptions options;
    return equilibrate(state, restrictions, options);
}

auto equilibrate(ChemicalState& state, const EquilibriumRestrictions& restrictions, const EquilibriumOptions& options) -> EquilibriumResult
{
    const ArrayXd b0 = state.componentAmounts();
    return solveClosedSystemEquilibrium(state, restrictions, options, b0);
}

auto equilibrate(ChemicalState& state, ArrayXdConstRef b0) -> EquilibriumResult
{
    EquilibriumOptions options;
    EquilibriumRestrictions restrictions(state.system());
    return equilibrate(state, restrictions, options, b0);
}

auto equilibrate(ChemicalState& state, const EquilibriumOptions& options, ArrayXdConstRef b0) -> EquilibriumResult
{
    EquilibriumRestrictions restrictions(state.system());
    return equilibrate(state, restrictions, options, b0);
}

auto equilibrate(ChemicalState& state, const EquilibriumRestrictions& restrictions, ArrayXdConstRef b0) -> EquilibriumResult
{
    EquilibriumOptions options;
    return equilibrate(state, restrictions, options, b0);
}

auto equilibrate(ChemicalState& state, const EquilibriumRestrictions& restrictions, const EquilibriumOptions& options, ArrayXdConstRef b0) -> EquilibriumResult
{
    return solveClosedSystemEquilibrium(state, restrictions, options, b0);
}

auto equilibrateWithGlobalizedSolidSolutionSplits(
    ChemicalState const& initialState,
    GlobalizedSolidSolutionEquilibriumRetryOptions options) -> GlobalizedSolidSolutionEquilibriumRetryResult
{
    const ArrayXd b0 = initialState.componentAmounts();
    return equilibrateWithGlobalizedSolidSolutionSplits(initialState, b0, std::move(options));
}

auto equilibrateWithGlobalizedSolidSolutionSplits(
    ChemicalState const& initialState,
    ArrayXdConstRef b0,
    GlobalizedSolidSolutionEquilibriumRetryOptions options) -> GlobalizedSolidSolutionEquilibriumRetryResult
{
    EquilibriumSpecs specs(initialState.system());
    specs.temperature();
    specs.pressure();

    EquilibriumConditions conditions(specs);
    conditions.temperature(initialState.temperature());
    conditions.pressure(initialState.pressure());
    conditions.setInitialComponentAmounts(b0);

    EquilibriumRestrictions restrictions(initialState.system());
    return equilibrateWithGlobalizedSolidSolutionSplits(initialState, specs, conditions, restrictions, std::move(options));
}

auto equilibrateWithGlobalizedSolidSolutionSplits(
    ChemicalState const& initialState,
    EquilibriumSpecs const& specs,
    EquilibriumConditions const& conditions,
    GlobalizedSolidSolutionEquilibriumRetryOptions options) -> GlobalizedSolidSolutionEquilibriumRetryResult
{
    EquilibriumRestrictions restrictions(initialState.system());
    return equilibrateWithGlobalizedSolidSolutionSplits(initialState, specs, conditions, restrictions, std::move(options));
}

auto equilibrateWithGlobalizedSolidSolutionSplits(
    ChemicalState const& initialState,
    EquilibriumSpecs const& specs,
    EquilibriumConditions const& conditions,
    EquilibriumRestrictions const& restrictions,
    GlobalizedSolidSolutionEquilibriumRetryOptions options) -> GlobalizedSolidSolutionEquilibriumRetryResult
{
    ChemicalState state(initialState);
    EquilibriumResult totalResult;
    Index numRebuilds = 0;
    Index numAcceptedSplitRetries = 0;
    Index numRejectedSplitRetries = 0;
    Index numRebuiltFallbackAttempts = 0;

    const auto solveRebuiltSystemWithFallback = [&](ChemicalState& rebuiltState) -> EquilibriumResult
    {
        const auto rebuiltSpecs = rebuildEquilibriumSpecsForSystem(specs, rebuiltState.system());
        const auto rebuiltConditions = rebuildEquilibriumConditionsForSystem(conditions, rebuiltSpecs);
        const auto rebuiltRestrictions = rebuildEquilibriumRestrictionsForSystem(restrictions, rebuiltState.system());

        auto stageOptions = options.equilibrium;
        auto rebuiltResult = solveEquilibrium(rebuiltState, rebuiltSpecs, rebuiltConditions, rebuiltRestrictions, stageOptions);

        if(rebuiltResult.succeeded() || !options.enableRebuiltFallbackSolve || options.rebuiltFallbackMaxRetries == 0)
            return rebuiltResult;

        auto accumulated = rebuiltResult;
        for(Index attempt = 0; attempt < options.rebuiltFallbackMaxRetries; ++attempt)
        {
            ++numRebuiltFallbackAttempts;
            auto fallbackOptions = options.equilibrium;
            if(attempt == 0 && options.rebuiltFallbackUseIdealFirst)
                fallbackOptions.use_ideal_activity_models = true;

            auto fallbackResult = solveEquilibrium(rebuiltState, rebuiltSpecs, rebuiltConditions, rebuiltRestrictions, fallbackOptions);
            accumulated += fallbackResult;
            if(fallbackResult.succeeded())
                break;
        }

        return accumulated;
    };

    for(;;)
    {
        const auto solvingSystem = state.system();
        const auto solvingOriginalSystem = state.system().id() == initialState.system().id();
        const auto result = solvingOriginalSystem
            ? solveEquilibrium(state, specs, conditions, restrictions, options.equilibrium)
            : solveRebuiltSystemWithFallback(state);

        if(!solvingOriginalSystem && result.succeeded())
            canonicalizeDuplicatedPhaseBranchAssignments(state, options.definitions);

        totalResult += result;

        if(result.failed() || options.definitions.empty() || numRebuilds >= options.maxRetries)
            return {state.system(), state, totalResult, numRebuilds, numAcceptedSplitRetries, numRejectedSplitRetries, numRebuiltFallbackAttempts};

        if(solvingOriginalSystem && options.enableSplitAcceptanceGate)
        {
            const auto accepted = shouldAcceptSplitRetry(result, state.props().extra(), options.definitions, options);
            if(!accepted)
            {
                ++numRejectedSplitRetries;
                return {state.system(), state, totalResult, numRebuilds, numAcceptedSplitRetries, numRejectedSplitRetries, numRebuiltFallbackAttempts};
            }
            ++numAcceptedSplitRetries;
        }

        const auto rebuiltSystem = ApplyGlobalizedSolidSolutionSplitRequests(
            solvingSystem,
            state.props().extra(),
            options.definitions);

        if(rebuiltSystem.id() == state.system().id())
            return {state.system(), state, totalResult, numRebuilds, numAcceptedSplitRetries, numRejectedSplitRetries, numRebuiltFallbackAttempts};

        state = rebuildStateForGlobalizedSolidSolutionSystem(state, rebuiltSystem, options.definitions);
        ++numRebuilds;
    }
}

} // namespace Reaktoro
