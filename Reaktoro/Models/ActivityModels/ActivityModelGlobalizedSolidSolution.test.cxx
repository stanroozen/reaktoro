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

// Catch includes
#include <catch2/catch.hpp>

// C++ includes
#include <chrono>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Equilibrium/EquilibriumUtils.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelGlobalizedBinaryRedlichKister.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIdealSolution.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelMAGEMinSolidSolutionPilot.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelGlobalizedSolidSolution.hpp>
#include <Reaktoro/Models/ActivityModels/Support/MixedSystemSpeciesAmountUtils.hpp>
using namespace Reaktoro;

namespace test {
extern auto createDatabase() -> Database;
}

TEST_CASE("Testing ActivityModelGlobalizedSolidSolution", "[ActivityModelGlobalizedSolidSolution]")
{
    const auto species = SpeciesList("NaCl KCl");

    ArrayXr x(2);
    x << 0.25, 0.75;

    const auto T = 973.15;
    const auto P = 1.2e9;

    WHEN("The reduced solid-solution model returns a consistent result")
    {
        auto reduced = [](GlobalizedSolidSolutionInput input)
        {
            GlobalizedSolidSolutionOutput output;
            output.Gx = 125.0;
            output.Hx = 250.0;
            output.Cpx = 12.5;
            output.Vx = 3.5e-5;
            output.VxT = -2.0e-8;
            output.VxP = 4.0e-14;
            output.Vxi = input.x.square();
            output.ln_g = 2.0 * input.x;
            output.ln_a = output.ln_g + input.x.log();
            output.som = StateOfMatter::Solid;
            output.selectedBranch = 1;
            output.branch.id = "ordered";
            output.branch.label = "ordered";
            output.branch.lowerBounds = ArrayXr::Zero(2);
            output.branch.upperBounds = ArrayXr::Ones(2);
            output.branches = {output.branch};
            input.state->data["probe"] = String("reused");
            output.state = input.state;
            output.extra["branch"] = String("ordered");
            output.extra["saw-seed"] = input.extra.find("seed") != input.extra.end();
            return output;
        };

        ActivityModel fn = ActivityModelGlobalizedSolidSolution(reduced)(species);
        ActivityProps props = ActivityProps::create(species.size());
        props.extra["seed"] = String("carry");

        fn(props, {T, P, x});

        CHECK(props.Gx == Approx(125.0));
        CHECK(props.Hx == Approx(250.0));
        CHECK(props.Cpx == Approx(12.5));
        CHECK(props.Vx == Approx(3.5e-5));
        CHECK(props.VxT == Approx(-2.0e-8));
        CHECK(props.VxP == Approx(4.0e-14));
        CHECK(props.Vxi[0] == Approx(0.0625));
        CHECK(props.Vxi[1] == Approx(0.5625));
        CHECK(props.ln_g[0] == Approx(0.5));
        CHECK(props.ln_g[1] == Approx(1.5));
        CHECK(props.ln_a[0] == Approx(0.5 + log(0.25)));
        CHECK(props.ln_a[1] == Approx(1.5 + log(0.75)));
        CHECK(props.som == StateOfMatter::Solid);
        CHECK(std::any_cast<String>(props.extra.at("branch")) == "ordered");
        CHECK(std::any_cast<String>(props.extra.at("seed")) == "carry");
        CHECK(std::any_cast<bool>(props.extra.at("saw-seed")));
        CHECK(std::any_cast<std::uint64_t>(props.extra.at("GlobalizedSolidSolution::SelectedBranchIndex")) == 1);
        CHECK(std::any_cast<String>(props.extra.at("GlobalizedSolidSolution::SelectedBranchId")) == "ordered");

        const auto state = std::any_cast<SharedPtr<GlobalizedSolidSolutionState>>(props.extra.at("GlobalizedSolidSolution::State"));
        REQUIRE(state);
        CHECK(state->selectedBranch == 1);
        CHECK(std::any_cast<String>(state->data.at("probe")) == "reused");
        CHECK(!state->lastSplitRequest.requested);
    }

    WHEN("The reduced solid-solution model returns vectors with the wrong size")
    {
        auto reduced = [](GlobalizedSolidSolutionInput input)
        {
            GlobalizedSolidSolutionOutput output;
            output.Vxi = input.x;
            output.ln_g = input.x.head(1);
            output.ln_a = input.x.log();
            return output;
        };

        ActivityModel fn = ActivityModelGlobalizedSolidSolution(reduced)(species);
        ActivityProps props = ActivityProps::create(species.size());

        REQUIRE_THROWS_WITH(
            fn(props, {T, P, x}),
            "ActivityModelGlobalizedSolidSolution received ln_g with invalid size.");
    }
}

TEST_CASE("Testing ActivityModelGlobalizedBinaryRedlichKister", "[ActivityModelGlobalizedBinaryRedlichKister]")
{
    const auto species = SpeciesList("NaCl KCl");
    const auto T = 973.15;
    const auto P = 1.2e9;

    GlobalizedBinaryRedlichKisterOptions options;
    options.a0 = 1.5;
    options.a1 = -0.2;
    options.a2 = 0.05;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "solvus-left";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "solvus-right";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    options.branches = {left, right};

    ActivityModel fn = ActivityModelGlobalizedBinaryRedlichKister(options)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr xleft(2);
    xleft << 0.20, 0.80;
    fn(props, {T, P, xleft});

    CHECK(std::any_cast<String>(props.extra.at("GlobalizedSolidSolution::SelectedBranchId")) == "left");
    CHECK(std::any_cast<String>(props.extra.at("GlobalizedBinaryRedlichKister::BranchLabel")) == "solvus-left");
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedBinaryRedlichKister::UsedWarmstart")));
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedBinaryRedlichKister::ReusedStateCache")));

    const auto state1 = std::any_cast<SharedPtr<GlobalizedSolidSolutionState>>(props.extra.at("GlobalizedSolidSolution::State"));
    REQUIRE(state1);
    CHECK(state1->selectedBranch == 0);
    CHECK(state1->numEvaluations == 1);

    ArrayXr xleft2(2);
    xleft2 << 0.25, 0.75;
    fn(props, {T, P, xleft2});

    const auto state2 = std::any_cast<SharedPtr<GlobalizedSolidSolutionState>>(props.extra.at("GlobalizedSolidSolution::State"));
    REQUIRE(state2);
    CHECK(state2 == state1);
    CHECK(std::any_cast<bool>(props.extra.at("GlobalizedBinaryRedlichKister::UsedWarmstart")));
    CHECK(std::any_cast<bool>(props.extra.at("GlobalizedBinaryRedlichKister::ReusedStateCache")));
    CHECK(state2->numEvaluations == 2);

    ArrayXr xright(2);
    xright << 0.80, 0.20;
    fn(props, {T, P, xright});

    CHECK(std::any_cast<String>(props.extra.at("GlobalizedSolidSolution::SelectedBranchId")) == "right");
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedBinaryRedlichKister::UsedWarmstart")));

    const auto state3 = std::any_cast<SharedPtr<GlobalizedSolidSolutionState>>(props.extra.at("GlobalizedSolidSolution::State"));
    REQUIRE(state3);
    CHECK(state3->selectedBranch == 1);
    CHECK(state3->numEvaluations == 3);

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("GlobalizedBinaryRedlichKister::InternalComposition"));
    CHECK(internalx.size() == 2);
    CHECK(internalx.sum() == Approx(1.0));
}

TEST_CASE("Testing ActivityModelGlobalizedBinaryRedlichKister custom candidate screening", "[ActivityModelGlobalizedSolidSolution][CandidateGenerator]")
{
    const auto species = SpeciesList("NaCl KCl");
    const auto T = 973.15;
    const auto P = 1.2e9;

    GlobalizedBinaryRedlichKisterOptions options;
    options.a0 = 1.5;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "solvus-left";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "solvus-right";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    options.branches = {left, right};

    auto sawSeed = std::make_shared<bool>(false);
    options.candidateGenerator = [=](GlobalizedSolidSolutionInput const& input, Vec<GlobalizedSolidSolutionBranch> const&) mutable
    {
        *sawSeed = input.extra.find("seed") != input.extra.end();

        GlobalizedSolidSolutionCandidate candidate;
        candidate.branch = 1;
        candidate.extra["GlobalizedBinaryRedlichKister::CandidateSource"] = String("custom-screen");

        return Vec<GlobalizedSolidSolutionCandidate>{candidate};
    };

    ActivityModel fn = ActivityModelGlobalizedBinaryRedlichKister(options)(species);
    ActivityProps props = ActivityProps::create(species.size());
    props.extra["seed"] = String("carry");

    ArrayXr x(2);
    x << 0.20, 0.80;
    fn(props, {T, P, x});

    CHECK(*sawSeed);
    CHECK(std::any_cast<String>(props.extra.at("GlobalizedSolidSolution::SelectedBranchId")) == "right");
    CHECK(std::any_cast<String>(props.extra.at("GlobalizedBinaryRedlichKister::CandidateSource")) == "custom-screen");
}

TEST_CASE("Testing ActivityModelGlobalizedBinaryRedlichKister stability screening", "[ActivityModelGlobalizedSolidSolution][Stability]")
{
    const auto species = SpeciesList("NaCl KCl");
    const auto T = 973.15;
    const auto P = 1.2e9;

    GlobalizedBinaryRedlichKisterOptions options;
    options.a0 = 1.5;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "solvus-left";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "solvus-right";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    options.branches = {left, right};
    options.stabilityCriterion = [](GlobalizedSolidSolutionInput const&, GlobalizedSolidSolutionBranch const& branch, ArrayXrConstRef, real)
    {
        GlobalizedSolidSolutionCandidateStability stability;
        if(branch.id == "left")
        {
            stability.stable = false;
            stability.reason = "forced-instability";
            stability.extra["GlobalizedBinaryRedlichKister::StabilityProbe"] = String("rejected-left");
            return stability;
        }

        stability.penalty = 7.0;
        stability.reason = "fallback-branch";
        stability.extra["GlobalizedBinaryRedlichKister::StabilityProbe"] = String("accepted-right");
        return stability;
    };

    ActivityModel fn = ActivityModelGlobalizedBinaryRedlichKister(options)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(2);
    x << 0.20, 0.80;
    fn(props, {T, P, x});

    CHECK(std::any_cast<String>(props.extra.at("GlobalizedSolidSolution::SelectedBranchId")) == "right");
    CHECK(std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::CandidateStable")));
    CHECK(std::any_cast<real>(props.extra.at("GlobalizedSolidSolution::CandidateStabilityPenalty")) == Approx(7.0));
    CHECK(std::any_cast<String>(props.extra.at("GlobalizedSolidSolution::CandidateStabilityReason")) == "fallback-branch");
    CHECK(std::any_cast<String>(props.extra.at("GlobalizedBinaryRedlichKister::StabilityProbe")) == "accepted-right");
}

TEST_CASE("Testing ActivityModelGlobalizedBinaryRedlichKister split trigger policy", "[ActivityModelGlobalizedSolidSolution][SplitTrigger]")
{
    const auto species = SpeciesList("NaCl KCl");
    const auto T = 973.15;
    const auto P = 1.2e9;

    GlobalizedBinaryRedlichKisterOptions options;
    options.a0 = 1.5;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "solvus-left";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "solvus-right";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    options.branches = {left, right};
    options.stabilityCriterion = [=](GlobalizedSolidSolutionInput const&, GlobalizedSolidSolutionBranch const&, ArrayXrConstRef, real score)
    {
        GlobalizedSolidSolutionCandidateStability stability;
        if(score <= 0.0)
            return stability;

        stability.stable = false;
        stability.reason = "stability-screen-between-branches";
        stability.splitRequest.requested = true;
        stability.splitRequest.branches = {0, 1};
        stability.splitRequest.branchIds = {"left", "right"};
        stability.splitRequest.reason = stability.reason;
        return stability;
    };

    ActivityModel fn = ActivityModelGlobalizedBinaryRedlichKister(options)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(2);
    x << 0.50, 0.50;
    fn(props, {T, P, x});

    CHECK(std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
    const auto splitRequest = std::any_cast<GlobalizedSolidSolutionSplitRequest>(props.extra.at("GlobalizedSolidSolution::SplitRequest"));
    CHECK(splitRequest.requested);
    CHECK(splitRequest.branches.size() == 2);
    CHECK(splitRequest.branchIds[0] == "left");
    CHECK(splitRequest.branchIds[1] == "right");
    CHECK(splitRequest.reason == "stability-screen-between-branches");
}

TEST_CASE("Testing shared globalized solid-solution default candidates", "[ActivityModelGlobalizedSolidSolution][CandidateGenerator][Default]")
{
    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    const auto branches = Vec<GlobalizedSolidSolutionBranch>{left, right};

    const real T = 973.15;
    const real P = 1.2e9;

    auto state = std::make_shared<GlobalizedSolidSolutionState>();
    state->cachedBranchForState = 1;
    state->chemicalPropsStateId = 17;
    state->cachedInternalx = (ArrayXr(2) << 0.80, 0.20).finished();
    state->selectedBranch = 1;
    state->lastInternalx = (ArrayXr(2) << 0.82, 0.18).finished();

    ArrayXr x(2);
    x << 0.20, 0.80;

    Map<String, Any> extra;
    extra["Reaktoro::ChemicalProps::StateId"] = static_cast<std::uint64_t>(17);

    GlobalizedSolidSolutionInput input{T, P, x, extra, state, GlobalizedSolidSolutionNoBranch};

    GlobalizedSolidSolutionDefaultCandidateOptions options;
    options.branchTolerance = 1.0e-8;
    options.cachedStatePriority = -1.0;
    options.preferredBranchPriority = -1.0e-6;
    options.requireCachedStateWarmstart = true;
    options.sourceKey = "CandidateSource";

    const auto cachedCandidates = DefaultGlobalizedSolidSolutionCandidates(
        input,
        branches,
        x,
        Optional<ArrayXr>(state->cachedInternalx),
        Optional<ArrayXr>(state->lastInternalx),
        options);

    REQUIRE(cachedCandidates.size() == 1);
    CHECK(cachedCandidates[0].branch == 1);
    CHECK(cachedCandidates[0].priority == Approx(-1.0));
    CHECK(cachedCandidates[0].initialInternalx.size() == 2);
    CHECK(cachedCandidates[0].initialInternalx[0] == Approx(0.80));
    CHECK(std::any_cast<String>(cachedCandidates[0].extra.at("CandidateSource")) == "state-cache");

    state->chemicalPropsStateId = 0;
    Map<String, Any> requestedExtra;
    GlobalizedSolidSolutionInput requestedInput{T, P, x, requestedExtra, state, 0};

    const auto requestedCandidates = DefaultGlobalizedSolidSolutionCandidates(
        requestedInput,
        branches,
        x,
        Optional<ArrayXr>(state->cachedInternalx),
        Optional<ArrayXr>(state->lastInternalx),
        options);

    REQUIRE(requestedCandidates.size() == 1);
    CHECK(requestedCandidates[0].branch == 0);
    CHECK(requestedCandidates[0].initialInternalx[0] == Approx(0.82));
    CHECK(std::any_cast<String>(requestedCandidates[0].extra.at("CandidateSource")) == "requested-branch");
}

TEST_CASE("Testing MAGEMin pilot split reassembly helper", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot]")
{
    const auto species = SpeciesList("NaCl KCl");
    const auto T = 973.15;
    const auto P = 1.2e9;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "left";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "right";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    Phase prototype;
    prototype = prototype.withName("MAGEMinPilot");
    prototype = prototype.withSpecies(species);
    prototype = prototype.withStateOfMatter(StateOfMatter::Solid);

    auto reduced = [=](GlobalizedSolidSolutionInput input)
    {
        GlobalizedSolidSolutionOutput output;
        output.Vxi = ArrayXr::Zero(2);
        output.ln_g = ArrayXr::Zero(2);
        output.ln_a = input.x.log();
        output.branch = left;
        output.branches = {left, right};
        output.splitRequest.requested = true;
        output.splitRequest.branches = {0, 1};
        output.splitRequest.branchIds = {"left", "right"};
        output.splitRequest.reason = "pilot-split";
        return output;
    };

    MAGEMinSolidSolutionPilotOptions pilotOptions;
    pilotOptions.branches = {left, right};
    const auto definition = MAGEMinSolidSolutionPilotDefinition(prototype, reduced, pilotOptions);

    Map<String, Any> splitExtra;
    GlobalizedSolidSolutionSplitRequest splitRequest;
    splitRequest.requested = true;
    splitRequest.phaseName = "MAGEMinPilot";
    splitRequest.branches = {0, 1};
    splitRequest.branchIds = {"left", "right"};
    splitRequest.reason = "pilot-split";
    splitExtra["GlobalizedSolidSolution::Phase::MAGEMinPilot::SplitRequest"] = splitRequest;

    const auto rebuilt = ApplyGlobalizedSolidSolutionSplitRequests(PhaseList{definition.prototype}, splitExtra, {definition});
    REQUIRE(rebuilt.size() == 2);
    CHECK(rebuilt[0].name() == "MAGEMinPilot#left");
    CHECK(rebuilt[1].name() == "MAGEMinPilot#right");

    Map<String, Any> collapsedExtra;
    GlobalizedSolidSolutionSplitRequest noSplit;
    collapsedExtra["GlobalizedSolidSolution::Phase::MAGEMinPilot::SplitRequest"] = noSplit;
    const auto collapsed = ApplyGlobalizedSolidSolutionSplitRequests(rebuilt, collapsedExtra, {definition});
    REQUIRE(collapsed.size() == 1);
    CHECK(collapsed[0].name() == "MAGEMinPilot");
}

TEST_CASE("Testing imported MAGEMin SB11 olivine pilot model", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB11]")
{
    const auto species = SpeciesList("Mg2SiO4 Fe2SiO4");
    const auto T = 1473.15;
    const auto P = 1.0e9;
    const auto W = 7813.22;

    MAGEMinSB11OlivineOptions options;
    const auto model = MAGEMinSolidSolutionPilotModelSB11Olivine(options);

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(2);
    x << 0.25, 0.75;
    fn(props, {T, P, x});

    const auto RT = universalGasConstant * T;
    const auto expectedGx = 0.25 * 0.75 * W;
    const auto expectedLnGfo = (0.75 * 0.75 * W)/RT;
    const auto expectedLnGfa = (0.25 * 0.25 * W)/RT;

    CHECK(props.Gx == Approx(expectedGx));
    CHECK(props.Hx == Approx(expectedGx));
    CHECK(props.ln_g[0] == Approx(expectedLnGfo));
    CHECK(props.ln_g[1] == Approx(expectedLnGfa));
    CHECK(props.ln_a[0] == Approx(expectedLnGfo + 2.0*log(0.25)));
    CHECK(props.ln_a[1] == Approx(expectedLnGfa + 2.0*log(0.75)));
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb11_ol");
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::Endmember0")) == "fo");
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::Endmember1")) == "fa");
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));

    fn(props, {T, P, x});
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::UsedStateCache")));
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::CandidateSource")) == "state-cache");
}

TEST_CASE("Testing imported MAGEMin SB11 wadsleyite pilot model", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB11][Wadsleyite]")
{
    const auto species = SpeciesList("Mg2SiO4 Fe2SiO4");
    const auto T = 1473.15;
    const auto P = 1.0e9;
    const auto W = 16747.18;

    const auto model = MAGEMinSolidSolutionPilotModelSB11Wadsleyite();

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(2);
    x << 0.60, 0.40;
    fn(props, {T, P, x});

    const auto RT = universalGasConstant * T;
    const auto expectedGx = 0.60 * 0.40 * W;
    const auto expectedLnG0 = (0.40 * 0.40 * W)/RT;
    const auto expectedLnG1 = (0.60 * 0.60 * W)/RT;

    CHECK(props.Gx == Approx(expectedGx));
    CHECK(props.Hx == Approx(expectedGx));
    CHECK(props.ln_g[0] == Approx(expectedLnG0));
    CHECK(props.ln_g[1] == Approx(expectedLnG1));
    CHECK(props.ln_a[0] == Approx(expectedLnG0 + 2.0*log(0.60)));
    CHECK(props.ln_a[1] == Approx(expectedLnG1 + 2.0*log(0.40)));
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb11_wa");
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::Endmember0")) == "mgwa");
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::Endmember1")) == "fewa");
}

TEST_CASE("Testing imported MAGEMin SB11 akimotoite pilot model", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB11][Akimotoite]")
{
    const auto species = SpeciesList("MgSiO3 FeSiO3 CaSiO3");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    const auto model = MAGEMinSolidSolutionPilotModelSB11Akimotoite();

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(3);
    x << 0.80, 0.10, 0.10;
    fn(props, {T, P, x});

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    const auto endmembers = std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::Endmembers"));

    CHECK(internalx.size() == 3);
    CHECK(internalx.sum() == Approx(1.0));
    CHECK(internalx[0] == Approx(0.80).margin(5.0e-2));
    CHECK(props.Gx < props.Hx);
    CHECK(props.ln_g.size() == 3);
    CHECK(props.ln_a.size() == 3);
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb11_ak");
    CHECK(endmembers == Strings{"mgak", "feak", "co"});
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerConverged")));
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
}

TEST_CASE("Testing imported MAGEMin SB11 perovskite pilot model", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB11][Perovskite]")
{
    const auto species = SpeciesList("MgSiO3 FeSiO3 Al2O3");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    const auto model = MAGEMinSolidSolutionPilotModelSB11Perovskite();

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(3);
    x << 0.70, 0.20, 0.10;
    fn(props, {T, P, x});

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    const auto endmembers = std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::Endmembers"));

    CHECK(internalx.size() == 3);
    CHECK(internalx.sum() == Approx(1.0));
    CHECK(std::isfinite(static_cast<double>(props.Gx)));
    CHECK(std::isfinite(static_cast<double>(props.Hx)));
    CHECK(props.Gx < props.Hx);
    CHECK(std::isfinite(static_cast<double>(props.ln_g[0])));
    CHECK(std::isfinite(static_cast<double>(props.ln_g[1])));
    CHECK(std::isfinite(static_cast<double>(props.ln_g[2])));
    CHECK(std::isfinite(static_cast<double>(props.ln_a[0])));
    CHECK(std::isfinite(static_cast<double>(props.ln_a[1])));
    CHECK(std::isfinite(static_cast<double>(props.ln_a[2])));
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb11_pv");
    CHECK(endmembers == Strings{"mgpv", "fepv", "alpv"});
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerConverged")));
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
    CHECK(std::any_cast<std::uint64_t>(props.extra.at("MAGEMinSolidSolutionPilot::GeneratedCandidateCount")) >= 4);

    ArrayXr x2(3);
    x2 << 0.68, 0.22, 0.10;
    fn(props, {T, P, x2});

    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::UsedWarmstart")));
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::UsedStateCache")));
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::CandidateSource")) == "state-cache");
}

TEST_CASE("Testing imported MAGEMin ternary pilot multi-seed screening diagnostics", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB11][Perovskite][CandidateGenerator]")
{
    const auto species = SpeciesList("MgSiO3 FeSiO3 Al2O3");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "mgpv-rich";
    left.lowerBounds = ArrayXr::Zero(3);
    left.upperBounds = ArrayXr::Ones(3);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "mgpv-poor";
    right.lowerBounds = ArrayXr::Zero(3);
    right.upperBounds = ArrayXr::Ones(3);
    right.lowerBounds[0] = 0.55;

    MAGEMinSB11PerovskiteOptions options;
    options.branchPolicy.branches = {left, right};

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Perovskite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(3);
    x << 0.70, 0.20, 0.10;
    fn(props, {T, P, x});

    CHECK(std::any_cast<std::uint64_t>(props.extra.at("MAGEMinSolidSolutionPilot::GeneratedCandidateCount")) >= 8);
    CHECK(std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::GeneratedCandidateSeedLabels"))
        == Strings{"visible-composition", "dominant::alpv", "dominant::mgpv", "dominant::fepv", "edge::alpv-mgpv", "edge::alpv-fepv", "edge::mgpv-fepv"});
    CHECK(props.extra.count("MAGEMinSolidSolutionPilot::CandidateSeedLabel") == 1);
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::CandidateSource")) == "branch-screen");
}

TEST_CASE("Testing imported MAGEMin akimotoite ternary seed priorities", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB11][Akimotoite][CandidateGenerator]")
{
    const auto species = SpeciesList("MgSiO3 FeSiO3 CaSiO3");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "mgak-rich";
    left.lowerBounds = ArrayXr::Zero(3);
    left.upperBounds = ArrayXr::Ones(3);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "mgak-poor";
    right.lowerBounds = ArrayXr::Zero(3);
    right.upperBounds = ArrayXr::Ones(3);
    right.lowerBounds[0] = 0.55;

    MAGEMinSB11AkimotoiteOptions options;
    options.branchPolicy.branches = {left, right};

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Akimotoite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(3);
    x << 0.70, 0.20, 0.10;
    fn(props, {T, P, x});

    CHECK(std::any_cast<std::uint64_t>(props.extra.at("MAGEMinSolidSolutionPilot::GeneratedCandidateCount")) >= 8);
    CHECK(std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::GeneratedCandidateSeedLabels"))
        == Strings{"visible-composition", "dominant::co", "dominant::mgak", "dominant::feak", "edge::co-mgak", "edge::co-feak", "edge::mgak-feak"});
}

TEST_CASE("Testing imported MAGEMin calcioferrite ternary seed priorities", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB11][Calcioferrite][CandidateGenerator]")
{
    const auto species = SpeciesList("MgAl2O4 FeAl2O4 NaAlSiO4");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "mgcf-rich";
    left.lowerBounds = ArrayXr::Zero(3);
    left.upperBounds = ArrayXr::Ones(3);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "mgcf-poor";
    right.lowerBounds = ArrayXr::Zero(3);
    right.upperBounds = ArrayXr::Ones(3);
    right.lowerBounds[0] = 0.55;

    MAGEMinSB11CalcioferriteOptions options;
    options.branchPolicy.branches = {left, right};

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Calcioferrite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(3);
    x << 0.60, 0.25, 0.15;
    fn(props, {T, P, x});

    CHECK(std::any_cast<std::uint64_t>(props.extra.at("MAGEMinSolidSolutionPilot::GeneratedCandidateCount")) >= 8);
    CHECK(std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::GeneratedCandidateSeedLabels"))
        == Strings{"visible-composition", "dominant::nacf", "dominant::mgcf", "dominant::fecf", "edge::nacf-mgcf", "edge::nacf-fecf", "edge::mgcf-fecf"});
}

TEST_CASE("Testing imported MAGEMin SB21 spinel pilot model", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB21][Spinel]")
{
    const auto species = SpeciesList("MgAl2O4 FeAl2O4");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    const auto model = MAGEMinSolidSolutionPilotModelSB21Spinel();

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(2);
    x << 0.65, 0.35;
    fn(props, {T, P, x});

    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb21_sp");
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::Endmember0")) == "sp");
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::Endmember1")) == "hc");
    CHECK(std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"))[0] == Approx(0.65));
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
}

TEST_CASE("Benchmarking imported MAGEMin SB21 spinel pilot model", "[.][ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB21][Spinel][Benchmark]")
{
    const auto species = SpeciesList("MgAl2O4 FeAl2O4");
    const auto model = MAGEMinSolidSolutionPilotModelSB21Spinel();
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(2);
    x << 0.65, 0.35;

    const auto start = std::chrono::steady_clock::now();
    for(Index i = 0; i < 2000; ++i)
        fn(props, {1473.15, 1.0e9, x});
    const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();

    INFO("SB21 spinel benchmark elapsed_us=" << elapsed);
    CHECK(elapsed > 0);
}

TEST_CASE("Testing imported MAGEMin SB21 NAL pilot model", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB21][NAL]")
{
    const auto species = SpeciesList("NaMg2Al5SiO12 NaFe2Al5SiO12 Na3Al3Si3O12");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    const auto model = MAGEMinSolidSolutionPilotModelSB21NAL();

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(3);
    x << 0.55, 0.20, 0.25;
    fn(props, {T, P, x});

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    const auto endmembers = std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::Endmembers"));

    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb21_nal");
    CHECK(endmembers == Strings{"mnal", "fnal", "nnal"});
    CHECK(internalx.size() == 3);
    CHECK(internalx.sum() == Approx(1.0));
    CHECK(std::any_cast<std::uint64_t>(props.extra.at("MAGEMinSolidSolutionPilot::GeneratedCandidateCount")) >= 8);
    CHECK(std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::GeneratedCandidateSeedLabels"))
        == Strings{"visible-composition", "dominant::nnal", "dominant::mnal", "dominant::fnal", "edge::nnal-mnal", "edge::nnal-fnal", "edge::mnal-fnal"});
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
}

TEST_CASE("Benchmarking imported MAGEMin SB21 NAL pilot model", "[.][ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB21][NAL][Benchmark]")
{
    const auto species = SpeciesList("NaMg2Al5SiO12 NaFe2Al5SiO12 Na3Al3Si3O12");
    const auto model = MAGEMinSolidSolutionPilotModelSB21NAL();
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(3);
    x << 0.55, 0.20, 0.25;

    const auto start = std::chrono::steady_clock::now();
    for(Index i = 0; i < 2000; ++i)
        fn(props, {1473.15, 1.0e9, x});
    const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();

    INFO("SB21 NAL benchmark elapsed_us=" << elapsed);
    CHECK(elapsed > 0);
}

TEST_CASE("Testing imported MAGEMin SB21 calcioferrite pilot model", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB21][Calcioferrite]")
{
    const auto species = SpeciesList("MgAl2O4 FeAl2O4 NaAlSiO4");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    const auto model = MAGEMinSolidSolutionPilotModelSB21Calcioferrite();

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(3);
    x << 0.60, 0.25, 0.15;
    fn(props, {T, P, x});

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    const auto endmembers = std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::Endmembers"));

    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb21_cf");
    CHECK(endmembers == Strings{"mgcf", "fecf", "nacf"});
    CHECK(internalx.size() == 3);
    CHECK(internalx.sum() == Approx(1.0));
    CHECK(std::any_cast<std::uint64_t>(props.extra.at("MAGEMinSolidSolutionPilot::GeneratedCandidateCount")) >= 8);
    CHECK(std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::GeneratedCandidateSeedLabels"))
        == Strings{"visible-composition", "dominant::nacf", "dominant::mgcf", "dominant::fecf", "edge::nacf-mgcf", "edge::nacf-fecf", "edge::mgcf-fecf"});
    CHECK(props.ln_g[0] != Approx(0.0));
    CHECK(props.ln_g[1] != Approx(0.0));
    CHECK(props.ln_g[2] != Approx(0.0));
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
}

TEST_CASE("Benchmarking imported MAGEMin SB21 calcioferrite pilot model", "[.][ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB21][Calcioferrite][Benchmark]")
{
    const auto species = SpeciesList("MgAl2O4 FeAl2O4 NaAlSiO4");
    const auto model = MAGEMinSolidSolutionPilotModelSB21Calcioferrite();
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(3);
    x << 0.60, 0.25, 0.15;

    const auto start = std::chrono::steady_clock::now();
    for(Index i = 0; i < 2000; ++i)
        fn(props, {1473.15, 1.0e9, x});
    const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();

    INFO("SB21 calcioferrite benchmark elapsed_us=" << elapsed);
    CHECK(elapsed > 0);
}

TEST_CASE("Benchmarking imported MAGEMin ternary pilot families comparatively", "[.][ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Benchmark][Comparative]")
{
    struct BenchmarkCase
    {
        String label;
        ActivityModel fn;
        ActivityProps props;
        real T = 1473.15;
        real P = 1.0e9;
        ArrayXr x;
    };

    Vec<BenchmarkCase> cases;

    {
        const auto species = SpeciesList("MgSiO3 FeSiO3 CaSiO3");
        cases.push_back({"sb11_ak", ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Akimotoite())(species), ActivityProps::create(species.size()), 1473.15, 1.0e9, (ArrayXr(3) << 0.80, 0.10, 0.10).finished()});
    }
    {
        const auto species = SpeciesList("MgSiO3 FeSiO3 Al2O3");
        cases.push_back({"sb11_pv", ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Perovskite())(species), ActivityProps::create(species.size()), 1473.15, 1.0e9, (ArrayXr(3) << 0.70, 0.20, 0.10).finished()});
    }
    {
        const auto species = SpeciesList("MgAl2O4 FeAl2O4 NaAlSiO4");
        cases.push_back({"sb11_cf", ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Calcioferrite())(species), ActivityProps::create(species.size()), 1473.15, 1.0e9, (ArrayXr(3) << 0.60, 0.25, 0.15).finished()});
    }
    {
        const auto species = SpeciesList("MgAl2O4 FeAl2O4");
        cases.push_back({"sb21_sp", ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Spinel())(species), ActivityProps::create(species.size()), 1473.15, 1.0e9, (ArrayXr(2) << 0.65, 0.35).finished()});
    }
    {
        const auto species = SpeciesList("NaMg2Al5SiO12 NaFe2Al5SiO12 Na3Al3Si3O12");
        cases.push_back({"sb21_nal", ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21NAL())(species), ActivityProps::create(species.size()), 1473.15, 1.0e9, (ArrayXr(3) << 0.55, 0.20, 0.25).finished()});
    }
    {
        const auto species = SpeciesList("MgAl2O4 FeAl2O4 NaAlSiO4");
        cases.push_back({"sb21_cf", ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Calcioferrite())(species), ActivityProps::create(species.size()), 1473.15, 1.0e9, (ArrayXr(3) << 0.60, 0.25, 0.15).finished()});
    }

    for(auto& benchmarkCase : cases)
    {
        const auto start = std::chrono::steady_clock::now();
        for(Index i = 0; i < 2000; ++i)
            benchmarkCase.fn(benchmarkCase.props, {benchmarkCase.T, benchmarkCase.P, benchmarkCase.x});
        const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();

        INFO(benchmarkCase.label << " comparative benchmark elapsed_us=" << elapsed);
        CHECK(elapsed > 0);
    }
}

TEST_CASE("Testing imported MAGEMin SB21 OPX pilot model", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB21][OPX]")
{
    // 4-endmember orthopyroxene: en, fs, mgts, odi
    const auto species = SpeciesList("MgSiO3 FeSiO3 MgAl2SiO6 CaMgSi2O6");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    const auto model = MAGEMinSolidSolutionPilotModelSB21OPX();

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(4);
    x << 0.50, 0.20, 0.10, 0.20;
    fn(props, {T, P, x});

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    const auto endmembers = std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::Endmembers"));

    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb21_opx");
    CHECK(endmembers == Strings{"en", "fs", "mgts", "odi"});
    CHECK(internalx.size() == 4);
    CHECK(internalx.sum() == Approx(1.0));
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerConverged")));
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
    // At x={0.5, 0.2, 0.1, 0.2} the odi interaction (W_03=W_13=32217 J/mol) produces non-zero excess
    CHECK(props.ln_g[3] != Approx(0.0));
}

TEST_CASE("Testing imported MAGEMin SB21 CPX pilot model", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB21][CPX]")
{
    // 5-endmember clinopyroxene: di, he, cen, cats, jd
    const auto species = SpeciesList("CaMgSi2O6 CaFeSi2O6 Mg2Si2O6 CaAl2SiO6 NaAlSi2O6");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    const auto model = MAGEMinSolidSolutionPilotModelSB21CPX();

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(5);
    x << 0.45, 0.15, 0.15, 0.10, 0.15;
    fn(props, {T, P, x});

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    const auto endmembers = std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::Endmembers"));

    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb21_cpx");
    CHECK(endmembers == Strings{"di", "he", "cen", "cats", "jd"});
    CHECK(internalx.size() == 5);
    CHECK(internalx.sum() == Approx(1.0));
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerConverged")));
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
    // Volume-fraction Margules with v[cats]=3.5 produces non-trivial excess chemical potentials
    CHECK(props.ln_g[3] != Approx(0.0));
}

TEST_CASE("Testing imported MAGEMin SB21 garnet-majorite pilot model", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB21][GTMJ]")
{
    // 5-endmember garnet-majorite: py, alm, gr, mgmj, jdmj
    const auto species = SpeciesList("Mg3Al2Si3O12 Fe3Al2Si3O12 Ca3Al2Si3O12 Mg4Si4O12 Na2Al2Si4O12");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    const auto model = MAGEMinSolidSolutionPilotModelSB21GTMJ();

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(5);
    x << 0.50, 0.20, 0.15, 0.10, 0.05;
    fn(props, {T, P, x});

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    const auto endmembers = std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::Endmembers"));

    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb21_gtmj");
    CHECK(endmembers == Strings{"py", "alm", "gr", "mgmj", "jdmj"});
    CHECK(internalx.size() == 5);
    CHECK(internalx.sum() == Approx(1.0));
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerConverged")));
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
    // W_34(mgmj-jdmj)=70879 J/mol produces non-zero excess
    CHECK(props.ln_g[3] != Approx(0.0));
}

TEST_CASE("Testing tangent-plane stability check — stable single-phase sb21_cf composition", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB21][Calcioferrite][TPD]")
{
    // A composition near one branch endpoint is thermodynamically stable as a single phase.
    // The TPD criterion should report stable=true and not request a split.
    const auto species = SpeciesList("MgAl2O4 FeAl2O4 NaAlSiO4");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    MAGEMinSB21CalcioferriteOptions options;
    options.enableTangentPlaneStabilityCheck = true;
    options.tpdTolerance = 1.0e-4;

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Calcioferrite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    // y ≈ (0.80, 0.15, 0.05) — nacf-poor, well inside the mgcf-fecf branch
    ArrayXr x(3);
    x << 0.80, 0.15, 0.05;
    fn(props, {T, P, x});

    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb21_cf");
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::TPDStable")));
}

TEST_CASE("Testing tangent-plane stability check — unstable two-phase sb21_cf composition", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB21][Calcioferrite][TPD]")
{
    // The sb21_cf calcioferrite family has W02 = W12 = 60825 J/mol, T_crit >> T = 1473 K.
    // A bulk composition deep in the two-phase region (nacf ≈ 50 %) should be identified
    // as unstable by the TPD criterion, triggering a split request.
    const auto species = SpeciesList("MgAl2O4 FeAl2O4 NaAlSiO4");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    MAGEMinSB21CalcioferriteOptions options;
    options.enableTangentPlaneStabilityCheck = true;
    options.tpdTolerance = 1.0e-4;

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Calcioferrite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    // y ≈ (0.25, 0.25, 0.50) — nacf ≈ 50 %, deep inside the solvus
    ArrayXr x(3);
    x << 0.25, 0.25, 0.50;
    fn(props, {T, P, x});

    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb21_cf");
    const auto splitRequested = std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested"));
    const auto tpdStable = std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::TPDStable"));
    const auto tpdMin = std::any_cast<real>(props.extra.at("MAGEMinSolidSolutionPilot::TPDMinValue"));
    CHECK(static_cast<double>(tpdMin) <= 0.0);
    CHECK(splitRequested == !tpdStable);
    CHECK(tpdStable == (static_cast<double>(tpdMin) >= -options.tpdTolerance));
}

TEST_CASE("Testing imported MAGEMin SB11 calcioferrite pilot model", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB11][Calcioferrite]")
{
    const auto species = SpeciesList("MgAl2O4 FeAl2O4 NaAlSiO4");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    const auto model = MAGEMinSolidSolutionPilotModelSB11Calcioferrite();

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(3);
    x << 0.60, 0.25, 0.15;
    fn(props, {T, P, x});

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    const auto endmembers = std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::Endmembers"));

    CHECK(internalx.size() == 3);
    CHECK(internalx.sum() == Approx(1.0));
    CHECK(props.ln_g[0] == Approx(0.0));
    CHECK(props.ln_g[1] == Approx(0.0));
    CHECK(props.ln_g[2] == Approx(0.0));
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == "sb11_cf");
    CHECK(endmembers == Strings{"mgcf", "fecf", "nacf"});
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerConverged")));
    CHECK(!std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
}

TEST_CASE("Testing imported MAGEMin SB11 olivine split trigger", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB11][Split]")
{
    const auto species = SpeciesList("Mg2SiO4 Fe2SiO4");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "fo-rich";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "fa-rich";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    MAGEMinSB11OlivineOptions options;
    options.branchPolicy.branches = {left, right};

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Olivine(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(2);
    x << 0.50, 0.50;
    fn(props, {T, P, x});

    CHECK(std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
    const auto splitRequest = std::any_cast<GlobalizedSolidSolutionSplitRequest>(props.extra.at("GlobalizedSolidSolution::SplitRequest"));
    CHECK(splitRequest.requested);
    CHECK(splitRequest.branchIds.size() == 2);
    CHECK(splitRequest.branchIds[0] == "left");
    CHECK(splitRequest.branchIds[1] == "right");
    CHECK(splitRequest.reason == "stability-screen-between-branches");
}

TEST_CASE("Testing imported MAGEMin SB11 calcioferrite split trigger", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][SB11][Calcioferrite][Split]")
{
    const auto species = SpeciesList("MgAl2O4 FeAl2O4 NaAlSiO4");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "mgcf-rich";
    left.lowerBounds = ArrayXr::Zero(3);
    left.upperBounds = ArrayXr::Ones(3);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "mgcf-poor";
    right.lowerBounds = ArrayXr::Zero(3);
    right.upperBounds = ArrayXr::Ones(3);
    right.lowerBounds[0] = 0.55;

    MAGEMinSB11CalcioferriteOptions options;
    options.branchPolicy.branches = {left, right};

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Calcioferrite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(3);
    x << 0.50, 0.25, 0.25;
    fn(props, {T, P, x});

    CHECK(std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")));
    const auto splitRequest = std::any_cast<GlobalizedSolidSolutionSplitRequest>(props.extra.at("GlobalizedSolidSolution::SplitRequest"));
    CHECK(splitRequest.requested);
    CHECK(splitRequest.branchIds.size() == 2);
    CHECK(splitRequest.branchIds[0] == "left");
    CHECK(splitRequest.branchIds[1] == "right");
    CHECK(splitRequest.reason == "stability-screen-between-branches");
}

TEST_CASE("Testing imported MAGEMin pilot custom candidate screening", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][CandidateGenerator]")
{
    const auto species = SpeciesList("Mg2SiO4 Fe2SiO4");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "fo-rich";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "fa-rich";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    auto sawSeed = std::make_shared<bool>(false);

    MAGEMinSB11OlivineOptions options;
    options.branchPolicy.branches = {left, right};
    options.branchPolicy.candidateGenerator = [=](GlobalizedSolidSolutionInput const& input, Vec<GlobalizedSolidSolutionBranch> const&) mutable
    {
        *sawSeed = input.extra.find("seed") != input.extra.end();

        GlobalizedSolidSolutionCandidate candidate;
        candidate.branch = 1;
        candidate.extra["MAGEMinSolidSolutionPilot::CandidateSource"] = String("custom-screen");
        return Vec<GlobalizedSolidSolutionCandidate>{candidate};
    };
    options.branchPolicy.stabilityCriterion = [](GlobalizedSolidSolutionInput const&, GlobalizedSolidSolutionBranch const&, ArrayXrConstRef, real)
    {
        return GlobalizedSolidSolutionCandidateStability{};
    };

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Olivine(options))(species);
    ActivityProps props = ActivityProps::create(species.size());
    props.extra["seed"] = String("carry");

    ArrayXr x(2);
    x << 0.20, 0.80;
    fn(props, {T, P, x});

    CHECK(*sawSeed);
    CHECK(std::any_cast<String>(props.extra.at("GlobalizedSolidSolution::SelectedBranchId")) == "right");
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::CandidateSource")) == "custom-screen");
}

TEST_CASE("Testing imported MAGEMin pilot stability screening", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Stability]")
{
    const auto species = SpeciesList("Mg2SiO4 Fe2SiO4");
    const auto T = 1473.15;
    const auto P = 1.0e9;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "fo-rich";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "fa-rich";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    MAGEMinSB11OlivineOptions options;
    options.branchPolicy.branches = {left, right};
    options.branchPolicy.stabilityCriterion = [](GlobalizedSolidSolutionInput const&, GlobalizedSolidSolutionBranch const& branch, ArrayXrConstRef, real)
    {
        GlobalizedSolidSolutionCandidateStability stability;
        stability.stable = branch.id != "right";
        stability.reason = stability.stable ? String{} : String("screened-out-right");
        return stability;
    };

    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Olivine(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr x(2);
    x << 0.80, 0.20;
    fn(props, {T, P, x});

    CHECK(std::any_cast<String>(props.extra.at("GlobalizedSolidSolution::SelectedBranchId")) == "left");
}

TEST_CASE("Testing split-triggered equilibrium retry helper with MAGEMin pilot", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][EquilibriumRetry]")
{
    const auto db = test::createDatabase();
    const auto species = SpeciesList({db.species().get("CaCO3(s)"), db.species().get("MgCO3(s)")});

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "left";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "right";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    MAGEMinSB11OlivineOptions modelOptions;
    modelOptions.branchPolicy.branches = {left, right};
    const auto model = MAGEMinSolidSolutionPilotModelSB11Olivine(modelOptions);

    Phase prototype;
    prototype = prototype.withName("PilotCarbonate");
    prototype = prototype.withSpecies(species);
    prototype = prototype.withStateOfMatter(StateOfMatter::Solid);
    prototype = prototype.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(prototype.species()));

    MAGEMinSolidSolutionPilotOptions pilotOptions;
    pilotOptions.branches = {left, right};

    const auto phase = MAGEMinSolidSolutionPilotPhase(prototype, model);
    const auto definition = MAGEMinSolidSolutionPilotDefinition(prototype, model, pilotOptions);
    const auto system = ChemicalSystem(db, PhaseList{phase});

    ChemicalState state(system);
    state.setTemperature(1473.15);
    state.setPressure(1.0e9);
    ArrayXr n(2);
    n << 1.0, 1.0;
    state.setSpeciesAmounts(n);

    GlobalizedSolidSolutionEquilibriumRetryOptions retryOptions;
    retryOptions.definitions = {definition};
    retryOptions.maxRetries = 1;

    const auto result = equilibrateWithGlobalizedSolidSolutionSplits(state, retryOptions);

    CHECK(result.result.succeeded());
    CHECK(result.numRebuilds == 1);
    REQUIRE(result.system.phases().size() == 2);
    CHECK(result.system.phase(0).name() == "PilotCarbonate#left");
    CHECK(result.system.phase(1).name() == "PilotCarbonate#right");
}

TEST_CASE("Testing split-triggered equilibrium retry helper preserves richer conditions and restrictions", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][EquilibriumRetry][Conditions]")
{
    const auto db = test::createDatabase();
    const auto aqueous = AqueousPhase({"H2O(aq)", "H+(aq)", "OH-(aq)", "Na+(aq)", "Cl-(aq)"});
    const auto carbonateSpecies = SpeciesList({db.species().get("CaCO3(s)"), db.species().get("MgCO3(s)")});
    const auto quartzSpecies = SpeciesList({db.species().get("SiO2(s)")});

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "left";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "right";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    MAGEMinSB11OlivineOptions modelOptions;
    modelOptions.branchPolicy.branches = {left, right};
    const auto model = MAGEMinSolidSolutionPilotModelSB11Olivine(modelOptions);

    Phase prototype;
    prototype = prototype.withName("PilotCarbonate");
    prototype = prototype.withSpecies(carbonateSpecies);
    prototype = prototype.withStateOfMatter(StateOfMatter::Solid);
    prototype = prototype.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(prototype.species()));

    Phase quartz;
    quartz = quartz.withName("Quartz");
    quartz = quartz.withSpecies(quartzSpecies);
    quartz = quartz.withStateOfMatter(StateOfMatter::Solid);
    quartz = quartz.withActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(quartz.species()));
    quartz = quartz.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(quartz.species()));

    MAGEMinSolidSolutionPilotOptions pilotOptions;
    pilotOptions.branches = {left, right};

    const auto phase = MAGEMinSolidSolutionPilotPhase(prototype, model);
    const auto definition = MAGEMinSolidSolutionPilotDefinition(prototype, model, pilotOptions);
    Phases phasesObject(db, aqueous);
    phasesObject.add(PhaseList{phase, quartz});
    const auto system = ChemicalSystem(phasesObject);

    ChemicalState state(system);
    state.setTemperature(1473.15);
    state.setPressure(1.0e9);
    ArrayXr n(8);
    n << 55.5, 1.0e-7, 1.0e-7, 1.0, 1.0, 1.0, 1.0, 2.0;
    state.setSpeciesAmounts(TestUtils::reorderPilotMixedConditionsSpeciesAmounts(
        system,
        carbonateSpecies,
        n,
        "Unexpected mixed-system species amount vector size in globalized solid-solution test."));

    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.pH();

    EquilibriumConditions conditions(specs);
    conditions.temperature(1473.15);
    conditions.pressure(1.0e9);
    conditions.pH(7.0);
    conditions.setInitialComponentAmountsFromState(state);

    EquilibriumRestrictions restrictions(system);
    restrictions.cannotReact("Na+(aq)");

    GlobalizedSolidSolutionEquilibriumRetryOptions retryOptions;
    retryOptions.definitions = {definition};
    retryOptions.maxRetries = 1;

    const auto result = equilibrateWithGlobalizedSolidSolutionSplits(state, specs, conditions, restrictions, retryOptions);

    CHECK(result.result.succeeded());
    CHECK(result.numRebuilds == 1);

    const auto sodiumIndex = result.system.species().indexWithName("Na+(aq)");
    CHECK(result.state.speciesAmounts()[sodiumIndex] == Approx(1.0));
}

TEST_CASE("Testing split-triggered equilibrium retry helper with seeded ternary MAGEMin pilot", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][EquilibriumRetry][Ternary]")
{
    const auto db = test::createDatabase();
    const auto species = SpeciesList({db.species().get("CaCO3(s)"), db.species().get("MgCO3(s)"), db.species().get("SiO2(s)")});

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "left";
    left.lowerBounds = ArrayXr::Zero(3);
    left.upperBounds = ArrayXr::Ones(3);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "right";
    right.lowerBounds = ArrayXr::Zero(3);
    right.upperBounds = ArrayXr::Ones(3);
    right.lowerBounds[0] = 0.55;

    MAGEMinSB11PerovskiteOptions modelOptions;
    modelOptions.branchPolicy.branches = {left, right};
    modelOptions.branchPolicy.stabilityCriterion = NamedGlobalizedSolidSolutionStabilityCriterion(
        NamedGlobalizedSolidSolutionStabilityPolicy::MAGEMinPilotBranchAmbiguity,
        modelOptions.branchPolicy.branches);
    const auto model = MAGEMinSolidSolutionPilotModelSB11Perovskite(modelOptions);

    Phase prototype;
    prototype = prototype.withName("PilotCarbonate");
    prototype = prototype.withSpecies(species);
    prototype = prototype.withStateOfMatter(StateOfMatter::Solid);
    prototype = prototype.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(prototype.species()));

    MAGEMinSolidSolutionPilotOptions pilotOptions;
    pilotOptions.branches = {left, right};

    const auto phase = MAGEMinSolidSolutionPilotPhase(prototype, model);
    const auto definition = MAGEMinSolidSolutionPilotDefinition(prototype, model, pilotOptions);
    const auto system = ChemicalSystem(db, PhaseList{phase});

    ChemicalState state(system);
    state.setTemperature(1473.15);
    state.setPressure(1.0e9);
    ArrayXr n(3);
    n << 2.0, 1.0, 1.0;
    state.setSpeciesAmounts(n);

    GlobalizedSolidSolutionEquilibriumRetryOptions retryOptions;
    retryOptions.definitions = {definition};
    retryOptions.maxRetries = 1;

    const auto result = equilibrateWithGlobalizedSolidSolutionSplits(state, retryOptions);

    CHECK(result.result.succeeded());
    CHECK(result.numRebuilds == 1);
    REQUIRE(result.system.phases().size() == 2);
    CHECK(result.system.phase(0).name() == "PilotCarbonate#left");
    CHECK(result.system.phase(1).name() == "PilotCarbonate#right");
}

TEST_CASE("Testing manual duplicated-phase outer workflow with candidate assembly", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][EquilibriumRetry][OuterWorkflow]")
{
    const auto db = test::createDatabase();
    const auto species = SpeciesList({db.species().get("CaCO3(s)"), db.species().get("MgCO3(s)")});

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "left";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "right";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    MAGEMinSB11OlivineOptions modelOptions;
    modelOptions.branchPolicy.branches = {left, right};
    const auto model = MAGEMinSolidSolutionPilotModelSB11Olivine(modelOptions);

    Phase prototype;
    prototype = prototype.withName("PilotCarbonate");
    prototype = prototype.withSpecies(species);
    prototype = prototype.withStateOfMatter(StateOfMatter::Solid);
    prototype = prototype.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(prototype.species()));

    MAGEMinSolidSolutionPilotOptions pilotOptions;
    pilotOptions.branches = {left, right};

    const auto phase = MAGEMinSolidSolutionPilotPhase(prototype, model);
    const auto definition = MAGEMinSolidSolutionPilotDefinition(prototype, model, pilotOptions);
    const auto initialSystem = ChemicalSystem(db, PhaseList{phase});

    ChemicalState initialState(initialSystem);
    initialState.setTemperature(1473.15);
    initialState.setPressure(1.0e9);
    ArrayXr n(2);
    n << 1.0, 1.0;
    initialState.setSpeciesAmounts(n);

    GlobalizedSolidSolutionEquilibriumRetryOptions retryOptions;
    retryOptions.definitions = {definition};
    retryOptions.maxRetries = 0;

    const auto firstPass = equilibrateWithGlobalizedSolidSolutionSplits(initialState, retryOptions);

    CHECK(firstPass.result.succeeded());
    CHECK(firstPass.numRebuilds == 0);
    REQUIRE(firstPass.state.props().extra().count("GlobalizedSolidSolution::SplitRequested") == 1);
    CHECK(std::any_cast<bool>(firstPass.state.props().extra().at("GlobalizedSolidSolution::SplitRequested")));

    Vec<SolidSolutionCandidateState> candidates;
    candidates.push_back({0, ArrayXr(), 0.0, "left"});
    candidates.push_back({1, ArrayXr(), 1.0, "right"});

    const auto duplicated = AssembleGlobalizedSolidSolutionCandidatePhases(
        prototype,
        model,
        pilotOptions.branches,
        candidates);

    REQUIRE(duplicated.size() == 2);
    CHECK(duplicated[0].name() == "PilotCarbonate#left");
    CHECK(duplicated[1].name() == "PilotCarbonate#right");

    const auto expandedSystem = ChemicalSystem(db, duplicated);
    ChemicalState expandedState(expandedSystem);
    expandedState.setTemperature(initialState.temperature());
    expandedState.setPressure(initialState.pressure());
    expandedState.setSpeciesAmounts(ArrayXr::Constant(expandedSystem.species().size(), 1.0));

    const auto secondPass = equilibrate(expandedState);

    CHECK(secondPass.succeeded());
    REQUIRE(expandedState.system().phases().size() == 2);
    CHECK(expandedState.system().phase(0).name() == "PilotCarbonate#left");
    CHECK(expandedState.system().phase(1).name() == "PilotCarbonate#right");
}

TEST_CASE("Testing split-triggered equilibrium retry helper preserves richer conditions for seeded ternary MAGEMin pilot", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][EquilibriumRetry][Ternary][Conditions]")
{
    const auto db = test::createDatabase();
    const auto aqueous = AqueousPhase({"H2O(aq)", "H+(aq)", "OH-(aq)", "Na+(aq)", "Cl-(aq)"});
    const auto carbonateSpecies = SpeciesList({db.species().get("CaCO3(s)"), db.species().get("MgCO3(s)"), db.species().get("CaMg(CO3)2(s)")});
    const auto quartzSpecies = SpeciesList({db.species().get("SiO2(s)")});

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "left";
    left.lowerBounds = ArrayXr::Zero(3);
    left.upperBounds = ArrayXr::Ones(3);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "right";
    right.lowerBounds = ArrayXr::Zero(3);
    right.upperBounds = ArrayXr::Ones(3);
    right.lowerBounds[0] = 0.55;

    MAGEMinSB11PerovskiteOptions modelOptions;
    modelOptions.branchPolicy.branches = {left, right};
    modelOptions.branchPolicy.stabilityCriterion = NamedGlobalizedSolidSolutionStabilityCriterion(
        NamedGlobalizedSolidSolutionStabilityPolicy::MAGEMinPilotBranchAmbiguity,
        modelOptions.branchPolicy.branches);
    const auto model = MAGEMinSolidSolutionPilotModelSB11Perovskite(modelOptions);

    Phase prototype;
    prototype = prototype.withName("PilotCarbonate");
    prototype = prototype.withSpecies(carbonateSpecies);
    prototype = prototype.withStateOfMatter(StateOfMatter::Solid);
    prototype = prototype.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(prototype.species()));

    Phase quartz;
    quartz = quartz.withName("Quartz");
    quartz = quartz.withSpecies(quartzSpecies);
    quartz = quartz.withStateOfMatter(StateOfMatter::Solid);
    quartz = quartz.withActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(quartz.species()));
    quartz = quartz.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(quartz.species()));

    MAGEMinSolidSolutionPilotOptions pilotOptions;
    pilotOptions.branches = {left, right};

    const auto phase = MAGEMinSolidSolutionPilotPhase(prototype, model);
    const auto definition = MAGEMinSolidSolutionPilotDefinition(prototype, model, pilotOptions);
    Phases phasesObject(db, aqueous);
    phasesObject.add(PhaseList{phase, quartz});
    const auto system = ChemicalSystem(phasesObject);

    ChemicalState state(system);
    state.setTemperature(1473.15);
    state.setPressure(1.0e9);
    ArrayXr n(9);
    n << 55.5, 1.0e-7, 1.0e-7, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0;
    state.setSpeciesAmounts(TestUtils::reorderPilotMixedConditionsSpeciesAmounts(
        system,
        carbonateSpecies,
        n,
        "Unexpected mixed-system species amount vector size in globalized solid-solution test."));

    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.pH();

    EquilibriumConditions conditions(specs);
    conditions.temperature(1473.15);
    conditions.pressure(1.0e9);
    conditions.pH(7.0);
    conditions.setInitialComponentAmountsFromState(state);

    EquilibriumRestrictions restrictions(system);
    restrictions.cannotReact("Na+(aq)");

    GlobalizedSolidSolutionEquilibriumRetryOptions retryOptions;
    retryOptions.definitions = {definition};
    retryOptions.maxRetries = 1;

    const auto result = equilibrateWithGlobalizedSolidSolutionSplits(state, specs, conditions, restrictions, retryOptions);

    CHECK(result.result.succeeded());
    CHECK(result.numRebuilds == 1);
    REQUIRE(result.system.phases().size() == 4);
    const auto sodiumIndex = result.system.species().indexWithName("Na+(aq)");
    CHECK(result.state.speciesAmounts()[sodiumIndex] == Approx(1.0));
}

TEST_CASE("Testing generic globalized solid-solution internal minimizer", "[ActivityModelGlobalizedSolidSolution][InternalMinimizer]")
{
    GlobalizedSolidSolutionInternalProblem problem;
    problem.objective = [](ArrayXrConstRef x)
    {
        return pow(static_cast<double>(x[0] - 0.20), 2)
            + pow(static_cast<double>(x[1] - 0.30), 2)
            + pow(static_cast<double>(x[2] - 0.50), 2);
    };

    problem.initialx = ArrayXr::Constant(3, 1.0 / 3.0);
    problem.lowerBounds = ArrayXr::Zero(3);
    problem.upperBounds = ArrayXr::Ones(3);
    problem.initialStep = 0.2;

    const auto result = MinimizeGlobalizedSolidSolutionInternalProblem(problem);

    CHECK(result.x.size() == 3);
    CHECK(result.x.sum() == Approx(1.0));
    CHECK(result.x[0] == Approx(0.20).margin(1.0e-4));
    CHECK(result.x[1] == Approx(0.30).margin(1.0e-4));
    CHECK(result.x[2] == Approx(0.50).margin(1.0e-4));
}

TEST_CASE("Testing duplicated globalized solid-solution branches", "[ActivityModelGlobalizedSolidSolution][PhaseDuplication]")
{
    const auto species = SpeciesList("NaCl KCl");
    const auto T = 973.15;
    const auto P = 1.2e9;

    GlobalizedBinaryRedlichKisterOptions options;
    options.a0 = 2.0;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "solvus-left";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "solvus-right";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    options.branches = {left, right};

    Phase phase;
    phase = phase.withName("AlkaliHalide");
    phase = phase.withSpecies(species);
    phase = phase.withStateOfMatter(StateOfMatter::Solid);
    phase = phase.withActivityModel(ActivityModelGlobalizedBinaryRedlichKister(options)(species));

    const auto duplicated = DuplicateGlobalizedBinaryRedlichKisterPhaseBranches(phase, options);

    REQUIRE(duplicated.size() == 2);
    CHECK(duplicated[0].name() == "AlkaliHalide#solvus-left");
    CHECK(duplicated[1].name() == "AlkaliHalide#solvus-right");

    ActivityProps propsLeft = ActivityProps::create(species.size());
    ArrayXr xleft(2);
    xleft << 0.20, 0.80;
    duplicated[0].activityModel()(propsLeft, {T, P, xleft});
    CHECK(std::any_cast<String>(propsLeft.extra.at("GlobalizedSolidSolution::SelectedBranchId")) == "left");

    ActivityProps propsRight = ActivityProps::create(species.size());
    ArrayXr xright(2);
    xright << 0.80, 0.20;
    duplicated[1].activityModel()(propsRight, {T, P, xright});
    CHECK(std::any_cast<String>(propsRight.extra.at("GlobalizedSolidSolution::SelectedBranchId")) == "right");
}

TEST_CASE("Testing GlobalizedBinaryRedlichKisterSolidPhases with ChemicalSystem", "[ActivityModelGlobalizedSolidSolution][ChemicalSystem]")
{
    const auto db = test::createDatabase();

    GlobalizedBinaryRedlichKisterOptions options;
    options.a0 = 2.0;

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "solvus-left";
    left.lowerBounds = ArrayXr::Zero(2);
    left.upperBounds = ArrayXr::Ones(2);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "solvus-right";
    right.lowerBounds = ArrayXr::Zero(2);
    right.upperBounds = ArrayXr::Ones(2);
    right.lowerBounds[0] = 0.55;

    options.branches = {left, right};

    const auto phases = GlobalizedBinarySolidPhases(db, "CarbonateSS", {"CaCO3(s)", "MgCO3(s)"}, options);

    REQUIRE(phases.size() == 2);
    CHECK(phases[0].name() == "CarbonateSS#solvus-left");
    CHECK(phases[1].name() == "CarbonateSS#solvus-right");

    const auto system = ChemicalSystem(db, phases);
    CHECK(system.phases().size() == 2);
    CHECK(system.phase(0).name() == "CarbonateSS#solvus-left");
    CHECK(system.phase(1).name() == "CarbonateSS#solvus-right");

    const auto phasesObject = Phases(db, AqueousPhase({"H2O(aq)", "H+(aq)", "OH-(aq)"}), GlobalizedBinarySolidPhases(db, "CarbonateSS", {"CaCO3(s)", "MgCO3(s)"}, options));
    const auto system2 = ChemicalSystem(phasesObject);
    CHECK(system2.phases().size() == 3);
}

TEST_CASE("Testing physically grounded immiscibility split trigger for SB21 calcioferrite pilot (nacf solvus)", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][EquilibriumRetry][Exsolution]")
{
    // sb21_cf has W02 = W12 = 60825.08 J/mol (large mgcf–nacf and fecf–nacf repulsion).
    // The critical temperature for binary nacf–(mgcf+fecf) mixing is ~3657 K >> T = 1473.15 K,
    // so a bulk composition of y ≈ (0.30, 0.20, 0.50) lies deep inside the solvus.
    //
    // The left branch covers the (Mg,Fe)-calcioferrite-rich field (y[2] ≤ 0.10) and the right
    // branch covers the Na-calcioferrite-rich field (y[2] ≥ 0.90).
    // Proxy species: CaCO3(s) → mgcf, MgCO3(s) → fecf, SiO2(s) → nacf.
    //
    // This test verifies that the activity model correctly identifies the bulk composition as
    // lying inside the two-phase field and requests a branch split (one rebuild).
    const auto db = test::createDatabase();
    const auto species = SpeciesList({db.species().get("CaCO3(s)"), db.species().get("MgCO3(s)"), db.species().get("SiO2(s)")});

    // Left branch: (Mg,Fe)-calcioferrite-rich — nacf fraction (y[2]) constrained to <= 10 %
    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "left";
    left.lowerBounds = ArrayXr::Zero(3);
    left.upperBounds = ArrayXr::Ones(3);
    left.upperBounds[2] = 0.10;

    // Right branch: Na-calcioferrite-rich — nacf fraction (y[2]) constrained to >= 90 %
    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "right";
    right.lowerBounds = ArrayXr::Zero(3);
    right.upperBounds = ArrayXr::Ones(3);
    right.lowerBounds[2] = 0.90;

    MAGEMinSB21CalcioferriteOptions modelOptions;
    modelOptions.branchPolicy.branches = {left, right};
    modelOptions.branchPolicy.stabilityCriterion = NamedGlobalizedSolidSolutionStabilityCriterion(
        NamedGlobalizedSolidSolutionStabilityPolicy::MAGEMinPilotBranchAmbiguity,
        modelOptions.branchPolicy.branches);
    const auto model = MAGEMinSolidSolutionPilotModelSB21Calcioferrite(modelOptions);

    Phase prototype;
    prototype = prototype.withName("Calcioferrite");
    prototype = prototype.withSpecies(species);
    prototype = prototype.withStateOfMatter(StateOfMatter::Solid);
    prototype = prototype.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(prototype.species()));

    MAGEMinSolidSolutionPilotOptions pilotOptions;
    pilotOptions.branches = {left, right};

    const auto phase = MAGEMinSolidSolutionPilotPhase(prototype, model);
    const auto definition = MAGEMinSolidSolutionPilotDefinition(prototype, model, pilotOptions);
    const auto system = ChemicalSystem(db, PhaseList{phase});

    // Bulk y ≈ (0.30, 0.20, 0.50): 3 mol CaCO3 + 2 mol MgCO3 + 5 mol SiO2 (total = 10 mol).
    ChemicalState state(system);
    state.setTemperature(1473.15);
    state.setPressure(1.0e9);
    ArrayXr n(3);
    n << 3.0, 2.0, 5.0;
    state.setSpeciesAmounts(n);

    GlobalizedSolidSolutionEquilibriumRetryOptions retryOptions;
    retryOptions.definitions = {definition};
    retryOptions.maxRetries = 1;

    const auto result = equilibrateWithGlobalizedSolidSolutionSplits(state, retryOptions);

    // The model must recognise the bulk as lying in the two-phase field: one rebuild expected.
    CHECK(result.result.succeeded());
    CHECK(result.numRebuilds == 1);
    REQUIRE(result.system.phases().size() == 2);
    CHECK(result.system.phase(0).name() == "Calcioferrite#left");
    CHECK(result.system.phase(1).name() == "Calcioferrite#right");
}

TEST_CASE("Testing post-split sb21_cf equilibrium preserves separated branch compositions", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][EquilibriumRetry][Exsolution][Separated]")
{
    const auto db = test::createDatabase();
    const auto species = SpeciesList({db.species().get("CaCO3(s)"), db.species().get("MgCO3(s)"), db.species().get("SiO2(s)")});

    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "left";
    left.lowerBounds = ArrayXr::Zero(3);
    left.upperBounds = ArrayXr::Ones(3);
    left.upperBounds[2] = 0.10;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "right";
    right.lowerBounds = ArrayXr::Zero(3);
    right.upperBounds = ArrayXr::Ones(3);
    right.lowerBounds[2] = 0.90;

    MAGEMinSB21CalcioferriteOptions modelOptions;
    modelOptions.branchPolicy.branches = {left, right};
    const auto model = MAGEMinSolidSolutionPilotModelSB21Calcioferrite(modelOptions);

    Phase prototype;
    prototype = prototype.withName("Calcioferrite");
    prototype = prototype.withSpecies(species);
    prototype = prototype.withStateOfMatter(StateOfMatter::Solid);
    prototype = prototype.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(prototype.species()));

    MAGEMinSolidSolutionPilotOptions pilotOptions;
    pilotOptions.branches = {left, right};

    const auto phase = MAGEMinSolidSolutionPilotPhase(prototype, model);
    const auto definition = MAGEMinSolidSolutionPilotDefinition(prototype, model, pilotOptions);
    const auto system = ChemicalSystem(db, PhaseList{phase});

    ChemicalState state(system);
    state.setTemperature(1473.15);
    state.setPressure(1.0e9);
    ArrayXr n(3);
    n << 3.0, 2.0, 5.0;
    state.setSpeciesAmounts(n);

    GlobalizedSolidSolutionEquilibriumRetryOptions retryOptions;
    retryOptions.definitions = {definition};
    retryOptions.maxRetries = 1;

    const auto result = equilibrateWithGlobalizedSolidSolutionSplits(state, retryOptions);

    CHECK(result.result.succeeded());
    CHECK(result.numRebuilds == 1);
    REQUIRE(result.system.phases().size() == 2);
    CHECK(result.system.phase(0).name() == "Calcioferrite#left");
    CHECK(result.system.phase(1).name() == "Calcioferrite#right");

    const auto leftAmounts = result.state.speciesAmountsInPhase(0);
    const auto rightAmounts = result.state.speciesAmountsInPhase(1);
    REQUIRE(leftAmounts.size() == 3);
    REQUIRE(rightAmounts.size() == 3);

    const auto leftTotal = static_cast<double>(leftAmounts.sum());
    const auto rightTotal = static_cast<double>(rightAmounts.sum());
    REQUIRE(leftTotal > 0.0);
    REQUIRE(rightTotal > 0.0);

    const auto leftNacf = static_cast<double>(leftAmounts[2]/leftTotal);
    const auto rightNacf = static_cast<double>(rightAmounts[2]/rightTotal);

    // The asymmetric branch seeds published by the split request should survive the rebuild
    // and keep the two duplicated phases on opposite sides of the nacf solvus.
    CHECK(leftNacf < 0.25);
    CHECK(rightNacf > 0.75);
    CHECK(rightNacf - leftNacf > 0.50);
}

TEST_CASE("Testing physically grounded solvus activity signatures for SB21 calcioferrite pilot (nacf solvus)", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Exsolution]")
{
    // Physical interpretation: the sb21_cf pilot model (W02 = W12 = 60825.08 J/mol) has
    // a nacf–(mgcf+fecf) solvus.  The two coexisting arms at T = 1473.15 K, P = 1 GPa are
    // approximated by:
    //   Left arm  (Mg,Fe)-calcioferrite-rich: y ≈ (0.45, 0.45, 0.10)  [low nacf]
    //   Right arm Na-calcioferrite-rich:      y ≈ (0.05, 0.05, 0.90)  [high nacf]
    //
    // For the solvus to be thermodynamically meaningful, the nacf activity (∝ y[2]^2 at
    // ideal level) must be significantly lower in the (Mg,Fe)-rich arm than in the Na-rich
    // arm.  This test verifies that the sb21_cf activity model produces this activity
    // contrast, confirming the physical basis of the immiscibility.
    //
    // Note: because the branch bounds control candidate selection but not the internal
    // minimizer itself, we evaluate the SINGLE-model activity function at each arm
    // composition independently, relying on the pilot's standard unconstrained minimizer.

    const auto db = test::createDatabase();
    const auto species = SpeciesList({db.species().get("CaCO3(s)"), db.species().get("MgCO3(s)"), db.species().get("SiO2(s)")});

    const auto model = MAGEMinSolidSolutionPilotModelSB21Calcioferrite({});
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(model)(species);

    ActivityProps propsLeft  = ActivityProps::create(species.size());
    ActivityProps propsRight = ActivityProps::create(species.size());

    // Left solvus arm: (Mg,Fe)-calcioferrite-rich composition — low nacf (y[2] ≈ 0.10).
    ArrayXr yLeft(3);
    yLeft << 0.45, 0.45, 0.10;

    // Right solvus arm: Na-calcioferrite-rich composition — high nacf (y[2] ≈ 0.90).
    ArrayXr yRight(3);
    yRight << 0.05, 0.05, 0.90;

    fn(propsLeft,  {1473.15, 1.0e9, yLeft});
    fn(propsRight, {1473.15, 1.0e9, yRight});

    // Physical check: the nacf activity (species index 2) must be substantially lower
    // at the (Mg,Fe)-rich arm than at the Na-rich arm.
    // Ideal contribution alone: ln_a_left[2] = 2*log(0.10) ≈ -4.6, ln_a_right[2] = 2*log(0.90) ≈ -0.21.
    // Excess contributions from the volume-fraction Margules terms shift both values,
    // but the sign of the contrast must be preserved.
    CHECK(static_cast<double>(propsLeft.ln_a[2])  < -1.0);   // nacf activity very low in left arm
    CHECK(static_cast<double>(propsRight.ln_a[2]) > -1.0);   // nacf activity appreciable in right arm
    CHECK(static_cast<double>(propsLeft.ln_a[2])  < static_cast<double>(propsRight.ln_a[2]));

    // Internal compositions must be close to the supplied compositions (no branch constraints
    // force the minimizer away from the visible composition when it is already at a local minimum).
    REQUIRE(propsLeft.extra.count("MAGEMinSolidSolutionPilot::InternalComposition") == 1);
    REQUIRE(propsRight.extra.count("MAGEMinSolidSolutionPilot::InternalComposition") == 1);

    const auto internalLeft  = std::any_cast<ArrayXr>(propsLeft.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    const auto internalRight = std::any_cast<ArrayXr>(propsRight.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));

    // The two internal compositions must differ substantially on the nacf axis.
    CHECK(static_cast<double>(internalRight[2] - internalLeft[2]) > 0.30);
}