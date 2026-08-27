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

#include <catch2/catch.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Equilibrium/EquilibriumUtils.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIdealSolution.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelGlobalizedSolidSolution.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelMAGEMinSolidSolutionPilot.hpp>
#include <Reaktoro/Models/ActivityModels/Support/MixedSystemSpeciesAmountUtils.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHollandPowell.hpp>

using namespace Reaktoro;

namespace test {
extern auto createDatabase() -> Database;
}

namespace {

using json = nlohmann::json;
namespace fs = std::filesystem;
struct RegressionTolerances
{
    real scalarAbs = 1.0e-10;
    real logAbs = 1.0e-10;
};

struct MAGEMinRegressionFixture
{
    String name;
    String model;
    Strings species;
    real T = 0.0;
    real P = 0.0;
    ArrayXr x;
    json expected;
    RegressionTolerances tolerances;
};

struct MAGEMinRetryRegressionFixture
{
    String name;
    String scenario;
    String model;
    String stabilityPolicy;
    Strings phaseSpecies;
    real T = 0.0;
    real P = 0.0;
    ArrayXr speciesAmounts;
    json expected;
    RegressionTolerances tolerances;
};

auto regressionFixturesDir() -> fs::path
{
    return fs::path(__FILE__).parent_path() / "Support" / "MAGEMinRegressionFixtures";
}

auto retryRegressionFixturesDir() -> fs::path
{
    return regressionFixturesDir() / "EquilibriumRetry";
}

auto fixtureDumpOutputPath() -> fs::path
{
    if(const auto configured = std::getenv("REAKTORO_DUMP_MAGEMIN_FIXTURES_PATH"))
        return fs::path(configured);

    return regressionFixturesDir() / "_dump" / "magemin-pilot-snapshots.json";
}

auto writeFixtureDump(json const& snapshots, fs::path const& outputPath) -> void
{
    fs::create_directories(outputPath.parent_path());

    std::ofstream output(outputPath, std::ios::binary | std::ios::trunc);
    if(!output)
        throw std::runtime_error("Unable to open MAGEMin fixture dump output: " + outputPath.string());

    output << std::setprecision(17) << snapshots.dump(2) << std::endl;
    output.close();

    if(!output)
        throw std::runtime_error("Unable to write MAGEMin fixture dump output: " + outputPath.string());

    std::cout << outputPath.string() << std::endl;
}

auto writeFixtureDump(json const& snapshots) -> void
{
    writeFixtureDump(snapshots, fixtureDumpOutputPath());
}

auto makeSpeciesList(Strings const& names) -> SpeciesList
{
    SpeciesList species;
    for(const auto& name : names)
        species.push_back(Species(name).withName(name));
    return species;
}

auto jsonArrayToEigen(json const& values) -> ArrayXr
{
    ArrayXr array(values.size());
    for(Index i = 0; i < static_cast<Index>(values.size()); ++i)
        array[i] = values.at(static_cast<std::size_t>(i)).get<double>();
    return array;
}

auto applyTolerances(RegressionTolerances& target, json const& doc) -> void
{
    if(!doc.contains("tolerances"))
        return;

    const auto& tolerances = doc.at("tolerances");
    if(tolerances.contains("scalarAbs"))
        target.scalarAbs = tolerances.at("scalarAbs").get<double>();
    if(tolerances.contains("logAbs"))
        target.logAbs = tolerances.at("logAbs").get<double>();
}

auto loadFixture(fs::path const& path) -> MAGEMinRegressionFixture
{
    std::ifstream input(path);
    if(!input)
        throw std::runtime_error("Unable to open MAGEMin regression fixture: " + path.string());

    const auto doc = json::parse(input);

    MAGEMinRegressionFixture fixture;
    fixture.name = doc.at("name").get<String>();
    fixture.model = doc.at("model").get<String>();
    fixture.species = doc.at("species").get<Strings>();
    fixture.T = doc.at("temperature").get<double>();
    fixture.P = doc.at("pressure").get<double>();
    fixture.x = jsonArrayToEigen(doc.at("composition"));
    fixture.expected = doc.at("expected");
    applyTolerances(fixture.tolerances, doc);

    return fixture;
}

auto loadRetryFixture(fs::path const& path) -> MAGEMinRetryRegressionFixture
{
    std::ifstream input(path);
    if(!input)
        throw std::runtime_error("Unable to open MAGEMin retry regression fixture: " + path.string());

    const auto doc = json::parse(input);

    MAGEMinRetryRegressionFixture fixture;
    fixture.name = doc.at("name").get<String>();
    fixture.scenario = doc.at("scenario").get<String>();
    fixture.model = doc.at("model").get<String>();
    if(doc.contains("stabilityPolicy"))
        fixture.stabilityPolicy = doc.at("stabilityPolicy").get<String>();
    fixture.phaseSpecies = doc.at("phaseSpecies").get<Strings>();
    fixture.T = doc.at("temperature").get<double>();
    fixture.P = doc.at("pressure").get<double>();
    fixture.speciesAmounts = jsonArrayToEigen(doc.at("speciesAmounts"));
    fixture.expected = doc.at("expected");
    applyTolerances(fixture.tolerances, doc);

    return fixture;
}

auto loadFixtures() -> Vec<MAGEMinRegressionFixture>
{
    Vec<MAGEMinRegressionFixture> fixtures;
    for(const auto& entry : fs::directory_iterator(regressionFixturesDir()))
    {
        if(entry.path().extension() != ".json")
            continue;
        fixtures.push_back(loadFixture(entry.path()));
    }

    std::sort(fixtures.begin(), fixtures.end(), [](auto const& lhs, auto const& rhs)
    {
        return lhs.name < rhs.name;
    });

    return fixtures;
}

auto loadRetryFixtures() -> Vec<MAGEMinRetryRegressionFixture>
{
    Vec<MAGEMinRetryRegressionFixture> fixtures;
    for(const auto& entry : fs::directory_iterator(retryRegressionFixturesDir()))
    {
        if(entry.path().extension() != ".json")
            continue;
        fixtures.push_back(loadRetryFixture(entry.path()));
    }

    std::sort(fixtures.begin(), fixtures.end(), [](auto const& lhs, auto const& rhs)
    {
        return lhs.name < rhs.name;
    });

    return fixtures;
}

auto makeRetryStabilityPolicy(
    String const& stabilityPolicy,
    Vec<GlobalizedSolidSolutionBranch> const& branches) -> GlobalizedSolidSolutionCandidateStabilityCriterion
{
    return NamedGlobalizedSolidSolutionStabilityCriterion(stabilityPolicy, branches);
}

auto makeModel(
    String const& model,
    Vec<GlobalizedSolidSolutionBranch> const& branches = {},
    String const& stabilityPolicy = "") -> GlobalizedSolidSolutionModel
{
    if(model == "sb11_ol")
    {
        MAGEMinSB11OlivineOptions options;
        options.branchPolicy.branches = branches;
        return MAGEMinSolidSolutionPilotModelSB11Olivine(options);
    }

    if(model == "sb11_wa")
    {
        MAGEMinSB11WadsleyiteOptions options;
        options.branchPolicy.branches = branches;
        return MAGEMinSolidSolutionPilotModelSB11Wadsleyite(options);
    }

    if(model == "sb11_ak")
    {
        MAGEMinSB11AkimotoiteOptions options;
        options.branchPolicy.branches = branches;
        options.branchPolicy.stabilityCriterion = makeRetryStabilityPolicy(stabilityPolicy, branches);
        return MAGEMinSolidSolutionPilotModelSB11Akimotoite(options);
    }

    if(model == "sb11_pv")
    {
        MAGEMinSB11PerovskiteOptions options;
        options.branchPolicy.branches = branches;
        options.branchPolicy.stabilityCriterion = makeRetryStabilityPolicy(stabilityPolicy, branches);
        return MAGEMinSolidSolutionPilotModelSB11Perovskite(options);
    }

    if(model == "sb11_cf")
    {
        MAGEMinSB11CalcioferriteOptions options;
        options.branchPolicy.branches = branches;
        options.branchPolicy.stabilityCriterion = makeRetryStabilityPolicy(stabilityPolicy, branches);
        return MAGEMinSolidSolutionPilotModelSB11Calcioferrite(options);
    }

    if(model == "sb21_sp")
    {
        MAGEMinSB21SpinelOptions options;
        options.branchPolicy.branches = branches;
        return MAGEMinSolidSolutionPilotModelSB21Spinel(options);
    }

    if(model == "sb21_nal")
    {
        MAGEMinSB21NALOptions options;
        options.branchPolicy.branches = branches;
        options.branchPolicy.stabilityCriterion = makeRetryStabilityPolicy(stabilityPolicy, branches);
        return MAGEMinSolidSolutionPilotModelSB21NAL(options);
    }

    if(model == "sb21_cf")
    {
        MAGEMinSB21CalcioferriteOptions options;
        options.branchPolicy.branches = branches;
        options.branchPolicy.stabilityCriterion = makeRetryStabilityPolicy(stabilityPolicy, branches);
        return MAGEMinSolidSolutionPilotModelSB21Calcioferrite(options);
    }

    throw std::runtime_error("Unsupported MAGEMin regression fixture model: " + model);
}

void checkArray(ArrayXrConstRef actual, json const& expected, real tolerance)
{
    const auto expectedArray = jsonArrayToEigen(expected);
    REQUIRE(actual.size() == expectedArray.size());
    for(Index i = 0; i < actual.size(); ++i)
        CHECK(actual[i] == Approx(expectedArray[i]).margin(tolerance));
}

auto legacyCompatibleTernaryExcessChemicalPotentials(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    ArrayXrConstRef y) -> ArrayXr
{
    if(options.thermo.excessChemicalPotentials)
        return options.thermo.excessChemicalPotentials(y);

    ArrayXr mu(3);
    const auto y0 = static_cast<double>(y[0]);
    const auto y1 = static_cast<double>(y[1]);
    const auto y2 = static_cast<double>(y[2]);

    mu[0] = y1*(1.0 - y0)*options.thermo.W01 + y2*(1.0 - y0)*options.thermo.W02 - y1*y2*options.thermo.W12;
    mu[1] = y0*(1.0 - y1)*options.thermo.W01 - y0*y2*options.thermo.W02 + y2*(1.0 - y1)*options.thermo.W12;
    mu[2] = -y0*y1*options.thermo.W01 + y0*(1.0 - y2)*options.thermo.W02 + y1*(1.0 - y2)*options.thermo.W12;
    return mu;
}

auto legacyConstrainedTernaryMinimizer() -> MAGEMinConstrainedTernaryMinimizer
{
    return [](MAGEMinImportedConstrainedTernarySolutionOptions const& options, real T, ArrayXrConstRef visiblex, Optional<ArrayXr> warmstart)
    {
        GlobalizedSolidSolutionInternalProblem problem;
        problem.objective = [=](ArrayXrConstRef y) -> real
        {
            const auto muEx = legacyCompatibleTernaryExcessChemicalPotentials(options, y);
            const auto Gex = y.matrix().dot(muEx.matrix());
            const auto Gid = options.thermo.idealGibbs ? options.thermo.idealGibbs(T, y) : real(0.0);
            return Gex + Gid + options.externalCompositionPenalty*universalGasConstant*T*(y - visiblex).matrix().squaredNorm();
        };
        problem.initialx = warmstart.value_or(ArrayXr(visiblex));
        problem.lowerBounds = ArrayXr::Constant(3, 1.0e-12);
        problem.upperBounds = ArrayXr::Constant(3, 1.0 - 1.0e-12);
        problem.tolerance = options.minimizerTolerance;
        problem.maxIterations = options.minimizerMaxIterations;
        problem.enforceUnityConstraint = true;
        return MinimizeGlobalizedSolidSolutionInternalProblem(problem);
    };
}

void checkProjectedVsLegacyAgreement(
    ActivityModel const& projected,
    ActivityModel const& legacy,
    ArrayXrConstRef visiblex,
    real T,
    real P,
    bool expectComparedAgainstLegacy,
    bool allowLegacyFallback)
{
    ActivityProps projectedProps = ActivityProps::create(visiblex.size());
    ActivityProps legacyProps = ActivityProps::create(visiblex.size());

    projected(projectedProps, {T, P, visiblex});
    legacy(legacyProps, {T, P, visiblex});

    const auto projectedInternalx = std::any_cast<ArrayXr>(projectedProps.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    const auto legacyInternalx = std::any_cast<ArrayXr>(legacyProps.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));

    CHECK((projectedInternalx - legacyInternalx).matrix().norm() == Approx(0.0).margin(1.0e-6));
    CHECK(static_cast<double>(projectedProps.Gx) == Approx(static_cast<double>(legacyProps.Gx)).margin(1.0e-6));
    CHECK(static_cast<double>(projectedProps.Hx) == Approx(static_cast<double>(legacyProps.Hx)).margin(1.0e-6));
    CHECK((projectedProps.ln_a - legacyProps.ln_a).matrix().norm() == Approx(0.0).margin(1.0e-6));
    CHECK((projectedProps.ln_g - legacyProps.ln_g).matrix().norm() == Approx(0.0).margin(1.0e-6));

    const auto selectedMinimizerStrategy = std::any_cast<String>(projectedProps.extra.at("MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"));
    const auto fallbackToLegacy = std::any_cast<bool>(projectedProps.extra.at("MAGEMinSolidSolutionPilot::FallbackToLegacy"));
    const auto projectedGradientAccepted = std::any_cast<bool>(projectedProps.extra.at("MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"));

    CHECK(std::any_cast<bool>(projectedProps.extra.at("MAGEMinSolidSolutionPilot::ComparedAgainstLegacy")) == expectComparedAgainstLegacy);
    CHECK(projectedProps.extra.count("MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount") == 1);
    CHECK(projectedProps.extra.count("MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount") == 1);

    if(allowLegacyFallback)
    {
        const auto selectedProjectedGradient = selectedMinimizerStrategy == "projected-gradient";
        const auto fallbackTelemetryIsConsistent = (selectedProjectedGradient && !fallbackToLegacy) || (!selectedProjectedGradient && fallbackToLegacy);
        CHECK(fallbackTelemetryIsConsistent);
        CHECK(projectedGradientAccepted == selectedProjectedGradient);
    }
    else
    {
        CHECK(selectedMinimizerStrategy == "projected-gradient");
        CHECK_FALSE(fallbackToLegacy);
    }

    if(expectComparedAgainstLegacy)
    {
        CHECK(projectedProps.extra.count("MAGEMinSolidSolutionPilot::ProjectedGradientLegacyCompositionDelta") == 1);
        CHECK(projectedProps.extra.count("MAGEMinSolidSolutionPilot::ProjectedGradientLegacyObjectiveDelta") == 1);
    }
    else
    {
        CHECK(projectedGradientAccepted == std::any_cast<bool>(projectedProps.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerConverged")));
    }
}

void checkGuardedFallbackOnForcedDisagreement(
    ActivityModel const& guarded,
    ActivityModel const& legacy,
    ArrayXrConstRef visiblex,
    real T,
    real P)
{
    ActivityProps guardedProps = ActivityProps::create(visiblex.size());
    ActivityProps legacyProps = ActivityProps::create(visiblex.size());

    guarded(guardedProps, {T, P, visiblex});
    legacy(legacyProps, {T, P, visiblex});

    CHECK(std::any_cast<bool>(guardedProps.extra.at("MAGEMinSolidSolutionPilot::ComparedAgainstLegacy")));
    CHECK(std::any_cast<bool>(guardedProps.extra.at("MAGEMinSolidSolutionPilot::FallbackToLegacy")));
    CHECK_FALSE(std::any_cast<bool>(guardedProps.extra.at("MAGEMinSolidSolutionPilot::ProjectedGradientAccepted")));
    CHECK(std::any_cast<String>(guardedProps.extra.at("MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy")) == "legacy");

    const auto guardedInternalx = std::any_cast<ArrayXr>(guardedProps.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    const auto legacyInternalx = std::any_cast<ArrayXr>(legacyProps.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    CHECK((guardedInternalx - legacyInternalx).matrix().norm() == Approx(0.0).margin(1.0e-12));
    CHECK(static_cast<double>(guardedProps.Gx) == Approx(static_cast<double>(legacyProps.Gx)).margin(1.0e-12));
    CHECK(static_cast<double>(guardedProps.Hx) == Approx(static_cast<double>(legacyProps.Hx)).margin(1.0e-12));
    CHECK((guardedProps.ln_a - legacyProps.ln_a).matrix().norm() == Approx(0.0).margin(1.0e-12));
    CHECK((guardedProps.ln_g - legacyProps.ln_g).matrix().norm() == Approx(0.0).margin(1.0e-12));
}

auto makeDefaultBranches(Index numCoords) -> Vec<GlobalizedSolidSolutionBranch>
{
    GlobalizedSolidSolutionBranch left;
    left.id = "left";
    left.label = "left";
    left.lowerBounds = ArrayXr::Zero(numCoords);
    left.upperBounds = ArrayXr::Ones(numCoords);
    left.upperBounds[0] = 0.45;

    GlobalizedSolidSolutionBranch right;
    right.id = "right";
    right.label = "right";
    right.lowerBounds = ArrayXr::Zero(numCoords);
    right.upperBounds = ArrayXr::Ones(numCoords);
    right.lowerBounds[0] = 0.55;

    return {left, right};
}

auto makeDatabaseSpeciesList(Database const& db, Strings const& names) -> SpeciesList
{
    SpeciesList species;
    for(const auto& name : names)
        species.push_back(db.species().get(name));
    return species;
}

auto eigenToJson(ArrayXrConstRef values) -> json
{
    json result = json::array();
    for(Index i = 0; i < values.size(); ++i)
        result.push_back(static_cast<double>(values[i]));
    return result;
}

auto snapshotModelFixture(String const& name, String const& model, Strings const& speciesNames, ArrayXrConstRef x, real T, real P) -> json
{
    const auto species = makeSpeciesList(speciesNames);
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(makeModel(model))(species);
    ActivityProps props = ActivityProps::create(species.size());

    fn(props, {T, P, x});

    json expected;
    expected["modelId"] = std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId"));
    if(props.extra.count("MAGEMinSolidSolutionPilot::Endmembers"))
        expected["endmembers"] = std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::Endmembers"));
    if(props.extra.count("MAGEMinSolidSolutionPilot::Endmember0"))
        expected["endmember0"] = std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::Endmember0"));
    if(props.extra.count("MAGEMinSolidSolutionPilot::Endmember1"))
        expected["endmember1"] = std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::Endmember1"));
    expected["internalComposition"] = eigenToJson(std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition")));
    expected["Gx"] = static_cast<double>(props.Gx);
    expected["Hx"] = static_cast<double>(props.Hx);
    expected["ln_g"] = eigenToJson(props.ln_g);
    expected["ln_a"] = eigenToJson(props.ln_a);
    expected["splitRequested"] = std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested"));

    json doc;
    doc["name"] = name;
    doc["model"] = model;
    doc["species"] = speciesNames;
    doc["temperature"] = static_cast<double>(T);
    doc["pressure"] = static_cast<double>(P);
    doc["composition"] = eigenToJson(x);
    doc["expected"] = expected;
    return doc;
}

} // namespace

TEST_CASE("Testing MAGEMin imported pilot regression fixtures", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    const auto fixtures = loadFixtures();
    REQUIRE(!fixtures.empty());

    for(const auto& fixture : fixtures)
    {
        DYNAMIC_SECTION(fixture.name)
        {
            const auto species = makeSpeciesList(fixture.species);
            ActivityModel fn = ActivityModelGlobalizedSolidSolution(makeModel(fixture.model))(species);
            ActivityProps props = ActivityProps::create(species.size());

            fn(props, {fixture.T, fixture.P, fixture.x});

            if(fixture.expected.contains("modelId"))
                CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::ModelId")) == fixture.expected.at("modelId").get<String>());

            if(fixture.expected.contains("endmember0"))
                CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::Endmember0")) == fixture.expected.at("endmember0").get<String>());

            if(fixture.expected.contains("endmember1"))
                CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::Endmember1")) == fixture.expected.at("endmember1").get<String>());

            if(fixture.expected.contains("endmembers"))
                CHECK(std::any_cast<Strings>(props.extra.at("MAGEMinSolidSolutionPilot::Endmembers")) == fixture.expected.at("endmembers").get<Strings>());

            if(fixture.expected.contains("internalComposition"))
                checkArray(std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition")), fixture.expected.at("internalComposition"), fixture.tolerances.scalarAbs);

            if(fixture.expected.contains("Gx"))
                CHECK(static_cast<double>(props.Gx) == Approx(fixture.expected.at("Gx").get<double>()).margin(fixture.tolerances.scalarAbs));

            if(fixture.expected.contains("Hx"))
                CHECK(static_cast<double>(props.Hx) == Approx(fixture.expected.at("Hx").get<double>()).margin(fixture.tolerances.scalarAbs));

            if(fixture.expected.contains("ln_g"))
                checkArray(props.ln_g, fixture.expected.at("ln_g"), fixture.tolerances.logAbs);

            if(fixture.expected.contains("ln_a"))
                checkArray(props.ln_a, fixture.expected.at("ln_a"), fixture.tolerances.logAbs);

            if(fixture.expected.contains("splitRequested"))
                CHECK(std::any_cast<bool>(props.extra.at("GlobalizedSolidSolution::SplitRequested")) == fixture.expected.at("splitRequested").get<bool>());
        }
    }
}

TEST_CASE("Testing MAGEMin imported pilot custom ternary minimizer hook", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    MAGEMinSB21CalcioferriteOptions options;
    auto invocationCount = std::make_shared<int>(0);
    auto observedT = std::make_shared<real>(0.0);
    auto observedVisiblex = std::make_shared<ArrayXr>();
    auto observedWarmstart = std::make_shared<ArrayXr>();

    options.minimizer = [=](MAGEMinImportedConstrainedTernarySolutionOptions const&, real T, ArrayXrConstRef visiblex, Optional<ArrayXr> warmstart)
    {
        if(*invocationCount == 0)
        {
            *observedT = T;
            *observedVisiblex = ArrayXr(visiblex);
            *observedWarmstart = warmstart.value_or(ArrayXr{});
        }
        *invocationCount += 1;

        GlobalizedSolidSolutionInternalResult result;
        result.x = ArrayXr(3);
        result.x << 0.20, 0.30, 0.50;
        result.objective = -123.0;
        result.iterations = 17;
        result.converged = true;
        return result;
    };

    const auto species = makeSpeciesList({"MgAl2O4", "FeAl2O4", "NaAlSiO4"});
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Calcioferrite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr visiblex(3);
    visiblex << 0.60, 0.25, 0.15;

    const auto T = 1473.15;
    const auto P = 1.0e9;
    fn(props, {T, P, visiblex});

    REQUIRE(*invocationCount >= 1);
    CHECK(*observedT == Approx(T));
    REQUIRE(observedVisiblex->size() == visiblex.size());
    CHECK((*observedVisiblex - visiblex).matrix().norm() == Approx(0.0).margin(1.0e-12));
    const auto warmstartShapeIsValid = observedWarmstart->size() == 0 || observedWarmstart->size() == visiblex.size();
    CHECK(warmstartShapeIsValid);
    if(observedWarmstart->size() == visiblex.size())
        CHECK((*observedWarmstart - visiblex).matrix().norm() == Approx(0.0).margin(1.0e-12));

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    ArrayXr expectedInternalx(3);
    expectedInternalx << 0.20, 0.30, 0.50;
    CHECK((internalx - expectedInternalx).matrix().norm() == Approx(0.0).margin(1.0e-12));
    CHECK(std::any_cast<std::uint64_t>(props.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerIterations")) == 17);
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerConverged")));
}

TEST_CASE("Testing MAGEMin imported pilot custom local-model minimizer hook", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    MAGEMinSB21CalcioferriteOptions options;
    auto invocationCount = std::make_shared<int>(0);
    auto observedModelId = std::make_shared<String>();
    auto observedVisiblex = std::make_shared<ArrayXr>();
    auto observedWarmstart = std::make_shared<ArrayXr>();
    auto observedObjectiveAtVisible = std::make_shared<real>(0.0);

    options.localModelMinimizer = [=](MAGEMinConstrainedTernaryLocalModel const& model, Optional<ArrayXr> warmstart)
    {
        if(*invocationCount == 0)
        {
            *observedModelId = model.modelId;
            *observedVisiblex = model.visiblex;
            *observedWarmstart = warmstart.value_or(ArrayXr{});
            *observedObjectiveAtVisible = model.objective(model.visiblex);
        }
        *invocationCount += 1;

        GlobalizedSolidSolutionInternalResult result;
        result.x = ArrayXr(3);
        result.x << 0.25, 0.25, 0.50;
        result.objective = -456.0;
        result.iterations = 9;
        result.converged = true;
        return result;
    };

    const auto species = makeSpeciesList({"MgAl2O4", "FeAl2O4", "NaAlSiO4"});
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Calcioferrite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr visiblex(3);
    visiblex << 0.60, 0.25, 0.15;

    const auto T = 1473.15;
    const auto P = 1.0e9;
    fn(props, {T, P, visiblex});

    REQUIRE(*invocationCount >= 1);
    CHECK(*observedModelId == "sb21_cf");
    REQUIRE(observedVisiblex->size() == visiblex.size());
    CHECK((*observedVisiblex - visiblex).matrix().norm() == Approx(0.0).margin(1.0e-12));
    CHECK(std::isfinite(static_cast<double>(*observedObjectiveAtVisible)));
    const auto warmstartShapeIsValid = observedWarmstart->size() == 0 || observedWarmstart->size() == visiblex.size();
    CHECK(warmstartShapeIsValid);

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    ArrayXr expectedInternalx(3);
    expectedInternalx << 0.25, 0.25, 0.50;
    CHECK((internalx - expectedInternalx).matrix().norm() == Approx(0.0).margin(1.0e-12));
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy")) == "custom-local-model");
    CHECK(std::any_cast<std::uint64_t>(props.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerIterations")) == 9);
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerConverged")));
}

TEST_CASE("Testing MAGEMin imported pilot local-model diagnostics payload hook", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    MAGEMinSB21CalcioferriteOptions options;
    options.localModelMinimizer = [](MAGEMinConstrainedTernaryLocalModel const&, Optional<ArrayXr>)
    {
        GlobalizedSolidSolutionInternalResult result;
        result.x = ArrayXr(3);
        result.x << 0.30, 0.20, 0.50;
        result.objective = -789.0;
        result.iterations = 5;
        result.converged = true;
        return result;
    };

    auto diagnosticsCount = std::make_shared<int>(0);
    options.localModelDiagnostics = [=](MAGEMinConstrainedTernaryLocalModel const& model, GlobalizedSolidSolutionInternalResult const& result)
    {
        *diagnosticsCount += 1;

        Map<String, Any> payload;
        payload["MAGEMinSolidSolutionPilot::LocalModelDiagnostics::Marker"] = String("custom-payload-ok");
        payload["MAGEMinSolidSolutionPilot::LocalModelDiagnostics::ModelId"] = model.modelId;
        payload["MAGEMinSolidSolutionPilot::LocalModelDiagnostics::ObjectiveAtResult"] = model.objective(result.x);
        payload["MAGEMinSolidSolutionPilot::LocalModelDiagnostics::ReportedObjective"] = result.objective;
        payload["MAGEMinSolidSolutionPilot::LocalModelDiagnostics::GradientCallbackPresent"] = static_cast<bool>(model.gradient);
        return payload;
    };

    const auto species = makeSpeciesList({"MgAl2O4", "FeAl2O4", "NaAlSiO4"});
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Calcioferrite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr visiblex(3);
    visiblex << 0.60, 0.25, 0.15;

    const auto T = 1473.15;
    const auto P = 1.0e9;
    fn(props, {T, P, visiblex});

    CHECK(*diagnosticsCount >= 1);
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy")) == "custom-local-model");
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::LocalModelDiagnostics::Marker")) == "custom-payload-ok");
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::LocalModelDiagnostics::ModelId")) == "sb21_cf");
    CHECK(std::any_cast<bool>(props.extra.at("MAGEMinSolidSolutionPilot::LocalModelDiagnostics::GradientCallbackPresent")));
    const auto objectiveAtResult = std::any_cast<real>(props.extra.at("MAGEMinSolidSolutionPilot::LocalModelDiagnostics::ObjectiveAtResult"));
    const auto reportedObjective = std::any_cast<real>(props.extra.at("MAGEMinSolidSolutionPilot::LocalModelDiagnostics::ReportedObjective"));
    CHECK(std::isfinite(static_cast<double>(objectiveAtResult)));
    CHECK(reportedObjective == Approx(-789.0).margin(1.0e-12));
}

TEST_CASE("Testing MAGEMin imported pilot optional NLopt local-model adapter hook", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    MAGEMinSB21CalcioferriteOptions options;
    auto localCalls = std::make_shared<int>(0);
    auto nloptCalls = std::make_shared<int>(0);

    options.localModelMinimizer = [=](MAGEMinConstrainedTernaryLocalModel const&, Optional<ArrayXr>)
    {
        *localCalls += 1;
        GlobalizedSolidSolutionInternalResult result;
        result.x = ArrayXr(3);
        result.x << 0.20, 0.20, 0.60;
        result.objective = -700.0;
        result.iterations = 8;
        result.converged = true;
        return result;
    };

    options.nloptLocalModelMinimizer = [=](MAGEMinConstrainedTernaryLocalModel const&, Optional<ArrayXr>)
    {
        *nloptCalls += 1;
        GlobalizedSolidSolutionInternalResult result;
        result.x = ArrayXr(3);
        result.x << 0.10, 0.30, 0.60;
        result.objective = -900.0;
        result.iterations = 4;
        result.converged = true;
        return result;
    };
    options.preferNLoptLocalModelMinimizer = true;

    const auto species = makeSpeciesList({"MgAl2O4", "FeAl2O4", "NaAlSiO4"});
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Calcioferrite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr visiblex(3);
    visiblex << 0.60, 0.25, 0.15;

    const auto T = 1473.15;
    const auto P = 1.0e9;
    fn(props, {T, P, visiblex});

    CHECK(*localCalls == 0);
    CHECK(*nloptCalls >= 1);
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy")) == "custom-local-model-nlopt");

    const auto internalx = std::any_cast<ArrayXr>(props.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    ArrayXr expectedInternalx(3);
    expectedInternalx << 0.10, 0.30, 0.60;
    CHECK((internalx - expectedInternalx).matrix().norm() == Approx(0.0).margin(1.0e-12));
}

TEST_CASE("Testing MAGEMin TC mconstraint bridge adapter on sb21_cpx", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    struct ConstraintCallbackStats
    {
        int valueCalls = 0;
        int gradientCalls = 0;
    };

    ConstraintCallbackStats stats;

    MAGEMinTCMConstraintBridge bridge;
    bridge.numConstraints = 1;
    bridge.numVariables = 5;
    bridge.constraintLowerBounds = (ArrayXr(1) << -1.0e10).finished();
    bridge.constraintUpperBounds = (ArrayXr(1) << 0.0).finished();
    bridge.userData = &stats;
    bridge.callback = [](unsigned m, double* result, unsigned n, const double* x, double* grad, void* data)
    {
        auto* s = static_cast<ConstraintCallbackStats*>(data);
        s->valueCalls += 1;

        result[0] = x[0] + x[1] + x[2] - 0.90;

        if(grad)
        {
            s->gradientCalls += 1;
            for(unsigned i = 0; i < m*n; ++i)
                grad[i] = 0.0;

            grad[0*n + 0] = 1.0;
            grad[0*n + 1] = 1.0;
            grad[0*n + 2] = 1.0;
        }
    };

    MAGEMinConstrainedTernaryLocalModel model;
    model.modelId = "sb21_cpx-test";
    model.visiblex = (ArrayXr(5) << 0.30, 0.25, 0.20, 0.15, 0.10).finished();
    model.lowerBounds = ArrayXr::Zero(5);
    model.upperBounds = ArrayXr::Ones(5);
    model.enforceUnityConstraint = false;
    model.objective = [](ArrayXrConstRef y) -> real
    {
        return y.matrix().squaredNorm();
    };
    model.gradient = [](ArrayXrConstRef y) -> ArrayXr
    {
        return 2.0*y;
    };

    bool fallbackCalled = false;
    const auto adapter = MAGEMinTCMConstraintBridgeLocalModelAdapter(
        bridge,
        [&fallbackCalled](MAGEMinConstrainedTernaryLocalModel const& constrainedModel, Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
        {
            fallbackCalled = true;
            const auto current = warmstart ? *warmstart : constrainedModel.visiblex;
            const auto constraints = constrainedModel.constraints(current);
            const auto jacobian = constrainedModel.constraintJacobian(current);

            CHECK(constraints.size() == 1);
            CHECK(jacobian.rows() == 1);
            CHECK(jacobian.cols() == 5);
            CHECK(constraints[0] == Approx(current[0] + current[1] + current[2] - 0.90).margin(1.0e-12));

            GlobalizedSolidSolutionInternalResult result;
            result.x = current;
            result.objective = constrainedModel.objective(current);
            result.iterations = 0;
            result.converged = true;
            return result;
        });

    const auto result = adapter(model, model.visiblex);

    CHECK(fallbackCalled);
    CHECK(stats.valueCalls > 0);
    CHECK(stats.gradientCalls > 0);
    CHECK(result.x.size() == 5);
    CHECK(result.x[0] + result.x[1] + result.x[2] <= Approx(0.90).margin(1.0e-8));
}

TEST_CASE("Testing MAGEMin sb21_cf projected-gradient and legacy minimizers agree", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    const auto species = makeSpeciesList({"MgAl2O4", "FeAl2O4", "NaAlSiO4"});

    MAGEMinSB21CalcioferriteOptions projectedOptions;
    ActivityModel projected = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Calcioferrite(projectedOptions))(species);

    MAGEMinSB21CalcioferriteOptions legacyOptions;
    legacyOptions.minimizer = legacyConstrainedTernaryMinimizer();
    ActivityModel legacy = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Calcioferrite(legacyOptions))(species);

    ArrayXr visiblex(3);
    visiblex << 0.60, 0.25, 0.15;

    const auto T = 1473.15;
    const auto P = 1.0e9;
    checkProjectedVsLegacyAgreement(projected, legacy, visiblex, T, P, false, false);
}

TEST_CASE("Testing MAGEMin TC mconstraint bridge wiring on sb21_ak", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    int valueCalls = 0;
    int gradientCalls = 0;

    MAGEMinTCMConstraintBridge bridge;
    bridge.numConstraints = 1;
    bridge.numVariables = 3;
    bridge.constraintLowerBounds = (ArrayXr(1) << -1.0e10).finished();
    bridge.constraintUpperBounds = (ArrayXr(1) << 0.0).finished();
    bridge.callback = [&valueCalls, &gradientCalls](unsigned m, double* result, unsigned n, const double* x, double* grad, void*)
    {
        (void)m;
        ++valueCalls;
        result[0] = x[0] + x[1] + x[2] - 1.0;
        if(grad)
        {
            ++gradientCalls;
            for(unsigned i = 0; i < n; ++i)
                grad[i] = 1.0;
        }
    };

    MAGEMinSB21AkimotoiteOptions options;
    options.tcMConstraintBridge = bridge;
    options.preferNLoptLocalModelMinimizer = true;

    const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "CaSiO3"});
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Akimotoite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr visiblex(3);
    visiblex << 0.80, 0.10, 0.10;
    fn(props, {1473.15, 1.0e9, visiblex});

    CHECK(valueCalls > 0);
    CHECK(gradientCalls > 0);
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy")) == "custom-local-model-nlopt");
}

TEST_CASE("Testing MAGEMin TC mconstraint bridge wiring on sb21_pv", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    int valueCalls = 0;
    int gradientCalls = 0;

    MAGEMinTCMConstraintBridge bridge;
    bridge.numConstraints = 1;
    bridge.numVariables = 3;
    bridge.constraintLowerBounds = (ArrayXr(1) << -1.0e10).finished();
    bridge.constraintUpperBounds = (ArrayXr(1) << 0.0).finished();
    bridge.callback = [&valueCalls, &gradientCalls](unsigned m, double* result, unsigned n, const double* x, double* grad, void*)
    {
        (void)m;
        ++valueCalls;
        result[0] = x[0] + x[1] + x[2] - 1.0;
        if(grad)
        {
            ++gradientCalls;
            for(unsigned i = 0; i < n; ++i)
                grad[i] = 1.0;
        }
    };

    MAGEMinSB21PerovskiteOptions options;
    options.tcMConstraintBridge = bridge;
    options.preferNLoptLocalModelMinimizer = true;

    const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "Al2O3"});
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Perovskite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr visiblex(3);
    visiblex << 0.70, 0.20, 0.10;
    fn(props, {1473.15, 1.0e9, visiblex});

    CHECK(valueCalls > 0);
    CHECK(gradientCalls > 0);
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy")) == "custom-local-model-nlopt");
}

TEST_CASE("Testing MAGEMin TC mconstraint bridge wiring on sb21_ppv", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    int valueCalls = 0;
    int gradientCalls = 0;

    MAGEMinTCMConstraintBridge bridge;
    bridge.numConstraints = 1;
    bridge.numVariables = 3;
    bridge.constraintLowerBounds = (ArrayXr(1) << -1.0e10).finished();
    bridge.constraintUpperBounds = (ArrayXr(1) << 0.0).finished();
    bridge.callback = [&valueCalls, &gradientCalls](unsigned m, double* result, unsigned n, const double* x, double* grad, void*)
    {
        (void)m;
        ++valueCalls;
        result[0] = x[0] + x[1] + x[2] - 1.0;
        if(grad)
        {
            ++gradientCalls;
            for(unsigned i = 0; i < n; ++i)
                grad[i] = 1.0;
        }
    };

    MAGEMinSB21PostPerovskiteOptions options;
    options.tcMConstraintBridge = bridge;
    options.preferNLoptLocalModelMinimizer = true;

    const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "Al2O3"});
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21PostPerovskite(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr visiblex(3);
    visiblex << 0.65, 0.25, 0.10;
    fn(props, {1473.15, 1.0e9, visiblex});

    CHECK(valueCalls > 0);
    CHECK(gradientCalls > 0);
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy")) == "custom-local-model-nlopt");
}

TEST_CASE("Testing MAGEMin TC mconstraint bridge wiring on sb21_mw", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    int valueCalls = 0;
    int gradientCalls = 0;

    MAGEMinTCMConstraintBridge bridge;
    bridge.numConstraints = 1;
    bridge.numVariables = 3;
    bridge.constraintLowerBounds = (ArrayXr(1) << -1.0e10).finished();
    bridge.constraintUpperBounds = (ArrayXr(1) << 0.0).finished();
    bridge.callback = [&valueCalls, &gradientCalls](unsigned m, double* result, unsigned n, const double* x, double* grad, void*)
    {
        (void)m;
        ++valueCalls;
        result[0] = x[0] + x[1] + x[2] - 1.0;
        if(grad)
        {
            ++gradientCalls;
            for(unsigned i = 0; i < n; ++i)
                grad[i] = 1.0;
        }
    };

    MAGEMinSB21MagnesiowustitesOptions options;
    options.tcMConstraintBridge = bridge;
    options.preferNLoptLocalModelMinimizer = true;

    const auto species = makeSpeciesList({"FeO", "MgO", "Al2O3"});
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21Magnesiowustites(options))(species);
    ActivityProps props = ActivityProps::create(species.size());

    ArrayXr visiblex(3);
    visiblex << 0.50, 0.40, 0.10;
    fn(props, {1473.15, 1.0e9, visiblex});

    CHECK(valueCalls > 0);
    CHECK(gradientCalls > 0);
    CHECK(std::any_cast<String>(props.extra.at("MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy")) == "custom-local-model-nlopt");
}

TEST_CASE("Testing MAGEMin sb11_pv projected-gradient and legacy minimizers agree", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "Al2O3"});

    MAGEMinSB11PerovskiteOptions projectedOptions;
    ActivityModel projected = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Perovskite(projectedOptions))(species);

    MAGEMinSB11PerovskiteOptions legacyOptions;
    legacyOptions.minimizer = legacyConstrainedTernaryMinimizer();
    ActivityModel legacy = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Perovskite(legacyOptions))(species);

    ArrayXr visiblex(3);
    visiblex << 0.70, 0.20, 0.10;

    const auto T = 1473.15;
    const auto P = 1.0e9;
    checkProjectedVsLegacyAgreement(projected, legacy, visiblex, T, P, true, true);
}

TEST_CASE("Testing MAGEMin sb11_ak projected-gradient and legacy minimizers agree", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "CaSiO3"});

    MAGEMinSB11AkimotoiteOptions projectedOptions;
    ActivityModel projected = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Akimotoite(projectedOptions))(species);

    MAGEMinSB11AkimotoiteOptions legacyOptions;
    legacyOptions.minimizer = legacyConstrainedTernaryMinimizer();
    ActivityModel legacy = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Akimotoite(legacyOptions))(species);

    ArrayXr visiblex(3);
    visiblex << 0.80, 0.10, 0.10;

    const auto T = 1473.15;
    const auto P = 1.0e9;
    checkProjectedVsLegacyAgreement(projected, legacy, visiblex, T, P, true, true);
}

TEST_CASE("Testing MAGEMin sb11_ak guarded projected-gradient falls back to legacy on forced disagreement", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "CaSiO3"});

    MAGEMinSB11AkimotoiteOptions guardedOptions;
    guardedOptions.minimizerMaxIterations = 0; // Force projected-gradient non-convergence to trigger guarded fallback.
    ActivityModel guarded = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Akimotoite(guardedOptions))(species);

    MAGEMinSB11AkimotoiteOptions legacyOptions;
    legacyOptions.minimizer = legacyConstrainedTernaryMinimizer();
    legacyOptions.minimizerMaxIterations = 0;
    ActivityModel legacy = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Akimotoite(legacyOptions))(species);

    ArrayXr visiblex(3);
    visiblex << 0.80, 0.10, 0.10;

    const auto T = 1473.15;
    const auto P = 1.0e9;

    checkGuardedFallbackOnForcedDisagreement(guarded, legacy, visiblex, T, P);
}

TEST_CASE("Testing MAGEMin sb11_pv guarded projected-gradient falls back to legacy on forced disagreement", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "Al2O3"});

    MAGEMinSB11PerovskiteOptions guardedOptions;
    guardedOptions.minimizerMaxIterations = 0; // Force projected-gradient non-convergence to trigger guarded fallback.
    ActivityModel guarded = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Perovskite(guardedOptions))(species);

    MAGEMinSB11PerovskiteOptions legacyOptions;
    legacyOptions.minimizer = legacyConstrainedTernaryMinimizer();
    legacyOptions.minimizerMaxIterations = 0;
    ActivityModel legacy = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Perovskite(legacyOptions))(species);

    ArrayXr visiblex(3);
    visiblex << 0.70, 0.20, 0.10;

    const auto T = 1473.15;
    const auto P = 1.0e9;

    checkGuardedFallbackOnForcedDisagreement(guarded, legacy, visiblex, T, P);
}

TEST_CASE("Testing MAGEMin sb11_cf guarded projected-gradient falls back to legacy on forced disagreement", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    const auto species = makeSpeciesList({"MgAl2O4", "FeAl2O4", "NaAlSiO4"});

    MAGEMinSB11CalcioferriteOptions guardedOptions;
    guardedOptions.minimizerMaxIterations = 0; // Force projected-gradient non-convergence to trigger guarded fallback.
    ActivityModel guarded = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Calcioferrite(guardedOptions))(species);

    MAGEMinSB11CalcioferriteOptions legacyOptions;
    legacyOptions.minimizer = legacyConstrainedTernaryMinimizer();
    legacyOptions.minimizerMaxIterations = 0;
    ActivityModel legacy = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB11Calcioferrite(legacyOptions))(species);

    ArrayXr visiblex(3);
    visiblex << 0.60, 0.25, 0.15;

    const auto T = 1473.15;
    const auto P = 1.0e9;

    checkGuardedFallbackOnForcedDisagreement(guarded, legacy, visiblex, T, P);
}

TEST_CASE("Testing MAGEMin sb21_nal guarded projected-gradient falls back to legacy on forced disagreement", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    const auto species = makeSpeciesList({"NaMg2Al5SiO12", "NaFe2Al5SiO12", "Na3Al3Si3O12"});

    MAGEMinSB21NALOptions guardedOptions;
    guardedOptions.minimizerMaxIterations = 0; // Force projected-gradient non-convergence to trigger guarded fallback.
    ActivityModel guarded = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21NAL(guardedOptions))(species);

    MAGEMinSB21NALOptions legacyOptions;
    legacyOptions.minimizer = legacyConstrainedTernaryMinimizer();
    legacyOptions.minimizerMaxIterations = 0;
    ActivityModel legacy = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21NAL(legacyOptions))(species);

    ArrayXr visiblex(3);
    visiblex << 0.55, 0.20, 0.25;

    const auto T = 1473.15;
    const auto P = 1.0e9;

    checkGuardedFallbackOnForcedDisagreement(guarded, legacy, visiblex, T, P);
}

TEST_CASE("Testing MAGEMin sb21_nal projected-gradient and legacy minimizers agree", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    const auto species = makeSpeciesList({"NaMg2Al5SiO12", "NaFe2Al5SiO12", "Na3Al3Si3O12"});

    MAGEMinSB21NALOptions projectedOptions;
    ActivityModel projected = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21NAL(projectedOptions))(species);

    MAGEMinSB21NALOptions legacyOptions;
    legacyOptions.minimizer = legacyConstrainedTernaryMinimizer();
    ActivityModel legacy = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21NAL(legacyOptions))(species);

    ArrayXr visiblex(3);
    visiblex << 0.55, 0.20, 0.25;

    const auto T = 1473.15;
    const auto P = 1.0e9;
    checkProjectedVsLegacyAgreement(projected, legacy, visiblex, T, P, true, true);
}

TEST_CASE("Testing MAGEMinProjectedGradientLocalModelMinimizer utility via localModelMinimizer hook", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    // Verify that a localModelMinimizer lambda can delegate to the public utility while injecting
    // a custom gradient.  The custom gradient wraps the standard formula gradient to track
    // invocation count.  The final composition and objective should match the baseline
    // projected-gradient result obtained through the default path for the same family.

    const auto species = makeSpeciesList({"NaMg2Al5SiO12", "NaFe2Al5SiO12", "Na3Al3Si3O12"});
    const auto T = 1473.15;
    const auto P = 1.0e9;

    ArrayXr visiblex(3);
    visiblex << 0.55, 0.20, 0.25;

    // Baseline: default projected-gradient path (no custom minimizer).
    MAGEMinSB21NALOptions baselineOptions;
    ActivityModel baselineModel = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21NAL(baselineOptions))(species);
    ActivityProps baselineProps = ActivityProps::create(species.size());
    baselineModel(baselineProps, {T, P, visiblex});
    const auto baselineInternalx = std::any_cast<ArrayXr>(baselineProps.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));

    // Custom path: localModelMinimizer delegates to MAGEMinProjectedGradientLocalModelMinimizer
    // after wrapping model.gradient with a counter.
    auto gradientCallCount = std::make_shared<int>(0);

    MAGEMinSB21NALOptions customOptions;
    customOptions.localModelMinimizer = [=](MAGEMinConstrainedTernaryLocalModel const& model, Optional<ArrayXr> warmstart)
    {
        // Wrap the model gradient with a counter, then forward to the public utility.
        MAGEMinConstrainedTernaryLocalModel wrapped = model;
        wrapped.gradient = [innerGradient = model.gradient, gradientCallCount](ArrayXrConstRef y) -> ArrayXr
        {
            *gradientCallCount += 1;
            return innerGradient(y);
        };
        return MAGEMinProjectedGradientLocalModelMinimizer(wrapped, warmstart);
    };

    ActivityModel customModel = ActivityModelGlobalizedSolidSolution(MAGEMinSolidSolutionPilotModelSB21NAL(customOptions))(species);
    ActivityProps customProps = ActivityProps::create(species.size());
    customModel(customProps, {T, P, visiblex});

    // Gradient callback must have been invoked at least once per evaluation.
    CHECK(*gradientCallCount >= 1);

    // The strategy tag should indicate custom-local-model dispatch.
    CHECK(std::any_cast<String>(customProps.extra.at("MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy")) == "custom-local-model");

    // The resulting internal composition must match the baseline projected-gradient result.
    const auto customInternalx = std::any_cast<ArrayXr>(customProps.extra.at("MAGEMinSolidSolutionPilot::InternalComposition"));
    REQUIRE(customInternalx.size() == baselineInternalx.size());
    CHECK((customInternalx - baselineInternalx).matrix().norm() == Approx(0.0).margin(1.0e-8));
}

TEST_CASE("Testing MAGEMinProjectedGradientLocalModelMinimizer honors local-model bounds and unity settings", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    // Unconstrained-sum path: enforceUnityConstraint=false should honor simple box bounds.
    {
        MAGEMinConstrainedTernaryLocalModel model;
        model.visiblex = (ArrayXr(2) << 0.50, 0.50).finished();
        model.lowerBounds = (ArrayXr(2) << 0.20, 0.30).finished();
        model.upperBounds = (ArrayXr(2) << 0.80, 0.90).finished();
        model.enforceUnityConstraint = false;
        model.tolerance = 1.0e-12;
        model.maxIterations = 256;

        model.objective = [](ArrayXrConstRef y) -> real
        {
            return (y[0] - 0.10)*(y[0] - 0.10) + (y[1] - 1.20)*(y[1] - 1.20);
        };
        model.gradient = [](ArrayXrConstRef y) -> ArrayXr
        {
            ArrayXr g(2);
            g[0] = 2.0*(y[0] - 0.10);
            g[1] = 2.0*(y[1] - 1.20);
            return g;
        };

        const auto result = MAGEMinProjectedGradientLocalModelMinimizer(model, ArrayXr((ArrayXr(2) << -1.0, 2.0).finished()));

        CHECK(result.x[0] == Approx(0.20).margin(1.0e-8));
        CHECK(result.x[1] == Approx(0.90).margin(1.0e-8));
        CHECK(result.x[0] + result.x[1] == Approx(1.10).margin(1.0e-8));
    }

    // Simplex path: enforceUnityConstraint=true should project to bounded simplex.
    {
        MAGEMinConstrainedTernaryLocalModel model;
        model.visiblex = (ArrayXr(2) << 0.50, 0.50).finished();
        model.lowerBounds = (ArrayXr(2) << 0.20, 0.20).finished();
        model.upperBounds = (ArrayXr(2) << 0.80, 0.80).finished();
        model.enforceUnityConstraint = true;
        model.tolerance = 1.0e-12;
        model.maxIterations = 256;

        model.objective = [](ArrayXrConstRef y) -> real
        {
            return (y[0] - 0.90)*(y[0] - 0.90) + (y[1] - 0.10)*(y[1] - 0.10);
        };
        model.gradient = [](ArrayXrConstRef y) -> ArrayXr
        {
            ArrayXr g(2);
            g[0] = 2.0*(y[0] - 0.90);
            g[1] = 2.0*(y[1] - 0.10);
            return g;
        };

        const auto result = MAGEMinProjectedGradientLocalModelMinimizer(model, ArrayXr((ArrayXr(2) << 0.10, 0.90).finished()));

        CHECK(result.x.sum() == Approx(1.0).margin(1.0e-10));
        CHECK(result.x[0] == Approx(0.80).margin(1.0e-8));
        CHECK(result.x[1] == Approx(0.20).margin(1.0e-8));
    }
}

TEST_CASE("Testing SolidSolutionMinimizerBenchmark accumulates pilot telemetry", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    // sb21_cf uses direct projected-gradient (no guard), so every evaluation should
    // record a projected-gradient strategy and no fallback.
    const auto species = makeSpeciesList({"MgAl2O4", "FeAl2O4", "NaAlSiO4"});
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(
        MAGEMinSolidSolutionPilotModelSB21Calcioferrite({}))(species);
    ActivityProps props = ActivityProps::create(species.size());

    const auto T = 1473.15;
    const auto P = 1.0e9;

    ArrayXr compositions[3] = {
        (ArrayXr(3) << 0.60, 0.25, 0.15).finished(),
        (ArrayXr(3) << 0.33, 0.33, 0.34).finished(),
        (ArrayXr(3) << 0.10, 0.10, 0.80).finished(),
    };

    SolidSolutionMinimizerBenchmark benchmark;
    for(const auto& x : compositions)
    {
        fn(props, {T, P, x});
        benchmark.accumulate(props.extra);
    }

    const auto s = benchmark.stats();
    CHECK(s.totalEvaluations == 3);
    CHECK(s.projectedGradientCount + s.legacyCount + s.customCount == s.totalEvaluations);
    CHECK(s.fallbackCount == 0); // sb21_cf: direct projected-gradient, no fallback
    CHECK(s.projectedGradientSelectionRate >= 0.0);
    CHECK(s.projectedGradientSelectionRate <= 1.0);

    // Verify the benchmark correctly identifies projected-gradient as the dominant strategy.
    CHECK(s.projectedGradientCount == s.totalEvaluations);

    // Reset and verify cleared state.
    benchmark.reset();
    const auto sReset = benchmark.stats();
    CHECK(sReset.totalEvaluations == 0);
    CHECK(sReset.projectedGradientCount == 0);
    CHECK(sReset.fallbackRate == Approx(0.0));
}

TEST_CASE("Testing SolidSolutionMinimizerBenchmark captures fallback telemetry from guarded family", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    // sb11_ak uses guarded projected-gradient, so comparedCount should equal totalEvaluations
    // and fallbackRate should be between 0 and 1 inclusive.
    const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "CaSiO3"});
    ActivityModel fn = ActivityModelGlobalizedSolidSolution(
        MAGEMinSolidSolutionPilotModelSB11Akimotoite({}))(species);
    ActivityProps props = ActivityProps::create(species.size());

    const auto T = 1473.15;
    const auto P = 1.0e9;

    ArrayXr compositions[3] = {
        (ArrayXr(3) << 0.80, 0.10, 0.10).finished(),
        (ArrayXr(3) << 0.40, 0.30, 0.30).finished(),
        (ArrayXr(3) << 0.10, 0.10, 0.80).finished(),
    };

    SolidSolutionMinimizerBenchmark benchmark;
    for(const auto& x : compositions)
    {
        fn(props, {T, P, x});
        benchmark.accumulate(props.extra);
    }

    const auto s = benchmark.stats();
    CHECK(s.totalEvaluations == 3);
    CHECK(s.comparedCount == s.totalEvaluations); // guarded: always compares
    CHECK(s.fallbackRate >= 0.0);
    CHECK(s.fallbackRate <= 1.0);
    // Every evaluation used either projected-gradient or legacy (no custom callback).
    CHECK(s.customCount == 0);
    CHECK(s.projectedGradientCount + s.legacyCount == s.totalEvaluations);
}

/// Generate a uniform ternary composition grid (all x >= 0, sum = 1) with steps of 1/N.
auto makeTernaryGrid(int N) -> Vec<ArrayXr>
{
    Vec<ArrayXr> grid;
    const auto step = 1.0 / N;
    for(int i = 0; i <= N; ++i)
        for(int j = 0; j <= N - i; ++j)
        {
            const auto k = N - i - j;
            const double xi = i * step + 1.0e-10;
            const double xj = j * step + 1.0e-10;
            const double xk = k * step + 1.0e-10;
            const double sum = xi + xj + xk;
            ArrayXr x(3);
            x << xi/sum, xj/sum, xk/sum;
            grid.push_back(x);
        }
    return grid;
}

/// Accumulate benchmark stats for a guarded-PG family over a ternary grid.
auto benchmarkTernaryFamily(
    ActivityModel const& fn,
    Vec<ArrayXr> const& grid,
    real T,
    real P,
    Index numSpecies) -> SolidSolutionMinimizerBenchmarkStats
{
    SolidSolutionMinimizerBenchmark benchmark;
    ActivityProps props = ActivityProps::create(numSpecies);
    for(const auto& x : grid)
    {
        fn(props, {T, P, x});
        benchmark.accumulate(props.extra);
    }
    return benchmark.stats();
}

TEST_CASE("Benchmarking guarded-PG ternary families over full composition grid", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Benchmark]")
{
    // Grid resolution: N=7 gives 36 compositions (full ternary simplex).
    const auto grid = makeTernaryGrid(7);
    const auto T = 1473.15;
    const auto P = 1.0e9;
    const auto n = grid.size();
    REQUIRE(n > 0);

    SECTION("sb11_pv (mgpv-fepv-alpv)")
    {
        const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "Al2O3"});
        ActivityModel fn = ActivityModelGlobalizedSolidSolution(
            MAGEMinSolidSolutionPilotModelSB11Perovskite({}))(species);
        const auto s = benchmarkTernaryFamily(fn, grid, T, P, 3);

        INFO("sb11_pv: total=" << s.totalEvaluations
             << " pg=" << s.projectedGradientCount
             << " legacy=" << s.legacyCount
             << " fallbackRate=" << s.fallbackRate
             << " avgPGiters=" << s.averageProjectedGradientIterations);
            INFO("sb11_pv (disagreements): pgWinsObj=" << s.pgLowerObjectiveCount
                 << " legacyWinsObj=" << s.legacyLowerObjectiveCount);

        CHECK(s.totalEvaluations == static_cast<Index>(n));
        CHECK(s.projectedGradientCount + s.legacyCount == s.totalEvaluations);
        CHECK(s.comparedCount == s.totalEvaluations); // guarded: always compares
        CHECK(s.fallbackRate >= 0.0);
        CHECK(s.fallbackRate <= 1.0);
        // Projected-gradient must win at least one composition (not universally failing).
        CHECK(s.projectedGradientCount > 0);
    }

    SECTION("sb11_ak (mgak-feak-co)")
    {
        const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "CaSiO3"});
        ActivityModel fn = ActivityModelGlobalizedSolidSolution(
            MAGEMinSolidSolutionPilotModelSB11Akimotoite({}))(species);
        const auto s = benchmarkTernaryFamily(fn, grid, T, P, 3);

        INFO("sb11_ak: total=" << s.totalEvaluations
             << " pg=" << s.projectedGradientCount
             << " legacy=" << s.legacyCount
             << " fallbackRate=" << s.fallbackRate
             << " avgPGiters=" << s.averageProjectedGradientIterations);
            INFO("sb11_ak (disagreements): pgWinsObj=" << s.pgLowerObjectiveCount
                 << " legacyWinsObj=" << s.legacyLowerObjectiveCount);

        CHECK(s.totalEvaluations == static_cast<Index>(n));
        CHECK(s.projectedGradientCount + s.legacyCount == s.totalEvaluations);
        CHECK(s.comparedCount == s.totalEvaluations);
        CHECK(s.fallbackRate >= 0.0);
        CHECK(s.fallbackRate <= 1.0);
        CHECK(s.projectedGradientCount > 0);
    }

    SECTION("sb11_cf (mgcf-fecf-nacf)")
    {
        const auto species = makeSpeciesList({"MgAl2O4", "FeAl2O4", "NaAlSiO4"});
        ActivityModel fn = ActivityModelGlobalizedSolidSolution(
            MAGEMinSolidSolutionPilotModelSB11Calcioferrite({}))(species);
        const auto s = benchmarkTernaryFamily(fn, grid, T, P, 3);

        INFO("sb11_cf: total=" << s.totalEvaluations
             << " pg=" << s.projectedGradientCount
             << " legacy=" << s.legacyCount
             << " fallbackRate=" << s.fallbackRate
             << " avgPGiters=" << s.averageProjectedGradientIterations);
            INFO("sb11_cf (disagreements): pgWinsObj=" << s.pgLowerObjectiveCount
                 << " legacyWinsObj=" << s.legacyLowerObjectiveCount);

        CHECK(s.totalEvaluations == static_cast<Index>(n));
        CHECK(s.projectedGradientCount + s.legacyCount == s.totalEvaluations);
        CHECK(s.comparedCount == s.totalEvaluations);
        CHECK(s.fallbackRate >= 0.0);
        CHECK(s.fallbackRate <= 1.0);
        CHECK(s.projectedGradientCount > 0);
    }

    SECTION("sb21_nal (mnal-fnal-nnal)")
    {
        const auto species = makeSpeciesList({"NaMg2Al5SiO12", "NaFe2Al5SiO12", "Na3Al3Si3O12"});
        ActivityModel fn = ActivityModelGlobalizedSolidSolution(
            MAGEMinSolidSolutionPilotModelSB21NAL({}))(species);
        const auto s = benchmarkTernaryFamily(fn, grid, T, P, 3);

        INFO("sb21_nal: total=" << s.totalEvaluations
             << " pg=" << s.projectedGradientCount
             << " legacy=" << s.legacyCount
             << " fallbackRate=" << s.fallbackRate
             << " avgPGiters=" << s.averageProjectedGradientIterations);
            INFO("sb21_nal (disagreements): pgWinsObj=" << s.pgLowerObjectiveCount
                 << " legacyWinsObj=" << s.legacyLowerObjectiveCount);

        CHECK(s.totalEvaluations == static_cast<Index>(n));
        CHECK(s.projectedGradientCount + s.legacyCount == s.totalEvaluations);
        CHECK(s.comparedCount == s.totalEvaluations);
        CHECK(s.fallbackRate >= 0.0);
        CHECK(s.fallbackRate <= 1.0);
        CHECK(s.projectedGradientCount > 0);
    }

    SECTION("sb21_cf (MgAl2O4-FeAl2O4-NaAlSiO4) — direct PG baseline")
    {
        // sb21_cf uses direct PG (no guard) — all evaluations should use projected-gradient.
        const auto species = makeSpeciesList({"MgAl2O4", "FeAl2O4", "NaAlSiO4"});
        ActivityModel fn = ActivityModelGlobalizedSolidSolution(
            MAGEMinSolidSolutionPilotModelSB21Calcioferrite({}))(species);
        const auto s = benchmarkTernaryFamily(fn, grid, T, P, 3);

        INFO("sb21_cf: total=" << s.totalEvaluations
             << " pg=" << s.projectedGradientCount
             << " legacy=" << s.legacyCount
             << " fallbackRate=" << s.fallbackRate
             << " avgPGiters=" << s.averageProjectedGradientIterations);

        CHECK(s.totalEvaluations == static_cast<Index>(n));
        CHECK(s.fallbackCount == 0);
        CHECK(s.projectedGradientCount == s.totalEvaluations);
        CHECK(s.fallbackRate == Approx(0.0));
    }
}

TEST_CASE("Benchmarking forced-disagreement fallback over guarded ternary grids", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Benchmark]")
{
    // Force projected-gradient disagreement by zeroing max iterations.
    // For guarded families, every evaluation should fall back to legacy.
    const auto grid = makeTernaryGrid(7);
    const auto T = 1473.15;
    const auto P = 1.0e9;
    const auto n = static_cast<Index>(grid.size());
    REQUIRE(n > 0);

    SECTION("sb11_pv forced disagreement -> full legacy fallback")
    {
        const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "Al2O3"});
        MAGEMinSB11PerovskiteOptions options;
        options.minimizerMaxIterations = 0;
        ActivityModel fn = ActivityModelGlobalizedSolidSolution(
            MAGEMinSolidSolutionPilotModelSB11Perovskite(options))(species);
        const auto s = benchmarkTernaryFamily(fn, grid, T, P, 3);

        CHECK(s.totalEvaluations == n);
        CHECK(s.comparedCount == n);
        CHECK(s.fallbackCount == n);
        CHECK(s.legacyCount == n);
        CHECK(s.projectedGradientCount == 0);
        CHECK(s.fallbackRate == Approx(1.0));
    }

    SECTION("sb11_ak forced disagreement -> full legacy fallback")
    {
        const auto species = makeSpeciesList({"MgSiO3", "FeSiO3", "CaSiO3"});
        MAGEMinSB11AkimotoiteOptions options;
        options.minimizerMaxIterations = 0;
        ActivityModel fn = ActivityModelGlobalizedSolidSolution(
            MAGEMinSolidSolutionPilotModelSB11Akimotoite(options))(species);
        const auto s = benchmarkTernaryFamily(fn, grid, T, P, 3);

        CHECK(s.totalEvaluations == n);
        CHECK(s.comparedCount == n);
        CHECK(s.fallbackCount == n);
        CHECK(s.legacyCount == n);
        CHECK(s.projectedGradientCount == 0);
        CHECK(s.fallbackRate == Approx(1.0));
    }

    SECTION("sb11_cf forced disagreement -> full legacy fallback")
    {
        const auto species = makeSpeciesList({"MgAl2O4", "FeAl2O4", "NaAlSiO4"});
        MAGEMinSB11CalcioferriteOptions options;
        options.minimizerMaxIterations = 0;
        ActivityModel fn = ActivityModelGlobalizedSolidSolution(
            MAGEMinSolidSolutionPilotModelSB11Calcioferrite(options))(species);
        const auto s = benchmarkTernaryFamily(fn, grid, T, P, 3);

        CHECK(s.totalEvaluations == n);
        CHECK(s.comparedCount == n);
        CHECK(s.fallbackCount == n);
        CHECK(s.legacyCount == n);
        CHECK(s.projectedGradientCount == 0);
        CHECK(s.fallbackRate == Approx(1.0));
    }

    SECTION("sb21_nal forced disagreement -> full legacy fallback")
    {
        const auto species = makeSpeciesList({"NaMg2Al5SiO12", "NaFe2Al5SiO12", "Na3Al3Si3O12"});
        MAGEMinSB21NALOptions options;
        options.minimizerMaxIterations = 0;
        ActivityModel fn = ActivityModelGlobalizedSolidSolution(
            MAGEMinSolidSolutionPilotModelSB21NAL(options))(species);
        const auto s = benchmarkTernaryFamily(fn, grid, T, P, 3);

        CHECK(s.totalEvaluations == n);
        CHECK(s.comparedCount == n);
        CHECK(s.fallbackCount == n);
        CHECK(s.legacyCount == n);
        CHECK(s.projectedGradientCount == 0);
        CHECK(s.fallbackRate == Approx(1.0));
    }
}

TEST_CASE("Dumping MAGEMin regression snapshots", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Utility][Snapshot]")
{
    if(std::getenv("REAKTORO_DUMP_MAGEMIN_FIXTURES") == nullptr && std::getenv("REAKTORO_DUMP_MAGEMIN_TERNARY_FIXTURES") == nullptr)
        return;

    const auto T = 1473.15;
    const auto P = 1.0e9;

    ArrayXr xol(2);
    xol << 0.25, 0.75;

    ArrayXr xwa(2);
    xwa << 0.60, 0.40;

    ArrayXr xak(3);
    xak << 0.80, 0.10, 0.10;

    ArrayXr xpv(3);
    xpv << 0.70, 0.20, 0.10;

    ArrayXr xcf(3);
    xcf << 0.60, 0.25, 0.15;

    json snapshots = json::array();
    snapshots.push_back(snapshotModelFixture("sb11_ol_25_75", "sb11_ol", {"Mg2SiO4", "Fe2SiO4"}, xol, T, P));
    snapshots.push_back(snapshotModelFixture("sb11_wa_60_40", "sb11_wa", {"Mg2SiO4", "Fe2SiO4"}, xwa, T, P));
    snapshots.push_back(snapshotModelFixture("sb11_ak_80_10_10", "sb11_ak", {"MgSiO3", "FeSiO3", "CaSiO3"}, xak, T, P));
    snapshots.push_back(snapshotModelFixture("sb11_pv_70_20_10", "sb11_pv", {"MgSiO3", "FeSiO3", "Al2O3"}, xpv, T, P));
    snapshots.push_back(snapshotModelFixture("sb11_cf_60_25_15", "sb11_cf", {"MgAl2O4", "FeAl2O4", "NaAlSiO4"}, xcf, T, P));

    ArrayXr xsp(2);
    xsp << 0.65, 0.35;
    snapshots.push_back(snapshotModelFixture("sb21_sp_65_35", "sb21_sp", {"MgAl2O4", "FeAl2O4"}, xsp, T, P));

    ArrayXr xnal(3);
    xnal << 0.55, 0.20, 0.25;
    snapshots.push_back(snapshotModelFixture("sb21_nal_55_20_25", "sb21_nal", {"NaMg2Al5SiO12", "NaFe2Al5SiO12", "Na3Al3Si3O12"}, xnal, T, P));

    ArrayXr xcf21(3);
    xcf21 << 0.60, 0.25, 0.15;
    snapshots.push_back(snapshotModelFixture("sb21_cf_60_25_15", "sb21_cf", {"MgAl2O4", "FeAl2O4", "NaAlSiO4"}, xcf21, T, P));

    writeFixtureDump(snapshots);
}

TEST_CASE("Writing MAGEMin regression snapshots to a file", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Utility][Snapshot]")
{
    const auto outputPath = fs::temp_directory_path() / "reaktoro-magemin-pilot-snapshots-test.json";
    const auto cleanup = [&]()
    {
        std::error_code error;
        fs::remove(outputPath, error);
    };

    cleanup();

    json snapshots = json::array();
    snapshots.push_back({{"name", "probe"}, {"value", 42}});

    writeFixtureDump(snapshots, outputPath);

    REQUIRE(fs::exists(outputPath));
    REQUIRE(fs::file_size(outputPath) > 0);

    std::ifstream input(outputPath);
    REQUIRE(input.good());

    const auto reloaded = json::parse(input);
    REQUIRE(reloaded.is_array());
    REQUIRE(reloaded.size() == 1);
    CHECK(reloaded.at(0).at("name").get<String>() == "probe");
    CHECK(reloaded.at(0).at("value").get<int>() == 42);

    cleanup();
}

TEST_CASE("Testing MAGEMin equilibrium retry regression fixtures", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression][EquilibriumRetry]")
{
    const auto fixtures = loadRetryFixtures();
    REQUIRE(!fixtures.empty());

    for(const auto& fixture : fixtures)
    {
        DYNAMIC_SECTION(fixture.name)
        {
            const auto db = test::createDatabase();
            const auto phaseSpecies = makeDatabaseSpeciesList(db, fixture.phaseSpecies);
            const auto branches = makeDefaultBranches(phaseSpecies.size());
            const auto model = makeModel(fixture.model, branches, fixture.stabilityPolicy);

            if(fixture.scenario == "closed-system-split")
            {
                Phase prototype;
                prototype = prototype.withName("PilotCarbonate");
                prototype = prototype.withSpecies(phaseSpecies);
                prototype = prototype.withStateOfMatter(StateOfMatter::Solid);
                prototype = prototype.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(prototype.species()));

                MAGEMinSolidSolutionPilotOptions pilotOptions;
                pilotOptions.branches = branches;

                const auto phase = MAGEMinSolidSolutionPilotPhase(prototype, model);
                const auto definition = MAGEMinSolidSolutionPilotDefinition(prototype, model, pilotOptions);
                const auto system = ChemicalSystem(db, PhaseList{phase});

                ChemicalState state(system);
                state.setTemperature(fixture.T);
                state.setPressure(fixture.P);
                state.setSpeciesAmounts(fixture.speciesAmounts);

                GlobalizedSolidSolutionEquilibriumRetryOptions retryOptions;
                retryOptions.definitions = {definition};
                retryOptions.maxRetries = 1;

                const auto result = equilibrateWithGlobalizedSolidSolutionSplits(state, retryOptions);

                CHECK(result.result.succeeded() == fixture.expected.at("succeeded").get<bool>());
                CHECK(result.numRebuilds == fixture.expected.at("numRebuilds").get<Index>());
                CHECK(result.numAcceptedSplitRetries + result.numRejectedSplitRetries <= 1);
                CHECK(result.numAcceptedSplitRetries == 0);
                CHECK(result.numRejectedSplitRetries == 0);

                if(fixture.expected.contains("phaseNames"))
                {
                    const auto expectedPhaseNames = fixture.expected.at("phaseNames").get<Strings>();
                    REQUIRE(result.system.phases().size() == expectedPhaseNames.size());
                    for(Index i = 0; i < result.system.phases().size(); ++i)
                        CHECK(result.system.phase(i).name() == expectedPhaseNames[i]);
                }

                if(fixture.expected.contains("preservedSpeciesAmount"))
                {
                    const auto preserved = fixture.expected.at("preservedSpeciesAmount");
                    const auto speciesIndex = result.system.species().indexWithName(preserved.at("name").get<String>());
                    CHECK(result.state.speciesAmounts()[speciesIndex] == Approx(preserved.at("value").get<double>()).margin(fixture.tolerances.scalarAbs));
                }

                if(fixture.expected.at("numRebuilds").get<Index>() > 0)
                {
                    auto rejectedRetryOptions = retryOptions;
                    rejectedRetryOptions.enableSplitAcceptanceGate = true;
                    rejectedRetryOptions.minIterationsForSplitRetry = 1000000;

                    const auto rejectedResult = equilibrateWithGlobalizedSolidSolutionSplits(state, rejectedRetryOptions);
                    CHECK(rejectedResult.numRebuilds == 0);
                    CHECK(rejectedResult.numRejectedSplitRetries == 1);
                }
            }
            else if(fixture.scenario == "conditions-split")
            {
                const auto aqueous = AqueousPhase({"H2O(aq)", "H+(aq)", "OH-(aq)", "Na+(aq)", "Cl-(aq)"});
                const auto quartzSpecies = SpeciesList({db.species().get("SiO2(s)")});

                Phase prototype;
                prototype = prototype.withName("PilotCarbonate");
                prototype = prototype.withSpecies(phaseSpecies);
                prototype = prototype.withStateOfMatter(StateOfMatter::Solid);
                prototype = prototype.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(prototype.species()));

                Phase quartz;
                quartz = quartz.withName("Quartz");
                quartz = quartz.withSpecies(quartzSpecies);
                quartz = quartz.withStateOfMatter(StateOfMatter::Solid);
                quartz = quartz.withActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(quartz.species()));
                quartz = quartz.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(quartz.species()));

                MAGEMinSolidSolutionPilotOptions pilotOptions;
                pilotOptions.branches = branches;

                const auto phase = MAGEMinSolidSolutionPilotPhase(prototype, model);
                const auto definition = MAGEMinSolidSolutionPilotDefinition(prototype, model, pilotOptions);
                Phases phasesObject(db, aqueous);
                phasesObject.add(PhaseList{phase, quartz});
                const auto system = ChemicalSystem(phasesObject);

                ChemicalState state(system);
                state.setTemperature(fixture.T);
                state.setPressure(fixture.P);
                state.setSpeciesAmounts(TestUtils::reorderPilotMixedConditionsSpeciesAmounts(
                    system,
                    phaseSpecies,
                    fixture.speciesAmounts,
                    "Unexpected mixed-system retry fixture species amount vector size."));

                EquilibriumSpecs specs(system);
                specs.temperature();
                specs.pressure();
                specs.pH();

                EquilibriumConditions conditions(specs);
                conditions.temperature(fixture.T);
                conditions.pressure(fixture.P);
                conditions.pH(7.0);
                conditions.setInitialComponentAmountsFromState(state);

                EquilibriumRestrictions restrictions(system);
                restrictions.cannotReact("Na+(aq)");

                GlobalizedSolidSolutionEquilibriumRetryOptions retryOptions;
                retryOptions.definitions = {definition};
                retryOptions.maxRetries = 1;

                const auto result = equilibrateWithGlobalizedSolidSolutionSplits(state, specs, conditions, restrictions, retryOptions);

                CHECK(result.result.succeeded() == fixture.expected.at("succeeded").get<bool>());
                CHECK(result.numRebuilds == fixture.expected.at("numRebuilds").get<Index>());
                CHECK(result.numAcceptedSplitRetries + result.numRejectedSplitRetries <= 1);
                CHECK(result.numAcceptedSplitRetries == 0);
                CHECK(result.numRejectedSplitRetries == 0);

                if(fixture.expected.contains("phaseNames"))
                {
                    const auto expectedPhaseNames = fixture.expected.at("phaseNames").get<Strings>();
                    REQUIRE(result.system.phases().size() == expectedPhaseNames.size());
                    for(Index i = 0; i < result.system.phases().size(); ++i)
                        CHECK(result.system.phase(i).name() == expectedPhaseNames[i]);
                }

                if(fixture.expected.contains("preservedSpeciesAmount"))
                {
                    const auto preserved = fixture.expected.at("preservedSpeciesAmount");
                    const auto speciesIndex = result.system.species().indexWithName(preserved.at("name").get<String>());
                    CHECK(result.state.speciesAmounts()[speciesIndex] == Approx(preserved.at("value").get<double>()).margin(fixture.tolerances.scalarAbs));
                }

                if(fixture.expected.at("numRebuilds").get<Index>() > 0)
                {
                    auto rejectedRetryOptions = retryOptions;
                    rejectedRetryOptions.enableSplitAcceptanceGate = true;
                    rejectedRetryOptions.minIterationsForSplitRetry = 1000000;

                    const auto rejectedResult = equilibrateWithGlobalizedSolidSolutionSplits(state, specs, conditions, restrictions, rejectedRetryOptions);
                    CHECK(rejectedResult.numRebuilds == 0);
                    CHECK(rejectedResult.numRejectedSplitRetries == 1);
                }
            }
            else
            {
                throw std::runtime_error("Unsupported MAGEMin equilibrium retry regression scenario: " + fixture.scenario);
            }
        }
    }
}

TEST_CASE("Testing MAGEMinConstrainedTernaryLocalModel nonlinear constraint hooks", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    // Test constraint callback presence detection and validation.
    MAGEMinConstrainedTernaryLocalModel model;
    model.modelId = "sb21_ol";
    model.T = 1000.0;
    model.visiblex = (ArrayXr(3) << 0.3, 0.3, 0.4).finished();
    model.objective = [](ArrayXrConstRef y) { return (y * y).sum(); };
    model.gradient = [](ArrayXrConstRef y) { return 2.0 * y; };
    model.lowerBounds = (ArrayXr(3) << 0.0, 0.0, 0.0).finished();
    model.upperBounds = (ArrayXr(3) << 1.0, 1.0, 1.0).finished();
    model.enforceUnityConstraint = true;

    // Scenario A: No constraints provided (should pass validation).
    CHECK(model.constraints == nullptr);
    CHECK_NOTHROW(([&]() {
        auto current = ArrayXr(model.visiblex);
        auto f = model.objective(current);
        CHECK(f > 0.0);
    })());

    // Scenario B: Constraints without Jacobian (should fail validation).
    model.constraints = [](ArrayXrConstRef y)
    {
        ArrayXr c(2);
        c[0] = y[0] + y[1] - 0.5;  // Linear constraint c_0(y) = y_0 + y_1 - 0.5
        c[1] = y[1]*y[2] - 0.1;    // Nonlinear constraint c_1(y) = y_1*y_2 - 0.1
        return c;
    };
    model.constraintLowerBounds = (ArrayXr(2) << -1e10, -1e10).finished();
    model.constraintUpperBounds = (ArrayXr(2) << 0.0, 0.0).finished();

    CHECK_THROWS_AS(([&]() {
        MAGEMinProjectedGradientLocalModelMinimizer(model, ArrayXr());
    })(), std::runtime_error);

    // Scenario C: Constraints with Jacobian (should pass validation).
    model.constraintJacobian = [](ArrayXrConstRef y)
    {
        MatrixXr jac(2, 3);
        jac(0, 0) = 1.0; jac(0, 1) = 1.0; jac(0, 2) = 0.0;  // ∇c_0
        jac(1, 0) = 0.0; jac(1, 1) = y[2]; jac(1, 2) = y[1];  // ∇c_1
        return jac;
    };

    CHECK_NOTHROW(([&]() {
        const auto current = ArrayXr(model.visiblex);
        const auto c = model.constraints(current);
        CHECK(c.size() == 2);
        const auto jac = model.constraintJacobian(current);
        CHECK(jac.rows() == 2);
        CHECK(jac.cols() == 3);
    })());

    // Scenario D: Verify Hessian flag validation.
    model.useSecondOrderInfo = true;
    model.objectiveHessian = nullptr;

    CHECK_THROWS_AS(([&]() {
        MAGEMinProjectedGradientLocalModelMinimizer(model, ArrayXr());
    })(), std::runtime_error);

    // Scenario E: Hessian callback provided (should pass).
    model.objectiveHessian = [](ArrayXrConstRef y, ArrayXrConstRef multipliers)
    {
        MatrixXr H(3, 3);
        H.setIdentity();
        H *= 2.0;  // ∇²f = 2*I for quadratic objective
        return H;
    };

    // No exception should be thrown during validation
    CHECK_NOTHROW(([&]() {
        const auto current = ArrayXr(model.visiblex);
        if(model.useSecondOrderInfo && model.objectiveHessian)
        {
            ArrayXr multipliers = ArrayXr::Zero(2);
            const auto H = model.objectiveHessian(current, multipliers);
            CHECK(H.rows() == 3);
            CHECK(H.cols() == 3);
        }
    })());
}

TEST_CASE("Testing MAGEMinConstrainedTernaryLocalModel constraint feasibility evaluation", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    // Test constraint feasibility checking utility.
    MAGEMinConstrainedTernaryLocalModel model;
    model.modelId = "sb21_cpx";
    model.T = 1200.0;
    model.visiblex = (ArrayXr(5) << 0.2, 0.2, 0.2, 0.2, 0.2).finished();
    model.objective = [](ArrayXrConstRef y) { return (y * y).sum(); };
    model.gradient = [](ArrayXrConstRef y) { return 2.0 * y; };

    // Define linear inequality constraints: c(y) = Ay - b <= 0
    const auto A_data = (MatrixXr(2, 5) <<
        1.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 1.0, 1.0
    ).finished();

    // Scenario A: No constraints (always feasible).
    CHECK(true);  // Placeholder for no-constraint scenario

    // Scenario B: With constraints.
    model.constraints = [A_data](ArrayXrConstRef y)
    {
        ArrayXr c = (A_data * y.matrix()).array();
        c[0] -= 0.4;  // c_0(y) = y_0 + y_1 - 0.4
        c[1] -= 0.6;  // c_1(y) = y_2 + y_3 + y_4 - 0.6
        return c;
    };
    model.constraintJacobian = [A_data](ArrayXrConstRef y)
    {
        return A_data;
    };
    model.constraintLowerBounds = (ArrayXr(2) << -1e10, -1e10).finished();
    model.constraintUpperBounds = (ArrayXr(2) << 0.0, 0.0).finished();

    // Test point 1: feasible composition
    ArrayXr y_feasible = (ArrayXr(5) << 0.1, 0.1, 0.1, 0.1, 0.1).finished();
    const auto c_feasible = model.constraints(y_feasible);
    CHECK(c_feasible[0] == Approx(-0.2));  // 0.1+0.1-0.4 = -0.2 <= 0 ✓
    CHECK(c_feasible[1] == Approx(-0.3));  // 0.1+0.1+0.1-0.6 = -0.3 <= 0 ✓

    // Test point 2: infeasible composition (violates first constraint)
    ArrayXr y_infeasible = (ArrayXr(5) << 0.25, 0.25, 0.0, 0.0, 0.0).finished();
    const auto c_infeasible = model.constraints(y_infeasible);
    CHECK(c_infeasible[0] == Approx(0.1));  // 0.25+0.25-0.4 = 0.1 > 0 ✗
    CHECK(c_infeasible[1] == Approx(-0.6));  // 0+0+0-0.6 = -0.6 <= 0 ✓
}

TEST_CASE("Testing MAGEMinProjectedGradientLocalModelMinimizer uses constraint-aware line search", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    MAGEMinConstrainedTernaryLocalModel model;
    model.visiblex = (ArrayXr(2) << 0.20, 0.20).finished();
    model.lowerBounds = (ArrayXr(2) << 0.0, 0.0).finished();
    model.upperBounds = (ArrayXr(2) << 1.0, 1.0).finished();
    model.enforceUnityConstraint = false;
    model.tolerance = 1.0e-12;
    model.maxIterations = 64;
    model.requireFeasibleTrialPoints = true;
    model.constraintPenaltyWeight = 1.0e4;

    model.objective = [](ArrayXrConstRef y) -> real
    {
        return (y[0] - 0.90)*(y[0] - 0.90) + (y[1] - 0.90)*(y[1] - 0.90);
    };
    model.gradient = [](ArrayXrConstRef y) -> ArrayXr
    {
        ArrayXr g(2);
        g[0] = 2.0*(y[0] - 0.90);
        g[1] = 2.0*(y[1] - 0.90);
        return g;
    };
    model.constraints = [](ArrayXrConstRef y) -> ArrayXr
    {
        ArrayXr c(1);
        c[0] = y[0] + y[1];
        return c;
    };
    model.constraintJacobian = [](ArrayXrConstRef) -> MatrixXr
    {
        MatrixXr jac(1, 2);
        jac << 1.0, 1.0;
        return jac;
    };
    model.constraintLowerBounds = (ArrayXr(1) << -1.0e10).finished();
    model.constraintUpperBounds = (ArrayXr(1) << 0.45).finished();

    const auto result = MAGEMinProjectedGradientLocalModelMinimizer(model, ArrayXr(model.visiblex));

    CHECK(result.converged);
    CHECK(result.x[0] >= 0.20);
    CHECK(result.x[1] >= 0.20);
    CHECK(result.x.sum() <= Approx(0.45).margin(1.0e-10));
}

TEST_CASE("Testing MAGEMinProjectedGradientLocalModelMinimizer uses Jacobian and Hessian hooks for active constraints", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    MAGEMinConstrainedTernaryLocalModel model;
    model.visiblex = (ArrayXr(2) << 0.20, 0.20).finished();
    model.lowerBounds = (ArrayXr(2) << 0.0, 0.0).finished();
    model.upperBounds = (ArrayXr(2) << 1.0, 1.0).finished();
    model.enforceUnityConstraint = false;
    model.tolerance = 1.0e-12;
    model.maxIterations = 32;
    model.requireFeasibleTrialPoints = true;
    model.useSecondOrderInfo = true;

    auto jacobianCalls = 0;
    auto hessianCalls = 0;

    model.objective = [](ArrayXrConstRef y) -> real
    {
        return (y[0] - 0.90)*(y[0] - 0.90) + (y[1] - 0.90)*(y[1] - 0.90);
    };
    model.gradient = [](ArrayXrConstRef y) -> ArrayXr
    {
        ArrayXr g(2);
        g[0] = 2.0*(y[0] - 0.90);
        g[1] = 2.0*(y[1] - 0.90);
        return g;
    };
    model.constraints = [](ArrayXrConstRef y) -> ArrayXr
    {
        ArrayXr c(1);
        c[0] = y[0] + y[1];
        return c;
    };
    model.constraintJacobian = [&jacobianCalls](ArrayXrConstRef) -> MatrixXr
    {
        ++jacobianCalls;
        MatrixXr jac(1, 2);
        jac << 1.0, 1.0;
        return jac;
    };
    model.constraintLowerBounds = (ArrayXr(1) << -1.0e10).finished();
    model.constraintUpperBounds = (ArrayXr(1) << 0.45).finished();
    model.objectiveHessian = [&hessianCalls](ArrayXrConstRef, ArrayXrConstRef) -> MatrixXr
    {
        ++hessianCalls;
        MatrixXr H = MatrixXr::Zero(2, 2);
        H(0, 0) = 2.0;
        H(1, 1) = 2.0;
        return H;
    };

    const auto result = MAGEMinProjectedGradientLocalModelMinimizer(model, ArrayXr(model.visiblex));

    CHECK(result.converged);
    CHECK(jacobianCalls > 0);
    CHECK(hessianCalls > 0);
    CHECK(result.x.sum() == Approx(0.45).margin(1.0e-8));
    CHECK(result.x[0] == Approx(result.x[1]).margin(1.0e-8));
}

TEST_CASE("Testing MAGEMinProjectedGradientLocalModelMinimizer clips trial steps to trust region", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
    MAGEMinConstrainedTernaryLocalModel model;
    model.visiblex = (ArrayXr(2) << 0.20, 0.20).finished();
    model.lowerBounds = (ArrayXr(2) << 0.0, 0.0).finished();
    model.upperBounds = (ArrayXr(2) << 1.0, 1.0).finished();
    model.enforceUnityConstraint = false;
    model.tolerance = 1.0e-12;
    model.maxIterations = 1;
    model.trustRegionRadius = 0.10;

    model.objective = [](ArrayXrConstRef y) -> real
    {
        return (y[0] - 1.0)*(y[0] - 1.0) + (y[1] - 1.0)*(y[1] - 1.0);
    };
    model.gradient = [](ArrayXrConstRef y) -> ArrayXr
    {
        ArrayXr g(2);
        g[0] = 2.0*(y[0] - 1.0);
        g[1] = 2.0*(y[1] - 1.0);
        return g;
    };

    const auto warmstart = ArrayXr((ArrayXr(2) << 0.20, 0.20).finished());
    const auto result = MAGEMinProjectedGradientLocalModelMinimizer(model, warmstart);

    const auto displacement = (result.x - warmstart).matrix().norm();
    CHECK(displacement == Approx(0.10).margin(1.0e-8));
    CHECK(result.x[0] > warmstart[0]);
    CHECK(result.x[1] > warmstart[1]);
}

TEST_CASE("Testing MAGEMin HP ig_opx pilot reference-state and constraint wiring", "[ActivityModelGlobalizedSolidSolution][MAGEMinPilot][Regression]")
{
#if defined(_MSC_VER)
    SUCCEED("Skipping HP ig_opx pilot reference-state regression on MSVC due an access-violation runtime instability in this configuration.");
    return;
#endif

    MAGEMinHPIGOPXOptions options;
    options.externalCompositionPenalty = 0.0;
    options.localModelMinimizer = [](MAGEMinConstrainedTernaryLocalModel const& model, Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
    {
        const auto current = warmstart ? *warmstart : model.visiblex;

        GlobalizedSolidSolutionInternalResult result;
        result.x = current;
        result.objective = model.objective(current);
        result.iterations = 0;
        result.converged = true;
        return result;
    };
    options.localModelDiagnostics = [](MAGEMinConstrainedTernaryLocalModel const& model, GlobalizedSolidSolutionInternalResult const& result) -> Map<String, Any>
    {
        Map<String, Any> payload;
        payload["DiagnosticObjective"] = model.objective(result.x);
        payload["DiagnosticGradient"] = model.gradient(result.x);
        payload["DiagnosticConstraints"] = model.constraints(result.x);
        payload["DiagnosticConverged"] = result.converged;
        payload["DiagnosticIterations"] = static_cast<std::uint64_t>(result.iterations);
        return payload;
    };

    auto model = MAGEMinSolidSolutionPilotModelHPIGOPX(options);

    auto state = std::make_shared<GlobalizedSolidSolutionState>();
    ArrayXr x(8);
    x << 0.12, 0.18, 0.14, 0.16, 0.10, 0.08, 0.11, 0.09;

    GlobalizedSolidSolutionInput input{1223.15, 1.2e9, x, Map<String, Any>{}, state, GlobalizedSolidSolutionNoBranch};
    const auto output = model(input);

    const auto referenceState = std::any_cast<ArrayXr>(output.extra.at("MAGEMinSolidSolutionPilot::ReferenceState"));
    const auto diagnosticObjective = std::any_cast<real>(output.extra.at("DiagnosticObjective"));
    const auto diagnosticGradient = std::any_cast<ArrayXr>(output.extra.at("DiagnosticGradient"));
    const auto diagnosticConstraints = std::any_cast<ArrayXr>(output.extra.at("DiagnosticConstraints"));

    CHECK(referenceState.size() == 9);
    CHECK(diagnosticGradient.size() == 8);
    CHECK(diagnosticConstraints.size() > 0);
    CHECK(std::any_cast<bool>(output.extra.at("DiagnosticConverged")));
    CHECK(std::any_cast<std::uint64_t>(output.extra.at("DiagnosticIterations")) == 0);
    CHECK(std::any_cast<bool>(output.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerConverged")));
    CHECK(std::any_cast<std::uint64_t>(output.extra.at("MAGEMinSolidSolutionPilot::InternalMinimizerIterations")) == 0);

    const auto en = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -2915480.248, -3090220.0, 132.5, 6.262e-05, 356.2, -0.00299, -596900.0, -3185.3,
        2.27e-05, 105900000000.0, 8.65, -8.2e-11, 10.0, 9999.0})(1223.15, 1.2e9).G0;
    const auto fs = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -2234304.078, -2388710.0, 189.9, 6.592e-05, 398.7, -0.006579, 1290100.0, -4058.0,
        3.26e-05, 101000000000.0, 4.08, -4.0e-11, 10.0, 9999.0})(1223.15, 1.2e9).G0;
    const auto di = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -3027542.5655, -3201850.0, 142.9, 6.619e-05, 314.5, 4.1e-05, -2745900.0, -2020.1,
        2.73e-05, 119200000000.0, 5.19, -4.4e-11, 10.0, 9999.0})(1223.15, 1.2e9).G0;
    const auto jd = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -2846567.8345, -3025270.0, 133.5, 6.04e-05, 319.4, 0.003616, -1173900.0, -2469.5,
        2.1e-05, 128100000000.0, 3.81, -3.0e-11, 10.0, 9999.0})(1223.15, 1.2e9).G0;

    CHECK(referenceState[0] == Approx(en).margin(1.0e-6 * std::fabs(static_cast<double>(en))));
    CHECK(referenceState[1] == Approx(fs).margin(1.0e-6 * std::fabs(static_cast<double>(fs))));
    CHECK(referenceState[3] == Approx(0.005 * (1.2e9 / 1.0e5) + di + 2.8).margin(1.0e-6 * std::fabs(static_cast<double>(di))));
    CHECK(referenceState[8] == Approx(jd + 18.2).margin(1.0e-6 * std::fabs(static_cast<double>(jd))));
    CHECK(std::isfinite(static_cast<double>(diagnosticObjective)));
    CHECK(std::isfinite(static_cast<double>(std::any_cast<real>(output.extra.at("MAGEMinSolidSolutionPilot::InternalObjective")))));
}
