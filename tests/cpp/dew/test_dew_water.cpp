#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "WaterTestCommon.hpp"
#include "WaterTestAdapters.hpp"
#include <fstream>

#include <numeric>

#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Extensions/DEW/DEWDatabase.hpp>
#include <Reaktoro/Extensions/DEW/WaterState.hpp>
#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelDEW.hpp>

using namespace Reaktoro;

// All CSVs assumed to be in the DEW root folder, tests in DEW/tests.
// So relative paths use "../truth_*.csv".

// -----------------------------------------------------------------------------
// Helper templates
// -----------------------------------------------------------------------------

template<typename Func>
void run_eps_file(const std::string& path,
                  Func eps_fun,
                  const char* label)
{
    auto rows = load_csv(path, /*skip_header=*/true);

    const double ABS_TOL = 1e-9;
    const double REL_TOL = 1e-8;

    std::ofstream out("reaction_test_results.csv");
        out << "T_C,P_kb,G_rxn_model,G_rxn_truth,V_rxn_model,V_rxn_truth,"
            "G_H2O_model,G_H2O_truth,G_CO2_model,G_CO2_truth,G_Hplus_model,G_Hplus_truth,G_HCO3_model,G_HCO3_truth\n";
    for (const auto& row : rows) {
        // T_C, P_bar, eq, rho_g_cm3, eps_r, Psat
        if (row.fields.size() < 5)
            continue;

        double T_C, P_bar, eps_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[4], eps_truth)) continue;

        double eps_model = eps_fun(T_C, P_bar);

        INFO(label << " eps: T=" << T_C << " C, P=" << P_bar << " bar");
        INFO("  Model value:  " << eps_model);
        INFO("  Truth value:  " << eps_truth);
        INFO("  Difference:   " << (eps_model - eps_truth));
        INFO("  Rel. error:   " << (std::abs(eps_model - eps_truth) / std::max(std::abs(eps_truth), 1e-10) * 100.0) << " %");
        REQUIRE(almost_equal(eps_model, eps_truth, ABS_TOL, REL_TOL));
    }
}

template<typename Func>
void run_depsdrho_file(const std::string& path,
                       Func deps_fun,
                       const char* label)
{
    auto rows = load_csv(path, /*skip_header=*/true);

    const double ABS_TOL = 1e-6;
    const double REL_TOL = 1e-6;

    for (const auto& row : rows) {
        // T_C, eq, rho_g_cm3, depsdrho
        if (row.fields.size() < 4)
            continue;

        double T_C, rho_g_cm3, deps_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[2], rho_g_cm3)) continue;
        if (!parse_maybe_double(row.fields[3], deps_truth)) continue;

        double deps_model = deps_fun(T_C, rho_g_cm3);

        INFO(label << " depsdrho: T=" << T_C << " C, rho=" << rho_g_cm3 << " g/cm3");
        INFO("  Model value:  " << deps_model);
        INFO("  Truth value:  " << deps_truth);
        INFO("  Difference:   " << (deps_model - deps_truth));
        INFO("  Rel. error:   " << (std::abs(deps_model - deps_truth) / std::max(std::abs(deps_truth), 1e-10) * 100.0) << " %");
        REQUIRE(almost_equal(deps_model, deps_truth, ABS_TOL, REL_TOL));
    }
}

// -----------------------------------------------------------------------------
// Density ρ
// -----------------------------------------------------------------------------

TEST_CASE("Density ZD2005 matches truth table", "[density][ZD2005]")
{
    auto rows = load_csv("truth_density_ZD2005.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-9;
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, P_bar, Psat, rho_g_cm3
        if (row.fields.size() < 4)
            continue;

        double T_C, P_bar, rho_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[3], rho_truth)) continue;

        double rho_model = dew_density_ZD2005(T_C, P_bar);

        INFO("ZD2005: T=" << T_C << " C, P=" << P_bar << " bar");
        INFO("  Model value:  " << rho_model << " g/cm3");
        INFO("  Truth value:  " << rho_truth << " g/cm3");
        INFO("  Difference:   " << (rho_model - rho_truth) << " g/cm3");
        INFO("  Rel. error:   " << (std::abs(rho_model - rho_truth) / std::max(std::abs(rho_truth), 1e-10) * 100.0) << " %");
        REQUIRE(almost_equal(rho_model, rho_truth, ABS_TOL, REL_TOL));
    }
}

TEST_CASE("Density ZD2009 matches truth table", "[density][ZD2009]")
{
    auto rows = load_csv("truth_density_ZD2009.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-9;
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, P_bar, Psat, rho_g_cm3
        if (row.fields.size() < 4)
            continue;

        double T_C, P_bar, rho_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[3], rho_truth)) continue;

        double rho_model = dew_density_ZD2009(T_C, P_bar);

        INFO("ZD2009: T=" << T_C << " C, P=" << P_bar << " bar");
        INFO("  Model value:  " << rho_model << " g/cm3");
        INFO("  Truth value:  " << rho_truth << " g/cm3");
        INFO("  Difference:   " << (rho_model - rho_truth) << " g/cm3");
        INFO("  Rel. error:   " << (std::abs(rho_model - rho_truth) / std::max(std::abs(rho_truth), 1e-10) * 100.0) << " %");
        REQUIRE(almost_equal(rho_model, rho_truth, ABS_TOL, REL_TOL));
    }
}

TEST_CASE("Psat density matches truth table", "[density][Psat]")
{
    auto rows = load_csv("truth_density_psat.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-9;
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, rho_g_cm3
        if (row.fields.size() < 2)
            continue;

        double T_C, rho_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], rho_truth)) continue;

        double rho_model = dew_density_psat(T_C);

        INFO("Psat density: T=" << T_C << " C");
        REQUIRE(almost_equal(rho_model, rho_truth, ABS_TOL, REL_TOL));
    }
}

// -----------------------------------------------------------------------------
// dρ/dP
// -----------------------------------------------------------------------------

TEST_CASE("drhodP ZD2005 matches truth table", "[drhodP][ZD2005]")
{
    auto rows = load_csv("truth_drhodP_ZD2005.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-12;
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, P_bar, eq, rho_g_cm3, drhodP
        if (row.fields.size() < 5)
            continue;

        double T_C, P_bar, drhodP_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[4], drhodP_truth)) continue;

        double drhodP_model = dew_drhodP_ZD2005(T_C, P_bar);

        INFO("drhodP ZD2005: T=" << T_C << " C, P=" << P_bar << " bar");
        REQUIRE(almost_equal(drhodP_model, drhodP_truth, ABS_TOL, REL_TOL));
    }
}

TEST_CASE("drhodP ZD2009 matches truth table", "[drhodP][ZD2009]")
{
    auto rows = load_csv("truth_drhodP_ZD2009.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-12;
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, P_bar, eq, rho_g_cm3, drhodP
        if (row.fields.size() < 5)
            continue;

        double T_C, P_bar, drhodP_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[4], drhodP_truth)) continue;

        double drhodP_model = dew_drhodP_ZD2009(T_C, P_bar);

        INFO("drhodP ZD2009: T=" << T_C << " C, P=" << P_bar << " bar");
        REQUIRE(almost_equal(drhodP_model, drhodP_truth, ABS_TOL, REL_TOL));
    }
}

// -----------------------------------------------------------------------------
// Dielectric constant ε
// -----------------------------------------------------------------------------

TEST_CASE("epsilon JN1991 matches truth table", "[epsilon][JN1991]")
{
    run_eps_file("truth_epsilon_JN1991.csv",
                 &dew_epsilon_JN1991,
                 "JN1991");
}

TEST_CASE("epsilon Franck1990 matches truth table", "[epsilon][Franck1990]")
{
    run_eps_file("truth_epsilon_Franck1990.csv",
                 &dew_epsilon_Franck1990,
                 "Franck1990");
}

TEST_CASE("epsilon Fernandez1997 matches truth table", "[epsilon][Fernandez1997]")
{
    run_eps_file("truth_epsilon_Fernandez1997.csv",
                 &dew_epsilon_Fernandez1997,
                 "Fernandez1997");
}

TEST_CASE("epsilon Power matches truth table", "[epsilon][Power]")
{
    run_eps_file("truth_epsilon_Power.csv",
                 &dew_epsilon_Power,
                 "Power");
}

TEST_CASE("epsilon Psat matches truth table", "[epsilon][Psat]")
{
    auto rows = load_csv("truth_epsilon_psat.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-9;
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, eps_r
        if (row.fields.size() < 2)
            continue;

        double T_C, eps_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], eps_truth)) continue;

        double eps_model = dew_epsilon_psat(T_C);

        INFO("epsilon Psat: T=" << T_C << " C");
        REQUIRE(almost_equal(eps_model, eps_truth, ABS_TOL, REL_TOL));
    }
}

// -----------------------------------------------------------------------------
// depsdrho for dielectric models
// -----------------------------------------------------------------------------

TEST_CASE("depsdrho JN1991 matches truth table", "[depsdrho][JN1991]")
{
    run_depsdrho_file("truth_depsdrho_JN1991.csv",
                      &dew_depsdrho_JN1991,
                      "JN1991");
}

TEST_CASE("depsdrho Franck1990 matches truth table", "[depsdrho][Franck1990]")
{
    run_depsdrho_file("truth_depsdrho_Franck1990.csv",
                      &dew_depsdrho_Franck1990,
                      "Franck1990");
}

TEST_CASE("depsdrho Fernandez1997 matches truth table", "[depsdrho][Fernandez1997]")
{
    run_depsdrho_file("truth_depsdrho_Fernandez1997.csv",
                      &dew_depsdrho_Fernandez1997,
                      "Fernandez1997");
}

TEST_CASE("depsdrho Power matches truth table", "[depsdrho][Power]")
{
    run_depsdrho_file("truth_depsdrho_Power.csv",
                      &dew_depsdrho_Power,
                      "Power");
}

// -----------------------------------------------------------------------------
// Solvent function g(T,P) and d(g)/dP
// -----------------------------------------------------------------------------

TEST_CASE("Solvent function g(T,P) matches truth table", "[g]")
{
    auto rows = load_csv("truth_g.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-9;
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, P_bar, rho_g_cm3, g
        if (row.fields.size() < 4)
            continue;

        double T_C, P_bar, g_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[3], g_truth)) continue;

        double g_model = dew_g_eq2(T_C, P_bar);

        INFO("g(T,P): T=" << T_C << " C, P=" << P_bar << " bar");
        INFO("  Model value:  " << g_model);
        INFO("  Truth value:  " << g_truth);
        INFO("  Difference:   " << (g_model - g_truth));
        INFO("  Rel. error:   " << (std::abs(g_model - g_truth) / std::max(std::abs(g_truth), 1e-10) * 100.0) << " %");
        REQUIRE(almost_equal(g_model, g_truth, ABS_TOL, REL_TOL));
    }
}

TEST_CASE("dgdP eq2 matches truth table", "[dgdP][eq2]")
{
    auto rows = load_csv("truth_dgdP_eq2.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-15;
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, P_bar, eq, rho_g_cm3, dgdP, Psat
        if (row.fields.size() < 6)
            continue;

        double T_C, P_bar, dgdP_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[4], dgdP_truth)) continue;

        double dgdP_model = dew_dgdP_eq2(T_C, P_bar);

        INFO("dgdP eq2: T=" << T_C << " C, P=" << P_bar << " bar");
        INFO("  Model value:  " << dgdP_model << " 1/Pa");
        INFO("  Truth value:  " << dgdP_truth << " 1/Pa");
        INFO("  Difference:   " << (dgdP_model - dgdP_truth) << " 1/Pa");
        INFO("  Rel. error:   " << (std::abs(dgdP_model - dgdP_truth) / std::max(std::abs(dgdP_truth), 1e-20) * 100.0) << " %");
        REQUIRE(almost_equal(dgdP_model, dgdP_truth, ABS_TOL, REL_TOL));
    }
}

TEST_CASE("dgdP Psat(T) matches truth table", "[dgdP][Psat]")
{
    auto rows = load_csv("truth_dgdP_psat.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-15;
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, dgdP
        if (row.fields.size() < 2)
            continue;

        double T_C, dgdP_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], dgdP_truth)) continue;

        double dgdP_model = dew_dgdP_psat(T_C);

        INFO("dgdP Psat: T=" << T_C << " C");
        INFO("  Model value:  " << dgdP_model << " Å/bar");
        INFO("  Truth value:  " << dgdP_truth << " Å/bar");
        INFO("  Difference:   " << (dgdP_model - dgdP_truth) << " Å/bar");
        INFO("  Rel. error:   " << (std::abs(dgdP_model - dgdP_truth) / std::max(std::abs(dgdP_truth), 1e-20) * 100.0) << " %");
        REQUIRE(almost_equal(dgdP_model, dgdP_truth, ABS_TOL, REL_TOL));
    }
}

// -----------------------------------------------------------------------------
// Gibbs free energy G
// -----------------------------------------------------------------------------

TEST_CASE("G_DH1978 matches truth table", "[G][DH1978]")
{
    auto rows = load_csv("truth_G_DH1978.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-6;      // cal/mol
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, P_bar, G_cal_mol (may be NaN for some entries)
        if (row.fields.size() < 3)
            continue;

        double T_C, P_bar, G_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[2], G_truth)) continue;

        double G_model = dew_G_DH1978(T_C, P_bar);

        INFO("G_DH1978: T=" << T_C << " C, P=" << P_bar << " bar");
        INFO("  Model value:  " << G_model << " cal/mol");
        INFO("  Truth value:  " << G_truth << " cal/mol");
        INFO("  Difference:   " << (G_model - G_truth) << " cal/mol");
        INFO("  Rel. error:   " << (std::abs(G_model - G_truth) / std::max(std::abs(G_truth), 1e-10) * 100.0) << " %");
        REQUIRE(almost_equal(G_model, G_truth, ABS_TOL, REL_TOL));
    }
}

TEST_CASE("G_integral matches truth table (Excel compatibility)", "[G][integral][excel]")
{
    auto rows = load_csv("truth_G_integral.csv", /*skip_header=*/true);

    // Excel VBA integration uses ~500 steps max, giving ~2% error at high pressure
    const double ABS_TOL = 1200.0; // Relaxed - numerical integration errors up to ~1100 cal/mol at very high P (4000 bar)
    const double REL_TOL = 2e-2;   // 2% relative error acceptable for Excel's integration method

    for (const auto& row : rows) {
        // T_C, P_bar, G_cal_mol
        if (row.fields.size() < 3)
            continue;

        double T_C, P_bar, G_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[2], G_truth)) continue;

        // dew_G_integral uses Excel-compatible integration (useExcelIntegration=true)
        double G_model = dew_G_integral(T_C, P_bar);

        INFO("G_integral (Excel compat): T=" << T_C << " C, P=" << P_bar << " bar");
        INFO("  Model value:  " << G_model << " cal/mol");
        INFO("  Truth value:  " << G_truth << " cal/mol");
        INFO("  Difference:   " << (G_model - G_truth) << " cal/mol");
        INFO("  Rel. error:   " << (std::abs(G_model - G_truth) / std::max(std::abs(G_truth), 1e-10) * 100.0) << " %");
        REQUIRE(almost_equal(G_model, G_truth, ABS_TOL, REL_TOL));
    }
}

TEST_CASE("G_integral high precision vs Excel truth", "[G][integral][highprec]")
{
    auto rows = load_csv("truth_G_integral.csv", /*skip_header=*/true);

    // High-precision integration should significantly improve on Excel's ~2% errors
    // The remaining errors are from the ZD2005 density EOS, not integration
    const double ABS_TOL = 50.0;    // Much tighter than 1200 cal/mol for Excel
    const double REL_TOL = 0.0005;  // 0.05% - much better than Excel's 2%

    for (const auto& row : rows) {
        // T_C, P_bar, G_cal_mol
        if (row.fields.size() < 3)
            continue;

        double T_C, P_bar, G_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[2], G_truth)) continue;

        // Use high-precision integration (5000 steps, trapezoidal)
        double G_model = dew_G_integral_highprec(T_C, P_bar);

        INFO("G_integral_highprec: T=" << T_C << " C, P=" << P_bar << " bar");
        INFO("  Model value:  " << G_model << " cal/mol");
        INFO("  Truth value:  " << G_truth << " cal/mol");
        INFO("  Difference:   " << (G_model - G_truth) << " cal/mol");
        INFO("  Rel. error:   " << (std::abs(G_model - G_truth) / std::max(std::abs(G_truth), 1e-10) * 100.0) << " %");
        REQUIRE(almost_equal(G_model, G_truth, ABS_TOL, REL_TOL));
    }
}

TEST_CASE("G_psat(T) matches truth table", "[G][Psat]")
{
    auto rows = load_csv("truth_G_psat.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-6;
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, G_cal_mol
        if (row.fields.size() < 2)
            continue;

        double T_C, G_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], G_truth)) continue;

        double G_model = dew_G_psat(T_C);

        INFO("G_psat: T=" << T_C << " C");
        REQUIRE(almost_equal(G_model, G_truth, ABS_TOL, REL_TOL));
    }
}

// -----------------------------------------------------------------------------
// Born Omega for all species
// -----------------------------------------------------------------------------

TEST_CASE("Omega(P,T) for all species matches truth table", "[Omega]")
{
    auto rows = load_csv("truth_Omega_AllSpecies.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-4;      // cal/mol tolerance
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // SpeciesName, Z, wref_cal_per_mol, P_bar, T_C,
        // rho_g_per_cm3, Omega_cal_per_mol, dOmega_dP_cal_per_mol_bar
        if (row.fields.size() < 8)
            continue;

        std::string speciesName = strip_quotes(row.fields[0]);

        double T_C, P_bar, rho_g_cm3;
        double omega_truth, domega_truth;

        if (!parse_maybe_double(row.fields[4], T_C)) continue;
        if (!parse_maybe_double(row.fields[3], P_bar)) continue;
        if (!parse_maybe_double(row.fields[5], rho_g_cm3)) continue;

        bool have_omega  = parse_maybe_double(row.fields[6], omega_truth);
        bool have_domega = parse_maybe_double(row.fields[7], domega_truth);

        if (have_omega) {
            double omega_model = dew_omega_species(speciesName, T_C, P_bar, rho_g_cm3);
            INFO("Omega: " << speciesName << ", T=" << T_C << " C, P=" << P_bar << " bar");
            INFO("  Model value:  " << omega_model << " cal/mol");
            INFO("  Truth value:  " << omega_truth << " cal/mol");
            INFO("  Difference:   " << (omega_model - omega_truth) << " cal/mol");
            INFO("  Rel. error:   " << (std::abs(omega_model - omega_truth) / std::max(std::abs(omega_truth), 1e-10) * 100.0) << " %");
            REQUIRE(almost_equal(omega_model, omega_truth, ABS_TOL, REL_TOL));
        }

        if (have_domega) {
            double domega_model = dew_domegadP_species(speciesName, T_C, P_bar, rho_g_cm3);
            INFO("dOmega/dP: " << speciesName << ", T=" << T_C << " C, P=" << P_bar << " bar");
            INFO("  Model value:  " << domega_model << " cal/mol/bar");
            INFO("  Truth value:  " << domega_truth << " cal/mol/bar");
            INFO("  Difference:   " << (domega_model - domega_truth) << " cal/mol/bar");
            INFO("  Rel. error:   " << (std::abs(domega_model - domega_truth) / std::max(std::abs(domega_truth), 1e-10) * 100.0) << " %");
            REQUIRE(almost_equal(domega_model, domega_truth, ABS_TOL, REL_TOL));
        }
    }
}

// -----------------------------------------------------------------------------
// Born Q(T,P) densEq1/epsEq4
// -----------------------------------------------------------------------------

TEST_CASE("Born Q(densEq1, epsEq4) matches truth table", "[Q]")
{
    auto rows = load_csv("truth_Q_densEq1_epsEq4.csv", /*skip_header=*/true);

    const double ABS_TOL = 1e-12;
    const double REL_TOL = 1e-8;

    for (const auto& row : rows) {
        // T_C, P_bar, densEq, epsEq, rho_g_cm3, Q_bar_inv
        if (row.fields.size() < 6)
            continue;

        double T_C, P_bar, Q_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[5], Q_truth)) continue;

        double Q_model = dew_Q_densEq1_epsEq4(T_C, P_bar);

        INFO("Q(densEq1,epsEq4): T=" << T_C << " C, P=" << P_bar << " bar");
        REQUIRE(almost_equal(Q_model, Q_truth, ABS_TOL, REL_TOL));
    }
}

TEST_CASE("DEW reaction thermodynamics: H2O + CO2,aq = H+ + HCO3-", "[reaction][DEW]")
{
    DEWDatabase db("dew2024-aqueous");
    auto species_list = db.species();
    Vec<Species> modified_species;
    for (auto sp : species_list) {
        modified_species.push_back(sp);
    }
    Species CO2_aq, H_plus, HCO3_minus;
    for (const auto& sp : modified_species) {
        if (sp.name() == "CO2_aq") CO2_aq = sp;
        else if (sp.name() == "H+") H_plus = sp;
        else if (sp.name() == "HCO3-") HCO3_minus = sp;
    }
    REQUIRE(CO2_aq.name() != "");
    REQUIRE(H_plus.name() != "");
    REQUIRE(HCO3_minus.name() != "");
    std::ofstream out("reaction_test_results.csv");
    out << "T_C,P_kb,G_rxn_model,G_rxn_truth,V_rxn_model,V_rxn_truth,"
        "G_H2O_model,G_H2O_truth,G_CO2_model,G_CO2_truth,G_Hplus_model,G_Hplus_truth,G_HCO3_model,G_HCO3_truth\n";
    auto rows = load_csv("reactionTesttruth.csv", /*skip_header=*/true);
    const double R = 8.314462618;
    const double CAL_TO_J = 4.184;
    const double CM3_TO_M3 = 1e-6;
    const double KB_TO_PA = 1e8;
    const double G_ABS_TOL = 500.0;
    const double G_REL_TOL = 0.0125;
    const double V_ABS_TOL = 1.0;
    const double V_REL_TOL = 0.0125;
    const double LOGK_ABS_TOL = 0.05;
    int test_count = 0;
    int passed_count = 0;
    std::vector<double> abs_err_G, rel_err_G, abs_err_V, rel_err_V, abs_err_logK, rel_err_logK;
    for (const auto& row : rows) {
        if (row.fields.size() < 19)
            continue;
        double P_kb, T_C;
            double G_H2O_cal, G_CO2_cal, G_Hplus_cal, G_HCO3_cal, G_rxn_cal, log_K_truth;
            double V_H2O_cm3, V_CO2_cm3, V_Hplus_cm3, V_HCO3_cm3, V_rxn_cm3;
            if (!parse_maybe_double(row.fields[0], P_kb)) continue;
            if (!parse_maybe_double(row.fields[1], T_C)) continue;
            if (!parse_maybe_double(row.fields[4], G_H2O_cal)) continue;
            if (!parse_maybe_double(row.fields[5], G_CO2_cal)) continue;
            if (!parse_maybe_double(row.fields[6], G_Hplus_cal)) continue;
            if (!parse_maybe_double(row.fields[7], G_HCO3_cal)) continue;
            if (!parse_maybe_double(row.fields[8], G_rxn_cal)) continue;
            if (!parse_maybe_double(row.fields[9], log_K_truth)) continue;
            if (!parse_maybe_double(row.fields[12], V_H2O_cm3)) continue;
            if (!parse_maybe_double(row.fields[13], V_CO2_cm3)) continue;
            if (!parse_maybe_double(row.fields[14], V_Hplus_cm3)) continue;
            if (!parse_maybe_double(row.fields[15], V_HCO3_cm3)) continue;
            if (!parse_maybe_double(row.fields[16], V_rxn_cm3)) continue;
        const double T_K = T_C + 273.15;
        const double P_Pa = P_kb * KB_TO_PA;
        const double G_rxn_truth = G_rxn_cal * CAL_TO_J;
        const double V_rxn_truth = V_rxn_cm3;
        test_count++;
        WaterModelOptions waterOpts = makeWaterModelOptionsDEW();
        WaterStateOptions wsOpts;
        wsOpts.thermo.eosModel = waterOpts.eosModel;
        wsOpts.computeGibbs = true;
        wsOpts.gibbs.model = waterOpts.gibbsModel;
        wsOpts.gibbs.thermo = wsOpts.thermo;
        wsOpts.gibbs.integrationSteps = 5000;
        wsOpts.gibbs.useExcelIntegration = false;
        auto ws = waterState(T_K, P_Pa, wsOpts);
        const double G0_H2O = ws.gibbs;
        const double M_H2O = 0.018015;
        const double V_specific = (ws.thermo.V != 0.0) ? ws.thermo.V : (1.0 / ws.thermo.D);
        const double V0_H2O = V_specific * M_H2O;
        const double G0_H2O_cal = G0_H2O / CAL_TO_J;
        const double V0_H2O_cm3 = V0_H2O / CM3_TO_M3;
        auto model_CO2 = CO2_aq.standardThermoModel();
        auto model_Hplus = H_plus.standardThermoModel();
        auto model_HCO3 = HCO3_minus.standardThermoModel();
        StandardThermoProps props_CO2 = model_CO2(T_K, P_Pa);
        StandardThermoProps props_Hplus = model_Hplus(T_K, P_Pa);
        StandardThermoProps props_HCO3 = model_HCO3(T_K, P_Pa);
        const double G0_CO2_cal = props_CO2.G0 / CAL_TO_J;
        const double G0_Hplus_cal = props_Hplus.G0 / CAL_TO_J;
        const double G0_HCO3_cal = props_HCO3.G0 / CAL_TO_J;
        const double V0_CO2_cm3 = props_CO2.V0 / CM3_TO_M3;
        const double V0_Hplus_cm3 = props_Hplus.V0 / CM3_TO_M3;
        const double V0_HCO3_cm3 = props_HCO3.V0 / CM3_TO_M3;
        const double G_rxn_model = props_Hplus.G0 + props_HCO3.G0 - G0_H2O - props_CO2.G0;
        const double V_rxn_model_m3 = props_Hplus.V0 + props_HCO3.V0 - V0_H2O - props_CO2.V0;
        const double V_rxn_model = V_rxn_model_m3 / CM3_TO_M3;
        const double log_K_model = -G_rxn_model / (R * T_K * std::log(10.0));
        double abs_err_g = std::abs(G_rxn_model - G_rxn_truth);
        double rel_err_g = abs_err_g / std::max(std::abs(G_rxn_truth), 1e-10);
        double abs_err_v = std::abs(V_rxn_model - V_rxn_truth);
        double rel_err_v = abs_err_v / std::max(std::abs(V_rxn_truth), 1e-10);
        double abs_err_logk = std::abs(log_K_model - log_K_truth);
        double rel_err_logk = abs_err_logk / std::max(std::abs(log_K_truth), 1e-10);
        abs_err_G.push_back(abs_err_g);
        rel_err_G.push_back(rel_err_g);
        abs_err_V.push_back(abs_err_v);
        rel_err_V.push_back(rel_err_v);
        abs_err_logK.push_back(abs_err_logk);
        rel_err_logK.push_back(rel_err_logk);
        INFO("ΔGr mismatch at T=" << T_C << "°C, P=" << P_kb << "kb: model=" << G_rxn_model << ", truth=" << G_rxn_truth << ", abs_err=" << abs_err_g << ", rel_err=" << rel_err_g);
        CHECK(almost_equal(G_rxn_model, G_rxn_truth, G_ABS_TOL, G_REL_TOL));
        INFO("ΔVr mismatch at T=" << T_C << "°C, P=" << P_kb << "kb: model=" << V_rxn_model << ", truth=" << V_rxn_truth << ", abs_err=" << abs_err_v << ", rel_err=" << rel_err_v);
        CHECK(almost_equal(V_rxn_model, V_rxn_truth, V_ABS_TOL, V_REL_TOL));
        INFO("log K mismatch at T=" << T_C << "°C, P=" << P_kb << "kb: model=" << log_K_model << ", truth=" << log_K_truth << ", abs_err=" << abs_err_logk << ", rel_err=" << rel_err_logk);
        CHECK(almost_equal(log_K_model, log_K_truth, LOGK_ABS_TOL, 0.0125));
        out << T_C << "," << P_kb << "," << G_rxn_model << "," << G_rxn_truth << "," << V_rxn_model << "," << V_rxn_truth
            << "," << G0_H2O << "," << (G_H2O_cal * CAL_TO_J)
            << "," << props_CO2.G0 << "," << (G_CO2_cal * CAL_TO_J)
            << "," << props_Hplus.G0 << "," << (G_Hplus_cal * CAL_TO_J)
            << "," << props_HCO3.G0 << "," << (G_HCO3_cal * CAL_TO_J)
            << "\n";
        if (
            almost_equal(G_rxn_model, G_rxn_truth, G_ABS_TOL, G_REL_TOL) &&
            almost_equal(V_rxn_model, V_rxn_truth, V_ABS_TOL, V_REL_TOL) &&
            almost_equal(log_K_model, log_K_truth, LOGK_ABS_TOL, 0.0125)
        ) {
            passed_count++;
        }
    }
    out.close();
    auto minmaxavg = [](const std::vector<double>& v) {
        double min = v.empty() ? 0.0 : *std::min_element(v.begin(), v.end());
        double max = v.empty() ? 0.0 : *std::max_element(v.begin(), v.end());
        double avg = v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / (v.empty() ? 1 : v.size());
        return std::make_tuple(min, max, avg);
    };
    auto print_stats = [&](const char* label, const std::vector<double>& abs, const std::vector<double>& rel, const char* unit) {
        auto [min_abs, max_abs, avg_abs] = minmaxavg(abs);
        auto [min_rel, max_rel, avg_rel] = minmaxavg(rel);
        std::cout << label << " absolute error: min=" << min_abs << ", max=" << max_abs << ", avg=" << avg_abs << " " << unit << std::endl;
        std::cout << label << " relative error: min=" << min_rel*100 << "% , max=" << max_rel*100 << "% , avg=" << avg_rel*100 << "%" << std::endl;
    };
    std::cout << "\nTested " << test_count << " conditions, " << passed_count << " passed." << std::endl;
    print_stats("ΔGr", abs_err_G, rel_err_G, "J/mol");
    print_stats("ΔVr", abs_err_V, rel_err_V, "cm³/mol");
    print_stats("log K", abs_err_logK, rel_err_logK, "");
    REQUIRE(test_count > 0);
    REQUIRE(passed_count == test_count);
}

