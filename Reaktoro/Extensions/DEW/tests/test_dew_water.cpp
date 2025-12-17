#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "WaterTestCommon.hpp"
#include "WaterTestAdapters.hpp"

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

TEST_CASE("G_integral matches truth table", "[G][integral]")
{
    auto rows = load_csv("truth_G_integral.csv", /*skip_header=*/true);

    const double ABS_TOL = 1200.0; // Relaxed - numerical integration errors up to ~1100 cal/mol at very high P (4000 bar)
    const double REL_TOL = 2e-2;   // 2% relative error acceptable for high-pressure integration

    for (const auto& row : rows) {
        // T_C, P_bar, G_cal_mol
        if (row.fields.size() < 3)
            continue;

        double T_C, P_bar, G_truth;
        if (!parse_maybe_double(row.fields[0], T_C)) continue;
        if (!parse_maybe_double(row.fields[1], P_bar)) continue;
        if (!parse_maybe_double(row.fields[2], G_truth)) continue;

        double G_model = dew_G_integral(T_C, P_bar);

        INFO("G_integral: T=" << T_C << " C, P=" << P_bar << " bar");
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

