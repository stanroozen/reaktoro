// Test program for Perple_X hybrid electrolyte machinery
// Validates dielectric constant, Debye-Hückel, and g-function calculations

#include "PerpleXElectrolyte.hpp"
#include "PerpleXFluidModel.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Reaktoro::PerpleX;

void printHeader(const std::string& title)
{
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

void testEpsh2o()
{
    printHeader("Test 1: Pure Water Dielectric Constant");

    // Test conditions
    std::vector<std::pair<double, double>> testCases = {
        {298.15, 18.0},   // 25°C, standard volume
        {373.15, 18.8},   // 100°C, expanded volume
        {473.15, 22.0},   // 200°C, high T
        {573.15, 28.0}    // 300°C, supercritical
    };

    std::cout << std::setw(12) << "T (K)"
              << std::setw(15) << "V (cm3/mol)"
              << std::setw(15) << "epsilon"
              << "\n" << std::string(42, '-') << "\n";

    for (const auto& [T, V] : testCases)
    {
        double v_jbar = V / 10.0;  // Convert to J/bar
        double eps = epsh2o(v_jbar, T);

        std::cout << std::setw(12) << std::fixed << std::setprecision(2) << T
                  << std::setw(15) << std::setprecision(2) << V
                  << std::setw(15) << std::setprecision(4) << eps
                  << "\n";
    }

    // Expected: ~78.5 at 25°C, decreasing with temperature
}

void testGfunc()
{
    printHeader("Test 2: Shock g-function (HKF)");

    // Test conditions (rho, P, T)
    std::vector<std::tuple<double, double, double>> testCases = {
        {0.997, 1.0, 298.15},      // Ambient water
        {0.95, 500.0, 473.15},     // Moderate P-T
        {0.80, 1000.0, 573.15},    // High P-T
        {0.60, 1500.0, 673.15},    // Very high P-T
        {0.35, 500.0, 623.15}      // Low density limit
    };

    std::cout << std::setw(12) << "rho (g/cm3)"
              << std::setw(12) << "P (bar)"
              << std::setw(12) << "T (K)"
              << std::setw(15) << "g (Angstrom)"
              << "\n" << std::string(51, '-') << "\n";

    for (const auto& [rho, P, T] : testCases)
    {
        double g = gfunc(rho, P, T);

        std::cout << std::setw(12) << std::fixed << std::setprecision(3) << rho
                  << std::setw(12) << std::setprecision(1) << P
                  << std::setw(12) << std::setprecision(2) << T
                  << std::setw(15) << std::setprecision(6) << g
                  << "\n";
    }

    // Expected: g decreases with density, becomes negative at low density
}

void testDebyeHuckel()
{
    printHeader("Test 3: Debye-Hückel Factor");

    // Test conditions (msol, vsolv, epsilon, T)
    std::vector<std::tuple<double, double, double, double>> testCases = {
        {1.0, 18.0, 78.47, 298.15},   // Standard aqueous
        {1.0, 18.8, 55.0, 373.15},    // High T water
        {1.0, 22.0, 35.0, 473.15},    // Very high T
        {2.0, 18.0, 78.47, 298.15}    // Concentrated solution
    };

    std::cout << std::setw(12) << "msol"
              << std::setw(12) << "V (cm3)"
              << std::setw(12) << "epsilon"
              << std::setw(12) << "T (K)"
              << std::setw(15) << "adh"
              << "\n" << std::string(63, '-') << "\n";

    for (const auto& [msol, vsolv, eps, T] : testCases)
    {
        double adh = debyeHuckel(msol, vsolv, eps, T);

        std::cout << std::setw(12) << std::fixed << std::setprecision(2) << msol
                  << std::setw(12) << std::setprecision(2) << vsolv
                  << std::setw(12) << std::setprecision(2) << eps
                  << std::setw(12) << std::setprecision(2) << T
                  << std::setw(15) << std::setprecision(6) << adh
                  << "\n";
    }

    // Expected: negative values, magnitude increases with T and msol
}

void testMixtureEps()
{
    printHeader("Test 4: Mixture Dielectric Constant");

    // Test H2O-CO2 mixture
    std::array<double, 19> vhyb{};
    std::array<double, 19> vf{};

    vhyb[0] = 18.0;   // H2O volume (cm3/mol)
    vhyb[1] = 32.5;   // CO2 volume (cm3/mol)

    std::vector<int> species = {1, 2};  // H2O, CO2

    std::cout << std::setw(12) << "X_H2O"
              << std::setw(12) << "X_CO2"
              << std::setw(15) << "epsilon"
              << "\n" << std::string(39, '-') << "\n";

    for (double xh2o = 0.0; xh2o <= 1.0; xh2o += 0.2)
    {
        double xco2 = 1.0 - xh2o;

        // Compute volume fractions
        double hyvol = xh2o * vhyb[0] + xco2 * vhyb[1];
        vf[0] = xh2o * vhyb[0] / hyvol;
        vf[1] = xco2 * vhyb[1] / hyvol;

        double eps = geteps(vhyb, vf, species, 2, 298.15);

        std::cout << std::setw(12) << std::fixed << std::setprecision(2) << xh2o
                  << std::setw(12) << std::setprecision(2) << xco2
                  << std::setw(15) << std::setprecision(4) << eps
                  << "\n";
    }

    // Expected: epsilon decreases from ~78 (pure H2O) to ~1.5 (pure CO2)
}

void testHybridSolventGibbs()
{
    printHeader("Test 5: Hybrid Solvent Gibbs Energy");

    std::array<double, 19> yf{};
    std::array<double, 19> ghybrid{};

    // Set up simple binary mixture
    ghybrid[0] = 1.0;  // Normalized fugacity ratio H2O
    ghybrid[1] = 0.95; // Normalized fugacity ratio CO2

    std::cout << std::setw(12) << "y_H2O"
              << std::setw(12) << "y_CO2"
              << std::setw(18) << "gsolv (J/mol)"
              << "\n" << std::string(42, '-') << "\n";

    double T = 298.15;

    for (double y1 = 0.2; y1 <= 1.0; y1 += 0.2)
    {
        yf[0] = y1;
        yf[1] = 1.0 - y1;

        std::vector<int> solventSpecies = {1, 2};
        std::array<double, 19> g0{};
        g0[0] = 1.0;
        g0[1] = 1.0;
        double gsolv = computeHybridSolventGibbs(yf, ghybrid, g0, solventSpecies, 2, T);

        std::cout << std::setw(12) << std::fixed << std::setprecision(2) << yf[0]
                  << std::setw(12) << std::setprecision(2) << yf[1]
                  << std::setw(18) << std::setprecision(2) << gsolv
                  << "\n";
    }
}

void testFullIntegration()
{
    printHeader("Test 6: Full Fluid Model Integration");

    // Test complete fluid model with electrolyte machinery
    PerpleXFluidOptions options;
    options.enableElectrolyte = true;
    options.hybridSpecies = {1, 2};  // H2O, CO2
    options.hybridOptions.water = HybridEosOptions::WaterEos::Mrk;
    options.hybridOptions.co2 = HybridEosOptions::CO2Eos::Mrk;

    std::array<double, 19> y{};
    y[0] = 0.8;  // H2O
    y[1] = 0.2;  // CO2

    std::vector<int> species = {1, 2};
    double P = 1000.0;  // bar
    double T = 473.15;  // K

    PerpleXFluidModel model;
    auto state = model.compute(species, y, P, T, options);

    std::cout << "\nFluid State at P=" << P << " bar, T=" << T << " K:\n";
    std::cout << "  Molar volume: " << state.vol << " cm3/mol\n";
    std::cout << "  Hybrid volume: " << state.hyvol << " cm3/mol\n";
    std::cout << "  Dielectric constant: " << state.dielectric.epsilon << "\n";
    std::cout << "  Debye-Hückel factor: " << state.dielectric.adh << "\n";
    std::cout << "  Shock g-function: " << state.dielectric.gf << " Angstrom\n";
    std::cout << "  Solvent Gibbs: " << state.gsolv << " J/mol\n";

    std::cout << "\n  Fugacity coefficients:\n";
    for (int i = 0; i < 2; ++i)
    {
        std::cout << "    Species " << species[i] << ": "
                  << std::exp(state.ln_f[species[i]-1]) / (P * y[species[i]-1])
                  << "\n";
    }
}

int main()
{
    try
    {
        std::cout << "\n" << std::string(60, '=');
        std::cout << "\nPerple_X Hybrid Electrolyte Machinery Tests";
        std::cout << "\n" << std::string(60, '=') << "\n";

        testEpsh2o();
        testGfunc();
        testDebyeHuckel();
        testMixtureEps();
        testHybridSolventGibbs();
        testFullIntegration();

        printHeader("All Tests Completed Successfully");

        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nERROR: " << e.what() << "\n";
        return 1;
    }
}
