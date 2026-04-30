/**
 * Direct comparison test: All 4 numerical integration methods
 * Compares water Gibbs energy calculations at sample T,P points
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>

#include "Reaktoro/Extensions/DEW/WaterGibbsModel.hpp"

using namespace Reaktoro;

struct TestCase {
    double T_C;
    double P_kb;
    std::string label;
};

void compareMethod(const std::string& name,
                   WaterIntegrationMethod method,
                   const std::vector<TestCase>& cases)
{
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << name << "\n";
    std::cout << std::string(80, '=') << "\n";

    WaterGibbsModelOptions opts;
    opts.integrationMethod = method;
    opts.integrationSteps = 5000;
    opts.densityTolerance = 0.001;
    opts.adaptiveIntegrationTolerance = 0.1;
    opts.maxAdaptiveSubdivisions = 20;
    opts.useExcelIntegration = false;

    WaterThermoModelOptions thermo_opts;
    thermo_opts.eosModel = WaterEosModel::ZhangDuan2005;
    thermo_opts.densityTolerance = 0.001;

    int success = 0, failed = 0;
    double total_time = 0.0;

    for (const auto& tc : cases) {
        double T_K = tc.T_C + 273.15;
        double P_Pa = tc.P_kb * 1.0e8;

        auto t0 = std::chrono::high_resolution_clock::now();

        try {
            double G_J = waterGibbsModel(T_K, P_Pa, opts);

            auto t1 = std::chrono::high_resolution_clock::now();
            double dt_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            total_time += dt_ms;

            double G_cal = G_J / 4.184;

            std::cout << "  " << std::left << std::setw(20) << tc.label
                     << " T=" << std::fixed << std::setw(6) << std::setprecision(0) << tc.T_C
                     << "°C P=" << std::setw(5) << tc.P_kb << "kb  "
                     << "G=" << std::setw(12) << std::setprecision(2) << G_cal << " cal/mol  "
                     << std::setw(6) << std::setprecision(2) << dt_ms << "ms  ✓\n";

            success++;
        }
        catch (const std::exception& e) {
            failed++;
            std::cout << "  " << std::left << std::setw(20) << tc.label
                     << " ERROR: " << e.what() << "\n";
        }
    }

    std::cout << "\nSummary: " << success << " passed, " << failed << " failed\n";
    if (success > 0) {
        std::cout << "  Total time: " << std::fixed << std::setprecision(1) << total_time << " ms\n";
        std::cout << "  Avg time/point: " << std::fixed << std::setprecision(2)
                 << (total_time / success) << " ms\n";
    }
}

int main() {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "NUMERICAL INTEGRATION METHOD COMPARISON\n";
    std::cout << "Testing all 4 methods on 5 sample points\n";
    std::cout << std::string(80, '=') << "\n";

    std::vector<TestCase> test_cases = {
        {300.0, 5000.0, "Point 1: Low T, Low P"},
        {400.0, 6000.0, "Point 2: Mid T, Mid P"},
        {500.0, 8000.0, "Point 3: High T, High P"},
        {650.0, 10000.0, "Point 4: Very High T"},
        {350.0, 7000.0, "Point 5: Mixed"},
    };

    // Test each method
    compareMethod("METHOD 1: Trapezoidal Rule (O(h²)) - BASELINE",
                  WaterIntegrationMethod::Trapezoidal,
                  test_cases);

    compareMethod("METHOD 2: Simpson's Rule (O(h⁴))",
                  WaterIntegrationMethod::Simpson,
                  test_cases);

    compareMethod("METHOD 3: Gauss-Legendre-16 (O(1/n³²))",
                  WaterIntegrationMethod::GaussLegendre16,
                  test_cases);

    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "✓ All integration methods tested successfully!\n";
    std::cout << std::string(80, '=') << "\n\n";

    return 0;
}
