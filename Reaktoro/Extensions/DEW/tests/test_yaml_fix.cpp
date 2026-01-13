// Quick test loading YAML directly to verify the fix
#include <Reaktoro/Reaktoro.hpp>
#include <Reaktoro/Extensions/DEW/WaterState.hpp>
#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>
#include <iostream>
#include <iomanip>

using namespace Reaktoro;

int main()
{
    const double CAL_TO_J = 4.184;
    const double CM3_TO_M3 = 1e-6;
    const double T_ref = 298.15;  // K
    const double P_ref = 1.0e5;   // Pa (1 bar)

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\n=== Testing YAML fix for a1 parameter ===" << std::endl;

    // Load database from YAML file directly
    const auto yaml_path = "C:/Users/stanroozen/Documents/Projects/reaktoro-dev/reaktoro/embedded/databases/DEW/dew2024-aqueous.yaml";

    Database db = Database::fromFile(yaml_path);

    auto CO2_aq = db.species().get("CO2(0)");

    // Get the HKF parameters
    auto model = CO2_aq.standardThermoModel();
    auto params = model.params();

    std::cout << "\nCO2,aq HKF parameters:" << std::endl;
    std::cout << "  a1 = " << std::scientific << params["HKF.a1"].asFloat() << " J/(mol·Pa) = m³/mol" << std::endl;
    std::cout << "  a2 = " << params["HKF.a2"].asFloat() << " J/mol" << std::endl;
    std::cout << "  a3 = " << params["HKF.a3"].asFloat() << " J·K/(mol·Pa)" << std::endl;
    std::cout << "  a4 = " << params["HKF.a4"].asFloat() << " J·K/mol" << std::endl;

    // Calculate V0 at 25°C, 1 bar
    auto props = model(T_ref, P_ref);
    double V0 = props.V0;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\nCalculated V° at 25°C, 1 bar:" << std::endl;
    std::cout << "  V° = " << V0 / CM3_TO_M3 << " cm³/mol" << std::endl;
    std::cout << "\nExpected from Excel: V° = 30.0 cm³/mol" << std::endl;

    double error_pct = std::abs((V0 / CM3_TO_M3) - 30.0) / 30.0 * 100.0;
    std::cout << "Error: " << error_pct << "%" << std::endl;

    if (error_pct < 5.0) {
        std::cout << "\n✓ SUCCESS: Volume calculation is correct!" << std::endl;
    } else {
        std::cout << "\n✗ FAILED: Volume still wrong!" << std::endl;
    }

    return 0;
}
