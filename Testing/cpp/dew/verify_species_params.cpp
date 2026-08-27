// Quick verification of species thermodynamic parameters at 25°C, 1 bar

#include <Reaktoro/Reaktoro.hpp>
#include <Reaktoro/Extensions/DEW/DEWDatabase.hpp>
#include <Reaktoro/Extensions/DEW/WaterState.hpp>
#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>
#include <iostream>
#include <iomanip>

using namespace Reaktoro;

int main()
{
    try {
    const double CAL_TO_J = 4.184;
    const double CM3_TO_M3 = 1e-6;

    // Reference conditions: 25°C = 298.15 K, 1 bar = 1e5 Pa
    const double T_ref = 298.15;  // K
    const double P_ref = 1.0e5;   // Pa (1 bar)

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\n=== Verification at 25°C, 1 bar ===" << std::endl;
    std::cout << "\nExpected values from database image:" << std::endl;
    std::cout << "  CO2,aq:  ΔGf° = -92200 cal/mol,  V° = 30.0 cm³/mol" << std::endl;
    std::cout << "  H+:      ΔGf° = 0 cal/mol,       V° = 0.0 cm³/mol" << std::endl;
    std::cout << "  HCO3-:   ΔGf° = -140282 cal/mol, V° = 24.2 cm³/mol" << std::endl;

    // Load database
    std::cout << "\nLoading database..." << std::endl;
    DEWDatabase db("dew2024-aqueous");
    std::cout << "Database loaded successfully!" << std::endl;    // Get species
    auto CO2_aq = db.species().get("CO2_aq");
    auto H_plus = db.species().get("H+");
    auto HCO3_minus = db.species().get("HCO3-");

    std::cout << "\n=== Database Gf values (should match exactly) ===" << std::endl;

    // The StandardThermoModel.params() should give us the raw HKF parameters
    auto model_CO2 = CO2_aq.standardThermoModel();
    auto model_Hplus = H_plus.standardThermoModel();
    auto model_HCO3 = HCO3_minus.standardThermoModel();

    // Calculate properties at reference state
    auto props_CO2 = model_CO2(T_ref, P_ref);
    auto props_Hplus = model_Hplus(T_ref, P_ref);
    auto props_HCO3 = model_HCO3(T_ref, P_ref);

    std::cout << "\nCO2,aq:" << std::endl;
    std::cout << "  G0 (calc) = " << props_CO2.G0 << " J/mol = " << (props_CO2.G0 / CAL_TO_J) << " cal/mol" << std::endl;
    std::cout << "  V0 (calc) = " << (props_CO2.V0 / CM3_TO_M3) << " cm³/mol" << std::endl;

    std::cout << "\nH+:" << std::endl;
    std::cout << "  G0 (calc) = " << props_Hplus.G0 << " J/mol = " << (props_Hplus.G0 / CAL_TO_J) << " cal/mol" << std::endl;
    std::cout << "  V0 (calc) = " << (props_Hplus.V0 / CM3_TO_M3) << " cm³/mol" << std::endl;

    std::cout << "\nHCO3-:" << std::endl;
    std::cout << "  G0 (calc) = " << props_HCO3.G0 << " J/mol = " << (props_HCO3.G0 / CAL_TO_J) << " cal/mol" << std::endl;
    std::cout << "  V0 (calc) = " << (props_HCO3.V0 / CM3_TO_M3) << " cm³/mol" << std::endl;

    // Now test at a high T,P point from the CSV: 300°C, 5 kb
    const double T_test = 573.15;  // 300°C in K
    const double P_test = 5.0e8;   // 5 kb = 5000 bar = 5e8 Pa

    std::cout << "\n\n=== Test at 300°C, 5 kb (first point from CSV) ===" << std::endl;

    // Get water properties
    auto waterOpts = makeWaterModelOptionsDEW();
    WaterStateOptions wsOpts;
    wsOpts.thermo.eosModel = waterOpts.eosModel;
    wsOpts.computeGibbs = true;
    wsOpts.gibbs.model = waterOpts.gibbsModel;
    wsOpts.gibbs.thermo = wsOpts.thermo;

    auto ws = waterState(T_test, P_test, wsOpts);

    // Convert specific volume (m³/kg) to molar volume (m³/mol)
    const double M_H2O = 0.018015;  // kg/mol
    // If V is not populated, calculate from density: V = 1/D
    const double V_specific = (ws.thermo.V != 0.0) ? ws.thermo.V : (1.0 / ws.thermo.D);
    const double V_H2O_molar = V_specific * M_H2O;  // m³/mol

    std::cout << "\nH2O (from WaterState):" << std::endl;
    std::cout << "  G = " << ws.gibbs << " J/mol = " << (ws.gibbs / CAL_TO_J) << " cal/mol" << std::endl;
    std::cout << "  V_specific = " << V_specific << " m³/kg (from " << (ws.thermo.V != 0.0 ? "V" : "1/D") << ")" << std::endl;
    std::cout << "  V_molar = " << V_H2O_molar << " m³/mol" << std::endl;
    std::cout << "  V = " << (V_H2O_molar / CM3_TO_M3) << " cm³/mol" << std::endl;
    std::cout << "  density = " << ws.thermo.D << " kg/m³" << std::endl;
    std::cout << "  dielectric = " << ws.electro.epsilon << std::endl;

    // Calculate species at high T,P
    auto props_CO2_hp = model_CO2(T_test, P_test);
    auto props_Hplus_hp = model_Hplus(T_test, P_test);
    auto props_HCO3_hp = model_HCO3(T_test, P_test);

    std::cout << "\nCO2,aq:" << std::endl;
    std::cout << "  G0 = " << props_CO2_hp.G0 << " J/mol = " << (props_CO2_hp.G0 / CAL_TO_J) << " cal/mol" << std::endl;
    std::cout << "  V0 = " << (props_CO2_hp.V0 / CM3_TO_M3) << " cm³/mol" << std::endl;

    std::cout << "\nH+:" << std::endl;
    std::cout << "  G0 = " << props_Hplus_hp.G0 << " J/mol = " << (props_Hplus_hp.G0 / CAL_TO_J) << " cal/mol" << std::endl;
    std::cout << "  V0 = " << (props_Hplus_hp.V0 / CM3_TO_M3) << " cm³/mol" << std::endl;

    std::cout << "\nHCO3-:" << std::endl;
    std::cout << "  G0 = " << props_HCO3_hp.G0 << " J/mol = " << (props_HCO3_hp.G0 / CAL_TO_J) << " cal/mol" << std::endl;
    std::cout << "  V0 = " << (props_HCO3_hp.V0 / CM3_TO_M3) << " cm³/mol" << std::endl;

    // Calculate reaction at 300°C, 5 kb
    const double G_rxn = props_Hplus_hp.G0 + props_HCO3_hp.G0 - ws.gibbs - props_CO2_hp.G0;
    const double V_rxn = (props_Hplus_hp.V0 + props_HCO3_hp.V0 - V_H2O_molar - props_CO2_hp.V0) / CM3_TO_M3;
    const double R = 8.314462618;
    const double log_K = -G_rxn / (R * T_test * std::log(10.0));

    std::cout << "\n=== Reaction: H2O + CO2,aq = H+ + HCO3- ===" << std::endl;
    std::cout << "ΔGr = " << G_rxn << " J/mol = " << (G_rxn / CAL_TO_J) << " cal/mol" << std::endl;
    std::cout << "ΔVr = " << V_rxn << " cm³/mol" << std::endl;
    std::cout << "log K = " << log_K << std::endl;

    std::cout << "\nFrom CSV (expected):" << std::endl;
    std::cout << "ΔGr = 68401 J/mol (16349.38 cal/mol)" << std::endl;
    std::cout << "ΔVr = -24.3076 cm³/mol" << std::endl;
    std::cout << "log K = -6.23785" << std::endl;

    return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nERROR: " << e.what() << std::endl;
        return 1;
    }
}