// Test DEW reaction thermodynamics: H2O + CO2,aq = H+ + HCO3-
// Regression test against reactionTesttruth.csv
//
// Units in CSV:
//   Pressure: kb (kilobar) -> multiply by 1e8 Pa = 1e5 bar = 1000 bar/kb
//   Temperature: °C -> convert to K by adding 273.15
//   ΔG°: cal/mol -> multiply by 4.184 to get J/mol
//   ΔV°: cm³/mol -> multiply by 1e-6 to get m³/mol
//
// The CSV contains:
//   - DeltaGo_H2O, DeltaGo_CO2_aq, DeltaGo_H+, DeltaGo_HCO3- (cal/mol)
//   - DeltaGro = DeltaGo_H+ + DeltaGo_HCO3- - DeltaGo_H2O - DeltaGo_CO2_aq (cal/mol)
//   - log_K = -DeltaGro / (R*T*ln(10))
//   - DeltaVo_H2O, DeltaVo_CO2_aq, DeltaVo_H+, DeltaVo_HCO3- (cm³/mol)
//   - DeltaVr = DeltaVo_H+ + DeltaVo_HCO3- - DeltaVo_H2O - DeltaVo_CO2_aq (cm³/mol)

#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "WaterTestCommon.hpp"

#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Extensions/DEW/DEWDatabase.hpp>
#include <Reaktoro/Extensions/DEW/WaterState.hpp>
#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelDEW.hpp>
#include <fstream>
#include <cstdio>

using namespace Reaktoro;

// Constants
const double R = 8.314462618;  // J/(mol·K) - gas constant
const double CAL_TO_J = 4.184;  // cal to J conversion
const double CM3_TO_M3 = 1e-6;  // cm³ to m³ conversion
const double KB_TO_PA = 1e8;    // kilobar to Pascal (1 kb = 1000 bar = 1e8 Pa)

TEST_CASE("DEW reaction thermodynamics: H2O + CO2,aq = H+ + HCO3-", "[reaction][DEW]")
{
    // Load DEW database with Excel-compatible settings for comparison
    // Note: Excel VBA uses 100 bar tolerance during Gibbs integration
    DEWDatabase db("dew2024-aqueous");
    auto species_list = db.species();

    // Modify species to use Excel VBA tolerance (100 bar) for exact comparison
    // This affects all water property calculations (ρ, ε, Q, ω) used in HKF model
    Vec<Species> modified_species;
    for (auto sp : species_list) {
        // Extract DEW parameters from existing model
        auto model = sp.standardThermoModel();
        // For now, we'll use the species as-is and just set water calc tolerance
        // The tolerance will be applied via WaterStateOptions below
        modified_species.push_back(sp);
    }

    // Find required species (note: water is not in the database, we'll use WaterState for it)
    Species CO2_aq, H_plus, HCO3_minus;

    for (const auto& sp : modified_species) {
        if (sp.name() == "CO2_aq") CO2_aq = sp;
        else if (sp.name() == "H+") H_plus = sp;
        else if (sp.name() == "HCO3-") HCO3_minus = sp;
    }

    REQUIRE(CO2_aq.name() != "");
    REQUIRE(H_plus.name() != "");
    REQUIRE(HCO3_minus.name() != "");

    INFO("Found species:");
    INFO("  CO2,aq: " << CO2_aq.name());
    INFO("  H+: " << H_plus.name());
    INFO("  HCO3-: " << HCO3_minus.name());    // Load test data
    auto rows = load_csv("reactionTesttruth.csv", /*skip_header=*/true);

    // Prepare CSV for per-point error logging (overwrite any existing file)
    std::remove("reaction_errors.csv");
    {
        std::ofstream _hdr("reaction_errors.csv");
        _hdr << "T_C,P_kb,abs_err_g,rel_err_g,abs_err_v,rel_err_v,abs_err_logk\n";
    }

    // Tolerances for thermodynamic properties (tightened per request)
    const double G_ABS_TOL = 50.0;     // J/mol absolute tolerance for Gibbs energy
    const double G_REL_TOL = 0.001;    // 0.1% relative tolerance
    const double V_ABS_TOL = 1e-3;     // cm³/mol absolute tolerance for volume
    const double V_REL_TOL = 0.001;    // 0.1% relative tolerance
    const double LOGK_ABS_TOL = 0.01;  // Absolute tolerance for log K
    const double LOGK_REL_TOL = 0.001; // Relative tolerance for log K (0.1%)

    int test_count = 0;
    int passed_count = 0;
    // Error statistics
    std::vector<double> abs_err_G, rel_err_G, abs_err_V, rel_err_V, abs_err_logK, rel_err_logK;

    for (const auto& row : rows) {
        // CSV columns:
        // 0: Pressure_kb
        // 1: Temp_C
        // 2: rhoH2O_gcm-3
        // 3: epsilon
        // 4: DeltaGo_H2O_calmol-1
        // 5: DeltaGo_CO2_aq_calmol-1
        // 6: DeltaGo_H+_calmol-1
        // 7: DeltaGo_HCO3-_calmol-1
        // 8: DeltaGro_calmol-1
        // 9: log_K
        // 10: Pressure_kb_2 (duplicate)
        // 11: Temp_C_2 (duplicate)
        // 12: DeltaVo_H2O_cm3mol-1
        // 13: DeltaVo_CO2_aq_cm3mol-1
        // 14: DeltaVo_H+_cm3mol-1
        // 15: DeltaVo_HCO3-_cm3mol-1
        // 16: DeltaVr_cm3mol-1
        // 17: Temp_C_3 (duplicate)
        // 18: Pressure_kb_3 (duplicate)

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

        // Convert units
        const double T_K = T_C + 273.15;
        const double P_Pa = P_kb * KB_TO_PA;

        // Convert truth values to SI
        const double G_rxn_truth = G_rxn_cal * CAL_TO_J;  // J/mol
        const double V_rxn_truth = V_rxn_cm3;  // Keep in cm³/mol for comparison

        test_count++;

        // Configure water model options for DEW
        WaterModelOptions waterOpts = makeWaterModelOptionsDEW();

        // Configure water state options
        WaterStateOptions wsOpts;
        wsOpts.thermo.eosModel = waterOpts.eosModel;
        wsOpts.computeGibbs = true;
        wsOpts.gibbs.model = waterOpts.gibbsModel;
        wsOpts.gibbs.thermo = wsOpts.thermo;
        // Use high-precision integration (5000 steps, not Excel mode)
        wsOpts.gibbs.integrationSteps = 5000;
        wsOpts.gibbs.useExcelIntegration = false;

        // Calculate water state at T, P
        auto ws = waterState(T_K, P_Pa, wsOpts);

        // Extract water Gibbs energy and volume
        const double G0_H2O = ws.gibbs;       // J/mol
        // Convert specific volume (m³/kg) to molar volume (m³/mol)
        const double M_H2O = 0.018015;  // kg/mol
        // If V is not populated, calculate from density: V = 1/D
        const double V_specific = (ws.thermo.V != 0.0) ? ws.thermo.V : (1.0 / ws.thermo.D);
        const double V0_H2O = V_specific * M_H2O;    // m³/mol

        // Convert for comparison
        const double G0_H2O_cal = G0_H2O / CAL_TO_J;  // cal/mol
        const double V0_H2O_cm3 = V0_H2O / CM3_TO_M3;  // cm³/mol

        // Get standard thermo models for each solute species
        auto model_CO2 = CO2_aq.standardThermoModel();
        auto model_Hplus = H_plus.standardThermoModel();
        auto model_HCO3 = HCO3_minus.standardThermoModel();

        // Calculate standard properties at T, P
        StandardThermoProps props_CO2 = model_CO2(T_K, P_Pa);
        StandardThermoProps props_Hplus = model_Hplus(T_K, P_Pa);
        StandardThermoProps props_HCO3 = model_HCO3(T_K, P_Pa);

        // Convert to cal/mol and cm³/mol for comparison with CSV
        const double G0_CO2_cal = props_CO2.G0 / CAL_TO_J;
        const double G0_Hplus_cal = props_Hplus.G0 / CAL_TO_J;
        const double G0_HCO3_cal = props_HCO3.G0 / CAL_TO_J;

        const double V0_CO2_cm3 = props_CO2.V0 / CM3_TO_M3;
        const double V0_Hplus_cm3 = props_Hplus.V0 / CM3_TO_M3;
        const double V0_HCO3_cm3 = props_HCO3.V0 / CM3_TO_M3;

        // Reaction: H2O + CO2(aq) = H+ + HCO3-
        // ΔGr = G(H+) + G(HCO3-) - G(H2O) - G(CO2,aq)
        const double G_rxn_model = props_Hplus.G0 + props_HCO3.G0 - G0_H2O - props_CO2.G0;

        // ΔVr = V(H+) + V(HCO3-) - V(H2O) - V(CO2,aq)
        const double V_rxn_model_m3 = props_Hplus.V0 + props_HCO3.V0 - V0_H2O - props_CO2.V0;
        const double V_rxn_model = V_rxn_model_m3 / CM3_TO_M3;  // Convert m³/mol to cm³/mol

        // Calculate log K from ΔGr
        // log K = -ΔGr / (R*T*ln(10))
        const double log_K_model = -G_rxn_model / (R * T_K * std::log(10.0));

        // Calculate errors
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

        // Log per-point errors and use CHECK so all points are evaluated
        {
            std::ofstream ferr("reaction_errors.csv", std::ios::app);
            ferr << T_C << "," << P_kb << "," << abs_err_g << "," << rel_err_g << "," << abs_err_v << "," << rel_err_v << "," << abs_err_logk << "\n";
        }

        CHECK(almost_equal(G_rxn_model, G_rxn_truth, G_ABS_TOL, G_REL_TOL));
        CHECK(almost_equal(V_rxn_model, V_rxn_truth, V_ABS_TOL, V_REL_TOL));
        CHECK(almost_equal(log_K_model, log_K_truth, LOGK_ABS_TOL, LOGK_REL_TOL));
        if (
            almost_equal(G_rxn_model, G_rxn_truth, G_ABS_TOL, G_REL_TOL) &&
            almost_equal(V_rxn_model, V_rxn_truth, V_ABS_TOL, V_REL_TOL) &&
            almost_equal(log_K_model, log_K_truth, LOGK_ABS_TOL, LOGK_REL_TOL)
        ) {
            passed_count++;
        }
        // ...existing code...
    }

    // Print error statistics
    auto minmaxavg = [](const std::vector<double>& v) {
        double min = v.empty() ? 0.0 : *std::min_element(v.begin(), v.end());
        double max = v.empty() ? 0.0 : *std::max_element(v.begin(), v.end());
        double avg = v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        return std::make_tuple(min, max, avg);
    };
    auto print_stats = [](const char* label, const std::vector<double>& abs, const std::vector<double>& rel, const char* unit) {
        auto [min_abs, max_abs, avg_abs] = minmaxavg(abs);
        auto [min_rel, max_rel, avg_rel] = minmaxavg(rel);
        std::cout << label << " absolute error: min=" << min_abs << ", max=" << max_abs << ", avg=" << avg_abs << " " << unit << std::endl;
        std::cout << label << " relative error: min=" << min_rel*100 << "%, max=" << max_rel*100 << "%, avg=" << avg_rel*100 << "%" << std::endl;
    };
    std::cout << "\nTested " << test_count << " conditions, " << passed_count << " passed." << std::endl;
    print_stats("ΔGr", abs_err_G, rel_err_G, "J/mol");
    print_stats("ΔVr", abs_err_V, rel_err_V, "cm³/mol");
    print_stats("log K", abs_err_logK, rel_err_logK, "");
    REQUIRE(test_count > 0);
    REQUIRE(passed_count == test_count);
}
