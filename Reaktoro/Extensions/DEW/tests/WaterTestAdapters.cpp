// WaterTestAdapters.cpp
//
// Adapter layer between your DEW water implementation and the regression
// tests defined in test_dew_water.cpp. All unit conversions and database
// lookup for aqueous species are handled here.

#include "WaterTestAdapters.hpp"

#include <cmath>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <algorithm>

#include <yaml-cpp/yaml.h>

// ----- Reaktoro / DEW headers (adjust paths as needed) ----------------

// Core types / options
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoProps.hpp>
#include <Reaktoro/Extensions/DEW/WaterThermoModel.hpp>
#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>

// EOS
#include <Reaktoro/Extensions/DEW/WaterEosZhangDuan2005.hpp>
#include <Reaktoro/Extensions/DEW/WaterEosZhangDuan2009.hpp>

// Psat polynomials
#include <Reaktoro/Extensions/DEW/WaterPsatPolynomialsDEW.hpp>

// Dielectric models and electro props
#include <Reaktoro/Extensions/DEW/WaterDielectricJohnsonNorton.hpp>
#include <Reaktoro/Extensions/DEW/WaterDielectricFranck1990.hpp>
#include <Reaktoro/Extensions/DEW/WaterDielectricFernandez1997.hpp>
#include <Reaktoro/Extensions/DEW/WaterDielectricPowerFunction.hpp>
#include <Reaktoro/Extensions/DEW/WaterElectroProps.hpp>

// Solvent function g and derivatives
#include <Reaktoro/Extensions/DEW/WaterSolventFunctionDEW.hpp>

// Gibbs free energy models
#include <Reaktoro/Extensions/DEW/WaterGibbsModel.hpp>

// Born Omega and dOmega/dP
#include <Reaktoro/Extensions/DEW/WaterBornOmegaDEW.hpp>

// ----------------------------------------------------------------------

using Reaktoro::real;
using Reaktoro::WaterThermoProps;
using Reaktoro::WaterThermoModelOptions;
using Reaktoro::WaterSolventFunctionOptions;
using Reaktoro::WaterGibbsModelOptions;
using Reaktoro::WaterGibbsModel;
using Reaktoro::WaterEosModel;
using Reaktoro::WaterElectroProps;
using Reaktoro::WaterBornOmegaOptions;

namespace {

// --------------------------- constants & helpers ---------------------------

constexpr double J_per_cal      = 4.184;
constexpr double cal_per_J      = 1.0 / J_per_cal;
constexpr double Pa_per_bar     = 1.0e5;
constexpr double bar_per_Pa     = 1.0e-5;
constexpr double kgm3_to_gcm3   = 1.0e-3;
constexpr double gcm3_to_kgm3   = 1.0e3;

inline double toKelvin(double T_C)     { return T_C + 273.15; }
inline double toPa(double P_bar)       { return P_bar * Pa_per_bar; }
inline double rho_si_to_gcm3(double r) { return r * kgm3_to_gcm3; }
inline double rho_gcm3_to_si(double r) { return r * gcm3_to_kgm3; }

// ------------------- Water EOS convenience wrappers -------------------

// “Equation 1” density: ZD2005 EOS
inline WaterThermoProps thermo_ZD2005(double T_K, double P_Pa)
{
    return Reaktoro::waterThermoPropsZhangDuan2005(T_K, P_Pa);
}

// “Equation 2” density: ZD2009 EOS
inline WaterThermoProps thermo_ZD2009(double T_K, double P_Pa)
{
    return Reaktoro::waterThermoPropsZhangDuan2009(T_K, P_Pa);
}

// ZD2009 via unified model (for dielectric, g, etc.)
inline WaterThermoProps thermo_model_ZD2009(double T_K, double P_Pa)
{
    WaterThermoModelOptions opt;
    opt.eosModel               = WaterEosModel::ZhangDuan2009;
    opt.usePsatPolynomials     = false;
    opt.psatRelativeTolerance  = 1e-3;
    return Reaktoro::waterThermoPropsModel(T_K, P_Pa, opt);
}

// ZD2005 via unified model (for Q with densEq1)
inline WaterThermoProps thermo_model_ZD2005(double T_K, double P_Pa)
{
    WaterThermoModelOptions opt;
    opt.eosModel               = WaterEosModel::ZhangDuan2005;
    opt.usePsatPolynomials     = false;
    opt.psatRelativeTolerance  = 1e-3;
    return Reaktoro::waterThermoPropsModel(T_K, P_Pa, opt);
}

// ------------------------- depsdrho helper thermo --------------------------
//
// We want (dε/dρ)_T in units matching your truth tables, not ε_P.
//
// In the dielectric modules, ε_P is computed as:
//
//    ε_P = (dε/dρ_g) * (dρ_g/dP_SI),
//
// with ρ_g in g/cm³ and P in Pa. But your truth tables for depsdrho_* give
// directly dε/dρ_g.
//
// If we construct a fake thermo state with:
//
//    wt.D  = ρ_SI  = ρ_gcm3 * 1000  [kg/m³]
//    wt.DP = 1000  [kg/m³/Pa],
//
// then:
//
//    dρ_g/dP_SI = (1/1000) * (dρ_SI/dP_SI) = wt.DP / 1000 = 1,
//    so ε_P = dε/dρ_g exactly.
//
// That’s what makeThermoForDepsdrho does.
//

inline WaterThermoProps makeThermoForDepsdrho(double rho_g_cm3)
{
    WaterThermoProps wt;
    wt.D  = rho_gcm3_to_si(rho_g_cm3); // kg/m³
    wt.DP = 1000.0;                    // so that dρ_g/dP_SI = 1
    return wt;
}

// ---------------------- Omega species DB lookup (DEW) ----------------------

struct SpeciesBornParams
{
    double Z;             // charge
    double wref_J_mol;    // wref in J/mol (unscaled)
    bool   isHydrogenLike;
};

inline std::string strip_quotes(std::string s)
{
    if (!s.empty() && (s.front() == '"' || s.front() == '\''))
        s.erase(s.begin());
    if (!s.empty() && (s.back() == '"' || s.back() == '\''))
        s.pop_back();
    return s;
}

// Path to dew2019-aqueous.yaml from embedded DEW DB.
// Define DEW_AQUEOUS_DB_PATH in your CMake to hard-wire a path.
inline std::string dewAqueousDbPath()
{
#ifdef DEW_AQUEOUS_DB_PATH
    return std::string(DEW_AQUEOUS_DB_PATH);
#else
    // Fallback: assume YAML is next to the binary / in cwd.
    return "dew2019-aqueous.yaml";
#endif
}

const std::unordered_map<std::string, SpeciesBornParams>& speciesBornTable()
{
    static std::unordered_map<std::string, SpeciesBornParams> table;
    static std::once_flag init_flag;

    std::call_once(init_flag, [](){
        const std::string path = dewAqueousDbPath();
        YAML::Node root = YAML::LoadFile(path);

        YAML::Node speciesNode = root["Species"];
        if (!speciesNode || !speciesNode.IsMap()) {
            throw std::runtime_error("dew2024-aqueous.yaml: missing or invalid 'Species' node");
        }

        for (auto it = speciesNode.begin(); it != speciesNode.end(); ++it) {
            const YAML::Node& spec = it->second;

            if (!spec["Name"] || !spec["Charge"])
                continue;

            std::string name = spec["Name"].as<std::string>(); // e.g. ACETATE,AQ
            double      Z    = spec["Charge"].as<double>();

            YAML::Node hkf = spec["StandardThermoModel"]["HKF"];
            if (!hkf || !hkf["wref"])
                continue;

            // HKF.wref in the YAML is provided in SI units (J/mol) by the
            // database generator; use it directly without additional scaling.
            double wref_param = hkf["wref"].as<double>();
            double wref_J_mol = wref_param;

            // Very simple heuristic for the hydrogen-like flag:
            bool isHydrogenLike = false;
            std::string upperName = name;
            std::transform(upperName.begin(), upperName.end(), upperName.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
            if (upperName == "H+,AQ" || upperName == "HYDROGEN-ION,AQ")
                isHydrogenLike = true;

            table.emplace(name, SpeciesBornParams{Z, wref_J_mol, isHydrogenLike});
        }

        if (table.empty()) {
            throw std::runtime_error("dew2024-aqueous.yaml: no species loaded for Born Omega");
        }
    });

    return table;
}

const SpeciesBornParams& getSpeciesParams(const std::string& csvNameRaw)
{
    std::string name = strip_quotes(csvNameRaw);
    const auto& table = speciesBornTable();
    auto it = table.find(name);
    if (it == table.end()) {
        throw std::runtime_error("Species not found in dew2024-aqueous.yaml: " + name);
    }
    return it->second;
}

} // end anonymous namespace

// ==========================================================================
// Density ρ
// ==========================================================================

double dew_density_ZD2005(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);

    WaterThermoProps wt = thermo_ZD2005(T_K, P_Pa);
    return rho_si_to_gcm3(wt.D); // kg/m³ → g/cm³
}

double dew_density_ZD2009(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);

    WaterThermoProps wt = thermo_ZD2009(T_K, P_Pa);
    return rho_si_to_gcm3(wt.D);
}

double dew_density_psat(double T_C)
{
    const double T_K = toKelvin(T_C);

    // You likely implemented this in WaterPsatPolynomialsDEW.cpp:
    // double waterPsatDensityDEW(double T_K);
    double rho_kg_m3 = Reaktoro::waterPsatDensityDEW(T_K);

    return rho_si_to_gcm3(rho_kg_m3);
}

// ==========================================================================
// dρ/dP in (g/cm³)/bar
// ==========================================================================
//
// WaterThermoProps.DP is in kg/m³/Pa. Thus:
//
//   dρ_g/dP_bar = dρ_SI/dP_Pa * (1e-3 g/kg) * (1e5 Pa/bar)
//               = DP * 100.
//

double dew_drhodP_ZD2005(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);

    WaterThermoProps wt = thermo_ZD2005(T_K, P_Pa);
    return wt.DP * 100.0;
}

double dew_drhodP_ZD2009(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);

    WaterThermoProps wt = thermo_ZD2009(T_K, P_Pa);
    return wt.DP * 100.0;
}

// ==========================================================================
// Dielectric ε(T,P)
// ==========================================================================

static WaterElectroProps eps_JN1991_props(double T_K, double P_Pa)
{
    WaterThermoProps wt = thermo_model_ZD2009(T_K, P_Pa);
    return Reaktoro::waterElectroPropsJohnsonNorton(T_K, P_Pa, wt);
}

static WaterElectroProps eps_Franck_props(double T_K, double P_Pa)
{
    WaterThermoProps wt = thermo_model_ZD2009(T_K, P_Pa);
    return Reaktoro::waterElectroPropsFranck1990(T_K, P_Pa, wt);
}

static WaterElectroProps eps_Fernandez_props(double T_K, double P_Pa)
{
    WaterThermoProps wt = thermo_model_ZD2009(T_K, P_Pa);
    return Reaktoro::waterElectroPropsFernandez1997(T_K, P_Pa, wt);
}

static WaterElectroProps eps_Power_props(double T_K, double P_Pa)
{
    WaterThermoProps wt = thermo_model_ZD2009(T_K, P_Pa);
    return Reaktoro::waterElectroPropsPowerFunction(T_K, P_Pa, wt);
}

double dew_epsilon_JN1991(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);
    return eps_JN1991_props(T_K, P_Pa).epsilon;
}

double dew_epsilon_Franck1990(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);
    return eps_Franck_props(T_K, P_Pa).epsilon;
}

double dew_epsilon_Fernandez1997(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);
    return eps_Fernandez_props(T_K, P_Pa).epsilon;
}

double dew_epsilon_Power(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);
    return eps_Power_props(T_K, P_Pa).epsilon;
}

double dew_epsilon_psat(double T_C)
{
    const double T_K = toKelvin(T_C);

    // Implemented in WaterPsatPolynomialsDEW:
    // double waterPsatEpsilonDEW(double T_K);
    return Reaktoro::waterPsatEpsilonDEW(T_K);
}

// ==========================================================================
// dε/dρ(T,ρ) for each dielectric model
// ==========================================================================
//
// We use makeThermoForDepsdrho to force ε_P = dε/dρ_g.
//

double dew_depsdrho_JN1991(double T_C, double rho_g_cm3)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(1000.0); // arbitrary

    WaterThermoProps wt = makeThermoForDepsdrho(rho_g_cm3);
    auto we = Reaktoro::waterElectroPropsJohnsonNorton(T_K, P_Pa, wt);
    return we.epsilonP;
}

double dew_depsdrho_Franck1990(double T_C, double rho_g_cm3)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(1000.0);

    WaterThermoProps wt = makeThermoForDepsdrho(rho_g_cm3);
    auto we = Reaktoro::waterElectroPropsFranck1990(T_K, P_Pa, wt);
    return we.epsilonP;
}

double dew_depsdrho_Fernandez1997(double T_C, double rho_g_cm3)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(1000.0);

    WaterThermoProps wt = makeThermoForDepsdrho(rho_g_cm3);
    auto we = Reaktoro::waterElectroPropsFernandez1997(T_K, P_Pa, wt);
    return we.epsilonP;
}

double dew_depsdrho_Power(double T_C, double rho_g_cm3)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(1000.0);

    WaterThermoProps wt = makeThermoForDepsdrho(rho_g_cm3);
    auto we = Reaktoro::waterElectroPropsPowerFunction(T_K, P_Pa, wt);
    return we.epsilonP;
}

// ==========================================================================
// Solvent function g(T,P) and d(g)/dP
// ==========================================================================

double dew_g_eq2(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);

    WaterThermoProps wt = thermo_model_ZD2009(T_K, P_Pa);

    WaterSolventFunctionOptions opt;
    opt.Psat = false;  // bulk P, eq2 branch

    // Implemented in WaterSolventFunctionDEW:
    // double waterSolventFunctionDEW(double T, double P, const WaterThermoProps&, const WaterSolventFunctionOptions&);
    double g = Reaktoro::waterSolventFunctionDEW(T_K, P_Pa, wt, opt);
    return g; // dimensionless
}

double dew_dgdP_eq2(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);

    WaterThermoProps wt = thermo_model_ZD2009(T_K, P_Pa);

    WaterSolventFunctionOptions opt;
    opt.Psat = false;

    double g      = Reaktoro::waterSolventFunctionDEW(T_K, P_Pa, wt, opt);
    double dgdP_SI = Reaktoro::waterSolventFunctionDgdP_DEW(T_K, P_Pa, wt, g, opt);

    // Truth uses 1/Pa, so no conversion.
    return dgdP_SI;
}

double dew_dgdP_psat(double T_C)
{
    const double T_K = toKelvin(T_C);

    // Implemented in WaterPsatPolynomialsDEW:
    // double waterPsatDgdPDEW(double T_K);
    double dgdP_A_per_bar = Reaktoro::waterPsatDgdPDEW(T_K);

    return dgdP_A_per_bar; // Å/bar
}

// ==========================================================================
// Gibbs free energy of water, in cal/mol
// ==========================================================================

double dew_G_DH1978(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);

    WaterGibbsModelOptions opt;
    opt.model              = WaterGibbsModel::DelaneyHelgeson1978;
    opt.usePsatPolynomials = false;

    double G_J_per_mol = Reaktoro::waterGibbsModel(T_K, P_Pa, opt);
    return G_J_per_mol * cal_per_J;
}

double dew_G_integral(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);

    WaterGibbsModelOptions opt;
    opt.model              = WaterGibbsModel::DewIntegral;
    opt.usePsatPolynomials = false;

    // Excel uses densityEquation=1 (ZD2005) for the integral calculation!
    // From Export_Gibbs_All: g = calculateGibbsOfWater(P, TT, 2, 1, False)
    // where equation=2 (integral method), densityEquation=1 (ZD2005)
    opt.thermo.eosModel              = WaterEosModel::ZhangDuan2005;
    opt.thermo.usePsatPolynomials    = false;
    opt.thermo.psatRelativeTolerance = 1e-3;

    // NOTE: ZD2005 implementation has hardcoded tolerance of 0.01 bar,
    // while Excel uses 100 bar for speed. This may cause small differences.

    double G_J_per_mol = Reaktoro::waterGibbsModel(T_K, P_Pa, opt);
    return G_J_per_mol * cal_per_J;
}

double dew_G_psat(double T_C)
{
    const double T_K = toKelvin(T_C);

    // Implemented in WaterPsatPolynomialsDEW:
    // double waterPsatGibbsDEW(double T_K);
    double G_J_per_mol = Reaktoro::waterPsatGibbsDEW(T_K);
    return G_J_per_mol * cal_per_J;
}

// ==========================================================================
// Born Omega and dOmega/dP for aqueous species
// ==========================================================================
//
// Uses dew2024-aqueous.yaml for Z and wref. The truth CSV values for Z and
// wref are *not* used as inputs, only for regression comparison.
//

double dew_omega_species(const std::string& speciesNameFromCsv,
                         double T_C, double P_bar, double rho_g_cm3)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);

    const SpeciesBornParams& p = getSpeciesParams(speciesNameFromCsv);

    WaterThermoProps wt;
    wt.D  = rho_gcm3_to_si(rho_g_cm3);
    wt.DP = 0.0; // not needed for Ω itself

    WaterBornOmegaOptions opt;
    opt.isHydrogenLike = p.isHydrogenLike;

    double omega_J_mol =
        Reaktoro::waterBornOmegaDEW(T_K, P_Pa, wt,
                                    p.wref_J_mol, p.Z, opt);

    // Convert J/mol → cal/mol for truth_Omega_AllSpecies.csv
    return omega_J_mol * cal_per_J;
}

double dew_domegadP_species(const std::string& speciesNameFromCsv,
                            double T_C, double P_bar, double rho_g_cm3)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);

    const SpeciesBornParams& p = getSpeciesParams(speciesNameFromCsv);

    // CRITICAL: Need wt.DP for dgdP calculation!
    // Excel uses densityEquation=1 (ZD2005) for drhodP in calculate_dgdP
    WaterThermoProps wt = thermo_model_ZD2005(T_K, P_Pa);

    // The density from CSV should match (within tolerance), but use
    // the CSV density to stay consistent with test data
    wt.D = rho_gcm3_to_si(rho_g_cm3);

    WaterBornOmegaOptions opt;
    opt.isHydrogenLike = p.isHydrogenLike;

    double domega_dP_J_per_mol_Pa =
        Reaktoro::waterBornDOmegaDP_DEW(T_K, P_Pa, wt,
                                        p.wref_J_mol, p.Z, opt);

    // Convert J/mol/Pa → cal/mol/bar:
    // J/mol/Pa → cal/mol/Pa: multiply by cal_per_J
    // cal/mol/Pa → cal/mol/bar: multiply by Pa_per_bar (since 1 bar = 100000 Pa)
    return domega_dP_J_per_mol_Pa * cal_per_J * Pa_per_bar;
}

// ==========================================================================
// Born Q(T,P): densEq1 (ZD2005) + epsEq4 = Power-law epsilon
// ==========================================================================
//
// DEW defines (dimensionless) Q as:
//   Q = (1 / ε²) * (∂ε/∂P)
//
// In your dielectric implementation, WaterElectroProps.bornQ is stored as
// ε_P / ε² with ε_P = (∂ε/∂P)_SI and P in Pa, so bornQ has units 1/Pa.
//
// The truth table truth_Q_densEq1_epsEq4.csv stores Q in 1/bar:
//
//   Q_bar_inv = bornQ * (Pa/bar) = bornQ * 1e5.
//

double dew_Q_densEq1_epsEq4(double T_C, double P_bar)
{
    const double T_K  = toKelvin(T_C);
    const double P_Pa = toPa(P_bar);

    // densEq1 = ZD2005 EOS via unified model:
    WaterThermoProps wt = thermo_model_ZD2005(T_K, P_Pa);

    // epsEq4 = Power-law dielectric model:
    WaterElectroProps we =
        Reaktoro::waterElectroPropsPowerFunction(T_K, P_Pa, wt);

    const double Q_1_per_Pa = we.bornQ;
    return Q_1_per_Pa * Pa_per_bar; // 1/Pa → 1/bar
}
