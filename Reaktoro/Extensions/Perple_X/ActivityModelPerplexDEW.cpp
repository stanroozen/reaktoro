#include "ActivityModelPerplexDEW.hpp"

// C++ includes
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

// Define M_PI if not defined (Windows)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXElectrolyte.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXHKF.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXMrkPure.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>

namespace Reaktoro {

using std::log;
using std::pow;
using std::sqrt;

namespace {

auto envVarEnabled(const char* name, bool defaultValue) -> bool
{
    const char* value = std::getenv(name);
    if(!value || !*value)
        return defaultValue;

    std::string text(value);
    for(char& ch : text)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));

    if(text == "0" || text == "false" || text == "off" || text == "no")
        return false;
    if(text == "1" || text == "true" || text == "on" || text == "yes")
        return true;

    return defaultValue;
}

// ============================================================================
// Physical constants
// ============================================================================
constexpr double c_NA   = 6.02214076e23;      // Avogadro constant  (mol⁻¹)
constexpr double c_e    = 1.602176634e-19;    // Elementary charge  (C)
constexpr double c_kB   = 1.380649e-23;       // Boltzmann constant (J/K)
constexpr double c_R    = c_NA * c_kB;        // Gas constant (J/(mol·K)) = 8.31446...
constexpr double c_eps0 = 8.8541878128e-12;   // Vacuum permittivity (F/m)
const     double c_ln10 = log(10.0);

// ============================================================================
// Per-species Born correction parameters (extracted at construction time)
// ============================================================================

/// Born solvation parameters needed for the runtime ε-mixing correction.
///
/// Perple_X recomputes G° for each species at every Newton step using the
/// MIXED-FLUID dielectric constant ε_mix (from slvnt1/Looyenga), so the
/// Born term ω(1/ε_mix − 1) automatically captures fluid-composition effects.
///
/// In Reaktoro, StandardThermoModelPerplexDEW evaluates G° once with pure-
/// water ε.  The activity model must therefore apply the residual:
///
///   Δln_g[i] = [ω_i(gf_mix)×(1/ε_mix − 1) − ω_i(gf_pure)×(1/ε_pure − 1)] / RT
///
/// This struct stores the static (construction-time) data needed to evaluate
/// calculateBornOmega(charge, bornRadius, omega0, gf) at runtime.
struct BornCorrParams
{
    double omega0     = 0.0;                                  // J/mol, Born ω₀ at ref (T,P)
    double bornRadius = PerpleX::HKFConstants::NEUTRAL_RADIUS; // Å, effective radius
    double charge     = 0.0;                                  // ionic charge
    bool   active     = false; // true only for species loaded with StandardThermoModelPerplexDEW
};

/// Reference water molality (mol/kg) - inverse of Mwater in kg/mol
constexpr double mwref = 55.51;

// ============================================================================
// Perple_X fluid species index table (1-based, COH-Fluid+ species)
// CO2=2, CO=3, CH4=4, H2=5, H2S=6, SO2=8, N2=10, NH3=11, C2H6=16, HF=17, HCl=18
// ============================================================================
const std::map<std::string, int> FLUID_SPEC_IDX = {
    {"CO2", 2}, {"CO", 3}, {"CH4", 4}, {"H2", 5}, {"H2S", 6},
    {"SO2", 8}, {"N2", 10}, {"NH3", 11}, {"C2H6", 16}, {"HF", 17}, {"HCl", 18}};

// ============================================================================
// Effective ionic radii (Å) — for ExtendedDH mode (Helgeson et al. 1981)
// ============================================================================
const std::map<std::string, real> IONIC_RADII = {
    {"H+"  , 3.08}, {"Fe+++", 3.46},
    {"Li+" , 1.64}, {"Al+++", 3.33},
    {"Na+" , 1.91}, {"Au+++", 3.72},
    {"K+"  , 2.27}, {"La+++", 3.96},
    {"Rb+" , 2.41}, {"Gd+++", 3.79},
    {"Cs+" , 2.61}, {"In+++", 3.63},
    {"NH4+", 2.31}, {"Ca+++", 3.44},
    {"Ag+" , 2.20}, {"F-"   , 1.33},
    {"Au+" , 2.31}, {"Cl-"  , 1.81},
    {"Cu+" , 1.90}, {"Br-"  , 1.96},
    {"Mg++", 2.54}, {"I-"   , 2.20},
    {"Sr++", 3.00}, {"OH-"  , 1.40},
    {"Ca++", 2.87}, {"HS-"  , 1.84},
    {"Ba++", 3.22}, {"NO3-" , 2.81},
    {"Pb++", 3.08}, {"HCO3-", 2.10},
    {"Zn++", 2.62}, {"HSO4-", 2.37},
    {"Cu++", 2.60}, {"ClO4-", 3.59},
    {"Cd++", 2.85}, {"ReO4-", 4.23},
    {"Hg++", 2.98}, {"SO4--", 3.15},
    {"Fe++", 2.62}, {"CO3--", 2.81},
    {"Mn++", 2.68}};

auto effectiveIonicRadius(const Species& sp) -> real
{
    for(const auto& [name, r] : IONIC_RADII)
        if(isAlternativeChargedSpeciesName(sp.name(), name))
            return r;
    const auto z = sp.charge();
    if(z == -1) return 1.81;
    if(z == -2) return 3.00;
    if(z == -3) return 4.20;
    if(z == +1) return 2.31;
    if(z == +2) return 2.80;
    if(z == +3) return 3.60;
    if(z == +4) return 4.50;
    if(z <  -3) return -z * 4.2 / 3.0;
    return z * 4.5 / 4.0;
}

// ============================================================================
// Debye-Hückel A and B (corrected SI formulas, both params use same coulomb_term)
// ============================================================================

/// A = sqrt(2π NA ρw) × (e²/4πε0εkBT)^(3/2) / ln10  ≈ 0.509 at 25°C, 1 bar
auto debyeHuckelParamA(real T, real rho_w, real epsilon) -> real
{
    if(rho_w <= 0.0 || epsilon <= 0.0 || T <= 0.0) return 0.0;
    const auto denom  = 4.0 * M_PI * c_eps0 * epsilon * c_kB * T;
    const auto ct     = (c_e * c_e) / denom;  // e²/(4πε0εkBT) [m]
    return sqrt(2.0 * M_PI * c_NA * rho_w) * pow(ct, 1.5) / c_ln10;
}

/// B = sqrt(8π NA ρw × e²/4πε0εkBT)  ≈ 3.28×10⁹ m⁻¹ at 25°C, 1 bar
auto debyeHuckelParamB(real T, real rho_w, real epsilon) -> real
{
    if(rho_w <= 0.0 || epsilon <= 0.0 || T <= 0.0) return 0.0;
    const auto denom = 4.0 * M_PI * c_eps0 * epsilon * c_kB * T;
    const auto ct    = (c_e * c_e) / denom;
    return sqrt(8.0 * M_PI * c_NA * rho_w * ct);
}

auto getExtraDouble(const Map<String, Any>& extra, const String& key, double& value) -> bool
{
    const auto it = extra.find(key);
    if(it == extra.end()) return false;
    if(const auto v = std::any_cast<double>(&it->second)) { value = *v; return true; }
    if(const auto v = std::any_cast<float>(&it->second))  { value = *v; return true; }
    if(const auto v = std::any_cast<int>(&it->second))    { value = static_cast<double>(*v); return true; }
    if(const auto v = std::any_cast<std::uint64_t>(&it->second)) { value = static_cast<double>(*v); return true; }
    return false;
}

auto getExtraStateId(const Map<String, Any>& extra, const String& key, std::uint64_t& value) -> bool
{
    const auto it = extra.find(key);
    if(it == extra.end()) return false;
    if(const auto v = std::any_cast<std::uint64_t>(&it->second)) { value = *v; return true; }
    if(const auto v = std::any_cast<Index>(&it->second))         { value = static_cast<std::uint64_t>(*v); return true; }
    if(const auto v = std::any_cast<int>(&it->second))           { value = static_cast<std::uint64_t>(*v); return true; }
    return false;
}

auto coupledFluidDebugEnabled() -> bool
{
    const char* flag = std::getenv("REAKTORO_PERPLEXDEW_DEBUG");
    return flag && std::string(flag) != "0";
}

// ============================================================================
// Mixed-fluid solvent state — Perple_X slvnt0 (pure H₂O) or slvnt1 (mixed)
// ============================================================================

/// Compute the solvent dielectric state consistent with Perple_X COH-Fluid+ GFSM.
///
/// Pure H₂O (no significant fluid co-solvents):
///   Uses slvnt0 path: ZD05 EOS + Fernandez/Sverjensky dielectric via epsh2o.
///
/// Mixed solvent (H₂O + CO₂, CH₄, H₂S, SO₂, H₂, CO, N₂, NH₃, HF, C₂H₆, HCl):
///   Uses slvnt1 path: exact mirror of the Perple_X slvnt1 subroutine.
///
///   Perple_X slvnt1 builds hybrid volumes as:
///     vhyb(i) = dvhy(i) + v_mrkmix(i)
///   where:
///     - v_mrkmix(i) = MRK MIXTURE partial molar volume (from mrkmix/ghybrid)
///     - dvhy(H2O)   = v_ZD05 − v_MRK_pure  (hybrid EOS correction for water)
///     - dvhy(other) = 0                     (no substitution for other species)
///
///   So:  vhyb(H2O)   = v_ZD05 + (v_MRK_mix − v_MRK_pure)
///        vhyb(other) = v_MRK_mix  (mixture PMV)
///
///   The volume fractions from these vhyb values enter the Looyenga mixing rule
///   for ε, and then adh = CDH × sqrt(10 msol/vsolv / (εT)³).
///
auto computeFluidSolventState(
    double Pbar,
    double T,
    double& waterVolume,
    const Vec<Index>& neutralIdx,
    const Vec<int>&   neutralPidx,
    const ArrayXr& m,
    Index iwater) -> PerpleX::DielectricState
{
    // Collect fluid co-solvent species with significant molality
    std::vector<int> fluid_species;
    const double m_h2o  = static_cast<double>(m[iwater]);  // ~55.51 mol/kg
    double m_total      = m_h2o;

    for(std::size_t k = 0; k < neutralIdx.size(); ++k)
    {
        const int pidx = neutralPidx[k];
        if(pidx <= 0) continue;
        const double mk = static_cast<double>(m[neutralIdx[k]]);
        if(mk > 1.0e-10 * m_h2o)
        {
            fluid_species.push_back(pidx);
            m_total += mk;
        }
    }

    // --- slvnt0: pure water ---------------------------------------------------
    if(fluid_species.empty())
    {
        auto pureState = PerpleX::getWaterSolventState(Pbar, T, waterVolume);
        // Correction is zero by construction (eps_mix = eps_pure); store refs for uniformity.
        pureState.epsilon_pure_water = pureState.epsilon;
        pureState.gf_pure_water      = pureState.gf;
        return pureState;
    }

    // --- slvnt1: mixed-fluid path (strict Perple_X slvnt1 equivalence) ---------
    // 1. Mole fractions (yf) for all active solvent species (including H₂O)
    std::array<double, 19> yf{};
    fluid_species.push_back(1);   // include H₂O (index 1)
    yf[1] = m_h2o / m_total;
    for(std::size_t k = 0; k < neutralIdx.size(); ++k)
    {
        const int pidx = neutralPidx[k];
        if(pidx <= 0) continue;
        const double mk = static_cast<double>(m[neutralIdx[k]]);
        if(mk > 1.0e-10 * m_h2o)
            yf[pidx] = mk / m_total;
    }

    // 2. H₂O volume from ZD05 — also capture pure-water solvent state for Born correction.
    //    epsilon_pure_water and gf_pure_water are used in ActivityModelPerplexDEW to compute
    //    ΔG°_Born = ω(gf_mix)×(1/ε_mix − 1) − ω(gf_pure)×(1/ε_pure − 1).
    auto pureWater = PerpleX::getWaterSolventState(Pbar, T, waterVolume);  // fills waterVolume

    // 3. MRK MIXTURE partial molar volumes for all species (Perple_X mrkmix/ghybrid)
    //    Perple_X slvnt1 uses v_mrkmix(i) after the mrkmix call inside ghybrid.
    auto mix = PerpleX::mrkMix(fluid_species, yf, Pbar, T, {});

    // 4. Pure MRK for ALL active solvent species — needed for the hybrid dvhy correction
    //    for H₂O and for the MRK non-ideal excess activity of each co-solvent:
    //      dvhy(H2O) = v_ZD05 − v_MRK_pure,   vhyb(H2O) = v_ZD05 + (v_MRK_mix − v_MRK_pure)
    //      ln_f_excess[j] = mix.ln_f[j] − mrkPureAll.ln_f[j]  (Perple_X ghybrid component)
    auto mrkPureAll = PerpleX::mrkPure(fluid_species, Pbar, T);

    // 5. Hybrid volumes: exact Perple_X slvnt1 formula vhyb(i) = dvhy(i)+v_mrkmix(i)
    std::array<double, 19> vhyb{};
    vhyb[1] = waterVolume + (mix.v[1] - mrkPureAll.v[1]);  // ZD05 + dvhy correction
    for(const int pidx : fluid_species)
        if(pidx != 1) vhyb[pidx] = mix.v[pidx];         // mixture PMV, dvhy=0

    // 6. Mixed ε (Looyenga), adh, gf via computeSolventState (= Perple_X slvnt1)
    auto solvent = PerpleX::computeSolventState(yf, vhyb, fluid_species,
                                                static_cast<int>(fluid_species.size()),
                                                Pbar, T);

    // 7. GFSM water activity: ln(f_H2O_mix / f_H2O_pure) via MRK fugacity ratio.
    //
    //    In the GFSM hybrid EOS (Perple_X ifug=39 / COH-Fluid+), the fugacity
    //    of H2O in the mixed fluid is:
    //      f_H2O_mix = g_ZD05_pure * (g_MRK_mix / g_MRK_pure) * P * y_H2O
    //    The pure-H2O reference is:
    //      f_H2O_pure = g_ZD05_pure * P
    //
    //    The ZD05/ZD09 pure-species correction g_ZD05_pure is IDENTICAL in
    //    numerator and denominator and cancels exactly.  Only the MRK
    //    cross-species mixing term survives:
    //      ln(a_H2O) = ln(f_H2O_mix) - ln(f_H2O_pure)
    //               = mix.ln_f[1] - mrkH2O.ln_f[1]
    //
    //    This replicates the Perple_X slvnt1 + aqact water-activity pathway
    //    for any GFSM co-solvent mixture (CO2, CH4, H2S, SO2, H2, CO, N2,
    //    NH3, HF, C2H6, HCl) and any combination thereof.
    //    Correct in the ideal-mixing limit (ln a_H2O → ln y_H2O as g_mix→g_pure)
    //    and at the pure-water limit (y_H2O=1 → ln a_H2O = 0).
    if(yf[1] > 0.0)
        solvent.ln_a_water = mix.ln_f[1] - mrkPureAll.ln_f[1];

    // Store pure-water reference for Born ε-correction in the activity model lambda.
    solvent.epsilon_pure_water = pureWater.epsilon;
    solvent.gf_pure_water      = pureWater.gf;

    // MRK non-ideal excess activity for each GFSM co-solvent species.
    // ln_f_excess[j] = ln(phi_mix[j]) - ln(phi_pure[j])  (Perple_X ghybrid pathway).
    // Applied in the lambda to neutral aqueous species that map to a GFSM fluid index
    // but have no HKF standard thermo model (pure GFSM treatment, no double-counting).
    for(const int pidx : fluid_species)
    {
        if(pidx == 1) continue;  // water activity handled separately via ln_a_water
        solvent.ln_f_excess[pidx] = mix.ln_f[pidx] - mrkPureAll.ln_f[pidx];
    }

    return solvent;
}

} // namespace

// ============================================================================
// Activity model implementation
// ============================================================================

auto activityModelPerplexDEW(
    const SpeciesList& species,
    const ActivityModelParamsPerplexDEW& options) -> ActivityModel
{
    const auto dhModel = options.dhModel;
    const auto failOnConflictingStandardState = options.errorOnConflictingStandardState;
    const auto requireCoupledFluidHandoff = options.requireCoupledGFSMHandoff;
    AqueousMixture mixture(species);

    const auto num_species          = mixture.species().size();
    const auto num_charged_species  = mixture.charged().size();
    const auto num_neutral_species  = mixture.neutral().size();
    const auto icharged_species     = mixture.indicesCharged();
    const auto ineutral_species     = mixture.indicesNeutral();
    const auto iwater               = mixture.indexWater();

    // Pre-compute Perple_X fluid index for each neutral species
    // (formula string, e.g. CO2(aq) → formula "CO2" → index 2)
    Vec<int> neutral_pidx;
    neutral_pidx.reserve(num_neutral_species);
    for(Index idx : ineutral_species)
    {
        const auto fstr = mixture.species(idx).formula().str();
        auto it = FLUID_SPEC_IDX.find(fstr);
        neutral_pidx.push_back(it != FLUID_SPEC_IDX.end() ? it->second : -1);
    }

    // Pre-compute z² for all charged species (Davies and ExtendedDH)
    Vec<double> z2_vec;
    z2_vec.reserve(num_charged_species);
    for(Index idx : icharged_species)
        z2_vec.push_back(static_cast<double>(mixture.species(idx).charge())
                       * static_cast<double>(mixture.species(idx).charge()));

    // Pre-compute ionic radii for ExtendedDH mode
    Vec<real> radii_vec;
    radii_vec.reserve(num_charged_species);
    for(Index idx : icharged_species)
        radii_vec.push_back(effectiveIonicRadius(mixture.species(idx)));

    // Pre-compute Born correction params for each species.
    // Populated only for species whose standard thermo model was created with
    // StandardThermoModelPerplexDEW (stores omega0, bornRadius, charge in params()).
    // Water is excluded — its activity is handled separately.
    Vec<BornCorrParams> born_vec(num_species);
    for(Index i = 0; i < num_species; ++i)
    {
        if(i == iwater) continue;
        const auto& p = mixture.species(i).standardThermoModel().params();
        if(p.exists("PerplexDEW") && p.exists("omega0"))
        {
            born_vec[i].omega0     = p["omega0"].asFloat();
            born_vec[i].bornRadius = p["bornRadius"].asFloat();
            born_vec[i].charge     = p["charge"].asFloat();
            born_vec[i].active     = true;
        }
    }

    // -------------------------------------------------------------------------
    // Pre-compute GFSM co-solvent flag.
    //
    // A "GFSM solvent" species is a neutral species that:
    //   • maps to a GFSM fluid index (CO2=2, CH4=4, H2=5, H2S=6, SO2=8, …), AND
    //   • has NO PerplexDEW HKF standard thermo model
    //
    // In Perple_X Aq_yuri these species are SOLVENTS (not solutes).  Their
    // contribution to chemical potential comes from ghybrid/MRK, and their
    // contribution to water activity comes from ln_a_water (MRK fugacity ratio).
    //
    // As a result they must be:
    //   (a) excluded from m_tot used in the Raoult water-activity term
    //       (adding them would double-count their effect on a_H2O)
    //   (b) excluded from the total-solvent molality reference correction
    //       (their activity is on a mole-fraction basis, not molal)
    // -------------------------------------------------------------------------
    Vec<bool> is_gfsm_solvent(num_species, false);
    for(Index k = 0; k < num_neutral_species; ++k)
    {
        const int pidx = neutral_pidx[k];
        if(pidx <= 0) continue;          // not a GFSM fluid species
        const Index i = ineutral_species[k];
        if(born_vec[i].active) continue; // has HKF model — solute, not pure GFSM
        if(i == iwater) continue;
        is_gfsm_solvent[i] = true;
    }

    // -------------------------------------------------------------------------
    // Conflict check: warn when a species has BOTH a PerplexDEW HKF model AND
    // maps to a GFSM fluid co-solvent index (CO2=2, CO=3, CH4=4, H2=5, H2S=6,
    // SO2=7, N2=9, C2H6=15, HF=16, HCl=17, ...).
    //
    // Perple_X treats these as mutually exclusive:
    //   • GFSM solvent path  — mole-fraction reference, MRK non-ideal excess
    //   • HKF/DEW solute path — molal reference, Born ε-correction
    //
    // Applying both simultaneously double-counts the chemical potential.
    // Remove the HKF entry and rely on the GFSM activity, or keep the HKF
    // entry and set the species' molality directly without GFSM.
    // -------------------------------------------------------------------------
    for(Index k = 0; k < num_neutral_species; ++k)
    {
        const int pidx = neutral_pidx[k];
        if(pidx <= 0) continue;  // not a GFSM fluid species
        const Index i = ineutral_species[k];
        if(!born_vec[i].active) continue;  // no HKF model, no conflict
        auto message = str(
            "PerplexDEW: species '", mixture.species(i).name(), "' is both a "
            "GFSM fluid co-solvent (Perple_X index ", pidx, ") and has a PerplexDEW "
            "HKF standard thermo model. These use incompatible standard states "
            "(mole-fraction GFSM vs. molal HKF) and will double-count the excess "
            "chemical potential. "
            "Solution A: remove the HKF species entry - activity comes from GFSM "
            "(mixture fugacity ratio, Perple_X ghybrid equivalent). "
            "Solution B: keep HKF only and exclude the species from the GFSM "
            "co-solvent list. "
            "To hard-fail on this conflict, set "
            "ActivityModelParamsPerplexDEW.errorOnConflictingStandardState=True.");

        if(failOnConflictingStandardState)
            errorif(true, message);
        warningif(true, message);
    }

    // -------------------------------------------------------------------------
    // Unmapped GFSM coupling diagnostic.
    //
    // A neutral non-water species that has NO HKF model and whose formula does
    // NOT match any FLUID_SPEC_IDX key (neutral_pidx[k] == -1) may still be a
    // Perple_X GFSM co-solvent stored under a non-standard formula string.
    // Common cause: database formula field has an aqueous suffix appended
    // (e.g. "CO2_aq" instead of "CO2"), so the exact-key lookup misses.
    //
    // The check strips a set of known suffixes and retests.  If the stripped
    // formula IS in FLUID_SPEC_IDX the user almost certainly has a naming
    // mismatch: GFSM coupling is inactive and the species gets activity = 1.
    // Controlled by ActivityModelParamsPerplexDEW.warnOnUnmappedGFSMCoupling.
    // -------------------------------------------------------------------------
    if(options.warnOnUnmappedGFSMCoupling)
    {
        static const std::vector<std::string> AQUEOUS_SUFFIXES =
            { ",aq", "_aq", "(aq)", "-aq", ".aq", ",AQ", "_AQ", "(AQ)" };

        for(Index k = 0; k < num_neutral_species; ++k)
        {
            if(neutral_pidx[k] > 0) continue;      // already coupled — no issue
            const Index i = ineutral_species[k];
            if(i == iwater) continue;
            if(born_vec[i].active) continue;        // HKF solute — intentional
            const std::string fstr = mixture.species(i).formula().str();
            std::string base = fstr;
            for(const auto& suf : AQUEOUS_SUFFIXES)
            {
                if(fstr.size() > suf.size() &&
                   fstr.compare(fstr.size() - suf.size(), suf.size(), suf) == 0)
                {
                    base = fstr.substr(0, fstr.size() - suf.size());
                    break;
                }
            }
            if(base == fstr) continue;  // no suffix stripped — no mismatch
            if(!FLUID_SPEC_IDX.count(base)) continue;
            warningif(true,
                "PerplexDEW: neutral species '", mixture.species(i).name(), "' "
                "has formula '", fstr, "' which did not match any Perple_X GFSM "
                "fluid co-solvent key. Stripping the aqueous suffix gives '", base,
                "' which IS a known GFSM species (FLUID_SPEC_IDX index ",
                FLUID_SPEC_IDX.at(base), "). GFSM mixed-fluid coupling for this "
                "species is INACTIVE — its activity coefficient will be 1 (no MRK "
                "non-ideal correction). Fix: change the database formula entry to "
                "'", base, "' so that formula().str() returns the bare molecular "
                "formula. To suppress this warning, set "
                "ActivityModelParamsPerplexDEW.warnOnUnmappedGFSMCoupling=False.");
        }
    }

    auto stateptr   = std::make_shared<AqueousMixtureState>();
    auto mixtureptr = std::make_shared<AqueousMixture>(mixture);

    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        const auto& [T, P, x] = args;
        auto const& state = *stateptr = mixture.state(T, P, x);

        props.som = StateOfMatter::Liquid;
        props.extra["AqueousMixtureState"] = stateptr;
        props.extra["AqueousMixture"]      = mixtureptr;

        const auto& I  = state.Is;
        const auto& m  = state.m;
        const auto& ms = state.ms;

        // Total solute molality (used for water activity — Raoult term).
        // GFSM co-solvent species are EXCLUDED: their effect on water activity
        // is already captured by solvent.ln_a_water (MRK fugacity ratio).
        // Including them here would double-count their contribution to a_H2O.
        real m_tot = 0.0;
        for(Index i = 0; i < num_species; ++i)
        {
            if(i == iwater) continue;
            if(is_gfsm_solvent[i]) continue;  // Perple_X GFSM solvent, not solute
            m_tot += m[i];
        }

        // --- Solvent state: slvnt0 (pure H₂O) or slvnt1 (mixed fluid) ----------
        const auto Pbar = P / 1.0e5;
        double waterVolume = 0.0;
        auto solvent = computeFluidSolventState(
            static_cast<double>(Pbar), static_cast<double>(T),
            waterVolume, ineutral_species, neutral_pidx, m, iwater);

        std::string waterActivitySource = "aqueous-neutral-fallback";
        bool usedCoupledFluidHandoff = false;

        // Coupled-fluid-first handoff: consume gas-phase GFSM signal if fresh.
        std::uint64_t stateid_expected = 0;
        std::uint64_t stateid_handoff = 0;
        const bool haveExpectedState = getExtraStateId(props.extra, "Reaktoro::ChemicalProps::StateId", stateid_expected);
        const bool haveHandoffState = getExtraStateId(props.extra, "PerplexGFSM::WaterActivity::StateId", stateid_handoff);

        if(haveExpectedState && haveHandoffState && stateid_expected == stateid_handoff)
        {
            double ln_ratio = 0.0;
            if(getExtraDouble(props.extra, "PerplexGFSM::WaterActivity::ln_f_ratio_h2o", ln_ratio))
            {
                solvent.ln_a_water = ln_ratio;
                usedCoupledFluidHandoff = true;
                waterActivitySource = "coupled-fluid";

                double v = 0.0;
                if(getExtraDouble(props.extra, "PerplexGFSM::WaterActivity::x_h2o_fluid", v) && v > 0.0)
                    solvent.x_H2O = v;
                if(getExtraDouble(props.extra, "PerplexGFSM::WaterActivity::msol_mix", v) && v > 0.0)
                    solvent.msol = v;
                if(getExtraDouble(props.extra, "PerplexGFSM::WaterActivity::vsolv_mix", v) && v > 0.0)
                    solvent.vsolv = v;
                if(getExtraDouble(props.extra, "PerplexGFSM::WaterActivity::epsilon_mix", v) && v > 0.0)
                    solvent.epsilon = v;
                if(getExtraDouble(props.extra, "PerplexGFSM::WaterActivity::epsilon_pure", v) && v > 0.0)
                    solvent.epsilon_pure_water = v;
                if(getExtraDouble(props.extra, "PerplexGFSM::WaterActivity::gf_mix", v))
                    solvent.gf = v;
                if(getExtraDouble(props.extra, "PerplexGFSM::WaterActivity::gf_pure", v))
                    solvent.gf_pure_water = v;

                // Keep Davies and Extended-DH paths on the same coupled-fluid state.
                solvent.adh = PerpleX::debyeHuckel(solvent.msol, solvent.vsolv, solvent.epsilon, static_cast<double>(T));
            }
        }

        if(requireCoupledFluidHandoff && !usedCoupledFluidHandoff)
        {
            errorif(true,
                "PerplexDEW strict coupled-fluid mode failed: "
                "requireCoupledGFSMHandoff=True but no fresh GFSM handoff was consumed. "
                "Diagnostics: haveExpectedStateId=", haveExpectedState ? "true" : "false",
                ", haveHandoffStateId=", haveHandoffState ? "true" : "false",
                ", expectedStateId=", stateid_expected,
                ", handoffStateId=", stateid_handoff,
                ". Ensure a PerplexGFSM gas phase is present and evaluated in the same ChemicalProps "
                "state before PerplexDEW, and that props.extra contains "
                "PerplexGFSM::WaterActivity::StateId and ::ln_f_ratio_h2o for the current state.");
        }

        props.extra["PerplexDEW::WaterActivity::Source"] = waterActivitySource;
        props.extra["PerplexDEW::WaterActivity::CoupledHandoffUsed"] = usedCoupledFluidHandoff;

        // Total water activity used for both water species and hydrous-equivalent
        // DEW corrections: ln(a_H2O) = ln(X_H2O,Raoult) + ln(f_H2O_mix/f_H2O_pure).
        const auto ln_X_water = log(mwref / (mwref + m_tot));
        const auto ln_a_h2o_gfsm = solvent.ln_a_water;
        const auto ln_a_water_total = ln_X_water + ln_a_h2o_gfsm;
        props.extra["PerplexDEW::WaterActivity::ln_X_water"] = static_cast<double>(ln_X_water);
        props.extra["PerplexDEW::WaterActivity::ln_a_h2o_gfsm"] = static_cast<double>(ln_a_h2o_gfsm);
        props.extra["PerplexDEW::WaterActivity::ln_a_h2o_total"] = static_cast<double>(ln_a_water_total);
        props.extra["PerplexDEW::WaterActivity::Enabled"] = true;

        const auto sqrt_I = sqrt(I);

        // -----------------------------------------------------------------------
        // Total-solvent molality scale factor (Perple_X slvnt2 consistency).
        //
        // Perple_X defines molalities as mol/kg_solvent_total:
        //   mo(k) = pa(k) / msol   where msol = Σ_j pa_j × M_j   (slvnt2)
        //
        // Here pa(k) for a dilute solute ≈ n_k / n_total_fluid, and
        // msol ≈ x_H2O^fluid × M_H2O + x_CO2^fluid × M_CO2 + ...
        //
        // Reaktoro's AqueousMixture uses mol/kg_H2O (pure-water denominator):
        //   m_reaktoro(k) = n_k / (M_H2O × n_H2O)
        //
        // Exact conversion (dilute-solute limit):
        //   mo(k) = m_reaktoro(k) × (M_H2O × x_H2O^fluid / msol)
        //
        // where x_H2O^fluid = n_H2O / (n_H2O + n_CO2 + …) = solvent.x_H2O = yf[1].
        //
        // Note: the simpler M_H2O/msol is missing the x_H2O factor and is
        // wrong by ~30% at X_CO2 = 0.3.
        //
        // msol_scale < 1.0 for mixed fluids, = 1.0 for pure water (slvnt0).
        // -----------------------------------------------------------------------
        constexpr double c_Mwater    = 18.01528e-3;  // kg/mol
        const double msol_scale      = c_Mwater * solvent.x_H2O / static_cast<double>(solvent.msol);
        const auto   sqrt_I_eff      = sqrt_I * std::sqrt(msol_scale); // sqrt(I_perp)
        const auto   I_eff           = I * msol_scale;                  // I_perp

        // --- Activity coefficients -------------------------------------------
        if(dhModel == ActivityDHModel::Davies)
        {
            // Perple_X Davies formula (aqact / slvnt2):
            //   ln(γᵢ) = zᵢ² × (adh × √I/(1+√I) + 0.2×I)
            // adh is negative (CDH is negative in Perple_X).
            // For mixed fluid, adh uses Looyenga-mixed ε and compositional ρ.
            const auto adh = solvent.adh;

            for(Index i = 0; i < num_charged_species; ++i)
            {
                const auto isp = icharged_species[i];
                if(ms[i] == 0.0) { props.ln_g[isp] = 0.0; continue; }
                props.ln_g[isp] = z2_vec[i] * (adh * sqrt_I_eff / (1.0 + sqrt_I_eff) + 0.2 * I_eff);
            }

            // Neutral species: z² = 0 → no DH correction (Perple_X slvnt2)
            for(Index i = 0; i < num_neutral_species; ++i)
                props.ln_g[ineutral_species[i]] = 0.0;
        }
        else  // ActivityDHModel::ExtendedDH
        {
            // Extended Debye-Hückel with ionic radii (corrected SI formulas):
            //   log₁₀(γᵢ) = −A zᵢ² √I / (1 + aᵢ B √I)
            // A and B are computed from the mixed-solvent ε and ρ
            // (same ε path as Davies, so both models are consistent).
            //
            // vsolv: cm³/mol, msol: kg/mol  ➡  rho_w: kg/m³
            const auto rho_w  = solvent.vsolv > 0.0
                              ? (solvent.msol / (solvent.vsolv * 1.0e-6))
                              : 0.0;
            const auto epsilon = solvent.epsilon;

            const auto A = debyeHuckelParamA(T, rho_w, epsilon);
            const auto B = debyeHuckelParamB(T, rho_w, epsilon);

            for(Index i = 0; i < num_charged_species; ++i)
            {
                const auto isp  = icharged_species[i];
                const auto z2   = z2_vec[i];
                const auto a_m  = static_cast<double>(radii_vec[i]) * 1.0e-10; // Å → m
                if(ms[i] == 0.0) { props.ln_g[isp] = 0.0; continue; }
                const auto lambda = 1.0 + a_m * B * sqrt_I_eff;
                props.ln_g[isp]   = -(A * z2 * sqrt_I_eff / lambda) * c_ln10;
            }

            // Neutral species: no DH correction
            for(Index i = 0; i < num_neutral_species; ++i)
                props.ln_g[ineutral_species[i]] = 0.0;
        }

        // Activities on molal scale (water-basis, mol/kg_H2O)
        props.ln_a = props.ln_g + m.log();

        // -----------------------------------------------------------------------
        // Total-solvent molality reference correction (Perple_X slvnt2 equivalent).
        //
        // Perple_X computes solute G as  G°_i + RT×ln(mo_i) + RT×lng0×z²
        // where mo_i = m_i × (M_H2O / msol_mix).  Applying the ln correction here
        // is equivalent and keeps all G° values on the standard water-molal scale
        // (consistent with HKF database parameterisation in DEW24).
        //
        //   Δln_a_i = ln(M_H2O / msol_mix) = ln(msol_scale)   for all HKF solutes
        //
        // Zero for pure water (msol_scale = 1).  Negative for mixed fluids.
        // Applied to charged species and HKF-neutral species.
        // GFSM co-solvent species are excluded (mole-fraction reference, not molal).
        // -----------------------------------------------------------------------
        if(msol_scale < 1.0 - 1.0e-9)  // mixed fluid — skip for pure water
        {
            const double delta_lnm = std::log(msol_scale);  // < 0 for mixed fluid
            for(Index i = 0; i < num_species; ++i)
            {
                if(i == iwater) continue;
                if(is_gfsm_solvent[i]) continue;  // mole-fraction reference, no correction
                props.ln_g[i] += delta_lnm;
                props.ln_a[i] += delta_lnm;
            }
        }

        // -----------------------------------------------------------------------
        // Born ε-mixing correction (Perple_X ghkf/slvnt1 equivalent).
        //
        // Perple_X recomputes G° for every species at every Newton step using the
        // MIXED-FLUID ε_mix (from Looyenga rule over H₂O+CO₂+…).  In Reaktoro,
        // StandardThermoModelPerplexDEW fixes ε = ε_pure_water.  The residual is:
        //
        //   ΔG°[i] = ω_i(gf_mix) × (1/ε_mix − 1) − ω_i(gf_pure) × (1/ε_pure − 1)
        //
        // Applied as Δln_a[i] = Δln_g[i] = ΔG°[i] / (RT).
        //
        // For pure water (slvnt0 path) ε_mix = ε_pure, so the block is a no-op.
        // For mixed fluids (e.g. H₂O+CO₂), ε_mix < ε_pure, and for SiO₂,aq
        // (ω₀ = 150,624 J/mol, neutral) this gives ~−0.05 to −0.3 in ln_a per
        // unit XCO₂, naturally reducing Si solubility without post-processing.
        // -----------------------------------------------------------------------
        const double eps_mix  = solvent.epsilon;
        const double eps_pure = solvent.epsilon_pure_water;

        if(eps_mix < eps_pure - 1.0e-8)  // non-trivial mixing: slvnt1 path active
        {
            const double gf_mix  = solvent.gf;
            const double gf_pure = solvent.gf_pure_water;
            const double RT      = c_R * static_cast<double>(T);

            for(Index i = 0; i < num_species; ++i)
            {
                if(i == iwater) continue;
                const auto& bp = born_vec[i];
                if(!bp.active) continue;

                const double w_mix  = PerpleX::calculateBornOmega(bp.charge, bp.bornRadius, bp.omega0, gf_mix);
                const double w_pure = PerpleX::calculateBornOmega(bp.charge, bp.bornRadius, bp.omega0, gf_pure);
                const double dG     = w_mix * (1.0/eps_mix - 1.0) - w_pure * (1.0/eps_pure - 1.0);

                props.ln_g[i] += dG / RT;
                props.ln_a[i] += dG / RT;
            }
        }

        // -----------------------------------------------------------------------
        // GFSM non-ideal excess activity for fluid co-solvent species
        // (Perple_X ghybrid/slvnt1 equivalent for neutral solutes without HKF).
        //
        // For species that:
        //   (a) map to a GFSM fluid index (CO2=2, CO=3, CH4=4, H2=5, H2S=6,
        //       SO2=7, N2=9, C2H6=15, HF=16, HCl=17, ...), AND
        //   (b) do NOT have a PerplexDEW HKF standard thermo model
        //
        // the activity coefficient is the MRK mixture excess term:
        //
        //   Δln_g[i] = ln(φ_mix[j]) − ln(φ_pure[j])
        //            = mix.ln_f[j] − mrkPure.ln_f[j]
        //
        // This matches Perple_X slvnt1 where co-solvent species are "solvent"
        // components and their excess μ comes entirely from the MRK equation.
        //
        // Species that DO have an HKF model get their G° corrected instead via
        // the Born ε-correction above — the two paths are mutually exclusive
        // (a construction-time warning flags accidental overlap).
        // -----------------------------------------------------------------------
        for(Index k = 0; k < num_neutral_species; ++k)
        {
            const int pidx = neutral_pidx[k];
            if(pidx <= 0) continue;              // not a GFSM fluid species
            const Index i = ineutral_species[k];
            if(born_vec[i].active) continue;     // has HKF model — Born path used instead
            if(i == iwater) continue;

            const double lnf_ex = solvent.ln_f_excess[pidx];
            if(lnf_ex == 0.0) continue;          // pure water path (slvnt0) or species absent

            props.ln_g[i] += lnf_ex;
            props.ln_a[i] += lnf_ex;
        }

        // Water activity: Raoult's law solute dilution + GFSM fluid-mixing correction.        //
        // solvent.ln_a_water = ln(f_H2O_mix / f_H2O_pure) from the MRK fugacity
        // ratio computed in the slvnt1 path above.  It is zero for pure water
        // (slvnt0 path) and negative whenever any GFSM co-solvent (CO2, CH4,
        // H2S, SO2, H2, CO, N2, NH3, HF, C2H6, HCl) is dissolved in the system.
        //
        // The combined expression:
        //   ln a_H2O = ln(m_ref / (m_ref + Σm_i))   ← solute dilution (Raoult)
        //            + ln(f_H2O_mix / f_H2O_pure)    ← GFSM fluid-mixing
        // is the Perple_X aqact/slvnt1 water activity for a COH+ fluid.
        // It generalises to any mixture supported by the GFSM solution model.
        if(coupledFluidDebugEnabled())
        {
            static int debugLines = 0;
            if(debugLines < 200)
            {
                std::cout
                    << "[PerplexDEW] source=" << waterActivitySource
                    << " lnX=" << static_cast<double>(ln_X_water)
                    << " ln_f_ratio=" << static_cast<double>(solvent.ln_a_water)
                    << " ln_a_total=" << static_cast<double>(ln_a_water_total)
                    << " xH2O=" << solvent.x_H2O
                    << " epsMix=" << solvent.epsilon
                    << " epsPure=" << solvent.epsilon_pure_water
                    << "\n";
                ++debugLines;
            }
        }

        props.ln_a[iwater] = ln_a_water_total;
        props.ln_g[iwater] = ln_a_water_total;
    };

    return fn;
}

auto ActivityModelPerplexDEW(const ActivityModelParamsPerplexDEW& params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& species) { return activityModelPerplexDEW(species, params); };
}

auto ActivityModelPerplexDEW(ActivityDHModel model) -> ActivityModelGenerator
{
    ActivityModelParamsPerplexDEW params;
    params.dhModel = model;
    return ActivityModelPerplexDEW(params);
}

} // namespace Reaktoro
