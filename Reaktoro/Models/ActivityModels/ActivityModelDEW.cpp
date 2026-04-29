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

#include "ActivityModelDEW.hpp"

// C++ includes
#include <cmath>
#include <map>

// Define M_PI if not defined (Windows)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Extensions/DEW/WaterState.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>

namespace Reaktoro {

using std::abs;
using std::log;
using std::log10;
using std::pow;
using std::sqrt;

namespace {

/// Physical and mathematical constants
const auto NA = 6.02214076e23;          // Avogadro's number (mol⁻¹)
const auto e = 1.602176634e-19;         // Elementary charge (C)
const auto kB = 1.380649e-23;           // Boltzmann constant (J/K)
const auto eps0 = 8.8541878128e-12;     // Vacuum permittivity (F/m)
const auto ln10 = log(10.0);            // Natural log of 10

/// Molar mass of water (kg/mol)
const auto Mwater = 0.018015;

/// Reference molar fraction of water at standard state
const auto mwref = 55.51;               // mol/kg at 25°C, 1 bar

/// Finite-water correction constant (DEW, Eq. 39)
const auto Cc_constant = 0.0180153;

/// Effective electrostatic radii of ionic species (Å, from Helgeson et al. 1981)
const std::map<std::string, real> effective_radii = {
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
    {"Mn++", 2.68}
};

/// Calculate the effective electrostatic radius of a species (in Ångströms)
auto effectiveIonicRadius(const Species& species) -> real
{
    // Try to find exact match in effective_radii table
    for(auto pair : effective_radii)
        if(isAlternativeChargedSpeciesName(species.name(), pair.first))
            return pair.second;

    // If not found, estimate based on charge (TOUGHREACT approach)
    const auto z = species.charge();

    if(z == -1) return 1.81;        // Based on Cl-
    if(z == -2) return 3.00;        // Based on CO3-- and SO4-- average
    if(z == -3) return 4.20;        // Estimation
    if(z == +1) return 2.31;        // Based on NH4+
    if(z == +2) return 2.80;        // Average of +2 species
    if(z == +3) return 3.60;        // Average of +3 species
    if(z == +4) return 4.50;        // Estimation
    if(z <  -3) return -z*4.2/3.0;  // Linear extrapolation
    return z*4.5/4.0;               // Linear extrapolation
}

/// Calculate Debye-Hückel parameter A(T,P) dynamically from DEW water properties
/// Formula: A = sqrt(2π NA ρw) × (e²/(4π ε₀ ε kB T))^(3/2) / ln10
/// Gives A ≈ 0.509 (mol/kg)^(-1/2) at 25°C, 1 bar (verified against literature)
auto debyeHuckelParamA(real T, real P, real rho_w, real epsilon) -> real
{
    if(rho_w <= 0.0 || epsilon <= 0.0 || T <= 0.0)
        return 0.0;

    // coulomb_term = e² / (4π ε₀ ε kB T)  [m]
    const auto denom = 4.0 * M_PI * eps0 * epsilon * kB * T;
    const auto coulomb_term = (e * e) / denom;

    // A = sqrt(2π NA ρw) × coulomb_term^(3/2) / ln10
    // rho_w in kg/m³ → consistent SI units throughout
    const auto A = sqrt(2.0 * M_PI * NA * rho_w) * pow(coulomb_term, 1.5) / ln10;

    return A;
}

/// Calculate Debye-Hückel parameter B(T,P) dynamically from DEW water properties
/// Formula: B = sqrt(8π NA ρw × e²/(4π ε₀ ε kB T)) = sqrt(2 NA ρw e²/(ε₀ ε kB T))
/// Gives B ≈ 3.28×10⁹ m⁻¹ (mol/kg)^(-1/2) at 25°C, 1 bar (verified against literature)
auto debyeHuckelParamB(real T, real P, real rho_w, real epsilon) -> real
{
    if(rho_w <= 0.0 || epsilon <= 0.0 || T <= 0.0)
        return 0.0;

    // coulomb_term = e² / (4π ε₀ ε kB T)  [m]  (shared factor with A)
    const auto denom = 4.0 * M_PI * eps0 * epsilon * kB * T;
    const auto coulomb_term = (e * e) / denom;

    // B = sqrt(8π NA ρw × coulomb_term)  [m⁻¹ (mol/kg)^(-1/2)]
    const auto B = sqrt(8.0 * M_PI * NA * rho_w * coulomb_term);

    return B;
}

} // namespace

auto activityModelDEW(const SpeciesList& species) -> ActivityModel
{
    // Create the aqueous mixture
    AqueousMixture mixture(species);

    // The number of species in the mixture
    const auto num_species = mixture.species().size();

    // The number of charged and neutral species
    const auto num_charged_species = mixture.charged().size();
    const auto num_neutral_species = mixture.neutral().size();

    // The indices of the charged and neutral species
    const auto icharged_species = mixture.indicesCharged();
    const auto ineutral_species = mixture.indicesNeutral();

    // The index of the water species
    const auto iwater = mixture.indexWater();

    // Extract effective ionic radii for all charged species
    Vec<real> effective_radii_vec;
    Vec<double> charges;

    for(Index idx_ion : icharged_species)
    {
        const Species& sp = mixture.species(idx_ion);
        effective_radii_vec.push_back(effectiveIonicRadius(sp));
        charges.push_back(sp.charge());
    }

    // Create shared pointers for caching state
    auto stateptr = std::make_shared<AqueousMixtureState>();
    auto mixtureptr = std::make_shared<AqueousMixture>(mixture);

    // Define the activity model function
    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // Extract T, P, composition
        const auto& [T, P, x] = args;

        // Evaluate the state of the aqueous mixture
        auto const& state = *stateptr = mixture.state(T, P, x);

        // Set the state of matter
        props.som = StateOfMatter::Liquid;

        // Export the aqueous mixture state and mixture via extra data
        props.extra["AqueousMixtureState"] = stateptr;
        props.extra["AqueousMixture"] = mixtureptr;

        // Get state variables
        const auto& I = state.Is;      // Ionic strength
        const auto& m = state.m;       // Molalities
        const auto& ms = state.ms;     // Stoichiometric molalities of charged species

        // ========== STEP 1: Compute total molality ==========
        real m_tot = 0.0;
        for(Index i = 0; i < num_species; ++i)
        {
            if(i != iwater)
                m_tot += m[i];
        }

        // ========== STEP 2: Retrieve water properties from DEW ==========
        // Get water state using ZD05 EOS + PowerFunction dielectric —
        // consistent with StandardThermoModelDEW defaults.
        WaterStateOptions opts;
        opts.thermo.eosModel    = WaterEosModel::ZhangDuan2005;
        opts.dielectric.primary = WaterDielectricPrimaryModel::PowerFunction;
        const auto water_state = waterState(T, P, opts);

        const auto rho_w = water_state.thermo.D;           // Water density (kg/m³)
        const auto epsilon = water_state.electro.epsilon;  // Dielectric constant

        // ========== STEP 3: Compute Debye-Hückel parameters A and B ==========
        const auto A = debyeHuckelParamA(T, P, rho_w, epsilon);
        const auto B = debyeHuckelParamB(T, P, rho_w, epsilon);

        // ========== STEP 4: Compute finite-water correction ==========
        const auto Cc = Cc_constant * log10(1.0 + m_tot);

        // ========== STEP 5: Calculate ionic strength square root ==========
        const auto sqrt_I = sqrt(I);

        // ========== STEP 6: Calculate activity coefficients for charged species ==========
        for(Index i = 0; i < num_charged_species; ++i)
        {
            const auto ispecies = icharged_species[i];
            const auto z = charges[i];
            const auto z2 = z * z;
            const auto a = effective_radii_vec[i] * 1e-10;  // Convert Å to m for consistency

            // Skip if molality is zero
            if(ms[i] == 0.0)
            {
                props.ln_g[ispecies] = 0.0;
                continue;
            }

            // ========== HKF/DEW EQUATION (Eq. 39 in Huang & Sverjensky 2019) ==========
            // log₁₀(γⱼ) = -(A*zⱼ²*√I) / (1 + aⱼ*B*√I) + bᶜ'ᵏ*I + Cᶜ

            // DEW assumption: bᶜ'ᵏ = 0 at deep-Earth conditions
            const auto b_extended = 0.0;  // Set to zero for DEW

            // Lambda factor
            const auto lambda = 1.0 + a * B * sqrt_I;

            // Debye-Hückel term
            const auto dh_term = -(A * z2 * sqrt_I) / lambda;

            // Total activity coefficient (log₁₀ scale)
            const auto log10_gamma = dh_term + b_extended * I + Cc;

            // Convert to natural log scale
            props.ln_g[ispecies] = log10_gamma * ln10;
        }

        // ========== STEP 7: Calculate activity coefficients for neutral species ==========
        for(Index i = 0; i < num_neutral_species; ++i)
        {
            const auto ispecies = ineutral_species[i];

            // For most neutrals: log₁₀(γₙ) = Cᶜ
            // (Special cases like CO₂, CH₄ would require parameterized bᶜ'ⁿ, omitted for brevity)
            const auto log10_gamma = Cc;

            props.ln_g[ispecies] = log10_gamma * ln10;
        }

        // ========== STEP 8: Calculate activities of solutes (molal scale) ==========
        props.ln_a = props.ln_g + m.log();

        // ========== STEP 9: Calculate activity of water (ideal mixing on molal scale) ==========
        // Water mole fraction: X_H2O = 55.51 / (55.51 + Σmᵢ)
        const auto X_water = mwref / (mwref + m_tot);
        const auto ln_X_water = log(X_water);

        props.ln_a[iwater] = ln_X_water;
        props.ln_g[iwater] = ln_X_water;  // For ideal mixing: γ_H2O ≈ 1 on molal scale
    };

    return fn;
}

auto ActivityModelDEW() -> ActivityModelGenerator
{
    return [](const SpeciesList& species) { return activityModelDEW(species); };
}

} // namespace Reaktoro
