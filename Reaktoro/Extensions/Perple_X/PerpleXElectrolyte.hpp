#pragma once

#include <array>
#include <vector>

namespace Reaktoro::PerpleX {

// Dielectric constant parameterization for pure molecular species
// Based on Harvey & Lemmon (2005) and Harvey & Mountain (2017)
struct DielectricParams
{
    // H&L 2005 non-polar: P/rho = A + A_mu/T + B*rho + C*rho^D
    // H&M 2017 polar: rho*(C_alpha + g(rho,T)*A_mu(T)/T)
    // Expanded form: 11 parameters per species
    double a0, a1, A_mu, b0, b1, c0, c1, D;
    double C_alpha, g0, g1;

    // Flag for polar (true) vs non-polar (false)
    bool isPolar = false;
};

// Dielectric constant calculation state
struct DielectricState
{
    double epsilon = 78.47;  // Reference value at 25°C, 1 bar for pure water
    double adh = 0.0;        // Debye-Hückel factor (ln activity coeff)
    double gf = 0.0;         // Shock g-function for HKF model
    double msol = 1.0;       // Molar mass of solvent (kg/mol)
    double vsolv = 18.0;     // Volume of solvent (cm3/mol)
    /// Mole fraction of H₂O among all fluid solvent species (= yf[1]).
    /// = 1.0 for pure water (slvnt0), < 1 for mixed fluid (slvnt1).
    /// Used to compute the correct Perple_X total-solvent molality scale factor:
    ///   msol_scale = M_H2O × x_H2O / msol
    /// rather than the simpler but incorrect M_H2O / msol.
    double x_H2O = 1.0;
    double hyvol = 0.0;      // Total hybrid volume (cm3/mol)

    /// ln(a_H2O) from GFSM fluid mixing, computed as MRK fugacity ratio
    /// ln(f_H2O_mix / f_H2O_pure) = mix.ln_f[1] - mrkPure.ln_f[1].
    /// The hybrid (ZD05/ZD09) pure-species correction cancels in the ratio;
    /// only the MRK cross-species mixing term contributes.
    /// Zero for pure water (slvnt0 path). Negative when co-solvents are present.
    /// Generalises to any GFSM mixture: CO2, CH4, H2S, SO2, H2, CO,
    /// N2, NH3, HF, C2H6, HCl — or any subset thereof.
    double ln_a_water = 0.0;

    /// Pure-water dielectric constant at the same (T, P).
    /// Equal to epsilon for slvnt0 (pure water). For slvnt1 (mixed fluid),
    /// this is ε_pure(H₂O) captured before Looyenga mixing, used as the
    /// reference for the Born ε-correction in ActivityModelPerplexDEW.
    double epsilon_pure_water = 78.47;

    /// Shock g-function for pure water at the same (T, P).
    /// Equal to gf for slvnt0. For slvnt1, this is gf from the pure-water
    /// ZD05 solvent state, used as the Born ε-correction reference.
    double gf_pure_water = 0.0;

    /// MRK non-ideal excess activity for each GFSM fluid co-solvent species.
    ///
    /// Index = Perple_X fluid species index (1=H2O, 2=CO2, 3=CO, 4=CH4, 5=H2,
    /// 6=H2S, 7=SO2, 9=N2, 10=NH3, 15=C2H6, 16=HF, 17=HCl).
    ///
    /// ln_f_excess[j] = mix.ln_f[j] - mrkPure.ln_f[j]
    ///
    /// This is the excess chemical potential of species j in the mixed fluid
    /// relative to its pure-species MRK reference, exactly replicating the
    /// Perple_X ghybrid/slvnt1 pathway.
    ///
    /// Applied as Δln_g[j] for neutral aqueous species that map to a GFSM
    /// fluid index but have NO HKF standard thermo model (e.g. CO2(aq) without
    /// a DEW entry). Zero for slvnt0 (pure water) and for H2O (index 1).
    std::array<double, 19> ln_f_excess{};  // zero-initialized
};

/// Calculate dielectric constant of pure water (Sverjensky 2014 / Fernandez et al. 1997)
/// @param v Molar volume in J/bar
/// @param temperatureK Temperature in Kelvin
/// @return Dielectric constant
double epsh2o(double v, double temperatureK);

/// Calculate dielectric constant for molecular mixture using Looyenga mixing rule
/// @param vhyb Array of hybrid partial molar volumes (cm3/mol)
/// @param vf Array of volume fractions
/// @param species Species indices (from Perple_X convention)
/// @param nSpecies Number of solvent species
/// @param temperatureK Temperature in Kelvin
/// @return Dielectric constant
double geteps(const std::array<double, 19>& vhyb,
              const std::array<double, 19>& vf,
              const std::vector<int>& species,
              int nSpecies,
              double temperatureK);

/// Calculate Shock et al. (1992) g-function for HKF aqueous model
/// @param rho Solvent density in g/cm3 (CGS units)
/// @param pressureBar Pressure in bar
/// @param temperatureK Temperature in Kelvin
/// @return g-function value in Angstrom
double gfunc(double rho, double pressureBar, double temperatureK);

/// Calculate Debye-Hückel factor for ionic activity coefficients
/// @param msol Molal concentration of solvent (mol/kg)
/// @param vsolv Volume of solvent (cm3/mol)
/// @param epsilon Dielectric constant
/// @param temperatureK Temperature in Kelvin
/// @return ln(activity coefficient) factor for z^2 term
double debyeHuckel(double msol, double vsolv, double epsilon, double temperatureK);

/// Compute hybrid solvent properties (dielectric, DH, g-function) from composition
/// @param yf Array of mole fractions
/// @param vhyb Array of hybrid partial molar volumes (cm3/mol)
/// @param species Species indices
/// @param nSpecies Number of solvent species
/// @param pressureBar Pressure in bar
/// @param temperatureK Temperature in Kelvin
/// @return DielectricState with computed properties
DielectricState computeSolventState(const std::array<double, 19>& yf,
                                    const std::array<double, 19>& vhyb,
                                    const std::vector<int>& species,
                                    int nSpecies,
                                    double pressureBar,
                                    double temperatureK);

/// Compute hybrid solvent Gibbs energy contribution
/// @param yf Array of mole fractions (all species)
/// @param ghybrid Array of hybrid Gibbs fugacity ratios
/// @param nSpecies Number of solvent species
/// @param temperatureK Temperature in Kelvin
/// @return Solvent Gibbs energy contribution (J/mol)
double computeHybridSolventGibbs(const std::array<double, 19>& yf,
                                 const std::array<double, 19>& gmix,
                                 const std::array<double, 19>& g0,
                                 const std::vector<int>& species,
                                 int nSpecies,
                                 double temperatureK);

/// Get default dielectric parameters for Perple_X species
/// @param speciesIndex Perple_X species index (1-18)
/// @return Dielectric parameters
DielectricParams getDefaultDielectricParams(int speciesIndex);

} // namespace Reaktoro::PerpleX
