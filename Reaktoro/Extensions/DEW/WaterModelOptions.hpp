#pragma once

#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

/// Equation of state model for pure water.
enum class WaterEosModel {
    WagnerPruss,     ///< IAPWS-95 style (via Helmholtz)
    HGK,             ///< Haar-Gallagher-Kell (if enabled)
    ZhangDuan2005,   ///< Zhang & Duan (2005), DEW-style
    ZhangDuan2009,   ///< Zhang & Duan (2009), DEW-style
};

/// Dielectric constant model for water.
enum class WaterDielectricModel {
    JohnsonNorton1991,  ///< Johnson & Norton (1991) (SUPCRT / DEW default)
    Franck1990,         ///< Franck et al. (1990)
    Fernandez1997,      ///< Fernandez et al. (1997)
    PowerFunction,      ///< Sverjensky-Harrison power law
};

/// Gibbs free energy model for water (required for thermodynamic calculations).
enum class WaterGibbsModel {
    DelaneyHelgeson1978,  ///< Polynomial fit (Delaney & Helgeson, 1978).
    DewIntegral,          ///< DEW-style integral of V(P,T) from 1 kb.
};

/// Born / solvation model for aqueous species (optional DEW compatibility).
enum class WaterBornModel {
    None,        ///< No DEW Shock-style omega model.
    Shock92Dew,  ///< Shock et al. (1992) / DEW omega(g) + domega/dP.
};

/// Global water model configuration for Reaktoro-DEW style usage.
struct WaterModelOptions
{
    /// EOS controlling water density & its derivatives.
    WaterEosModel eosModel = WaterEosModel::ZhangDuan2005;

    /// Dielectric model controlling epsilon(T, P, rho).
    WaterDielectricModel dielectricModel = WaterDielectricModel::PowerFunction;

    /// Gibbs free energy model (required for species thermodynamics).
    WaterGibbsModel gibbsModel = WaterGibbsModel::DewIntegral;

    /// Optional Born omega model for solvation (if needed).
    WaterBornModel bornModel = WaterBornModel::None;

    /// If true, use DEW-style Psat polynomials in a small neighborhood of Psat(T).
    bool usePsatPolynomials = false;

    /// Relative tolerance |P - Psat| / Psat below which Psat polynomials are applied.
    /// Only used when usePsatPolynomials == true.
    real psatRelTol = 1e-3;

    /// Density calculation tolerance [bar] for Zhang & Duan EOS.
    /// Default 0.001 bar gives high accuracy.
    double densityTolerance = 0.001;
};

/// Convenience: Construct a WaterModelOptions corresponding to
/// "classic DEW-style" behavior.
auto makeWaterModelOptionsDEW() -> WaterModelOptions;

} // namespace Reaktoro
