#include <Reaktoro/Extensions/DEW/WaterElectroModel.hpp>

#include <cmath>

#include <Reaktoro/Extensions/DEW/WaterDielectricJohnsonNorton.hpp>
#include <Reaktoro/Extensions/DEW/WaterDielectricFranck1990.hpp>
#include <Reaktoro/Extensions/DEW/WaterDielectricFernandez1997.hpp>
#include <Reaktoro/Extensions/DEW/WaterDielectricPowerFunction.hpp>

#include <Reaktoro/Extensions/DEW/WaterPsatPolynomialsDEW.hpp>
#include <Reaktoro/Extensions/DEW/WaterHelmholtzPropsWagnerPruss.hpp>
#include <Reaktoro/Extensions/DEW/WaterUtils.hpp> // for saturation helpers if present

namespace Reaktoro {
namespace {

using std::abs;

// Compute Psat(T) via Wagner–Pruss (IAPWS-95) helper.
// Adjust the function name if your existing API differs.
inline double saturationPressure(double T_K)
{
    return waterSaturationPressureWagnerPruss(T_K); // [Pa]
}

inline bool isNearPsat(double T_K,
                       double P_Pa,
                       double relTol)
{
    if (!std::isfinite(T_K) || !std::isfinite(P_Pa) || relTol <= 0.0)
        return false;

    const double Psat = saturationPressure(T_K);
    if (!(Psat > 0.0))
        return false;

    const double diff = std::abs(P_Pa - Psat);
    return diff <= relTol * Psat;
}

/// Apply DEW Psat polynomials to override epsilon & Q near saturation.
///
/// Matches the Excel `Psat = True` semantics:
///   - Epsilon along Psat(T) given by Psat epsilon polynomial.
///   - Q along Psat(T) given by Psat Q polynomial.
///   - Then epsilonP is reconstructed from Q via:
///         Q = (1/epsilon^2) * (d epsilon / dP)
///       => d epsilon / dP = Q * epsilon^2
///
/// We do *not* fabricate higher-order derivatives. Unset fields are left
/// as in the underlying model output to avoid inventing behavior.
inline void maybeApplyPsatOverride(real T,
                                   real P,
                                   const WaterElectroModelOptions& opt,
                                   WaterElectroProps& we)
{
    if (!opt.usePsatPolynomials)
        return;

    if (!isNearPsat(T, P, opt.psatRelativeTolerance))
        return;

    // Epsilon along Psat(T)
    const double eps_p = waterPsatEpsilonDEW(T);

    // Born Q along Psat(T) [1/Pa]
    const double Q_p = waterPsatBornQDEW(T);

    if (eps_p <= 0.0)
        return; // avoid nonsense; keep underlying model values

    // Override epsilon
    we.epsilon = eps_p;

    // Born Z = -1/epsilon
    we.bornZ = -1.0 / eps_p;

    // Override bornQ with Psat polynomial:
    we.bornQ = Q_p;

    // From definition Q = (1/ε²) * ε_P  =>  ε_P = Q * ε²
    const double eps2 = eps_p * eps_p;
    we.epsilonP = Q_p * eps2;

    // We deliberately do NOT force epsilonT, epsilonTT, etc. here.
    // They remain whatever the underlying model produced (often 0.0),
    // which matches DEW: Psat polynomials are special-case fits,
    // not a full thermodynamic surface.
}

} // namespace

//------------------------------------------------------------------------------
// Public selector
//------------------------------------------------------------------------------

auto waterElectroPropsModel(real T,
                            real P,
                            const WaterThermoProps& wt,
                            const WaterElectroModelOptions& opt)
    -> WaterElectroProps
{
    WaterElectroProps we;

    // 1) Base model selection: exactly one of the modular implementations.
    switch (opt.model)
    {
        case WaterDielectricModel::JohnsonNorton1991:
        {
            we = waterElectroPropsJohnsonNorton(T, P, wt);
            break;
        }

        case WaterDielectricModel::Franck1990:
        {
            we = waterElectroPropsFranck1990(T, P, wt);
            break;
        }

        case WaterDielectricModel::Fernandez1997:
        {
            we = waterElectroPropsFernandez1997(T, P, wt);
            break;
        }

        case WaterDielectricModel::PowerFunction:
        {
            we = waterElectroPropsPowerFunction(T, P, wt);
            break;
        }
    }

    // 2) Optional DEW Psat polynomial override (Psat=True semantics).
    maybeApplyPsatOverride(T, P, opt, we);

    return we;
}

} // namespace Reaktoro
