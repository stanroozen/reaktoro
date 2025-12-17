#include <Reaktoro/Extensions/DEW/WaterDielectricModel.hpp>

#include <cmath>

#include <Reaktoro/Extensions/DEW/WaterDielectricJohnsonNorton.hpp>
#include <Reaktoro/Extensions/DEW/WaterDielectricFranck1990.hpp>
#include <Reaktoro/Extensions/DEW/WaterDielectricFernandez1997.hpp>
#include <Reaktoro/Extensions/DEW/WaterDielectricPowerFunction.hpp>
#include <Reaktoro/Extensions/DEW/WaterPsatPolynomialsDEW.hpp> // module 4
#include <Reaktoro/Extensions/DEW/WaterUtils.hpp>
#include <Reaktoro/Extensions/DEW/WaterUtils.hpp>

namespace Reaktoro {
namespace {

inline real saturationPressure(real T)
{
    // From Wagner–Pruss or your chosen Psat implementation.
    return waterSaturationPressureWagnerPruss(T); // [Pa]
}

inline bool isNearPsat(real T,
                       real P,
                       real reltol)
{
    if (reltol <= 0.0)
        return false;

    const real Psat = saturationPressure(T);
    if (!(Psat > 0.0))
        return false;

    const real diff = autodiff::detail::abs(P - Psat);
    return diff <= reltol * Psat;
}

// Apply DEW Psat(T) epsilon/Q overrides to a base WaterElectroProps.
inline void applyPsatOverrides(real T,
                               const WaterDielectricModelOptions& opts,
                               WaterElectroProps& we)
{
    // epsilon(T) along Psat from Excel polynomial:
    const real eps_psat = waterPsatEpsilonDEW(T);

    we.epsilon = eps_psat;
    we.bornZ   = (eps_psat != 0.0) ? -1.0 / eps_psat : real(0.0);

    if (opts.overrideQWithPsatFit)
    {
        // Q(T) along Psat from Excel polynomial:
        const real Q_psat = waterPsatBornQDEW(T); // 1/bar in Excel → convert inside impl to 1/Pa
        we.bornQ = Q_psat;
    }

    // Along Psat(T), Excel polynomials define only T-dependence.
    // Leave remaining derivatives at 0.0 unless you have fitted forms.
    we.epsilonT  = 0.0;
    we.epsilonP  = 0.0;
    we.epsilonTT = 0.0;
    we.epsilonTP = 0.0;
    we.epsilonPP = 0.0;

    // X, N, U are second derivatives of Z; set to 0 here unless you add
    // consistent Psat fits for them.
    we.bornX = 0.0;
    we.bornN = 0.0;
    we.bornU = 0.0;
}

inline WaterElectroProps evaluatePrimaryModel(real T,
                                              real P,
                                              const WaterThermoProps& wt,
                                              WaterDielectricPrimaryModel model)
{
    switch(model)
    {
        case WaterDielectricPrimaryModel::JohnsonNorton1991:
            return waterElectroPropsJohnsonNorton(T, P, wt);

        case WaterDielectricPrimaryModel::Franck1990:
            return waterElectroPropsFranck1990(T, P, wt);

        case WaterDielectricPrimaryModel::Fernandez1997:
            return waterElectroPropsFernandez1997(T, P, wt);

        case WaterDielectricPrimaryModel::PowerFunction:
            return waterElectroPropsPowerFunction(T, P, wt);
    }

    // Fallback: Johnson & Norton (robust over wide range)
    return waterElectroPropsJohnsonNorton(T, P, wt);
}

} // namespace

//----------------------------------------------------------------------------//
// Public API
//----------------------------------------------------------------------------//

auto waterElectroPropsModel(real T,
                            real P,
                            const WaterThermoProps& wt,
                            const WaterDielectricModelOptions& opts)
    -> WaterElectroProps
{
    // Step 1: evaluate the chosen primary dielectric model
    WaterElectroProps we = evaluatePrimaryModel(T, P, wt, opts.primary);

    // Step 2: optionally override with Psat(T) polynomials
    switch(opts.psatMode)
    {
        case WaterDielectricPsatMode::None:
            // do nothing; base model stands
            break;

        case WaterDielectricPsatMode::UsePsatWhenNear:
            if (isNearPsat(T, P, opts.psatRelativeTolerance))
                applyPsatOverrides(T, opts, we);
            break;

        case WaterDielectricPsatMode::ForcePsat:
            applyPsatOverrides(T, opts, we);
            break;
    }

    return we;
}

} // namespace Reaktoro
