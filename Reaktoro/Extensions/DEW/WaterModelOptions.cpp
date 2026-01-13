#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>

namespace Reaktoro {

auto makeWaterModelOptionsDEW() -> WaterModelOptions
{
    WaterModelOptions opt;

    // DEW canonical choices:
    opt.eosModel        = WaterEosModel::ZhangDuan2005;
    opt.dielectricModel = WaterDielectricModel::PowerFunction;

    // Use DEW-style ∫V dP integral for G(H2O)
    opt.gibbsModel      = WaterGibbsModel::DewIntegral;

    // Use Shock et al. (1992) / DEW omega(g) model for ions.
    opt.bornModel       = WaterBornModel::Shock92Dew;

    // Enable DEW Psat polynomials near saturation.
    opt.usePsatPolynomials = true;
    opt.psatRelTol         = 1e-3; // can tune if needed

    return opt;
}

} // namespace Reaktoro
