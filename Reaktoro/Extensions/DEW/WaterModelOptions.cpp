#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>

namespace Reaktoro {

auto makeWaterModelOptionsDEW() -> WaterModelOptions
{
    WaterModelOptions opt;

    // DEW canonical choices (based on the Excel/VBA implementation):
    opt.eosModel        = WaterEosModel::ZhangDuan2005;
    opt.dielectricModel = WaterDielectricModel::JohnsonNorton1991;

    // Use Delaney & Helgeson (1978) G(H2O), as in DEW option 1.
    opt.gibbsModel      = WaterGibbsModel::DelaneyHelgeson1978;

    // Use Shock et al. (1992) / DEW omega(g) model for ions.
    opt.bornModel       = WaterBornModel::Shock92Dew;

    // Enable DEW Psat polynomials near saturation.
    opt.usePsatPolynomials = true;
    opt.psatRelTol         = 1e-3; // can tune if needed

    return opt;
}

} // namespace Reaktoro
