#pragma once

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

/// Activity coefficient model for ActivityModelPerplexDEW.
enum class ActivityDHModel
{
    /// Davies approximation matching Perple_X GFSM (slvnt2 / aqact) exactly:
    ///   ln(γᵢ) = zᵢ² × (adh × √I/(1+√I) + 0.2×I)
    /// where adh = CDH/√(V_H₂O/10 × (εT)³)  (CDH = −5661800.47810 for pure H₂O,
    /// Looyenga-mixed ε for multi-component solvents).
    Davies,

    /// Extended Debye-Hückel with ionic radii (corrected SI formulas):
    ///   log₁₀(γᵢ) = −A zᵢ² √I / (1 + aᵢ B √I)
    /// A and B computed from the mixed-solvent ε and ρ
    /// (ZD05 + Looyenga/Fernandez-Sverjensky, same as Davies path).
    ExtendedDH
};

struct ActivityModelParamsPerplexDEW
{
    ActivityDHModel model = ActivityDHModel::Davies;
    // Hydrous-equivalent correction is optional and disabled by default.
    bool enableHydrousSpeciesCorrection = false;
};

/// Return an aqueous activity model using Perple_X solvent internals.
///
/// Solvent state selection follows Perple_X GFSM (COH-Fluid+, model type 39):
///  - Pure H₂O  → slvnt0 path: ZD05 EOS + Fernandez/Sverjensky dielectric.
///  - Mixed fluid (H₂O + CO₂, CH₄, H₂S, SO₂, H₂, CO, N₂, NH₃, HF, C₂H₆, HCl)
///               → slvnt1 path: MRK pure volumes + Looyenga mixing rule for ε.
///
/// @param params  PerplexDEW activity-model options.
auto ActivityModelPerplexDEW(const ActivityModelParamsPerplexDEW& params)
    -> ActivityModelGenerator;

/// Convenience overload using only the Debye-Hückel model selector.
auto ActivityModelPerplexDEW(ActivityDHModel model = ActivityDHModel::Davies)
    -> ActivityModelGenerator;

} // namespace Reaktoro
