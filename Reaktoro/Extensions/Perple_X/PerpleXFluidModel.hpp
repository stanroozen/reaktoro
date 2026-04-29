#pragma once

#include <array>
#include <vector>

#include "PerpleXHybridEos.hpp"
#include "PerpleXMrkMixture.hpp"
#include "PerpleXMrkPure.hpp"
#include "PerpleXElectrolyte.hpp"

namespace Reaktoro::PerpleX {

/// ============================================================================
/// DEPRECATED: Use PerpleXGFSMModel instead for GFSM (ifug=39) calculations
/// ============================================================================
///
/// This class provides a generic interface to the older fluid model system.
/// For new code, prefer PerpleXGFSMModel which explicitly implements the
/// Generic Fluid Solution Model (GFSM, Perple_X ifug=39).
///
/// NOTE ON MODEL DISTINCTION:
/// - PerpleXFluidModel: Legacy generic interface (may support multiple models)
/// - PerpleXGFSMModel: Explicit GFSM implementation (ifug=39, speciation space)
/// - Binary models (ifug=0-5): Fixed binary composition space (NOT in this class)
///

struct PerpleXFluidOptions
{
    bool useLowTMrk = false;
    bool enableElectrolyte = false;  // Enable solution model coupling
    MrkMixOptions mixOptions{};
    HybridEosOptions hybridOptions{};
    std::vector<int> hybridSpecies;
};

struct PerpleXFluidState
{
    std::array<double, 19> ln_f{};
    std::array<double, 19> g{};
    std::array<double, 19> v{};
    std::array<double, 19> gh{};
    std::array<double, 19> vh{};
    std::array<double, 19> vhyb{};
    std::array<double, 19> vf{};
    double vol = 0.0;
    double hyvol = 0.0;
    DielectricState dielectric{};
    double gsolv = 0.0;  // Hybrid solvent Gibbs contribution
};

class PerpleXFluidModel
{
public:
    PerpleXFluidState compute(const std::vector<int>& species,
                              const std::array<double, 19>& y,
                              double pressureBar,
                              double temperatureK,
                              const PerpleXFluidOptions& options,
                              MrkRootState* rootState = nullptr) const;
};

} // namespace Reaktoro::PerpleX
