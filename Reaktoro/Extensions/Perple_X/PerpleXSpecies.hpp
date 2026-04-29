#pragma once

#include <string>
#include <string_view>

namespace Reaktoro::PerpleX {

/// Species indices as used by Perple_X MRK routines.
enum class Species : int {
    H2O  = 1,
    CO2  = 2,
    CO   = 3,
    CH4  = 4,
    H2   = 5,
    H2S  = 6,
    O2   = 7,
    SO2  = 8,
    N2   = 10,
    NH3  = 11,
    C2H6 = 16,
    HF   = 17,
    HCl  = 18,
};

constexpr int speciesCount() noexcept { return 13; }

constexpr int toIndex(Species s) noexcept { return static_cast<int>(s); }

inline std::string_view name(Species s)
{
    switch(s)
    {
    case Species::H2O:  return "H2O";
    case Species::CO2:  return "CO2";
    case Species::CO:   return "CO";
    case Species::CH4:  return "CH4";
    case Species::H2:   return "H2";
    case Species::H2S:  return "H2S";
    case Species::O2:   return "O2";
    case Species::SO2:  return "SO2";
    case Species::N2:   return "N2";
    case Species::NH3:  return "NH3";
    case Species::C2H6: return "C2H6";
    case Species::HF:   return "HF";
    case Species::HCl:  return "HCl";
    }

    return "";
}

} // namespace Reaktoro::PerpleX
