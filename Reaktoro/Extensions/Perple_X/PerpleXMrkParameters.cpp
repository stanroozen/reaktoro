#include "PerpleXMrkParameters.hpp"

#include <cmath>

namespace Reaktoro::PerpleX {
namespace {

constexpr std::array<double, 19> ark = {
    0.0,
    0.0,                // 1 H2O (overridden)
    0.0,                // 2 CO2 (overridden)
    16.98e6,            // 3 CO
    32.154811e6,        // 4 CH4
    2.821e5,            // 5 H2 (overridden)
    89.0e6,             // 6 H2S
    174026.0e2,         // 7 O2
    133.1e6,            // 8 SO2
    130.0e6,            // 9 COS
    136.0e5,            // 10 N2
    631.0e5,            // 11 NH3
    174026.0e2,         // 12 O
    424441664.6,        // 13 SiO (placeholder; overridden below)
    7373939618.0,       // 14 SiO2 (overridden)
    3767833334.0,       // 15 Si (overridden)
    98774720.4,         // 16 C2H6
    7284049.7,          // 17 HF
    7284049.7           // 18 HCl (HF placeholder)
};

constexpr std::array<double, 19> brk = {
    0.0,
    14.6,               // 1 H2O
    29.7,               // 2 CO2
    27.38,              // 3 CO
    29.681,             // 4 CH4
    15.7699,            // 5 H2 (overridden)
    29.94,              // 6 H2S
    22.07,              // 7 O2
    37.4,               // 8 SO2
    43.0,               // 9 COS
    23.42,              // 10 N2
    18.84,              // 11 NH3
    22.07,              // 12 O
    23.81,              // 13 SiO
    25.83798814,        // 14 SiO2
    10.35788774,        // 15 Si
    45.139,             // 16 C2H6
    17.93096733,        // 17 HF
    17.93096733         // 18 HCl (HF placeholder)
};

} // namespace

MrkParameters mrkParameters(double T)
{
    MrkParameters params{};

    for(int i = 1; i <= 18; ++i)
    {
        params.b[i] = brk[i];

        if(i == 1)
        {
            if(T > 300.0)
            {
                params.a[i] = 0.1452535403e8
                    + T * (306893.3587 + T * (-307.9995871
                    + T * (0.09226256008 - 0.2930106337e-5 * T)));
            }
            else
            {
                params.b[i] = 16.0;
                params.a[i] = 127354240.0;
            }
        }
        else if(i == 2)
        {
            params.a[i] = 92935540.0 + T * (-82130.73 + 21.29 * T);
        }
        else if(i == 5)
        {
            params.b[i] = 12.81508162;
            params.a[i] = 0.391950132949994654e8
                + T * (-0.881231157499978144e5)
                + (T * T) * 0.890185987380923081e2
                + (T * T * T) * (-0.286881183333320412e-1);
        }
        else if(i == 14)
        {
            const double dT = T - 1999.0;
            params.a[i] = (
                -0.370631646315402e9 - 88784.52
                + 0.710713269453173e8 * std::log(T)
                - 0.468778070702675e7 / T
                + (0.194790021605110e4 * std::sqrt(T)
                   - 0.110935131465938e6
                   - 0.120230245951606e2 * T) * T
            ) * 1.0e2
            + 32300.0 * dT + 14.25 * dT * dT;
        }
        else if(i == 15)
        {
            const double dT = T - 1687.0;
            params.a[i] = (
                0.131596431388077e7
                - ((0.380259023635694e-1 * T
                    + 0.124090483523393e4) * T
                    + 0.170392520137105e7) * std::sqrt(T)
                + 0.151371320806448e6 / std::sqrt(T)
                + 0.427563259532326e7 * std::log(T)
                + (0.108181901455347e2 * T + 0.711400073165747e5) * T
                + 17737.22
                - 50.5 * dT
                - 2.04e-2 * dT * dT
            ) * 1.0e2;
        }
        else
        {
            params.a[i] = ark[i];
        }

        if(params.a[i] < 0.0)
            params.a[i] = 1.0;
    }

    // SiO is derived from SiO2 in Perple_X
    params.a[13] = ark[14] / 20.0;
    params.b[13] = brk[13];

    return params;
}

} // namespace Reaktoro::PerpleX
