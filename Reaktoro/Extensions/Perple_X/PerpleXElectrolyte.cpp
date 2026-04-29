#include "PerpleXElectrolyte.hpp"

#include <Reaktoro/Common/Exception.hpp>

#include <atomic>
#include <cmath>
#include <stdexcept>

namespace Reaktoro::PerpleX {

namespace {
    constexpr double R = 8.31451;           // J/(mol·K)
    constexpr double CDH = -42182668.74;    // Debye-Hückel constant
    constexpr double EPSILON0 = 78.244;     // Reference dielectric at 25°C, 1 bar (SUPCRT92/DEW, J&N1991 @ 997 kg/m³)
    constexpr double TR = 273.16;           // Reference temperature (K)
    constexpr int MAX_GFUNC_WARNINGS = 10;  // Mirror Perple_X-style warning throttling
    std::atomic<int> gfuncWarningCount{0};

    // Dielectric parameter table for Perple_X species
    // Index: 0=H2O, 1=CO2, 2=CO, 3=CH4, 4=H2, 5=H2S, 6=O2, 7=SO2, 8=COS, 9=N2,
    //        10=NH3, 15=C2H6, 16=HF, 17=HCl
    constexpr DielectricParams SPECIES_PARAMS[19] = {
        // 0 - H2O (handled separately via epsh2o)
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false},
        // 1 - CO2 (H&L 2005)
        {7.3455, 3.35e-3, 0, 83.93, 145.1, -578.8, -1012, 1.55, 0, 0, 0, false},
        // 2 - CO (approximated by O2)
        {3.9578, 6.5e-3, 0, 0.575, 1.028, -8.96, -5.15, 1.5, 0, 0, 0, false},
        // 3 - CH4 (H&L 2005)
        {6.5443, 1.33e-2, 0, 8.4578, 3.7196, -352.97, -100.65, 2, 0, 0, 0, false},
        // 4 - H2 (H&L 2005)
        {2.0306, 5.6e-3, 0, 0.181, 0.021, -7.4, 0, 2, 0, 0, 0, false},
        // 5 - H2S (H&M 2017, polar)
        {1.18, 5829.059676, 9.232464738, -0.01213537391, 0.9, -453374.7482, 3.5, 1.241, -0.241, -16.61833221, 0.5, true},
        // 6 - O2 (H&L 2005)
        {3.9578, 6.5e-3, 0, 0.575, 1.028, -8.96, -5.15, 1.5, 0, 0, 0, false},
        // 7 - SO2 (H&M 2017, polar)
        {2.516, 16242.2847, 10.31715322, -0.00225289526, 0.98, -44.03397284, 1.2, 1.335, 0.335, -16.19171204, 0.75, true},
        // 8 - COS (approximated by CO2)
        {7.3455, 3.35e-3, 0, 83.93, 145.1, -578.8, -1012, 1.55, 0, 0, 0, false},
        // 9 - N2 (H&L 2005)
        {4.3872, 2.26e-3, 0, 2.206, 1.135, -169, -35.83, 2.1, 0, 0, 0, false},
        // 10 - NH3 (not parameterized)
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false},
        // 11-14: Si-O species (not parameterized)
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false},
        // 15 - C2H6 (H&L 2005)
        {11.1552, 0.0112, 0, 36.759, 23.639, -808.03, -378.84, 1.75, 0, 0, 0, false},
        // 16 - HF (approximated by H2S)
        {1.18, 5829.059676, 9.232464738, -0.01213537391, 0.9, -453374.7482, 3.5, 1.241, -0.241, -16.61833221, 0.5, true},
        // 17 - HCl (approximated by H2S)
        {1.18, 5829.059676, 9.232464738, -0.01213537391, 0.9, -453374.7482, 3.5, 1.241, -0.241, -16.61833221, 0.5, true},
        // 18 - placeholder
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false}
    };

    // Molar masses (kg/mol) for Perple_X fluid species indices (1..18).
    // Unused entries are zero.
    constexpr double SPECIES_MW_KG_PER_MOL[19] = {
        0.0,            // 0 unused
        18.01528e-3,    // 1  H2O
        44.0095e-3,     // 2  CO2
        28.0101e-3,     // 3  CO
        16.04246e-3,    // 4  CH4
        2.01588e-3,     // 5  H2
        34.0809e-3,     // 6  H2S
        31.9988e-3,     // 7  O2
        64.0638e-3,     // 8  SO2
        60.0751e-3,     // 9  COS
        28.0134e-3,     // 10 N2
        17.03052e-3,    // 11 NH3
        16.00e-3,       // 12 O
        44.0849e-3,     // 13 SiO
        60.0843e-3,     // 14 SiO2
        28.0855e-3,     // 15 Si
        30.06904e-3,    // 16 C2H6
        20.00634e-3,    // 17 HF
        36.46094e-3     // 18 HCl
    };

    inline double psat2(double temperatureK)
    {
        // Same correlation as Perple_X (rlib.f/PerpleXPureEos.cpp), returns saturation P in bar.
        static constexpr double a[8] = {
            -7.8889166, 2.5514255, -6.716169, 33.239495,
            -105.38479, 174.35319, -148.39348, 48.631602
        };

        if(temperatureK <= 314.0)
            return std::exp(6.3573118 - 8858.843 / temperatureK + 607.56335 / std::pow(temperatureK, 0.6));

        double v = temperatureK / 647.25;
        double w = std::abs(1.0 - v);
        double wsq = std::sqrt(w);
        double ff = 0.0;

        for(int i = 0; i < 8; ++i)
        {
            ff += a[i] * w;
            w *= wsq;
        }

        return 220.93 * std::exp(ff / v);
    }
}

double epsh2o(double v, double temperatureK)
{
    // Sverjenski 2014 / Fernandez et al. 1997 dielectric constant for pure water
    // v is molar volume in J/bar
    // Based on rlib.f lines 3117-3139

    // Convert J/bar to cm3/mol (1 J/bar = 10 cm3)
    double vcm3 = v * 10.0;

    double sqrtt = (temperatureK >= 273.15) ? std::sqrt(temperatureK - 273.15) : 0.0;

    double eps = std::exp(-8.016651e-5 * temperatureK + 4.769870482 - 0.06871618 * sqrtt)
               * std::pow(18.01526833 / vcm3, -1.576377e-3 * temperatureK + 1.185462878 + 0.06810288 * sqrtt);

    return eps;
}

double geteps(const std::array<double, 19>& vhyb,
              const std::array<double, 19>& vf,
              const std::vector<int>& species,
              int nSpecies,
              double temperatureK)
{
    // Based on rlib.f lines 11851-12000 (geteps subroutine)

    double trt = temperatureK / TR - 1.0;
    double epsln = 0.0;

    // Process all non-water species.
    for (int i = 0; i < nSpecies; ++i)
    {
        const int j = species[i];
        if (j <= 0 || j >= 19 || j == 1)
            continue;

        const auto& po = SPECIES_PARAMS[j - 1];

        // Density in 1/cm3
        double rho = 1.0 / vhyb[j];

        double eps;

        if (!po.isPolar)
        {
            // Non-polar: Eq 5 of Harvey & Lemmon 2005
            // Polarization/rho
            eps = po.a0 + po.a1 * trt + (po.b0 + po.b1 * trt) * rho
                + (po.c0 + po.c1 * trt) * std::pow(rho, po.D);

            // Invert Clausius-Mossotti relation for dielectric constant
            eps = (2.0 * eps * rho + 1.0) / (1.0 - rho * eps);
        }
        else
        {
            // Polar: Harvey & Mountain 2017
            // Struct field ↔ Fortran po(j,k) mapping:
            //   a0=po(1), a1=po(2)=A_mu coeff, A_mu=po(3)=C_alpha_outer,
            //   b0=po(4), b1=po(5), c0=po(6), c1=po(7),
            //   D=po(8), C_alpha=po(9)=Kirkwood g mult, g0=po(10)=exp coeff, g1=po(11)=rho exp
            //
            // Fortran:
            //   eps = rho*(po(3) + po(2)*(po(1)*exp(po(4)*T^po(5))*(1-exp(po(6)*rho^po(7)))+1)
            //             *(po(8)+po(9)*exp(po(10)*rho^po(11)))^2/T)
            eps = rho * (po.A_mu + po.a1
                * (po.a0 * std::exp(po.b0 * std::pow(temperatureK, po.b1))
                * (1.0 - std::exp(po.c0 * std::pow(rho, po.c1))) + 1.0)
                * std::pow(po.D + po.C_alpha * std::exp(po.g0 * std::pow(rho, po.g1)), 2) / temperatureK);

            // Invert Kirkwood relation
            eps = 2.25 * eps + 0.25
                + std::sqrt((5.0625 * eps + 1.125) * eps + 0.5625);
        }

        // Looyenga mixing rule with volume fractions
        epsln += vf[j] * std::pow(eps, 1.0/3.0);
    }

    // Add water, if present, using epsh2o.
    bool hasWater = false;
    for (int i = 0; i < nSpecies; ++i)
    {
        if (species[i] == 1)
        {
            hasWater = true;
            break;
        }
    }
    if (hasWater)
    {
        // Convert vhyb from cm3/mol to J/bar
        double v_jbar = vhyb[1] / 10.0;
        epsln += vf[1] * std::pow(epsh2o(v_jbar, temperatureK), 1.0/3.0);
    }

    return std::pow(epsln, 3.0);
}
double gfunc(double rho, double pressureBar, double temperatureK)
{
    // Shock et al. (1992) g-function for HKF aqueous model
    // Based on rlib.f lines 3039-3109
    // rho is CGS solvent density (g/cm3)

    double g = 0.0;

    if (rho > 1.0)
    {
        // Region III: rho = 1 g/cm3, g = 0
        return 0.0;
    }

    // Region I function
    g = ((-6.557892e-6 * temperatureK + 9.3295764e-3) * temperatureK - 4.096745422)
      * std::pow(1.0 - rho, (1.268348e-5 * temperatureK - 1.767275512e-2) * temperatureK + 9.98834792);

    if (temperatureK > 428.15 && pressureBar < 1000.0)
    {
        // Add region II perturbation term
        double tf = (temperatureK / 300.0 - 1.427166667);

        g -= (std::pow(tf, 4.8) + 0.366666e-15 * std::pow(tf, 16))
           * ((((5.01799e-14 * pressureBar - 5.0224e-11) * pressureBar - 1.504074e-7) * pressureBar
               + 2.507672e-4) * pressureBar - 0.1003157);
    }

    // Physical bounds check (Perple_X rlib.f behavior): warn and zero g outside valid region.
    const bool outOfBounds =
        rho < 0.35
        || (temperatureK > 623.15 && pressureBar < 500.0)
        || (temperatureK <= 623.15 && pressureBar < psat2(temperatureK));
    if(outOfBounds)
    {
        const int seen = gfuncWarningCount.load();
        if(seen < MAX_GFUNC_WARNINGS)
        {
            warningif(true,
                "PerpleXElectrolyte::gfunc: T=", temperatureK, " K, P=", pressureBar,
                " bar is beyond HKF g-function validity limits; returning g=0 to match Perple_X.");
            const int now = gfuncWarningCount.fetch_add(1) + 1;
            if(now == MAX_GFUNC_WARNINGS)
            {
                warningif(true,
                    "PerpleXElectrolyte::gfunc: reached warning limit (", MAX_GFUNC_WARNINGS,
                    "); further out-of-bounds warnings are suppressed.");
            }
        }
        return 0.0;
    }

    return g;
}

double debyeHuckel(double msol, double vsolv, double epsilon, double temperatureK)
{
    // Debye-Hückel factor calculation
    // Based on rlib.f lines 11724-11725
    // adh = cdh * sqrt(10*msol/vsolv/(epsln*t)^3)

    if (epsilon <= 0.0 || temperatureK <= 0.0)
    {
        throw std::runtime_error("Invalid parameters for Debye-Hückel calculation");
    }

    double denominator = std::pow(epsilon * temperatureK, 3);
    double adh = CDH * std::sqrt(10.0 * msol / vsolv / denominator);

    return adh;
}

DielectricParams getDefaultDielectricParams(int speciesIndex)
{
    if (speciesIndex < 1 || speciesIndex > 18)
    {
        return DielectricParams{};
    }

    return SPECIES_PARAMS[speciesIndex - 1];
}

DielectricState computeSolventState(const std::array<double, 19>& yf,
                                    const std::array<double, 19>& vhyb,
                                    const std::vector<int>& species,
                                    int nSpecies,
                                    double pressureBar,
                                    double temperatureK)
{
    // Based on rlib.f lines 11640-11756 (slvnt1 subroutine)

    DielectricState state;

    // Compute total hybrid volume and volume fractions
    state.hyvol = 0.0;
    double ysum = 0.0;
    for (int i = 0; i < nSpecies; ++i)
    {
        const int j = species[i];
        if (j > 0 && j < 19)
        {
            ysum += yf[j];
            state.hyvol += yf[j] * vhyb[j];
        }
    }

    // Compute volume fractions
    std::array<double, 19> vf{};
    if (state.hyvol > 0.0)
    {
        for (int i = 0; i < nSpecies; ++i)
        {
            const int j = species[i];
            if (j > 0 && j < 19)
            {
                vf[j] = yf[j] * vhyb[j] / state.hyvol;
            }
        }
    }

    // Compute dielectric constant using volume fractions
    state.epsilon = geteps(vhyb, vf, species, nSpecies, temperatureK);

    // Compute solvent properties for DH and g-function.
    // Match Perple_X behavior: solvent mass is composition-weighted over solvent species.
    state.msol = 0.0;
    for (int i = 0; i < nSpecies; ++i)
    {
        const int j = species[i];
        if (j > 0 && j < 19)
            state.msol += yf[j] * SPECIES_MW_KG_PER_MOL[j];
    }

    // Store H₂O mole fraction among solvents (yf[1]) for msol_scale correction.
    // Used in ActivityModelPerplexDEW to compute the correct Perple_X total-solvent
    // molality scale:  msol_scale = M_H2O × x_H2O / msol  (not just M_H2O / msol).
    state.x_H2O = yf[1];

    // Perple_X uses vsolv = ysum * hyvol.
    state.vsolv = ysum * state.hyvol;

    // Compute solvent density for g-function (CGS: g/cm3)
    double rho_solvent = (state.vsolv > 0.0) ? (state.msol * 1000.0 / state.vsolv) : 0.0;

    // Compute g-function
    state.gf = gfunc(rho_solvent, pressureBar, temperatureK);

    // Compute Debye-Hückel factor
    if (state.vsolv > 0.0)
        state.adh = debyeHuckel(state.msol, state.vsolv, state.epsilon, temperatureK);

    return state;
}

double computeHybridSolventGibbs(const std::array<double, 19>& yf,
                                 const std::array<double, 19>& gmix,
                                 const std::array<double, 19>& g0,
                                 const std::vector<int>& species,
                                 int nSpecies,
                                 double temperatureK)
{
    // Based on Perple_X slvnt1 + ghybrid:
    // ghybrid(x) = R*T*sum_i x_i*log(x_i*gmix_i/g0_i)
    // gsolv correction = ysum*(ghybrid(ysolv) + R*T*log(ysum))

    constexpr double R = 8.31451;  // J/(mol·K)

    // Compute normalized solvent fractions
    double ysum = 0.0;
    for (int i = 0; i < nSpecies; ++i)
    {
        const int j = species[i];
        if (j > 0 && j < 19)
            ysum += yf[j];
    }

    if (ysum <= 0.0)
    {
        return 0.0;
    }

    double gtemp = 0.0;
    for (int i = 0; i < nSpecies; ++i)
    {
        const int j = species[i];
        if (j <= 0 || j >= 19)
            continue;

        const double x = yf[j] / ysum;
        if (x <= 0.0)
            continue;

        // Guard against undefined logs in degenerate states.
        if (gmix[j] <= 0.0 || g0[j] <= 0.0)
            continue;

        gtemp += x * std::log(x * gmix[j] / g0[j]);
    }

    const double ghyb = R * temperatureK * gtemp;
    return ysum * (ghyb + R * temperatureK * std::log(ysum));
}

} // namespace Reaktoro::PerpleX

