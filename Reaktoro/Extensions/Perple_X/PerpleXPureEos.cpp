#include "PerpleXPureEos.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "PerpleXMrkPure.hpp"

namespace Reaktoro::PerpleX {

/// ============================================================================
/// PURE EOS IMPLEMENTATIONS FOR GFSM (EXPLICIT SPECIATION SPACE MODEL)
/// ============================================================================
///
/// GFSM SPECIATION SPACE (Explicit, NOT composition space):
/// ==========================================================
///
/// GFSM computes properties in EXPLICIT SPECIATION SPACE:
/// - Input: All 12 mole fractions (Xn_CO2, Xn_H2O, Xn_CH4, ... Xn_HCl)
/// - Process: Evaluate pure EOS for each species independently
/// - Output: Thermodynamic properties (fugacity, volume, Gibbs energy)
///
/// Pure EOS Options in GFSM:
///
/// For H2O (7 total alternatives):
/// - MRK: Modified Redlich-Kwong (always available)
/// - HSMRK: Hybrid with hydrogen bonding effects
/// - CORK: Coupled Oscillator RK formulation
/// - PSEOS: Perturbation scaled EOS
/// - Haar: Wagner-Haar pure water correlation
/// - ZhangDuan05: Zhang-Duan 2005 formulation
/// - ZhangDuan09: Zhang-Duan 2009 formulation
///
/// For CO2 (6 total alternatives):
/// - MRK, HSMRK, CORK, BRMRK, PSEOS, ZhangDuan09
///
/// For CH4 (3 total alternatives):
/// - MRK, HSMRK, ZhangDuan09
///
/// For other 9 species (H2S, SO2, H2, CO, N2, NH3, HF, C2H6, HCl):
/// - FIXED to MRK (no alternatives available)
///
/// How GFSM uses these:
/// 1. Start with MRK foundation (all 12 species via MRK mixing rules)
/// 2. REPLACE H2O, CO2, CH4 pure EOS with user-selected alternatives
/// 3. Other 9 species remain on MRK
/// 4. Result: EXPLICIT hybrid model in speciation space
namespace {

constexpr double rkR = 83.1441; // cm3·bar/(mol·K)
constexpr double zdR = 8.31441; // Unscaled gas constant used in Perple_X ZD routines
constexpr double haarVolumeCompat = 1.2053110886584423;
constexpr double zd09Co2VolumeCompat = 0.7636447435629266;
constexpr double zd09Ch4VolumeCompat = 0.7450000000000000;

inline double fugp(double rtt, double b, double yz, double c, double d, double e, double v)
{
    const double y = b / 4.0 / v;
    const double vpb = v + b;
    const double dvbv = std::log(vpb / v);
    const double y1 = 1.0 - y;

    return ((4.0 - 3.0 * y) * y + (2.0 - y) * 2.0 * y / y1) / (y1 * y1)
        + (d * (dvbv / b + (4.0 * y + 2.0) / vpb - 3.0 / v)
           - c * (dvbv + b / vpb)
           + e * ((4.0 / b - 2.0 / v) / v - dvbv / b / b
                  + ((2.0 * y - 1.5) / v - 3.0 / b) / vpb)) / rtt / b
        - std::log(yz);
}

inline bool nurap(double b, double c, double d, double e,
                  double& yz, double& vi,
                  double temperatureK, double pressureBar,
                  const PerpleXPureEosOptions& options,
                  double r = rkR)
{
    const double t12 = std::sqrt(temperatureK);
    const double s1 = r * temperatureK * t12;
    const double s2 = b * s1;
    const double s3 = t12 * pressureBar * b;

    const double p0 = -256.0 * s1;
    const double b2 = b * b;

    const double q0 = 256.0 * t12 * pressureBar;
    const double q1 = 256.0 * (s3 - s1);
    const double q2 = (-160.0 * s3 - 512.0 * s1) * b + 256.0 * c;
    const double q3 = (-80.0 * s3 + p0) * b2 + 256.0 * d;
    const double q4 = ((65.0 * s3 + 8.0 * s1) * b - 160.0 * c) * b2 + 256.0 * e;
    const double q5 = -b2 * (((14.0 * s3 - 15.0 * s1) * b - 80.0 * c) * b + 160.0 * d);
    const double q6 = b2 * ((((s3 + 6.0 * s1) * b - 15.0 * c) * b + 80.0 * d) * b - 160.0 * e);
    const double q7 = b * b * b * (((-s1 * b + c) * b - 15.0 * d) * b + 80.0 * e);
    const double q8 = b2 * b2 * (-15.0 * e + d * b);
    const double q9 = e * b2 * b2 * b;

    const double p1 = 512.0 * c - 768.0 * s2;
    const double p2 = (-832.0 * s2 - 256.0 * c) * b + 768.0 * d;
    const double p3 = ((-368.0 * s2 - 64.0 * c) * b - 256.0 * d) * b + 1024.0 * e;
    const double p4 = -b * (((33.0 * s2 - 64.0 * c) * b + 224.0 * d) * b + 256.0 * e);
    const double p5 = 2.0 * b2 * (b * ((s2 - c) * 7.0 * b + 72.0 * d) - 192.0 * e);
    const double p6 = -b * b * b * (b * ((s2 - c) * b + 29.0 * d) - 224.0 * e);
    const double p7 = 2.0 * b2 * b2 * (d * b - 22.0 * e);
    const double p8 = 3.0 * q9;

    for(int k = 0; k <= options.maxIter; ++k)
    {
        const double num = ((((((((((q0 * vi + q1) * vi + q2) * vi + q3) * vi + q4) * vi + q5) * vi
                           + q6) * vi + q7) * vi + q8) * vi + q9) * vi);
        const double den = (((((((p0 * vi + p1) * vi + p2) * vi + p3) * vi + p4) * vi + p5) * vi + p6) * vi + p7) * vi + p8;
        const double cor = num / den;
        vi = vi + cor;

        if(std::abs(cor / vi) < options.tol)
        {
            yz = vi * pressureBar / r / temperatureK;
            return true;
        }

        if(vi < 0.0)
            return false;
    }

    return false;
}

inline double psat2(double temperatureK)
{
    static constexpr double a[8] = {
        -7.8889166, 2.5514255, -6.716169, 33.239495,
        -105.38479, 174.35319, -148.39348, 48.631602
    };

    if(temperatureK <= 314.0)
    {
        return std::exp(6.3573118 - 8858.843 / temperatureK + 607.56335 / std::pow(temperatureK, 0.6));
    }

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

inline double aideal(double tr, double rt)
{
    static constexpr double ci[18] = {
        1.9730271018e1, 2.09662681977e1, -4.83429455355e-1,
        6.05743189245, 2.256023885e1, -9.87532442,
        -4.3135538513, 4.58155781e-1, -4.7754901883e-2,
        4.1238460633e-3, -2.7929052852e-4, 1.4481695261e-5,
        -5.6473658748e-7, 1.6200446e-8, -3.303822796e-10,
        4.51916067368e-12, -3.70734122708e-14, 1.37546068238e-16
    };

    double w = std::pow(tr, -3.0);
    double aid = 1.0 + (ci[0] / tr + ci[1]) * std::log(tr);

    for(int i = 2; i < 18; ++i)
    {
        aid += ci[i] * w;
        w *= tr;
    }

    return -rt * aid;
}

inline double br_a(double v, double a1, double a2, double a3)
{
    const double a3v = a3 / v;
    const double a3v3 = a3v * a3v * a3v;
    return a1 * (a3v3 - a3v3 * a3v3) + a2;
}

inline double br_b(double v, double a3, double b1, double b2)
{
    return (std::log(v / a3) + b1) / b2;
}

inline void brvol(double pressureBar, double temperatureK, double& vol,
                  const PerpleXPureEosOptions& options)
{
    static constexpr double rbar = 83.143;
    static constexpr double a1 = 6.566e7;
    static constexpr double a2 = 7.276e7;
    static constexpr double a3 = 37.3;

    const double rt = rbar * temperatureK;
    const double t12 = std::sqrt(temperatureK);

    double v = vol;
    double dv = 5.0e-5;

    for(int it = 0; it <= 50; ++it)
    {
        double b1 = 7.352629;
        double b2 = 0.241413;

        if(v <= 47.22)
        {
            b1 = 1.856669;
            b2 = 0.0637935;
        }
        else if(v < 180.0)
        {
            b1 = 11.707864;
            b2 = 0.363955;
        }

        const double b = br_b(v, a3, b1, b2);
        const double a = br_a(v, a1, a2, a3);

        const double vp = v + dv;
        const double bp = br_b(vp, a3, b1, b2);
        const double ap = br_a(vp, a1, a2, a3);

        const double fv = rt / (v - b) - a / (v * (v + b) * t12) - pressureBar;
        const double dpdv = (fv - (rt / (vp - bp) - ap / (vp * (vp + bp) * t12) - pressureBar)) / dv;

        const double corr = fv / dpdv;
        v = v + corr;

        if(std::abs(corr) < 1.0e-3)
        {
            vol = v;
            return;
        }
    }

    vol = v;
}

inline double vdpdv(double v, double temperatureK)
{
    static constexpr double rbar = 83.143;
    static constexpr double a1 = 6.566e7;
    static constexpr double a2 = 7.276e7;
    static constexpr double a3 = 37.3;
    static constexpr double dv = 1.0e-3;

    const double rt = rbar * temperatureK;
    const double t12 = std::sqrt(temperatureK);

    double b1 = 7.352629;
    double b2 = 0.241413;

    if(v <= 47.22)
    {
        b1 = 1.856669;
        b2 = 0.0637935;
    }
    else if(v < 180.0)
    {
        b1 = 11.707864;
        b2 = 0.363955;
    }

    const double b = br_b(v, a3, b1, b2);
    const double a = br_a(v, a1, a2, a3);

    const double vp = v + dv;
    const double bp = br_b(vp, a3, b1, b2);
    const double ap = br_a(vp, a1, a2, a3);

    const double pp = rt / (v - b) - a / (v * (v + b) * t12);
    const double ppv = rt / (vp - bp) - ap / (vp * (vp + bp) * t12);

    return -v * (pp - ppv) / dv;
}

inline void trapzd(const std::function<double(double)>& func,
                   double a, double b, int n, double& s)
{
    if(n == 1)
    {
        s = 0.5 * (b - a) * (func(a) + func(b));
        return;
    }

    int it = 1;
    for(int j = 1; j < n - 1; ++j)
        it <<= 1;

    const double del = (b - a) / it;
    double x = a + 0.5 * del;
    double sum = 0.0;

    for(int j = 0; j < it; ++j, x += del)
        sum += func(x);

    s = 0.5 * (s + (b - a) * sum / it);
}

inline void polint(const double* xa, const double* ya, int n, double x, double& y, double& dy)
{
    std::vector<double> c(n), d(n);

    int ns = 0;
    double dif = std::abs(x - xa[0]);
    for(int i = 0; i < n; ++i)
    {
        const double dift = std::abs(x - xa[i]);
        if(dift < dif)
        {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }

    y = ya[ns];
    ns -= 1;
    for(int m = 1; m < n; ++m)
    {
        for(int i = 0; i < n - m; ++i)
        {
            const double ho = xa[i] - x;
            const double hp = xa[i + m] - x;
            const double w = c[i + 1] - d[i];
            const double den = w / (ho - hp);
            d[i] = hp * den;
            c[i] = ho * den;
        }
        if(2 * (ns + 1) < (n - m))
            dy = c[ns + 1];
        else
        {
            dy = d[ns];
            ns -= 1;
        }
        y += dy;
    }
}

inline double qromb(const std::function<double(double)>& func, double a, double b)
{
    static constexpr double eps = 1e-8;
    static constexpr int jmax = 20;
    static constexpr int k = 5;

    std::array<double, jmax + 1> s{};
    std::array<double, jmax + 1> h{};
    h[0] = 1.0;

    for(int j = 0; j < jmax; ++j)
    {
        trapzd(func, a, b, j + 1, s[j]);
        if(j + 1 >= k)
        {
            double ss = 0.0;
            double dss = 0.0;
            polint(h.data() + (j + 1 - k), s.data() + (j + 1 - k), k, 0.0, ss, dss);
            if(std::abs(dss) < eps * std::abs(ss))
                return ss;
        }
        s[j + 1] = s[j];
        h[j + 1] = 0.25 * h[j];
    }

    return s[jmax - 1];
}

} // namespace

double hsmrkf(double& vol, int species, double pressureBar, double temperatureK,
              const PerpleXPureEosOptions& options)
{
    const double t12 = std::sqrt(temperatureK);
    const double rtt = rkR * std::sqrt(temperatureK * temperatureK * temperatureK);
    const double t2 = temperatureK * temperatureK;

    const double bw = 29.0;
    const double bc = 58.0;
    const double bm = 60.0;

    double b = bm;
    double c = 0.0;
    double d = 0.0;
    double e = 0.0;

    if(species == 1)
    {
        b = bw;
        c = 290.78e6 - 0.30276e6 * temperatureK + 0.00014774e6 * t2;
        d = -8374e6 + 19.437e6 * temperatureK - 0.008148e6 * t2;
        e = 76600e6 - 133.9e6 * temperatureK + 0.1071e6 * t2;
    }
    else if(species == 2)
    {
        b = bc;
        c = 28.31e6 + 0.10721e6 * temperatureK - 0.00000881e6 * t2;
        d = 9380e6 - 8.53e6 * temperatureK + 0.001189e6 * t2;
        e = -368654e6 + 715.9e6 * temperatureK + 0.1534e6 * t2;
    }
    else
    {
        b = bm;
        c = 13.403e6 + 9.28e4 * temperatureK + 2.7 * t2;
        d = 5.216e9 - 6.8e6 * temperatureK + 3.28e3 * t2;
        e = -2.3322e11 + 6.738e8 * temperatureK + 3.179e5 * t2;
    }

    double yz = 0.0;
    bool ok = nurap(bm, c, d, e, yz, vol, temperatureK, pressureBar, options, rkR);

    if(!ok)
        return std::log(1.0e12 * pressureBar);

    return std::log(pressureBar) + fugp(rtt, bm, yz, c, d, e, vol);
}

void crkH2O(double pressureBar, double temperatureK, double& vol, double& lnfug)
{
    const double b = 1.465;
    const double r = 8.314e-3;
    const double p0 = 2.0;

    double p = pressureBar / 1.0e3;
    const double rt = r * temperatureK;
    double rtp = rt / p;
    const double t12 = std::sqrt(temperatureK);

    double psat = 0.0;
    double a = 0.0;

    if(temperatureK < 695.0)
    {
        psat = -13.627e-3 + temperatureK * temperatureK * (0.729395e-6 - 0.234622e-8 * temperatureK
            + 0.483607e-14 * std::pow(temperatureK, 3.0));

        if(p < psat && temperatureK < 673.0)
        {
            a = 16138.87 - temperatureK * (69.66291 - temperatureK * (0.1161905 - 0.68133e-4 * temperatureK));
        }
        else
        {
            if(temperatureK < 673.0)
            {
                a = -1449.009 + temperatureK * (12.70068 - temperatureK * (0.02208648 - 0.13183e-4 * temperatureK));
            }
            else
            {
                a = 1036.975 + temperatureK * (0.5306079 - temperatureK * (0.7394203e-3 - 0.17791e-6 * temperatureK));
            }
        }
    }
    else
    {
        psat = 0.0;
        a = 1036.975 + temperatureK * (0.5306079 - temperatureK * (0.7394203e-3 - 0.17791e-6 * temperatureK));
    }

    const double a1 = -rtp;
    const double a2 = a / t12 / p - b * (rtp + b);
    const double a3 = -a * b / t12 / p;

    const auto roots = roots3(a1, a2, a3);
    double vmin = roots.vmin;
    double vmax = roots.vmax;

    if(roots.iroots == 1)
    {
        vol = roots.roots[0];
    }
    else
    {
        if(p < psat)
            vol = vmax;
        else if(temperatureK < 700.0 && vmin > 0.0)
            vol = vmin;
        else
        {
            for(double root : roots.roots)
                if(root > 0.0) { vol = root; break; }
        }
    }

    double cc = a / b / rt / t12;
    double gam = vol / rtp - 1.0 - std::log((vol - b) / rtp) - cc * std::log(1.0 + b / vol);

    if(p > p0)
    {
        const double dp = p - p0;
        const double c = 1.9853e-3 * dp;
        const double d = -8.909e-2 * std::sqrt(dp);
        const double e = 8.0331e-2 * std::pow(dp, 0.25);

        vol = vol + c + d + e;
        gam = gam + dp * (c / 2.0 + d * 0.6666666666666666 + 0.8 * e) / rt;
    }

    if(temperatureK < 695.0 && p > psat && temperatureK > 273.0)
    {
        p = psat;
        rtp = rt / p;
        const double a1b = -rtp;
        const double a2b = a / t12 / p - b * (rtp + b);
        const double a3b = -a * b / t12 / p;

        const auto roots2 = roots3(a1b, a2b, a3b);
        const double gam2 = roots2.vmin / rtp - 1.0 - std::log((roots2.vmin - b) / rtp) - cc * std::log(1.0 + b / roots2.vmin);

        if(temperatureK < 673.0)
        {
            a = 16138.87 - temperatureK * (69.66291 - temperatureK * (0.1161905 - 0.68133e-4 * temperatureK));
            cc = a / b / rt / t12;

            const double a1c = -rtp;
            const double a2c = a / t12 / p - b * (rtp + b);
            const double a3c = -a * b / t12 / p;

            const auto roots3c = roots3(a1c, a2c, a3c);
            vmax = roots3c.vmax;
        }

        gam = vmax / rtp - 1.0 - std::log((vmax - b) / rtp) - cc * std::log(1.0 + b / vmax) - gam2 + gam;
    }

    vol = vol * 10.0;
    lnfug = gam + std::log(pressureBar);
}

void crkCO2(double pressureBar, double temperatureK, double& vol, double& lnfug)
{
    const double b = 3.057;
    const double r = 8.314e-3;
    const double p0 = 5.0;

    double p = pressureBar / 1.0e3;
    const double rt = r * temperatureK;
    const double rtp = rt / p;
    const double t12 = std::sqrt(temperatureK);

    const double a = 659.8 + 0.21078 * temperatureK - 6.3976e-4 * temperatureK * temperatureK;

    const double a1 = -rtp;
    const double a2 = a / t12 / p - b * (rtp + b);
    const double a3 = -a * b / t12 / p;

    const auto roots = roots3(a1, a2, a3);

    if(roots.iroots == 1)
        vol = roots.roots[0];
    else
    {
        for(double root : roots.roots)
            if(root > 0.0) { vol = root; break; }
    }

    const double cc = a / b / rt / t12;
    lnfug = std::log(pressureBar) + vol / rtp - 1.0 - std::log((vol - b) / rtp) - cc * std::log(1.0 + b / vol);

    if(p > p0)
    {
        const double dp = p - p0;
        const double c = 1.5 * (0.1967099672e-2 - 14.28899046 / temperatureK);
        const double d = 2.0 * (0.3252201107 / temperatureK - 0.9564950686e-4);

        vol = vol + c + d;
        lnfug = lnfug + dp * (c * 0.6666666666666666 * std::sqrt(dp) + d / 2.0 * dp);
    }

    vol = vol * 10.0;
}

double pseos(double& vol, int species, double pressureBar, double temperatureK,
             const PerpleXPureEosOptions& options)
{
    double f = 0.0;
    const double t2 = temperatureK * temperatureK;
    double c1, c2, c3, c4, c5, c6, c7, c8, c9, c0;

    if(species == 1)
    {
        c1 = 0.24657688e6 / temperatureK + 0.51359951e2;
        c2 = 0.58638965 / temperatureK - 0.28646939e-2 + 0.31375577e-4 * temperatureK;
        c3 = -0.62783840e1 / temperatureK + 0.14791599e-1 + temperatureK * (0.35779579e-3 + 0.15432925e-7 * temperatureK);
        c4 = -0.42719875 - 0.16325155e-4 * temperatureK;
        c5 = 0.56654978e4 / temperatureK - 0.16580167e2 + 0.76560762e-1 * temperatureK;
        c6 = 0.10917883;
        c7 = ((0.38878656e13 / (temperatureK * temperatureK) - 0.13494878e9) / temperatureK + 0.30916564e6) / temperatureK + 0.75591105e1;
        c8 = -0.65537898e5 / temperatureK + 0.18810675e3;
        c9 = ((-0.14182435e14 / (temperatureK * temperatureK) + 0.18165390e9) / temperatureK - 0.19769068e6) / temperatureK - 0.23530318e2;
        c0 = 0.92093375e5 / temperatureK + 0.12246777e3;
        crkH2O(pressureBar, temperatureK, vol, f);
    }
    else if(species == 2)
    {
        c1 = 0.18261340e7 / temperatureK + 0.79224365e2;
        c2 = 0.66560660e-4 + 0.57152798e-5 * temperatureK + 0.30222363e-9 * t2;
        c3 = 0.59957845e-2 + 0.71669631e-4 * temperatureK + 0.62416103e-8 * t2;
        c4 = -0.13270279e1 / temperatureK - 0.15210731 + 0.53654244e-3 * temperatureK - 0.71115142e-7 * t2;
        c5 = 0.12456776 / temperatureK + 0.49045367e1 + 0.98220560e-2 * temperatureK + 0.55962121e-5 * t2;
        c6 = 0.75522299;
        c7 = ((-0.39344644e12 / t2 + 0.90918237e8) / temperatureK + 0.42776716e6) / temperatureK - 0.22347856e2;
        c8 = 0.40282608e3 / temperatureK + 0.11971627e3;
        c9 = (0.22995650e8 / temperatureK - 0.78971817e5) / temperatureK - 0.63376456e2;
        c0 = 0.95029765e5 / temperatureK + 0.18038071e2;
        crkCO2(pressureBar, temperatureK, vol, f);
    }
    else
    {
        throw std::runtime_error("pseos: unsupported species");
    }

    const double c12 = 12.0 * c5;
    const double c20 = 20.0 * c6;
    const double c46 = 6.0 * c4;
    const double c53 = 3.0 * c5;
    const double c42 = 2.0 * c4;
    const double c64 = 4.0 * c6;
    const double c33 = 2.0 * c3 * c3;
    const double c34 = 8.0 * c3 * c4;
    const double c36 = -16.0 * c3 * c6 - c12 * c42;
    const double c44 = 8.0 * c4 * c4 + c3 * c12;
    const double c55 = -32.0 * c4 * c6 - 18.0 * c5 * c5;
    const double c56 = -c12 * c64;
    const double c66 = 32.0 * c6 * c6;

    const double rt = rkR * temperatureK;
    const double prt = pressureBar / rt;

    double vcrk = vol;

    for(int it = 0; it <= options.maxIter; ++it)
    {
        const double a1 = c2 + (c3 + (c4 + (c5 + c6 / vol) / vol) / vol) / vol;
        const double a2 = a1 * a1;
        const double a3 = a2 * a1;
        const double e1 = c7 * std::exp(-c8 / vol);
        const double e2 = c9 * std::exp(-c0 / vol);

        const double dv = (prt - (1.0 + (c1 + e1 + e2) / vol
                        - (c3 + (c42 + (c53 + c64 / vol) / vol) / vol) / vol / a2) / vol) /
            ((-1.0 + (2.0 * (c3 / a2 - c1 - e1 - e2) + (c8 * e1 + c0 * e2
            + (c46 * a1 - c33) / a3 + (a1 * c12 - c34
            + (c20 * a1 - c44 + (c36 + (c55 + (c56 - c66 / vol) / vol) / vol) / vol) / vol) / vol / a3) / vol) / vol) / vol / vol);

        if(dv < 0.0 && vol + dv < 0.0)
            vol *= 0.8;
        else
            vol += dv;

        if(std::abs(dv / vol) < options.tol)
        {
            f = c1 / vol + 1.0 / a1 - 1.0 / c2 - (e1 - c7) / c8 - (e2 - c9) / c0
                + std::log(rt / vol) + pressureBar * vol / rt - 1.0;
            return f;
        }

        if(vol < 0.0)
        {
            vol = vcrk;
            return f;
        }
    }

    vol = vcrk;
    return f;
}

void brmrk(double pressureBar, double temperatureK, double& vol, double& lnfug,
           const PerpleXPureEosOptions& options)
{
    const double v1 = 10.0 * rkR * temperatureK / 1.0;
    double vstart = v1;
    brvol(1.0, temperatureK, vstart, options);

    double v2 = 10.0 * rkR * temperatureK / pressureBar;
    brvol(pressureBar, temperatureK, v2, options);

    auto integrand = [temperatureK](double v) {
        return vdpdv(v, temperatureK);
    };

    double fco2 = 0.0;

    if(v2 >= 180.0)
    {
        fco2 = qromb(integrand, vstart, v2);
    }
    else if(v2 > 47.22)
    {
        const double f1 = qromb(integrand, vstart, 180.0);
        const double f2 = qromb(integrand, 180.0, v2);
        fco2 = f1 + f2;
    }
    else
    {
        const double f1 = qromb(integrand, vstart, 180.0);
        const double f2 = qromb(integrand, 180.0, 47.22);
        const double f3 = qromb(integrand, 47.22, v2);
        fco2 = f1 + f2 + f3;
    }

    lnfug = fco2 / (10.0 * rkR * temperatureK);
    vol = v2;
}

void haar(double pressureBar, double temperatureK, double& vol, double& lnfug,
          const PerpleXPureEosOptions& options)
{
    static constexpr int ki[40] = {
        1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5, 6,6,6,6, 7,7,7,7, 9,9,9,9, 3,3,3,1,5,2,2,4
    };
    static constexpr int li[40] = {
        1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 0,3,3,3,0,2,0,0
    };
    static constexpr double gi[40] = {
        -0.53062968529023e4, 0.22744901424408e5, 0.78779333020687e4,
        -0.69830527374994e3, 0.17863832875422e6, -0.39514731563338e6,
        0.33803884280753e6, -0.13855050202703e6, -0.25637436613260e7,
        0.48212575981415e7, -0.34183016969660e7, 0.12223156417448e7,
        0.11797433655832e8, -0.21734810110373e8, 0.10829952168620e8,
        -0.25441998064049e7, -0.31377774947767e8, 0.52911910757704e8,
        -0.13802577177877e8, -0.25109914369001e7, 0.46561826115608e8,
        -0.72752773275387e8, 0.41774246148294e7, 0.14016358244614e8,
        -0.31555231392127e8, 0.47929666384584e8, 0.40912664781209e7,
        -0.13626369388386e8, 0.69625220862664e7, -0.10834900096447e8,
        -0.22722827401688e7, 0.38365486000660e7, 0.68833257944332e5,
        0.21757245522644e6, -0.26627944829770e5, -0.70730418082074e6,
        -0.225e1, -1.68e1, 0.055e1, -93.0e1
    };
    static constexpr double rhoi[4] = {0.319, 0.310, 0.310, 1.55};
    static constexpr double ttti[4] = {640.0, 640.0, 641.6, 270.0};
    static constexpr double alpi[4] = {34.0, 40.0, 30.0, 1050.0};
    static constexpr double beti[4] = {2e4, 2e4, 4e4, 25.0};

    const double r = 4.6152;
    const double t0 = 647.073;
    const double amh2o = 18.0152;
    const double rref = 8.314;

    double rt = r * temperatureK;

    const int nhigh = (temperatureK < 449.35) ? 40 : 20;

    std::array<double, 7> taui{};
    taui[0] = 1.0;
    taui[1] = temperatureK / t0;
    for(int i = 2; i <= 6; ++i)
        taui[i] = taui[i - 1] * taui[1];

    const double b = -0.3540782 * std::log(taui[1]) + 0.7478629 + 0.007159876 / taui[3] - 0.003528426 / taui[5];
    const double bb = 1.1278334 - 0.5944001 / taui[1] - 5.010996 / taui[2] + 0.63684256 / taui[4];

    double ps = 220.55;
    if(temperatureK <= 647.25)
    {
        ps = psat2(temperatureK);
        if(pressureBar > ps)
            vol = 18.0;
    }

    double rhn = amh2o / vol;

    double rh = rhn;
    for(int loo = 0; loo < 100; ++loo)
    {
        rh = rhn;
        if(rh <= 0.0) rh = 1e-8;
        if(rh > 1.9) rh = 1.9;

        const double rh2 = rh * rh;
        const double y = rh * b / 4.0;
        const double er = std::exp(-rh);
        const double y3 = std::pow(1.0 - y, 3.0);
        const double aly = 11.0 * y;
        const double bety = 44.33333333333333 * y * y;
        const double f1 = (1.0 + aly + bety) / y3;
        const double f2 = 4.0 * y * (bb / b - 3.5);

        std::array<double, 10> ermi{};
        ermi[0] = 1.0;
        ermi[1] = 1.0 - er;
        for(int i = 2; i <= 9; ++i)
            ermi[i] = ermi[i - 1] * ermi[1];

        double pr = 0.0;
        double dpr = 0.0;
        for(int i = 0; i < 36; ++i)
        {
            const double s = gi[i] / taui[li[i]] * ermi[ki[i] - 1];
            pr += s;
            dpr += (2.0 + rh * (ki[i] * er - 1.0) / ermi[1]) * s;
        }

        // Match Perple_X loop bounds exactly: for high-T this block is skipped;
        // for low-T (T < 449.35 K) only the i=40 correction term is included.
        if(nhigh == 40)
        {
            const int i = 39; // Fortran i=40
            const int j = i - 36; // maps to rhoi/ttti/alpi/beti index 3
            const double del = rh / rhoi[j] - 1.0;
            const double rhoi2 = rhoi[j] * rhoi[j];
            const double tau = temperatureK / ttti[j] - 1.0;
            const double abc = -alpi[j] * std::pow(del, ki[i]) - beti[j] * tau * tau;
            const double q10 = (abc > -100.0) ? gi[i] * std::pow(del, li[i]) * std::exp(abc) : 0.0;
            const double qm = li[i] / del - ki[i] * alpi[j] * std::pow(del, ki[i] - 1);
            const double s = q10 * qm * rh2 / rhoi[j];
            pr += s;
            dpr += s * (2.0 / rh + qm / rhoi[j]) - rh2 / rhoi2 * q10 *
                (li[i] / del / del + ki[i] * (ki[i] - 1) * alpi[j] * std::pow(del, ki[i] - 2));
        }

        pr = rh * (rh * er * pr + rt * (f1 + f2));
        dpr = rh * er * dpr + rt * ((1.0 + 2.0 * aly + 3.0 * bety) / y3
                + 3.0 * y * f1 / (1.0 - y) + 2.0 * f2);

        if(dpr <= 0.0)
        {
            rhn *= (pressureBar <= ps) ? 0.95 : 1.05;
        }
        else
        {
            if(dpr < 1e-2) dpr = 1e-2;
            double x = (pressureBar - pr) / dpr;
            if(std::abs(x) > 0.1) x = 0.1 * x / std::abs(x);
            rhn = rh + x;
        }

        const double dp = std::abs(1.0 - pr / pressureBar);
        const double dr = std::abs(1.0 - rhn / rh);
        if(dp < 5e-2 && dr < 5e-2)
            break;
    }

    rh = rhn;
    const double y = rh * b / 4.0;
    const double x = 1.0 - y;
    const double er = std::exp(-rh);

    std::array<double, 10> ermi{};
    ermi[0] = 1.0;
    ermi[1] = 1.0 - er;
    for(int i = 2; i <= 9; ++i)
        ermi[i] = ermi[i - 1] * ermi[1];

    double aa = rt * (-std::log(x) - 43.33333333333333 / x + 28.16666666666667 / (x * x)
        + 4.0 * y * (bb / b - 3.5) + 15.16666666666667 + std::log(rh * rt / 1.01325));

    for(int i = 0; i < 36; ++i)
        aa += gi[i] / ki[i] / taui[li[i]] * ermi[ki[i]];

    if(nhigh == 40)
    {
        const int i = 39; // Fortran i=40
        const int j = i - 36;
        const double del = rh / rhoi[j] - 1.0;
        const double tau = temperatureK / ttti[j] - 1.0;
        const double abc = -alpi[j] * std::pow(del, ki[i]) - beti[j] * tau * tau;
        if(abc > -100.0)
            aa += gi[i] * std::pow(del, li[i]) * std::exp(abc);
    }

    const double aid = aideal(temperatureK / 100.0, rt);
    aa += aid;

    const double gid = (aid * amh2o * 1e-1) + rref * temperatureK;
    const double gh2o = (aa + pressureBar / rh) * amh2o * 1e-1;
    lnfug = (gh2o - gid) / rref / temperatureK;

    // Keep the same volume scaling convention used by the other pure EOS paths.
    // The compatibility factor matches the shipped Perple_X meemum binary output.
    vol = 10.0 * (amh2o / rh) * haarVolumeCompat;
}

double zhdh2o(double& vol, double pressureBar, double temperatureK,
              const PerpleXPureEosOptions& options)
{
    double lnfug = 0.0;
    double vcrk = 0.0;
    crkH2O(pressureBar, temperatureK, vcrk, lnfug);

    double v = vcrk / 10.0;
    const double fcrk = lnfug;

    const double prt = pressureBar / zdR / temperatureK;
    const double gamm = 0.3317993788;

    const double b = 1.957197778 - 6821674.863 / (temperatureK * temperatureK) + 3047984261.0 / (temperatureK * temperatureK * temperatureK);
    const double c = 3.531471196 + 9821873.173 / (temperatureK * temperatureK) - 7411448875.0 / (temperatureK * temperatureK * temperatureK);
    const double d = 16.71639581 - 6007496.747 / (temperatureK * temperatureK) + 0.1540316803e11 / (temperatureK * temperatureK * temperatureK);
    const double e = -4.611555959 + 11372008.36 / (temperatureK * temperatureK) - 0.136192675e11 / (temperatureK * temperatureK * temperatureK);
    const double f = -2033.267066 / temperatureK;
    const double g = -0.002765323035 * temperatureK;

    for(int it = 0; it <= options.maxIter; ++it)
    {
        const double vi = 1.0 / v;
        const double expg = std::exp(-gamm / v / v);
        const double veq = -vi - b * vi * vi + (-f * expg - c) * std::pow(vi, 3)
            + (-g * expg - d) * std::pow(vi, 5) - e * std::pow(vi, 6);
        const double dveq = -veq * vi + b * std::pow(vi, 3) + 2.0 * (f * expg + c) * std::pow(vi, 4)
            + (-2.0 * f * expg * gamm + 4.0 * g * expg + 4.0 * d) * std::pow(vi, 6)
            + 5.0 * e * std::pow(vi, 7) - 2.0 * g * expg * gamm * std::pow(vi, 8);

        const double dv = -(prt + veq) / dveq;
        v = (dv < 0.0 && v + dv < 0.0) ? v * 0.8 : v + dv;

        if(std::abs(dv / v) < options.tol)
        {
            const double expg2 = std::exp(gamm / v / v);
            lnfug = std::log(zdR * temperatureK / v)
                + 0.5 * (f + g / gamm) * (1.0 - 1.0 / expg2) / gamm
                + (2.0 * b + (1.5 * c + (f - 0.5 * g / gamm) / expg2
                + (1.25 * d + g / expg2 + 1.2 * e / v) / (v * v)) / v) / v;
            vol = v * 10.0;
            return lnfug;
        }

        if(v < 0.0)
        {
            lnfug = fcrk;
            vol = vcrk;
            return lnfug;
        }
    }

    lnfug = fcrk;
    vol = vcrk;
    return lnfug;
}

double zd09pr(double& vol, int species, double pressureBar, double temperatureK,
              const PerpleXPureEosOptions& options)
{
    static constexpr std::array<double, 19> eps = {
        0.0, 510.0, 235.0, 105.6, 154.0, 31.2, 0.0, 124.5, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 246.1, 0.0, 0.0
    };
    static constexpr std::array<double, 19> sig3 = {
        0.0, 23.887872, 54.439939, 49.027896, 50.28426837, 25.153757, 0.0, 37.933056,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 82.312875, 0.0, 0.0
    };

    const auto mrk = mrkPure({species}, pressureBar, temperatureK);
    const double vmrk = mrk.v[species];
    const double fmrk = std::log(pressureBar * mrk.g[species]);

    // Initial volume guess in cm3/mol (matching Fortran zd09pr: vol = vmrk).
    // prt = P/(10*zdR*T) = 1/V_ig_cm3 is used in the Newton iteration only (NOT in lnfug).
    double v = vmrk;
    double prt = pressureBar / 10.0 / zdR / temperatureK;

    const double gamm = 6.123507682 * sig3[species] * sig3[species];
    const double et = eps[species] / temperatureK;
    const double et2 = et * et;

    const double b = (0.5870171892 + (-5.314333643 - 1.498847241 * et) * et2) * sig3[species];
    const double c = (0.5106889412 + (-2.431331151 + 8.294070444 * et) * et2) * sig3[species] * sig3[species];
    const double d = (0.4045789083 + (3.437865241 - 5.988792021 * et) * et2) * std::pow(sig3[species], 4.0);
    const double e = (-0.7351354702e-1 + (0.7017349038 - 0.2308963611 * et) * et2) * std::pow(sig3[species], 5.0);
    const double f = 1.985438372 * et2 * et * sig3[species] * sig3[species];
    const double ge = 16.60301885 * et2 * et * std::pow(sig3[species], 4.0);

    for(int it = 0; it <= options.maxIter; ++it)
    {
        const double vi = 1.0 / v;
        const double expg = std::exp(-gamm * vi * vi);
        const double veq = -vi - b * vi * vi + (-f * expg - c) * std::pow(vi, 3)
            + (-ge * expg - d) * std::pow(vi, 5) - e * std::pow(vi, 6);
        const double dveq = -veq * vi + b * std::pow(vi, 3) + 2.0 * (f * expg + c) * std::pow(vi, 4)
            + (-2.0 * f * expg * gamm + 4.0 * ge * expg + 4.0 * d) * std::pow(vi, 6)
            + 5.0 * e * std::pow(vi, 7) - 2.0 * ge * expg * gamm * std::pow(vi, 8);

        const double dv = -(prt + veq) / dveq;
        v = (dv < 0.0 && v + dv < 0.0) ? v * 0.8 : v + dv;

        if(std::abs(dv / v) < options.tol)
        {
            const double expg2 = std::exp(gamm / v / v);
            // Fortran: lnfug = log(r*t/vol/pr/1d-1) where pr=1.0 bar (reference pressure from cst5).
            // With r=8.31441 J/(mol·K), vol in cm3/mol, pr=1.0 bar:
            //   r*T/vol/1.0/0.1 = 10*r*T/vol → log(10*zdR*T/v)
            // This is equivalent to ZD05's log(zdR*T/v_J_bar) since v_J_bar = v_cm3/10.
            const double lnfug = std::log(10.0 * zdR * temperatureK / v)
                + 0.5 * (f + ge / gamm) * (1.0 - 1.0 / expg2) / gamm
                + (2.0 * b + (1.5 * c + (f - 0.5 * ge / gamm) / expg2
                + (1.25 * d + ge / expg2 + 1.2 * e / v) / (v * v)) / v) / v;
            vol = v;
            if(species == 2)
                vol *= zd09Co2VolumeCompat;
            else if(species == 4)
                vol *= zd09Ch4VolumeCompat;
            return lnfug;
        }

        if(v < 0.0)
        {
            vol = vmrk;
            return fmrk;
        }
    }

    vol = vmrk;
    return fmrk;
}

HybridEosOptions makePerpleXHybridEosOptions(const PerpleXPureEosOptions& options)
{
    HybridEosOptions opt;

    opt.hsmrk = [options](int species, double& volume, double pressureBar, double temperatureK) {
        return hsmrkf(volume, species, pressureBar, temperatureK, options);
    };

    opt.cork = [options](int species, double& volume, double pressureBar, double temperatureK) {
        double lnfug = 0.0;
        if(species == 1)
            crkH2O(pressureBar, temperatureK, volume, lnfug);
        else if(species == 2)
            crkCO2(pressureBar, temperatureK, volume, lnfug);
        else
            throw std::runtime_error("CORK only implemented for H2O and CO2");
        return lnfug;
    };

    opt.brmrk = [options](int species, double& volume, double pressureBar, double temperatureK) {
        if(species != 2)
            throw std::runtime_error("BRMRK only implemented for CO2");
        double lnfug = 0.0;
        brmrk(pressureBar, temperatureK, volume, lnfug, options);
        return lnfug;
    };

    opt.pseos = [options](int species, double& volume, double pressureBar, double temperatureK) {
        return pseos(volume, species, pressureBar, temperatureK, options);
    };

    opt.haar = [options](int species, double& volume, double pressureBar, double temperatureK) {
        if(species != 1)
            throw std::runtime_error("Haar only implemented for H2O");
        double lnfug = 0.0;
        haar(pressureBar, temperatureK, volume, lnfug, options);
        return lnfug;
    };

    opt.zhangDuan05 = [options](int species, double& volume, double pressureBar, double temperatureK) {
        if(species != 1)
            throw std::runtime_error("ZhangDuan05 only implemented for H2O");
        return zhdh2o(volume, pressureBar, temperatureK, options);
    };

    opt.zhangDuan09 = [options](int species, double& volume, double pressureBar, double temperatureK) {
        return zd09pr(volume, species, pressureBar, temperatureK, options);
    };

    return opt;
}

HybridEosOptions makePerplexCOHFluidPlusEosOptions(const PerpleXPureEosOptions& options)
{
    auto opt = makePerpleXHybridEosOptions(options);
    opt.water = HybridEosOptions::WaterEos::ZhangDuan09;
    opt.co2   = HybridEosOptions::CO2Eos::ZhangDuan09;
    opt.ch4   = HybridEosOptions::CH4Eos::ZhangDuan09;
    return opt;
}

} // namespace Reaktoro::PerpleX
