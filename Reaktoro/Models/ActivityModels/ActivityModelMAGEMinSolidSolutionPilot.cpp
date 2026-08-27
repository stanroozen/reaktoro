// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "ActivityModelMAGEMinSolidSolutionPilot.hpp"

// C++ includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

// Third-party includes
#include <Eigen/Dense>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHollandPowell.hpp>

#ifdef _MSC_VER
// Workaround for MSVC C2668: 'fpclassify' ambiguous overload resolution when
// Eigen's isFinite/hasNaN member functions are instantiated with autodiff::real1st
// scalar types. Adding overloads in autodiff::detail enables ADL resolution.
namespace autodiff { namespace detail {
template<size_t N, typename T> inline int  fpclassify(const Real<N,T>& x) noexcept { return std::fpclassify(static_cast<double>(x)); }
template<size_t N, typename T> inline bool isinf     (const Real<N,T>& x) noexcept { return std::isinf    (static_cast<double>(x)); }
template<size_t N, typename T> inline bool isnan     (const Real<N,T>& x) noexcept { return std::isnan    (static_cast<double>(x)); }
template<size_t N, typename T> inline bool isfinite  (const Real<N,T>& x) noexcept { return std::isfinite (static_cast<double>(x)); }
}} // namespace autodiff::detail
#endif

namespace Reaktoro {

namespace {

using std::log;

constexpr auto CompositionFloor = 1.0e-12;

constexpr auto CandidateSeedTolerance = 1.0e-8;

constexpr auto ProjectedGradientArmijo = 1.0e-4;

constexpr auto ProjectedGradientAgreementTolerance = 1.0e-6;

constexpr auto BranchStabilityObjectiveTolerance = 1.0e-3;

constexpr auto BranchStabilitySeedGapTolerance = 5.0e-2;

constexpr auto NativeGradientFiniteDifferenceStep = 1.0e-7;

constexpr auto SplitCandidateStatesKey = "MAGEMinSolidSolutionPilot::SplitCandidates";

constexpr auto SplitCandidateObjectiveGapKey = "MAGEMinSolidSolutionPilot::CompetingStableBranchObjectiveGap";

constexpr auto SplitCandidateCountKey = "MAGEMinSolidSolutionPilot::CompetingStableBranchCount";

constexpr auto InternalObjectiveKey = "MAGEMinSolidSolutionPilot::InternalObjective";

constexpr auto PrecomputedCandidateIndexKey = "MAGEMinSolidSolutionPilot::PrecomputedCandidateIndex";

constexpr auto BuiltinLegacyMinimizerStrategy = "legacy";

constexpr auto BuiltinProjectedGradientMinimizerStrategy = "projected-gradient";

struct ConstrainedTernaryMinimizationOutcome
{
    GlobalizedSolidSolutionInternalResult result;
    Map<String, Any> extra;
};

auto sb11OlivineThermo() -> MAGEMinImportedBinarySolutionThermoModel
{
    MAGEMinImportedBinarySolutionThermoModel thermo;
    thermo.modelId = "sb11_ol";
    thermo.endmember0 = "fo";
    thermo.endmember1 = "fa";
    thermo.W = 7813.22;
    thermo.idealSiteMultiplicity = 2.0;
    return thermo;
}

auto sb11WadsleyiteThermo() -> MAGEMinImportedBinarySolutionThermoModel
{
    MAGEMinImportedBinarySolutionThermoModel thermo;
    thermo.modelId = "sb11_wa";
    thermo.endmember0 = "mgwa";
    thermo.endmember1 = "fewa";
    thermo.W = 16747.18;
    thermo.idealSiteMultiplicity = 2.0;
    return thermo;
}

auto sb21SpinelThermo() -> MAGEMinImportedBinarySolutionThermoModel
{
    MAGEMinImportedBinarySolutionThermoModel thermo;
    thermo.modelId = "sb21_sp";
    thermo.endmember0 = "sp";
    thermo.endmember1 = "hc";
    thermo.W = -533.21;
    thermo.idealSiteMultiplicity = 2.0;
    return thermo;
}

auto sb21NALThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb21_nal";
    thermo.endmembers = {"mnal", "fnal", "nnal"};
    thermo.W01 = 0.0;
    thermo.W02 = -60781.47;
    thermo.W12 = -60781.47;
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        const auto siteA = (5.0/6.0)*y[0] + (5.0/6.0)*y[1] + 0.5*y[2];
        const auto siteB = (1.0/6.0)*y[0] + (1.0/6.0)*y[1] + 0.5*y[2];
        const auto sum = y[0] + y[1] + y[2];

        return universalGasConstant * T * (
            2.0*y[0]*log(y[0])
            + sum*log(sum)
            + 6.0*siteA*log(siteA)
            + 6.0*siteB*log(siteB)
            + 2.0*y[1]*log(y[1])
            + 2.0*y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        const auto siteA = (5.0/6.0)*y[0] + (5.0/6.0)*y[1] + 0.5*y[2];
        const auto siteB = (1.0/6.0)*y[0] + (1.0/6.0)*y[1] + 0.5*y[2];
        const auto sum = y[0] + y[1] + y[2];

        ArrayXr ln_a(3);
        ln_a[0] = 2.0*log(y[0]) + log(sum) + 5.0*log(siteA) + log(siteB);
        ln_a[1] = log(sum) + 5.0*log(siteA) + log(siteB) + 2.0*log(y[1]);
        ln_a[2] = log(sum) + 3.0*log(siteA) + 3.0*log(siteB) + 2.0*log(y[2]);
        return ln_a;
    };
    return thermo;
}

auto sb21OPXThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb21_opx";
    thermo.endmembers = {"en", "fs", "mgts", "odi"};
    // Symmetric regular Margules — pairs (en-fs, en-mgts, en-odi, fs-mgts, fs-odi, mgts-odi)
    thermo.excessChemicalPotentials = [](ArrayXrConstRef y)
    {
        constexpr Index n = 4;
        constexpr double W[] = {0.0, 0.0, 32217.44, 0.0, 32217.44, 48370.41};
        ArrayXr mu(n);
        for(Index i = 0; i < n; ++i)
        {
            real Gex = 0.0;
            Index it = 0;
            for(Index j = 0; j < n; ++j)
            {
                const auto tmp = static_cast<double>(i == j ? 1.0 : 0.0) - static_cast<double>(y[j]);
                for(Index k = j + 1; k < n; ++k)
                {
                    const auto delta = static_cast<double>(i == k ? 1.0 : 0.0) - static_cast<double>(y[k]);
                    Gex -= tmp * delta * W[it++];
                }
            }
            mu[i] = Gex;
        }
        return mu;
    };
    // Site-mixing entropy from MAGEMin Sconfig for sb21_opx
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        return universalGasConstant * T * (
            (y[0]+y[3])*log(y[0]+y[3])
            + (y[0]+y[2])*log(y[0]+y[2])
            + 2.0*y[1]*log(y[1])
            + y[2]*log(y[2])
            + y[3]*log(y[3]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        ArrayXr ln_a(4);
        ln_a[0] = log(y[0]+y[3]) + log(y[0]+y[2]);
        ln_a[1] = 2.0*log(y[1]);
        ln_a[2] = log(y[2]) + log(y[0]+y[2]);
        ln_a[3] = log(y[0]+y[3]) + log(y[3]);
        return ln_a;
    };
    return thermo;
}

auto sb21CPXThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb21_cpx";
    thermo.endmembers = {"di", "he", "cen", "cats", "jd"};
    // Volume-fraction asymmetric Margules — v = {1, 1, 1, 3.5, 1}
    // pairs (di-he, di-cen, di-cats, di-jd, he-cen, he-cats, he-jd, cen-cats, cen-jd, cats-jd)
    thermo.excessChemicalPotentials = [](ArrayXrConstRef y)
    {
        constexpr Index n = 5;
        const ArrayXr volumes = (ArrayXr(n) << 1.0, 1.0, 1.0, 3.5, 1.0).finished();
        constexpr double W[] = {0.0, 24740.0, 26000.0, 24300.0, 24740.0, 26000.0, 24300.0, 60132.81, 46046.07, 10000.0};

        const auto sumv = static_cast<double>(y.matrix().dot(volumes.matrix()));
        ArrayXr phi(n);
        for(Index i = 0; i < n; ++i)
            phi[i] = static_cast<double>(y[i]) * static_cast<double>(volumes[i]) / sumv;

        ArrayXr mu(n);
        for(Index i = 0; i < n; ++i)
        {
            real Gex = 0.0;
            Index it = 0;
            for(Index j = 0; j < n; ++j)
            {
                const auto tmp = static_cast<double>(i == j ? 1.0 : 0.0) - static_cast<double>(phi[j]);
                for(Index k = j + 1; k < n; ++k)
                {
                    const auto delta = static_cast<double>(i == k ? 1.0 : 0.0) - static_cast<double>(phi[k]);
                    Gex -= tmp * delta * (W[it++] * 2.0 * static_cast<double>(volumes[i]) / (static_cast<double>(volumes[j]) + static_cast<double>(volumes[k])));
                }
            }
            mu[i] = Gex;
        }
        return mu;
    };
    // Site-mixing entropy from MAGEMin Sconfig for sb21_cpx
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        const auto s1 = y[0]+y[1]+y[3];
        const auto s2 = y[0]+y[1]+y[2]+0.5*y[3]+y[4];
        const auto s3 = y[0]+y[2];
        const auto s4 = y[3]+y[4];
        return universalGasConstant * T * (
            s1*log(s1) + 2.0*s2*log(s2) + s3*log(s3)
            + y[1]*log(y[1]) + y[2]*log(y[2])
            + y[3]*log(0.5*y[3]) + s4*log(s4) + y[4]*log(y[4]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        const auto s1 = y[0]+y[1]+y[3];
        const auto s2 = y[0]+y[1]+y[2]+0.5*y[3]+y[4];
        const auto s3 = y[0]+y[2];
        const auto s4 = y[3]+y[4];
        ArrayXr ln_a(5);
        ln_a[0] = log(s1) + 2.0*log(s2) + log(s3);
        ln_a[1] = log(s1) + 2.0*log(s2) + log(y[1]);
        ln_a[2] = log(y[2]) + 2.0*log(s2) + log(s3);
        ln_a[3] = log(0.5*y[3]) + log(s1) + log(s2) + log(s4);
        ln_a[4] = 2.0*log(s2) + log(s4) + log(y[4]);
        return ln_a;
    };
    return thermo;
}

auto sb21CPXTCMConstraintBridge() -> MAGEMinTCMConstraintBridge
{
    MAGEMinTCMConstraintBridge bridge;
    bridge.numConstraints = 9;
    bridge.numVariables = 5;
    bridge.constraintLowerBounds = ArrayXr::Constant(9, -1.0e12);
    bridge.constraintUpperBounds = ArrayXr::Zero(9);
    bridge.callback = [](unsigned m, double* result, unsigned n, const double* x, double* grad, void*)
    {
        if(m != 9 || n != 5)
            throw std::runtime_error("sb21_cpx TC mconstraint bridge expects m=9 and n=5.");

        // Direct transcription of MAGEMin TC callback `cpx_mtl_c` from NLopt_opt_function.c.
        result[0] = (-x[0]*x[1] - x[0]*x[3] + x[0] - x[1]*x[4] + x[1] - x[3]*x[4] + x[3] + x[4] - 1.0);
        result[1] = (x[0]*x[1] + x[0]*x[3] - x[0] + x[1]*x[4] + x[3]*x[4] - x[4]);
        result[2] = (-x[1] - x[3]);
        result[3] = (x[0]*x[2] + x[1]*x[4] - x[2] + x[3]*x[4] - x[4]);
        result[4] = (-x[0]*x[2] - x[1]*x[4] - x[3]*x[4] + x[4]);
        result[5] = (x[2] + x[3] - 1.0);
        result[6] = (-x[3]);
        result[7] = (0.5*x[1] - 1.0);
        result[8] = (-0.5*x[1]);

        if(grad)
        {
            grad[0] = -x[1] - x[3] + 1.0;
            grad[1] = -x[0] - x[4] + 1.0;
            grad[2] = 0.0;
            grad[3] = -x[0] - x[4] + 1.0;
            grad[4] = -x[1] - x[3] + 1.0;
            grad[5] = x[1] + x[3] - 1.0;
            grad[6] = x[0] + x[4];
            grad[7] = 0.0;
            grad[8] = x[0] + x[4];
            grad[9] = x[1] + x[3] - 1.0;
            grad[10] = 0.0;
            grad[11] = -1.0;
            grad[12] = 0.0;
            grad[13] = -1.0;
            grad[14] = 0.0;
            grad[15] = x[2];
            grad[16] = x[4];
            grad[17] = x[0] - 1.0;
            grad[18] = x[4];
            grad[19] = x[1] + x[3] - 1.0;
            grad[20] = -x[2];
            grad[21] = -x[4];
            grad[22] = -x[0];
            grad[23] = -x[4];
            grad[24] = -x[1] - x[3] + 1.0;
            grad[25] = 0.0;
            grad[26] = 0.0;
            grad[27] = 1.0;
            grad[28] = 1.0;
            grad[29] = 0.0;
            grad[30] = 0.0;
            grad[31] = 0.0;
            grad[32] = 0.0;
            grad[33] = -1.0;
            grad[34] = 0.0;
            grad[35] = 0.0;
            grad[36] = 0.5;
            grad[37] = 0.0;
            grad[38] = 0.0;
            grad[39] = 0.0;
            grad[40] = 0.0;
            grad[41] = -0.5;
            grad[42] = 0.0;
            grad[43] = 0.0;
            grad[44] = 0.0;
        }
    };
    return bridge;
}

auto normalizedVisibleFromNative(ArrayXr visible, Index expectedSize) -> ArrayXr
{
    if(visible.size() != expectedSize)
        throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: nativeToVisible returned an unexpected vector size.");

    visible = visible.max(CompositionFloor);
    const auto total = static_cast<double>(visible.sum());
    if(total <= 0.0)
        visible = ArrayXr::Constant(expectedSize, 1.0/static_cast<double>(expectedSize));
    else visible /= total;

    return visible;
}

auto sb21AkimotoiteNativeToVisible(ArrayXrConstRef x) -> ArrayXr
{
    if(x.size() != 2)
        throw std::runtime_error("sb21_ak native map expects two TC variables.");

    ArrayXr p(3);
    p[0] = x[1];
    p[1] = x[0]*x[1] - x[0] - x[1] + 1.0;
    p[2] = -x[0]*x[1] + x[0];
    return normalizedVisibleFromNative(std::move(p), 3);
}

auto sb21AkimotoiteVisibleToNative(ArrayXrConstRef y) -> ArrayXr
{
    if(y.size() < 3)
        throw std::runtime_error("sb21_ak visible map expects at least three visible coordinates.");

    auto visible = normalizedVisibleFromNative(ArrayXr(y.head(3)), 3);
    const auto x1 = std::clamp(static_cast<double>(visible[0]), CompositionFloor, 1.0 - CompositionFloor);
    const auto denom = std::max(CompositionFloor, 1.0 - x1);
    const auto x0 = std::clamp(static_cast<double>(visible[2] / denom), CompositionFloor, 1.0 - CompositionFloor);

    ArrayXr native(2);
    native << x0, x1;
    return native;
}

auto sb21PerovskiteNativeToVisible(ArrayXrConstRef x) -> ArrayXr
{
    if(x.size() != 4)
        throw std::runtime_error("sb21_pv native map expects four TC variables.");

    ArrayXr p(3);
    p[0] = x[0]*x[1] + x[0]*x[2] + x[0]*x[3] - x[0] - x[1] - x[2] - x[3] + 1.0;
    p[1] = -x[0]*x[1] - x[0]*x[2] - x[0]*x[3] + x[0];
    p[2] = x[2];
    return normalizedVisibleFromNative(std::move(p), 3);
}

auto sb21PerovskiteVisibleToNative(ArrayXrConstRef y) -> ArrayXr
{
    if(y.size() < 3)
        throw std::runtime_error("sb21_pv visible map expects at least three visible coordinates.");

    auto visible = normalizedVisibleFromNative(ArrayXr(y.head(3)), 3);
    const auto denom = std::max(CompositionFloor, static_cast<double>(visible[0] + visible[1]));
    const auto x0 = std::clamp(static_cast<double>(visible[1] / denom), CompositionFloor, 1.0 - CompositionFloor);
    const auto x2 = std::clamp(static_cast<double>(visible[2]), CompositionFloor, 1.0 - CompositionFloor);

    ArrayXr native(4);
    native << x0, CompositionFloor, x2, CompositionFloor;
    return native;
}

auto sb21MagnesiowustiteNativeToVisible(ArrayXrConstRef x) -> ArrayXr
{
    if(x.size() != 1)
        throw std::runtime_error("sb21_mw native map expects one TC variable.");

    ArrayXr p(3);
    p[0] = 1.0 - x[0];
    p[1] = x[0];
    p[2] = CompositionFloor;
    return normalizedVisibleFromNative(std::move(p), 3);
}

auto sb21MagnesiowustiteVisibleToNative(ArrayXrConstRef y) -> ArrayXr
{
    if(y.size() < 3)
        throw std::runtime_error("sb21_mw visible map expects at least three visible coordinates.");

    auto visible = normalizedVisibleFromNative(ArrayXr(y.head(3)), 3);
    const auto denom = std::max(CompositionFloor, static_cast<double>(visible[0] + visible[1]));

    ArrayXr native(1);
    native[0] = std::clamp(static_cast<double>(visible[1] / denom), CompositionFloor, 1.0 - CompositionFloor);
    return native;
}

auto sb21AkimotoiteTCMConstraintBridge() -> MAGEMinTCMConstraintBridge
{
    MAGEMinTCMConstraintBridge bridge;
    bridge.numConstraints = 5;
    bridge.numVariables = 2;
    bridge.constraintLowerBounds = ArrayXr::Constant(5, -1.0e12);
    bridge.constraintUpperBounds = ArrayXr::Zero(5);
    bridge.variableLowerBounds = ArrayXr::Constant(2, CompositionFloor);
    bridge.variableUpperBounds = ArrayXr::Constant(2, 1.0 - CompositionFloor);
    bridge.enforceUnityConstraint = false;
    bridge.nativeToVisible = [](ArrayXrConstRef x) { return sb21AkimotoiteNativeToVisible(x); };
    bridge.visibleToNative = [](ArrayXrConstRef y) { return sb21AkimotoiteVisibleToNative(y); };
    bridge.callback = [](unsigned m, double* result, unsigned n, const double* x, double* grad, void*)
    {
        if(m != 5 || n != 2)
            throw std::runtime_error("sb21_ak TC mconstraint bridge expects m=5 and n=2.");

        // Direct transcription of MAGEMin TC callback `aki_mtl_c`.
        result[0] = (-x[1]);
        result[1] = (-x[0]*x[1] + x[0] + x[1] - 1.0);
        result[2] = (x[0]*x[1] - x[0]);
        result[3] = (-x[1]);
        result[4] = (x[1] - 1.0);

        if(grad)
        {
            grad[0] = 0.0;
            grad[1] = -1.0;
            grad[2] = 1.0 - x[1];
            grad[3] = 1.0 - x[0];
            grad[4] = x[1] - 1.0;
            grad[5] = x[0];
            grad[6] = 0.0;
            grad[7] = -1.0;
            grad[8] = 0.0;
            grad[9] = 1.0;
        }
    };
    return bridge;
}

auto sb21PerovskiteTCMConstraintBridge() -> MAGEMinTCMConstraintBridge
{
    MAGEMinTCMConstraintBridge bridge;
    bridge.numConstraints = 7;
    bridge.numVariables = 4;
    bridge.constraintLowerBounds = ArrayXr::Constant(7, -1.0e12);
    bridge.constraintUpperBounds = ArrayXr::Zero(7);
    bridge.variableLowerBounds = ArrayXr::Constant(4, CompositionFloor);
    bridge.variableUpperBounds = ArrayXr::Constant(4, 1.0 - CompositionFloor);
    bridge.enforceUnityConstraint = false;
    bridge.nativeToVisible = [](ArrayXrConstRef x) { return sb21PerovskiteNativeToVisible(x); };
    bridge.visibleToNative = [](ArrayXrConstRef y) { return sb21PerovskiteVisibleToNative(y); };
    bridge.callback = [](unsigned m, double* result, unsigned n, const double* x, double* grad, void*)
    {
        if(m != 7 || n != 4)
            throw std::runtime_error("sb21_pv TC mconstraint bridge expects m=7 and n=4.");

        // Direct transcription of MAGEMin TC callback `mpv_mtl_c`.
        result[0] = (-x[2]);
        result[1] = (-x[0]*x[1] - x[0]*x[2] - x[0]*x[3] + x[0] + x[1] + x[2] + x[3] - 1.0);
        result[2] = (x[0]*x[1] + x[0]*x[2] + x[0]*x[3] - x[0]);
        result[3] = (-0.5*x[3]);
        result[4] = (-x[1] - 0.5*x[3]);
        result[5] = (-x[1]);
        result[6] = (x[1] - 1.0);

        if(grad)
        {
            grad[0] = 0.0;
            grad[1] = 0.0;
            grad[2] = -1.0;
            grad[3] = 0.0;
            grad[4] = -x[1] - x[2] - x[3] + 1.0;
            grad[5] = 1.0 - x[0];
            grad[6] = 1.0 - x[0];
            grad[7] = 1.0 - x[0];
            grad[8] = x[1] + x[2] + x[3] - 1.0;
            grad[9] = x[0];
            grad[10] = x[0];
            grad[11] = x[0];
            grad[12] = 0.0;
            grad[13] = 0.0;
            grad[14] = 0.0;
            grad[15] = -0.50;
            grad[16] = 0.0;
            grad[17] = -1.0;
            grad[18] = 0.0;
            grad[19] = -0.50;
            grad[20] = 0.0;
            grad[21] = -1.0;
            grad[22] = 0.0;
            grad[23] = 0.0;
            grad[24] = 0.0;
            grad[25] = 1.0;
            grad[26] = 0.0;
            grad[27] = 0.0;
        }
    };
    return bridge;
}

auto sb21PostPerovskiteTCMConstraintBridge() -> MAGEMinTCMConstraintBridge
{
    auto bridge = sb21PerovskiteTCMConstraintBridge();
    // `cpv_mtl_c` is algebraically identical to `mpv_mtl_c` in MAGEMin TC callbacks.
    return bridge;
}

auto sb21MagnesiowustitesTCMConstraintBridge() -> MAGEMinTCMConstraintBridge
{
    MAGEMinTCMConstraintBridge bridge;
    bridge.numConstraints = 2;
    bridge.numVariables = 1;
    bridge.constraintLowerBounds = ArrayXr::Constant(2, -1.0e12);
    bridge.constraintUpperBounds = ArrayXr::Zero(2);
    bridge.variableLowerBounds = ArrayXr::Constant(1, CompositionFloor);
    bridge.variableUpperBounds = ArrayXr::Constant(1, 1.0 - CompositionFloor);
    bridge.enforceUnityConstraint = false;
    bridge.nativeToVisible = [](ArrayXrConstRef x) { return sb21MagnesiowustiteNativeToVisible(x); };
    bridge.visibleToNative = [](ArrayXrConstRef y) { return sb21MagnesiowustiteVisibleToNative(y); };
    bridge.callback = [](unsigned m, double* result, unsigned n, const double* x, double* grad, void*)
    {
        if(m != 2 || n != 1)
            throw std::runtime_error("sb21_mw TC mconstraint bridge expects m=2 and n=1.");

        // Direct transcription of MAGEMin TC callback `fp_mtl_c`.
        result[0] = (x[0] - 1.0);
        result[1] = (-x[0]);

        if(grad)
        {
            grad[0] = 1.0;
            grad[1] = -1.0;
        }
    };
    return bridge;
}

auto hpIGOPXLowerBounds() -> ArrayXr
{
    return (ArrayXr(8) <<
        CompositionFloor,
        CompositionFloor,
        CompositionFloor,
        -1.0 + CompositionFloor,
        CompositionFloor,
        CompositionFloor,
        CompositionFloor,
        CompositionFloor).finished();
}

auto hpIGOPXUpperBounds() -> ArrayXr
{
    return (ArrayXr(8) <<
        1.0 - CompositionFloor,
        2.0 - CompositionFloor,
        1.0 - CompositionFloor,
        1.0 - CompositionFloor,
        1.0 - CompositionFloor,
        1.0 - CompositionFloor,
        1.0 - CompositionFloor,
        1.0 - CompositionFloor).finished();
}

auto hpIGOPXTCMConstraintBridge() -> MAGEMinTCMConstraintBridge
{
    MAGEMinTCMConstraintBridge bridge;
    bridge.numConstraints = 12;
    bridge.numVariables = 8;
    bridge.constraintLowerBounds = ArrayXr::Constant(12, -1.0e12);
    bridge.constraintUpperBounds = ArrayXr::Zero(12);
    bridge.callback = [](unsigned m, double* result, unsigned n, const double* x, double* grad, void*)
    {
        if(m != 12 || n != 8)
            throw std::runtime_error("ig_opx TC mconstraint bridge expects m=12 and n=8.");

        // Direct transcription of MAGEMin TC callback `opx_ig_c`.
        result[0] = (-x[7]*x[3] - x[7]*x[0] + x[7] + x[3]*x[5] - x[3]*x[1] + x[3] + x[5]*x[0] - x[5] - x[0]*x[1] + x[0] + x[1] - 1.0);
        result[1] = (x[7]*x[3] + x[7]*x[0] - x[3]*x[5] + x[3]*x[1] - x[3] - x[5]*x[0] + x[0]*x[1] - x[0]);
        result[2] = (x[6] + x[4] - x[7] + 2.0*x[5] - x[1]);
        result[3] = (-x[4]);
        result[4] = (-x[6]);
        result[5] = (-x[5]);
        result[6] = (-x[2]*x[0] + x[2] + x[7]*x[3] - x[7]*x[0] + x[7] - x[3]*x[5] + x[3]*x[1] - x[3] + x[0] - 1.0);
        result[7] = (x[2]*x[0] - x[7]*x[3] + x[7]*x[0] + x[3]*x[5] - x[3]*x[1] + x[3] - x[0]);
        result[8] = (-x[2]);
        result[9] = (-x[7]);
        result[10] = (0.5*x[1] - 1.0);
        result[11] = (-0.5*x[1]);

        if(grad)
        {
            grad[0] = -x[7] + x[5] - x[1] + 1.0;
            grad[1] = -x[3] - x[0] + 1.0;
            grad[2] = 0.0;
            grad[3] = -x[7] + x[5] - x[1] + 1.0;
            grad[4] = 0.0;
            grad[5] = x[3] + x[0] - 1.0;
            grad[6] = 0.0;
            grad[7] = -x[3] - x[0] + 1.0;
            grad[8] = x[7] - x[5] + x[1] - 1.0;
            grad[9] = x[3] + x[0];
            grad[10] = 0.0;
            grad[11] = x[7] - x[5] + x[1] - 1.0;
            grad[12] = 0.0;
            grad[13] = -x[3] - x[0];
            grad[14] = 0.0;
            grad[15] = x[3] + x[0];
            grad[16] = 0.0;
            grad[17] = -1.0;
            grad[18] = 0.0;
            grad[19] = 0.0;
            grad[20] = 1.0;
            grad[21] = 2.0;
            grad[22] = 1.0;
            grad[23] = -1.0;
            grad[24] = 0.0;
            grad[25] = 0.0;
            grad[26] = 0.0;
            grad[27] = 0.0;
            grad[28] = -1.0;
            grad[29] = 0.0;
            grad[30] = 0.0;
            grad[31] = 0.0;
            grad[32] = 0.0;
            grad[33] = 0.0;
            grad[34] = 0.0;
            grad[35] = 0.0;
            grad[36] = 0.0;
            grad[37] = 0.0;
            grad[38] = -1.0;
            grad[39] = 0.0;
            grad[40] = 0.0;
            grad[41] = 0.0;
            grad[42] = 0.0;
            grad[43] = 0.0;
            grad[44] = 0.0;
            grad[45] = -1.0;
            grad[46] = 0.0;
            grad[47] = 0.0;
            grad[48] = -x[2] - x[7] + 1.0;
            grad[49] = x[3];
            grad[50] = 1.0 - x[0];
            grad[51] = x[7] - x[5] + x[1] - 1.0;
            grad[52] = 0.0;
            grad[53] = -x[3];
            grad[54] = 0.0;
            grad[55] = x[3] - x[0] + 1.0;
            grad[56] = x[2] + x[7] - 1.0;
            grad[57] = -x[3];
            grad[58] = x[0];
            grad[59] = -x[7] + x[5] - x[1] + 1.0;
            grad[60] = 0.0;
            grad[61] = x[3];
            grad[62] = 0.0;
            grad[63] = -x[3] + x[0];
            grad[64] = 0.0;
            grad[65] = 0.0;
            grad[66] = -1.0;
            grad[67] = 0.0;
            grad[68] = 0.0;
            grad[69] = 0.0;
            grad[70] = 0.0;
            grad[71] = 0.0;
            grad[72] = 0.0;
            grad[73] = 0.0;
            grad[74] = 0.0;
            grad[75] = 0.0;
            grad[76] = 0.0;
            grad[77] = 0.0;
            grad[78] = 0.0;
            grad[79] = -1.0;
            grad[80] = 0.0;
            grad[81] = 0.5;
            grad[82] = 0.0;
            grad[83] = 0.0;
            grad[84] = 0.0;
            grad[85] = 0.0;
            grad[86] = 0.0;
            grad[87] = 0.0;
            grad[88] = 0.0;
            grad[89] = -0.5;
            grad[90] = 0.0;
            grad[91] = 0.0;
            grad[92] = 0.0;
            grad[93] = 0.0;
            grad[94] = 0.0;
            grad[95] = 0.0;
        }
    };
    return bridge;
}

auto hpIGOPXEndmemberFractions(ArrayXrConstRef x) -> ArrayXr
{
    ArrayXr p(9);
    p[0] = -x[3]*x[7] + x[3]*x[5] - x[3]*x[1] + x[3] + x[2]*x[0] - x[2] + x[7]*x[0] - x[7] - x[0] - x[1] + 1.0;
    p[1] = -x[3]*x[7] + x[3]*x[5] - x[3]*x[1] + x[3] - x[7]*x[0] + x[5]*x[0] - x[0]*x[1] + x[0];
    p[2] = 2.0*x[3]*x[7] - 2.0*x[3]*x[5] + 2.0*x[3]*x[1] - 2.0*x[3] - x[2]*x[0] - x[5]*x[0] + x[0]*x[1];
    p[3] = x[2];
    p[4] = -x[6] - x[4] - 2.0*x[5] + x[1];
    p[5] = x[6];
    p[6] = 2.0*x[5];
    p[7] = x[4];
    p[8] = x[7];
    return p;
}

auto hpIGOPXSiteFractions(ArrayXrConstRef x) -> ArrayXr
{
    ArrayXr sf(12);
    sf[0] = x[0]*x[1] - x[0]*x[5] + x[0]*x[7] - x[0] + x[1]*x[3] - x[1] - x[3]*x[5] + x[3]*x[7] - x[3] + x[5] - x[7] + 1.0;
    sf[1] = -x[0]*x[1] + x[0]*x[5] - x[0]*x[7] + x[0] - x[1]*x[3] + x[3]*x[5] - x[3]*x[7] + x[3];
    sf[2] = x[1] - x[4] - 2.0*x[5] - x[6] + x[7];
    sf[3] = x[4];
    sf[4] = x[6];
    sf[5] = x[5];
    sf[6] = x[0]*x[2] + x[0]*x[7] - x[0] - x[1]*x[3] - x[2] + x[3]*x[5] - x[3]*x[7] + x[3] - x[7] + 1.0;
    sf[7] = -x[0]*x[2] - x[0]*x[7] + x[0] + x[1]*x[3] - x[3]*x[5] + x[3]*x[7] - x[3];
    sf[8] = x[2];
    sf[9] = x[7];
    sf[10] = 1.0 - 0.5*x[1];
    sf[11] = 0.5*x[1];
    return sf;
}

auto hpIGOPXEnThermo() -> const StandardThermoModel&
{
    static const auto model = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -2915480.248, -3090220.0, 132.5, 6.262e-05, 356.2, -0.00299, -596900.0, -3185.3,
        2.27e-05, 105900000000.0, 8.65, -8.2e-11, 10.0, 9999.0});
    return model;
}

auto hpIGOPXFsThermo() -> const StandardThermoModel&
{
    static const auto model = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -2234304.078, -2388710.0, 189.9, 6.592e-05, 398.7, -0.006579, 1290100.0, -4058.0,
        3.26e-05, 101000000000.0, 4.08, -4.0e-11, 10.0, 9999.0});
    return model;
}

auto hpIGOPXDiThermo() -> const StandardThermoModel&
{
    static const auto model = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -3027542.5655, -3201850.0, 142.9, 6.619e-05, 314.5, 4.1e-05, -2745900.0, -2020.1,
        2.73e-05, 119200000000.0, 5.19, -4.4e-11, 10.0, 9999.0});
    return model;
}

auto hpIGOPXMgtsThermo() -> const StandardThermoModel&
{
    static const auto model = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -3019929.6615, -3196670.0, 131.0, 6.05e-05, 371.4, -0.004082, -398400.0, -3547.1,
        2.17e-05, 102800000000.0, 8.55, -8.3e-11, 10.0, 9999.0});
    return model;
}

auto hpIGOPXCatsThermo() -> const StandardThermoModel&
{
    static const auto model = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -3131893.819, -3310110.0, 135.0, 6.356e-05, 347.6, -0.006974, -1781600.0, -2757.5,
        2.08e-05, 119200000000.0, 5.19, -4.4e-11, 10.0, 9999.0});
    return model;
}

auto hpIGOPXJdThermo() -> const StandardThermoModel&
{
    static const auto model = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -2846567.8345, -3025270.0, 133.5, 6.04e-05, 319.4, 0.003616, -1173900.0, -2469.5,
        2.1e-05, 128100000000.0, 3.81, -3.0e-11, 10.0, 9999.0});
    return model;
}

auto hpIGOPXKosThermo() -> const StandardThermoModel&
{
    static const auto model = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -2581405.5095, -2746840.0, 149.65, 6.309e-05, 309.2, 0.005419, -664600.0, -2176.6,
        1.94e-05, 130800000000.0, 3.0, -2.3e-11, 10.0, 9999.0});
    return model;
}

auto hpIGOPXCorThermo() -> const StandardThermoModel&
{
    static const auto model = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -1581770.16, -1675270.0, 50.9, 2.558e-05, 139.5, 0.00589, -2460600.0, -589.2,
        1.8e-05, 254000000000.0, 4.34, -1.7e-11, 5.0, 9999.0});
    return model;
}

auto hpIGOPXRuThermo() -> const StandardThermoModel&
{
    static const auto model = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -889103.8605, -944360.0, 50.5, 1.882e-05, 90.4, 0.0029, 0.0, -623.8,
        2.24e-05, 222000000000.0, 4.24, -1.9e-11, 3.0, 9999.0});
    return model;
}

auto hpIGOPXPerThermo() -> const StandardThermoModel&
{
    static const auto model = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -569097.243, -601530.0, 26.5, 1.125e-05, 60.5, 0.000362, -535800.0, -299.2,
        3.11e-05, 161600000000.0, 3.95, -2.4e-11, 2.0, 9999.0});
    return model;
}

auto hpIGOPXAcmThermo() -> const StandardThermoModel&
{
    static const auto model = StandardThermoModelHollandPowell(StandardThermoModelParamsHollandPowell{
        -2416108.22, -2583430.0, 170.6, 6.459e-05, 307.1, 0.016758, -1685500.0, -2125.8,
        2.11e-05, 106000000000.0, 4.08, -3.8e-11, 10.0, 9999.0});
    return model;
}

auto hpIGOPXReferenceState(real T, real P) -> ArrayXr
{
    const auto en = hpIGOPXEnThermo()(T, P).G0;
    const auto fs = hpIGOPXFsThermo()(T, P).G0;
    const auto di = hpIGOPXDiThermo()(T, P).G0;
    const auto mgts = hpIGOPXMgtsThermo()(T, P).G0;
    const auto cats = hpIGOPXCatsThermo()(T, P).G0;
    const auto jd = hpIGOPXJdThermo()(T, P).G0;
    const auto kos = hpIGOPXKosThermo()(T, P).G0;
    const auto cor = hpIGOPXCorThermo()(T, P).G0;
    const auto per = hpIGOPXPerThermo()(T, P).G0;
    const auto ru = hpIGOPXRuThermo()(T, P).G0;
    const auto acm = hpIGOPXAcmThermo()(T, P).G0;

    const auto Pbar = P / 1.0e5;

    ArrayXr gb(9);
    gb[0] = en;
    gb[1] = fs;
    gb[2] = 0.5*en + 0.5*fs - 6.6;
    gb[3] = 0.005*Pbar + di + 2.8;
    gb[4] = mgts;
    gb[5] = 0.14*Pbar - di + cats + en - jd + kos - 6.0;
    gb[6] = 0.37*Pbar - 0.0051*T - 0.5*cor + mgts + 0.5*per + 0.5*ru - 3.91;
    gb[7] = acm - jd + mgts + 3.0;
    gb[8] = jd + 18.2;
    return gb;
}

auto hpIGOPXObjective(real T, real P, ArrayXrConstRef visiblex, ArrayXrConstRef x, real externalCompositionPenalty, Fn<ArrayXr(real, real)> const& referenceState) -> real
{
    constexpr std::array<double, 36> W = {
        7.00, 3.50, 29.0, 12.5, 8.00, 6.00, 8.00, 35.0, 4.50,
        23.0, 11.0, 10.0, 7.00, 10.0, 35.0, 19.0, 15.0, 12.0,
        8.00, 12.0, 35.0, 75.5, 20.0, 40.0, 20.0, 35.0, 2.00,
        10.0, 2.00, 7.00, 6.00, 2.00, -11.0, 6.00, 20.0, -11.0};
    constexpr std::array<double, 9> v = {1.00, 1.00, 1.00, 1.20, 1.00, 1.00, 1.00, 1.00, 1.20};

    const auto p = hpIGOPXEndmemberFractions(x);
    const auto sf = hpIGOPXSiteFractions(x);
    const auto gb = referenceState(T, P);

    const auto sumv = static_cast<double>(
        p[0]*v[0] + p[1]*v[1] + p[2]*v[2] + p[3]*v[3] + p[4]*v[4] +
        p[5]*v[5] + p[6]*v[6] + p[7]*v[7] + p[8]*v[8]);
    if(!(sumv > 0.0))
        return std::numeric_limits<double>::infinity();

    ArrayXr phi(9);
    for(Index i = 0; i < 9; ++i)
        phi[i] = p[i] * v[static_cast<std::size_t>(i)] / sumv;

    ArrayXr muGex(9);
    for(Index i = 0; i < 9; ++i)
    {
        real Gex = 0.0;
        Index it = 0;
        for(Index j = 0; j < 8; ++j)
        {
            const auto eyeij = (i == j) ? 1.0 : 0.0;
            const auto tmp = eyeij - static_cast<double>(phi[j]);
            for(Index k = j + 1; k < 9; ++k)
            {
                const auto eyeik = (k < 8 && i == k) ? 1.0 : 0.0;
                Gex -= tmp*(eyeik - static_cast<double>(phi[k]))*
                    (W[static_cast<std::size_t>(it)] * 2.0 * v[static_cast<std::size_t>(i)] /
                        (v[static_cast<std::size_t>(j)] + v[static_cast<std::size_t>(k)]));
                ++it;
            }
        }
        muGex[i] = Gex;
    }

    const auto safe = [](real value) -> real { return std::max(value, static_cast<real>(CompositionFloor)); };
    const auto ln2 = static_cast<real>(std::log(2.0));
    const auto lnsqrt2 = static_cast<real>(0.5 * std::log(2.0));

    ArrayXr mu(9);
    mu[0] = gb[0] + universalGasConstant*T*(log(safe(sf[0])) + 0.5*log(safe(sf[10])) + log(safe(sf[6]))) + muGex[0];
    mu[1] = gb[1] + universalGasConstant*T*(0.5*log(safe(sf[10])) + log(safe(sf[1])) + log(safe(sf[7]))) + muGex[1];
    mu[2] = gb[2] + universalGasConstant*T*(log(safe(sf[0])) + 0.5*log(safe(sf[10])) + log(safe(sf[7]))) + muGex[2];
    mu[3] = gb[3] + universalGasConstant*T*(log(safe(sf[0])) + 0.5*log(safe(sf[10])) + log(safe(sf[8]))) + muGex[3];
    mu[4] = gb[4] + universalGasConstant*T*(lnsqrt2 + 0.25*log(safe(sf[10])) + 0.25*log(safe(sf[11])) + log(safe(sf[2])) + log(safe(sf[6]))) + muGex[4];
    mu[5] = gb[5] + universalGasConstant*T*(lnsqrt2 + 0.25*log(safe(sf[10])) + 0.25*log(safe(sf[11])) + log(safe(sf[4])) + log(safe(sf[6]))) + muGex[5];
    mu[6] = gb[6] + universalGasConstant*T*(ln2 + lnsqrt2 + 0.5*log(safe(sf[0])) + 0.25*log(safe(sf[10])) + 0.25*log(safe(sf[11])) + 0.5*log(safe(sf[5])) + log(safe(sf[6]))) + muGex[6];
    mu[7] = gb[7] + universalGasConstant*T*(lnsqrt2 + 0.25*log(safe(sf[10])) + 0.25*log(safe(sf[11])) + log(safe(sf[3])) + log(safe(sf[6]))) + muGex[7];
    mu[8] = gb[8] + universalGasConstant*T*(0.5*log(safe(sf[10])) + log(safe(sf[2])) + log(safe(sf[9]))) + muGex[8];

    real objective = p.matrix().dot(mu.matrix());
    objective += externalCompositionPenalty*universalGasConstant*T*(x - visiblex).matrix().squaredNorm();
    return objective;
}

auto hpIGOPXObjectiveGradientFiniteDifference(
    real T,
    real P,
    ArrayXrConstRef visiblex,
    ArrayXrConstRef x,
    ArrayXrConstRef lowerBounds,
    ArrayXrConstRef upperBounds,
    real externalCompositionPenalty,
    Fn<ArrayXr(real, real)> const& referenceState) -> ArrayXr
{
    ArrayXr gradient(x.size());
    for(Index i = 0; i < x.size(); ++i)
    {
        const auto xi = static_cast<double>(x[i]);
        const auto h = std::max(1.0e-8, 1.0e-7*std::max(1.0, std::abs(xi)));

        ArrayXr xp = x;
        ArrayXr xm = x;
        xp[i] = std::clamp(xi + h, static_cast<double>(lowerBounds[i]), static_cast<double>(upperBounds[i]));
        xm[i] = std::clamp(xi - h, static_cast<double>(lowerBounds[i]), static_cast<double>(upperBounds[i]));

        const auto dx = static_cast<double>(xp[i] - xm[i]);
        if(dx <= 0.0)
        {
            gradient[i] = 0.0;
            continue;
        }

        const auto fp = hpIGOPXObjective(T, P, visiblex, xp, externalCompositionPenalty, referenceState);
        const auto fm = hpIGOPXObjective(T, P, visiblex, xm, externalCompositionPenalty, referenceState);
        gradient[i] = (fp - fm) / dx;
    }
    return gradient;
}

auto sb21GTMJThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb21_gtmj";
    thermo.endmembers = {"py", "alm", "gr", "mgmj", "jdmj"};
    // Symmetric regular Margules — pairs (py-alm, py-gr, py-mgmj, py-jdmj,
    //   alm-gr, alm-mgmj, alm-jdmj, gr-mgmj, gr-jdmj, mgmj-jdmj).
    // W[1]=W[py-gr] and W[7]=W[gr-mgmj] have a small P-dependent correction
    // (+1.03e-5 * P_bar) that is neglected here (pilot-level approximation).
    thermo.excessChemicalPotentials = [](ArrayXrConstRef y)
    {
        constexpr Index n = 5;
        constexpr double W[] = {0.0, 21117.58, 22672.42, 22672.42, 21117.58, 22672.42, 22672.42, 60718.2, 60718.2, 70879.14};
        ArrayXr mu(n);
        for(Index i = 0; i < n; ++i)
        {
            real Gex = 0.0;
            Index it = 0;
            for(Index j = 0; j < n; ++j)
            {
                const auto tmp = static_cast<double>(i == j ? 1.0 : 0.0) - static_cast<double>(y[j]);
                for(Index k = j + 1; k < n; ++k)
                {
                    const auto delta = static_cast<double>(i == k ? 1.0 : 0.0) - static_cast<double>(y[k]);
                    Gex -= tmp * delta * W[it++];
                }
            }
            mu[i] = Gex;
        }
        return mu;
    };
    // Site-mixing entropy from MAGEMin Sconfig for sb21_gtmj
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        const auto s1 = y[0]+y[1]+y[2];
        const auto s2 = y[0]+y[1]+y[2]+y[4];
        const auto s3 = y[0]+y[3];
        const auto s4 = y[3]+y[4];
        return universalGasConstant * T * (
            s1*log(s1) + s2*log(s2) + 3.0*s3*log(s3)
            + 3.0*y[1]*log(y[1]) + 3.0*y[2]*log(y[2])
            + y[3]*log(y[3]) + s4*log(s4)
            + 2.0*y[4]*log(2.0/3.0*y[4]) + y[4]*log(1.0/3.0*y[4]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        const auto s1 = y[0]+y[1]+y[2];
        const auto s2 = y[0]+y[1]+y[2]+y[4];
        const auto s3 = y[0]+y[3];
        const auto s4 = y[3]+y[4];
        ArrayXr ln_a(5);
        ln_a[0] = log(s1) + log(s2) + 3.0*log(s3);
        ln_a[1] = log(s1) + log(s2) + 3.0*log(y[1]);
        ln_a[2] = log(s1) + log(s2) + 3.0*log(y[2]);
        ln_a[3] = 3.0*log(s3) + log(y[3]) + log(s4);
        ln_a[4] = 2.0*log(2.0/3.0*y[4]) + log(1.0/3.0*y[4]) + log(s2) + log(s4);
        return ln_a;
    };
    return thermo;
}

auto sb21CalcioferriteThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb21_cf";
    thermo.endmembers = {"mgcf", "fecf", "nacf"};
    thermo.W01 = 0.0;
    thermo.W02 = 60825.08;
    thermo.W12 = 60825.08;
    thermo.excessChemicalPotentials = [](ArrayXrConstRef y)
    {
        const ArrayXr volumes = (ArrayXr(3) << 1.0, 1.0, 4.4532).finished();
        const auto sumv = y.matrix().dot(volumes.matrix());

        ArrayXr phi(3);
        phi[0] = y[0]*volumes[0]/sumv;
        phi[1] = y[1]*volumes[1]/sumv;
        phi[2] = y[2]*volumes[2]/sumv;

        const ArrayXr interactions = (ArrayXr(3) << 0.0, 60825.08, 60825.08).finished();

        ArrayXr mu(3);
        for(Index i = 0; i < 3; ++i)
        {
            real Gex = 0.0;
            Index interaction = 0;
            for(Index j = 0; j < 3; ++j)
            {
                const auto tmp = ((i == j) ? 1.0 : 0.0) - static_cast<double>(phi[j]);
                for(Index k = j + 1; k < 3; ++k)
                {
                    const auto delta = ((i == k) ? 1.0 : 0.0) - static_cast<double>(phi[k]);
                    Gex -= tmp*delta*(interactions[interaction]*2.0*volumes[i]/(volumes[j] + volumes[k]));
                    ++interaction;
                }
            }
            mu[i] = Gex;
        }

        return mu;
    };
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        return universalGasConstant * T * (
            y[0]*log(y[0])
            + (y[0] + y[1])*log(y[0] + y[1])
            + y[1]*log(y[1])
            + 2.0*y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        ArrayXr ln_a(3);
        ln_a[0] = log(y[0]) + log(y[0] + y[1]);
        ln_a[1] = log(y[0] + y[1]) + log(y[1]);
        ln_a[2] = 2.0*log(y[2]);
        return ln_a;
    };
    return thermo;
}

auto sb11AkimotoiteThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb11_ak";
    thermo.endmembers = {"mgak", "feak", "co"};
    thermo.W01 = 0.0;
    thermo.W02 = 66000.0;
    thermo.W12 = 0.0;
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        return universalGasConstant * T * (
            2.0*y[0]*log(y[0])
            + (y[1] + y[2])*log(y[1] + y[2])
            + y[1]*log(y[1])
            + y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        ArrayXr ln_a(3);
        ln_a[0] = 2.0*log(y[0]);
        ln_a[1] = log(y[1] + y[2]) + log(y[1]);
        ln_a[2] = log(y[1] + y[2]) + log(y[2]);
        return ln_a;
    };
    return thermo;
}

auto sb11PerovskiteThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb11_pv";
    thermo.endmembers = {"mgpv", "fepv", "alpv"};
    thermo.W01 = 0.0;
    thermo.W02 = 116000.0;
    thermo.W12 = 0.0;
    thermo.excessChemicalPotentials = [](ArrayXrConstRef y)
    {
        const auto y0 = static_cast<double>(y[0]);
        const auto y1 = static_cast<double>(y[1]);
        const auto y2 = static_cast<double>(y[2]);
        const auto v0 = 1.0;
        const auto v1 = 1.0;
        const auto v2 = 0.39;
        const auto sumv = y0*v0 + y1*v1 + y2*v2;

        ArrayXr phi(3);
        phi[0] = y0*v0/sumv;
        phi[1] = y1*v1/sumv;
        phi[2] = y2*v2/sumv;

        ArrayXr mu(3);
        mu[0] = 2.0*116000.0*v0*phi[2]*(1.0 - phi[0])/(v0 + v2);
        mu[1] = -2.0*116000.0*v1*phi[0]*phi[2]/(v0 + v2);
        mu[2] = 2.0*116000.0*v2*phi[0]*(1.0 - phi[2])/(v0 + v2);
        return mu;
    };
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        return universalGasConstant * T * (
            2.0*y[0]*log(y[0])
            + (y[1] + y[2])*log(y[1] + y[2])
            + y[1]*log(y[1])
            + y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        ArrayXr ln_a(3);
        ln_a[0] = 2.0*log(y[0]);
        ln_a[1] = log(y[1] + y[2]) + log(y[1]);
        ln_a[2] = log(y[1] + y[2]) + log(y[2]);
        return ln_a;
    };
    return thermo;
}

auto sb11CalcioferriteThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb11_cf";
    thermo.endmembers = {"mgcf", "fecf", "nacf"};
    thermo.W01 = 0.0;
    thermo.W02 = 0.0;
    thermo.W12 = 0.0;
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        return universalGasConstant * T * (
            y[0]*log(y[0])
            + (y[0] + y[2])*log(y[0] + y[2])
            + 2.0*y[1]*log(y[1])
            + y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        ArrayXr ln_a(3);
        ln_a[0] = log(y[0]) + log(y[0] + y[2]);
        ln_a[1] = 2.0*log(y[1]);
        ln_a[2] = log(y[2]) + log(y[0] + y[2]);
        return ln_a;
    };
    return thermo;
}

auto normalizedMAGEMinPilotBranches(Vec<GlobalizedSolidSolutionBranch> branches, Index numCoords) -> Vec<GlobalizedSolidSolutionBranch>
{
    return NormalizeGlobalizedSolidSolutionBranches(
        std::move(branches),
        numCoords,
        "MAGEMin pilot branches must define bounds matching the visible composition size.");
}

// ---------------------------------------------------------------------------
// SB21 thermo helper functions (binary phases)
// ---------------------------------------------------------------------------

auto sb21PLGThermo() -> MAGEMinImportedBinarySolutionThermoModel
{
    MAGEMinImportedBinarySolutionThermoModel thermo;
    thermo.modelId = "sb21_plg";
    thermo.endmember0 = "an";
    thermo.endmember1 = "ab";
    thermo.W = 13000.0;
    thermo.idealSiteMultiplicity = 1.0;
    return thermo;
}

auto sb21OlivineThermo() -> MAGEMinImportedBinarySolutionThermoModel
{
    MAGEMinImportedBinarySolutionThermoModel thermo;
    thermo.modelId = "sb21_ol";
    thermo.endmember0 = "fo";
    thermo.endmember1 = "fa";
    thermo.W = 4694.66;
    thermo.idealSiteMultiplicity = 2.0;
    return thermo;
}

auto sb21WadsleyiteThermo() -> MAGEMinImportedBinarySolutionThermoModel
{
    MAGEMinImportedBinarySolutionThermoModel thermo;
    thermo.modelId = "sb21_wa";
    thermo.endmember0 = "mgwa";
    thermo.endmember1 = "fewa";
    thermo.W = 13202.38;
    thermo.idealSiteMultiplicity = 2.0;
    return thermo;
}

auto sb21RingwooditeThermo() -> MAGEMinImportedBinarySolutionThermoModel
{
    MAGEMinImportedBinarySolutionThermoModel thermo;
    thermo.modelId = "sb21_ri";
    thermo.endmember0 = "mgri";
    thermo.endmember1 = "feri";
    thermo.W = 7600.74;
    thermo.idealSiteMultiplicity = 2.0;
    return thermo;
}

auto sb21HPCPXThermo() -> MAGEMinImportedBinarySolutionThermoModel
{
    MAGEMinImportedBinarySolutionThermoModel thermo;
    thermo.modelId = "sb21_hpcpx";
    thermo.endmember0 = "hpcen";
    thermo.endmember1 = "hpcfs";
    thermo.W = 0.0;
    thermo.idealSiteMultiplicity = 2.0;
    return thermo;
}

// ---------------------------------------------------------------------------
// SB21 thermo helper functions (ternary phases)
// ---------------------------------------------------------------------------

auto sb21AkimotoiteThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb21_ak";
    thermo.endmembers = {"mgak", "feak", "co"};
    // Symmetric regular Margules — pairs (mgak-feak, mgak-co, feak-co)
    thermo.excessChemicalPotentials = [](ArrayXrConstRef y)
    {
        constexpr Index n = 3;
        constexpr double W[] = {0.0, 59348.69, 59348.69};
        ArrayXr mu(n);
        for(Index i = 0; i < n; ++i)
        {
            real Gex = 0.0;
            Index it = 0;
            for(Index j = 0; j < n; ++j)
            {
                const auto tmp = static_cast<double>(i == j ? 1.0 : 0.0) - static_cast<double>(y[j]);
                for(Index k = j + 1; k < n; ++k)
                {
                    const auto delta = static_cast<double>(i == k ? 1.0 : 0.0) - static_cast<double>(y[k]);
                    Gex -= tmp * delta * W[it++];
                }
            }
            mu[i] = Gex;
        }
        return mu;
    };
    // Site-mixing entropy: Sconfig/RT = p0*ln(p0) + (p0+p1)*ln(p0+p1) + p1*ln(p1) + 2*p2*ln(p2)
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        const auto s01 = y[0] + y[1];
        return universalGasConstant * T * (
            y[0]*log(y[0]) + s01*log(s01) + y[1]*log(y[1]) + 2.0*y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        const auto s01 = y[0] + y[1];
        ArrayXr ln_a(3);
        ln_a[0] = log(y[0]) + log(s01);
        ln_a[1] = log(s01) + log(y[1]);
        ln_a[2] = 2.0*log(y[2]);
        return ln_a;
    };
    return thermo;
}

auto sb21PerovskiteThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb21_pv";
    thermo.endmembers = {"mgpv", "fepv", "alpv"};
    // Symmetric regular Margules — pairs (mgpv-fepv, mgpv-alpv, fepv-alpv)
    thermo.excessChemicalPotentials = [](ArrayXrConstRef y)
    {
        constexpr Index n = 3;
        constexpr double W[] = {-11396.17, 34979.87, 0.0};
        ArrayXr mu(n);
        for(Index i = 0; i < n; ++i)
        {
            real Gex = 0.0;
            Index it = 0;
            for(Index j = 0; j < n; ++j)
            {
                const auto tmp = static_cast<double>(i == j ? 1.0 : 0.0) - static_cast<double>(y[j]);
                for(Index k = j + 1; k < n; ++k)
                {
                    const auto delta = static_cast<double>(i == k ? 1.0 : 0.0) - static_cast<double>(y[k]);
                    Gex -= tmp * delta * W[it++];
                }
            }
            mu[i] = Gex;
        }
        return mu;
    };
    // Sconfig/RT = p0*ln(p0) + (p0+p1)*ln(p0+p1) + p1*ln(p1) + 2*p2*ln(p2)
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        const auto s01 = y[0] + y[1];
        return universalGasConstant * T * (
            y[0]*log(y[0]) + s01*log(s01) + y[1]*log(y[1]) + 2.0*y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        const auto s01 = y[0] + y[1];
        ArrayXr ln_a(3);
        ln_a[0] = log(y[0]) + log(s01);
        ln_a[1] = log(s01) + log(y[1]);
        ln_a[2] = 2.0*log(y[2]);
        return ln_a;
    };
    return thermo;
}

auto sb21PostPerovskiteThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb21_ppv";
    thermo.endmembers = {"mppv", "fppv", "appv"};
    // Symmetric regular Margules — pairs (mppv-fppv, mppv-appv, fppv-appv)
    thermo.excessChemicalPotentials = [](ArrayXrConstRef y)
    {
        constexpr Index n = 3;
        constexpr double W[] = {-10955.49, 34979.87, 34979.87};
        ArrayXr mu(n);
        for(Index i = 0; i < n; ++i)
        {
            real Gex = 0.0;
            Index it = 0;
            for(Index j = 0; j < n; ++j)
            {
                const auto tmp = static_cast<double>(i == j ? 1.0 : 0.0) - static_cast<double>(y[j]);
                for(Index k = j + 1; k < n; ++k)
                {
                    const auto delta = static_cast<double>(i == k ? 1.0 : 0.0) - static_cast<double>(y[k]);
                    Gex -= tmp * delta * W[it++];
                }
            }
            mu[i] = Gex;
        }
        return mu;
    };
    // Sconfig/RT = p0*ln(p0) + (p0+p1)*ln(p0+p1) + p1*ln(p1) + 2*p2*ln(p2)
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        const auto s01 = y[0] + y[1];
        return universalGasConstant * T * (
            y[0]*log(y[0]) + s01*log(s01) + y[1]*log(y[1]) + 2.0*y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        const auto s01 = y[0] + y[1];
        ArrayXr ln_a(3);
        ln_a[0] = log(y[0]) + log(s01);
        ln_a[1] = log(s01) + log(y[1]);
        ln_a[2] = 2.0*log(y[2]);
        return ln_a;
    };
    return thermo;
}

auto sb21MagnesiowustitesThermo() -> MAGEMinImportedConstrainedTernarySolutionThermoModel
{
    MAGEMinImportedConstrainedTernarySolutionThermoModel thermo;
    thermo.modelId = "sb21_mw";
    thermo.endmembers = {"pe", "wu", "anao"};
    // Symmetric regular Margules — pairs (pe-wu, pe-anao, wu-anao)
    // W[0] = 44000 + 4.4e-6*P in MAGEMin (P in MPa); P-correction negligible at pilot level.
    thermo.excessChemicalPotentials = [](ArrayXrConstRef y)
    {
        constexpr Index n = 3;
        constexpr double W[] = {44000.0, 120000.0, 120000.0};
        ArrayXr mu(n);
        for(Index i = 0; i < n; ++i)
        {
            real Gex = 0.0;
            Index it = 0;
            for(Index j = 0; j < n; ++j)
            {
                const auto tmp = static_cast<double>(i == j ? 1.0 : 0.0) - static_cast<double>(y[j]);
                for(Index k = j + 1; k < n; ++k)
                {
                    const auto delta = static_cast<double>(i == k ? 1.0 : 0.0) - static_cast<double>(y[k]);
                    Gex -= tmp * delta * W[it++];
                }
            }
            mu[i] = Gex;
        }
        return mu;
    };
    // Sconfig/RT = 4*(p0*ln(p0) + p1*ln(p1) + p2*ln(p2))
    thermo.idealGibbs = [](real T, ArrayXrConstRef y)
    {
        return universalGasConstant * T * 4.0 * (
            y[0]*log(y[0]) + y[1]*log(y[1]) + y[2]*log(y[2]));
    };
    thermo.idealLnActivities = [](ArrayXrConstRef y)
    {
        ArrayXr ln_a(3);
        ln_a[0] = 4.0*log(y[0]);
        ln_a[1] = 4.0*log(y[1]);
        ln_a[2] = 4.0*log(y[2]);
        return ln_a;
    };
    return thermo;
}



auto projectSeedToBranch(ArrayXr seed, GlobalizedSolidSolutionBranch const& branch) -> ArrayXr
{
    if(seed.size() != branch.lowerBounds.size() || seed.size() != branch.upperBounds.size())
        return seed;

    for(Index i = 0; i < seed.size(); ++i)
        seed[i] = std::clamp(static_cast<double>(seed[i]), static_cast<double>(branch.lowerBounds[i]), static_cast<double>(branch.upperBounds[i]));

    for(Index iter = 0; iter < seed.size() * 8; ++iter)
    {
        const auto residual = 1.0 - static_cast<double>(seed.sum());
        if(std::abs(residual) <= CandidateSeedTolerance)
            break;

        real capacity = 0.0;
        for(Index i = 0; i < seed.size(); ++i)
        {
            const auto slack = residual > 0.0 ? branch.upperBounds[i] - seed[i] : seed[i] - branch.lowerBounds[i];
            if(slack > CandidateSeedTolerance)
                capacity += slack;
        }

        if(capacity <= CandidateSeedTolerance)
            break;

        for(Index i = 0; i < seed.size(); ++i)
        {
            const auto slack = residual > 0.0 ? branch.upperBounds[i] - seed[i] : seed[i] - branch.lowerBounds[i];
            if(slack <= CandidateSeedTolerance)
                continue;

            const auto delta = residual * (slack / capacity);
            seed[i] = std::clamp(
                static_cast<double>(seed[i] + delta),
                static_cast<double>(branch.lowerBounds[i]),
                static_cast<double>(branch.upperBounds[i]));
        }
    }

    seed = seed.max(CompositionFloor);
    seed /= seed.sum();
    return seed;
}

auto seedsEquivalent(ArrayXrConstRef lhs, ArrayXrConstRef rhs, real tolerance) -> bool
{
    return lhs.size() == rhs.size() && (lhs - rhs).cwiseAbs().maxCoeff() <= tolerance;
}

auto appendDistinctCandidate(
    Vec<GlobalizedSolidSolutionCandidate>& candidates,
    GlobalizedSolidSolutionCandidate candidate,
    String const& seedLabel,
    real extraPriority = 0.0) -> void
{
    if(candidate.initialInternalx.size() == 0)
        return;

    for(const auto& existing : candidates)
    {
        if(existing.branch != candidate.branch)
            continue;
        if(existing.initialInternalx.size() == candidate.initialInternalx.size()
            && seedsEquivalent(existing.initialInternalx, candidate.initialInternalx, CandidateSeedTolerance))
            return;
    }

    candidate.priority += extraPriority;
    candidate.extra["MAGEMinSolidSolutionPilot::CandidateSeedLabel"] = seedLabel;
    candidates.push_back(std::move(candidate));
}

struct MAGEMinTernaryCandidateSeedSpec
{
    ArrayXr seed;
    String label;
    real priority = 0.0;
};

auto dominantEndmemberSeed(Index size, Index dominant) -> ArrayXr
{
    ArrayXr seed = ArrayXr::Constant(size, CompositionFloor);
    seed[dominant] = 1.0 - (size - 1)*CompositionFloor;
    return seed;
}

auto binaryEdgeMidpointSeed(Index size, Index lhs, Index rhs) -> ArrayXr
{
    ArrayXr seed = ArrayXr::Constant(size, CompositionFloor);
    const auto retained = 1.0 - (size - 2)*CompositionFloor;
    seed[lhs] = 0.5*retained;
    seed[rhs] = 0.5*retained;
    return seed;
}

auto normalizedDominantEndmemberOrder(
    MAGEMinStructuredTernaryProposalOptions const& options,
    Index size) -> Indices
{
    Indices order;
    order.reserve(size);

    for(const auto index : options.dominantEndmemberOrder)
    {
        if(index >= size)
            continue;
        if(std::find(order.begin(), order.end(), index) == order.end())
            order.push_back(index);
    }

    for(Index i = 0; i < size; ++i)
    {
        if(std::find(order.begin(), order.end(), i) == order.end())
            order.push_back(i);
    }

    return order;
}

auto defaultTernaryCandidateSeedSpecs(
    MAGEMinImportedConstrainedTernarySolutionThermoModel const& thermo,
    MAGEMinStructuredTernaryProposalOptions const& options,
    ArrayXrConstRef visiblex) -> Vec<MAGEMinTernaryCandidateSeedSpec>
{
    Vec<MAGEMinTernaryCandidateSeedSpec> specs;

    if(options.includeVisibleCompositionSeed)
        specs.push_back({ArrayXr(visiblex), "visible-composition", options.visibleCompositionPriority});

    const auto dominantOrder = normalizedDominantEndmemberOrder(options, visiblex.size());

    auto appendSeed = [&](ArrayXr seed, String label, real priority)
    {
        specs.push_back({std::move(seed), std::move(label), priority});
    };

    if(options.includeDominantEndmemberSeeds)
    {
        for(Index i = 0; i < dominantOrder.size(); ++i)
        {
            const auto dominant = dominantOrder[i];
            appendSeed(
                dominantEndmemberSeed(visiblex.size(), dominant),
                "dominant::" + thermo.endmembers[dominant],
                options.dominantEndmemberPriority + i*options.dominantEndmemberPriorityStep);
        }
    }

    if(options.includeBinaryEdgeMidpointSeeds)
    {
        Index pairRank = 0;
        for(Index i = 0; i < dominantOrder.size(); ++i)
        {
            for(Index j = i + 1; j < dominantOrder.size(); ++j)
            {
                const auto lhs = dominantOrder[i];
                const auto rhs = dominantOrder[j];
                appendSeed(
                    binaryEdgeMidpointSeed(visiblex.size(), lhs, rhs),
                    "edge::" + thermo.endmembers[lhs] + "-" + thermo.endmembers[rhs],
                    options.binaryEdgePriority - pairRank*1.0e-7);
                pairRank += 1;
            }
        }
    }

    return specs;
}

auto augmentDefaultConstrainedTernaryCandidates(
    MAGEMinImportedConstrainedTernarySolutionThermoModel const& thermo,
    MAGEMinStructuredTernaryProposalOptions const& proposalOptions,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    ArrayXrConstRef visiblex,
    Vec<GlobalizedSolidSolutionCandidate> candidates) -> Vec<GlobalizedSolidSolutionCandidate>
{
    if(candidates.empty())
        return candidates;

    const auto sourceIt = candidates[0].extra.find("MAGEMinSolidSolutionPilot::CandidateSource");
    const auto source = sourceIt != candidates[0].extra.end()
        ? std::any_cast<String>(&sourceIt->second)
        : nullptr;

    if(!(source && (*source == "branch-screen" || *source == "requested-branch")))
        return candidates;

    const auto seedSpecs = defaultTernaryCandidateSeedSpecs(thermo, proposalOptions, visiblex);
    Strings generatedSeedLabels;
    generatedSeedLabels.reserve(seedSpecs.size());
    for(const auto& spec : seedSpecs)
        generatedSeedLabels.push_back(spec.label);

    Vec<GlobalizedSolidSolutionCandidate> augmented = candidates;
    for(const auto& candidate : candidates)
    {
        const auto& branch = branches[candidate.branch];
        for(const auto& spec : seedSpecs)
        {
            auto seeded = candidate;
            seeded.initialInternalx = projectSeedToBranch(spec.seed, branch);
            appendDistinctCandidate(augmented, std::move(seeded), spec.label, spec.priority);
        }
    }

    for(auto& candidate : augmented)
    {
        candidate.extra["MAGEMinSolidSolutionPilot::GeneratedCandidateCount"] = static_cast<std::uint64_t>(augmented.size());
        candidate.extra["MAGEMinSolidSolutionPilot::GeneratedCandidateSeedLabels"] = generatedSeedLabels;
    }

    return augmented;
}

auto defaultCandidates(
    MAGEMinSolidSolutionPilotBranchPolicyOptions const& options,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionInput const& input,
    ArrayXrConstRef seedx) -> Vec<GlobalizedSolidSolutionCandidate>
{
    Optional<ArrayXr> cachedWarmstart = std::nullopt;
    Optional<ArrayXr> branchWarmstart = std::nullopt;
    if(input.state && input.state->cachedInternalx.size() == seedx.size())
        cachedWarmstart = input.state->cachedInternalx;
    if(input.state && input.state->lastInternalx.size() == seedx.size())
        branchWarmstart = input.state->lastInternalx;

    GlobalizedSolidSolutionDefaultCandidateOptions candidateOptions;
    candidateOptions.branchTolerance = options.branchTolerance;
    candidateOptions.cachedStatePriority = -1.0;
    candidateOptions.preferredBranchPriority = -options.branchScoreHysteresis;
    candidateOptions.requireCachedStateWarmstart = true;
    candidateOptions.sourceKey = "MAGEMinSolidSolutionPilot::CandidateSource";
    candidateOptions.invalidRequestedBranchMessage = "Requested MAGEMin imported solid-solution pilot branch is out of range.";

    // During duplicated-phase requested-branch evaluations, disable cached-state reuse so
    // requested branch selection is not short-circuited by stale cached branch state.
    // Keep branch warmstart enabled to preserve solver stability.
    if(input.requestedBranch != GlobalizedSolidSolutionNoBranch)
    {
        cachedWarmstart = std::nullopt;
    }

    ArrayXr initialSeed = cachedWarmstart.has_value()
        ? ArrayXr(*cachedWarmstart)
        : ArrayXr(seedx);

    auto candidates = DefaultGlobalizedSolidSolutionCandidates(
        input,
        branches,
        initialSeed,
        cachedWarmstart,
        branchWarmstart,
        candidateOptions);

    for(auto& candidate : candidates)
        candidate.extra["MAGEMinSolidSolutionPilot::GeneratedCandidateCount"] = static_cast<std::uint64_t>(candidates.size());

    return candidates;
}

auto branchBoundsEqual(ArrayXrConstRef lhs, ArrayXrConstRef rhs, real tolerance) -> bool
{
    if(lhs.size() != rhs.size())
        return false;
    if(lhs.size() == 0)
        return true;
    return (lhs - rhs).cwiseAbs().maxCoeff() <= tolerance;
}

auto findBranchIndex(
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionBranch const& branch,
    real tolerance) -> Index
{
    for(Index i = 0; i < branches.size(); ++i)
    {
        if(!branch.id.empty() && branches[i].id == branch.id)
            return i;

        if(branchBoundsEqual(branches[i].lowerBounds, branch.lowerBounds, tolerance)
            && branchBoundsEqual(branches[i].upperBounds, branch.upperBounds, tolerance))
            return i;
    }

    return GlobalizedSolidSolutionNoBranch;
}

auto defaultPilotStabilityCriterion(
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    real branchTolerance,
    String const& splitViolationKey) -> GlobalizedSolidSolutionCandidateStabilityCriterion
{
    GlobalizedSolidSolutionBranchAmbiguityStabilityOptions options;
    options.branchTolerance = branchTolerance;
    options.violationKey = splitViolationKey;
    return NamedGlobalizedSolidSolutionStabilityCriterion(
        NamedGlobalizedSolidSolutionStabilityPolicy::MAGEMinPilotBranchAmbiguity,
        branches,
        options);
}

struct PrecomputedConstrainedTernaryCandidateEvaluation
{
    GlobalizedSolidSolutionBranchSelection selection;
};

auto branchCandidateLabel(GlobalizedSolidSolutionBranch const& branch, Index branchIndex) -> String
{
    if(!branch.label.empty())
        return branch.label;
    if(!branch.id.empty())
        return branch.id;
    return std::to_string(branchIndex);
}

auto candidateReusedState(GlobalizedSolidSolutionCandidate const& candidate) -> bool
{
    const auto it = candidate.extra.find("MAGEMinSolidSolutionPilot::CandidateSource");
    if(it == candidate.extra.end())
        return false;
    if(const auto source = std::any_cast<String>(&it->second))
        return (*source == "state-cache");
    return false;
}

auto solveConstrainedTernaryInternalProblem(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult;

auto solveConstrainedTernaryInternalProblemWithDiagnostics(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> ConstrainedTernaryMinimizationOutcome;

auto candidateInternalObjective(GlobalizedSolidSolutionBranchSelection const& selection) -> real
{
    const auto it = selection.extra.find(InternalObjectiveKey);
    if(it == selection.extra.end())
        return std::numeric_limits<double>::infinity();
    if(const auto objective = std::any_cast<real>(&it->second))
        return *objective;
    return std::numeric_limits<double>::infinity();
}

auto formatArray(ArrayXrConstRef values) -> String
{
    std::ostringstream out;
    out << "[";
    for(Index i = 0; i < values.size(); ++i)
    {
        if(i)
            out << ", ";
        out << static_cast<double>(values[i]);
    }
    out << "]";
    return out.str();
}

auto makeConstrainedTernarySelection(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    GlobalizedSolidSolutionCandidate const& candidate,
    GlobalizedSolidSolutionInput const& evaluationInput,
    Vec<GlobalizedSolidSolutionBranch> const& evaluationBranches,
    ArrayXrConstRef visiblex,
    real branchTolerance) -> GlobalizedSolidSolutionBranchSelection
{
    const Optional<ArrayXr> warmstart = candidate.initialInternalx.size() == visiblex.size()
        ? Optional<ArrayXr>(candidate.initialInternalx)
        : std::nullopt;

    const auto minimized = solveConstrainedTernaryInternalProblemWithDiagnostics(options, evaluationInput.T, visiblex, warmstart);
    ArrayXr internalx = minimized.result.x;
    internalx = internalx.max(CompositionFloor);
    internalx /= internalx.sum();
    auto requestedSeedProjectionApplied = false;

    if(evaluationInput.requestedBranch != GlobalizedSolidSolutionNoBranch
        && candidate.branch == evaluationInput.requestedBranch
        && candidate.branch < evaluationBranches.size())
    {
        auto requestedSeed = warmstart.has_value()
            ? projectSeedToBranch(ArrayXr(*warmstart), evaluationBranches[candidate.branch])
            : projectSeedToBranch(ArrayXr(visiblex), evaluationBranches[candidate.branch]);

        if(requestedSeed.size() == internalx.size())
        {
            const auto minimizedViolation = GlobalizedSolidSolutionBranchViolation(internalx, evaluationBranches[candidate.branch], branchTolerance);
            const auto seedViolation = GlobalizedSolidSolutionBranchViolation(requestedSeed, evaluationBranches[candidate.branch], branchTolerance);
            if(seedViolation + branchTolerance < minimizedViolation)
            {
                internalx = requestedSeed;
                requestedSeedProjectionApplied = true;
            }
        }
    }

    GlobalizedSolidSolutionBranchSelection selection;
    selection.branch = candidate.branch;
    selection.score = GlobalizedSolidSolutionBranchViolation(internalx, evaluationBranches[candidate.branch], branchTolerance);
    selection.internalx = internalx;
    selection.usedWarmstart = warmstart.has_value();
    selection.reusedState = candidateReusedState(candidate);
    selection.extra["MAGEMinSolidSolutionPilot::InternalMinimizerIterations"] = static_cast<std::uint64_t>(minimized.result.iterations);
    selection.extra["MAGEMinSolidSolutionPilot::InternalMinimizerConverged"] = minimized.result.converged;
    selection.extra[InternalObjectiveKey] = minimized.result.objective;
    selection.extra["MAGEMinSolidSolutionPilot::RequestedBranchProjectedSeedApplied"] = requestedSeedProjectionApplied;

    for(const auto& [key, value] : minimized.extra)
        selection.extra[key] = value;

    for(const auto& [key, value] : candidate.extra)
        selection.extra[key] = value;

    const auto diagnosticsEnabled = []() -> bool
    {
        const auto env = std::getenv("REAKTORO_MAGEMIN_PILOT_DIAG");
        return env && String(env) == "1";
    }();

    if(diagnosticsEnabled && evaluationInput.requestedBranch != GlobalizedSolidSolutionNoBranch)
    {
        const auto seedUsed = warmstart.has_value() ? *warmstart : ArrayXr(visiblex);
        std::cerr
            << "[MAGEMinPilotDiag] model=" << options.thermo.modelId
            << " requestedBranch=" << evaluationInput.requestedBranch
            << " candidateBranch=" << candidate.branch
            << " selectedBranch=" << selection.branch
            << " seed=" << formatArray(seedUsed)
            << " internalx=" << formatArray(selection.internalx)
            << " objective=" << static_cast<double>(minimized.result.objective)
            << " projectedSeedApplied=" << (requestedSeedProjectionApplied ? 1 : 0)
            << "\n";
    }

    return selection;
}

auto constrainedTernarySplitCandidates(
    Vec<PrecomputedConstrainedTernaryCandidateEvaluation> const& evaluations,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionInput const& input,
    real branchTolerance) -> Vec<SolidSolutionCandidateState>
{
    struct BranchBestCandidate
    {
        bool valid = false;
        real objective = std::numeric_limits<double>::infinity();
        ArrayXr internalx;
    };

    Vec<BranchBestCandidate> best(branches.size());
    for(const auto& evaluation : evaluations)
    {
        const auto branchIndex = evaluation.selection.branch;
        if(branchIndex == GlobalizedSolidSolutionNoBranch || branchIndex >= branches.size())
            continue;

        const auto branchViolation = GlobalizedSolidSolutionBranchViolation(
            evaluation.selection.internalx,
            branches[branchIndex],
            branchTolerance);
        if(branchViolation > branchTolerance)
            continue;

        const auto objective = candidateInternalObjective(evaluation.selection);
        if(!best[branchIndex].valid || objective < best[branchIndex].objective)
        {
            best[branchIndex].valid = true;
            best[branchIndex].objective = objective;
            best[branchIndex].internalx = evaluation.selection.internalx;
        }
    }

    auto bestObjective = std::numeric_limits<double>::infinity();
    for(const auto& candidate : best)
        if(candidate.valid)
            bestObjective = std::min(bestObjective, static_cast<double>(candidate.objective));

    if(!std::isfinite(bestObjective))
        return {};

    Vec<SolidSolutionCandidateState> candidates;
    ArrayXr referenceSeed;
    for(Index i = 0; i < branches.size(); ++i)
    {
        if(!best[i].valid)
            continue;

        const auto objectiveGap = static_cast<double>(best[i].objective - bestObjective);
        if(objectiveGap > BranchStabilityObjectiveTolerance)
            continue;

        if(referenceSeed.size() == 0)
        {
            referenceSeed = best[i].internalx;
        }
        else if((best[i].internalx - referenceSeed).matrix().norm() <= BranchStabilitySeedGapTolerance)
        {
            continue;
        }

        candidates.push_back({i, best[i].internalx, objectiveGap, branchCandidateLabel(branches[i], i)});
    }

    if(candidates.size() < 2)
        return {};

    auto visibleBetweenBranches = false;
    for(const auto& candidate : candidates)
    {
        const auto branchViolation = GlobalizedSolidSolutionBranchViolation(input.x, branches[candidate.branch], branchTolerance);
        if(branchViolation > branchTolerance)
        {
            visibleBetweenBranches = true;
            break;
        }
    }

    return visibleBetweenBranches ? candidates : Vec<SolidSolutionCandidateState>{};
}

auto constrainedTernarySplitRequest(
    GlobalizedSolidSolutionInput const& input,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    Indices const& candidateBranches,
    Vec<SolidSolutionCandidateState> const& candidates,
    Index triggeringBranch,
    real branchTolerance,
    String const& splitViolationKey) -> GlobalizedSolidSolutionSplitRequest
{
    // Build the working seed list. Prefer minimizer-found seeds; if none are
    // within branch bounds (bulk in spinodal / two-phase region), synthesize seeds
    // by projecting the bulk composition onto each branch's feasible simplex.
    Vec<SolidSolutionCandidateState> workingCandidates = candidates;

    if(workingCandidates.size() < 2)
    {
        // Determine which branch indices to generate seeds for.
        // Use explicitly nominated candidateBranches if available; otherwise all branches.
        workingCandidates.clear();
        const auto n = input.x.size();
        const auto numBranches = static_cast<Index>(branches.size());
        const bool useCandidateBranches = static_cast<Index>(candidateBranches.size()) >= 2;

        for(Index bi = 0; bi < (useCandidateBranches ? static_cast<Index>(candidateBranches.size()) : numBranches); ++bi)
        {
            const auto branchIndex = useCandidateBranches ? candidateBranches[bi] : bi;
            if(branchIndex >= static_cast<Index>(branches.size()))
                continue;
            const auto& branch = branches[branchIndex];
            auto seed = projectSeedToBranch(ArrayXr(input.x), branch);

            const auto seedSum = static_cast<double>(seed.sum());
            if(seedSum <= 0.0)
                continue;
            seed /= seedSum;
            SolidSolutionCandidateState cs;
            cs.branch = branchIndex;
            cs.seedx = seed;
            cs.priority = 0.0;
            cs.label = branchCandidateLabel(branch, branchIndex) + "-projected";
            workingCandidates.push_back(std::move(cs));
        }
    }

    if(workingCandidates.size() < 2)
        return {};

    auto splitRequest = DefaultGlobalizedSolidSolutionSplitRequest(
        input,
        branches,
        triggeringBranch,
        branchTolerance,
        splitViolationKey);

    // Collect branch indices from workingCandidates (may differ from candidateBranches
    // when seeds were synthesized from all branches).
    Indices workingBranchIndices;
    workingBranchIndices.reserve(workingCandidates.size());
    for(const auto& wc : workingCandidates)
        workingBranchIndices.push_back(wc.branch);

    splitRequest.requested = true;
    splitRequest.triggeringBranch = triggeringBranch;
    splitRequest.branches = workingBranchIndices;
    splitRequest.branchIds.clear();
    for(const auto branchIndex : workingBranchIndices)
        if(branchIndex < static_cast<Index>(branches.size()))
            splitRequest.branchIds.push_back(branches[branchIndex].id);
    splitRequest.reason = "branch-stability-between-branches";
    splitRequest.extra[SplitCandidateStatesKey] = workingCandidates;
    splitRequest.extra[SplitCandidateCountKey] = static_cast<std::uint64_t>(workingCandidates.size());

    auto objectiveGap = 0.0;
    for(const auto& candidate : workingCandidates)
        objectiveGap = std::max(objectiveGap, static_cast<double>(candidate.priority));
    splitRequest.extra[SplitCandidateObjectiveGapKey] = objectiveGap;

    return splitRequest;
}


auto selectBranch(
    MAGEMinSolidSolutionPilotBranchPolicyOptions const& options,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionInput const& input,
    ArrayXrConstRef internalx,
    String const& emptyCandidatesMessage,
    String const& invalidBranchMessage,
    String const& rejectedCandidatesMessage,
    String const& splitViolationKey) -> GlobalizedSolidSolutionBranchSelection
{
    const auto defaultGenerator = [=](GlobalizedSolidSolutionInput const& screeningInput, Vec<GlobalizedSolidSolutionBranch> const& screeningBranches)
    {
        return defaultCandidates(options, screeningBranches, screeningInput, internalx);
    };

    const auto stabilityCriterion = options.stabilityCriterion
        ? options.stabilityCriterion
        : defaultPilotStabilityCriterion(branches, options.branchTolerance, splitViolationKey);

    return ComposeGlobalizedSolidSolutionBranch(
        branches,
        input,
        options.candidateGenerator,
        defaultGenerator,
        stabilityCriterion,
        [=](GlobalizedSolidSolutionCandidate const& candidate, GlobalizedSolidSolutionInput const&, Vec<GlobalizedSolidSolutionBranch> const& evaluationBranches)
        {
            GlobalizedSolidSolutionBranchSelection selection;
            selection.branch = candidate.branch;
            selection.score = GlobalizedSolidSolutionBranchViolation(internalx, evaluationBranches[candidate.branch], options.branchTolerance);
            selection.reusedState = candidateReusedState(candidate);
            return selection;
        },
        {},
        emptyCandidatesMessage,
        invalidBranchMessage,
        rejectedCandidatesMessage);
}

auto regularTernaryExcessChemicalPotentials(
    MAGEMinImportedConstrainedTernarySolutionThermoModel const& thermo,
    ArrayXrConstRef y) -> ArrayXr;

/// Compute the pure (no penalty) Gibbs mixing energy at composition y.
auto constrainedTernaryPureGibbs(
    MAGEMinImportedConstrainedTernarySolutionThermoModel const& thermo,
    real T,
    ArrayXrConstRef y) -> real
{
    const auto muEx = regularTernaryExcessChemicalPotentials(thermo, y);
    const auto Gex = y.matrix().dot(muEx.matrix());
    const auto Gid = thermo.idealGibbs ? thermo.idealGibbs(T, y) : real(0.0);
    return Gex + Gid;
}

/// Compute the gradient of the pure Gibbs mixing energy at composition y (d/dy G_mix, no penalty).
auto constrainedTernaryPureGibbsGradient(
    MAGEMinImportedConstrainedTernarySolutionThermoModel const& thermo,
    real T,
    ArrayXrConstRef y) -> ArrayXr
{
    const auto RT = universalGasConstant * T;
    ArrayXr grad = regularTernaryExcessChemicalPotentials(thermo, y);
    if(thermo.idealLnActivities)
        grad += RT * thermo.idealLnActivities(y);
    return grad;
}

/// Build a tangent-plane distance stability criterion from precomputed branch-local evaluations.
///
/// For each candidate composition ξ, the criterion evaluates the tangent-plane distance
///   TPD(ξ') = G_mix(ξ') - G_mix(ξ) - ∇G_mix(ξ)·(ξ' - ξ)
/// at every competing branch's best local minimum ξ'. A negative TPD (more negative than
/// -tpdTolerance * R * T) indicates instability against splitting toward that branch.
auto constrainedTernaryTangentPlaneStabilityCriterion(
    MAGEMinImportedConstrainedTernarySolutionThermoModel const& thermo,
    Vec<PrecomputedConstrainedTernaryCandidateEvaluation> const& evaluations,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    real T,
    real tpdTolerance,
    real branchTolerance,
    String const& splitViolationKey) -> GlobalizedSolidSolutionCandidateStabilityCriterion
{
    // Pre-collect best internal composition per branch from evaluations.
    struct BranchBest {
        bool valid = false;
        ArrayXr internalx;
        real objective = std::numeric_limits<double>::infinity();
    };
    Vec<BranchBest> best(branches.size());
    for(const auto& ev : evaluations)
    {
        const auto bi = ev.selection.branch;
        if(bi == GlobalizedSolidSolutionNoBranch || bi >= branches.size())
            continue;
        const auto branchViolation = GlobalizedSolidSolutionBranchViolation(
            ev.selection.internalx, branches[bi], branchTolerance);
        if(branchViolation > branchTolerance)
            continue;
        const auto obj = candidateInternalObjective(ev.selection);
        if(!best[bi].valid || obj < best[bi].objective)
        {
            best[bi].valid = true;
            best[bi].objective = obj;
            best[bi].internalx = ev.selection.internalx;
        }
    }

    const auto RT = universalGasConstant * T;
    const auto threshold = -tpdTolerance * RT;

    return [=](GlobalizedSolidSolutionInput const& input,
               GlobalizedSolidSolutionBranch const& /*branch*/,
               ArrayXrConstRef xi,
               real /*score*/) -> GlobalizedSolidSolutionCandidateStability
    {
        if(xi.size() == 0)
            return {};

        const auto Gxi = constrainedTernaryPureGibbs(thermo, T, xi);
        const auto gradGxi = constrainedTernaryPureGibbsGradient(thermo, T, xi);

        double minTPD = std::numeric_limits<double>::infinity();
        Indices unstableBranches;
        Indices allCompetingBranches;

        for(Index bi = 0; bi < static_cast<Index>(best.size()); ++bi)
        {
            if(!best[bi].valid)
                continue;
            const auto& xip = best[bi].internalx;
            if(xip.size() != xi.size())
                continue;

            // Compositions are essentially the same — skip self-check.
            if((xip - xi).matrix().norm() < branchTolerance)
                continue;

            const auto Gxip = constrainedTernaryPureGibbs(thermo, T, xip);
            const auto tpd = static_cast<double>(Gxip - Gxi - gradGxi.matrix().dot((xip - xi).matrix()));

            allCompetingBranches.push_back(bi);
            if(tpd < static_cast<double>(threshold))
            {
                unstableBranches.push_back(bi);
                minTPD = std::min(minTPD, tpd);
            }
        }

        if(unstableBranches.empty())
        {
            GlobalizedSolidSolutionCandidateStability stability;
            stability.extra["MAGEMinSolidSolutionPilot::TPDStable"] = true;
            stability.extra["MAGEMinSolidSolutionPilot::TPDMinValue"] =
                allCompetingBranches.empty() ? real(0.0) : real(minTPD);
            return stability;
        }

        // Collect all unstable + current branch indices for the split request.
        // (The split request builder needs ≥ 2 branches.)
        Vec<SolidSolutionCandidateState> splitCandidates;
        for(const auto bi : unstableBranches)
        {
            SolidSolutionCandidateState cs;
            cs.branch = bi;
            cs.seedx = best[bi].internalx;
            cs.priority = 0.0;
            cs.label = branchCandidateLabel(branches[bi], bi);
            splitCandidates.push_back(cs);
        }
        // Also add the candidate's own branch so we have ≥ 2 entries.
        {
            const auto triggerBranch = GlobalizedSolidSolutionBranchViolation(xi, branches.empty() ? GlobalizedSolidSolutionBranch{} : branches[0], branchTolerance) < branchTolerance ? Index(0) : GlobalizedSolidSolutionNoBranch;
            // Find the branch that xi belongs to.
            Index xiBranch = GlobalizedSolidSolutionNoBranch;
            for(Index bi = 0; bi < static_cast<Index>(branches.size()); ++bi)
            {
                if(GlobalizedSolidSolutionBranchViolation(xi, branches[bi], branchTolerance) < branchTolerance)
                {
                    xiBranch = bi;
                    break;
                }
            }
            (void)triggerBranch;
            const auto alreadyPresent = std::any_of(splitCandidates.begin(), splitCandidates.end(),
                [xiBranch](const SolidSolutionCandidateState& c){ return c.branch == xiBranch; });
            if(!alreadyPresent && xiBranch != GlobalizedSolidSolutionNoBranch)
            {
                SolidSolutionCandidateState cs;
                cs.branch = xiBranch;
                cs.seedx = ArrayXr(xi);
                cs.priority = 0.0;
                cs.label = branchCandidateLabel(branches[xiBranch], xiBranch);
                splitCandidates.push_back(cs);
            }
        }

        Indices splitBranchIndices;
        Strings splitBranchIds;
        for(const auto& sc : splitCandidates)
        {
            splitBranchIndices.push_back(sc.branch);
            if(sc.branch < static_cast<Index>(branches.size()))
                splitBranchIds.push_back(branches[sc.branch].id);
        }

        GlobalizedSolidSolutionSplitRequest splitRequest;
        splitRequest.requested = true;
        splitRequest.triggeringBranch = splitCandidates.empty() ? GlobalizedSolidSolutionNoBranch : splitCandidates.front().branch;
        splitRequest.branches = splitBranchIndices;
        splitRequest.branchIds = splitBranchIds;
        splitRequest.reason = "tangent-plane-instability";
        splitRequest.extra[SplitCandidateStatesKey] = splitCandidates;
        splitRequest.extra[SplitCandidateCountKey] = static_cast<std::uint64_t>(splitCandidates.size());
        splitRequest.extra["MAGEMinSolidSolutionPilot::TPDMinValue"] = real(minTPD);
        splitRequest.extra["MAGEMinSolidSolutionPilot::TPDStable"] = false;

        GlobalizedSolidSolutionCandidateStability stability;
        stability.stable = false;
        stability.reason = "tangent-plane-instability";
        stability.splitRequest = splitRequest;
        stability.extra["MAGEMinSolidSolutionPilot::TPDStable"] = false;
        stability.extra["MAGEMinSolidSolutionPilot::TPDMinValue"] = real(minTPD);
        return stability;
    };
}

auto selectConstrainedTernaryBranch(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    MAGEMinSolidSolutionPilotBranchPolicyOptions const& branchPolicy,
    Vec<GlobalizedSolidSolutionBranch> const& branches,
    GlobalizedSolidSolutionInput const& input,
    ArrayXrConstRef visiblex,
    String const& emptyCandidatesMessage,
    String const& invalidBranchMessage,
    String const& rejectedCandidatesMessage,
    String const& splitViolationKey) -> GlobalizedSolidSolutionBranchSelection
{
    const auto defaultGenerator = [=](GlobalizedSolidSolutionInput const& screeningInput, Vec<GlobalizedSolidSolutionBranch> const& screeningBranches)
    {
        return augmentDefaultConstrainedTernaryCandidates(
            options.thermo,
            options.proposals,
            screeningBranches,
            visiblex,
            defaultCandidates(branchPolicy, screeningBranches, screeningInput, visiblex));
    };

    const auto generator = branchPolicy.candidateGenerator
        ? branchPolicy.candidateGenerator
        : defaultGenerator;
    auto candidates = generator(input, branches);
    if(candidates.empty())
        throw std::runtime_error(emptyCandidatesMessage);

    Vec<PrecomputedConstrainedTernaryCandidateEvaluation> evaluations(candidates.size());
    for(Index i = 0; i < candidates.size(); ++i)
    {
        candidates[i].extra[PrecomputedCandidateIndexKey] = static_cast<std::uint64_t>(i);
        evaluations[i].selection = makeConstrainedTernarySelection(
            options,
            candidates[i],
            input,
            branches,
            visiblex,
            branchPolicy.branchTolerance);
    }

    const auto splitCandidates = constrainedTernarySplitCandidates(evaluations, branches, input, branchPolicy.branchTolerance);
    Indices splitCandidateBranches;
    splitCandidateBranches.reserve(splitCandidates.size());
    for(const auto& candidate : splitCandidates)
        splitCandidateBranches.push_back(candidate.branch);

    GlobalizedSolidSolutionCandidateStabilityCriterion stabilityCriterion;
    if(branchPolicy.stabilityCriterion)
    {
        stabilityCriterion = branchPolicy.stabilityCriterion;
    }
    else if(options.enableTangentPlaneStabilityCheck)
    {
        // Use the thermodynamically grounded tangent-plane distance criterion.
        stabilityCriterion = constrainedTernaryTangentPlaneStabilityCriterion(
            options.thermo, evaluations, branches,
            static_cast<double>(input.T), options.tpdTolerance,
            branchPolicy.branchTolerance, splitViolationKey);
    }
    else
    {
        stabilityCriterion = defaultPilotStabilityCriterion(branches, branchPolicy.branchTolerance, splitViolationKey);
    }

    return ComposeGlobalizedSolidSolutionBranch(
        branches,
        input,
        [candidates](GlobalizedSolidSolutionInput const&, Vec<GlobalizedSolidSolutionBranch> const&) { return candidates; },
        defaultGenerator,
        stabilityCriterion,
        [evaluations](GlobalizedSolidSolutionCandidate const& candidate, GlobalizedSolidSolutionInput const&, Vec<GlobalizedSolidSolutionBranch> const&)
        {
            const auto it = candidate.extra.find(PrecomputedCandidateIndexKey);
            if(it == candidate.extra.end())
                throw std::runtime_error("Precomputed constrained ternary candidate is missing its evaluation index.");
            const auto* index = std::any_cast<std::uint64_t>(&it->second);
            if(!index || *index >= evaluations.size())
                throw std::runtime_error("Precomputed constrained ternary candidate has an invalid evaluation index.");
            return evaluations[static_cast<Index>(*index)].selection;
        },
        [=](GlobalizedSolidSolutionInput const& splitInput, Vec<GlobalizedSolidSolutionBranch> const& splitBranches, Index triggeringBranch)
        {
            return constrainedTernarySplitRequest(
                splitInput,
                splitBranches,
                splitCandidateBranches,
                splitCandidates,
                triggeringBranch,
                branchPolicy.branchTolerance,
                splitViolationKey);
        },
        emptyCandidatesMessage,
        invalidBranchMessage,
        rejectedCandidatesMessage);
}

auto regularTernaryExcessChemicalPotentials(MAGEMinImportedConstrainedTernarySolutionThermoModel const& thermo, ArrayXrConstRef y) -> ArrayXr
{
    if(thermo.excessChemicalPotentials)
        return thermo.excessChemicalPotentials(y);

    ArrayXr mu(3);
    const auto y0 = static_cast<double>(y[0]);
    const auto y1 = static_cast<double>(y[1]);
    const auto y2 = static_cast<double>(y[2]);

    mu[0] = y1*(1.0 - y0)*thermo.W01 + y2*(1.0 - y0)*thermo.W02 - y1*y2*thermo.W12;
    mu[1] = y0*(1.0 - y1)*thermo.W01 - y0*y2*thermo.W02 + y2*(1.0 - y1)*thermo.W12;
    mu[2] = -y0*y1*thermo.W01 + y0*(1.0 - y2)*thermo.W02 + y1*(1.0 - y2)*thermo.W12;
    return mu;
}

auto normalizedInternalComposition(ArrayXr y) -> ArrayXr
{
    if(y.size() == 0)
        return y;

    for(Index iter = 0; iter < 8; ++iter)
    {
        y = y.max(CompositionFloor);
        const auto sum = static_cast<double>(y.sum());
        if(sum <= 0.0)
            y = ArrayXr::Constant(y.size(), 1.0/static_cast<double>(y.size()));
        else y /= sum;
    }

    return y;
}

auto hasBounds(ArrayXrConstRef lowerBounds, ArrayXrConstRef upperBounds, Index size) -> bool
{
    const auto hasLower = lowerBounds.size() != 0;
    const auto hasUpper = upperBounds.size() != 0;
    if(hasLower != hasUpper)
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: lowerBounds and upperBounds must be both empty or both populated.");
    if(!hasLower)
        return false;

    if(lowerBounds.size() != size || upperBounds.size() != size)
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: bounds size must match composition size.");

    for(Index i = 0; i < size; ++i)
        if(lowerBounds[i] > upperBounds[i])
            throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: lowerBounds must not exceed upperBounds.");

    return true;
}

auto projectToBoundedSimplex(ArrayXrConstRef v, ArrayXrConstRef lowerBounds, ArrayXrConstRef upperBounds) -> ArrayXr
{
    const auto sumLower = static_cast<double>(lowerBounds.sum());
    const auto sumUpper = static_cast<double>(upperBounds.sum());
    if(sumLower > 1.0 || sumUpper < 1.0)
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: bounded simplex is infeasible (sum(lowerBounds) <= 1 <= sum(upperBounds) violated).");

    auto lo = std::numeric_limits<double>::infinity();
    auto hi = -std::numeric_limits<double>::infinity();
    for(Index i = 0; i < v.size(); ++i)
    {
        lo = std::min(lo, static_cast<double>(v[i] - upperBounds[i]));
        hi = std::max(hi, static_cast<double>(v[i] - lowerBounds[i]));
    }

    auto projected = ArrayXr(v.size());
    for(Index iter = 0; iter < 80; ++iter)
    {
        const auto lambda = 0.5*(lo + hi);
        auto sum = 0.0;
        for(Index i = 0; i < v.size(); ++i)
        {
            const auto value = std::clamp(
                static_cast<double>(v[i] - lambda),
                static_cast<double>(lowerBounds[i]),
                static_cast<double>(upperBounds[i]));
            projected[i] = value;
            sum += value;
        }

        if(sum > 1.0)
            lo = lambda;
        else
            hi = lambda;
    }

    return projected;
}

auto applyLocalModelConstraints(
    MAGEMinConstrainedTernaryLocalModel const& model,
    ArrayXr y) -> ArrayXr
{
    if(y.size() == 0)
        return y;

    const auto bounded = hasBounds(model.lowerBounds, model.upperBounds, y.size());

    if(model.enforceUnityConstraint)
    {
        ArrayXr lowerBounds = ArrayXr::Constant(y.size(), CompositionFloor);
        ArrayXr upperBounds = ArrayXr::Constant(y.size(), real(1.0) - CompositionFloor);
        if(bounded)
        {
            lowerBounds = model.lowerBounds.max(CompositionFloor);
            upperBounds = model.upperBounds.min(real(1.0) - CompositionFloor);
        }

        for(Index i = 0; i < y.size(); ++i)
        {
            if(lowerBounds[i] > upperBounds[i])
                throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: constrained bounds are invalid after applying composition floor.");
        }

        return projectToBoundedSimplex(y, lowerBounds, upperBounds);
    }

    if(bounded)
        return y.max(model.lowerBounds).min(model.upperBounds);

    return y;
}

/// Check if the local model has nonlinear constraint callbacks defined.
auto hasNonlinearConstraints(MAGEMinConstrainedTernaryLocalModel const& model) -> bool
{
    return static_cast<bool>(model.constraints);
}

/// Validate constraint callback consistency and bounds.
auto validateConstraintCallbacks(MAGEMinConstrainedTernaryLocalModel const& model) -> void
{
    if(model.constraintPenaltyWeight < 0.0)
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: constraintPenaltyWeight must be non-negative.");
    if(model.constraintBarrierWeight < 0.0)
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: constraintBarrierWeight must be non-negative.");
    if(model.trustRegionRadius < 0.0)
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: trustRegionRadius must be non-negative.");

    if(!hasNonlinearConstraints(model))
        return;  // No constraints to validate

    // If constraints callback exists, Jacobian should also be provided for solver robustness
    if(model.constraints && !model.constraintJacobian)
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: "
            "constraintJacobian callback must be provided if constraints callback is set.");

    // Ensure constraint bounds are consistent
    if(model.constraintLowerBounds.size() != model.constraintUpperBounds.size())
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: "
            "constraintLowerBounds and constraintUpperBounds must have the same size.");

    if(model.constraintLowerBounds.size() == 0)
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: constraint bounds must be populated when constraints callback is set.");

    // If Hessian is provided and useSecondOrderInfo is true, ensure consistency
    if(model.useSecondOrderInfo && !model.objectiveHessian)
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: "
            "useSecondOrderInfo=true requires objectiveHessian callback to be defined.");
}

/// Evaluate nonlinear constraint feasibility at given composition.
auto evaluateConstraintFeasibility(
    MAGEMinConstrainedTernaryLocalModel const& model,
    ArrayXrConstRef y) -> bool
{
    if(!hasNonlinearConstraints(model))
        return true;  // No constraints = feasible

    const auto c = model.constraints(y);
    for(Index i = 0; i < c.size(); ++i)
    {
        if(c[i] > model.constraintUpperBounds[i] + 1.0e-12 ||
           c[i] < model.constraintLowerBounds[i] - 1.0e-12)
            return false;  // Constraint violation detected
    }
    return true;
}

auto nonlinearConstraintViolation(
    MAGEMinConstrainedTernaryLocalModel const& model,
    ArrayXrConstRef y) -> real
{
    if(!hasNonlinearConstraints(model))
        return 0.0;

    const auto c = model.constraints(y);
    if(c.size() != model.constraintLowerBounds.size())
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: constraint callback size must match constraint bounds size.");

    real violation = 0.0;
    for(Index i = 0; i < c.size(); ++i)
    {
        const auto lowerResidual = std::max(real(0.0), model.constraintLowerBounds[i] - c[i]);
        const auto upperResidual = std::max(real(0.0), c[i] - model.constraintUpperBounds[i]);
        violation += lowerResidual*lowerResidual + upperResidual*upperResidual;
    }
    return violation;
}

auto nonlinearConstraintBarrier(
    MAGEMinConstrainedTernaryLocalModel const& model,
    ArrayXrConstRef y) -> real
{
    if(!hasNonlinearConstraints(model) || model.constraintBarrierWeight <= 0.0)
        return 0.0;

    const auto c = model.constraints(y);
    if(c.size() != model.constraintLowerBounds.size())
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: constraint callback size must match constraint bounds size.");

    real barrier = 0.0;
    for(Index i = 0; i < c.size(); ++i)
    {
        const auto lowerSlack = c[i] - model.constraintLowerBounds[i];
        const auto upperSlack = model.constraintUpperBounds[i] - c[i];
        if(lowerSlack <= 0.0 || upperSlack <= 0.0)
            return std::numeric_limits<double>::infinity();
        barrier -= model.constraintBarrierWeight*(log(lowerSlack) + log(upperSlack));
    }
    return barrier;
}

auto localModelMerit(
    MAGEMinConstrainedTernaryLocalModel const& model,
    ArrayXrConstRef y) -> real
{
    const auto objective = model.objective(y);
    const auto violation = nonlinearConstraintViolation(model, y);
    const auto barrier = nonlinearConstraintBarrier(model, y);
    return objective + model.constraintPenaltyWeight*violation + barrier;
}

auto applyTrustRegion(
    MAGEMinConstrainedTernaryLocalModel const& model,
    ArrayXrConstRef step) -> ArrayXr
{
    if(model.trustRegionRadius <= 0.0)
        return ArrayXr(step);

    const auto norm = step.matrix().norm();
    if(norm <= model.trustRegionRadius || norm == 0.0)
        return ArrayXr(step);

    return (step * (model.trustRegionRadius / norm)).eval();
}

auto validateTCMConstraintBridge(MAGEMinTCMConstraintBridge const& bridge, Index modelSize) -> void
{
    if(bridge.numConstraints == 0)
        throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: bridge.numConstraints must be greater than zero.");
    if(bridge.numVariables == 0)
        throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: bridge.numVariables must be greater than zero.");
    if(!bridge.nativeToVisible && static_cast<Index>(bridge.numVariables) != modelSize)
        throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: bridge.numVariables must match local-model composition size.");
    if(!bridge.callback)
        throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: bridge.callback must be populated.");
    if(bridge.constraintLowerBounds.size() != bridge.constraintUpperBounds.size())
        throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: bridge constraint bounds must have matching sizes.");
    if(bridge.constraintLowerBounds.size() != static_cast<Index>(bridge.numConstraints))
        throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: bridge constraint bounds size must match bridge.numConstraints.");
    if((bridge.variableLowerBounds.size() == 0) != (bridge.variableUpperBounds.size() == 0))
        throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: native variable bounds must both be empty or both be populated.");
    if(bridge.variableLowerBounds.size() && bridge.variableLowerBounds.size() != static_cast<Index>(bridge.numVariables))
        throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: native variable bounds size must match bridge.numVariables.");
}

auto evaluateTCMConstraintValues(
    MAGEMinTCMConstraintBridge const& bridge,
    ArrayXrConstRef y) -> ArrayXr
{
    std::vector<double> x(bridge.numVariables, 0.0);
    for(unsigned j = 0; j < bridge.numVariables; ++j)
        x[j] = static_cast<double>(y[static_cast<Index>(j)]);

    std::vector<double> values(bridge.numConstraints, 0.0);
    bridge.callback(
        bridge.numConstraints,
        values.data(),
        bridge.numVariables,
        x.data(),
        nullptr,
        bridge.userData);

    ArrayXr c(static_cast<Index>(bridge.numConstraints));
    for(unsigned i = 0; i < bridge.numConstraints; ++i)
        c[static_cast<Index>(i)] = values[i];
    return c;
}

auto evaluateTCMConstraintJacobian(
    MAGEMinTCMConstraintBridge const& bridge,
    ArrayXrConstRef y) -> MatrixXr
{
    std::vector<double> x(bridge.numVariables, 0.0);
    for(unsigned j = 0; j < bridge.numVariables; ++j)
        x[j] = static_cast<double>(y[static_cast<Index>(j)]);

    std::vector<double> values(bridge.numConstraints, 0.0);
    std::vector<double> gradient(bridge.numConstraints*bridge.numVariables, 0.0);
    bridge.callback(
        bridge.numConstraints,
        values.data(),
        bridge.numVariables,
        x.data(),
        gradient.data(),
        bridge.userData);

    MatrixXr J(static_cast<Index>(bridge.numConstraints), static_cast<Index>(bridge.numVariables));
    for(unsigned i = 0; i < bridge.numConstraints; ++i)
        for(unsigned j = 0; j < bridge.numVariables; ++j)
            J(static_cast<Index>(i), static_cast<Index>(j)) = gradient[i*bridge.numVariables + j];
    return J;
}

auto finiteDifferenceNativeGradient(
    Fn<real(ArrayXrConstRef)> const& objective,
    ArrayXrConstRef nativex) -> ArrayXr
{
    ArrayXr g(nativex.size());
    for(Index i = 0; i < nativex.size(); ++i)
    {
        const auto xi = static_cast<double>(nativex[i]);
        const auto h = NativeGradientFiniteDifferenceStep*std::max(1.0, std::abs(xi));

        auto xp = ArrayXr(nativex);
        auto xm = ArrayXr(nativex);
        xp[i] = xi + h;
        xm[i] = xi - h;

        g[i] = (objective(xp) - objective(xm))/(2.0*h);
    }
    return g;
}

auto finiteDifferenceMapJacobian(
    Fn<ArrayXr(ArrayXrConstRef)> const& mapNativeToVisible,
    ArrayXrConstRef nativex,
    Index visibleSize) -> MatrixXr
{
    MatrixXr J = MatrixXr::Zero(visibleSize, nativex.size());
    for(Index j = 0; j < nativex.size(); ++j)
    {
        const auto xj = static_cast<double>(nativex[j]);
        const auto h = NativeGradientFiniteDifferenceStep*std::max(1.0, std::abs(xj));

        auto xp = ArrayXr(nativex);
        auto xm = ArrayXr(nativex);
        xp[j] = xj + h;
        xm[j] = xj - h;

        const auto yp = mapNativeToVisible(xp);
        const auto ym = mapNativeToVisible(xm);
        if(yp.size() != visibleSize || ym.size() != visibleSize)
            throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: nativeToVisible returned an unexpected vector size during Jacobian estimation.");

        J.col(j) = ((yp - ym)/(2.0*h)).matrix();
    }
    return J;
}

struct ActiveConstraintLinearization
{
    MatrixXr A;
    ArrayXr b;
    Vec<Index> constraintIndices;
    Vec<real> constraintSigns;
};

auto buildActiveConstraintLinearization(
    MAGEMinConstrainedTernaryLocalModel const& model,
    ArrayXrConstRef y) -> ActiveConstraintLinearization
{
    ActiveConstraintLinearization linearization;

    const auto n = y.size();
    const auto activeTolerance = std::max(1.0e-8, 10.0*static_cast<double>(model.tolerance));

    if(model.enforceUnityConstraint)
    {
        linearization.constraintIndices.push_back(-1);
        linearization.constraintSigns.push_back(1.0);
    }

    if(hasNonlinearConstraints(model))
    {
        const auto c = model.constraints(y);
        const auto J = model.constraintJacobian(y);

        if(c.size() != model.constraintLowerBounds.size())
            throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: constraint callback size must match constraint bounds size.");
        if(J.rows() != c.size() || J.cols() != n)
            throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: constraintJacobian dimensions must be m x n with m equal to the constraint count.");

        for(Index i = 0; i < c.size(); ++i)
        {
            const auto upperSlack = model.constraintUpperBounds[i] - c[i];
            if(upperSlack <= activeTolerance || c[i] > model.constraintUpperBounds[i])
            {
                linearization.constraintIndices.push_back(i);
                linearization.constraintSigns.push_back(1.0);
            }

            const auto lowerSlack = c[i] - model.constraintLowerBounds[i];
            if(lowerSlack <= activeTolerance || c[i] < model.constraintLowerBounds[i])
            {
                linearization.constraintIndices.push_back(i);
                linearization.constraintSigns.push_back(-1.0);
            }
        }

        linearization.A = MatrixXr::Zero(linearization.constraintIndices.size(), n);
        linearization.b = ArrayXr::Zero(linearization.constraintIndices.size());

        Index row = 0;
        if(model.enforceUnityConstraint)
        {
            linearization.A.row(row).setOnes();
            linearization.b[row] = 1.0 - y.sum();
            ++row;
        }

        for(Index k = model.enforceUnityConstraint ? 1 : 0; k < linearization.constraintIndices.size(); ++k, ++row)
        {
            const auto i = linearization.constraintIndices[k];
            const auto sign = linearization.constraintSigns[k];
            linearization.A.row(row) = sign*J.row(i);
            linearization.b[row] = sign > 0.0
                ? model.constraintUpperBounds[i] - c[i]
                : c[i] - model.constraintLowerBounds[i];
        }

        return linearization;
    }

    linearization.A = MatrixXr::Zero(linearization.constraintIndices.size(), n);
    linearization.b = ArrayXr::Zero(linearization.constraintIndices.size());

    if(model.enforceUnityConstraint)
    {
        linearization.A.row(0).setOnes();
        linearization.b[0] = 1.0 - y.sum();
    }

    return linearization;
}

auto regularizeLocalModelHessian(MatrixXr H) -> MatrixXr
{
    H = 0.5*(H + H.transpose()).eval();
    const auto diagonalScale = std::max(1.0, static_cast<double>(H.diagonal().cwiseAbs().maxCoeff()));
    H.diagonal().array() += 1.0e-8*diagonalScale;
    return H;
}

auto mapConstraintMultipliers(
    ActiveConstraintLinearization const& linearization,
    ArrayXrConstRef activeMultipliers,
    Index nonlinearConstraintCount) -> ArrayXr
{
    ArrayXr multipliers = ArrayXr::Zero(nonlinearConstraintCount);
    for(Index row = 0; row < activeMultipliers.size(); ++row)
    {
        const auto constraintIndex = linearization.constraintIndices[row];
        if(constraintIndex < 0)
            continue;

        multipliers[constraintIndex] = multipliers[constraintIndex] + linearization.constraintSigns[row]*activeMultipliers[row];
    }
    return multipliers;
}

auto constrainedSteepestDescentDirection(
    ArrayXrConstRef gradient,
    ActiveConstraintLinearization const& linearization) -> ArrayXr
{
    auto direction = (-gradient).eval();

    if(linearization.A.rows() == 0)
        return direction;

    MatrixXr normal = (linearization.A*linearization.A.transpose()).eval();
    normal.diagonal().array() += 1.0e-10;

    const auto residual = (linearization.A*direction.matrix()).array() - linearization.b;
    const auto correctionMultipliers = normal.partialPivLu().solve(residual.matrix()).array();
    direction -= (linearization.A.transpose()*correctionMultipliers.matrix()).array();
    return direction;
}

auto computeConstraintAwareSearchDirection(
    MAGEMinConstrainedTernaryLocalModel const& model,
    ArrayXrConstRef y,
    ArrayXrConstRef gradient,
    ArrayXrConstRef multiplierEstimate) -> std::pair<ArrayXr, ArrayXr>
{
    const auto linearization = buildActiveConstraintLinearization(model, y);
    const auto nonlinearConstraintCount = model.constraintLowerBounds.size();

    if(model.useSecondOrderInfo && model.objectiveHessian)
    {
        MatrixXr H = regularizeLocalModelHessian(model.objectiveHessian(y, multiplierEstimate));
        if(H.rows() != y.size() || H.cols() != y.size())
            throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: objectiveHessian must return an n x n matrix.");

        if(linearization.A.rows() == 0)
        {
            const auto Hd = H.template cast<double>();
            const auto rhsd = (-gradient).matrix().template cast<double>();
            const auto directiond = Hd.partialPivLu().solve(rhsd);
            const auto direction = directiond.array().template cast<real>();
            return std::make_pair(direction, ArrayXr::Zero(nonlinearConstraintCount));
        }

        const auto n = y.size();
        const auto m = linearization.A.rows();
        MatrixXr kkt = MatrixXr::Zero(n + m, n + m);
        kkt.topLeftCorner(n, n) = H;
        kkt.topRightCorner(n, m) = linearization.A.transpose();
        kkt.bottomLeftCorner(m, n) = linearization.A;

        ArrayXr rhs = ArrayXr::Zero(n + m);
        rhs.head(n) = -gradient;
        rhs.tail(m) = linearization.b;

        const auto kktd = kkt.template cast<double>();
        const auto rhsd = rhs.matrix().template cast<double>();
        const auto solutiond = kktd.partialPivLu().solve(rhsd);
        const auto direction = solutiond.head(n).array().template cast<real>();
        const auto activeMultipliers = solutiond.tail(m).array().template cast<real>();
        return std::make_pair(direction, mapConstraintMultipliers(linearization, activeMultipliers, nonlinearConstraintCount));
    }

    return std::make_pair(
        constrainedSteepestDescentDirection(gradient, linearization),
        ArrayXr::Zero(nonlinearConstraintCount)
    );
}

auto constrainedTernaryObjective(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    ArrayXrConstRef y) -> real
{
    const auto muEx = regularTernaryExcessChemicalPotentials(options.thermo, y);
    const auto Gex = y.matrix().dot(muEx.matrix());
    const auto Gid = options.thermo.idealGibbs ? options.thermo.idealGibbs(T, y) : real(0.0);
    return Gex + Gid + options.externalCompositionPenalty*universalGasConstant*T*(y - x).matrix().squaredNorm();
}

auto constrainedTernaryGradient(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    ArrayXrConstRef y) -> ArrayXr
{
    const auto RT = universalGasConstant*T;
    ArrayXr gradient = regularTernaryExcessChemicalPotentials(options.thermo, y);
    if(options.thermo.idealLnActivities)
        gradient += RT*options.thermo.idealLnActivities(y);
    gradient += 2.0*options.externalCompositionPenalty*RT*(y - x);
    return gradient;
}

auto makeConstrainedTernaryLocalModel(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x) -> MAGEMinConstrainedTernaryLocalModel
{
    MAGEMinConstrainedTernaryLocalModel model;
    model.modelId = options.thermo.modelId;
    model.T = T;
    model.visiblex = ArrayXr(x);
    model.objective = [=](ArrayXrConstRef y) -> real { return constrainedTernaryObjective(options, T, x, y); };
    model.gradient = [=](ArrayXrConstRef y) -> ArrayXr { return constrainedTernaryGradient(options, T, x, y); };
    model.lowerBounds = ArrayXr::Constant(x.size(), CompositionFloor);
    model.upperBounds = ArrayXr::Constant(x.size(), 1.0 - CompositionFloor);
    model.enforceUnityConstraint = true;
    model.tolerance = options.minimizerTolerance;
    model.maxIterations = options.minimizerMaxIterations;
    return model;
}

auto projectedGradientConstrainedTernaryMinimizer(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
{
    auto current = normalizedInternalComposition(warmstart.value_or(ArrayXr(x)));
    auto currentObjective = constrainedTernaryObjective(options, T, x, current);

    Index iterations = 0;
    auto converged = false;

    for(; iterations < options.minimizerMaxIterations; ++iterations)
    {
        const auto gradient = constrainedTernaryGradient(options, T, x, current);
        auto projectedGradient = gradient.array() - gradient.mean();

        if(projectedGradient.cwiseAbs().maxCoeff() <= options.minimizerTolerance)
        {
            converged = true;
            break;
        }

        auto step = 0.25;
        auto accepted = false;
        const auto referenceSlope = projectedGradient.matrix().squaredNorm();

        for(Index backtrack = 0; backtrack < 20; ++backtrack)
        {
            auto trial = normalizedInternalComposition(current - step*projectedGradient);
            const auto trialObjective = constrainedTernaryObjective(options, T, x, trial);

            if(trialObjective <= currentObjective - ProjectedGradientArmijo*step*referenceSlope)
            {
                current = std::move(trial);
                currentObjective = trialObjective;
                accepted = true;
                break;
            }

            step *= 0.5;
        }

        if(!accepted)
        {
            converged = true; // Armijo failure: stuck at local minimum
            break;
        }
    }

    return {
        current,
        currentObjective,
        iterations,
        converged,
    };
}

auto defaultConstrainedTernaryMinimizer(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult;

auto constrainedTernaryBuiltinMinimizerOutcome(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> ConstrainedTernaryMinimizationOutcome
{
    ConstrainedTernaryMinimizationOutcome outcome;

    if(options.defaultMinimizerStrategy != BuiltinProjectedGradientMinimizerStrategy)
    {
        outcome.result = defaultConstrainedTernaryMinimizer(options, T, x, warmstart);
        outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String(BuiltinLegacyMinimizerStrategy);
        outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount"] = static_cast<std::uint64_t>(0);
        outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(outcome.result.iterations);
        return outcome;
    }

    const auto projected = projectedGradientConstrainedTernaryMinimizer(options, T, x, warmstart);
    outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount"] = static_cast<std::uint64_t>(projected.iterations);
    outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientConverged"] = projected.converged;

    if(!options.compareProjectedGradientAgainstLegacy)
    {
        outcome.result = projected;
        outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String(BuiltinProjectedGradientMinimizerStrategy);
        outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = projected.converged;
        outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(0);
        return outcome;
    }

    const auto legacy = defaultConstrainedTernaryMinimizer(options, T, x, warmstart);
    const auto sameShape = projected.x.size() == legacy.x.size();
    const auto compositionDelta = sameShape
        ? static_cast<double>((projected.x - legacy.x).cwiseAbs().maxCoeff())
        : std::numeric_limits<double>::infinity();
    const auto objectiveDelta = std::abs(static_cast<double>(projected.objective - legacy.objective));
    const auto agreement = projected.converged
        && sameShape
        && compositionDelta <= ProjectedGradientAgreementTolerance
        && objectiveDelta <= ProjectedGradientAgreementTolerance;
    const auto fallbackToLegacy = options.fallbackToLegacyOnProjectedGradientDisagreement && !agreement;

    outcome.result = fallbackToLegacy ? legacy : projected;
    outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String(fallbackToLegacy ? BuiltinLegacyMinimizerStrategy : BuiltinProjectedGradientMinimizerStrategy);
    outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = true;
    outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = fallbackToLegacy;
    outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = agreement;
    outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(legacy.iterations);
    outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientLegacyCompositionDelta"] = compositionDelta;
    outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientLegacyObjectiveDelta"] = objectiveDelta;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientHasLowerObjective"] = bool(projected.objective < legacy.objective);
    return outcome;
}

auto defaultConstrainedTernaryMinimizer(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
{
    GlobalizedSolidSolutionInternalProblem problem;
    problem.objective = [=](ArrayXrConstRef y) -> real { return constrainedTernaryObjective(options, T, x, y); };
    problem.initialx = warmstart.value_or(x);
    problem.lowerBounds = ArrayXr::Constant(x.size(), CompositionFloor);
    problem.upperBounds = ArrayXr::Constant(x.size(), 1.0 - CompositionFloor);
    problem.tolerance = options.minimizerTolerance;
    problem.maxIterations = options.minimizerMaxIterations;
    problem.enforceUnityConstraint = true;
    return MinimizeGlobalizedSolidSolutionInternalProblem(problem);
}

auto solveConstrainedTernaryInternalProblem(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
{
    if(options.preferNLoptLocalModelMinimizer && options.nloptLocalModelMinimizer)
        return options.nloptLocalModelMinimizer(makeConstrainedTernaryLocalModel(options, T, x), warmstart);

    if(options.localModelMinimizer)
        return options.localModelMinimizer(makeConstrainedTernaryLocalModel(options, T, x), warmstart);

    if(options.nloptLocalModelMinimizer)
        return options.nloptLocalModelMinimizer(makeConstrainedTernaryLocalModel(options, T, x), warmstart);

    if(options.minimizer)
        return options.minimizer(options, T, x, warmstart);

    return defaultConstrainedTernaryMinimizer(options, T, x, warmstart);
}

auto solveConstrainedTernaryInternalProblemWithDiagnostics(
    MAGEMinImportedConstrainedTernarySolutionOptions const& options,
    real T,
    ArrayXrConstRef x,
    Optional<ArrayXr> warmstart) -> ConstrainedTernaryMinimizationOutcome
{
    if(options.preferNLoptLocalModelMinimizer && options.nloptLocalModelMinimizer)
    {
        ConstrainedTernaryMinimizationOutcome outcome;
        const auto model = makeConstrainedTernaryLocalModel(options, T, x);
        outcome.result = options.nloptLocalModelMinimizer(model, warmstart);
        outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String("custom-local-model-nlopt");
        outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount"] = static_cast<std::uint64_t>(0);
        outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(0);

        if(options.localModelDiagnostics)
        {
            const auto payload = options.localModelDiagnostics(model, outcome.result);
            for(const auto& [key, value] : payload)
                outcome.extra[key] = value;
        }

        return outcome;
    }

    if(options.localModelMinimizer)
    {
        ConstrainedTernaryMinimizationOutcome outcome;
        const auto model = makeConstrainedTernaryLocalModel(options, T, x);
        outcome.result = options.localModelMinimizer(model, warmstart);
        outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String("custom-local-model");
        outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount"] = static_cast<std::uint64_t>(0);
        outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(0);

        if(options.localModelDiagnostics)
        {
            const auto payload = options.localModelDiagnostics(model, outcome.result);
            for(const auto& [key, value] : payload)
                outcome.extra[key] = value;
        }

        return outcome;
    }

    if(options.nloptLocalModelMinimizer)
    {
        ConstrainedTernaryMinimizationOutcome outcome;
        const auto model = makeConstrainedTernaryLocalModel(options, T, x);
        outcome.result = options.nloptLocalModelMinimizer(model, warmstart);
        outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String("custom-local-model-nlopt");
        outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount"] = static_cast<std::uint64_t>(0);
        outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(0);

        if(options.localModelDiagnostics)
        {
            const auto payload = options.localModelDiagnostics(model, outcome.result);
            for(const auto& [key, value] : payload)
                outcome.extra[key] = value;
        }

        return outcome;
    }

    if(options.minimizer)
    {
        ConstrainedTernaryMinimizationOutcome outcome;
        outcome.result = options.minimizer(options, T, x, warmstart);
        outcome.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String("custom");
        outcome.extra["MAGEMinSolidSolutionPilot::ComparedAgainstLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::FallbackToLegacy"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientAccepted"] = false;
        outcome.extra["MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount"] = static_cast<std::uint64_t>(0);
        outcome.extra["MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount"] = static_cast<std::uint64_t>(0);
        return outcome;
    }

    return constrainedTernaryBuiltinMinimizerOutcome(options, T, x, warmstart);
}

} // namespace

auto MAGEMinProjectedGradientLocalModelMinimizer(
    MAGEMinConstrainedTernaryLocalModel const& model,
    Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
{
    if(!model.gradient)
        throw std::runtime_error("MAGEMinProjectedGradientLocalModelMinimizer: model.gradient callback must be populated.");

    // Validate nonlinear constraint callbacks if provided
    validateConstraintCallbacks(model);

    auto current = applyLocalModelConstraints(model, warmstart.value_or(ArrayXr(model.visiblex)));
    auto currentObjective = model.objective(current);
    auto currentMerit = localModelMerit(model, current);
    ArrayXr constraintMultipliers = ArrayXr::Zero(model.constraintLowerBounds.size());

    Index iterations = 0;
    auto converged = false;

    for(; iterations < model.maxIterations; ++iterations)
    {
        const auto gradient = model.gradient(current);
        ArrayXr fallbackDirection = gradient;
        if(model.enforceUnityConstraint)
            fallbackDirection = (gradient.array() - gradient.mean()).eval();

        auto searchDirection = (-fallbackDirection).eval();
        if(hasNonlinearConstraints(model))
        {
            const auto stepData = computeConstraintAwareSearchDirection(
                model,
                current,
                gradient,
                constraintMultipliers);
            searchDirection = stepData.first;
            constraintMultipliers = stepData.second;

            if(!searchDirection.cast<double>().matrix().allFinite() || searchDirection.matrix().norm() == 0.0)
                searchDirection = (-fallbackDirection).eval();
        }

        const auto violation = nonlinearConstraintViolation(model, current);
        const auto directionNorm = searchDirection.matrix().norm();
        if((hasNonlinearConstraints(model) && directionNorm <= model.tolerance && violation <= model.tolerance) ||
           (!hasNonlinearConstraints(model) && fallbackDirection.cwiseAbs().maxCoeff() <= model.tolerance))
        {
            converged = true;
            break;
        }

        auto step = 0.25;
        auto accepted = false;
        auto referenceSlope = static_cast<double>((-gradient.matrix().dot(searchDirection.matrix())));
        if(!(referenceSlope > 0.0))
            referenceSlope = std::max(1.0e-16, static_cast<double>(searchDirection.matrix().squaredNorm()));

        for(Index backtrack = 0; backtrack < 20; ++backtrack)
        {
            auto trialStep = applyTrustRegion(model, (step*searchDirection).eval());
            auto trial = applyLocalModelConstraints(model, current + trialStep);
            if(model.requireFeasibleTrialPoints && !evaluateConstraintFeasibility(model, trial))
            {
                step *= 0.5;
                continue;
            }

            const auto trialObjective = model.objective(trial);
            const auto trialMerit = localModelMerit(model, trial);

            if(trialMerit <= currentMerit - ProjectedGradientArmijo*step*referenceSlope)
            {
                current = std::move(trial);
                currentObjective = trialObjective;
                currentMerit = trialMerit;
                accepted = true;
                break;
            }

            step *= 0.5;
        }

        if(!accepted)
        {
            converged = true; // Armijo failure: stuck at local minimum
            break;
        }
    }

    return {
        current,
        currentObjective,
        iterations,
        converged,
    };
}

auto MAGEMinTCMConstraintBridgeLocalModelAdapter(
    MAGEMinTCMConstraintBridge bridge,
    MAGEMinConstrainedTernaryLocalModelMinimizer fallback) -> MAGEMinConstrainedTernaryLocalModelMinimizer
{
    return [bridge = std::move(bridge), fallback = std::move(fallback)](
        MAGEMinConstrainedTernaryLocalModel const& model,
        Optional<ArrayXr> warmstart) -> GlobalizedSolidSolutionInternalResult
    {
        validateTCMConstraintBridge(bridge, model.visiblex.size());

        const auto mapNativeToVisible = [&](ArrayXrConstRef nativex) -> ArrayXr
        {
            if(bridge.nativeToVisible)
            {
                auto mapped = bridge.nativeToVisible(nativex);
                if(mapped.size() == model.visiblex.size())
                    return normalizedVisibleFromNative(std::move(mapped), model.visiblex.size());

                if(mapped.size() <= 0 || mapped.size() > model.visiblex.size())
                    throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: nativeToVisible returned an unexpected vector size.");

                ArrayXr expanded = ArrayXr::Constant(model.visiblex.size(), CompositionFloor);
                expanded.head(mapped.size()) = mapped;
                return normalizedVisibleFromNative(std::move(expanded), model.visiblex.size());
            }

            if(nativex.size() != model.visiblex.size())
                throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: native and visible variable sizes differ but no nativeToVisible mapping was provided.");

            return ArrayXr(nativex);
        };

        Optional<ArrayXr> nativeWarmstart = std::nullopt;
        if(warmstart)
        {
            if(warmstart->size() == static_cast<Index>(bridge.numVariables))
                nativeWarmstart = *warmstart;
            else if(bridge.visibleToNative && warmstart->size() == model.visiblex.size())
                nativeWarmstart = bridge.visibleToNative(*warmstart);
            else if(warmstart->size() == model.visiblex.size() && static_cast<Index>(bridge.numVariables) == model.visiblex.size())
                nativeWarmstart = *warmstart;
        }

        if(!nativeWarmstart && bridge.visibleToNative)
            nativeWarmstart = bridge.visibleToNative(model.visiblex);

        auto constrainedModel = model;
        constrainedModel.visiblex = nativeWarmstart
            ? ArrayXr(*nativeWarmstart)
            : ArrayXr::Constant(static_cast<Index>(bridge.numVariables), 1.0/static_cast<double>(bridge.numVariables));

        constrainedModel.objective = [model, mapNativeToVisible](ArrayXrConstRef nativex) -> real
        {
            return model.objective(mapNativeToVisible(nativex));
        };

        constrainedModel.gradient = [bridge, model, mapNativeToVisible](ArrayXrConstRef nativex) -> ArrayXr
        {
            const auto wrappedObjective = [&](ArrayXrConstRef xn) -> real
            {
                return model.objective(mapNativeToVisible(xn));
            };

            if(!model.gradient)
                return finiteDifferenceNativeGradient(wrappedObjective, nativex);

            const auto visiblex = mapNativeToVisible(nativex);
            const auto visibleGradient = model.gradient(visiblex);

            if(visibleGradient.size() != visiblex.size())
                throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: visible gradient size must match nativeToVisible output size.");

            if(bridge.nativeToVisibleJacobian)
            {
                const auto J = bridge.nativeToVisibleJacobian(nativex);
                if(J.rows() != visibleGradient.size() || J.cols() != nativex.size())
                    throw std::runtime_error("MAGEMinTCMConstraintBridgeLocalModelAdapter: nativeToVisibleJacobian dimensions must be visible_size x numVariables.");
                return (J.transpose()*visibleGradient.matrix()).array();
            }

            if(bridge.nativeToVisible)
            {
                const auto J = finiteDifferenceMapJacobian(mapNativeToVisible, nativex, visibleGradient.size());
                return (J.transpose()*visibleGradient.matrix()).array();
            }

            if(visibleGradient.size() == nativex.size())
                return visibleGradient;

            return finiteDifferenceNativeGradient(wrappedObjective, nativex);
        };

        constrainedModel.constraints = [bridge](ArrayXrConstRef nativex) -> ArrayXr
        {
            return evaluateTCMConstraintValues(bridge, nativex);
        };
        constrainedModel.constraintJacobian = [bridge](ArrayXrConstRef nativex) -> MatrixXr
        {
            return evaluateTCMConstraintJacobian(bridge, nativex);
        };
        constrainedModel.constraintLowerBounds = bridge.constraintLowerBounds;
        constrainedModel.constraintUpperBounds = bridge.constraintUpperBounds;
        if(bridge.variableLowerBounds.size())
            constrainedModel.lowerBounds = bridge.variableLowerBounds;
        else constrainedModel.lowerBounds = ArrayXr::Constant(static_cast<Index>(bridge.numVariables), CompositionFloor);

        if(bridge.variableUpperBounds.size())
            constrainedModel.upperBounds = bridge.variableUpperBounds;
        else constrainedModel.upperBounds = ArrayXr::Constant(static_cast<Index>(bridge.numVariables), 1.0 - CompositionFloor);

        const auto defaultUnityConstraint = bridge.nativeToVisible ? false : model.enforceUnityConstraint;
        constrainedModel.enforceUnityConstraint = bridge.enforceUnityConstraint.value_or(defaultUnityConstraint);
        constrainedModel.requireFeasibleTrialPoints = model.requireFeasibleTrialPoints || static_cast<bool>(bridge.nativeToVisible);

        const auto minimizer = fallback ? fallback : MAGEMinProjectedGradientLocalModelMinimizer;
        return minimizer(constrainedModel, nativeWarmstart);
    };
}

auto MAGEMinSolidSolutionPilotModelImportedBinary(
    MAGEMinImportedBinarySolutionOptions options) -> GlobalizedSolidSolutionModel
{
    const auto thermo = options.thermo;
    auto branchPolicy = options.branchPolicy;
    branchPolicy.branches = normalizedMAGEMinPilotBranches(branchPolicy.branches, 2);

    return [=](GlobalizedSolidSolutionInput input)
    {
        if(input.x.size() != 2)
            throw std::runtime_error("MAGEMin imported binary pilot model requires exactly two species.");

        auto state = input.state ? input.state : std::make_shared<GlobalizedSolidSolutionState>();

        ArrayXr internalx = input.x;
        internalx[0] = std::clamp(static_cast<double>(internalx[0]), CompositionFloor, 1.0 - CompositionFloor);
        internalx[1] = 1.0 - internalx[0];

        const auto selected = selectBranch(
            branchPolicy,
            branchPolicy.branches,
            input,
            internalx,
            "MAGEMin imported binary pilot candidate generator returned no candidates.",
            "MAGEMin imported binary pilot candidate generator returned an invalid branch.",
            "MAGEMin imported binary pilot stability screen rejected all branch candidates.",
            "MAGEMinSolidSolutionPilot::SplitViolation");
        const auto selectedBranch = selected.branch;

        const auto RT = universalGasConstant * input.T;
        const auto x0 = static_cast<double>(internalx[0]);
        const auto x1 = static_cast<double>(internalx[1]);
        const auto muEx0 = x1*x1*thermo.W;
        const auto muEx1 = x0*x0*thermo.W;

        GlobalizedSolidSolutionOutput output;
        output.branches = branchPolicy.branches;
        output.selectedBranch = selectedBranch;
        output.branch = branchPolicy.branches[selectedBranch];
        output.state = state;
        output.som = StateOfMatter::Solid;
        output.Vxi = ArrayXr::Zero(2);
        output.ln_g = ArrayXr::Zero(2);
        output.ln_a = ArrayXr::Zero(2);

        output.ln_g[0] = muEx0/RT;
        output.ln_g[1] = muEx1/RT;
        output.ln_a = output.ln_g + thermo.idealSiteMultiplicity*log(internalx);

        output.Gx = x0*x1*thermo.W;
        output.Hx = output.Gx;
        output.splitRequest = selected.splitRequest;

        const auto stateid = CurrentGlobalizedSolidSolutionChemicalPropsStateId(input.extra);
        state->chemicalPropsStateId = stateid;
        state->selectedBranch = selectedBranch;
        state->cachedBranchForState = selectedBranch;
        state->cachedInternalx = internalx;
        state->numEvaluations += 1;
        state->lastT = input.T;
        state->lastP = input.P;
        state->lastx = input.x;
        state->lastInternalx = internalx;
        state->lastSplitRequest = output.splitRequest;
        state->data["MAGEMinSolidSolutionPilot::ModelId"] = thermo.modelId;
        state->data["MAGEMinSolidSolutionPilot::Endmember0"] = thermo.endmember0;
        state->data["MAGEMinSolidSolutionPilot::Endmember1"] = thermo.endmember1;

        output.extra["MAGEMinSolidSolutionPilot::ModelId"] = thermo.modelId;
        output.extra["MAGEMinSolidSolutionPilot::Endmember0"] = thermo.endmember0;
        output.extra["MAGEMinSolidSolutionPilot::Endmember1"] = thermo.endmember1;
        output.extra["MAGEMinSolidSolutionPilot::InternalComposition"] = internalx;
        output.extra["MAGEMinSolidSolutionPilot::BinaryInteractionParameterW"] = thermo.W;
        output.extra["MAGEMinSolidSolutionPilot::IdealSiteMultiplicity"] = thermo.idealSiteMultiplicity;
        output.extra["MAGEMinSolidSolutionPilot::UsedStateCache"] = selected.reusedState;
        for(const auto& [key, value] : selected.extra)
            output.extra[key] = value;
        output.extra["GlobalizedSolidSolution::SplitRequested"] = output.splitRequest.requested;
        output.extra["GlobalizedSolidSolution::SplitRequest"] = output.splitRequest;
        if(output.splitRequest.requested)
            output.extra["GlobalizedSolidSolution::SplitReason"] = output.splitRequest.reason;

        if(input.requestedBranch != GlobalizedSolidSolutionNoBranch)
        {
            std::cerr
                << "[MAGEMinPilotDiag] model=" << thermo.modelId
                << " requestedBranch=" << input.requestedBranch
                << " selectedBranch=" << output.selectedBranch
                << " finalInternalx=" << formatArray(internalx)
                << " splitRequested=" << (output.splitRequest.requested ? 1 : 0)
                << "\n";
        }

        return output;
    };
}

auto MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(
    MAGEMinImportedConstrainedTernarySolutionOptions options) -> GlobalizedSolidSolutionModel
{
    const auto thermo = options.thermo;
    const auto numEndmembers = static_cast<Index>(thermo.endmembers.size());
    auto branchPolicy = options.branchPolicy;
    branchPolicy.branches = normalizedMAGEMinPilotBranches(branchPolicy.branches, numEndmembers);

    return [=](GlobalizedSolidSolutionInput input)
    {
        if(input.x.size() != numEndmembers)
            throw std::runtime_error("MAGEMin imported constrained ternary pilot model: species count (" + std::to_string(input.x.size()) + ") does not match endmember count (" + std::to_string(numEndmembers) + ").");

        auto state = input.state ? input.state : std::make_shared<GlobalizedSolidSolutionState>();

        ArrayXr visiblex = input.x;
        visiblex = visiblex.max(CompositionFloor);
        visiblex /= visiblex.sum();

        const auto selected = selectConstrainedTernaryBranch(
            options,
            branchPolicy,
            branchPolicy.branches,
            input,
            visiblex,
            "MAGEMin imported constrained ternary pilot candidate generator returned no candidates.",
            "MAGEMin imported constrained ternary pilot candidate generator returned an invalid branch.",
            "MAGEMin imported constrained ternary pilot stability screen rejected all branch candidates.",
            "MAGEMinSolidSolutionPilot::SplitViolation");

        ArrayXr internalx = selected.internalx;
        internalx = internalx.max(CompositionFloor);
        internalx /= internalx.sum();

        const auto muEx = regularTernaryExcessChemicalPotentials(thermo, internalx);
        const auto RT = universalGasConstant * input.T;
        const auto Gex = internalx.matrix().dot(muEx.matrix());
        const auto Gid = thermo.idealGibbs ? thermo.idealGibbs(input.T, internalx) : real(0.0);
        const auto idealLnA = thermo.idealLnActivities ? thermo.idealLnActivities(internalx) : ArrayXr::Zero(3);

        GlobalizedSolidSolutionOutput output;
        output.branches = branchPolicy.branches;
        output.selectedBranch = selected.branch;
        output.branch = branchPolicy.branches[selected.branch];
        output.state = state;
        output.som = StateOfMatter::Solid;
        output.Vxi = ArrayXr::Zero(numEndmembers);
        output.ln_g = muEx/RT;
        output.ln_a = output.ln_g + idealLnA;
        output.Gx = Gex + Gid;
        output.Hx = Gex;
        output.splitRequest = selected.splitRequest;

        const auto stateid = CurrentGlobalizedSolidSolutionChemicalPropsStateId(input.extra);
        state->chemicalPropsStateId = stateid;
        state->selectedBranch = selected.branch;
        state->cachedBranchForState = selected.branch;
        state->cachedInternalx = internalx;
        state->numEvaluations += 1;
        state->lastT = input.T;
        state->lastP = input.P;
        state->lastx = input.x;
        state->lastInternalx = internalx;
        state->lastSplitRequest = output.splitRequest;
        state->data["MAGEMinSolidSolutionPilot::ModelId"] = thermo.modelId;
        state->data["MAGEMinSolidSolutionPilot::Endmembers"] = thermo.endmembers;

        output.extra["MAGEMinSolidSolutionPilot::ModelId"] = thermo.modelId;
        output.extra["MAGEMinSolidSolutionPilot::Endmembers"] = thermo.endmembers;
        output.extra["MAGEMinSolidSolutionPilot::InternalComposition"] = internalx;
        output.extra["MAGEMinSolidSolutionPilot::UsedStateCache"] = selected.reusedState;
        output.extra["MAGEMinSolidSolutionPilot::UsedWarmstart"] = selected.usedWarmstart;
        output.extra["MAGEMinSolidSolutionPilot::BinaryInteraction01"] = thermo.W01;
        output.extra["MAGEMinSolidSolutionPilot::BinaryInteraction02"] = thermo.W02;
        output.extra["MAGEMinSolidSolutionPilot::BinaryInteraction12"] = thermo.W12;
        for(const auto& [key, value] : selected.extra)
            output.extra[key] = value;
        output.extra["GlobalizedSolidSolution::SplitRequested"] = output.splitRequest.requested;
        output.extra["GlobalizedSolidSolution::SplitRequest"] = output.splitRequest;
        if(output.splitRequest.requested)
            output.extra["GlobalizedSolidSolution::SplitReason"] = output.splitRequest.reason;

        return output;
    };
}

auto MAGEMinSolidSolutionPilotModelSB11Olivine(
    MAGEMinSB11OlivineOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedBinarySolutionOptions imported;
    imported.thermo = sb11OlivineThermo();
    imported.branchPolicy = options.branchPolicy;
    return MAGEMinSolidSolutionPilotModelImportedBinary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB11Wadsleyite(
    MAGEMinSB11WadsleyiteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedBinarySolutionOptions imported;
    imported.thermo = sb11WadsleyiteThermo();
    imported.branchPolicy = options.branchPolicy;
    return MAGEMinSolidSolutionPilotModelImportedBinary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB11Akimotoite(
    MAGEMinSB11AkimotoiteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb11AkimotoiteThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1};
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;
    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.compareProjectedGradientAgainstLegacy = true;
    imported.fallbackToLegacyOnProjectedGradientDisagreement = true;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB11Perovskite(
    MAGEMinSB11PerovskiteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb11PerovskiteThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1};
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;
    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.compareProjectedGradientAgainstLegacy = true;
    imported.fallbackToLegacyOnProjectedGradientDisagreement = true;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB11Calcioferrite(
    MAGEMinSB11CalcioferriteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb11CalcioferriteThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1};
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;
    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.compareProjectedGradientAgainstLegacy = true;
    imported.fallbackToLegacyOnProjectedGradientDisagreement = true;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21Spinel(
    MAGEMinSB21SpinelOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedBinarySolutionOptions imported;
    imported.thermo = sb21SpinelThermo();
    imported.branchPolicy = options.branchPolicy;
    return MAGEMinSolidSolutionPilotModelImportedBinary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21NAL(
    MAGEMinSB21NALOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb21NALThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1}; // nnal-first: largest W interactions
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;
    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.compareProjectedGradientAgainstLegacy = true;
    imported.fallbackToLegacyOnProjectedGradientDisagreement = true;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21Calcioferrite(
    MAGEMinSB21CalcioferriteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb21CalcioferriteThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1};
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;
    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    imported.enableTangentPlaneStabilityCheck = options.enableTangentPlaneStabilityCheck;
    imported.tpdTolerance = options.tpdTolerance;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21OPX(
    MAGEMinSB21OPXOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb21OPXThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {3, 2, 0, 1}; // odi-first: largest W interactions
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;

    if(!imported.nloptLocalModelMinimizer)
    {
        const auto bridge = options.tcMConstraintBridge ? *options.tcMConstraintBridge : sb21AkimotoiteTCMConstraintBridge();
        imported.nloptLocalModelMinimizer = MAGEMinTCMConstraintBridgeLocalModelAdapter(bridge);
        if(!options.localModelMinimizer)
            imported.preferNLoptLocalModelMinimizer = true;
    }

    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.compareProjectedGradientAgainstLegacy = true;
    imported.fallbackToLegacyOnProjectedGradientDisagreement = true;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21CPX(
    MAGEMinSB21CPXOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb21CPXThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {3, 2, 4, 0, 1}; // cats-first: largest W interactions
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;

    if(!imported.nloptLocalModelMinimizer)
    {
        const auto bridge = options.tcMConstraintBridge ? *options.tcMConstraintBridge : sb21CPXTCMConstraintBridge();
        imported.nloptLocalModelMinimizer = MAGEMinTCMConstraintBridgeLocalModelAdapter(bridge);
        if(!options.localModelMinimizer)
            imported.preferNLoptLocalModelMinimizer = true;
    }

    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.compareProjectedGradientAgainstLegacy = true;
    imported.fallbackToLegacyOnProjectedGradientDisagreement = true;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21GTMJ(
    MAGEMinSB21GTMJOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb21GTMJThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {3, 4, 2, 0, 1}; // mgmj-first: largest W interactions
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;

    if(!imported.nloptLocalModelMinimizer)
    {
        const auto bridge = options.tcMConstraintBridge ? *options.tcMConstraintBridge : sb21PerovskiteTCMConstraintBridge();
        imported.nloptLocalModelMinimizer = MAGEMinTCMConstraintBridgeLocalModelAdapter(bridge);
        if(!options.localModelMinimizer)
            imported.preferNLoptLocalModelMinimizer = true;
    }

    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.compareProjectedGradientAgainstLegacy = true;
    imported.fallbackToLegacyOnProjectedGradientDisagreement = true;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21PLG(
    MAGEMinSB21PLGOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedBinarySolutionOptions imported;
    imported.thermo = sb21PLGThermo();
    imported.branchPolicy = options.branchPolicy;
    return MAGEMinSolidSolutionPilotModelImportedBinary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21Olivine(
    MAGEMinSB21OlivineOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedBinarySolutionOptions imported;
    imported.thermo = sb21OlivineThermo();
    imported.branchPolicy = options.branchPolicy;
    return MAGEMinSolidSolutionPilotModelImportedBinary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21Wadsleyite(
    MAGEMinSB21WadsleyiteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedBinarySolutionOptions imported;
    imported.thermo = sb21WadsleyiteThermo();
    imported.branchPolicy = options.branchPolicy;
    return MAGEMinSolidSolutionPilotModelImportedBinary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21Ringwoodite(
    MAGEMinSB21RingwooditeOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedBinarySolutionOptions imported;
    imported.thermo = sb21RingwooditeThermo();
    imported.branchPolicy = options.branchPolicy;
    return MAGEMinSolidSolutionPilotModelImportedBinary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21HPCPX(
    MAGEMinSB21HPCPXOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedBinarySolutionOptions imported;
    imported.thermo = sb21HPCPXThermo();
    imported.branchPolicy = options.branchPolicy;
    return MAGEMinSolidSolutionPilotModelImportedBinary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21Akimotoite(
    MAGEMinSB21AkimotoiteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb21AkimotoiteThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1}; // co-first: largest W interactions
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;

    if(!imported.nloptLocalModelMinimizer)
    {
        const auto bridge = options.tcMConstraintBridge ? *options.tcMConstraintBridge : sb21PostPerovskiteTCMConstraintBridge();
        imported.nloptLocalModelMinimizer = MAGEMinTCMConstraintBridgeLocalModelAdapter(bridge);
        if(!options.localModelMinimizer)
            imported.preferNLoptLocalModelMinimizer = true;
    }

    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21Perovskite(
    MAGEMinSB21PerovskiteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb21PerovskiteThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1}; // alpv-first: largest W interactions
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;

    if(!imported.nloptLocalModelMinimizer)
    {
        const auto bridge = options.tcMConstraintBridge ? *options.tcMConstraintBridge : sb21MagnesiowustitesTCMConstraintBridge();
        imported.nloptLocalModelMinimizer = MAGEMinTCMConstraintBridgeLocalModelAdapter(bridge);
        if(!options.localModelMinimizer)
            imported.preferNLoptLocalModelMinimizer = true;
    }

    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21PostPerovskite(
    MAGEMinSB21PostPerovskiteOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb21PostPerovskiteThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {2, 0, 1}; // appv-first: largest W interactions
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;

    if(!imported.nloptLocalModelMinimizer && options.tcMConstraintBridge)
    {
        imported.nloptLocalModelMinimizer = MAGEMinTCMConstraintBridgeLocalModelAdapter(*options.tcMConstraintBridge);
        if(!options.localModelMinimizer)
            imported.preferNLoptLocalModelMinimizer = true;
    }

    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelSB21Magnesiowustites(
    MAGEMinSB21MagnesiowustitesOptions options) -> GlobalizedSolidSolutionModel
{
    MAGEMinImportedConstrainedTernarySolutionOptions imported;
    imported.thermo = sb21MagnesiowustitesThermo();
    imported.branchPolicy = options.branchPolicy;
    imported.proposals.dominantEndmemberOrder = {0, 1, 2}; // pe-first
    imported.minimizer = options.minimizer;
    imported.localModelMinimizer = options.localModelMinimizer;
    imported.nloptLocalModelMinimizer = options.nloptLocalModelMinimizer;
    imported.preferNLoptLocalModelMinimizer = options.preferNLoptLocalModelMinimizer;

    if(!imported.nloptLocalModelMinimizer && options.tcMConstraintBridge)
    {
        imported.nloptLocalModelMinimizer = MAGEMinTCMConstraintBridgeLocalModelAdapter(*options.tcMConstraintBridge);
        if(!options.localModelMinimizer)
            imported.preferNLoptLocalModelMinimizer = true;
    }

    imported.localModelDiagnostics = options.localModelDiagnostics;
    imported.defaultMinimizerStrategy = BuiltinProjectedGradientMinimizerStrategy;
    imported.externalCompositionPenalty = options.externalCompositionPenalty;
    imported.minimizerTolerance = options.minimizerTolerance;
    imported.minimizerMaxIterations = options.minimizerMaxIterations;
    return MAGEMinSolidSolutionPilotModelImportedConstrainedTernary(imported);
}

auto MAGEMinSolidSolutionPilotModelHPIGOPX(
    MAGEMinHPIGOPXOptions options) -> GlobalizedSolidSolutionModel
{
    auto branchPolicy = options.branchPolicy;
    branchPolicy.branches = normalizedMAGEMinPilotBranches(branchPolicy.branches, 8);

    const auto lowerBounds = hpIGOPXLowerBounds();
    const auto upperBounds = hpIGOPXUpperBounds();
    const auto bridge = options.tcMConstraintBridge ? *options.tcMConstraintBridge : hpIGOPXTCMConstraintBridge();
    const auto referenceState = options.referenceState ? options.referenceState : hpIGOPXReferenceState;

    return [=](GlobalizedSolidSolutionInput input)
    {
        if(input.x.size() != 8)
            throw std::runtime_error("MAGEMin HP ig_opx xeos-native pilot model requires exactly eight coordinates.");

        auto state = input.state ? input.state : std::make_shared<GlobalizedSolidSolutionState>();

        ArrayXr visiblex = input.x;
        for(Index i = 0; i < visiblex.size(); ++i)
            visiblex[i] = std::clamp(static_cast<double>(visiblex[i]), static_cast<double>(lowerBounds[i]), static_cast<double>(upperBounds[i]));

        const auto selected = selectBranch(
            branchPolicy,
            branchPolicy.branches,
            input,
            visiblex,
            "MAGEMin HP ig_opx candidate generator returned no candidates.",
            "MAGEMin HP ig_opx candidate generator returned an invalid branch.",
            "MAGEMin HP ig_opx stability screen rejected all branch candidates.",
            "MAGEMinSolidSolutionPilot::SplitViolation");

        MAGEMinConstrainedTernaryLocalModel localModel;
        localModel.modelId = "ig_opx";
        localModel.T = input.T;
        localModel.visiblex = visiblex;
        localModel.objective = [=](ArrayXrConstRef y) -> real
        {
            return hpIGOPXObjective(input.T, input.P, visiblex, y, options.externalCompositionPenalty, referenceState);
        };
        localModel.gradient = [=](ArrayXrConstRef y) -> ArrayXr
        {
            return hpIGOPXObjectiveGradientFiniteDifference(
                input.T,
                input.P,
                visiblex,
                y,
                lowerBounds,
                upperBounds,
                options.externalCompositionPenalty,
                referenceState);
        };
        localModel.lowerBounds = lowerBounds;
        localModel.upperBounds = upperBounds;
        localModel.enforceUnityConstraint = false;
        localModel.tolerance = options.minimizerTolerance;
        localModel.maxIterations = options.minimizerMaxIterations;

        validateTCMConstraintBridge(bridge, localModel.visiblex.size());
        localModel.constraints = [bridge](ArrayXrConstRef y) -> ArrayXr
        {
            return evaluateTCMConstraintValues(bridge, y);
        };
        localModel.constraintJacobian = [bridge](ArrayXrConstRef y) -> MatrixXr
        {
            return evaluateTCMConstraintJacobian(bridge, y);
        };
        localModel.constraintLowerBounds = bridge.constraintLowerBounds;
        localModel.constraintUpperBounds = bridge.constraintUpperBounds;

        Optional<ArrayXr> warmstart = std::nullopt;
        if(input.state && input.state->lastInternalx.size() == visiblex.size())
            warmstart = input.state->lastInternalx;
        else warmstart = visiblex;

        GlobalizedSolidSolutionInternalResult internalResult;
        if(options.preferNLoptLocalModelMinimizer && options.nloptLocalModelMinimizer)
            internalResult = options.nloptLocalModelMinimizer(localModel, warmstart);
        else if(options.localModelMinimizer)
            internalResult = options.localModelMinimizer(localModel, warmstart);
        else if(options.nloptLocalModelMinimizer)
            internalResult = options.nloptLocalModelMinimizer(localModel, warmstart);
        else
            internalResult = MAGEMinProjectedGradientLocalModelMinimizer(localModel, warmstart);

        ArrayXr internalx = internalResult.x;
        if(internalx.size() != visiblex.size())
            internalx = visiblex;
        for(Index i = 0; i < internalx.size(); ++i)
            internalx[i] = std::clamp(static_cast<double>(internalx[i]), static_cast<double>(lowerBounds[i]), static_cast<double>(upperBounds[i]));

        const auto gradient = localModel.gradient(internalx);
        const auto RT = universalGasConstant * input.T;

        GlobalizedSolidSolutionOutput output;
        output.branches = branchPolicy.branches;
        output.selectedBranch = selected.branch;
        output.branch = branchPolicy.branches[output.selectedBranch];
        output.state = state;
        output.som = StateOfMatter::Solid;
        output.Vxi = ArrayXr::Zero(internalx.size());
        output.ln_a = gradient / RT;
        output.ln_g = output.ln_a;
        output.Gx = internalResult.objective;
        output.Hx = output.Gx;
        output.splitRequest = selected.splitRequest;

        const auto stateid = CurrentGlobalizedSolidSolutionChemicalPropsStateId(input.extra);
        state->chemicalPropsStateId = stateid;
        state->selectedBranch = output.selectedBranch;
        state->cachedBranchForState = output.selectedBranch;
        state->cachedInternalx = internalx;
        state->numEvaluations += 1;
        state->lastT = input.T;
        state->lastP = input.P;
        state->lastx = input.x;
        state->lastInternalx = internalx;
        state->lastSplitRequest = output.splitRequest;

        output.extra["MAGEMinSolidSolutionPilot::ModelId"] = String("ig_opx");
        output.extra["MAGEMinSolidSolutionPilot::CoordinateSystem"] = String("xeos");
        output.extra["MAGEMinSolidSolutionPilot::InternalComposition"] = internalx;
        output.extra["MAGEMinSolidSolutionPilot::InternalMinimizerIterations"] = static_cast<std::uint64_t>(internalResult.iterations);
        output.extra["MAGEMinSolidSolutionPilot::InternalMinimizerConverged"] = internalResult.converged;
        output.extra["MAGEMinSolidSolutionPilot::InternalObjective"] = internalResult.objective;
        output.extra["MAGEMinSolidSolutionPilot::ReferenceState"] = referenceState(input.T, input.P);

        if(options.preferNLoptLocalModelMinimizer && options.nloptLocalModelMinimizer)
            output.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String("custom-local-model-nlopt");
        else if(options.localModelMinimizer)
            output.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String("custom-local-model");
        else if(options.nloptLocalModelMinimizer)
            output.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String("custom-local-model-nlopt");
        else
            output.extra["MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy"] = String(BuiltinProjectedGradientMinimizerStrategy);

        if(options.localModelDiagnostics)
        {
            const auto payload = options.localModelDiagnostics(localModel, internalResult);
            for(const auto& [key, value] : payload)
                output.extra[key] = value;
        }

        for(const auto& [key, value] : selected.extra)
            output.extra[key] = value;

        output.extra["GlobalizedSolidSolution::SplitRequested"] = output.splitRequest.requested;
        output.extra["GlobalizedSolidSolution::SplitRequest"] = output.splitRequest;
        if(output.splitRequest.requested)
            output.extra["GlobalizedSolidSolution::SplitReason"] = output.splitRequest.reason;

        return output;
    };
}

auto MAGEMinSolidSolutionPilotPhase(
    Phase const& phase,
    GlobalizedSolidSolutionModel model) -> Phase
{
    auto pilot = phase.clone();
    pilot = pilot.withActivityModel(ActivityModelGlobalizedSolidSolution(model, phase.name())(phase.species()));
    return pilot;
}

auto MAGEMinSolidSolutionPilotDefinition(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    MAGEMinSolidSolutionPilotOptions options) -> GlobalizedSolidSolutionPhaseDefinition
{
    return MakeGlobalizedSolidSolutionPhaseDefinition(
        MAGEMinSolidSolutionPilotPhase(phase, model),
        model,
        options.branches,
        options.suffixSeparator);
}

auto MAGEMinSolidSolutionPilotPhases(
    Phase const& phase,
    GlobalizedSolidSolutionModel model,
    MAGEMinSolidSolutionPilotOptions options) -> PhaseList
{
    return DuplicateGlobalizedSolidSolutionPhaseBranches(
        MAGEMinSolidSolutionPilotPhase(phase, model),
        model,
        options.branches,
        options.suffixSeparator);
}

} // namespace Reaktoro