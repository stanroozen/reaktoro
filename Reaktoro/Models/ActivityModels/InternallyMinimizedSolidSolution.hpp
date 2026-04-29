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

#pragma once

// C++ includes
#include <any>
#include <cstdint>

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

//=============================================================================
// Core internally-minimized solid-solution types
//=============================================================================

/// Input for evaluating an internally minimized solid solution.
///
/// This is the phase-level problem input: given a visible (outer) mol-fraction
/// composition and a thermodynamic state, compute the internally optimal
/// site-fraction or endmember-coordinate composition and the resulting reduced
/// thermodynamics (Gx, ln_g, ln_a).
struct InternalSolidSolutionInput
{
    /// Temperature (in K).
    real T = 0.0;

    /// Pressure (in Pa).
    real P = 0.0;

    /// Externally visible mol-fraction composition as seen by the outer Reaktoro solver.
    ArrayXr visiblex;

    /// Optional warm-start internal composition from the previous evaluation.
    /// When provided, the internal minimizer should use it as an initial guess.
    Optional<ArrayXr> warmstart;
};

/// Thermodynamic output and diagnostics from an internally minimized solid solution.
///
/// The returned Gx, ln_g, ln_a must be derived from the envelope-theorem form of
/// the internally minimized free energy so that the outer Reaktoro solver sees
/// consistent reduced thermodynamics.
struct InternalSolidSolutionResult
{
    /// The minimized internal composition (site fractions or endmember coordinates).
    ArrayXr internalx;

    /// Corrective molar Gibbs energy after internal minimization (in J/mol).
    real Gx = 0.0;

    /// Corrective molar enthalpy after internal minimization (in J/mol).
    real Hx = 0.0;

    /// Activity coefficients in natural log.
    ArrayXr ln_g;

    /// Activities in natural log.
    ArrayXr ln_a;

    /// Index of the branch selected during this evaluation.
    /// The sentinel value static_cast<Index>(-1) means no branch was selected.
    Index branchIndex = static_cast<Index>(-1);

    /// String identifier of the selected branch, if the model tracks branches.
    String branchId;

    /// Number of internal minimizer iterations performed.
    Index minimizerIterations = 0;

    /// Whether the internal minimizer converged within the requested tolerance.
    bool minimizerConverged = false;

    /// Strategy used by the internal minimizer.
    /// Typical values: "legacy" (bounded coordinate search), "projected-gradient", "custom".
    String minimizerStrategy;

    /// Whether the internal result was compared against the legacy bounded search.
    bool comparedAgainstLegacy = false;

    /// Whether the evaluation fell back from the projected-gradient path to the legacy path.
    bool fallbackToLegacy = false;

    /// Whether the projected-gradient result was accepted (within agreement tolerance of legacy).
    bool projectedGradientAccepted = false;

    /// Number of projected-gradient iterations (set when the projected-gradient path ran).
    Index projectedGradientIterations = 0;

    /// Number of legacy minimizer iterations (set when the legacy path ran or was compared).
    Index legacyIterations = 0;
};

/// Model type for an internally minimized single solid solution.
///
/// Given an InternalSolidSolutionInput, returns the reduced thermodynamics and
/// diagnostics. This interface is separate from the MAGEMin pilot wrapper and from
/// the globalized seam, so it can be reused with different databases and families.
using InternalSolidSolutionModel = Fn<InternalSolidSolutionResult(InternalSolidSolutionInput const&)>;

//=============================================================================
// Candidate-generation interface for multi-phase assembly
//=============================================================================

/// One proposed candidate state to be explored before branch-local refinement.
///
/// In a multi-candidate phase assembly the outer Reaktoro problem contains
/// one phase instance per candidate. Each candidate is assigned a branch and
/// an optional seed composition that biases the first evaluation toward a
/// specific region of the internal composition space.
struct SolidSolutionCandidateState
{
    /// Index of the branch this candidate belongs to.
    /// The sentinel value static_cast<Index>(-1) means the candidate is not branch-constrained.
    Index branch = static_cast<Index>(-1);

    /// Initial internal composition seed for this candidate.
    /// When non-empty, passed as the warm-start on the first evaluation.
    ArrayXr seedx;

    /// Relative priority for ordering candidates; lower values are preferred.
    real priority = 0.0;

    /// Human-readable label appended to the candidate phase name in multi-phase assembly.
    /// When empty, the branch label/id or branch index is used instead.
    String label;
};

/// Generator that proposes candidate states at a reference thermodynamic condition.
///
/// Used by multi-candidate phase assembly to enumerate competing internal states
/// (different branches, endmember corners, or immiscible compositions) before
/// the outer Reaktoro solver refines the assembled equilibrium. Each call should
/// return a small ordered set of candidates; the assembly creates one Reaktoro
/// phase per candidate so the outer minimization can choose among them.
using SolidSolutionCandidateGenerator = Fn<Vec<SolidSolutionCandidateState>(
    real T, real P, ArrayXrConstRef referencex)>;

//=============================================================================
// Telemetry extraction and benchmarking
//=============================================================================

/// Telemetry extracted from one activity model evaluation of a MAGEMin pilot model.
///
/// Populated by `ExtractSolidSolutionMinimizerTelemetry` from the extra map that
/// the MAGEMin pilot activity model populates after each evaluation.
struct SolidSolutionMinimizerTelemetry
{
    /// Selected minimizer strategy for this evaluation.
    /// Typical values: "legacy", "projected-gradient", "custom", "".
    String strategy;

    /// Whether the evaluation compared projected-gradient against the legacy path.
    bool comparedAgainstLegacy = false;

    /// Whether the evaluation fell back to the legacy minimizer.
    bool fallbackToLegacy = false;

    /// Whether the projected-gradient result was accepted.
    bool projectedGradientAccepted = false;

    /// Number of projected-gradient iterations for this evaluation.
    Index projectedGradientIterations = 0;

    /// Number of legacy minimizer iterations for this evaluation.
    Index legacyIterations = 0;

    /// Max-composition delta between projected-gradient and legacy (only set when compared).
    double compositionDelta = 0.0;

    /// Objective delta between projected-gradient and legacy (only set when compared).
    double objectiveDelta = 0.0;

    /// Whether the evaluation reused a cached state from a prior evaluation.
        /// Whether projected-gradient had a lower objective than legacy (only set when compared and disagreed).
        bool pgHasLowerObjective = false;

        /// Whether the evaluation reused a cached state from a prior evaluation.
        bool usedStateCache = false;

    /// Number of internal minimizer iterations as emitted in InternalMinimizerIterations.
    Index internalMinimizerIterations = 0;

    /// Whether the internal minimizer converged.
    bool internalMinimizerConverged = false;
};

/// Aggregated benchmark statistics from accumulated MAGEMin pilot telemetry.
struct SolidSolutionMinimizerBenchmarkStats
{
    /// Total number of accumulated evaluations.
    Index totalEvaluations = 0;

    /// Number of evaluations where projected-gradient was selected.
    Index projectedGradientCount = 0;

    /// Number of evaluations where the legacy minimizer was selected.
    Index legacyCount = 0;

    /// Number of evaluations where a custom minimizer was used.
    Index customCount = 0;

    /// Number of evaluations where projected-gradient was compared against legacy.
    Index comparedCount = 0;

    /// Number of evaluations where the path fell back from projected-gradient to legacy.
    Index fallbackCount = 0;

    /// Number of evaluations that reused a cached internal state.
    Index stateCacheHits = 0;

    /// Average projected-gradient iteration count over evaluations that ran projected-gradient.
    real averageProjectedGradientIterations = 0.0;

    /// Average legacy minimizer iteration count over evaluations that ran the legacy solver.
    real averageLegacyIterations = 0.0;

    /// Fraction of compared evaluations that fell back to legacy (fallbackCount / comparedCount).
    real fallbackRate = 0.0;

    /// Fraction of total evaluations where projected-gradient was selected.
    real projectedGradientSelectionRate = 0.0;

        /// When PG and legacy disagree: number of fallback cases where PG had the lower objective.
        Index pgLowerObjectiveCount = 0;

        /// When PG and legacy disagree: number of fallback cases where legacy had the lower objective.
        Index legacyLowerObjectiveCount = 0;
};

/// Extract MAGEMin pilot minimizer telemetry from an activity model extra map.
/// Returns a default-constructed telemetry struct if the expected keys are absent.
inline auto ExtractSolidSolutionMinimizerTelemetry(
    Map<String, Any> const& extra) -> SolidSolutionMinimizerTelemetry
{
    SolidSolutionMinimizerTelemetry t;

    // Helper: try to read a typed value from the map.
    auto tryGet = [&](String const& key, auto& target)
    {
        using T = std::decay_t<decltype(target)>;
        auto it = extra.find(key);
        if(it != extra.end())
            if(auto p = std::any_cast<T>(&it->second))
                target = *p;
    };

    tryGet("MAGEMinSolidSolutionPilot::SelectedMinimizerStrategy",   t.strategy);
    tryGet("MAGEMinSolidSolutionPilot::ComparedAgainstLegacy",       t.comparedAgainstLegacy);
    tryGet("MAGEMinSolidSolutionPilot::FallbackToLegacy",            t.fallbackToLegacy);
    tryGet("MAGEMinSolidSolutionPilot::ProjectedGradientAccepted",   t.projectedGradientAccepted);
    tryGet("MAGEMinSolidSolutionPilot::UsedStateCache",              t.usedStateCache);
    tryGet("MAGEMinSolidSolutionPilot::InternalMinimizerConverged",  t.internalMinimizerConverged);
    tryGet("MAGEMinSolidSolutionPilot::ProjectedGradientLegacyCompositionDelta", t.compositionDelta);
    tryGet("MAGEMinSolidSolutionPilot::ProjectedGradientLegacyObjectiveDelta",   t.objectiveDelta);
    tryGet("MAGEMinSolidSolutionPilot::ProjectedGradientHasLowerObjective",       t.pgHasLowerObjective);

    // uint64_t keys (emitted as std::uint64_t by the pilot).
    auto tryGetU64 = [&](String const& key, Index& target)
    {
        auto it = extra.find(key);
        if(it != extra.end())
            if(auto p = std::any_cast<std::uint64_t>(&it->second))
                target = static_cast<Index>(*p);
    };

    tryGetU64("MAGEMinSolidSolutionPilot::ProjectedGradientIterationCount", t.projectedGradientIterations);
    tryGetU64("MAGEMinSolidSolutionPilot::LegacyMinimizerIterationCount",   t.legacyIterations);
    tryGetU64("MAGEMinSolidSolutionPilot::InternalMinimizerIterations",     t.internalMinimizerIterations);

    return t;
}

/// Benchmark accumulator for MAGEMin pilot minimizer telemetry.
///
/// Call `accumulate` after each activity model evaluation to track statistics
/// over many evaluations. Use `stats()` to retrieve the aggregated results and
/// `reset()` to clear the accumulator for a new measurement window.
class SolidSolutionMinimizerBenchmark
{
public:
    /// Accumulate telemetry from one activity model evaluation via its extra map.
    void accumulate(Map<String, Any> const& extra)
    {
        accumulate(ExtractSolidSolutionMinimizerTelemetry(extra));
    }

    /// Accumulate a pre-extracted telemetry struct.
    void accumulate(SolidSolutionMinimizerTelemetry const& t)
    {
        totalEvaluations_ += 1;

        if(t.strategy == "projected-gradient")
            projectedGradientCount_ += 1;
        else if(t.strategy == "legacy")
            legacyCount_ += 1;
        else if(!t.strategy.empty())
            customCount_ += 1;

        if(t.comparedAgainstLegacy)
        {
            comparedCount_ += 1;
            if(t.fallbackToLegacy)
                fallbackCount_ += 1;
                if(t.fallbackToLegacy)
                {
                    if(t.pgHasLowerObjective)
                        pgLowerObjectiveCount_ += 1;
                    else
                        legacyLowerObjectiveCount_ += 1;
                }
        }

        if(t.usedStateCache)
            stateCacheHits_ += 1;

        if(t.projectedGradientIterations > 0)
        {
            sumProjectedGradientIterations_ += t.projectedGradientIterations;
            projectedGradientIterationSamples_ += 1;
        }

        if(t.legacyIterations > 0)
        {
            sumLegacyIterations_ += t.legacyIterations;
            legacyIterationSamples_ += 1;
        }
    }

    /// Reset all accumulated statistics.
    void reset()
    {
        totalEvaluations_ = 0;
        projectedGradientCount_ = 0;
        legacyCount_ = 0;
        customCount_ = 0;
        comparedCount_ = 0;
        fallbackCount_ = 0;
        stateCacheHits_ = 0;
        sumProjectedGradientIterations_ = 0;
        projectedGradientIterationSamples_ = 0;
        sumLegacyIterations_ = 0;
        legacyIterationSamples_ = 0;
            pgLowerObjectiveCount_ = 0;
            legacyLowerObjectiveCount_ = 0;
    }

    /// Return current aggregated benchmark statistics.
    auto stats() const -> SolidSolutionMinimizerBenchmarkStats
    {
        SolidSolutionMinimizerBenchmarkStats s;
        s.totalEvaluations             = totalEvaluations_;
        s.projectedGradientCount       = projectedGradientCount_;
        s.legacyCount                  = legacyCount_;
        s.customCount                  = customCount_;
        s.comparedCount                = comparedCount_;
        s.fallbackCount                = fallbackCount_;
        s.stateCacheHits               = stateCacheHits_;

        s.averageProjectedGradientIterations = (projectedGradientIterationSamples_ > 0)
            ? real(sumProjectedGradientIterations_) / real(projectedGradientIterationSamples_)
            : real(0.0);

        s.averageLegacyIterations = (legacyIterationSamples_ > 0)
            ? real(sumLegacyIterations_) / real(legacyIterationSamples_)
            : real(0.0);

        s.fallbackRate = (comparedCount_ > 0)
            ? real(fallbackCount_) / real(comparedCount_)
            : real(0.0);

        s.projectedGradientSelectionRate = (totalEvaluations_ > 0)
            ? real(projectedGradientCount_) / real(totalEvaluations_)
            : real(0.0);
            s.pgLowerObjectiveCount     = pgLowerObjectiveCount_;
            s.legacyLowerObjectiveCount = legacyLowerObjectiveCount_;

        return s;
    }

private:
    Index totalEvaluations_                  = 0;
    Index projectedGradientCount_            = 0;
    Index legacyCount_                       = 0;
    Index customCount_                       = 0;
    Index comparedCount_                     = 0;
    Index fallbackCount_                     = 0;
    Index stateCacheHits_                    = 0;
    Index sumProjectedGradientIterations_    = 0;
    Index projectedGradientIterationSamples_ = 0;
    Index sumLegacyIterations_               = 0;
    Index legacyIterationSamples_            = 0;
        Index pgLowerObjectiveCount_             = 0;
        Index legacyLowerObjectiveCount_         = 0;
};

} // namespace Reaktoro
