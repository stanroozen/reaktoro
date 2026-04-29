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

// Reaktoro includes
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelGlobalizedSolidSolution.hpp>

namespace Reaktoro {

/// Options for a branch-aware binary Redlich-Kister solid solution.
struct GlobalizedBinaryRedlichKisterOptions
{
    /// The Redlich-Kister parameter a0.
    real a0 = 0.0;

    /// The Redlich-Kister parameter a1.
    real a1 = 0.0;

    /// The Redlich-Kister parameter a2.
    real a2 = 0.0;

    /// The admissible composition branches of the solution.
    Vec<GlobalizedSolidSolutionBranch> branches;

    /// Tolerance used when checking if an external composition lies in a branch.
    real branchTolerance = 1.0e-12;

    /// Penalty added to disabled branches to keep them out of the outer minimization.
    real inactiveBranchPenalty = 1.0e6;

    /// Quadratic penalty that couples the minimized internal branch coordinate to the external composition.
    real externalCompositionPenalty = 25.0;

    /// Half-width of the score hysteresis used to prefer the persisted branch when branch scores are similar.
    real branchScoreHysteresis = 1.0e-8;

    /// Tolerance of the internal one-dimensional branch minimizer.
    real minimizerTolerance = 1.0e-10;

    /// Maximum number of iterations for the internal one-dimensional branch minimizer.
    Index minimizerMaxIterations = 128;

    /// Optional candidate generator used to screen branches before local refinement.
    GlobalizedSolidSolutionCandidateGenerator candidateGenerator;

    /// Optional branch-local stability screen evaluated after local refinement and before branch selection.
    GlobalizedSolidSolutionCandidateStabilityCriterion stabilityCriterion;
};

/// Construct a reduced binary Redlich-Kister solid-solution model.
auto GlobalizedBinaryRedlichKisterModel(GlobalizedBinaryRedlichKisterOptions options) -> GlobalizedSolidSolutionModel;

/// Construct an activity-model generator backed by the reduced binary Redlich-Kister model.
auto ActivityModelGlobalizedBinaryRedlichKister(GlobalizedBinaryRedlichKisterOptions options) -> ActivityModelGenerator;

/// Duplicate a phase into one phase per branch of a binary Redlich-Kister solid solution.
auto DuplicateGlobalizedBinaryRedlichKisterPhaseBranches(
    Phase const& phase,
    GlobalizedBinaryRedlichKisterOptions options,
    String suffixSeparator = "#") -> PhaseList;

/// Construct duplicated solid phases for a binary Redlich-Kister solid solution directly from a database.
auto GlobalizedBinaryRedlichKisterSolidPhases(
    Database const& db,
    String name,
    Strings const& species,
    GlobalizedBinaryRedlichKisterOptions options,
    String suffixSeparator = "#") -> PhaseList;

} // namespace Reaktoro