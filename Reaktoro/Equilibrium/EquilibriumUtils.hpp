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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelGlobalizedSolidSolution.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class EquilibriumConditions;
class EquilibriumRestrictions;
class EquilibriumSpecs;
struct EquilibriumOptions;
struct EquilibriumResult;

/// Options controlling split-triggered equilibrium retries for globalized solid-solution phases.
struct GlobalizedSolidSolutionEquilibriumRetryOptions
{
	/// Definitions describing which phase families can request duplication.
	Vec<GlobalizedSolidSolutionPhaseDefinition> definitions;

	/// Equilibrium options applied to each solve attempt.
	EquilibriumOptions equilibrium;

	/// Maximum number of rebuild-and-retry passes after the initial solve.
	Index maxRetries = 1;

	/// Enable pre-rebuild split acceptance gating based on solve effort and split diagnostics.
	bool enableSplitAcceptanceGate = false;

	/// Minimum number of iterations in the initial solve before split retry is considered.
	Index minIterationsForSplitRetry = 30;

	/// Minimum required objective-gap evidence for accepting a split request when available.
	real splitImprovementTolerance = 1.0e-8;

	/// Enable fallback solve passes for rebuilt systems when the normal rebuilt solve fails.
	bool enableRebuiltFallbackSolve = true;

	/// Run an ideal-activity-model solve first in rebuilt fallback passes.
	bool rebuiltFallbackUseIdealFirst = true;

	/// Maximum number of fallback solve passes attempted per rebuilt system.
	Index rebuiltFallbackMaxRetries = 1;
};

/// Result returned by split-triggered equilibrium retries.
struct GlobalizedSolidSolutionEquilibriumRetryResult
{
	/// The final rebuilt chemical system used by the returned state.
	ChemicalSystem system;

	/// The final chemical state after all retries.
	ChemicalState state;

	/// The accumulated equilibrium result over all solve attempts.
	EquilibriumResult result;

	/// Number of times the system was rebuilt and retried.
	Index numRebuilds = 0;

	/// Number of split-retry requests accepted by the gating policy.
	Index numAcceptedSplitRetries = 0;

	/// Number of split-retry requests rejected by the gating policy.
	Index numRejectedSplitRetries = 0;

	/// Number of rebuilt fallback solve attempts performed.
	Index numRebuiltFallbackAttempts = 0;
};

/// Perform a chemical equilibrium calculation on a given chemical state.
/// The calculation is performed with fixed temperature and pressure obtained
/// from the chemical state, and the chemical system is closed, so chemical
/// elements and electric charge are conserved.
///@{
auto equilibrate(ChemicalState& state) -> EquilibriumResult;
auto equilibrate(ChemicalState& state, const EquilibriumOptions& options) -> EquilibriumResult;
auto equilibrate(ChemicalState& state, const EquilibriumRestrictions& restrictions) -> EquilibriumResult;
auto equilibrate(ChemicalState& state, const EquilibriumRestrictions& restrictions, const EquilibriumOptions& options) -> EquilibriumResult;

auto equilibrate(ChemicalState& state, ArrayXdConstRef b0) -> EquilibriumResult;
auto equilibrate(ChemicalState& state, const EquilibriumOptions& options, ArrayXdConstRef b0) -> EquilibriumResult;
auto equilibrate(ChemicalState& state, const EquilibriumRestrictions& restrictions, ArrayXdConstRef b0) -> EquilibriumResult;
auto equilibrate(ChemicalState& state, const EquilibriumRestrictions& restrictions, const EquilibriumOptions& options, ArrayXdConstRef b0) -> EquilibriumResult;
///@}

/// Perform a closed-system equilibrium calculation with bounded split-trigger retries.
///
/// The helper keeps Reaktoro's usual fixed-temperature/fixed-pressure solve path, but if one of the
/// supplied globalized solid-solution phase families requests duplication through the phase-scoped
/// diagnostics stored in `ChemicalProps::extra`, the chemical system is rebuilt and the solve is retried.
auto equilibrateWithGlobalizedSolidSolutionSplits(
	ChemicalState const& initialState,
	GlobalizedSolidSolutionEquilibriumRetryOptions options = {}) -> GlobalizedSolidSolutionEquilibriumRetryResult;

/// Perform a closed-system equilibrium calculation with bounded split-trigger retries and explicit component totals.
auto equilibrateWithGlobalizedSolidSolutionSplits(
	ChemicalState const& initialState,
	ArrayXdConstRef b0,
	GlobalizedSolidSolutionEquilibriumRetryOptions options = {}) -> GlobalizedSolidSolutionEquilibriumRetryResult;

/// Perform a split-triggered equilibrium calculation while preserving a richer equilibrium specification template.
auto equilibrateWithGlobalizedSolidSolutionSplits(
	ChemicalState const& initialState,
	EquilibriumSpecs const& specs,
	EquilibriumConditions const& conditions,
	GlobalizedSolidSolutionEquilibriumRetryOptions options = {}) -> GlobalizedSolidSolutionEquilibriumRetryResult;

/// Perform a split-triggered equilibrium calculation while preserving equilibrium conditions and restrictions across rebuilt systems.
auto equilibrateWithGlobalizedSolidSolutionSplits(
	ChemicalState const& initialState,
	EquilibriumSpecs const& specs,
	EquilibriumConditions const& conditions,
	EquilibriumRestrictions const& restrictions,
	GlobalizedSolidSolutionEquilibriumRetryOptions options = {}) -> GlobalizedSolidSolutionEquilibriumRetryResult;

} // namespace Reaktoro
