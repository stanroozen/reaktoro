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
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {

class ChemicalSystem;

/// Options controlling sweep execution.
struct EquilibriumSweepOptions
{
    /// Number of worker threads.
    /// If <= 1, solve points serially with continuation support.
    Index num_threads = 1;

    /// Reuse previous converged state as initial guess in serial sweeps.
    bool continuation = true;
};

/// Result of a sweep operation.
struct EquilibriumSweepResult
{
    /// The solved chemical state per sweep point.
    std::vector<ChemicalState> states;

    /// The equilibrium solve result per sweep point.
    std::vector<EquilibriumResult> results;

    /// Return number of sweep points.
    auto size() const -> Index;

    /// Return number of successful sweep points.
    auto succeededCount() const -> Index;

    /// Return the total dissolved element molality array (mol/kg water) for each sweep point.
    auto elementMolalityArray(StringOrIndex const& element) const -> ArrayXd;
};

/// Result of a 2D sweep operation over x and y variables.
struct EquilibriumSweepGridResult
{
    /// Grid values for x-axis variable.
    ArrayXd xvalues;

    /// Grid values for y-axis variable.
    ArrayXd yvalues;

    /// The solved chemical state per grid point in row-major x/y order.
    std::vector<ChemicalState> states;

    /// The equilibrium solve result per grid point in row-major x/y order.
    std::vector<EquilibriumResult> results;

    /// Return number of x-grid points.
    auto sizeX() const -> Index;

    /// Return number of y-grid points.
    auto sizeY() const -> Index;

    /// Return total number of grid points.
    auto size() const -> Index;

    /// Return number of successful grid points.
    auto succeededCount() const -> Index;

    /// Return the log10 activity grid for a species.
    auto logActivityGrid(StringOrIndex const& species) const -> ArrayXXd;

    /// Return the predominance grid index for a list of species names.
    /// The returned values are the 0-based index into the input species list.
    auto predominantSpeciesGrid(Strings const& species) const -> ArrayXXd;

    /// Return the saturation index grid for a non-aqueous species.
    auto saturationIndexGrid(StringOrIndex const& species) const -> ArrayXXd;

    /// Return the total dissolved element molality grid (mol/kg water).
    auto elementMolalityGrid(StringOrIndex const& element) const -> ArrayXXd;
};

/// Used for batch equilibrium calculations over multiple input points.
class EquilibriumSweepSolver
{
public:
    /// Construct an EquilibriumSweepSolver object with given chemical system.
    explicit EquilibriumSweepSolver(ChemicalSystem const& system);

    /// Construct an EquilibriumSweepSolver object with given equilibrium specifications.
    explicit EquilibriumSweepSolver(EquilibriumSpecs const& specs);

    /// Set the underlying equilibrium options.
    auto setOptions(EquilibriumOptions const& options) -> void;

    /// Set sweep execution options.
    auto setSweepOptions(EquilibriumSweepOptions const& options) -> void;

    /// Get the underlying equilibrium options.
    auto options() const -> EquilibriumOptions const&;

    /// Get the sweep execution options.
    auto sweepOptions() const -> EquilibriumSweepOptions const&;

    /// Sweep temperature and pressure with one solve per point.
    auto sweepTP(
        ChemicalState const& initial,
        ArrayXdConstRef const& temperatures,
        ArrayXdConstRef const& pressures,
        String const& temperature_unit,
        String const& pressure_unit
    ) -> EquilibriumSweepResult;

    /// Sweep an arbitrary input variable from EquilibriumConditions input names.
    auto sweepInput(
        ChemicalState const& initial,
        String const& input,
        ArrayXdConstRef const& values
    ) -> EquilibriumSweepResult;

    /// Sweep pH values with one solve per point.
    auto sweepPH(
        ChemicalState const& initial,
        ArrayXdConstRef const& values
    ) -> EquilibriumSweepResult;

    /// Sweep Eh values with one solve per point.
    auto sweepEh(
        ChemicalState const& initial,
        ArrayXdConstRef const& values,
        String const& unit
    ) -> EquilibriumSweepResult;

    /// Sweep pH and Eh values over a 2D grid.
    auto sweepPHEhGrid(
        ChemicalState const& initial,
        ArrayXdConstRef const& pH_values,
        ArrayXdConstRef const& Eh_values,
        String const& Eh_unit
    ) -> EquilibriumSweepGridResult;

    /// Sweep log10 activity of two species over a 2D grid.
    /// The caller's EquilibriumSpecs must have lgActivity(speciesX) and
    /// lgActivity(speciesY) declared before constructing the solver.
    auto sweepLgActivityGrid(
        ChemicalState const& initial,
        String const& speciesX,
        ArrayXdConstRef const& lgaX_values,
        String const& speciesY,
        ArrayXdConstRef const& lgaY_values
    ) -> EquilibriumSweepGridResult;

    /// Sweep temperature and pressure over a 2D grid.
    auto sweepTPGrid(
        ChemicalState const& initial,
        ArrayXdConstRef const& temperature_values,
        String const& temperature_unit,
        ArrayXdConstRef const& pressure_values,
        String const& pressure_unit
    ) -> EquilibriumSweepGridResult;

    /// Sweep log10 fugacity of O2 and pH over a 2D grid.
    /// The caller's EquilibriumSpecs must have fugacity("O2") and pH declared.
    /// @param logfO2_values  log10 of O2 fugacity values (dimensionless log10 bar).
    /// @param fug_unit       Unit of the fugacity values passed to conditions.fugacity() (e.g. "bar").
    auto sweepLogfO2pHGrid(
        ChemicalState const& initial,
        ArrayXdConstRef const& logfO2_values,
        String const& fug_unit,
        ArrayXdConstRef const& pH_values
    ) -> EquilibriumSweepGridResult;

private:
    EquilibriumSpecs m_specs;
    EquilibriumOptions m_options;
    EquilibriumSweepOptions m_sweep_options;
};

} // namespace Reaktoro
