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
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

class Database;

/// One stoichiometric coefficient term linking target element and species.
struct ElementStoichiometryTerm
{
    String element;
    String species;
    double coefficient = 0.0;
};

/// Result type for interpolation-based residual matching.
struct InterpolationResidualResult
{
    ArrayXd calculated;
    ArrayXd residuals;
};

/// Quantile band result over Monte Carlo samples.
struct UncertaintyBandResult
{
    ArrayXd lower;
    ArrayXd median;
    ArrayXd upper;
};

/// Apply per-entity Gibbs energy shifts (J/mol) to a mineral JSON database text.
auto perturbMineralDatabaseJSON(
    String const& base_json,
    Strings const& entities,
    ArrayXdConstRef const& shifts_j_per_mol
) -> String;

/// Apply per-entity Gibbs energy shifts and return sampled database in memory.
auto perturbMineralDatabase(
    String const& base_json,
    Strings const& entities,
    ArrayXdConstRef const& shifts_j_per_mol
) -> Database;

/// Apply per-sample Gibbs shifts (rows=samples, cols=entities) and return all sampled databases.
auto perturbMineralDatabases(
    String const& base_json,
    Strings const& entities,
    MatrixXdConstRef const& shifts_j_per_mol_samples,
    Index num_threads=1
) -> Vec<Database>;

/// Return aqueous species names whose formula elements are all in `allowed_elements`.
auto aqueousSpeciesNamesWithAllowedElements(
    Database const& database,
    Strings const& allowed_elements,
    Strings const& excluded_species={}
) -> Strings;

/// Return stoichiometry terms (element, species, coefficient) for selected species/target elements.
auto elementStoichiometryTerms(
    Database const& database,
    Strings const& species_names,
    Strings const& target_elements
) -> Vec<ElementStoichiometryTerm>;

/// Return an interpolated curve value using finite in-range points only.
auto interpolateCurveValue(
    ArrayXdConstRef const& x,
    ArrayXdConstRef const& y,
    double x_query,
    double atol=1e-8
) -> double;

/// Compute residuals for query points using interpolation of a reference curve.
auto computeResidualsInterpolated(
    ArrayXdConstRef const& curve_x,
    ArrayXdConstRef const& curve_y,
    ArrayXdConstRef const& query_x,
    ArrayXdConstRef const& query_y,
    double atol=1e-8
) -> InterpolationResidualResult;

/// Compute central confidence band (lower/median/upper) per column over sample matrix.
auto computeUncertaintyBand(
    ArrayXXdConstRef const& samples,
    double ci_percent=95.0
) -> UncertaintyBandResult;

} // namespace Reaktoro
