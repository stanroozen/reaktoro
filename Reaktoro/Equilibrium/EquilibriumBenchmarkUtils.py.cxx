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

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Equilibrium/EquilibriumBenchmarkUtils.hpp>
using namespace Reaktoro;

void exportEquilibriumBenchmarkUtils(py::module& m)
{
    py::class_<ElementStoichiometryTerm>(m, "ElementStoichiometryTerm")
        .def(py::init<>())
        .def_readwrite("element", &ElementStoichiometryTerm::element)
        .def_readwrite("species", &ElementStoichiometryTerm::species)
        .def_readwrite("coefficient", &ElementStoichiometryTerm::coefficient)
        ;

    py::class_<InterpolationResidualResult>(m, "InterpolationResidualResult")
        .def(py::init<>())
        .def_readwrite("calculated", &InterpolationResidualResult::calculated)
        .def_readwrite("residuals", &InterpolationResidualResult::residuals)
        ;

    py::class_<UncertaintyBandResult>(m, "UncertaintyBandResult")
        .def(py::init<>())
        .def_readwrite("lower", &UncertaintyBandResult::lower)
        .def_readwrite("median", &UncertaintyBandResult::median)
        .def_readwrite("upper", &UncertaintyBandResult::upper)
        ;

    m.def("perturbMineralDatabaseJSON", &perturbMineralDatabaseJSON,
        py::arg("base_json"),
        py::arg("entities"),
        py::arg("shifts_j_per_mol"));

    m.def("perturbMineralDatabase", &perturbMineralDatabase,
        py::arg("base_json"),
        py::arg("entities"),
        py::arg("shifts_j_per_mol"));

    m.def("perturbMineralDatabases", &perturbMineralDatabases,
        py::arg("base_json"),
        py::arg("entities"),
        py::arg("shifts_j_per_mol_samples"),
        py::arg("num_threads") = 1);

    m.def("aqueousSpeciesNamesWithAllowedElements", &aqueousSpeciesNamesWithAllowedElements,
        py::arg("database"),
        py::arg("allowed_elements"),
        py::arg("excluded_species") = Strings{});

    m.def("elementStoichiometryTerms", &elementStoichiometryTerms,
        py::arg("database"),
        py::arg("species_names"),
        py::arg("target_elements"));

    m.def("interpolateCurveValue", &interpolateCurveValue,
        py::arg("x"),
        py::arg("y"),
        py::arg("x_query"),
        py::arg("atol") = 1e-8);

    m.def("computeResidualsInterpolated", &computeResidualsInterpolated,
        py::arg("curve_x"),
        py::arg("curve_y"),
        py::arg("query_x"),
        py::arg("query_y"),
        py::arg("atol") = 1e-8);

    m.def("computeUncertaintyBand", &computeUncertaintyBand,
        py::arg("samples"),
        py::arg("ci_percent") = 95.0);
}
