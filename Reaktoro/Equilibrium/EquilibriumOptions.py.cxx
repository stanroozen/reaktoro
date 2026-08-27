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
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
using namespace Reaktoro;

void exportEquilibriumOptions(py::module& m)
{
    // Export GibbsHessian enum
    py::enum_<GibbsHessian>(m, "GibbsHessian")
        .value("Exact", GibbsHessian::Exact, "The Hessian of the Gibbs energy function is fully exact.")
        .value("PartiallyExact", GibbsHessian::PartiallyExact, "The Hessian of the Gibbs energy function is partially exact, partially approximated.")
        .value("Approx", GibbsHessian::Approx, "The Hessian of the Gibbs energy function is approximated using ideal thermodynamic models.")
        .value("ApproxDiagonal", GibbsHessian::ApproxDiagonal, "The Hessian of the Gibbs energy function is a diagonal matrix approximation using ideal thermodynamic models.")
        ;

    // Export EquilibriumOptions class
    py::class_<EquilibriumOptions>(m, "EquilibriumOptions")
        .def(py::init<>())
        .def_readwrite("optima", &EquilibriumOptions::optima)
        .def_readwrite("epsilon", &EquilibriumOptions::epsilon)
        .def_readwrite("logarithm_barrier_factor", &EquilibriumOptions::logarithm_barrier_factor)
        .def_readwrite("warmstart", &EquilibriumOptions::warmstart)
        .def_readwrite("use_ideal_activity_models", &EquilibriumOptions::use_ideal_activity_models)
        .def_readwrite("hessian", &EquilibriumOptions::hessian)
        ;
}
