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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSweepSolver.hpp>
using namespace Reaktoro;

void exportEquilibriumSweepSolver(py::module& m)
{
    py::class_<EquilibriumSweepOptions>(m, "EquilibriumSweepOptions")
        .def(py::init<>())
        .def_readwrite("num_threads", &EquilibriumSweepOptions::num_threads, "Number of worker threads. If <= 1, points are solved serially.")
        .def_readwrite("continuation", &EquilibriumSweepOptions::continuation, "Reuse previous converged state as initial guess in serial sweeps.")
        ;

    py::class_<EquilibriumSweepResult>(m, "EquilibriumSweepResult")
        .def(py::init<>())
        .def("size", &EquilibriumSweepResult::size, "Return number of sweep points.")
        .def("succeededCount", &EquilibriumSweepResult::succeededCount, "Return number of successful sweep points.")
        .def("elementMolalityArray", &EquilibriumSweepResult::elementMolalityArray,
            py::arg("element"),
            "Return total dissolved element molality (mol/kg water) for each sweep point.")
        .def_readwrite("states", &EquilibriumSweepResult::states)
        .def_readwrite("results", &EquilibriumSweepResult::results)
        ;

    py::class_<EquilibriumSweepGridResult>(m, "EquilibriumSweepGridResult")
        .def(py::init<>())
        .def("sizeX", &EquilibriumSweepGridResult::sizeX, "Return number of x-grid points.")
        .def("sizeY", &EquilibriumSweepGridResult::sizeY, "Return number of y-grid points.")
        .def("size", &EquilibriumSweepGridResult::size, "Return total number of grid points.")
        .def("succeededCount", &EquilibriumSweepGridResult::succeededCount, "Return number of successful grid points.")
        .def("logActivityGrid", &EquilibriumSweepGridResult::logActivityGrid,
            py::arg("species"),
            "Return 2D log10 activity grid for a species.")
        .def("predominantSpeciesGrid", &EquilibriumSweepGridResult::predominantSpeciesGrid,
            py::arg("species"),
            "Return 2D predominance grid with 0-based indices into the provided species list.")
        .def("saturationIndexGrid", &EquilibriumSweepGridResult::saturationIndexGrid,
            py::arg("species"),
            "Return 2D saturation index grid for a non-aqueous species.")
        .def("elementMolalityGrid", &EquilibriumSweepGridResult::elementMolalityGrid,
            py::arg("element"),
            "Return 2D total dissolved element molality grid (mol/kg water).")
        .def_readwrite("xvalues", &EquilibriumSweepGridResult::xvalues)
        .def_readwrite("yvalues", &EquilibriumSweepGridResult::yvalues)
        .def_readwrite("states", &EquilibriumSweepGridResult::states)
        .def_readwrite("results", &EquilibriumSweepGridResult::results)
        ;

    py::class_<EquilibriumSweepSolver>(m, "EquilibriumSweepSolver")
        .def(py::init<ChemicalSystem const&>())
        .def(py::init<EquilibriumSpecs const&>())

        .def("setOptions", &EquilibriumSweepSolver::setOptions)
        .def("setSweepOptions", &EquilibriumSweepSolver::setSweepOptions)

        .def("sweepTP", &EquilibriumSweepSolver::sweepTP,
            py::arg("initial"),
            py::arg("temperatures"),
            py::arg("pressures"),
            py::arg("temperature_unit") = "K",
            py::arg("pressure_unit") = "Pa")

        .def("sweepInput", &EquilibriumSweepSolver::sweepInput,
            py::arg("initial"),
            py::arg("input"),
            py::arg("values"))

        .def("sweepPH", &EquilibriumSweepSolver::sweepPH,
            py::arg("initial"),
            py::arg("values"))

        .def("sweepEh", &EquilibriumSweepSolver::sweepEh,
            py::arg("initial"),
            py::arg("values"),
            py::arg("unit") = "V")

        .def("sweepPHEhGrid", &EquilibriumSweepSolver::sweepPHEhGrid,
            py::arg("initial"),
            py::arg("pH_values"),
            py::arg("Eh_values"),
            py::arg("Eh_unit") = "V")

        .def("sweepLgActivityGrid", &EquilibriumSweepSolver::sweepLgActivityGrid,
            py::arg("initial"),
            py::arg("speciesX"),
            py::arg("lgaX_values"),
            py::arg("speciesY"),
            py::arg("lgaY_values"),
            "Sweep log10 activity of two species over a 2D grid. "
            "EquilibriumSpecs must have lgActivity(speciesX) and lgActivity(speciesY) declared.")

        .def("sweepTPGrid", &EquilibriumSweepSolver::sweepTPGrid,
            py::arg("initial"),
            py::arg("temperature_values"),
            py::arg("temperature_unit") = "K",
            py::arg("pressure_values"),
            py::arg("pressure_unit") = "Pa",
            "Sweep temperature and pressure over a 2D grid.")

        .def("sweepLogfO2pHGrid", &EquilibriumSweepSolver::sweepLogfO2pHGrid,
            py::arg("initial"),
            py::arg("logfO2_values"),
            py::arg("fug_unit") = "bar",
            py::arg("pH_values"),
            "Sweep log10 O2 fugacity and pH over a 2D grid. "
            "EquilibriumSpecs must have fugacity(\"O2\") and pH declared.")
        ;
}
