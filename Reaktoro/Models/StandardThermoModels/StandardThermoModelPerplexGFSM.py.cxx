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
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelPerplexGFSM.hpp>
using namespace Reaktoro;

void exportStandardThermoModelPerplexGFSM(py::module& m)
{
    py::class_<StandardThermoModelParamsPerplexGFSM>(m, "StandardThermoModelParamsPerplexGFSM")
        .def(py::init<>())
        .def_readwrite("speciesIndex", &StandardThermoModelParamsPerplexGFSM::speciesIndex)
        .def_readwrite("G0", &StandardThermoModelParamsPerplexGFSM::G0)
        .def_readwrite("H0", &StandardThermoModelParamsPerplexGFSM::H0)
        .def_readwrite("V0", &StandardThermoModelParamsPerplexGFSM::V0)
        .def_readwrite("hybridEosOptions", &StandardThermoModelParamsPerplexGFSM::hybridEosOptions)
        .def_readwrite("mrkMixOptions", &StandardThermoModelParamsPerplexGFSM::mrkMixOptions)
        .def_readwrite("pureEosOptions", &StandardThermoModelParamsPerplexGFSM::pureEosOptions)
        .def_readwrite("useLowTMrk", &StandardThermoModelParamsPerplexGFSM::useLowTMrk)
        .def_readwrite("Tmax", &StandardThermoModelParamsPerplexGFSM::Tmax)
        ;

    m.def("StandardThermoModelPerplexGFSM", &StandardThermoModelPerplexGFSM,
        py::arg("params"),
        "Return a standard thermodynamic model for a Perple_X GFSM fluid species.");
}