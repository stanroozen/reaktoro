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
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelPerplexDEW.hpp>
using namespace Reaktoro;

void exportStandardThermoModelPerplexDEW(py::module& m)
{
    py::class_<StandardThermoModelParamsPerplexDEW>(m, "StandardThermoModelParamsPerplexDEW")
        .def(py::init<>())
        .def_readwrite("Gf", &StandardThermoModelParamsPerplexDEW::Gf)
        .def_readwrite("Hf", &StandardThermoModelParamsPerplexDEW::Hf)
        .def_readwrite("Sr", &StandardThermoModelParamsPerplexDEW::Sr)
        .def_readwrite("a1", &StandardThermoModelParamsPerplexDEW::a1)
        .def_readwrite("a2", &StandardThermoModelParamsPerplexDEW::a2)
        .def_readwrite("a3", &StandardThermoModelParamsPerplexDEW::a3)
        .def_readwrite("a4", &StandardThermoModelParamsPerplexDEW::a4)
        .def_readwrite("c1", &StandardThermoModelParamsPerplexDEW::c1)
        .def_readwrite("c2", &StandardThermoModelParamsPerplexDEW::c2)
        .def_readwrite("wref", &StandardThermoModelParamsPerplexDEW::wref)
        .def_readwrite("charge", &StandardThermoModelParamsPerplexDEW::charge)
        .def_readwrite("Tmax", &StandardThermoModelParamsPerplexDEW::Tmax)
        ;

    m.def("StandardThermoModelPerplexDEW", &StandardThermoModelPerplexDEW);
}
