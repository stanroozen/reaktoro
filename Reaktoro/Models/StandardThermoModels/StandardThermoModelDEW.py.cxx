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
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelDEW.hpp>
using namespace Reaktoro;

void exportStandardThermoModelDEW(py::module& m)
{
    py::class_<StandardThermoModelParamsDEW>(m, "StandardThermoModelParamsDEW")
        .def(py::init<>())
        .def_readwrite("Gf",            &StandardThermoModelParamsDEW::Gf)
        .def_readwrite("Hf",            &StandardThermoModelParamsDEW::Hf)
        .def_readwrite("Sr",            &StandardThermoModelParamsDEW::Sr)
        .def_readwrite("a1",            &StandardThermoModelParamsDEW::a1)
        .def_readwrite("a2",            &StandardThermoModelParamsDEW::a2)
        .def_readwrite("a3",            &StandardThermoModelParamsDEW::a3)
        .def_readwrite("a4",            &StandardThermoModelParamsDEW::a4)
        .def_readwrite("c1",            &StandardThermoModelParamsDEW::c1)
        .def_readwrite("c2",            &StandardThermoModelParamsDEW::c2)
        .def_readwrite("wref",          &StandardThermoModelParamsDEW::wref)
        .def_readwrite("charge",        &StandardThermoModelParamsDEW::charge)
        .def_readwrite("Tmax",          &StandardThermoModelParamsDEW::Tmax)
        .def_readwrite("waterOptions",  &StandardThermoModelParamsDEW::waterOptions)
        ;

    m.def("StandardThermoModelDEW", StandardThermoModelDEW);
}
