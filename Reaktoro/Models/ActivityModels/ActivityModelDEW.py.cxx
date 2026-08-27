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
#include <Reaktoro/Models/ActivityModels/ActivityModelDEW.hpp>
using namespace Reaktoro;

void exportActivityModelDEW(py::module& m)
{
    py::enum_<ActivityDHModel>(m, "ActivityDHModel")
        .value("ExtendedDH", ActivityDHModel::ExtendedDH)
        .value("Davies",     ActivityDHModel::Davies)
        .export_values();

    py::class_<ActivityModelParamsDEW>(m, "ActivityModelParamsDEW")
        .def(py::init<>())
        .def_readwrite("dhModel",      &ActivityModelParamsDEW::dhModel)
        .def_readwrite("waterOptions", &ActivityModelParamsDEW::waterOptions)
        .def_readwrite("bExtended",    &ActivityModelParamsDEW::bExtended)
        ;

    m.def("ActivityModelDEW", py::overload_cast<>(&ActivityModelDEW),
        "Return the activity model for aqueous electrolyte phases based on the DEW (Deep Earth Water) model.\n"
        "\n"
        "This model implements the Helgeson–Kirkham–Flowers extended Debye–Hückel formulation\n"
        "with electrostatic coefficients (A and B parameters) evaluated dynamically from water\n"
        "density and dielectric constant using the DEW equation of state, valid up to 1000°C\n"
        "and 60 kbar.\n"
        "\n"
        "Returns\n"
        "-------\n"
        "ActivityModelGenerator\n"
        "    A generator function that creates an ActivityModel for aqueous phases."
    );

    m.def("ActivityModelDEW", py::overload_cast<const ActivityModelParamsDEW&>(&ActivityModelDEW),
        py::arg("params"),
        "Return the activity model for aqueous electrolyte phases based on the DEW (Deep Earth Water) model using explicit parameters.");
}


