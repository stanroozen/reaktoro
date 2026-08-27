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
#include <Reaktoro/Models/ActivityModels/ActivityModelPerplexDEW.hpp>
using namespace Reaktoro;

void exportActivityModelPerplexDEW(py::module& m)
{
    // ActivityDHModel is registered by exportActivityModelDEW (called first).
    // No re-registration needed here.

    py::class_<ActivityModelParamsPerplexDEW>(m, "ActivityModelParamsPerplexDEW")
        .def(py::init<>())
        .def_readwrite("dhModel", &ActivityModelParamsPerplexDEW::dhModel)
        .def_readwrite("errorOnConflictingStandardState", &ActivityModelParamsPerplexDEW::errorOnConflictingStandardState)
        .def_readwrite("warnOnUnmappedGFSMCoupling", &ActivityModelParamsPerplexDEW::warnOnUnmappedGFSMCoupling)
        .def_readwrite("requireCoupledGFSMHandoff", &ActivityModelParamsPerplexDEW::requireCoupledGFSMHandoff)
        ;

    m.def("ActivityModelPerplexDEW",
        py::overload_cast<const ActivityModelParamsPerplexDEW&>(&ActivityModelPerplexDEW),
        py::arg("params"),
        "Return the activity model for aqueous electrolyte phases using "
        "Perple_X-linked HKF solvent/electrolyte internals.");

    m.def("ActivityModelPerplexDEW",
        [](ActivityDHModel model) { return ActivityModelPerplexDEW(model); },
        py::arg("model") = ActivityDHModel::Davies,
        "Return the activity model for aqueous electrolyte phases using "
        "Perple_X-linked HKF solvent/electrolyte internals.");
}
