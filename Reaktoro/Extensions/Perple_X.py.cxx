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
#include <Reaktoro/Extensions/Perple_X/PerpleXHybridEos.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXMrkMixture.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXPureEos.hpp>
using namespace Reaktoro;

void exportExtensionPerpleX(py::module& m)
{
    py::enum_<PerpleX::HybridEosOptions::WaterEos>(m, "PerpleXWaterEos")
        .value("Mrk", PerpleX::HybridEosOptions::WaterEos::Mrk)
        .value("Hsmrk", PerpleX::HybridEosOptions::WaterEos::Hsmrk)
        .value("Cork", PerpleX::HybridEosOptions::WaterEos::Cork)
        .value("Pseos", PerpleX::HybridEosOptions::WaterEos::Pseos)
        .value("Haar", PerpleX::HybridEosOptions::WaterEos::Haar)
        .value("ZhangDuan05", PerpleX::HybridEosOptions::WaterEos::ZhangDuan05)
        .value("ZhangDuan09", PerpleX::HybridEosOptions::WaterEos::ZhangDuan09)
        ;

    py::enum_<PerpleX::HybridEosOptions::CO2Eos>(m, "PerpleXCO2Eos")
        .value("Mrk", PerpleX::HybridEosOptions::CO2Eos::Mrk)
        .value("Hsmrk", PerpleX::HybridEosOptions::CO2Eos::Hsmrk)
        .value("Cork", PerpleX::HybridEosOptions::CO2Eos::Cork)
        .value("Brmrk", PerpleX::HybridEosOptions::CO2Eos::Brmrk)
        .value("Pseos", PerpleX::HybridEosOptions::CO2Eos::Pseos)
        .value("ZhangDuan09", PerpleX::HybridEosOptions::CO2Eos::ZhangDuan09)
        ;

    py::enum_<PerpleX::HybridEosOptions::CH4Eos>(m, "PerpleXCH4Eos")
        .value("Mrk", PerpleX::HybridEosOptions::CH4Eos::Mrk)
        .value("Hsmrk", PerpleX::HybridEosOptions::CH4Eos::Hsmrk)
        .value("ZhangDuan09", PerpleX::HybridEosOptions::CH4Eos::ZhangDuan09)
        ;

    py::class_<PerpleX::PerpleXPureEosOptions>(m, "PerpleXPureEosOptions")
        .def(py::init<>())
        .def_readwrite("maxIter", &PerpleX::PerpleXPureEosOptions::maxIter)
        .def_readwrite("maxWarn", &PerpleX::PerpleXPureEosOptions::maxWarn)
        .def_readwrite("tol", &PerpleX::PerpleXPureEosOptions::tol)
        ;

    py::class_<PerpleX::MrkMixOptions>(m, "PerpleXMrkMixOptions")
        .def(py::init<>())
        .def_readwrite("iavg", &PerpleX::MrkMixOptions::iavg)
        .def_readwrite("minY", &PerpleX::MrkMixOptions::minY)
        ;

    py::class_<PerpleX::HybridEosOptions>(m, "PerpleXHybridEosOptions")
        .def(py::init([]() { return PerpleX::makePerpleXHybridEosOptions(); }))
        .def_readwrite("water", &PerpleX::HybridEosOptions::water)
        .def_readwrite("co2", &PerpleX::HybridEosOptions::co2)
        .def_readwrite("ch4", &PerpleX::HybridEosOptions::ch4)
        ;

    m.def("makePerpleXHybridEosOptions", &PerpleX::makePerpleXHybridEosOptions,
        py::arg("options") = PerpleX::PerpleXPureEosOptions{},
        "Create HybridEosOptions with Perple_X pure-EOS callbacks configured.");

    m.def("makePerplexCOHFluidPlusEosOptions", &PerpleX::makePerplexCOHFluidPlusEosOptions,
        py::arg("options") = PerpleX::PerpleXPureEosOptions{},
        "Create HybridEosOptions with ZhangDuan09 for H2O, CO2, and CH4, matching the Perple_X COH-Fluid+ configuration (iopt(25)=iopt(26)=iopt(27)=7).");
}