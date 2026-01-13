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
#include <Reaktoro/Extensions/DEW/WaterModelOptions.hpp>
using namespace Reaktoro;

void exportDEWWaterModels(py::module& m)
{
    // Export WaterEosModel enum
    py::enum_<WaterEosModel>(m, "WaterEosModel")
        .value("WagnerPruss", WaterEosModel::WagnerPruss)
        .value("HGK", WaterEosModel::HGK)
        .value("ZhangDuan2005", WaterEosModel::ZhangDuan2005)
        .value("ZhangDuan2009", WaterEosModel::ZhangDuan2009)
        ;

    // Export WaterDielectricModel enum
    py::enum_<WaterDielectricModel>(m, "WaterDielectricModel")
        .value("JohnsonNorton1991", WaterDielectricModel::JohnsonNorton1991)
        .value("Franck1990", WaterDielectricModel::Franck1990)
        .value("Fernandez1997", WaterDielectricModel::Fernandez1997)
        .value("PowerFunction", WaterDielectricModel::PowerFunction)
        ;

    // Export WaterGibbsModel enum
    py::enum_<WaterGibbsModel>(m, "WaterGibbsModel")
        .value("DelaneyHelgeson1978", WaterGibbsModel::DelaneyHelgeson1978)
        .value("DewIntegral", WaterGibbsModel::DewIntegral)
        ;

    // Export WaterBornModel enum
    py::enum_<WaterBornModel>(m, "WaterBornModel")
        .value("None", WaterBornModel::None)
        .value("Shock92Dew", WaterBornModel::Shock92Dew)
        ;

    // Export WaterModelOptions struct
    py::class_<WaterModelOptions>(m, "WaterModelOptions")
        .def(py::init<>())
        .def_readwrite("eosModel", &WaterModelOptions::eosModel)
        .def_readwrite("dielectricModel", &WaterModelOptions::dielectricModel)
        .def_readwrite("gibbsModel", &WaterModelOptions::gibbsModel)
        .def_readwrite("bornModel", &WaterModelOptions::bornModel)
        .def_readwrite("usePsatPolynomials", &WaterModelOptions::usePsatPolynomials)
        .def_readwrite("psatRelTol", &WaterModelOptions::psatRelTol)
        ;

    // Export makeWaterModelOptionsDEW function
    m.def("makeWaterModelOptionsDEW", &makeWaterModelOptionsDEW,
        "Create WaterModelOptions with DEW default settings");
}
