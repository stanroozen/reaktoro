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
#include <Reaktoro/Models/ActivityModels/ActivityModelPerplexGFSM.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXHybridEos.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXMrkMixture.hpp>
#include <Reaktoro/Extensions/Perple_X/PerpleXPureEos.hpp>
using namespace Reaktoro;
using namespace Reaktoro::PerpleX;

void exportActivityModelPerplexGFSM(py::module& m)
{
    // Some build configurations invoke this exporter twice. Guard against
    // re-registration of classes/functions in the same Python module.
    static bool exported = false;
    if(exported)
        return;
    exported = true;

    // -----------------------------------------------------------------------
    // PerpleXPureEosOptions — convergence / iteration control for pure EOS
    // -----------------------------------------------------------------------
    // PerpleXPureEosOptions is also exposed by Perple_X.py.cxx. Guard the
    // registration so importing the extension does not fail on duplicate types.
    if(!py::hasattr(m, "PerpleXPureEosOptions")) {
        py::class_<PerpleXPureEosOptions>(m, "PerpleXPureEosOptions")
            .def(py::init<>())
            .def_readwrite("maxIter", &PerpleXPureEosOptions::maxIter,
                "Maximum iterations for root-finding in pure EOS.")
            .def_readwrite("maxWarn", &PerpleXPureEosOptions::maxWarn,
                "Maximum warnings printed during pure EOS computation.")
            .def_readwrite("tol", &PerpleXPureEosOptions::tol,
                "Convergence tolerance for nonlinear root-finding.")
            ;
    }

    // -----------------------------------------------------------------------
    // MrkMixOptions — MRK mixture averaging / floor options
    // -----------------------------------------------------------------------
    if(py::hasattr(m, "PerpleXMrkMixOptions")) {
        m.attr("MrkMixOptions") = m.attr("PerpleXMrkMixOptions");
    } else {
        py::class_<MrkMixOptions>(m, "MrkMixOptions")
            .def(py::init<>())
            .def_readwrite("iavg", &MrkMixOptions::iavg,
                "Interaction-parameter averaging scheme: 1=geometric (default), "
                "2=arithmetic, else=harmonic.")
            .def_readwrite("minY", &MrkMixOptions::minY,
                "Floor for mole fractions in log expressions (avoids log(0)).")
            ;
    }

    // -----------------------------------------------------------------------
    // HybridEosOptions — per-species pure EOS selection for H2O, CO2, CH4
    // -----------------------------------------------------------------------
    if(py::hasattr(m, "PerpleXHybridEosOptions")) {
        m.attr("HybridEosOptions") = m.attr("PerpleXHybridEosOptions");
        if(py::hasattr(m, "PerpleXWaterEos")) m.attr("WaterEos") = m.attr("PerpleXWaterEos");
        if(py::hasattr(m, "PerpleXCO2Eos")) m.attr("CO2Eos") = m.attr("PerpleXCO2Eos");
        if(py::hasattr(m, "PerpleXCH4Eos")) m.attr("CH4Eos") = m.attr("PerpleXCH4Eos");
    } else {
        py::class_<HybridEosOptions> cls_heo(m, "HybridEosOptions");

        py::enum_<HybridEosOptions::WaterEos>(cls_heo, "WaterEos")
            .value("Mrk",        HybridEosOptions::WaterEos::Mrk,      "Modified Redlich-Kwong (fast, default)")
            .value("Hsmrk",      HybridEosOptions::WaterEos::Hsmrk,    "Hard-sphere MRK (Kerrick & Jacobs)")
            .value("Cork",       HybridEosOptions::WaterEos::Cork,      "Holland & Powell CORK (high P/T)")
            .value("Pseos",      HybridEosOptions::WaterEos::Pseos,     "Pseff pseudo-Einstein")
            .value("Haar",       HybridEosOptions::WaterEos::Haar,      "Haar et al. 1979 (most accurate)")
            .value("ZhangDuan05", HybridEosOptions::WaterEos::ZhangDuan05, "Zhang & Duan 2005 correlation")
            .value("ZhangDuan09", HybridEosOptions::WaterEos::ZhangDuan09, "Zhang & Duan 2009 correlation")
            .export_values()
            ;

        py::enum_<HybridEosOptions::CO2Eos>(cls_heo, "CO2Eos")
            .value("Mrk",        HybridEosOptions::CO2Eos::Mrk,       "Modified Redlich-Kwong (default)")
            .value("Hsmrk",      HybridEosOptions::CO2Eos::Hsmrk,     "Hard-sphere MRK")
            .value("Cork",       HybridEosOptions::CO2Eos::Cork,       "Holland & Powell CORK")
            .value("Brmrk",      HybridEosOptions::CO2Eos::Brmrk,     "Bottinga & Richet MRK variant")
            .value("Pseos",      HybridEosOptions::CO2Eos::Pseos,      "Pseff model")
            .value("ZhangDuan09", HybridEosOptions::CO2Eos::ZhangDuan09, "Zhang & Duan 2009")
            .export_values()
            ;

        py::enum_<HybridEosOptions::CH4Eos>(cls_heo, "CH4Eos")
            .value("Mrk",        HybridEosOptions::CH4Eos::Mrk,       "Modified Redlich-Kwong (default)")
            .value("Hsmrk",      HybridEosOptions::CH4Eos::Hsmrk,     "Hard-sphere MRK")
            .value("ZhangDuan09", HybridEosOptions::CH4Eos::ZhangDuan09, "Zhang & Duan 2009")
            .export_values()
            ;

        cls_heo
            .def(py::init<>())
            .def_readwrite("water", &HybridEosOptions::water,
                "Pure EOS selection for H2O (HybridEosOptions.WaterEos).")
            .def_readwrite("co2",   &HybridEosOptions::co2,
                "Pure EOS selection for CO2 (HybridEosOptions.CO2Eos).")
            .def_readwrite("ch4",   &HybridEosOptions::ch4,
                "Pure EOS selection for CH4 (HybridEosOptions.CH4Eos).")
            ;
    }

    // Convenience factory: default COH-Fluid+ options (ZhangDuan09 for all three).
    m.def("makePerplexCOHFluidPlusEosOptions", &makePerplexCOHFluidPlusEosOptions,
        "Return HybridEosOptions pre-configured for Perple_X COH-Fluid+ "
        "(ZhangDuan09 for H2O, CO2, CH4 — matches iopt(25)=7, iopt(26)=7, iopt(27)=7).");

    // -----------------------------------------------------------------------
    // ActivityModelParamsPerplexGFSM
    // -----------------------------------------------------------------------
    py::class_<ActivityModelParamsPerplexGFSM>(m, "ActivityModelParamsPerplexGFSM")
        .def(py::init<>())
        .def_readwrite("hybridEosOptions",  &ActivityModelParamsPerplexGFSM::hybridEosOptions,
            "Per-species pure EOS selection for H2O, CO2, CH4 (HybridEosOptions).")
        .def_readwrite("mrkMixOptions",     &ActivityModelParamsPerplexGFSM::mrkMixOptions,
            "MRK mixture averaging/floor options (MrkMixOptions).")
        .def_readwrite("useLowTMrk",        &ActivityModelParamsPerplexGFSM::useLowTMrk,
            "Use low-temperature MRK variant.")
        .def_readwrite("enableElectrolyte", &ActivityModelParamsPerplexGFSM::enableElectrolyte,
            "Enable electrolyte solvent properties in GFSM evaluation.")
        .def_readwrite("pureEosOptions",    &ActivityModelParamsPerplexGFSM::pureEosOptions,
            "Convergence / iteration control for pure EOS solvers (PerpleXPureEosOptions).")
        ;

    m.def("ActivityModelPerplexGFSM", py::overload_cast<const ActivityModelParamsPerplexGFSM&>(&ActivityModelPerplexGFSM),
        py::arg("params") = ActivityModelParamsPerplexGFSM{},
        "Return the activity model for Perple_X GFSM (ifug=39) fluid phases.");
}