"""Computational smoke test against the real Reaktoro forward model and the
actual DEW17HP622_Zn_2025 database. This is an integration check that the
full thermoinvert stack (forward model, database overlay, numerical
sensitivity) runs correctly against genuine Reaktoro equilibrium solves; it
is NOT a scientific calibration or a claim about literature-consistent
enthalpy values.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))

from thermoinvert.data import Experiment, ThermodynamicParameter
from thermoinvert.forward_model import ReaktoroForwardModel
from thermoinvert.sensitivity import identify_sensitivity, numerical_sensitivity

DB_PATH = (
    Path(__file__).parents[1]
    / "embedded/databases/perplex/DEW17HP622_Zn_2025-reaktoro.json"
)

# A restricted mineral set keeps each real equilibrium solve fast; this is a
# subset of the full Vazante mineral assemblage used elsewhere in the repo.
AQUEOUS_SPECIES = [
    "H2O",
    "H+",
    "OH-",
    "Al+3",
    "Ca+2",
    "Mg+2",
    "Fe+2",
    "K+",
    "Mn+2",
    "Zn2+",
    "SiO2,aq",
    "HS-",
    "HCO3-",
    "CO3-2",
]
MINERAL_NAMES = ["cc", "mag", "Znc"]

VAZANTE_BULK_COMPOSITION = {
    "Al": 0.02275,
    "Ca": 0.50823,
    "Fe": 0.03507,
    "K": 0.00155,
    "Mg": 0.48879,
    "Mn": 0.001917,
    "Si": 0.03246,
    "Zn": 0.02286,
    "S": 0.000312,
    "C": 1.0407,
}


def _vazante_experiment() -> Experiment:
    return Experiment(
        id="vazante-500C-2000bar",
        pressure_bar=2000.0,
        temperature_c=500.0,
        bulk_composition=VAZANTE_BULK_COMPOSITION,
        observed_phases=("cc", "mag"),
        metadata={"water_moles": 55.51},
    )


def test_real_forward_model_solves_vazante_bulk_composition():
    forward_model = ReaktoroForwardModel(DB_PATH, AQUEOUS_SPECIES, MINERAL_NAMES)
    prediction = forward_model.predict(_vazante_experiment())
    assert prediction.valid
    assert prediction.stable_phases


def test_numerical_sensitivity_runs_against_real_equilibrium_solves():
    forward_model = ReaktoroForwardModel(DB_PATH, AQUEOUS_SPECIES, MINERAL_NAMES)
    experiment = _vazante_experiment()
    parameters = [
        ThermodynamicParameter("dH_cc", "cc", prior_sigma_j_mol=2000.0),
    ]
    result = numerical_sensitivity(forward_model, [experiment], parameters)
    assert result.jacobian.shape[1] == 1
    assert result.jacobian.shape[0] == len(result.observable_names)
    identifiability = identify_sensitivity(result.jacobian, ["dH_cc"])
    assert identifiability.rank in (0, 1)
