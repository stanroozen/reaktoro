#!/usr/bin/env python
"""Check quartz Gibbs free energy in SUPCRTBL database."""

import sys
import os
import autodiff

# Try to import Reaktoro; fall back to local extension module if not installed
try:
    from reaktoro import *
except ModuleNotFoundError:
    # Add local build path for reaktoro4py if available
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    ROOT_DIR = os.path.dirname(SCRIPT_DIR)
    PYD_DIR = os.path.join(ROOT_DIR, "build", "Reaktoro", "Release")
    if os.path.isdir(PYD_DIR) and PYD_DIR not in sys.path:
        sys.path.insert(0, PYD_DIR)
    try:
        from reaktoro4py import *

        print("Using local reaktoro4py extension from build.")
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "Could not import Reaktoro. Install the 'reaktoro' package or ensure reaktoro4py is on PYTHONPATH."
        ) from e

# Load SUPCRTBL database
print("\nLoading SUPCRTBL database...")
db = SupcrtDatabase("supcrtbl")

# Get quartz species
quartz = db.species("Quartz")
print(f"\nQuartz Species Information:")
print(f"  Name: {quartz.name()}")
print(f"  Formula: {quartz.formula()}")

# Check what methods are available
print(f"\nGetting standard thermodynamic properties...")

# Get the thermo model
thermo_model = quartz.standardThermoModel()
print(f"Thermo Model: {thermo_model}")

# Evaluate at reference state (25C, 1 bar)
T_ref = autodiff.real(298.15)  # 25C in K
P_ref = autodiff.real(1.0e5)  # 1 bar in Pa

props_ref = quartz.standardThermoProps(T_ref, P_ref)
print(f"\nStandard Thermodynamic Properties at 25Â°C, 1 bar:")
print(f"  G0 (Gibbs Free Energy): {props_ref.G0.val() / 1000:.2f} kJ/mol")
print(f"  H0 (Enthalpy): {props_ref.H0.val() / 1000:.2f} kJ/mol")

# Test at higher temperature/pressure
T_test = autodiff.real(373.15)  # 100C
P_test = autodiff.real(100.0e5)  # 100 bar in Pa

props_test = quartz.standardThermoProps(T_test, P_test)
print(f"\nStandard Thermodynamic Properties at 100Â°C, 100 bar:")
print(f"  G0 (Gibbs Free Energy): {props_test.G0.val() / 1000:.2f} kJ/mol")
print(f"  H0 (Enthalpy): {props_test.H0.val() / 1000:.2f} kJ/mol")

print("\nDone.")

