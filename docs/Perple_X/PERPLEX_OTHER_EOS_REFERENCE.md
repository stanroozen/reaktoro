# Perple_X Non-GFSM Pure EoS Reference

This document lists pure equations of state (ifug options) present in Perple_X that are **not** used by the GFSM hybrid model in the Reaktoro extension. These are documented here to keep the main GFSM reference focused on MRK + hybrid substitutions for H2O/CO2/CH4.

## Non-GFSM Pure EoS Inventory

| ifug | Name | Species | Status in Reaktoro | Reference |
|------|------|---------|-------------------|-----------|
| 7 | GFSEOS | H₂O | ❌ Not exposed | Gottschalk (2007) |
| 8-10 | PRSV, PREOS, PRKEOS | Various | 🔄 Documented | Peng-Robinson variants |
| 11 | H2O (low-P correlation) | H₂O | 🔄 Documented | Simple correlation |
| 12 | SOAVE | Mixed | 🔄 Documented | Soave (1972) |
| 13 | TOOP | H₂O-CO₂ | 🔄 Documented | Toop-Samis mixing |
| 14-20 | Various | Mixed | 🔄 Documented | Special formulations |
| 21+ | VTEOS | Various | 🔄 Documented | Volume-temperature EoS |

## Scope Notes

- These models are listed in Perple_X data and code (ifug options) but are **not** used by the GFSM hybrid substitution path in the Reaktoro implementation.
- This document is a placeholder for future expansion if these models are implemented or exposed in Reaktoro.
- For the GFSM hybrid pipeline, see [Reaktoro/Extensions/Perple_X/PERPLEX_COMPLETE_REFERENCE.md](Reaktoro/Extensions/Perple_X/PERPLEX_COMPLETE_REFERENCE.md).

## Where These Are Referenced

In the main reference, these entries were removed from the GFSM-relevant table to keep the scope focused. If you want these models promoted into the GFSM workflow, we can add:

1. Pure-EoS wrappers in PerpleXPureEos.
2. Hybrid selection options in PerpleXHybridEos.
3. Validation tests against Perple_X outputs.
