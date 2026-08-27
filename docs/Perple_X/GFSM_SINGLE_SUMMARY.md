# GFSM Single Summary (Reaktoro Perple_X)

**Date**: March 2026
**Scope**: GFSM (Perple_X type-39 path) and its regression architecture

---

## 1) What GFSM Is

GFSM (Generic Fluid Solution Model) is an explicit speciation-space fluid model.
Users pass a composition over the GFSM species basis and Reaktoro computes
mixture properties directly using MRK as the foundation with optional hybrid
pure-EOS substitution for H2O, CO2, and CH4.

Key distinction:
- Composition-space legacy workflows tabulate as a function of bulk sectioning variables.
- GFSM workflow evaluates directly from explicit species composition at (P, T).

---

## 2) Current GFSM Computational Pipeline

1. Input normalization and species validation against the GFSM-allowed set.
2. MRK baseline evaluation for active species.
3. Hybrid pure-EOS substitution for H2O/CO2/CH4 when configured.
4. Final activity/fugacity and volume assembly.
5. Optional electrolyte/HKF side calculations where tests require them.

---

## 3) EOS Coverage in the Active Regression Gate

- H2O anchors: mrk, hsmrk, cork, pseos, haar, zd05, zd09 (7).
- CO2 anchors: mrk, hsmrk, cork, pseos, zd09 (5).
- CH4 anchors: mrk, hsmrk, zd09 (3).

Note:
- CO2 brmrk is intentionally excluded from the regression gate because pure-CO2
  behavior is unstable for the tested operating range in the current reference
  generation flow.

---

## 4) Architecture Layers and Files

Core model and orchestration:
- PerpleXGFSMModel.hpp/cpp
- PerpleXFluidModel.hpp/cpp

MRK core:
- PerpleXMrkMixture.hpp/cpp
- PerpleXMrkPure.hpp/cpp
- PerpleXMrkParameters.hpp/cpp

Hybrid EOS dispatch and pure-EOS implementations:
- PerpleXHybridEos.hpp/cpp
- PerpleXPureEos.hpp/cpp

Electrolyte and HKF parity components:
- PerpleXElectrolyte.hpp/cpp
- PerpleXHKF.hpp/cpp

Regression harness and workflow:
- test_regression.cpp
- run_regression.py
- generate_gfsm_meemum_references.py
- test/gfsm_regression_matrix.csv
- test/gfsm_case_matrix.csv

---

## 5) Regression Architecture Snapshot

- Generator stage: `generate_gfsm_meemum_references.py` regenerates GFSM CSV references.
- Matrix stage: `test/gfsm_regression_matrix.csv` and `test/gfsm_case_matrix.csv` define
  required coverage and explicit stress/redox cases.
- Execution stage: `test_regression.cpp` performs EOS gate checks, matrix volume parity,
  DEW/HKF parity checks, and integration-flow assertions.
- Orchestration stage: `run_regression.py` wires generation and execution into one command,
  including strict gating via `--require-full-gfsm-pure`.

Current strict coverage footprint in the gate is 51 GFSM reference CSV cases:
- 15 pure anchors (H2O/CO2/CH4).
- 27 binary mixtures (H2O+CO2, H2O+CH4, CO2+CH4).
- 3 ternary (H2O+CO2+CH4).
- 6 O2/redox stress cases.

---

## 6) Current Validation Status

Latest recorded regression report:
- 11 tests passed, 0 failed.
- 116 assertions passed, 0 failed.
- 51 row-level assertions passed, 0 failed.

Integration-flow scenarios now include explicit high-component checks with
minimum active-species requirements aligned to the current generated dataset.

---

## 7) Summary

GFSM in Reaktoro is now documented and tested as a matrix-driven, generator-backed
architecture with explicit coverage gates and CI-friendly pass/fail metrics.
