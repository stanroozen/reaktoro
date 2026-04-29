# Complete Perple_X Implementation Test Checklist
## Ensuring Reaktoro ≡ Perple_X

---

## 1. PURE SPECIES EOS TESTS

### 1.1 Pure Water (H2O)
- [ ] **HSMRK at 0.1 GPa, 523 K** (dielectric, properties)
- [ ] **HSMRK at 0.5 GPa, 523 K** (dielectric, properties)
- [ ] **HSMRK at 1.0 GPa, 523 K** (dielectric, properties)
- [ ] **HSMRK at 2.0 GPa, 523 K** (extreme pressure)
- [ ] **HSMRK at 0.001 GPa, 373 K** (low pressure)
- [ ] **MRK pure H2O** (comparison with HSMRK)
- **Properties to validate:**
  - Molar volume (cm³/mol)
  - Gibbs energy (J/mol)
  - Fugacity coefficient (dimensionless)
  - Dielectric constant ε (unitless)

### 1.2 Pure CO2
- [ ] **CORK at 0.1 GPa, 523 K**
- [ ] **CORK at 0.5 GPa, 523 K**
- [ ] **CORK at 1.0 GPa, 523 K**
- [ ] **CORK at 2.0 GPa, 523 K**
- [ ] **MRK pure CO2** (comparison with CORK)
- **Properties to validate:**
  - Molar volume (cm³/mol)
  - Gibbs energy (J/mol)
  - Fugacity coefficient (dimensionless)

### 1.3 Other Pure Species (if supported)
- [ ] **NaCl** (if electrolyte model available)
- [ ] **SiO2** (if applicable)
- [ ] **Other minerals/species in model**

---

## 2. BINARY MIXTURE TESTS

### 2.1 H2O-CO2 (MRK Model)
**Composition series (constant P, T):**
- [ ] **Test at P=0.1 GPa, T=523 K, XCO2 = 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0**
- [ ] **Test at P=0.5 GPa, T=523 K, XCO2 = 0.0, 0.1, 0.2, ..., 1.0**
- [ ] **Test at P=1.0 GPa, T=523 K, XCO2 = 0.0, 0.1, 0.2, ..., 1.0**

**Pressure-temperature grid (fixed composition):**
- [ ] **XCO2 = 0.5, P ∈ [0.1, 0.5, 1.0, 2.0] GPa, T ∈ [373, 473, 523, 573] K**

**Properties to validate for all:**
- Molar volume (cm³/mol)
- Partial molar volumes: V̄(H2O), V̄(CO2)
- Fugacity coefficients: ln f(H2O), ln f(CO2)
- Total Gibbs energy
- Activity coefficients (if applicable)

### 2.2 H2O-CO2 (Hybrid EOS Model)
**Same composition series as 2.1:**
- [ ] **P=0.1 GPa, T=523 K, XCO2 = 0.0, 0.1, 0.2, ..., 1.0**
- [ ] **P=0.5 GPa, T=523 K, XCO2 = 0.0, 0.1, 0.2, ..., 1.0**
- [ ] **P=1.0 GPa, T=523 K, XCO2 = 0.0, 0.1, 0.2, ..., 1.0**

**Properties to validate:**
- Molar volume
- Partial molar volumes
- Fugacity coefficients
- Hybrid solvent contribution (gsolv)
- Comparison with MRK model

### 2.3 H2O-NaCl (if electrolyte available)
- [ ] **XNaCl = 0.0, 0.01, 0.05, 0.1, 0.2, at P=0.1 GPa, T=523 K**
- [ ] **Salting-out effects on H2O and CO2 (if ternary)**
- **Properties:**
  - Osmotic coefficient
  - Activity of water
  - Ion activity products

### 2.4 Other Binary Systems (if applicable)
- [ ] **H2O-N2** (if applicable)
- [ ] **H2O-CH4** (if applicable)
- [ ] **CO2-CH4** (if applicable)

---

## 3. TERNARY & HIGHER-ORDER SYSTEMS

### 3.1 H2O-CO2-NaCl (if available)
- [ ] **Variable salt concentration with CO2 partitioning**
- [ ] **Salting-out coefficient validation**
- [ ] **Activity coefficients**

### 3.2 Other Ternary Systems
- [ ] **Any other supported 3+ component systems**

---

## 4. THERMODYNAMIC PROPERTY DERIVATIVE TESTS

### 4.1 Temperature Derivatives
- [ ] **∂V/∂T|P,x** (thermal expansion)
- [ ] **∂G/∂T|P,x** (entropy from Gibbs)**
- [ ] **∂ln(f)/∂T|P,x** (fugacity temperature dependence)

### 4.2 Pressure Derivatives
- [ ] **∂V/∂P|T,x** (isothermal compressibility)
- [ ] **∂G/∂P|T,x** (molar volume from Gibbs)
- [ ] **∂ln(f)/∂P|T,x** (fugacity pressure dependence)

### 4.3 Composition Derivatives
- [ ] **∂V/∂xi|P,T** (partial molar volumes)
- [ ] **∂G/∂xi|P,T** (chemical potentials)
- [ ] **∂ln(f)/∂xi|P,T** (composition-fugacity coupling)

---

## 5. DIELECTRIC CONSTANT TESTS (epsh2o)

### 5.1 Pure Water Dielectric
- [ ] **ε(P, T) across grid: P ∈ [0.1, 0.5, 1.0, 2.0] GPa, T ∈ [373, 473, 523, 573] K**
- [ ] **Comparison with tabulated/literature values**
- [ ] **Temperature and pressure dependencies**

### 5.2 Mixed System Dielectric
- [ ] **ε(P, T, XCO2) for H2O-CO2 at key compositions**
- [ ] **Low-ε (non-polar) fluid mixing rules**

---

## 6. SHOCK G-FUNCTION TESTS (gf)

### 6.1 Pure H2O g-function
- [ ] **g-function at 10 key (P,T) points**
- [ ] **Sech² shape validation**
- [ ] **Derivative validation: dg/dT, dg/dP**

### 6.2 H2O-CO2 g-function
- [ ] **Mixed system g-function at key compositions**
- [ ] **Composition-dependence of g-function width**

---

## 7. DEBYE-HÜCKEL FACTOR TESTS (adh)

### 7.1 Pure Water Debye-Hückel
- [ ] **Ã(P,T) across grid: P ∈ [0.1, 0.5, 1.0, 2.0] GPa, T ∈ [373, 473, 523, 573] K**
- [ ] **Temperature and pressure dependencies**
- [ ] **Comparison with standard equations**

### 7.2 Ionic Strength Dependence
- [ ] **Debye-Hückel factor as function of ionic strength (if available)**

---

## 8. BORN OMEGA (ω) TESTS

### 8.1 Pure Ion Born Omega
- [ ] **ω(P, T) for representative ions at key conditions**
- [ ] **Charge-dependence validation**
- [ ] **Pressure and temperature dependencies**

### 8.2 Temperature/Pressure Derivatives
- [ ] **dω/dT|P** (temperature dependence)
- [ ] **dω/dP|T** (pressure dependence)

---

## 9. HKF MODEL TESTS (Helgeson-Kirkham-Flowers)

### 9.1 Single Ion Properties
- [ ] **Na⁺ at multiple (P, T) points**
- [ ] **Cl⁻ at multiple (P, T) points**
- [ ] **Other major ions: K⁺, Ca²⁺, Mg²⁺, SO₄²⁻, etc. (if available)**

**Properties to validate:**
- Gibbs energy (J/mol)
- Enthalpy (J/mol)
- Entropy (J/mol/K)
- Heat capacity (J/mol/K)
- Partial molar volume (cm³/mol)

### 9.2 Ion Pairs
- [ ] **NaCl ion pairing (if available)**
- [ ] **CaCl₂, MgCl₂ (if available)**

### 9.3 HKF Parameter Validation
- [ ] **∂(HKF properties)/∂T|P** (derivatives)
- [ ] **∂(HKF properties)/∂P|T** (derivatives)

---

## 10. ELECTROLYTE/AQUEOUS MODEL COUPLING

### 10.1 Electrolyte-Solvent Coupling
- [ ] **Activity coefficient of H2O near saturation (if coupled)**
- [ ] **Temperature and pressure effects on solvation**
- [ ] **Ionic strength effects on dielectric constant (if available)**

### 10.2 Mixed Solvent-Electrolyte
- [ ] **H2O-CO2-NaCl ternary at multiple compositions**
- [ ] **Salting-out of CO2 with varying NaCl**
- [ ] **Non-ideal mixing rules validation**

---

## 11. EXTREME & EDGE CASE TESTS

### 11.1 Low-Pressure Conditions
- [ ] **P = 0.001 GPa (1 bar), T = 373 K** (near saturation)
- [ ] **P = 0.01 GPa, T = 473 K**
- [ ] **Fugacity behavior near atmospheric**

### 11.2 High-Pressure Conditions
- [ ] **P = 2.0 GPa, T = 523 K** (extreme compression)
- [ ] **P = 5.0 GPa, T = 573 K** (if model supports)
- [ ] **Volume behavior at high pressure**

### 11.3 Low-Temperature Conditions
- [ ] **T = 273 K, P = 0.1 GPa** (cold subsurface)
- [ ] **T = 373 K, P = 0.5 GPa** (moderate cold)

### 11.4 High-Temperature Conditions
- [ ] **T = 773 K, P = 1.0 GPa** (crustal temperatures)
- [ ] **T = 873 K, P = 2.0 GPa** (mantle-like)

### 11.5 Supercritical Fluid Conditions
- [ ] **H2O above critical point: T > 647 K, P > 22 MPa**
- [ ] **CO2 above critical point: T > 304 K, P > 7.4 MPa**
- [ ] **Thermodynamic singularity behavior**

### 11.6 Two-Phase Region Edge Cases
- [ ] **Near H2O saturation curve**
- [ ] **Near CO2 liquid-gas boundary**
- [ ] **Narrow two-phase window testing (if applicable)**

---

## 12. NUMERICAL STABILITY & CONVERGENCE

### 12.1 Root-Finding Convergence
- [ ] **Iteration count for MRK root at various compositions**
- [ ] **Convergence at mixture extremes (X → 0 or 1)**
- [ ] **Convergence at high pressure (compression behavior)**

### 12.2 Analytical vs. Numerical Derivatives
- [ ] **Compare ∂V/∂T (analytical) vs. numerical finite difference**
- [ ] **Compare ∂G/∂P (analytical) vs. numerical finite difference**
- [ ] **Derivative accuracy at singular points**

### 12.3 Numerical Precision
- [ ] **Machine epsilon effects on fugacity coefficients**
- [ ] **Rounding error accumulation in long calculation chains**
- [ ] **Consistency of forward/backward calculations (e.g., V from G and vice versa)**

---

## 13. MODEL VARIANT TESTS

### 13.1 Low-Temperature MRK Option
- [ ] **useLowTMrk = true vs. false at T < 500 K**
- [ ] **Comparison with Perple_X low-T option**

### 13.2 Hybrid EOS Variants
- [ ] **Different hybrid option flags impact on results**
- [ ] **Comparison of hybrid vs. pure EOS at different conditions**

### 13.3 Electrolyte Coupling On/Off
- [ ] **enableElectrolyte = false (pure fluid)**
- [ ] **enableElectrolyte = true (with solvation)**
- [ ] **Impact on H2O and CO2 properties**

---

## 14. REFERENCE DATA MATRIX

**Minimum comprehensive test grid (if all tests cannot be done):**

### Core Grid: P-T-X Space
```
Pressures: [0.1, 0.5, 1.0, 2.0] GPa
Temperatures: [373, 473, 523, 573, 673] K
Compositions (H2O-CO2): [0.0, 0.25, 0.5, 0.75, 1.0]
Models: [MRK, Hybrid, CORK, HSMRK]
```

**This creates: 4 × 5 × 5 × 4 = 400 reference points minimum**

### Components Grid: P-T for Properties
```
Pressures: [0.1, 0.5, 1.0, 2.0] GPa
Temperatures: [373, 473, 523, 573, 673] K
For: [epsh2o, g-function, adh, ω]
```

**This creates: 4 × 5 × 4 = 80 reference points**

**Total minimum: ~480 validated reference points**

---

## 15. VALIDATION ACCEPTANCE CRITERIA

### Strict Tolerances (Core Properties)
- **Fugacity coefficient:** ≤ 1e-6 relative error
- **Molar volume:** ≤ 1e-6 relative error
- **Gibbs energy:** ≤ 100 J/mol absolute error

### Moderate Tolerances (Electrolyte/Solvation)
- **Dielectric constant:** ≤ 1e-4 absolute error
- **Debye-Hückel factor:** ≤ 1e-6 relative error
- **Born omega:** ≤ 1e-6 relative error
- **g-function:** ≤ 1e-8 Ångstrom absolute error

### Loose Tolerances (Derivatives/Extreme Conditions)
- **Numerical derivatives:** ≤ 1e-4 relative error
- **High-pressure volumes:** ≤ 1e-3 relative error (numerical sensitivity)

---

## 16. TESTING WORKFLOW

1. **Phase 1: Generate Perple_X References**
   - Use `generate_reference_matrix.py` to generate MRK/Hybrid/CORK/HSMRK references
   - Use `generate_reference_components.py` for dielectric/g-function/DH/Born/HKF
   - Save all as CSV files in `test/`

2. **Phase 2: Implement Reaktoro Tests**
   - Create test functions for each category above
   - Link to reference CSV files
   - Implement property comparisons

3. **Phase 3: Run Full Suite**
   - Execute all tests
   - Generate detailed mismatch report
   - Identify systematic vs. random errors

4. **Phase 4: Resolve Discrepancies**
   - Identify which components are correct/incorrect
   - Trace through source code (PerpleX Fortran vs. Reaktoro C++)
   - Adjust Reaktoro implementation as needed

5. **Phase 5: Final Validation**
   - Re-run full suite
   - Document all passing tests
   - Archive reference data and test results

---

## Current Status

**Passing (6/12):**
- [x] Pure H2O MRK
- [x] Pure CO2 MRK
- [x] Pure CO2 CORK
- [x] H2O-CO2 Binary (MRK)
- [x] H2O-CO2 Composition Series (MRK)
- [x] P-T Grid Test

**Failing (3/12):**
- [ ] Pure H2O HSMRK (volume mismatch)
- [ ] Hybrid MRK Binary (fugacity & volume mismatch)
- [ ] Hybrid MRK Composition Series (volume mismatches)

**Not Yet Implemented/Run:**
- [ ] HKF Sodium ion tests
- [ ] Born omega tests
- [ ] Ternary systems
- [ ] Extreme P-T cases
- [ ] Electrolyte coupling tests
- [ ] Derivative validation tests

---

## Next Actions

1. **Debug HSMRK volume discrepancy** → trace calculation in PerpleXPureEos.cpp
2. **Debug Hybrid EOS volume/fugacity** → review hybrid mixing rules in PerpleXHybridEos.cpp
3. **Activate HKF and Born tests** → ensure reference data generated
4. **Expand to ternary systems** → add H2O-CO2-NaCl tests
5. **Add derivative tests** → validate ∂G/∂P, ∂V/∂T numerically

