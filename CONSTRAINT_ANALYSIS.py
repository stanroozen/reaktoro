#!/usr/bin/env python
"""
ANALYSIS: Constraint Types for Willemite Solubility Control

This script evaluates available constraint options and their thermodynamic rigor.
"""

analysis = """
╔════════════════════════════════════════════════════════════════════════════╗
║        CONSTRAINT TYPE COMPARISON FOR WILLEMITE SOLUBILITY CONTROL        ║
╚════════════════════════════════════════════════════════════════════════════╝

AVAILABLE WORKING CONSTRAINTS IN REAKTORO:
──────────────────────────────────────────────────────────────────────────────

FOR SOLUTE/AQUEOUS SPECIES (SiO2,aq, HS⁻, etc.):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1. ACTIVITY-BASED (lgActivity, lnActivity, activity)  ✓ CURRENT CHOICE
   ─────────────────────────────────────────────────────────────
   • Constraint: lg(a[species]) = X or a[species] = 10^X
   • Units: dimensionless (activity is a normalized measure)
   • Thermodynamic basis: Appears in all equilibrium constants
     - K_eq = ∏(a_i^ν_i)  [fundamental]
   • Meaning: Controls concentration relative to reference state
     - a > 1 → supersaturation
     - a < 1 → undersaturation
     - a = 1 → equilibrium with reference mineral

   ✓ Pros:
     - Direct: activity is what equilibrium constants use
     - Practical: intuitive scaling with concentration
     - Easy: no special character parsing issues
     - Standard: used in most geochemical codes (PHREEQC, etc.)

   ✗ Cons:
     - Normalized to reference state (not absolute)
     - For SiO2,aq: depends on chosen reference pressure/temperature
     - For HS⁻: mixes redox state + sulfur speciation (less clean)


2. GIBBS ENERGY / CHEMICAL POTENTIAL  ✗ WANTED BUT BLOCKED
   ──────────────────────────────────────────────────────────
   • Constraint: G(species) = X or μ(species) = X  [at given T, P]
   • Units: J/mol
   • Thermodynamic basis: Fundamental equation μ = μ° + RT·ln(a)
   • Meaning: Controls absolute equilibrium potential
     - More negative → more stable/favorable
     - Directly shows driving force for dissolution

   ✓ Pros:
     - Fundamental: no reference state dependence
     - Universal: same units across all systems
     - Rigorous: G° is standard thermodynamic property
     - Intuitive for driving force: ΔG = G_reactants - G_products

   ✗ Cons:
     - API limitation: chemicalPotential() can't parse "SiO2,aq" (has comma)
     - Less standard in geochemistry (most use activities)
     - CURRENT WORKAROUND: convert μ → activity, then use lgActivity
       [mathematically identical, but roundabout]


3. FUGACITY (for gases)  ✓ AVAILABLE FOR O2
   ────────────────────────────────────────────
   • Constraint: f(gas) = X
   • Units: bar (or atm)
   • Thermodynamic basis: f = φ·P where φ is fugacity coefficient
   • Meaning: "Effective pressure" accounting for non-ideality
     - f = P for ideal gases
     - f < P if molecules attract (common)
     - f > P if molecules repel (rare)

   ✓ Pros:
     - Most practical for gases: directly related to pressure
     - Already used in current script for O2
     - Very standard: PHREEQC, PerpleX both use this

   ✗ Cons:
     - Can't use for aqueous species (not designed for it)
     - Less meaningful for dissolved species


4. REDOX POTENTIAL (Eh or pE)  ✓ AVAILABLE, BETTER FOR HS⁻
   ──────────────────────────────────────────────────────────
   • Constraint: Eh = X (volts) or pE = X (dimensionless)
   • Units: Eh in mV or volts; pE = -log(e⁻ activity)
   • Thermodynamic basis: O2/e⁻ pair redox equation
   • Meaning: Controls oxidation state of the system
     - More positive Eh → oxidized (SO4²⁻, O2)
     - More negative Eh → reduced (H2S, HS⁻, Fe²⁺)

   ✓ Pros:
     - FUNDAMENTAL for redox-sensitive systems
     - Controls entire redox state (not just HS⁻ alone)
     - Affects speciation: SO4²⁻ ↔ HS⁻ ↔ H2S automatically
     - Cleaner than constraining HS⁻ activity (which mixes speciation)

   ✗ Cons:
     - Doesn't directly control SiO2,aq
     - Less intuitive than activity (need to know reference electrode)


5. ENTHALPY / ENTROPY  ✓ AVAILABLE BUT UNUSUAL
   ────────────────────────────────────────────
   • Constraint: H(system) = X or S(system) = X
   • Thermodynamic basis: H = U + PV, S = fundamental property
   • Meaning: Controls energy content or disorder

   ✗ Cons:
     - Rarely used for equilibrium constraints
     - Not standard in Reaktoro tutorials
     - Doesn't directly control solubility


╔════════════════════════════════════════════════════════════════════════════╗
║                        RECOMMENDATION MATRIX                              ║
╚════════════════════════════════════════════════════════════════════════════╝

GOAL: Control Zn solubility in willemite as function of:
      (1) SiO2 availability, (2) Redox state, (3) Oxygen fugacity

┌─────────────────────────────────────────────────────────────────────────┐
│ CURRENT APPROACH (Activity-based):                                      │
├─────────────────────────────────────────────────────────────────────────┤
│ SiO2,aq:  constraint on lg(a[SiO2,aq]) via μ conversion ✓ GOOD         │
│ HS⁻:      constraint on lg(a[HS⁻])  via μ conversion   ~ ACCEPTABLE    │
│ O2:       constraint on fugacity(O2)                    ✓ GOOD          │
│                                                                         │
│ Score: PRACTICAL & WORKING (EXIT=0)                                   │
│ Rigor: SUFFICIENT (activity is standard, but roundabout)               │
└─────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────┐
│ RECOMMENDED IMPROVEMENT (Hybrid):                                       │
├─────────────────────────────────────────────────────────────────────────┤
│ SiO2,aq:  constraint on lg(a[SiO2,aq])               ✓ DIRECT APPROACH │
│ Eh/pE:    constraint on Eh (or pE)                  ✓ BETTER RIGOR    │
│ O2:       constraint on fugacity(O2)                ✓ SAME             │
│                                                                         │
│ Why better:                                                            │
│ • Eh/pE directly controls H2S ↔ SO4²⁻ speciation                       │
│ • Not just constraining HS⁻ activity (which mixes in speciation)       │
│ • Standard geochemical approach (HSC, EQ3/6, PHREEQC)                  │
│ • More intuitive: know redox potential, get redox state                │
│                                                                         │
│ Score: MORE THERMODYNAMICALLY CLEAN                                    │
│ Effort: LOW (one function change)                                      │
└─────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────┐
│ MOST FUNDAMENTAL APPROACH (Gibbs Energy):                               │
├─────────────────────────────────────────────────────────────────────────┤
│ SiO2,aq:  constraint on G(SiO2,aq) in J/mol          ✗ API BLOCKED    │
│ [All]:    constraint via μ converted to activity     ✓ WORKAROUND OK  │
│                                                                         │
│ Status: Currently using this (indirectly via activity conversion)      │
│ Why not direct: chemicalPotential() can't parse "SiO2,aq" with comma   │
│ Rigor: HIGHEST (G is absolute, no reference state)                    │
│ Practicality: LOWER (less standard, more obscure)                      │
└─────────────────────────────────────────────────────────────────────────┘

╔════════════════════════════════════════════════════════════════════════════╗
║                             FINAL VERDICT                                 ║
╚════════════════════════════════════════════════════════════════════════════╝

YOUR CURRENT IMPLEMENTATION (Activity via μ conversion):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

✓ Is this the best way?
  DEPENDS ON YOUR GOAL:

  → If goal = "most thermodynamically rigorous":
    YES, because it controls absolute μ (fundamental)
    (Better than pure activity, worse than direct chemicalPotential)

  → If goal = "most practical & interpretable":
    PARTIALLY. Activity is more practical, but current approach is fine.

  → If goal = "explore redox sensitivity":
    NO. You should use Eh/pE instead of HS⁻ activity.


SUGGESTED ALTERNATIVES (in order of recommendation):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1. SIMPLIFY (Use pure activity constraints):
   • Current: lg(a[SiO2,aq]) ← computed from μ
   • Better:  lg(a[SiO2,aq]) ← specified directly
   • Reason:  Simpler code, more standard, same physics
   • Change:  Delete μ → activity conversion, use lgActivity directly
   • Benefit: Removes indirection, clearer intent


2. IMPROVE REDOX CONTROL (Replace HS⁻ activity with Eh):
   • Current: lg(a[HS⁻]) ← computed from μ
   • Better:  Eh (volts) ← specified directly
   • Reason:  Controls redox state properly, cleaner speciation
   • Change:  compute_chemical_potential_sensitivity(..., "HS-", ...)
              →  compute_redox_sensitivity(..., Eh_range_mv, ...)
   • Benefit: More rigorous for S redox chemistry


3. USE GIBBS ENERGY (if you really want absolute μ):
   • Current: lg(a[SiO2,aq]) ← converted from μ
   • Better:  gibbsEnergy() constraint works in Reaktoro
   • Reason:  More fundamental, no reference state
   • Change:  Implement compute_gibbs_energy_sensitivity(...)
   • Benefit: Most thermodynamically transparent


BOTTOM LINE:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Your current approach (activity via μ conversion) is:
  ✓ Mathematically correct
  ✓ Thermodynamically rigorous
  ✓ Working (EXIT=0, convergence verified)
  ~ But roundabout (could simplify)

RECOMMENDATION:
  Keep as-is IF: you want maximal thermodynamic transparency
  Simplify to:  direct activity constraints (cleaner code)
  Enhance to:   activity(SiO2) + Eh/pE (redox) (better physics)
  Switch to:    gibbsEnergy() + activity (most fundamental)
"""

print(analysis)
