#pragma once

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelDEW.hpp>  // ActivityDHModel shared enum

/// ============================================================
/// ActivityModelPerplexDEW — capability matrix
/// ============================================================
///
/// Supported workflows
///   ✔  Pure aqueous chemistry with Perple_X-derived (DEW24HP622) database
///   ✔  T/P constraint: temperature() + pressure()
///   ✔  pH constraint: specs.pH()
///   ✔  Fugacity constraint: specs.fugacity()
///   ✔  EquilibriumSensitivity (dndw)
///   ✔  KineticsSolver (model-agnostic kinetics path)
///   ✔  Operator-splitting reactive transport (Python level)
///   ✔  Davies Debye-Hückel (dhModel = Davies, default; matches Perple_X aqact)
///   ✔  Extended Debye-Hückel with ion-size a (dhModel = ExtendedDH)
///   ✔  Perple_X GFSM gas-phase water-activity handoff (StateId protocol)
///   ✔  Optional fail-fast mode for missing/stale GFSM handoff
///   ✔  GFSM standard-state conflict detection (errorOnConflictingStandardState)
///   ✔  Unmapped GFSM coupling diagnostic (warnOnUnmappedGFSMCoupling)
///   ✔  Exports AqueousMixtureState into props.extra
///
/// Not supported in this model
///   ✘  Configurable water submodels (EOS, dielectric, Gibbs, Born) via waterOptions
///         (Perple_X uses ZD05 + Looyenga mixing internally; not switchable)
///   ✘  Extended DH b-dot term (no bExtended field)
///
/// Constructor signatures (C++ and Python)
///   ActivityModelPerplexDEW()
///   ActivityModelPerplexDEW(ActivityModelParamsPerplexDEW const& params)
///   ActivityModelPerplexDEW(ActivityDHModel model)
///
/// Debye-Hückel default
///   dhModel = ActivityDHModel::Davies
///
/// See also: ActivityModelDEW.hpp for the DEW-spreadsheet-backed variant
///           that declares the shared ActivityDHModel enum.
/// ============================================================

namespace Reaktoro {

struct ActivityModelParamsPerplexDEW
{
    /// Debye-Hückel variant. Default: Davies (matches Perple_X GFSM aqact exactly).
    /// Set to ActivityDHModel::ExtendedDH to use ionic-radius-resolved DH.
    ActivityDHModel dhModel = ActivityDHModel::Davies;

    /// When true, throw an exception if a neutral species is both a Perple_X
    /// GFSM fluid co-solvent (via formula lookup) and has a PerplexDEW HKF
    /// standard thermo model — double-counting the chemical potential.
    /// Default (false): emit a warning to stderr instead of throwing.
    bool errorOnConflictingStandardState = false;

    /// When true (default), emit a warning to stderr if a neutral species has
    /// a formula that does not match any Perple_X GFSM fluid key but becomes a
    /// match after stripping a common aqueous suffix (e.g., "_aq", ",aq",
    /// "(aq)").  This catches formula-string mismatches where GFSM coupling
    /// would be silently skipped.
    bool warnOnUnmappedGFSMCoupling = true;

    /// When true, PerplexDEW requires a fresh GFSM coupled-fluid handoff on
    /// every evaluation and throws if it must fall back to the aqueous-only
    /// solvent path. Useful for workflows that depend on strict gas-aqueous
    /// coupling and should never silently degrade.
    /// Default (false): allow fallback and continue.
    bool requireCoupledGFSMHandoff = false;
};

/// Return an aqueous activity model using Perple_X solvent internals.
///
/// Solvent state selection follows Perple_X GFSM (COH-Fluid+, model type 39):
///  - Pure H₂O  → slvnt0 path: ZD05 EOS + Fernandez/Sverjensky dielectric.
///  - Mixed fluid (H₂O + CO₂, CH₄, H₂S, SO₂, H₂, CO, N₂, NH₃, HF, C₂H₆, HCl)
///               → slvnt1 path: MRK pure volumes + Looyenga mixing rule for ε.
///
/// @param params  PerplexDEW activity-model options.
auto ActivityModelPerplexDEW(const ActivityModelParamsPerplexDEW& params)
    -> ActivityModelGenerator;

/// Convenience overload using only the Debye-Hückel model selector.
auto ActivityModelPerplexDEW(ActivityDHModel model = ActivityDHModel::Davies)
    -> ActivityModelGenerator;

} // namespace Reaktoro
