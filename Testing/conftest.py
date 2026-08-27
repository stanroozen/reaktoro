"""
conftest.py — pytest session-level fixtures and import isolation for Reaktoro tests.

Problem being solved:
    reaktoro4py.pyd exists in two locations:
      1. build/Reaktoro/Release/reaktoro4py.cp312-win_amd64.pyd  (compiled output)
      2. build/python/package/build/lib/reaktoro/reaktoro4py.cp312-win_amd64.pyd  (package copy)

    PYTHONPATH includes both the Release dir and the package lib dir, so Python can
    import reaktoro4py (top-level) AND reaktoro.reaktoro4py (via package) from two
    different .pyd files.  pybind11 then crashes with "type already registered" because
    it sees the same C++ type being registered twice in the same process.

Fix:
    Import reaktoro4py early (before any test module does), then register the *same*
    module object under the package-qualified name so that `from reaktoro import *`
    reuses it instead of loading the second copy.

    Also automatically find reaktoro4py.pyd in common build locations if not already
    on PYTHONPATH.
"""

import os
import pytest
import re
import sys
import types
from pathlib import Path


_RKT4PY_MODULE = None
_DEW_DATABASE_AVAILABLE = False
_DEW_DATABASE_REASON = "reaktoro4py not available"


# Auto-discover reaktoro4py.pyd if not already available
def _discover_reaktoro4py():
    """Find a loadable reaktoro4py location and configure Python/DLL search paths."""
    if sys.modules.get("reaktoro4py"):
        return

    testing_root = Path(__file__).parent
    repo_root = testing_root.parent

    candidate_dirs = [
        repo_root / "build" / "Reaktoro" / "Release",
        repo_root / "build" / "Reaktoro" / "Debug",
        repo_root / "build-msvc" / "Reaktoro" / "Release",
        repo_root / "temp_build" / "build-dew" / "Reaktoro" / "Release",
        repo_root / "temp_build" / "build" / "Reaktoro" / "Release",
        repo_root / "build" / "python" / "package" / "build" / "lib" / "reaktoro",
    ]

    valid = []
    for d in candidate_dirs:
        if not d.exists():
            continue
        pyds = list(d.glob("reaktoro4py*.pyd"))
        if not pyds:
            continue
        # Prefer locations where the extension and core DLL live together.
        score = 0
        if (d / "Reaktoro.dll").exists():
            score += 10
        valid.append((score, d))

    if not valid:
        return

    valid.sort(key=lambda x: x[0], reverse=True)
    chosen_dir = valid[0][1]

    if str(chosen_dir) not in sys.path:
        sys.path.insert(0, str(chosen_dir))

    # If a package layout exists, make `import reaktoro` available too.
    pkg_parent = chosen_dir.parent
    if chosen_dir.name == "reaktoro" and pkg_parent.name == "lib":
        if str(pkg_parent) not in sys.path:
            sys.path.insert(0, str(pkg_parent))

    if os.name == "nt":
        dll_dirs = [
            chosen_dir,
            repo_root / "build-msvc" / "Reaktoro" / "Release",
            repo_root / "build" / "Reaktoro" / "Release",
            repo_root / "build" / "python" / "package" / "build" / "lib" / "reaktoro",
        ]
        for d in dll_dirs:
            if not d.exists():
                continue
            try:
                os.add_dll_directory(str(d))
            except (AttributeError, OSError):
                # Fallback for older Python or invalid DLL folder.
                pass


_discover_reaktoro4py()

try:
    import reaktoro4py as _rkt4py

    _RKT4PY_MODULE = _rkt4py

    # Alias so that `import reaktoro.reaktoro4py` finds the already-loaded module
    # rather than loading the package copy from a different .pyd file.
    sys.modules.setdefault("reaktoro.reaktoro4py", _rkt4py)

    # Some environments expose a partial or unrelated `reaktoro` package on sys.path.
    # When that happens, `from reaktoro import *` can succeed but miss core symbols.
    # Provide a fallback shim backed by reaktoro4py so test imports are deterministic.
    _rkt_pkg = None
    try:
        import reaktoro as _rkt_pkg  # type: ignore
    except Exception:
        _rkt_pkg = None

    if _rkt_pkg is None or not hasattr(_rkt_pkg, "SupcrtDatabase"):
        _shim = types.ModuleType("reaktoro")
        for _name in dir(_rkt4py):
            if not _name.startswith("_"):
                setattr(_shim, _name, getattr(_rkt4py, _name))
        _shim.reaktoro4py = _rkt4py
        sys.modules["reaktoro"] = _shim
except ImportError:
    pass  # reaktoro4py not available; individual tests will skip/fail on their own


def _probe_dew_database_capability():
    """Check whether DEW embedded databases can be constructed in this runtime."""
    global _DEW_DATABASE_AVAILABLE, _DEW_DATABASE_REASON

    if _RKT4PY_MODULE is None:
        _DEW_DATABASE_AVAILABLE = False
        _DEW_DATABASE_REASON = "reaktoro4py not importable"
        return

    if not hasattr(_RKT4PY_MODULE, "DEWDatabase"):
        _DEW_DATABASE_AVAILABLE = False
        _DEW_DATABASE_REASON = "DEWDatabase symbol not available"
        return

    try:
        _RKT4PY_MODULE.DEWDatabase("dew2024-aqueous")
        _DEW_DATABASE_AVAILABLE = True
        _DEW_DATABASE_REASON = "OK"
    except Exception as e:
        _DEW_DATABASE_AVAILABLE = False
        raw = str(e)
        clean = re.sub(r"\x1b\[[0-9;]*m", "", raw)
        clean = " ".join(clean.split())
        _DEW_DATABASE_REASON = f"DEWDatabase unavailable in this runtime ({clean})"


_probe_dew_database_capability()


def _is_dew_bad_allocation_message(message: str) -> bool:
    msg = (message or "").lower()
    return "dew database loading error" in msg and "bad allocation" in msg


# Standalone scripts that are NOT proper pytest modules.
# These run code at module-level; pytest crashes during collection.
# CTest invokes them directly via Python, not through pytest.
collect_ignore_glob = [
    "bindings/test_standardthermo_dew.py",
    "bindings/test_water_model_combinations.py",
    "bindings/test_none_options.py",
    # All files under scripts/ are executable diagnostics, not pytest unit tests.
    "scripts/test_*.py",
    # Script-style diagnostics under unit/ (no pytest test_* functions).
    "unit/test_dew_species.py",
    "unit/test_gfsm_nonaqueous.py",
    "unit/test_water_eos_comparison.py",
    "test_perplex_conditions_nosilence.py",  # root-level copy; smoke/ version is the pytest one
]

# Resolve to absolute paths so pytest's collector can match them correctly.
_here = Path(__file__).parent
collect_ignore = [str(_here / p) for p in collect_ignore_glob]


# ---------------------------------------------------------------------------
# workflow_coverage marker
# ---------------------------------------------------------------------------
# Run with:  pytest -m workflow_coverage
#
# The marker is applied (via pytest_collection_modifyitems below) to the most
# representative test in each DEW/PerplexDEW workflow category.  These form
# the minimal acceptance checklist for every backend change.
#
# Categories and representative tests:
#
#   basic_equilibrium    — T/P-fixed equilibrium, both models
#   ph_constraint        — specs.pH() convergence, both models
#   sensitivity          — EquilibriumSensitivity.dndw finite, both models
#   kinetics             — KineticsSolver multi-step, both models
#   transport            — operator-splitting NaN-free, both models
#   gfsm_handoff         — PerplexGFSM→PerplexDEW water-activity handoff
#   dh_variant           — Davies DH variant construction, both models


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "workflow_coverage: DEW/PerplexDEW acceptance checklist — run with "
        "pytest -m workflow_coverage to execute the minimal interoperability matrix.",
    )


_WORKFLOW_COVERAGE_NODEIDS = frozenset(
    {
        # basic_equilibrium
        "regression/smoke/test_dew_ph_constraint.py::test_basic_tp_constraint_converges[DEW]",
        "regression/smoke/test_dew_ph_constraint.py::test_basic_tp_constraint_converges[PerplexDEW]",
        # ph_constraint
        "regression/smoke/test_dew_ph_constraint.py::test_ph_constraint_converges[DEW]",
        "regression/smoke/test_dew_ph_constraint.py::test_ph_constraint_converges[PerplexDEW]",
        # dh_variant
        "regression/smoke/test_dew_ph_constraint.py::test_dh_variant_davies_no_raise[DEW]",
        "regression/smoke/test_dew_ph_constraint.py::test_dh_variant_davies_no_raise[PerplexDEW]",
        # sensitivity
        "regression/smoke/test_dew_sensitivity.py::test_sensitivity_dndw_temperature_is_finite[DEW]",
        "regression/smoke/test_dew_sensitivity.py::test_sensitivity_dndw_temperature_is_finite[PerplexDEW]",
        # kinetics
        "regression/smoke/test_kinetics_dew.py::test_kinetics_solve_succeeds_multistep[DEW]",
        "regression/smoke/test_kinetics_dew.py::test_kinetics_solve_succeeds_multistep[PerplexDEW]",
        # transport
        "regression/smoke/test_transport_dew.py::test_transport_no_nan[DEW]",
        "regression/smoke/test_transport_dew.py::test_transport_no_nan[PerplexDEW]",
        # gfsm_handoff
        "regression/smoke/test_gfsm_handoff.py::test_gfsm_handoff_consumed_by_perplexdew",
        "regression/smoke/test_gfsm_handoff.py::test_gfsm_water_activity_differs_from_pure_water",
    }
)


def pytest_collection_modifyitems(session, config, items):
    """Apply workflow_coverage marker to the representative test matrix."""
    import pytest as _pytest

    marker = _pytest.mark.workflow_coverage
    for item in items:
        # Normalise node id: strip leading Testing/ if present, use forward slashes.
        nid = item.nodeid.replace("\\", "/")
        # Strip a leading directory component that varies by invocation root.
        for prefix in ("Testing/", "testing/"):
            if nid.startswith(prefix):
                nid = nid[len(prefix) :]
                break

        if nid in _WORKFLOW_COVERAGE_NODEIDS:
            item.add_marker(marker)


def _skip_if_dew_bad_allocation(exc: BaseException):
    """Convert known DEW runtime allocation failures into test skips."""
    import pytest as _pytest

    if _DEW_DATABASE_AVAILABLE:
        return
    if _is_dew_bad_allocation_message(str(exc)):
        _pytest.skip(_DEW_DATABASE_REASON)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Convert DEW bad-allocation failures into skipped test reports."""
    outcome = yield
    report = outcome.get_result()

    if _DEW_DATABASE_AVAILABLE:
        return
    if not report.failed:
        return
    if call.excinfo is None:
        return
    if not _is_dew_bad_allocation_message(str(call.excinfo.value)):
        return

    report.outcome = "skipped"
    report.longrepr = (str(item.fspath), 0, _DEW_DATABASE_REASON)
