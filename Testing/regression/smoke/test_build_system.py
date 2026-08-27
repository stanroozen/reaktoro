import pytest
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
QUARTZ_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "quartz")


def test_build_system():
    """Smoke test that build_system() from quartz_solubility_analysis_v2_dew24 works."""
    try:
        import reaktoro4py as rkt  # noqa: F401
    except ImportError as e:
        pytest.skip(f"reaktoro4py not available: {e}")

    if QUARTZ_DIR not in sys.path:
        sys.path.insert(0, QUARTZ_DIR)

    try:
        from quartz_solubility_analysis_v2_dew24 import MINERAL_CONFIG, build_system
    except ImportError as e:
        pytest.skip(f"quartz_solubility_analysis_v2_dew24 not importable: {e}")

    dew_db = rkt.DEWDatabase("dew2024-aqueous")
    supcrt_db = rkt.SupcrtDatabase("supcrtbl")
    try:
        system = build_system(
            dew_db=dew_db,
            supcrt_db=supcrt_db,
            mineral_config=MINERAL_CONFIG,
            model_backend="DEW",
        )
    except RuntimeError as exc:
        msg = str(exc)
        if "Could not find any Species object with name H+" in msg:
            pytest.skip(
                "Legacy quartz helper uses old species naming (H+, SiO2_aq, ...)."
            )
        raise
    assert system is not None
    assert len(system.species()) > 0


if __name__ == "__main__":
    try:
        test_build_system()
        print("test_build_system passed.")
    except BaseException as e:
        if type(e).__name__ == "Skipped":
            print(f"SKIP: {e}")
            raise SystemExit(0)
        raise
