import sys
import os
import pytest

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)


def test_run_quartz_dew():
    """Regression: quartz solubility vs experimental data."""
    try:
        import reaktoro4py  # noqa: F401
    except ImportError as e:
        pytest.skip(f"reaktoro4py not available: {e}")

    from unittest.mock import patch
    import quartz_solubility_analysis_v2_dew24 as _m
    with patch("sys.argv", [__file__, "--backend", "DEW"]):
        _m.main()


if __name__ == "__main__":
    sys.argv = [sys.argv[0], "--backend", "DEW"]
    import quartz_solubility_analysis_v2_dew24 as _m
    _m.main()
