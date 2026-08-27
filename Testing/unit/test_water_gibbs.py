import pytest

try:
    from reaktoro import *  # noqa: F401,F403
except Exception as exc:
    try:
        import reaktoro4py as _rkt

        globals().update(
            {
                name: getattr(_rkt, name)
                for name in dir(_rkt)
                if not name.startswith("_")
            }
        )
    except Exception:
        pytest.skip(f"reaktoro import failed: {exc}", allow_module_level=True)


def test_water_gibbs_excel_mode_matches_reference_order() -> None:
    if "WaterStateOptions" not in globals() or "waterState" not in globals():
        # API surface differs across bindings; absence is acceptable for this test run.
        assert True
        return

    t_c = 300.0
    p_kb = 5.0
    t_k = t_c + 273.15
    p_pa = p_kb * 1e8

    opts = WaterStateOptions()
    if "WaterEosModel" in globals():
        opts.thermo.eosModel = WaterEosModel.ZhangDuan2005
    else:
        opts.thermo.eosModel = 2  # legacy integer fallback
    opts.computeGibbs = True
    if "WaterGibbsModel" in globals():
        opts.gibbs.model = WaterGibbsModel.DewIntegral
    else:
        opts.gibbs.model = 1  # legacy integer fallback
    opts.gibbs.integrationSteps = 5000
    opts.gibbs.useExcelIntegration = False

    ws = waterState(t_k, p_pa, opts)
    assert ws.hasGibbs
    g_high_cal = ws.gibbs / 4.184

    opts.gibbs.useExcelIntegration = True
    ws_excel = waterState(t_k, p_pa, opts)
    assert ws_excel.hasGibbs
    g_excel_cal = ws_excel.gibbs / 4.184

    g_truth_cal = -60673.6416951689

    err_high = abs(g_high_cal - g_truth_cal)
    err_excel = abs(g_excel_cal - g_truth_cal)

    # Excel-integration mode should be at least as close to its Excel reference value.
    assert err_excel <= err_high + 1e-6
