"""
Test DEW water model combinations.
Verifies that StandardThermoModelDEW can be configured with different
water EOS, dielectric, Gibbs-integration, and Born model combinations.
"""

import pytest


def _rkt():
    """Import reaktoro4py, skipping the test if unavailable."""
    try:
        import reaktoro4py as rkt

        return rkt
    except ImportError as e:
        pytest.skip(f"reaktoro4py not available: {e}")


def test_dew_database_loads():
    """DEW2019 aqueous database loads and has expected species count."""
    rkt = _rkt()
    db = rkt.DEWDatabase("dew2019-aqueous")
    assert len(db.species()) > 0, "DEW database should have species"


def test_dew_water_model_default():
    """DEW preset: ZhangDuan2005 + PowerFunction + DewIntegral + Shock92Dew."""
    rkt = _rkt()
    params = rkt.StandardThermoModelParamsDEW()
    params.waterOptions.eosModel = rkt.WaterEosModel.ZhangDuan2005
    params.waterOptions.dielectricModel = rkt.WaterDielectricModel.PowerFunction
    params.waterOptions.gibbsModel = rkt.WaterGibbsModel.DewIntegral
    params.waterOptions.bornModel = rkt.WaterBornModel.Shock92Dew
    model = rkt.StandardThermoModelDEW(params)
    assert model is not None


def test_dew_water_model_alternative():
    """Alternative: ZhangDuan2009 + JohnsonNorton1991 + DelaneyHelgeson1978."""
    rkt = _rkt()
    params = rkt.StandardThermoModelParamsDEW()
    params.waterOptions.eosModel = rkt.WaterEosModel.ZhangDuan2009
    params.waterOptions.dielectricModel = rkt.WaterDielectricModel.JohnsonNorton1991
    params.waterOptions.gibbsModel = rkt.WaterGibbsModel.DelaneyHelgeson1978
    params.waterOptions.bornModel = rkt.WaterBornModel.Shock92Dew
    model = rkt.StandardThermoModelDEW(params)
    assert model is not None


if __name__ == "__main__":
    test_dew_database_loads()
    test_dew_water_model_default()
    test_dew_water_model_alternative()
    print("All smoke tests passed.")
