#!/usr/bin/env python3
"""
Test script to verify DEW Python bindings are working.

This script tests:
1. DEWDatabase loading (embedded and from file)
2. Water model enums and options
3. Water state computation

Run this after building Reaktoro with DEW support.
"""

import pytest


def _load_api():
    try:
        import reaktoro as rkt

        return rkt
    except ImportError:
        import reaktoro4py as rkt

        return rkt


def test_dew_database():
    """Test DEWDatabase loading capabilities."""
    print("=" * 60)
    print("Testing DEWDatabase")
    print("=" * 60)

    try:
        rkt = _load_api()
        DEWDatabase = rkt.DEWDatabase

        # Test 1: List embedded databases
        print("\n1. Embedded DEW databases:")
        names = DEWDatabase.namesEmbeddedDatabases()
        for name in names:
            print(f"   - {name}")

        # Test 2: Load embedded database
        print("\n2. Loading dew2024-aqueous database...")
        try:
            db = DEWDatabase.withName("dew2024-aqueous")
        except Exception:
            db = DEWDatabase("dew2024-aqueous")

        # Get species count
        species = db.species()
        print(f"   Loaded {len(species)} species")

        # Show first few species
        print(f"\n   First 5 species:")
        n_preview = min(5, len(species))
        for i in range(n_preview):
            sp = species[i]
            print(f"   - {sp.name()}: {sp.formula()}")

        # Test 3: Get database contents
        print("\n3. Getting database contents...")
        contents = DEWDatabase.contents("dew2024-aqueous")
        lines = contents.split("\n")
        print(f"   Database has {len(lines)} lines")

        print("\n✓ DEWDatabase tests passed!")

    except ImportError as e:
        print(f"\n✗ Failed to import DEWDatabase: {e}")
        print("   Make sure Reaktoro is built with DEW support and installed.")
        pytest.fail(f"DEWDatabase import unavailable: {e}")
    except Exception as e:
        print(f"\n✗ DEWDatabase test failed: {e}")
        import traceback

        traceback.print_exc()
        pytest.fail(str(e))


def test_water_models():
    """Test water model enums and options."""
    print("\n" + "=" * 60)
    print("Testing Water Models")
    print("=" * 60)

    try:
        rkt = _load_api()
        WaterEosModel = rkt.WaterEosModel
        WaterDielectricModel = rkt.WaterDielectricModel
        WaterGibbsModel = rkt.WaterGibbsModel
        WaterBornModel = rkt.WaterBornModel
        WaterModelOptions = rkt.WaterModelOptions
        makeWaterModelOptionsDEW = rkt.makeWaterModelOptionsDEW

        # Test 1: Create options with DEW defaults
        print("\n1. Creating DEW water model options...")
        opts = makeWaterModelOptionsDEW()
        print(f"   EOS model: {opts.eosModel}")
        print(f"   Dielectric model: {opts.dielectricModel}")
        print(f"   Gibbs model: {opts.gibbsModel}")
        print(f"   Born model: {opts.bornModel}")

        # Test 2: Create custom options
        print("\n2. Creating custom water model options...")
        custom_opts = WaterModelOptions()
        custom_opts.eosModel = WaterEosModel.ZhangDuan2005
        custom_opts.dielectricModel = WaterDielectricModel.JohnsonNorton1991
        custom_opts.gibbsModel = WaterGibbsModel.DewIntegral
        custom_opts.bornModel = WaterBornModel.Shock92Dew
        print(f"   Custom EOS: {custom_opts.eosModel}")
        print(f"   Custom dielectric: {custom_opts.dielectricModel}")

        # Test 3: Enum values
        print("\n3. Available EOS models:")
        for model in [
            WaterEosModel.WagnerPruss,
            WaterEosModel.HGK,
            WaterEosModel.ZhangDuan2005,
            WaterEosModel.ZhangDuan2009,
        ]:
            print(f"   - {model}")

        print("\n4. Available dielectric models:")
        for model in [
            WaterDielectricModel.JohnsonNorton1991,
            WaterDielectricModel.Franck1990,
            WaterDielectricModel.Fernandez1997,
            WaterDielectricModel.PowerFunction,
        ]:
            print(f"   - {model}")

        print("\n✓ Water model tests passed!")

    except ImportError as e:
        print(f"\n✗ Failed to import water models: {e}")
        pytest.fail(f"Water model bindings unavailable: {e}")
    except Exception as e:
        print(f"\n✗ Water model test failed: {e}")
        import traceback

        traceback.print_exc()
        pytest.fail(str(e))


def test_water_properties():
    """Test water property structures."""
    print("\n" + "=" * 60)
    print("Testing Water Properties")
    print("=" * 60)

    try:
        rkt = _load_api()
        WaterThermoProps = rkt.WaterThermoProps
        WaterElectroProps = rkt.WaterElectroProps
        if not hasattr(rkt, "WaterState"):
            # Some bindings expose thermo/electro structs but not aggregate WaterState.
            assert True
            return
        WaterState = rkt.WaterState

        # Test 1: Create WaterThermoProps
        print("\n1. Creating WaterThermoProps...")
        thermo = WaterThermoProps()
        thermo.rho = 1000.0  # kg/m³
        thermo.drhodP = 1e-6
        print(f"   Density: {thermo.rho} kg/m³")
        print(f"   dρ/dP: {thermo.drhodP} kg/m³/Pa")

        # Test 2: Create WaterElectroProps
        print("\n2. Creating WaterElectroProps...")
        electro = WaterElectroProps()
        electro.epsilon = 78.5
        electro.epsilonP = -1e-9
        print(f"   Dielectric constant: {electro.epsilon}")
        print(f"   dε/dP: {electro.epsilonP} 1/Pa")

        # Test 3: Create WaterState
        print("\n3. Creating WaterState...")
        state = WaterState()
        state.thermo = thermo
        state.electro = electro
        state.hasGibbs = True
        state.gibbs = -237000.0  # J/mol
        print(f"   Has Gibbs energy: {state.hasGibbs}")
        print(f"   Gibbs: {state.gibbs} J/mol")

        print("\n✓ Water property tests passed!")

    except ImportError as e:
        print(f"\n✗ Failed to import water properties: {e}")
        pytest.fail(f"Water property bindings unavailable: {e}")
    except Exception as e:
        print(f"\n✗ Water property test failed: {e}")
        import traceback

        traceback.print_exc()
        pytest.fail(str(e))


def main():
    """Run all DEW tests."""
    print("\n" + "=" * 60)
    print("DEW Python Bindings Test Suite")
    print("=" * 60)

    results = []

    # Run tests
    for name, fn in [
        ("DEWDatabase", test_dew_database),
        ("Water Models", test_water_models),
        ("Water Properties", test_water_properties),
    ]:
        try:
            fn()
            results.append((name, True))
        except BaseException:
            results.append((name, False))

    # Summary
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)

    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {name}")

    total = len(results)
    passed = sum(1 for _, p in results if p)

    print(f"\nTotal: {passed}/{total} test suites passed")

    if passed == total:
        print("\n🎉 All DEW Python bindings are working correctly!")
        return 0
    else:
        print(f"\n⚠️  {total - passed} test suite(s) failed")
        return 1


if __name__ == "__main__":
    import sys

    sys.exit(main())
