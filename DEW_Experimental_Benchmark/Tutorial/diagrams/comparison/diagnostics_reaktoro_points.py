import csv
import os
import sys

import numpy as np
import autodiff as ad

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
BASE = os.path.abspath(os.path.dirname(__file__))

PYD_DIR = os.path.join(REPO, "build", "python", "package", "build", "lib", "reaktoro")
if PYD_DIR not in sys.path:
    sys.path.insert(0, PYD_DIR)
os.add_dll_directory(PYD_DIR)

import reaktoro4py as rkt


def nearest_idx(vals, x):
    vals = np.asarray(vals, dtype=float)
    return int(np.argmin(np.abs(vals - float(x))))


def write_case(case_name, grid, xvals, yvals, species, out_path):
    pred = np.asarray(grid.predominantSpeciesGrid(species), dtype=float)

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "case",
                "pH",
                "Eh_V",
                "pred_idx",
                "pred_species_reaktoro",
                "top1_species_reaktoro",
                "top1_log_activity",
                "top2_species_reaktoro",
                "top2_log_activity",
            ]
        )

        states = list(grid.states)
        ny = len(yvals)

        for pH in xvals:
            ix = nearest_idx(grid.xvalues, pH)
            for Eh in yvals:
                iy = nearest_idx(grid.yvalues, Eh)

                idx_val = pred[ix, iy]
                pred_name = ""
                if (
                    np.isfinite(idx_val)
                    and int(idx_val) >= 0
                    and int(idx_val) < len(species)
                ):
                    pred_name = species[int(idx_val)]

                flat = ix * ny + iy
                state = states[flat]
                props = rkt.ChemicalProps(state)

                vals = []
                for sp in species:
                    try:
                        v = float(props.speciesActivityLg(sp))
                    except Exception:
                        v = float("nan")
                    vals.append((sp, v))

                vals = [x for x in vals if np.isfinite(x[1])]
                vals.sort(key=lambda z: z[1], reverse=True)

                top1_sp = vals[0][0] if len(vals) > 0 else ""
                top1_v = vals[0][1] if len(vals) > 0 else float("nan")
                top2_sp = vals[1][0] if len(vals) > 1 else ""
                top2_v = vals[1][1] if len(vals) > 1 else float("nan")

                # If field extraction returns NaN, fall back to top-activity species at that point.
                if not pred_name and top1_sp:
                    pred_name = top1_sp

                w.writerow(
                    [
                        case_name,
                        float(pH),
                        float(Eh),
                        "" if not np.isfinite(idx_val) else int(idx_val),
                        pred_name,
                        top1_sp,
                        top1_v,
                        top2_sp,
                        top2_v,
                    ]
                )


def main():
    db = rkt.SupcrtDatabase("supcrtbl")

    # --- Fe Pourbaix (CHNOSZ-like) ---
    aq_p = [
        "H2O(aq)",
        "H+",
        "OH-",
        "Fe+2",
        "Fe+3",
        "FeO(aq)",
        "FeO+",
        "FeO2-",
        "FeOH+",
        "FeOH+2",
        "HFeO2(aq)",
        "HFeO2-",
    ]
    min_p = ["Ferropericlase", "Goethite", "Hematite", "Iron", "Magnetite"]
    species_p = [
        "Fe+2",
        "Fe+3",
        "FeO(aq)",
        "FeO+",
        "FeO2-",
        "FeOH+",
        "FeOH+2",
        "HFeO2(aq)",
        "HFeO2-",
        "Ferropericlase",
        "Goethite",
        "Hematite",
        "Iron",
        "Magnetite",
    ]

    sys_p = rkt.ChemicalSystem(
        db,
        rkt.AqueousPhase(rkt.speciate(aq_p)),
        rkt.MineralPhases(rkt.StringList(min_p)),
    )
    specs_p = rkt.EquilibriumSpecs(sys_p)
    specs_p.temperature()
    specs_p.pressure()
    specs_p.pH()
    specs_p.Eh()
    st_p = rkt.ChemicalState(sys_p)
    st_p.set("H2O(aq)", ad.real(1.0), "kg")
    st_p.set("Fe+2", ad.real(1e-6), "mol")
    solver_p = rkt.EquilibriumSweepSolver(specs_p)

    pH_pts_p = np.array([-2.0, 1.0, 4.0, 7.0, 10.0, 13.0, 16.0], dtype=float)
    Eh_pts_p = np.array([-2.0, -1.0, 0.0, 1.0, 2.0], dtype=float)
    grid_p = solver_p.sweepPHEhGrid(st_p, pH_pts_p, Eh_pts_p, "V")

    out_p = os.path.join(BASE, "diagnostics_reaktoro_pourbaix_points.csv")
    write_case("Fe_Pourbaix", grid_p, pH_pts_p, Eh_pts_p, species_p, out_p)

    # --- Fe Mosaic (CHNOSZ-like) ---
    aq_m = [
        "H2O(aq)",
        "H+",
        "OH-",
        "Fe+2",
        "Fe+3",
        "HFeO2-",
        "SO4-2",
        "HSO4-",
        "HS-",
        "H2S(aq)",
        "CO3-2",
        "HCO3-",
        "CO2(aq)",
    ]
    min_m = ["Pyrite", "Pyrrhotite,trot", "Siderite", "Hematite", "Magnetite"]
    species_m = [
        "Fe+2",
        "Fe+3",
        "HFeO2-",
        "Pyrite",
        "Pyrrhotite,trot",
        "Siderite",
        "Hematite",
        "Magnetite",
    ]

    sys_m = rkt.ChemicalSystem(
        db,
        rkt.AqueousPhase(rkt.speciate(aq_m)),
        rkt.MineralPhases(rkt.StringList(min_m)),
    )
    specs_m = rkt.EquilibriumSpecs(sys_m)
    specs_m.temperature()
    specs_m.pressure()
    specs_m.pH()
    specs_m.Eh()
    st_m = rkt.ChemicalState(sys_m)
    st_m.set("H2O(aq)", ad.real(1.0), "kg")
    st_m.set("Fe+2", ad.real(1e-6), "mol")
    st_m.set("SO4-2", ad.real(1e-6), "mol")
    st_m.set("CO3-2", ad.real(1.0), "mol")
    solver_m = rkt.EquilibriumSweepSolver(specs_m)

    pH_pts_m = np.array([0.0, 2.0, 5.0, 8.0, 11.0, 14.0], dtype=float)
    Eh_pts_m = np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float)
    grid_m = solver_m.sweepPHEhGrid(st_m, pH_pts_m, Eh_pts_m, "V")

    out_m = os.path.join(BASE, "diagnostics_reaktoro_mosaic_points.csv")
    write_case("Fe_Mosaic", grid_m, pH_pts_m, Eh_pts_m, species_m, out_m)

    print("Wrote:", out_p)
    print("Wrote:", out_m)


if __name__ == "__main__":
    main()
