import csv
import os
import re

BASE = os.path.dirname(os.path.abspath(__file__))

CASES = {
    "Fe_Pourbaix": {
        "reaktoro": [
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
        ],
        "chnosz": [
            "Fe+2",
            "Fe+3",
            "FeOH+",
            "FeOH+2",
            "HFeO2-",
            "HFeO2",
            "FeO+",
            "FeO2-",
            "hematite",
            "magnetite",
            "goethite",
            "iron",
        ],
    },
    "Fe_Mosaic": {
        "reaktoro": [
            "Fe+2",
            "Fe+3",
            "HFeO2-",
            "Pyrite",
            "Pyrrhotite,trot",
            "Siderite",
            "Hematite",
            "Magnetite",
        ],
        "chnosz": [
            "Fe+2",
            "Fe+3",
            "HFeO2-",
            "pyrite",
            "pyrrhotite",
            "siderite",
            "hematite",
            "magnetite",
        ],
    },
}


def canon_exact(name: str) -> str:
    return name.strip().lower()


def canon_approx(name: str) -> str:
    n = name.strip().lower()
    n = n.replace("(aq)", "")
    n = n.replace(",trot", "")
    n = re.sub(r"[^a-z0-9+-]", "", n)
    return n


def build_mapping(case_name: str, rk_list, ch_list):
    ch_exact = {canon_exact(x): x for x in ch_list}
    ch_approx = {canon_approx(x): x for x in ch_list}

    rows = []
    used_ch = set()

    for rk in rk_list:
        rk_e = canon_exact(rk)
        rk_a = canon_approx(rk)

        if rk_e in ch_exact:
            ch = ch_exact[rk_e]
            status = "exact"
            reason = "case-insensitive exact name match"
        elif rk_a in ch_approx:
            ch = ch_approx[rk_a]
            status = "approx"
            reason = "normalized alias (aq polymorph/comma/punctuation removed)"
        else:
            ch = ""
            status = "missing_in_chnosz"
            reason = "no CHNOSZ species with equivalent normalized key"

        if ch:
            used_ch.add(ch)

        rows.append(
            {
                "case": case_name,
                "reaktoro_name": rk,
                "chnosz_name": ch,
                "status": status,
                "reason": reason,
            }
        )

    for ch in ch_list:
        if ch not in used_ch:
            rows.append(
                {
                    "case": case_name,
                    "reaktoro_name": "",
                    "chnosz_name": ch,
                    "status": "missing_in_reaktoro",
                    "reason": "present in CHNOSZ list but unmatched from Reaktoro side",
                }
            )

    return rows


def write_outputs(rows):
    csv_path = os.path.join(BASE, "species_mapping_table.csv")
    md_path = os.path.join(BASE, "species_mapping_table.md")

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["case", "reaktoro_name", "chnosz_name", "status", "reason"],
        )
        w.writeheader()
        w.writerows(rows)

    by_case = {}
    for r in rows:
        by_case.setdefault(r["case"], []).append(r)

    lines = ["# Species Mapping: Reaktoro vs CHNOSZ", ""]
    for case, recs in by_case.items():
        lines.append(f"## {case}")
        counts = {
            "exact": sum(1 for x in recs if x["status"] == "exact"),
            "approx": sum(1 for x in recs if x["status"] == "approx"),
            "missing_in_chnosz": sum(
                1 for x in recs if x["status"] == "missing_in_chnosz"
            ),
            "missing_in_reaktoro": sum(
                1 for x in recs if x["status"] == "missing_in_reaktoro"
            ),
        }
        lines.append(
            f"- exact: {counts['exact']}, approx: {counts['approx']}, missing_in_chnosz: {counts['missing_in_chnosz']}, missing_in_reaktoro: {counts['missing_in_reaktoro']}"
        )
        lines.append("")
        lines.append("| Reaktoro | CHNOSZ | Status | Reason |")
        lines.append("|---|---|---|---|")
        for r in recs:
            lines.append(
                f"| {r['reaktoro_name']} | {r['chnosz_name']} | {r['status']} | {r['reason']} |"
            )
        lines.append("")

    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return csv_path, md_path


def main():
    all_rows = []
    for case_name, d in CASES.items():
        all_rows.extend(build_mapping(case_name, d["reaktoro"], d["chnosz"]))

    csv_path, md_path = write_outputs(all_rows)
    print("Wrote:", csv_path)
    print("Wrote:", md_path)


if __name__ == "__main__":
    main()
