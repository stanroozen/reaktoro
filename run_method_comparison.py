#!/usr/bin/env python3
"""Test all 4 integration methods systematically"""

import subprocess, re
from pathlib import Path
import time


class Tester:
    def __init__(self):
        self.hpp = Path("Reaktoro/Extensions/DEW/WaterGibbsModel.hpp")
        self.exe = Path("build-msvc/Reaktoro/Release/reaktoro-cpptests.exe")

    def set_method(self, name, enum_val):
        with open(self.hpp, "r", encoding="utf-8") as f:
            c = f.read()
        new_c = re.sub(r"(integrationMethod = )\w+;", f"\\1{enum_val};", c)
        with open(self.hpp, "w", encoding="utf-8") as f:
            f.write(new_c)
        print(f"[+] {name}")
        return True

    def rebuild(self):
        r = subprocess.run(
            [
                "cmake",
                "--build",
                "build-msvc",
                "--config",
                "Release",
                "--target",
                "reaktoro-cpptests",
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )
        ok = r.returncode == 0
        print("[+] Build" if ok else "[x] Build FAIL")
        return ok

    def test(self):
        r = subprocess.run(
            [str(self.exe), "[dew][reaction]"],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=".",
        )
        out = r.stdout + r.stderr
        res = {}

        if "180 passed" in out:
            res["passed"] = 180
            res["status"] = "PASS"

        m = re.search(r"avg=([0-9.]+)\s*J/mol", out)
        res["error"] = float(m.group(1)) if m else 0

        for line in out.split("\n"):
            if "Tested" in line or "passed" in line:
                print(f"  {line.strip()}")

        return res

    def run_all(self):
        methods = [
            ("Trapezoidal", "Trapezoidal"),
            ("Simpson", "Simpson"),
            ("Gauss-Legendre-16", "GaussLegendre16"),
        ]

        print("\nMETHOD COMPARISON\n" + "=" * 60)
        results = []

        for name, enum_val in methods:
            print(f"\nTesting: {name}")
            print("-" * 60)
            self.set_method(name, enum_val)
            if self.rebuild():
                r = self.test()
                r["name"] = name
                results.append(r)
            time.sleep(1)

        print("\n" + "=" * 60)
        print("RESULTS")
        print("=" * 60)
        for r in results:
            print(
                f"{r['name']:<20} | {r.get('status', '?'):<6} | Error: {r.get('error', 0):.4f} J/mol"
            )
        print("=" * 60)

        # Restore
        with open(self.hpp, "r", encoding="utf-8") as f:
            c = f.read()
        with open(self.hpp, "w", encoding="utf-8") as f:
            f.write(re.sub(r"(integrationMethod = )\w+;", r"\1Trapezoidal;", c))
        print("[+] Restored")


if __name__ == "__main__":
    Tester().run_all()
