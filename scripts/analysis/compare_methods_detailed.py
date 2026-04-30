#!/usr/bin/env python3
"""
Detailed Performance Comparison: All 4 Integration Methods
Measures timing, accuracy, and convergence characteristics
"""

import subprocess, re, time
from pathlib import Path


class DetailedComparison:
    def __init__(self):
        self.hpp = Path("Reaktoro/Extensions/DEW/WaterGibbsModel.hpp")
        self.exe = Path("build-msvc/Reaktoro/Release/reaktoro-cpptests.exe")

    def set_method(self, name, enum_val):
        with open(self.hpp, "r", encoding="utf-8") as f:
            c = f.read()
        new_c = re.sub(r"(integrationMethod = )\w+;", f"\\1{enum_val};", c)
        with open(self.hpp, "w", encoding="utf-8") as f:
            f.write(new_c)

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
        return r.returncode == 0

    def test(self):
        t0 = time.time()
        r = subprocess.run(
            [str(self.exe), "[dew][reaction]"],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=".",
        )
        elapsed = time.time() - t0

        out = r.stdout + r.stderr
        res = {
            "passed": 0,
            "status": "FAIL",
            "error_min": 0,
            "error_max": 0,
            "error_avg": 0,
            "error_rel": 0,
            "time_sec": elapsed,
        }

        if "180 passed" in out:
            res["passed"] = 180
            res["status"] = "PASS"

        # Extract error min/max/avg
        m = re.search(r"min=([0-9.]+),\s*max=([0-9.]+),\s*avg=([0-9.]+)", out)
        if m:
            res["error_min"] = float(m.group(1))
            res["error_max"] = float(m.group(2))
            res["error_avg"] = float(m.group(3))

        # Extract relative error
        m = re.search(r"relative error:.*?avg=([0-9.]+)%", out)
        if m:
            res["error_rel"] = float(m.group(1))

        return res

    def compare_all(self):
        methods = [
            ("Trapezoidal (O(h2))", "Trapezoidal", "5000 fixed steps"),
            ("Simpson (O(h4))", "Simpson", "5000 steps, parabolic"),
            (
                "Gauss-Legendre-16 (O(1/n32))",
                "GaussLegendre16",
                "312 segments x 16 nodes",
            ),
            ("Adaptive Simpson", "AdaptiveSimpson", "Auto-subdivide to 0.1 J/mol"),
        ]

        print("\n" + "=" * 100)
        print("DETAILED PERFORMANCE COMPARISON - ALL 4 INTEGRATION METHODS")
        print("=" * 100)
        print()

        results = []

        for name, enum_val, desc in methods:
            print(f"\nTesting: {name}")
            print(f"  Description: {desc}")
            print("-" * 100)

            self.set_method(name, enum_val)

            if self.rebuild():
                print("  [+] Build OK, running tests...")
                r = self.test()
                r["name"] = name
                results.append(r)
                print(f"  [+] Time: {r['time_sec']:.2f}s")
                print(
                    f"  [+] Error: min={r['error_min']:.2f}, max={r['error_max']:.2f}, avg={r['error_avg']:.4f} J/mol"
                )
                print(f"  [+] Relative error: {r['error_rel']:.4f}%")
            else:
                print("  [x] Build failed")

            time.sleep(1)

        # Detailed comparison table
        print("\n" + "=" * 100)
        print("COMPARISON TABLE")
        print("=" * 100)
        print()
        print(
            f"{'Method':<30} {'Status':<8} {'Min Err':<12} {'Max Err':<12} {'Avg Err':<12} {'Rel %':<10} {'Time (s)':<10}"
        )
        print("-" * 100)

        for r in results:
            print(
                f"{r['name']:<30} {r['status']:<8} {r['error_min']:>10.2f} J/mol {r['error_max']:>10.2f} J/mol {r['error_avg']:>10.4f} J/mol {r['error_rel']:>8.4f}% {r['time_sec']:>8.2f}s"
            )

        # Analysis
        print("\n" + "=" * 100)
        print("ANALYSIS")
        print("=" * 100)
        print()

        if len(results) > 1:
            baseline = results[0]
            print(
                f"Baseline (Trapezoidal): {baseline['error_avg']:.4f} J/mol, {baseline['time_sec']:.2f}s\n"
            )

            for i, r in enumerate(results[1:], 1):
                err_improvement = (
                    (
                        (baseline["error_avg"] - r["error_avg"])
                        / baseline["error_avg"]
                        * 100
                    )
                    if baseline["error_avg"] > 0
                    else 0
                )
                time_ratio = r["time_sec"] / baseline["time_sec"]
                efficiency = (
                    (err_improvement / (time_ratio * 100)) if time_ratio > 0 else 0
                )

                print(f"{r['name']}:")
                print(f"  Error improvement: {err_improvement:+.1f}%")
                print(f"  Relative time: {time_ratio:.2f}x")
                if abs(err_improvement) > 0.1:
                    print(f"  Efficiency (improvement/cost): {efficiency:.3f}")
                print()

        print("=" * 100)

        # Restore default
        with open(self.hpp, "r", encoding="utf-8") as f:
            c = f.read()
        with open(self.hpp, "w", encoding="utf-8") as f:
            f.write(re.sub(r"(integrationMethod = )\w+;", r"\1Trapezoidal;", c))
        print("\n[+] Restored default to Trapezoidal")


if __name__ == "__main__":
    comp = DetailedComparison()
    comp.compare_all()
