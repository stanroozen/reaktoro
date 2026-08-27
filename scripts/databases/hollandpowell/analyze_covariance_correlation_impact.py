#!/usr/bin/env python3
"""
Analyze covariance matrices exported by convert_tc_ds_to_reaktoro_json.py.

This script quantifies how correlations change generalized variance by comparing:
- Full covariance C
- Diagonal-only covariance D = diag(C)

Because tc-ds covariance matrices can be singular (det(C)=0), the script also
computes regularized comparisons using Cr = C + lambda*I and Dr = D + lambda*I.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import List, Tuple


def unpack_upper_triangle(packed: List[float], n: int) -> List[List[float]]:
    m = [[0.0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i, n):
            v = float(packed[k])
            m[i][j] = v
            m[j][i] = v
            k += 1
    return m


def lu_logdet_sign(a: List[List[float]], eps: float = 1e-300) -> Tuple[float, float]:
    """
    Compute sign(det(A)) and log(abs(det(A))) with partial-pivoting LU.
    Returns (0.0, -inf) for singular matrices.
    """
    n = len(a)
    sign = 1.0
    logabsdet = 0.0

    for k in range(n):
        piv = k
        piv_abs = abs(a[k][k])
        for r in range(k + 1, n):
            ar = abs(a[r][k])
            if ar > piv_abs:
                piv_abs = ar
                piv = r

        if piv_abs <= eps:
            return 0.0, float("-inf")

        if piv != k:
            a[k], a[piv] = a[piv], a[k]
            sign *= -1.0

        pivot = a[k][k]
        if pivot < 0.0:
            sign *= -1.0
        logabsdet += math.log(abs(pivot))

        inv_pivot = 1.0 / pivot
        for i in range(k + 1, n):
            factor = a[i][k] * inv_pivot
            if factor == 0.0:
                continue
            a[i][k] = 0.0
            row_i = a[i]
            row_k = a[k]
            for j in range(k + 1, n):
                row_i[j] -= factor * row_k[j]

    return sign, logabsdet


def correlation_stats(
    c: List[List[float]], diag: List[float]
) -> Tuple[float, float, float, int]:
    vals = []
    near_one = 0
    n = len(c)

    for i in range(n):
        di = diag[i]
        for j in range(i + 1, n):
            dj = diag[j]
            denom = math.sqrt(di * dj) if di > 0.0 and dj > 0.0 else 0.0
            corr = (c[i][j] / denom) if denom > 0.0 else 0.0
            a = abs(corr)
            vals.append(a)
            if a >= 0.999:
                near_one += 1

    vals.sort()
    m = len(vals)
    mean_abs = sum(vals) / m if m else 0.0
    p95 = vals[int(0.95 * (m - 1))] if m else 0.0
    max_abs = vals[-1] if m else 0.0
    return mean_abs, p95, max_abs, near_one


def analyze_covariance_file(path: Path, rel_lambdas: List[float]) -> str:
    data = json.loads(path.read_text(encoding="utf-8"))

    n = int(data["NumEntities"])
    packed = data["PackedUpperTriangle"]
    expected = n * (n + 1) // 2
    if len(packed) < expected:
        raise ValueError(
            f"{path.name}: packed length {len(packed)} < expected {expected}"
        )

    c = unpack_upper_triangle(packed[:expected], n)
    diag = [c[i][i] for i in range(n)]

    trace_c = sum(diag)
    pos = sum(1 for d in diag if d > 0.0)
    zero = sum(1 for d in diag if d == 0.0)
    neg = sum(1 for d in diag if d < 0.0)

    sign_raw, logdet_raw = lu_logdet_sign([row[:] for row in c])

    mean_abs_corr, p95_abs_corr, max_abs_corr, near_one = correlation_stats(c, diag)

    pos_vals = [d for d in diag if d > 0.0]
    mean_pos_var = sum(pos_vals) / len(pos_vals) if pos_vals else 0.0

    lines = []
    lines.append(f"=== {path.name} ===")
    lines.append(f"n = {n}")
    lines.append(f"trace(C) = {trace_c:.12g}")
    lines.append(f"diag counts: >0={pos}, =0={zero}, <0={neg}")
    lines.append(f"raw sign(det(C)) = {sign_raw}")
    lines.append(f"raw log(det(C)) = {logdet_raw:.12g}")
    lines.append(
        "corr stats: "
        f"mean|corr|={mean_abs_corr:.6f}, p95|corr|={p95_abs_corr:.6f}, "
        f"max|corr|={max_abs_corr:.6f}, |corr|>=0.999 pairs={near_one}"
    )

    if mean_pos_var <= 0.0:
        lines.append("No positive diagonal variances; cannot regularize.")
        return "\n".join(lines)

    lines.append("regularized determinant comparison vs diagonal-only baseline:")

    for rel in rel_lambdas:
        lam = rel * mean_pos_var
        cr = [row[:] for row in c]
        for i in range(n):
            cr[i][i] += lam

        sign_cr, logdet_cr = lu_logdet_sign(cr)

        logdet_dr = 0.0
        for d in diag:
            logdet_dr += math.log(max(d + lam, 1e-300))

        if sign_cr <= 0.0 or not math.isfinite(logdet_cr):
            lines.append(f"  rel={rel:.0e}: sign(det(Cr))={sign_cr}, cannot compare")
            continue

        log_ratio = logdet_cr - logdet_dr
        det_ratio = math.exp(log_ratio)
        vol_ratio = math.exp(0.5 * log_ratio)
        per_dim_std = math.exp(0.5 * log_ratio / n)
        orders = -log_ratio / math.log(10.0)

        lines.append(
            f"  rel={rel:.0e}: "
            f"log(det(Cr)/det(Dr))={log_ratio:.6f}, "
            f"det_ratio={det_ratio:.3e}, "
            f"det_reduction_orders={orders:.2f}, "
            f"volume_ratio={vol_ratio:.3e}, "
            f"per_dim_std_factor={per_dim_std:.5f}"
        )

    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Analyze covariance correlation impact (trace vs determinant, regularized comparisons)"
    )
    parser.add_argument(
        "files",
        nargs="*",
        help="Covariance JSON files. Default: tc-ds62-covariance.json and tc-ds636-covariance.json in this folder.",
    )
    parser.add_argument(
        "--regularization",
        nargs="*",
        type=float,
        default=[1e-12, 1e-10, 1e-8, 1e-6],
        help="Relative lambda values, where lambda = rel * mean(positive diagonal variance).",
    )
    args = parser.parse_args()

    here = Path(__file__).resolve().parent

    if args.files:
        paths = [
            Path(p) if Path(p).is_absolute() else (Path.cwd() / p).resolve()
            for p in args.files
        ]
    else:
        paths = [
            here / "tc-ds62-covariance.json",
            here / "tc-ds636-covariance.json",
        ]

    for path in paths:
        if not path.exists():
            print(f"Missing file: {path}")
            continue
        print(analyze_covariance_file(path, args.regularization))
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
