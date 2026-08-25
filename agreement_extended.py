#!/usr/bin/env python3
"""
Bland-Altman Agreement — Extended Statistics
Proportional-bias regression (difference on mean), Lin's CCC, Gwet's AC1,
repeatability coefficient, and bootstrap CIs for bias/LoA.

Zero-dependency. Author: Dr. Abu Suraih Sakhri. License: MIT.
"""
import argparse
import json
import math
import random
import statistics
import sys
from dataclasses import dataclass
from typing import Dict, Any, List, Sequence, Tuple


@dataclass
class MethodPair:
    subject_id: str
    method_a: float
    method_b: float


def basic_bland_altman(pairs: List[MethodPair]) -> Dict[str, Any]:
    diffs = [p.method_a - p.method_b for p in pairs]
    means = [(p.method_a + p.method_b) / 2 for p in pairs]
    n = len(diffs)
    bias = statistics.mean(diffs)
    sd = statistics.stdev(diffs) if n > 1 else 0.0
    loa = (bias - 1.96 * sd, bias + 1.96 * sd)
    return {"n": n, "bias": round(bias, 4), "sd_diff": round(sd, 4),
            "lower_loa": round(loa[0], 4), "upper_loa": round(loa[1], 4),
            "mean_of_means": round(statistics.mean(means), 4)}


def proportional_bias_regression(pairs: List[MethodPair]) -> Dict[str, Any]:
    """Regress difference on average: d = a + b*avg; test H0: b = 0 (t-test)."""
    xs = [(p.method_a + p.method_b) / 2 for p in pairs]
    ys = [p.method_a - p.method_b for p in pairs]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx if sxx else 0.0
    intercept = my - slope * mx
    resid = [y - (intercept + slope * x) for x, y in zip(xs, ys)]
    s_err = math.sqrt(sum(r * r for r in resid) / max(n - 2, 1))
    se_slope = s_err / math.sqrt(sxx) if sxx else float("inf")
    t = slope / se_slope if se_slope else 0.0
    # two-sided p via t-distribution survival approximation (normal for n>30,
    # exact small-sample via incomplete beta is overkill here -> report t and df)
    return {
        "intercept": round(intercept, 4),
        "slope": round(slope, 5),
        "se_slope": f"{se_slope:.3e}",
        "t_statistic": round(t, 3),
        "df": n - 2,
        "proportional_bias_detected": abs(t) > _tcrit_95(n - 2),
        "interpretation": ("difference scales with measurement level: use percentage "
                           "limits of agreement" if abs(t) > _tcrit_95(n - 2)
                           else "no significant proportional bias"),
    }


def _tcrit_95(df: int) -> float:
    """Approximate two-sided 95% critical t."""
    if df >= 30:
        return 1.96
    table = {2: 12.71, 3: 4.30, 4: 3.18, 5: 2.78, 6: 2.57, 8: 2.36, 10: 2.23,
             15: 2.13, 20: 2.09, 25: 2.06}
    keys = sorted(table.keys())
    for k in keys:
        if df <= k:
            return table[k]
    return 2.09


def lins_ccc(pairs: List[MethodPair]) -> Dict[str, Any]:
    """Lin's concordance correlation coefficient."""
    xa = [p.method_a for p in pairs]
    xb = [p.method_b for p in pairs]
    n = len(xa)
    ma, mb = sum(xa) / n, sum(xb) / n
    va = sum((x - ma) ** 2 for x in xa) / n
    vb = sum((x - mb) ** 2 for x in xb) / n
    cov = sum((a - ma) * (b - mb) for a, b in zip(xa, xb)) / n
    pearson = cov / math.sqrt(va * vb) if va and vb else 0.0
    ccc = (2 * cov) / ((ma - mb) ** 2 + va + vb) if (va + vb) else 0.0
    return {
        "ccc": round(ccc, 4),
        "pearson_r": round(pearson, 4),
        "precision_vs_accuracy": {
            "pearson_measures_precision_only": round(pearson, 4),
            "ccc_penalizes_scale_and_offset_error": round(ccc, 4),
            "accuracy_shift_penalty": round(abs(ma - mb), 4),
        },
        "interpretation": _ccc_interpretation(ccc),
    }


def _ccc_interpretation(c: float) -> str:
    if c >= 0.99: return "almost perfect"
    if c >= 0.95: return "substantial"
    if c >= 0.90: return "moderate"
    if c >= 0.80: return "fair agreement; review LoA against clinical tolerance"
    return "poor"


def gwet_ac1(labels_a: Sequence[str], labels_b: Sequence[str]) -> Dict[str, Any]:
    """Gwet's AC1 first-order agreement coefficient (robust to prevalence)."""
    if len(labels_a) != len(labels_b) or not labels_a:
        return {"error": "label length mismatch"}
    cats = sorted(set(labels_a) | set(labels_b))
    k = len(cats)
    n = len(labels_a)
    pa = sum(1 for a, b in zip(labels_a, labels_b) if a == b) / n
    marg = {c: (sum(1 for x in labels_a if x == c) +
                sum(1 for x in labels_b if x == c)) / (2 * n) for c in cats}
    pe = sum(p * (1 - p) for p in marg.values()) / max(k - 1, 1)
    ac1 = (pa - pe) / (1 - pe) if pe < 1 else 0.0
    kappa = None
    pk = sum(marg[c] ** 2 for c in cats)
    if pk < 1:
        kappa = (pa - pk) / (1 - pk)
    return {"observed_agreement": round(pa, 4), "chance_expected_ac1": round(pe, 4),
            "gwet_ac1": round(ac1, 4), "cohens_kappa_for_comparison":
                round(kappa, 4) if kappa is not None else None}


def repeatability_coefficient(repeated_subjects: Dict[str, List[float]]) -> float:
    """Within-subject SD across replicates -> r = 1.96*sqrt(2)*Sw."""
    ssq = []
    for _, vals in repeated_subjects.items():
        if len(vals) >= 2:
            m = sum(vals) / len(vals)
            ssq.append(sum((v - m) ** 2 for v in vals) / (len(vals) - 1))
    sw = math.sqrt(sum(ssq) / len(ssq)) if ssq else 0.0
    return round(1.96 * math.sqrt(2.0) * sw, 4)


def bootstrap_intervals(pairs: List[MethodPair], iterations: int = 2000,
                        seed: int = 42) -> Dict[str, Any]:
    rng = random.Random(seed)
    biases, lows, highs = [], [], []
    n = len(pairs)
    for _ in range(iterations):
        sample = [pairs[rng.randrange(n)] for _ in range(n)]
        diffs = [p.method_a - p.method_b for p in sample]
        b = statistics.mean(diffs)
        s = statistics.stdev(diffs) if n > 1 else 0.0
        biases.append(b)
        lows.append(b - 1.96 * s)
        highs.append(b + 1.96 * s)

    def pct(v: List[float], q: float) -> float:
        s = sorted(v)
        return round(s[min(len(s) - 1, int(q * len(s)))], 4)

    return {
        "iterations": iterations,
        "bias_95ci_percentile": [pct(biases, 0.025), pct(biases, 0.975)],
        "loa_lower_95ci": [pct(lows, 0.025), pct(lows, 0.975)],
        "loa_upper_95ci": [pct(highs, 0.025), pct(highs, 0.975)],
        "method": "percentile bootstrap",
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Extended Bland-Altman demo")
    ap.add_argument("--iters", type=int, default=2000)
    args = ap.parse_args()

    # synthetic method-comparison with proportional bias: B under-reads at high values
    pairs = []
    for i, true_val in enumerate([4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26]):
        a = true_val + (0.15 * i % 3) - 0.1
        b = true_val - 0.02 * true_val + 0.08
        pairs.append(MethodPair(f"S{i+1}", a, b))

    print(json.dumps(basic_bland_altman(pairs), indent=2))
    print(json.dumps(proportional_bias_regression(pairs), indent=2))
    print(json.dumps(lins_ccc(pairs), indent=2))
    print(json.dumps(gwet_ac1(["pos", "neg", "pos", "neg", "pos"],
                              ["pos", "neg", "pos", "pos", "pos"]), indent=2))
    reps = {"A": [5.1, 5.4], "B": [7.0, 6.8], "C": [9.2, 9.9]}
    print(json.dumps({"repeatability_coefficient": repeatability_coefficient(reps)}, indent=2))
    print(json.dumps(bootstrap_intervals(pairs, iterations=args.iters), indent=2))
