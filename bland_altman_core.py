#!/usr/bin/env python3
"""
Bland-Altman Method Comparison & Agreement Analyzer

Real implementation of:
  - Mean difference (bias): d̄ = Σ(di)/n where di = method1 - method2
  - Standard deviation of differences: SD = √(Σ(di-d̄)²/(n-1))
  - Limits of agreement: d̄ ± 1.96×SD
  - 95% CI for bias and limits
  - Percentage of points within limits
  - Proportional bias detection (regression of differences on means)
  - Intraclass Correlation Coefficient (ICC)

Pure Python stdlib — no external dependencies.
"""

import math
import csv
import json
import sys
from typing import List, Dict, Any, Optional, Tuple


# ---------------------------------------------------------------------------
# Core Bland-Altman analysis
# ---------------------------------------------------------------------------

def bland_altman(
    method1: List[float],
    method2: List[float],
) -> Dict[str, Any]:
    """
    Perform Bland-Altman analysis on two measurement methods.

    Parameters
    ----------
    method1 : list of float
        Measurements from method 1.
    method2 : list of float
        Measurements from method 2.

    Returns
    -------
    dict with keys:
        n               : number of paired observations
        mean_difference  : d̄ (bias)
        sd_differences   : SD of differences
        limits_of_agreement : (lower, upper) = d̄ ± 1.96×SD
        ci_bias          : 95% CI for bias (lower, upper)
        ci_lower_limit   : 95% CI for lower limit of agreement
        ci_upper_limit   : 95% CI for upper limit of agreement
        pct_within       : percentage of points within limits
        proportional_bias : dict with regression results
        differences      : list of individual differences
        means            : list of individual means
    """
    if len(method1) != len(method2):
        raise ValueError("method1 and method2 must have the same length")
    n = len(method1)
    if n < 2:
        raise ValueError("At least 2 paired observations are required")

    # Compute differences and means
    differences = [m1 - m2 for m1, m2 in zip(method1, method2)]
    means = [(m1 + m2) / 2.0 for m1, m2 in zip(method1, method2)]

    # Mean difference (bias)
    mean_diff = sum(differences) / n

    # Standard deviation of differences
    sd_diff = math.sqrt(sum((d - mean_diff) ** 2 for d in differences) / (n - 1))

    # Limits of agreement
    lower_limit = mean_diff - 1.96 * sd_diff
    upper_limit = mean_diff + 1.96 * sd_diff

    # 95% CI for bias
    se_bias = sd_diff / math.sqrt(n)
    ci_bias_lower = mean_diff - 1.96 * se_bias
    ci_bias_upper = mean_diff + 1.96 * se_bias

    # 95% CI for limits of agreement
    # SE of limits = SD * sqrt(3/n) (approximate)
    se_limit = sd_diff * math.sqrt(3.0 / n)
    ci_lower_lo = lower_limit - 1.96 * se_limit
    ci_lower_hi = lower_limit + 1.96 * se_limit
    ci_upper_lo = upper_limit - 1.96 * se_limit
    ci_upper_hi = upper_limit + 1.96 * se_limit

    # Percentage within limits
    within = sum(1 for d in differences if lower_limit <= d <= upper_limit)
    pct_within = (within / n) * 100.0

    # Proportional bias: regression of differences on means
    pb = _proportional_bias_test(means, differences)

    return {
        "n": n,
        "mean_difference": round(mean_diff, 6),
        "sd_differences": round(sd_diff, 6),
        "limits_of_agreement": (round(lower_limit, 6), round(upper_limit, 6)),
        "ci_bias": (round(ci_bias_lower, 6), round(ci_bias_upper, 6)),
        "ci_lower_limit": (round(ci_lower_lo, 6), round(ci_lower_hi, 6)),
        "ci_upper_limit": (round(ci_upper_lo, 6), round(ci_upper_hi, 6)),
        "pct_within": round(pct_within, 2),
        "proportional_bias": pb,
        "differences": differences,
        "means": means,
    }


def limits_of_agreement(method1: List[float], method2: List[float]) -> Tuple[float, float]:
    """Return the 95% limits of agreement."""
    result = bland_altman(method1, method2)
    return result["limits_of_agreement"]


def mean_bias(method1: List[float], method2: List[float]) -> float:
    """Return the mean difference (bias)."""
    result = bland_altman(method1, method2)
    return result["mean_difference"]


# ---------------------------------------------------------------------------
# Intraclass Correlation Coefficient (ICC)
# ---------------------------------------------------------------------------

def intraclass_correlation(
    method1: List[float],
    method2: List[float],
    model: str = "two_way_random",
) -> Dict[str, Any]:
    """
    Compute the Intraclass Correlation Coefficient (ICC).

    Uses a two-way random effects model (ICC(2,1) - single measures,
    absolute agreement).

    ICC = (MSR - MSE) / (MSR + (k-1)*MSE + k*(MSB-MSE)/n)

    Where:
    - MSR = mean square for rows (subjects)
    - MSE = mean square for error
    - MSB = mean square for columns (methods)
    - k = number of methods (2)
    - n = number of subjects

    For the simpler case with 2 methods:
    ICC = (BMS - WMS) / (BMS + WMS)
    where BMS = between-subject MS, WMS = within-subject MS
    """
    if len(method1) != len(method2):
        raise ValueError("method1 and method2 must have the same length")
    n = len(method1)
    if n < 2:
        raise ValueError("At least 2 paired observations are required")

    k = 2  # number of methods

    # Grand mean
    all_vals = method1 + method2
    grand_mean = sum(all_vals) / (2 * n)

    # Subject means
    subject_means = [(m1 + m2) / 2.0 for m1, m2 in zip(method1, method2)]

    # Between-subject sum of squares
    ss_between = k * sum((sm - grand_mean) ** 2 for sm in subject_means)

    # Within-subject sum of squares
    ss_within = sum((m1 - sm) ** 2 + (m2 - sm) ** 2
                    for m1, m2, sm in zip(method1, method2, subject_means))

    # Method sum of squares
    mean_m1 = sum(method1) / n
    mean_m2 = sum(method2) / n
    ss_method = n * ((mean_m1 - grand_mean) ** 2 + (mean_m2 - grand_mean) ** 2)

    # Error sum of squares (within - method)
    ss_error = ss_within - ss_method

    # Mean squares
    df_between = n - 1
    df_within = n * (k - 1)
    df_method = k - 1
    df_error = (n - 1) * (k - 1)

    ms_between = ss_between / df_between if df_between > 0 else 0
    ms_within = ss_within / df_within if df_within > 0 else 0
    ms_error = ss_error / df_error if df_error > 0 else 0

    # ICC(2,1) - Two-way random, single measures, absolute agreement
    if (ms_between + ms_error) == 0:
        icc = 0.0
    else:
        icc = (ms_between - ms_error) / (ms_between + ms_error)

    # 95% CI for ICC (approximate using F-distribution)
    # F = MS_between / MS_error
    if ms_error > 0:
        f_stat = ms_between / ms_error
    else:
        f_stat = float('inf')

    # Approximate CI
    alpha = 0.05
    # Lower bound
    if f_stat > 0:
        fl = f_stat / _f_inv(1 - alpha / 2, df_between, df_error) if df_error > 0 else f_stat
        fu = f_stat * _f_inv(1 - alpha / 2, df_error, df_between) if df_between > 0 else f_stat
    else:
        fl = 0.0
        fu = 0.0

    icc_lower = (fl - 1) / (fl + 1) if fl > 0 else 0.0
    icc_upper = (fu - 1) / (fu + 1) if fu > 0 else 1.0

    return {
        "icc": round(icc, 6),
        "ci_lower": round(max(-1, icc_lower), 6),
        "ci_upper": round(min(1, icc_upper), 6),
        "ms_between": round(ms_between, 6),
        "ms_within": round(ms_within, 6),
        "ms_error": round(ms_error, 6),
        "f_statistic": round(f_stat, 4),
        "n_subjects": n,
        "n_methods": k,
    }


# ---------------------------------------------------------------------------
# Proportional bias detection
# ---------------------------------------------------------------------------

def _proportional_bias_test(
    means: List[float],
    differences: List[float],
) -> Dict[str, Any]:
    """
    Test for proportional bias by regressing differences on means.

    If the slope is significantly different from zero, there is
    proportional bias (the methods disagree more as the measured
    value increases).
    """
    n = len(means)
    if n < 3:
        return {"slope": 0.0, "intercept": 0.0, "r_squared": 0.0, "p_value": 1.0, "proportional_bias": False}

    # Simple linear regression: differences = a + b * means
    mean_x = sum(means) / n
    mean_y = sum(differences) / n

    ss_xy = sum((means[i] - mean_x) * (differences[i] - mean_y) for i in range(n))
    ss_xx = sum((means[i] - mean_x) ** 2 for i in range(n))
    ss_yy = sum((differences[i] - mean_y) ** 2 for i in range(n))

    if ss_xx == 0:
        return {"slope": 0.0, "intercept": mean_y, "r_squared": 0.0, "p_value": 1.0, "proportional_bias": False}

    slope = ss_xy / ss_xx
    intercept = mean_y - slope * mean_x

    # R-squared
    if ss_yy == 0:
        r_squared = 0.0
    else:
        r_squared = (ss_xy ** 2) / (ss_xx * ss_yy)

    # t-test for slope significance
    if n < 4:
        return {"slope": slope, "intercept": intercept, "r_squared": r_squared, "p_value": 1.0, "proportional_bias": False}

    # Residual standard error
    ss_res = ss_yy - slope * ss_xy
    mse = ss_res / (n - 2)
    se_slope = math.sqrt(mse / ss_xx) if ss_xx > 0 and mse > 0 else 0.0

    t_stat = slope / se_slope if se_slope > 0 else 0.0
    # p-value from t-distribution (two-tailed)
    p_value = 2.0 * _t_sf(abs(t_stat), n - 2)

    return {
        "slope": round(slope, 6),
        "intercept": round(intercept, 6),
        "r_squared": round(r_squared, 6),
        "t_statistic": round(t_stat, 4),
        "p_value": round(p_value, 6),
        "proportional_bias": p_value < 0.05,
    }


# ---------------------------------------------------------------------------
# CSV batch processing
# ---------------------------------------------------------------------------

def process_csv(input_path: str, output_path: str) -> Dict[str, Any]:
    """Process a CSV file for Bland-Altman analysis."""
    with open(input_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    # Find method columns
    method1_col = _find_col(fieldnames, ["method1", "m1", "test", "new", "device1", "v1"])
    method2_col = _find_col(fieldnames, ["method2", "m2", "reference", "gold", "device2", "v2"])

    method1 = [float(r[method1_col]) for r in rows]
    method2 = [float(r[method2_col]) for r in rows]

    result = bland_altman(method1, method2)
    icc_result = intraclass_correlation(method1, method2)

    # Write detailed output
    out_fields = ["subject", "method1", "method2", "difference", "mean", "within_limits"]
    out_rows = []
    lower, upper = result["limits_of_agreement"]
    for i, r in enumerate(rows):
        out_rows.append({
            "subject": str(i),
            "method1": str(method1[i]),
            "method2": str(method2[i]),
            "difference": str(round(result["differences"][i], 6)),
            "mean": str(round(result["means"][i], 6)),
            "within_limits": str(lower <= result["differences"][i] <= upper),
        })

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    return {**result, "icc": icc_result}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_col(fieldnames: List[str], candidates: List[str]) -> str:
    lower_map = {c.lower(): c for c in fieldnames}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return fieldnames[0]


def _t_sf(t: float, df: int) -> float:
    """Survival function of t-distribution."""
    if df <= 0:
        return 1.0
    if df > 30:
        return 2.0 * (1.0 - _norm_cdf(abs(t)))
    x = df / (df + t * t)
    return _inc_beta(df / 2.0, 0.5, x)


def _norm_cdf(x: float) -> float:
    """Standard normal CDF."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _f_inv(p: float, df1: int, df2: int) -> float:
    """
    Approximate inverse F-distribution (quantile function).
    Uses the normal approximation for large df.
    """
    if df1 <= 0 or df2 <= 0:
        return 1.0

    # For large df, use normal approximation
    # F ≈ exp(2 * z * sqrt(1/a)) where a = min(df1, df2)
    z = _norm_ppf(p)
    a = min(df1, df2)

    # Simple approximation
    if a > 30:
        return math.exp(2 * z / math.sqrt(a))

    # Use a rough iterative approach
    # For moderate df, use the relationship with chi-square
    # This is a simplified approximation
    return 1.0 + 2.0 * z / math.sqrt(a) + z * z / a


def _norm_ppf(p: float) -> float:
    """Inverse standard normal CDF."""
    if p <= 0:
        return float('-inf')
    if p >= 1:
        return float('inf')
    if p == 0.5:
        return 0.0

    if p < 0.5:
        return -_rational_approx(math.sqrt(-2.0 * math.log(p)))
    else:
        return _rational_approx(math.sqrt(-2.0 * math.log(1.0 - p)))


def _rational_approx(t: float) -> float:
    c0, c1, c2 = 2.515517, 0.802853, 0.010328
    d1, d2, d3 = 1.432788, 0.189269, 0.001308
    return t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t)


def _inc_beta(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta function."""
    if x < 0 or x > 1:
        return 0.0
    if x == 0:
        return 0.0
    if x == 1:
        return 1.0

    max_iter = 200
    eps = 1e-12

    lbeta = _log_gamma(a) + _log_gamma(b) - _log_gamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1.0 - x) * b - lbeta)

    if x < (a + 1) / (a + b + 2):
        cf = _beta_cf(a, b, x, max_iter, eps)
        return front * cf / a
    else:
        cf = _beta_cf(b, a, 1.0 - x, max_iter, eps)
        return 1.0 - front * cf / b


def _beta_cf(a: float, b: float, x: float, max_iter: int, eps: float) -> float:
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0

    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d

    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c

        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta

        if abs(delta - 1.0) < eps:
            break

    return h


def _log_gamma(x: float) -> float:
    if x < 0.5:
        return math.log(math.pi / math.sin(math.pi * x)) - _log_gamma(1.0 - x)
    x -= 1.0
    g = 7
    c = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ]
    total = c[0]
    for i in range(1, g + 2):
        total += c[i] / (x + i)
    st = x + g + 0.5
    return 0.5 * math.log(2.0 * math.pi) + (x + 0.5) * math.log(st) - st + math.log(total)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def summary(method1: List[float], method2: List[float]) -> str:
    """Return a formatted summary of the Bland-Altman analysis."""
    result = bland_altman(method1, method2)
    icc = intraclass_correlation(method1, method2)

    lines = []
    lines.append("Bland-Altman Agreement Analysis")
    lines.append("=" * 60)
    lines.append(f"  N pairs:           {result['n']}")
    lines.append(f"  Mean difference:   {result['mean_difference']:.4f}")
    lines.append(f"  SD of differences: {result['sd_differences']:.4f}")
    lines.append("")
    lo, hi = result["limits_of_agreement"]
    lines.append(f"  Limits of Agreement:")
    lines.append(f"    Lower: {lo:.4f}  95% CI: [{result['ci_lower_limit'][0]:.4f}, {result['ci_lower_limit'][1]:.4f}]")
    lines.append(f"    Upper: {hi:.4f}  95% CI: [{result['ci_upper_limit'][0]:.4f}, {result['ci_upper_limit'][1]:.4f}]")
    lines.append("")
    lines.append(f"  Bias 95% CI:       [{result['ci_bias'][0]:.4f}, {result['ci_bias'][1]:.4f}]")
    lines.append(f"  Within limits:     {result['pct_within']:.1f}%")
    lines.append("")
    lines.append(f"  ICC:               {icc['icc']:.4f}  95% CI: [{icc['ci_lower']:.4f}, {icc['ci_upper']:.4f}]")
    lines.append("")
    pb = result["proportional_bias"]
    lines.append(f"  Proportional Bias:")
    lines.append(f"    Slope:     {pb['slope']:.6f}")
    lines.append(f"    R-squared: {pb['r_squared']:.6f}")
    lines.append(f"    p-value:   {pb['p_value']:.6f}")
    lines.append(f"    Detected:  {'Yes' if pb['proportional_bias'] else 'No'}")

    return "\n".join(lines)
