"""Bland-Altman Method Comparison: agreement analysis with bias, limits, and regression."""
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
import math


@dataclass
class MethodPair:
    sample_id: str
    method_a: float
    method_b: float


class BlandAltmanAnalyzer:
    """Bland-Altman analysis for method comparison in clinical assays."""

    def __init__(self):
        self._pairs: List[MethodPair] = []

    def add_pair(self, pair: MethodPair) -> None:
        self._pairs.append(pair)

    def analyze(self) -> Dict[str, Any]:
        """Full Bland-Altman analysis."""
        if len(self._pairs) < 2:
            return {"status": "insufficient_data", "min_required": 2}

        means = [(p.method_a + p.method_b) / 2 for p in self._pairs]
        diffs = [p.method_a - p.method_b for p in self._pairs]

        n = len(diffs)
        mean_diff = sum(diffs) / n
        std_diff = math.sqrt(sum((d - mean_diff) ** 2 for d in diffs) / max(n - 1, 1))

        ucl = mean_diff + 1.96 * std_diff
        lcl = mean_diff - 1.96 * std_diff

        within_limits = sum(1 for d in diffs if lcl <= d <= ucl)
        percent_within = within_limits / n * 100

        proportional_bias = self._test_proportional_bias(means, diffs)
        constant_bias = abs(mean_diff) > 0.1 * std_diff

        return {
            "num_pairs": n,
            "bias": round(mean_diff, 4),
            "bias_std": round(std_diff, 4),
            "upper_limit_of_agreement": round(ucl, 4),
            "lower_limit_of_agreement": round(lcl, 4),
            "percent_within_limits": round(percent_within, 1),
            "proportional_bias_detected": proportional_bias,
            "constant_bias_detected": constant_bias,
            "acceptable_agreement": percent_within >= 95 and not constant_bias,
        }

    def _test_proportional_bias(self, means: List[float], diffs: List[float]) -> bool:
        n = len(means)
        if n < 3:
            return False
        mean_m = sum(means) / n
        mean_d = sum(diffs) / n
        cov = sum((m - mean_m) * (d - mean_d) for m, d in zip(means, diffs)) / max(n - 1, 1)
        var_m = sum((m - mean_m) ** 2 for m in means) / max(n - 1, 1)
        if var_m == 0:
            return False
        slope = cov / var_m
        return abs(slope) > 0.1

    def trend_analysis(self) -> Dict[str, Any]:
        """Analyze trends in differences across the measurement range."""
        if len(self._pairs) < 3:
            return {"status": "insufficient_data"}

        means = [(p.method_a + p.method_b) / 2 for p in self._pairs]
        diffs = [p.method_a - p.method_b for p in self._pairs]

        sorted_pairs = sorted(zip(means, diffs), key=lambda x: x[0])
        third = len(sorted_pairs) // 3

        low_mean = sum(d for _, d in sorted_pairs[:third]) / max(third, 1)
        mid_mean = sum(d for _, d in sorted_pairs[third:2*third]) / max(third, 1)
        high_mean = sum(d for _, d in sorted_pairs[2*third:]) / max(len(sorted_pairs) - 2*third, 1)

        trend = "stable"
        if high_mean - low_mean > 0.5 * abs(low_mean if low_mean != 0 else 1):
            trend = "increasing_difference"
        elif low_mean - high_mean > 0.5 * abs(high_mean if high_mean != 0 else 1):
            trend = "decreasing_difference"

        return {
            "low_range_bias": round(low_mean, 4),
            "mid_range_bias": round(mid_mean, 4),
            "high_range_bias": round(high_mean, 4),
            "trend": trend,
        }
