"""Bland-Altman Statistical Tests: regression on differences, variance ratio, and agreement indices."""
from typing import Dict, Any, List
import math


class AgreementStatistics:
    """Statistical tests for method agreement beyond basic Bland-Altman."""

    def concordance_correlation_coefficient(
        self, method_a: List[float], method_b: List[float]
    ) -> Dict[str, Any]:
        """Lin's CCC - measures both precision and accuracy."""
        n = len(method_a)
        if n < 2:
            return {"ccc": 0.0, "status": "insufficient_data"}

        mean_a = sum(method_a) / n
        mean_b = sum(method_b) / n
        var_a = sum((x - mean_a) ** 2 for x in method_a) / (n - 1)
        var_b = sum((x - mean_b) ** 2 for x in method_b) / (n - 1)
        cov_ab = sum((a - mean_a) * (b - mean_b) for a, b in zip(method_a, method_b)) / (n - 1)

        s2_prime = (var_a + var_b) / 2
        m_prime = (mean_a + mean_b) / 2
        d_prime_sq = (mean_a - mean_b) ** 2

        ccc = (cov_ab + d_prime_sq * (-1)) / (s2_prime + d_prime_sq) if (s2_prime + d_prime_sq) > 0 else 0
        ccc = max(-1, min(1, ccc))

        if ccc >= 0.9:
            quality = "excellent"
        elif ccc >= 0.75:
            quality = "good"
        elif ccc >= 0.5:
            quality = "moderate"
        else:
            quality = "poor"

        return {
            "ccc": round(ccc, 4),
            "quality": quality,
            "mean_difference": round(mean_a - mean_b, 4),
            "precision": round(math.sqrt(var_a * var_b) / s2_prime, 4) if s2_prime > 0 else 0,
            "accuracy_shift": round(d_prime_sq, 4),
        }

    def total_deviation_index(
        self, method_a: List[float], method_b: List[float], tolerance: float = 1.0
    ) -> Dict[str, Any]:
        """TDI - proportion of pairs within tolerance."""
        n = len(method_a)
        diffs = [abs(a - b) for a, b in zip(method_a, method_b)]
        within = sum(1 for d in diffs if d <= tolerance)
        tdi = within / max(n, 1)

        return {
            "tdi": round(tdi, 4),
            "tolerance": tolerance,
            "within_tolerance": within,
            "total_pairs": n,
            "percent_agreement": round(tdi * 100, 1),
        }

    def variance_ratio_test(
        self, method_a: List[float], method_b: List[float]
    ) -> Dict[str, Any]:
        """F-test for equality of variances."""
        n = len(method_a)
        if n < 2:
            return {"status": "insufficient_data"}

        mean_a = sum(method_a) / n
        mean_b = sum(method_b) / n
        var_a = sum((x - mean_a) ** 2 for x in method_a) / (n - 1)
        var_b = sum((x - mean_b) ** 2 for x in method_b) / (n - 1)

        f_stat = var_a / var_b if var_b > 0 else float('inf')

        if f_stat > 1:
            df_n, df_d = n - 1, n - 1
        else:
            f_stat = 1 / f_stat if f_stat > 0 else 0
            df_n, df_d = n - 1, n - 1

        ratio_larger = max(var_a, var_b)
        ratio_smaller = min(var_a, var_b)

        return {
            "f_statistic": round(f_stat, 4),
            "variance_method_a": round(var_a, 4),
            "variance_method_b": round(var_b, 4),
            "variance_ratio": round(ratio_larger / max(ratio_smaller, 1e-10), 4),
            "equal_variances": 0.5 < f_stat < 2.0,
        }

    def paired_t_test(
        self, method_a: List[float], method_b: List[float]
    ) -> Dict[str, Any]:
        """Paired t-test for systematic bias."""
        n = len(method_a)
        if n < 2:
            return {"status": "insufficient_data"}

        diffs = [a - b for a, b in zip(method_a, method_b)]
        mean_d = sum(diffs) / n
        var_d = sum((d - mean_d) ** 2 for d in diffs) / (n - 1)
        se = math.sqrt(var_d / n)

        t_stat = mean_d / se if se > 0 else 0

        return {
            "t_statistic": round(t_stat, 4),
            "mean_difference": round(mean_d, 4),
            "std_difference": round(math.sqrt(var_d), 4),
            "standard_error": round(se, 4),
            "systematic_bias": abs(t_stat) > 2.0,
        }
