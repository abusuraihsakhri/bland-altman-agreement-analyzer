#!/usr/bin/env python3
"""
Tests for Bland-Altman Agreement Analyzer.
"""
import math
import os
import sys
import tempfile

import pytest

from bland_altman_core import (
    bland_altman,
    limits_of_agreement,
    mean_bias,
    intraclass_correlation,
    summary,
    process_csv,
    _norm_cdf,
    _norm_ppf,
)


# ---------------------------------------------------------------------------
# Core Bland-Altman analysis
# ---------------------------------------------------------------------------

class TestBlandAltmanBasic:
    def test_identical_methods(self):
        """Identical methods should have zero bias and narrow limits."""
        m1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        m2 = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = bland_altman(m1, m2)
        assert abs(result["mean_difference"]) < 1e-10
        assert abs(result["sd_differences"]) < 1e-10
        assert result["pct_within"] == 100.0

    def test_constant_offset(self):
        """Constant offset should give bias equal to offset."""
        m1 = [10.0, 20.0, 30.0, 40.0, 50.0]
        m2 = [9.0, 19.0, 29.0, 39.0, 49.0]
        result = bland_altman(m1, m2)
        assert abs(result["mean_difference"] - 1.0) < 1e-10
        assert abs(result["sd_differences"]) < 1e-10

    def test_limits_symmetric(self):
        """Limits should be symmetric around the mean difference."""
        m1 = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
        m2 = [9.5, 11.8, 13.5, 15.9, 17.2, 19.8]
        result = bland_altman(m1, m2)
        lo, hi = result["limits_of_agreement"]
        mean_diff = result["mean_difference"]
        assert abs((hi - mean_diff) - (mean_diff - lo)) < 1e-10

    def test_n_subjects(self):
        """N should match input length."""
        m1 = [1.0, 2.0, 3.0]
        m2 = [1.1, 2.1, 3.1]
        result = bland_altman(m1, m2)
        assert result["n"] == 3

    def test_mismatched_lengths_raises(self):
        """Mismatched lengths should raise ValueError."""
        with pytest.raises(ValueError):
            bland_altman([1, 2], [1])

    def test_too_few_raises(self):
        """Less than 2 pairs should raise ValueError."""
        with pytest.raises(ValueError):
            bland_altman([1.0], [1.1])


# ---------------------------------------------------------------------------
# Mean difference (bias)
# ---------------------------------------------------------------------------

class TestMeanBias:
    def test_bias_calculation(self):
        """Bias should be mean of differences."""
        m1 = [10.0, 12.0, 14.0, 16.0, 18.0]
        m2 = [9.0, 11.0, 13.0, 15.0, 17.0]
        result = bland_altman(m1, m2)
        # All differences are 1.0
        assert abs(result["mean_difference"] - 1.0) < 1e-10

    def test_bias_convenience_function(self):
        """Test the convenience function."""
        m1 = [10.0, 12.0, 14.0]
        m2 = [9.0, 11.0, 13.0]
        assert abs(mean_bias(m1, m2) - 1.0) < 1e-10


# ---------------------------------------------------------------------------
# Standard deviation of differences
# ---------------------------------------------------------------------------

class TestSDDifferences:
    def test_sd_calculation(self):
        """SD should be computed correctly."""
        m1 = [10.0, 12.0, 14.0, 16.0, 18.0]
        m2 = [9.0, 11.0, 13.0, 15.0, 17.0]
        result = bland_altman(m1, m2)
        # All differences are 1.0, so SD = 0
        assert abs(result["sd_differences"]) < 1e-10

    def test_sd_with_variation(self):
        """SD should capture variation in differences."""
        m1 = [10.0, 12.0, 14.0, 16.0, 18.0]
        m2 = [9.0, 11.5, 13.0, 15.5, 17.0]
        result = bland_altman(m1, m2)
        # Differences: 1.0, 0.5, 1.0, 0.5, 1.0
        # Mean = 0.8, SD = sqrt(((0.2)^2 + (-0.3)^2 + (0.2)^2 + (-0.3)^2 + (0.2)^2)/4)
        expected_sd = math.sqrt((0.04 + 0.09 + 0.04 + 0.09 + 0.04) / 4)
        assert abs(result["sd_differences"] - expected_sd) < 1e-6


# ---------------------------------------------------------------------------
# Limits of agreement
# ---------------------------------------------------------------------------

class TestLimitsOfAgreement:
    def test_loa_calculation(self):
        """LoA should be d̄ ± 1.96×SD."""
        m1 = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
        m2 = [9.5, 11.8, 13.5, 15.9, 17.2, 19.8]
        result = bland_altman(m1, m2)
        lo, hi = result["limits_of_agreement"]
        mean_diff = result["mean_difference"]
        sd = result["sd_differences"]
        assert abs(lo - (mean_diff - 1.96 * sd)) < 1e-6
        assert abs(hi - (mean_diff + 1.96 * sd)) < 1e-6

    def test_loa_convenience_function(self):
        """Test the convenience function."""
        m1 = [10.0, 12.0, 14.0]
        m2 = [9.0, 11.5, 13.0]
        lo, hi = limits_of_agreement(m1, m2)
        assert lo < hi


# ---------------------------------------------------------------------------
# Confidence intervals
# ---------------------------------------------------------------------------

class TestConfidenceIntervals:
    def test_ci_bias_contains_mean(self):
        """CI for bias should contain the mean difference."""
        m1 = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
        m2 = [9.5, 11.8, 13.5, 15.9, 17.2, 19.8]
        result = bland_altman(m1, m2)
        ci_lo, ci_hi = result["ci_bias"]
        assert ci_lo <= result["mean_difference"] <= ci_hi

    def test_ci_limits_ordered(self):
        """CI limits should be properly ordered."""
        m1 = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
        m2 = [9.5, 11.8, 13.5, 15.9, 17.2, 19.8]
        result = bland_altman(m1, m2)
        assert result["ci_lower_limit"][0] <= result["ci_lower_limit"][1]
        assert result["ci_upper_limit"][0] <= result["ci_upper_limit"][1]


# ---------------------------------------------------------------------------
# Percentage within limits
# ---------------------------------------------------------------------------

class TestPctWithin:
    def test_all_within(self):
        """Identical methods: 100% within limits."""
        m1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        m2 = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = bland_altman(m1, m2)
        assert result["pct_within"] == 100.0

    def test_pct_range(self):
        """Percentage should be between 0 and 100."""
        m1 = [10.0, 12.0, 14.0, 16.0, 18.0]
        m2 = [9.0, 11.0, 13.0, 15.0, 17.0]
        result = bland_altman(m1, m2)
        assert 0.0 <= result["pct_within"] <= 100.0


# ---------------------------------------------------------------------------
# Proportional bias
# ---------------------------------------------------------------------------

class TestProportionalBias:
    def test_no_proportional_bias(self):
        """Constant offset should have no proportional bias."""
        m1 = [10.0, 20.0, 30.0, 40.0, 50.0]
        m2 = [9.0, 19.0, 29.0, 39.0, 49.0]
        result = bland_altman(m1, m2)
        assert not result["proportional_bias"]["proportional_bias"]

    def test_proportional_bias_structure(self):
        """Proportional bias result should have correct keys."""
        m1 = [10.0, 12.0, 14.0, 16.0, 18.0]
        m2 = [9.0, 11.0, 13.0, 15.0, 17.0]
        result = bland_altman(m1, m2)
        pb = result["proportional_bias"]
        assert "slope" in pb
        assert "intercept" in pb
        assert "r_squared" in pb
        assert "p_value" in pb
        assert "proportional_bias" in pb


# ---------------------------------------------------------------------------
# Intraclass Correlation Coefficient
# ---------------------------------------------------------------------------

class TestICC:
    def test_perfect_agreement(self):
        """Perfect agreement should give ICC = 1."""
        m1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        m2 = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = intraclass_correlation(m1, m2)
        assert abs(result["icc"] - 1.0) < 1e-6

    def test_icc_range(self):
        """ICC should be between -1 and 1."""
        m1 = [10.0, 12.0, 14.0, 16.0, 18.0]
        m2 = [9.5, 11.8, 13.5, 15.9, 17.2]
        result = intraclass_correlation(m1, m2)
        assert -1.0 <= result["icc"] <= 1.0

    def test_icc_structure(self):
        """ICC result should have correct keys."""
        m1 = [10.0, 12.0, 14.0]
        m2 = [9.0, 11.0, 13.0]
        result = intraclass_correlation(m1, m2)
        assert "icc" in result
        assert "ci_lower" in result
        assert "ci_upper" in result
        assert "f_statistic" in result

    def test_icc_mismatched_raises(self):
        """Mismatched lengths should raise ValueError."""
        with pytest.raises(ValueError):
            intraclass_correlation([1, 2], [1])


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

class TestSummary:
    def test_summary_output(self):
        """Summary should return formatted string."""
        m1 = [10.0, 12.0, 14.0, 16.0, 18.0]
        m2 = [9.5, 11.8, 13.5, 15.9, 17.2]
        s = summary(m1, m2)
        assert isinstance(s, str)
        assert "Bland-Altman" in s
        assert "ICC" in s


# ---------------------------------------------------------------------------
# CSV processing
# ---------------------------------------------------------------------------

class TestCSVProcessing:
    def test_batch_process(self, tmp_path):
        """Batch CSV processing should produce output."""
        csv_in = tmp_path / "input.csv"
        csv_out = tmp_path / "output.csv"
        csv_in.write_text(
            "subject,method1,method2\n1,14.5,14.2\n2,8.0,7.8\n3,28.4,27.9\n4,12.3,12.1\n5,19.7,19.5\n",
            encoding="utf-8"
        )
        result = process_csv(str(csv_in), str(csv_out))
        assert csv_out.exists()
        assert result["n"] == 5


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------

class TestStatHelpers:
    def test_norm_cdf(self):
        """Standard normal CDF at known values."""
        assert abs(_norm_cdf(0) - 0.5) < 1e-10
        assert abs(_norm_cdf(1.96) - 0.975) < 0.001

    def test_norm_ppf(self):
        """Inverse normal CDF at known values."""
        assert abs(_norm_ppf(0.5)) < 1e-10
        assert abs(_norm_ppf(0.975) - 1.96) < 0.01


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_large_dataset(self):
        """Performance test with larger dataset."""
        import random
        random.seed(42)
        n = 200
        m1 = [random.uniform(5, 50) for _ in range(n)]
        m2 = [x + random.gauss(0, 0.5) for x in m1]
        result = bland_altman(m1, m2)
        assert result["n"] == n
        assert result["pct_within"] > 90.0  # Should be ~95%

    def test_small_dataset(self):
        """Minimum dataset (2 pairs)."""
        m1 = [10.0, 20.0]
        m2 = [9.0, 19.0]
        result = bland_altman(m1, m2)
        assert result["n"] == 2
        assert abs(result["mean_difference"] - 1.0) < 1e-10
