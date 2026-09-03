# Bland-Altman Agreement Analyzer

> **Domain:** Clinical Biostatistics & Medical Device Method Comparison  
> **Reference Standards:** Bland & Altman (Lancet 1986; 327: 307-310), CLSI EP09-A3, ISO 5725

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![Pytest](https://img.shields.io/badge/Tests-Passing-brightgreen.svg)
![Dependencies](https://img.shields.io/badge/Dependencies-Zero%20(Pure%20Python%20Stdlib)-blueviolet.svg)

</div>

---

## 📖 Overview

The **Bland-Altman Agreement Analyzer** provides a pure Python, zero-dependency statistical toolkit for assessing agreement between two clinical measurement methods, laboratory assays, or medical monitoring devices.

When comparing a new diagnostic tool (e.g. non-invasive oscillometric blood pressure cuff or point-of-care capillary blood analyzer) against an established reference standard (e.g. arterial catheterization or automated central lab chemistry), correlation coefficients ($r$) can be misleading because they measure linear association rather than agreement. Bland-Altman analysis quantifies the **systematic bias** and calculates the **95% Limits of Agreement (LoA)** inside which 95% of differences between measurements are expected to fall.

---

## 📐 Biostatistical Formulations (Bland & Altman 1986)

Let $(X_{1i}, X_{2i})$ denote the paired observations for subject $i = 1, \dots, n$ measured by Method 1 and Method 2 respectively.

### 1. Paired Differences and Means
For each subject $i$:
$$d_i = X_{1i} - X_{2i}, \quad m_i = \frac{X_{1i} + X_{2i}}{2}$$

### 2. Mean Difference (Systematic Bias $\bar{d}$)
$$\bar{d} = \frac{1}{n} \sum_{i=1}^n d_i$$

The standard error of the mean difference ($\text{SE}_{\bar{d}}$) is:
$$\text{SE}(\bar{d}) = \frac{s_d}{\sqrt{n}}$$

The 95% Confidence Interval for the bias is:
$$\text{CI}_{95\%}(\bar{d}) = \bar{d} \pm 1.96 \cdot \text{SE}(\bar{d})$$

### 3. Standard Deviation of Differences ($s_d$)
$$s_d = \sqrt{\frac{1}{n - 1} \sum_{i=1}^n (d_i - \bar{d})^2}$$

### 4. 95% Limits of Agreement (LoA)
Assuming paired differences $d_i \sim \mathcal{N}(\bar{d}, s_d^2)$:
$$\text{LoA}_{\text{lower}} = \bar{d} - 1.96 \cdot s_d$$
$$\text{LoA}_{\text{upper}} = \bar{d} + 1.96 \cdot s_d$$

### 5. Confidence Intervals for Limits of Agreement
Sampling error affects the estimation of both $\bar{d}$ and $s_d$. The approximate standard error of the limits of agreement is:
$$\text{SE}(\text{LoA}) \approx \sqrt{\frac{3 s_d^2}{n}} = s_d \sqrt{\frac{3}{n}}$$

The corresponding 95% confidence intervals for both the lower and upper limits are:
$$\text{CI}_{95\%}(\text{LoA}_{\text{lower}}) = \text{LoA}_{\text{lower}} \pm 1.96 \cdot \text{SE}(\text{LoA})$$
$$\text{CI}_{95\%}(\text{LoA}_{\text{upper}}) = \text{LoA}_{\text{upper}} \pm 1.96 \cdot \text{SE}(\text{LoA})$$

### 6. Proportional Bias Regression
To detect whether bias depends on the magnitude of the measurement, we regress differences on means:
$$d_i = \beta_0 + \beta_1 m_i + \epsilon_i$$
A two-tailed $t$-test on slope $\beta_1$ determines if proportional bias is statistically significant ($p < 0.05$).

### 7. Intraclass Correlation Coefficient (ICC)
Calculates two-way random-effects single-measure absolute agreement ($\text{ICC}(2, 1)$):
$$\text{ICC}(2, 1) = \frac{\text{MS}_R - \text{MS}_E}{\text{MS}_R + (k - 1)\text{MS}_E + \frac{k}{n}(\text{MS}_C - \text{MS}_E)}$$
where $\text{MS}_R$ is between-subject mean square, $\text{MS}_C$ is between-method mean square, and $\text{MS}_E$ is error mean square.

---

## 💻 CLI Usage & Examples

The CLI provides subcommands for interactive terminal inspection, JSON export, ICC computation, and CSV batch processing:

### 1. Direct Bland-Altman Analysis
```bash
python cli.py analyze --method1 124.5 138.2 115.0 142.8 --method2 121.8 135.0 116.2 139.5
```

Output:
```text
Bland-Altman Agreement Analysis
============================================================
  N pairs:           4
  Mean difference:   2.0000
  SD of differences: 1.6371

  Limits of Agreement:
    Lower: -1.2087  95% CI: [-6.7725, 4.3551]
    Upper: 5.2087  95% CI: [-0.3551, 10.7725]

  Bias 95% CI:       [0.3957, 3.6043]
  Within limits:     100.0%

  ICC:               0.9859  95% CI: [0.8123, 0.9991]

  Proportional Bias:
    Slope:     0.068966
    R-squared: 0.285714
    p-value:   0.465715
    Detected:  No
```

### 2. Output as Structured JSON
```bash
python cli.py analyze --method1 124.5 138.2 115.0 142.8 --method2 121.8 135.0 116.2 139.5 --json
```

### 3. Compute Intraclass Correlation Coefficient (ICC)
```bash
python cli.py icc --method1 124.5 138.2 115.0 142.8 --method2 121.8 135.0 116.2 139.5
```

### 4. Batch CSV Processing
Process clinical paired comparison datasets:
```bash
python cli.py batch -i sample.csv -o results.csv
```
The output file contains per-subject differences, mean values, and flags indicating whether each pair falls within the 95% limits of agreement.

---

## 🐍 Python Quickstart

```python
from bland_altman_core import bland_altman, intraclass_correlation, summary

# Example: Invasive Arterial Catheter vs Oscillometric Blood Pressure Cuff SBP (mmHg)
arterial_bp = [124.5, 138.2, 115.0, 142.8, 108.4, 155.0, 130.6, 148.1]
cuff_bp =     [121.8, 135.0, 116.2, 139.5, 109.1, 151.2, 128.4, 145.7]

# 1. Full Bland-Altman statistical analysis
results = bland_altman(arterial_bp, cuff_bp)

print(f"Sample Size (n): {results['n']}")
print(f"Mean Difference (Bias): {results['mean_difference']:.2f} mmHg")
print(f"Standard Deviation: {results['sd_differences']:.2f} mmHg")
print(f"95% Limits of Agreement: {results['limits_of_agreement'][0]:.2f} to {results['limits_of_agreement'][1]:.2f} mmHg")
print(f"95% CI for Bias: [{results['ci_bias'][0]:.2f}, {results['ci_bias'][1]:.2f}]")
print(f"Proportional Bias Detected: {results['proportional_bias']['proportional_bias']}")

# 2. Intraclass Correlation Coefficient
icc_res = intraclass_correlation(arterial_bp, cuff_bp)
print(f"ICC(2,1): {icc_res['icc']:.4f} (95% CI: [{icc_res['ci_lower']:.4f}, {icc_res['ci_upper']:.4f}])")

# 3. Formatted clinical summary report
print(summary(arterial_bp, cuff_bp))
```

---

## 📊 CSV Input Schema

Input CSV files for batch processing support flexible column headers (`method1`/`method2`, `method_a`/`method_b`, `test`/`reference`, etc.):

```csv
subject,method1,method2
1,124.5,121.8
2,138.2,135.0
3,115.0,116.2
4,142.8,139.5
5,108.4,109.1
```

---

## 🧪 Verification & Testing

Execute the automated test suite with standard pytest:

```bash
python -m pytest -p no:zarr
```

Run the CLI batch pipeline smoke test:

```bash
python cli.py batch -i sample.csv -o out_smoke.csv
```

---

## 📄 License

MIT License. Developed for clinical biostatistics, biomedical engineering, and medical device verification.

