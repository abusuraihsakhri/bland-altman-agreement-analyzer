# Bland-Altman Agreement Analyzer

A pure-Python (stdlib-only) implementation of Bland-Altman analysis for method comparison studies, including limits of agreement, confidence intervals, proportional bias detection, and intraclass correlation coefficient.

## Features

- **Mean difference (bias)**: d̄ = Σ(di)/n where di = method1 − method2
- **Standard deviation of differences**: SD = √(Σ(di−d̄)²/(n−1))
- **95% Limits of agreement**: d̄ ± 1.96×SD
- **95% Confidence intervals**: For bias and limits of agreement
- **Percentage within limits**: Proportion of observations within LoA
- **Proportional bias detection**: Regression of differences on means
- **Intraclass Correlation Coefficient (ICC)**: Two-way random model
- **CSV batch processing**

## Formulas

### Mean Difference (Bias)
d̄ = Σ(di) / n, where di = method1_i − method2_i

### Limits of Agreement
Lower = d̄ − 1.96 × SD
Upper = d̄ + 1.96 × SD

### 95% CI for Bias
CI = d̄ ± 1.96 × (SD / √n)

### 95% CI for Limits
SE_limit = SD × √(3/n)
CI = limit ± 1.96 × SE_limit

### ICC (Two-way Random, Single Measures)
ICC = (MSB − MSE) / (MSB + MSE)

## Usage

### Command Line

```bash
# Full Bland-Altman analysis
python cli.py analyze --method1 10 12 14 16 18 --method2 9.5 11.8 13.5 15.9 17.2

# ICC only
python cli.py icc --method1 10 12 14 16 18 --method2 9.5 11.8 13.5 15.9 17.2

# Batch CSV processing
python cli.py batch --input sample.csv --output results.csv
```

### Python API

```python
from bland_altman_core import bland_altman, intraclass_correlation

# Full analysis
result = bland_altman(
    method1=[14.5, 8.0, 28.4, 12.3, 19.7],
    method2=[14.2, 7.8, 27.9, 12.1, 19.5]
)
print(f"Bias: {result['mean_difference']:.4f}")
print(f"Limits: {result['limits_of_agreement']}")
print(f"Within limits: {result['pct_within']:.1f}%")

# ICC
icc = intraclass_correlation(method1, method2)
print(f"ICC: {icc['icc']:.4f}")
```

## CSV Format

| subject | method1 | method2 |
|---------|---------|---------|
| 1       | 14.5    | 14.2    |
| 2       | 8.0     | 7.8     |

## Testing

```bash
python -m pytest test_bland_altman_core.py -v
```

## License

MIT
