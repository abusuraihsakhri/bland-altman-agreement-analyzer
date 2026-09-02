# Bland Altman Agreement Analyzer

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Bland-Altman Agreement — Extended Statistics
Proportional-bias regression (difference on mean), Lin's CCC, Gwet's AC1,
repeatability coefficient, and bootstrap CIs for bias/LoA.

Zero-dependency. Author: Dr. Abu Suraih Sakhri. License: MIT.

Bland-Altman Statistical Tests: regression on differences, variance ratio, and agreement indices.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`MethodPair`** — dedicated module for method pair evaluation and state verification.
- **`AgreementStatistics`**: Statistical tests for method agreement beyond basic Bland-Altman.
- **`MethodPair`** — dedicated module for method pair evaluation and state verification.
- **`BlandAltmanAnalyzer`**: Bland-Altman analysis for method comparison in clinical assays.

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --method1 <value> --method2 <value> --input <value> --output <value>
```

### Parameter Reference
- `--method1`: Specifies input measurement or parameter value.
- `--method2`: Specifies input measurement or parameter value.
- `--input`: Specifies input measurement or parameter value.
- `--output`: Specifies input measurement or parameter value.
- `---`: Specifies input measurement or parameter value.
- `--json`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `subject` | Parameter / observation metric | Required |
| `method1` | Parameter / observation metric | Required |
| `method2` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t bland-altman-agreement-analyzer .
docker run -p 8000:8000 bland-altman-agreement-analyzer
```
