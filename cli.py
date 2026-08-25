#!/usr/bin/env python3
"""
CLI for Bland-Altman Agreement Analyzer.

Usage:
  python cli.py analyze --method1 10 12 14 16 18 --method2 9.5 11.8 13.5 15.9 17.2
  python cli.py icc --method1 10 12 14 16 18 --method2 9.5 11.8 13.5 15.9 17.2
  python cli.py batch --input sample.csv --output results.csv
"""

import argparse
import json
import sys

from bland_altman_core import bland_altman, intraclass_correlation, summary, process_csv


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="bland-altman",
        description="Bland-Altman Agreement Analyzer",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # --- analyze ---
    p_analyze = sub.add_parser("analyze", help="Run Bland-Altman analysis")
    p_analyze.add_argument("--method1", type=float, nargs="+", required=True, help="Method 1 measurements")
    p_analyze.add_argument("--method2", type=float, nargs="+", required=True, help="Method 2 measurements")
    p_analyze.add_argument("--json", action="store_true", help="Output as JSON")

    # --- icc ---
    p_icc = sub.add_parser("icc", help="Compute Intraclass Correlation Coefficient")
    p_icc.add_argument("--method1", type=float, nargs="+", required=True)
    p_icc.add_argument("--method2", type=float, nargs="+", required=True)

    # --- batch ---
    p_batch = sub.add_parser("batch", help="Process CSV file")
    p_batch.add_argument("-i", "--input", required=True, help="Input CSV path")
    p_batch.add_argument("-o", "--output", default="results.csv", help="Output CSV path")

    args = parser.parse_args(argv)

    if args.command == "analyze":
        if args.json:
            result = bland_altman(args.method1, args.method2)
            # Remove non-serializable lists
            result.pop("differences", None)
            result.pop("means", None)
            print(json.dumps(result, indent=2))
        else:
            print(summary(args.method1, args.method2))
        return 0

    elif args.command == "icc":
        result = intraclass_correlation(args.method1, args.method2)
        print("Intraclass Correlation Coefficient")
        print("=" * 40)
        print(f"  ICC:        {result['icc']:.4f}")
        print(f"  95% CI:     [{result['ci_lower']:.4f}, {result['ci_upper']:.4f}]")
        print(f"  F-statistic: {result['f_statistic']:.4f}")
        return 0

    elif args.command == "batch":
        result = process_csv(args.input, args.output)
        print(f"Processed {result['n']} pairs -> {args.output}")
        return 0

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
