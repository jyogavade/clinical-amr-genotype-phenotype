"""Compare each selected gene family with submitted meropenem R/S results."""

from pathlib import Path
import argparse
import math

import pandas as pd


DEFAULT_COHORT = Path("work/data/processed/isolates_analysis_ready.csv")
DEFAULT_CALLS = Path("work/data/processed/selected_gene_call_details.csv")
DEFAULT_OUTPUT = Path("work/data/processed/gene_family_comparison_summary.csv")


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    den = 1 + z * z / n
    center = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return center - half, center + half


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate one 2x2 table per selected gene family.")
    parser.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    parser.add_argument("--calls", type=Path, default=DEFAULT_CALLS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    cohort = pd.read_csv(args.cohort, encoding="utf-8-sig")
    calls = pd.read_csv(args.calls, encoding="utf-8-sig")
    selected = cohort[cohort["meropenem_phenotype"].isin(["R", "S"])].copy()
    families = list(calls["gene_family"].dropna().unique())
    rows = []
    for family in families:
        family_complete = set(calls.loc[
            calls["gene_family"].eq(family) & calls["NCBI_call_status"].eq("COMPLETE"),
            "Isolate",
        ])
        positive = selected["Isolate"].isin(family_complete)
        resistant = selected["meropenem_phenotype"].eq("R")
        susceptible = selected["meropenem_phenotype"].eq("S")
        tp = int((resistant & positive).sum())
        fn = int((resistant & ~positive).sum())
        fp = int((susceptible & positive).sum())
        tn = int((susceptible & ~positive).sum())
        metrics = [
            ("Sensitivity", tp, tp + fn),
            ("Specificity", tn, tn + fp),
            ("PPV", tp, tp + fp),
            ("NPV", tn, tn + fn),
            ("Overall agreement", tp + tn, tp + fn + fp + tn),
        ]
        row = {"gene_family": family, "TP": tp, "FN": fn, "FP": fp, "TN": tn,
               "R_total": tp + fn, "S_total": tn + fp}
        for name, k, n in metrics:
            estimate = k / n if n else math.nan
            low, high = wilson(k, n)
            row[name] = estimate
            row[f"{name}_ci95_low"] = low
            row[f"{name}_ci95_high"] = high
        rows.append(row)

    result = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False, encoding="utf-8-sig")
    print(f"Saved: {args.output.resolve()}")
    print("Per-family counts and selected estimates:")
    for row in rows:
        print(f"\n{row['gene_family']}: TP={row['TP']}, FN={row['FN']}, FP={row['FP']}, TN={row['TN']}")
        for name in ["Sensitivity", "Specificity", "PPV", "NPV"]:
            print(f"  {name}: {row[name]:.3f} ({row[name + '_ci95_low']:.3f}–{row[name + '_ci95_high']:.3f})")
    print("\nThese are exploratory within-dataset comparisons, not clinical test performance.")


if __name__ == "__main__":
    main()
