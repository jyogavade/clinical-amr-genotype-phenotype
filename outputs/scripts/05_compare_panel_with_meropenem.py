"""Calculate an exploratory 2x2 comparison using complete gene calls only."""

from pathlib import Path
import argparse
import math

import pandas as pd


DEFAULT_INPUT = Path("work/data/processed/isolates_analysis_ready.csv")
DEFAULT_CALLS = Path("work/data/processed/selected_gene_call_details.csv")
DEFAULT_OUTPUT = Path("work/data/processed/panel_vs_meropenem_summary.csv")


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """95% Wilson score interval for a binomial proportion."""
    if total == 0:
        return (math.nan, math.nan)
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return (center - half, center + half)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare complete selected gene-family annotations with meropenem R/S."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--calls", type=Path, default=DEFAULT_CALLS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    cohort = pd.read_csv(args.input, encoding="utf-8-sig")
    calls = pd.read_csv(args.calls, encoding="utf-8-sig")
    required = {"Isolate", "meropenem_phenotype"}
    missing = required.difference(cohort.columns)
    if missing:
        raise ValueError(f"Required columns missing from cohort: {sorted(missing)}")
    if not {"Isolate", "NCBI_call_status"}.issubset(calls.columns):
        raise ValueError("Gene-call detail file is missing required columns.")

    # Only COMPLETE calls count as a positive gene-panel result in this primary comparison.
    complete_ids = set(calls.loc[calls["NCBI_call_status"].eq("COMPLETE"), "Isolate"])
    phenotypes = cohort[["Isolate", "meropenem_phenotype"]].drop_duplicates()
    if phenotypes["Isolate"].duplicated().any():
        raise ValueError("An isolate has conflicting phenotype entries in the detail file.")
    phenotypes = phenotypes[phenotypes["meropenem_phenotype"].isin(["R", "S"])].copy()
    phenotypes["panel_positive"] = phenotypes["Isolate"].isin(complete_ids)

    tp = int(((phenotypes["meropenem_phenotype"] == "R") & phenotypes["panel_positive"]).sum())
    fn = int(((phenotypes["meropenem_phenotype"] == "R") & ~phenotypes["panel_positive"]).sum())
    fp = int(((phenotypes["meropenem_phenotype"] == "S") & phenotypes["panel_positive"]).sum())
    tn = int(((phenotypes["meropenem_phenotype"] == "S") & ~phenotypes["panel_positive"]).sum())

    measures = [
        ("Sensitivity", tp, tp + fn),
        ("Specificity", tn, tn + fp),
        ("Positive predictive value (PPV)", tp, tp + fp),
        ("Negative predictive value (NPV)", tn, tn + fn),
        ("Overall agreement", tp + tn, tp + fn + fp + tn),
    ]
    rows = [
        {"measure": "True positive (R and panel positive)", "numerator": tp, "denominator": tp + fn, "estimate": ""},
        {"measure": "False negative (R and panel negative)", "numerator": fn, "denominator": tp + fn, "estimate": ""},
        {"measure": "False positive (S and panel positive)", "numerator": fp, "denominator": tn + fp, "estimate": ""},
        {"measure": "True negative (S and panel negative)", "numerator": tn, "denominator": tn + fp, "estimate": ""},
    ]
    print(f"R/S isolates analyzed: {len(phenotypes)} (intermediate excluded: {int(cohort.loc[cohort['meropenem_phenotype'].eq('I'), 'Isolate'].nunique())})")
    print("2x2 table (rows = meropenem phenotype, columns = gene panel):")
    print(f"                 Panel positive   Panel negative")
    print(f"Resistant (R)          {tp:>4}              {fn:>4}")
    print(f"Susceptible (S)        {fp:>4}              {tn:>4}")
    print("\nExploratory estimates with 95% Wilson confidence intervals:")
    for name, successes, total in measures:
        estimate = successes / total if total else math.nan
        low, high = wilson(successes, total)
        rows.append({"measure": name, "numerator": successes, "denominator": total,
                     "estimate": estimate, "ci95_low": low, "ci95_high": high})
        print(f"  {name}: {estimate:.3f} ({low:.3f}–{high:.3f}; {successes}/{total})")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False, encoding="utf-8-sig")
    print(f"\nSummary saved: {args.output.resolve()}")
    print("Interpret only as agreement within this selected public dataset, not as clinical diagnostic performance.")


if __name__ == "__main__":
    main()
