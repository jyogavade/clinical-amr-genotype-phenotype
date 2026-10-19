"""Create an SVG 2x2 figure and a plain-language results report."""

from pathlib import Path
import argparse

import pandas as pd


DEFAULT_SUMMARY = Path("work/data/processed/panel_vs_meropenem_summary.csv")
DEFAULT_FIGURE = Path("outputs/figures/panel_vs_meropenem.svg")
DEFAULT_REPORT = Path("outputs/RESULTS_SUMMARY.md")


def find_count(table: pd.DataFrame, name: str) -> int:
    row = table.loc[table["measure"].eq(name)]
    if row.empty:
        raise ValueError(f"Missing count row: {name}")
    return int(row.iloc[0]["numerator"])


def main() -> None:
    parser = argparse.ArgumentParser(description="Make the 2x2 figure and a short report.")
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    table = pd.read_csv(args.summary, encoding="utf-8-sig")
    tp = find_count(table, "True positive (R and panel positive)")
    fn = find_count(table, "False negative (R and panel negative)")
    fp = find_count(table, "False positive (S and panel positive)")
    tn = find_count(table, "True negative (S and panel negative)")

    values = [[tp, fn], [fp, tn]]
    colors = [["#dcefe5", "#fff0d6"], ["#fff0d6", "#dcefe5"]]
    svg = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="540" viewBox="0 0 900 540">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#192b27}.head{font-size:22px;font-weight:bold}.label{font-size:17px;font-weight:bold}.num{font-size:36px;font-weight:bold}.tag{font-size:14px}</style>',
        '<text x="450" y="48" text-anchor="middle" class="head">NCBI gene panel vs. meropenem result</text>',
        '<text x="500" y="100" text-anchor="middle" class="label">Panel positive</text>',
        '<text x="730" y="100" text-anchor="middle" class="label">Panel negative</text>',
        '<text x="30" y="220" class="label">Resistant (R)</text>',
        '<text x="30" y="400" class="label">Susceptible (S)</text>',
    ]
    for r in range(2):
        for c in range(2):
            x, y = 390 + c * 230, 135 + r * 180
            svg.append(f'<rect x="{x}" y="{y}" width="210" height="155" rx="4" fill="{colors[r][c]}" stroke="#49615a" stroke-width="2"/>')
            label = "MATCH" if r == c else "MISMATCH"
            svg.append(f'<text x="{x + 105}" y="{y + 78}" text-anchor="middle" class="num">{values[r][c]}</text>')
            svg.append(f'<text x="{x + 105}" y="{y + 118}" text-anchor="middle" class="tag">{label}</text>')
    svg.extend([
        '<text x="450" y="500" text-anchor="middle" font-family="Arial,sans-serif" font-size="13">Positive = at least one selected gene call marked COMPLETE; 55 intermediate isolates excluded.</text>',
        '</svg>',
    ])
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    args.figure.write_text("\n".join(svg), encoding="utf-8")

    report = f"""# Preliminary results: selected carbapenemase panel and meropenem

## What was compared

The comparison includes {tp + fn + fp + tn:,} public clinical *Klebsiella pneumoniae* isolate records with a submitted meropenem interpretation of resistant (R) or susceptible (S). The 55 intermediate (I) records were set aside. A panel-positive result means at least one selected `blaKPC`, `blaNDM`, or selected OXA-48-like gene annotation was marked `COMPLETE` by NCBI. Each isolate is counted once even if it has multiple panel genes.

![2x2 comparison of the selected NCBI gene panel and meropenem phenotype](figures/panel_vs_meropenem.svg)

## Counts

| Submitted meropenem result | Panel positive | Panel negative | Total |
|---|---:|---:|---:|
| Resistant (R) | {tp} | {fn} | {tp + fn} |
| Susceptible (S) | {fp} | {tn} | {fp + tn} |
| Total | {tp + fp} | {fn + tn} | {tp + fn + fp + tn} |

## Exploratory estimates

The 95% intervals below are Wilson score intervals calculated from the displayed counts.

| Measure | Estimate (95% interval) |
|---|---|
"""
    for name in ["Sensitivity", "Specificity", "Positive predictive value (PPV)",
                 "Negative predictive value (NPV)", "Overall agreement"]:
        row = table.loc[table["measure"].eq(name)].iloc[0]
        report += f"| {name} | {float(row['estimate']):.3f} ({float(row['ci95_low']):.3f}–{float(row['ci95_high']):.3f}) |\n"
    report += """

## How to read this

Within this particular exported dataset, 582 resistant isolates had a complete call from the selected panel, while 81 did not. Also, 46 susceptible isolates had a complete panel call. The panel therefore does not perfectly separate the two submitted AST categories in this dataset.

This is an exploratory comparison of public, submitter-provided data. It is not a clinical diagnostic test evaluation, does not show that a gene caused an isolate's phenotype, and should not guide treatment. Isolates in public databases are selected rather than a random sample; AST methods, breakpoints, assemblies, and gene-call quality may vary. The gene panel is limited, and other resistance mechanisms may be present. Intermediate isolates were excluded from this R/S comparison, not discarded from the project.

The browser reported 436 SNP clusters for 1,306 records, so some isolates may be related. The Wilson intervals shown here treat records as independent and do not adjust for clustering; they may be too narrow. Interpret them as descriptive uncertainty intervals, not as population-level inference.

The discordant gene-frequency table is exploratory. Frequently annotated genes such as `emrD` or `fosA` appear in both discordant groups and do not by themselves explain the mismatch. Review individual NCBI call status, AST metadata, and candidate alternative mechanisms before making biological claims.

## Reproduce

Run scripts 04, 05, and 08 first, followed by:

```powershell
python outputs/scripts/07_make_results_figure.py
```

The underlying comparison summary is `work/data/processed/panel_vs_meropenem_summary.csv`.
"""
    family_path = Path("work/data/processed/gene_family_comparison_summary.csv")
    if family_path.exists():
        family_table = pd.read_csv(family_path, encoding="utf-8-sig")
        report += "\n## Exploratory gene-family comparisons\n\n"
        report += "Each family was compared separately using the same R/S cohort and complete-call rule. Families can overlap within an isolate, so these estimates are not independent and should not be treated as head-to-head statistical comparisons.\n\n"
        report += "| Family | TP | FN | FP | TN | Sensitivity (95% interval) | Specificity (95% interval) |\n"
        report += "|---|---:|---:|---:|---:|---|---|\n"
        for _, row in family_table.iterrows():
            report += (f"| {row['gene_family']} | {int(row['TP'])} | {int(row['FN'])} | "
                       f"{int(row['FP'])} | {int(row['TN'])} | "
                       f"{row['Sensitivity']:.3f} ({row['Sensitivity_ci95_low']:.3f}–{row['Sensitivity_ci95_high']:.3f}) | "
                       f"{row['Specificity']:.3f} ({row['Specificity_ci95_low']:.3f}–{row['Specificity_ci95_high']:.3f}) |\n")
        report += "\nThese family-specific counts are exploratory descriptions of the same selected public dataset.\n"
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report, encoding="utf-8")
    print(f"Figure created: {args.figure.resolve()}")
    print(f"Results report created: {args.report.resolve()}")


if __name__ == "__main__":
    main()
