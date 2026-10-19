# Preliminary results: selected carbapenemase panel and meropenem

## What was compared

The comparison includes 1,251 public clinical *Klebsiella pneumoniae* isolate records with a submitted meropenem interpretation of resistant (R) or susceptible (S). The 55 intermediate (I) records were set aside. A panel-positive result means at least one selected `blaKPC`, `blaNDM`, or selected OXA-48-like gene annotation was marked `COMPLETE` by NCBI. Each isolate is counted once even if it has multiple panel genes.

![2x2 comparison of the selected NCBI gene panel and meropenem phenotype](figures/panel_vs_meropenem.svg)

## Counts

| Submitted meropenem result | Panel positive | Panel negative | Total |
|---|---:|---:|---:|
| Resistant (R) | 582 | 81 | 663 |
| Susceptible (S) | 46 | 542 | 588 |
| Total | 628 | 623 | 1251 |

## Exploratory estimates

The 95% intervals below are Wilson score intervals calculated from the displayed counts.

| Measure | Estimate (95% interval) |
|---|---|
| Sensitivity | 0.878 (0.851–0.901) |
| Specificity | 0.922 (0.897–0.941) |
| Positive predictive value (PPV) | 0.927 (0.904–0.945) |
| Negative predictive value (NPV) | 0.870 (0.841–0.894) |
| Overall agreement | 0.898 (0.881–0.914) |


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

## Exploratory gene-family comparisons

Each family was compared separately using the same R/S cohort and complete-call rule. Families can overlap within an isolate, so these estimates are not independent and should not be treated as head-to-head statistical comparisons.

| Family | TP | FN | FP | TN | Sensitivity (95% interval) | Specificity (95% interval) |
|---|---:|---:|---:|---:|---|---|
| KPC | 287 | 376 | 4 | 584 | 0.433 (0.396–0.471) | 0.993 (0.983–0.997) |
| NDM | 214 | 449 | 5 | 583 | 0.323 (0.288–0.359) | 0.991 (0.980–0.996) |
| OXA-48-like selected variants | 164 | 499 | 37 | 551 | 0.247 (0.216–0.282) | 0.937 (0.914–0.954) |

These family-specific counts are exploratory descriptions of the same selected public dataset.
