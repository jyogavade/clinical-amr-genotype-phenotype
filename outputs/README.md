# Meropenem resistance and carbapenemase annotations in *Klebsiella pneumoniae*

## Project at a glance

This beginner portfolio project compares NCBI AMRFinderPlus gene annotations with submitted meropenem antimicrobial susceptibility testing (AST) results for public clinical *Klebsiella pneumoniae* isolates. The question is whether a selected carbapenemase gene panel agrees with the submitted resistant or susceptible category in this downloaded dataset.

This is an exploratory analysis of public records. It is not a clinical diagnostic evaluation and must not be used to guide treatment.

## Current results

The NCBI export downloaded on 17 October 2026 contained 1,306 records: 663 resistant (R), 588 susceptible (S), and 55 intermediate (I). All were labeled clinical in the exported isolation-type field. The R/S comparison includes 1,251 records; the 55 I records are kept separate and excluded from that comparison.

For the main comparison, a panel-positive isolate has at least one **COMPLETE** NCBI call from the selected `blaKPC`, `blaNDM`, or selected OXA-48-like variants. The 2×2 table is:

| Submitted meropenem result | Panel positive | Panel negative | Total |
|---|---:|---:|---:|
| Resistant (R) | 582 | 81 | 663 |
| Susceptible (S) | 46 | 542 | 588 |

The combined panel had 87.8% sensitivity and 92.2% specificity in these records. Full estimates, denominators, Wilson intervals, family-level comparisons, and limitations are in [RESULTS_SUMMARY.md](RESULTS_SUMMARY.md). A figure is available as [SVG](figures/panel_vs_meropenem.svg); the Jupyter notebook can create a high-resolution PNG and another SVG.

The browser displayed 436 SNP clusters among 1,306 records. Some records may therefore be related. The reported Wilson intervals do not adjust for cluster relatedness and may be too narrow.

## Data source and cohort

- Source: [NCBI Pathogen Detection Isolates Browser](https://www.ncbi.nlm.nih.gov/pathogens/isolates/)
- Retrieval date: 17 October 2026
- Search query:

  ```text
  taxgroup_name:"Klebsiella pneumoniae" AND epi_type:clinical AND (AST_phenotypes:meropenem=R OR AST_phenotypes:meropenem=S OR AST_phenotypes:meropenem=I) AND AMR_genotypes:*
  ```

- Export: metadata CSV with 1,306 rows and 18 columns
- Local raw file: `work/data/raw/isolates.csv` (kept unchanged; not included as a project deliverable)

NCBI BigQuery and Google Cloud are not required for this workflow. No raw sequencing reads are downloaded or analyzed.

## Reproduce the analysis

Use Python 3.9 or newer and install `pandas`. The Jupyter figure notebook also needs `matplotlib`. In a fresh project checkout, download the metadata CSV using the query above and save it as `work/data/raw/isolates.csv`. Keep the source export unchanged. Run commands from the project folder using a Python 3 interpreter. The root [README](../README.md) shows how to create a virtual environment on Windows.

Run these scripts in order:

```powershell
python outputs/scripts/01_audit_metadata.py
python outputs/scripts/02_prepare_analysis_table.py
python outputs/scripts/03_flag_carbapenemase_genes.py
python outputs/scripts/04_review_gene_call_status.py
python outputs/scripts/05_compare_panel_with_meropenem.py
python outputs/scripts/06_review_discordant_isolates.py
python outputs/scripts/08_compare_gene_families.py
python outputs/scripts/07_make_results_figure.py
```

Each script reads a prior file and writes a separate output under `work/data/processed/`; the raw input remains unchanged. Script 07 reads the family summary if script 08 has already produced it.

### Optional Jupyter figure

Open [plot_results_in_jupyter.ipynb](plot_results_in_jupyter.ipynb) in Jupyter and run its cells from top to bottom. The notebook reads `panel_vs_meropenem_summary.csv` and saves `panel_vs_meropenem_jupyter.png` (300 dpi) and `.svg` in `outputs/figures/`. It requires `pandas` and `matplotlib` in the Jupyter kernel.

## What the scripts do

1. `01_audit_metadata.py` checks row counts, identifiers, labels, meropenem results, and missing fields.
2. `02_prepare_analysis_table.py` extracts the `R`, `S`, or `I` label into a dedicated column.
3. `03_flag_carbapenemase_genes.py` makes exploratory family flags using any matching annotation; it does not require a complete call.
4. `04_review_gene_call_status.py` lists selected gene annotations and their NCBI status.
5. `05_compare_panel_with_meropenem.py` builds the main 2×2 table. It includes all R and S records, excludes I, and requires at least one selected **COMPLETE** call for panel-positive.
6. `06_review_discordant_isolates.py` lists the 81 resistant/panel-negative and 46 susceptible/panel-positive records and summarizes other gene annotations for review.
7. `08_compare_gene_families.py` calculates a separate exploratory comparison for each family. A single isolate can carry more than one family.
8. `07_make_results_figure.py` produces the SVG figure and results summary.

## Interpretation and limitations

- NCBI AST results are submitter-reported. Test methods and breakpoint standards may vary.
- Public database isolates are selected, not a random sample of all infections.
- The selected gene panel is incomplete. A panel-negative result does not mean that no resistance mechanism exists.
- A detected gene annotation does not prove expression or explain a specific AST result.
- Gene-family comparisons overlap and are not independent head-to-head tests.
- The 436 SNP clusters suggest potential relatedness. Confidence intervals are ordinary Wilson intervals and are not cluster-adjusted.
- Dataset contents can change; the query and retrieval date must accompany any reported results.

## Project files

- [BEGINNER_GUIDE.md](BEGINNER_GUIDE.md) — walkthrough for a new learner
- [COHORT_FEASIBILITY.md](COHORT_FEASIBILITY.md) — query, retrieval, and initial data audit
- [RESULTS_SUMMARY.md](RESULTS_SUMMARY.md) — counts, estimates, figure, and limitations
- `scripts/` — analysis scripts
- `figures/` — figure output
- `plot_results_in_jupyter.ipynb` — interactive figure notebook

## CV description

Built a reproducible Python workflow to compare NCBI AMRFinderPlus carbapenemase annotations with submitter-reported meropenem AST results in 1,306 public clinical *K. pneumoniae* isolate records. Audited metadata, defined a complete-call gene panel, calculated exploratory agreement statistics, reviewed discordant records, and documented data and clustering limitations.

Use this description only after you can explain the analysis and its limitations.
