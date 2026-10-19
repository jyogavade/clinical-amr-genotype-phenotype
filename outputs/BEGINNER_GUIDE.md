# Beginner guide: the AMR gene and meropenem project

This guide explains what the project does, where its files are, and how to repeat the analysis. You do not need Google Cloud or NCBI BigQuery. We use the NCBI Isolates Browser, a downloaded metadata table, Python scripts, and an optional Jupyter notebook.

## 1. The project question

We compare two kinds of information about the same bacterial isolate:

- **Phenotype:** the laboratory antimicrobial susceptibility test (AST) result. For meropenem, the export reports `R` (resistant), `S` (susceptible), or `I` (intermediate).
- **Genotype:** resistance-associated gene annotations reported by NCBI AMRFinderPlus.

The question is: **in this downloaded group of public *Klebsiella pneumoniae* records, how often does a selected carbapenemase-gene panel agree with the submitted meropenem AST result?**

This is a learning project using public records. It does not prove that a gene caused a result, validate a clinical test, or guide anyone's treatment.

## 2. The project folders

The project folder is the main folder containing this README and the `outputs` directory. On your own computer, its exact location depends on where you save or clone the project. For example:

```text
C:\path\to\clinical-amr-genotype-phenotype
```

Important locations:

- `work/data/raw/isolates.csv` — the unchanged NCBI download.
- `work/data/processed/` — working copies and analysis tables created by scripts.
- `outputs/scripts/` — Python instructions that perform each step.
- `outputs/figures/` — figures.
- `outputs/RESULTS_SUMMARY.md` — current results and limitations.

Keep the original CSV unchanged. The analysis scripts read it and create separate files.

## 3. How to run a Python script on your computer

Use **PowerShell**, not the File Explorer search box. You can open PowerShell with **Windows key + R**, type `powershell`, then press Enter.

In a new PowerShell window, run these two setup lines, one at a time:

```powershell
Set-Location 'C:\path\to\clinical-amr-genotype-phenotype'
$PY = 'C:\path\to\Python3\python.exe'
```

Replace each example path with the actual folder and Python 3 program path on your computer. The first line changes PowerShell to the project folder. The second stores the Python program path in a short variable called `$PY`. The `$PY` setting lasts only in that PowerShell window; repeat the setup lines in a new one. For a fresh installation, see the [main README](../README.md) for creating a Python 3 environment.

Run a script like this:

```powershell
& $PY outputs/scripts/01_audit_metadata.py
```

The command prompt should start with the project folder name. A line like `PS C:\Users\YourName>` means you still need to run `Set-Location` first. Do not type the prompt itself.

## 4. What we did, in order

### Step 1 — Choose a cohort in NCBI

We used the [NCBI Pathogen Detection Isolates Browser](https://www.ncbi.nlm.nih.gov/pathogens/isolates/) with this search:

```text
taxgroup_name:"Klebsiella pneumoniae" AND epi_type:clinical AND (AST_phenotypes:meropenem=R OR AST_phenotypes:meropenem=S OR AST_phenotypes:meropenem=I) AND AMR_genotypes:*
```

The export was downloaded on 17 October 2026. It contained 1,306 records: 663 `R`, 588 `S`, and 55 `I`. The browser showed 436 SNP clusters. A cluster is a group of related records; it is not the same as one isolate.

### Step 2 — Check the CSV

Run:

```powershell
& $PY outputs/scripts/01_audit_metadata.py
```

The script checks the number of rows, repeated isolate IDs, organism and source labels, meropenem categories, and missing optional information. The initial check found 1,306 unique isolate accessions and no duplicate accession values.

### Step 3 — Make a tidy working copy

Run:

```powershell
& $PY outputs/scripts/02_prepare_analysis_table.py
```

This creates `work/data/processed/isolates_analysis_ready.csv`. It adds a separate `meropenem_phenotype` column containing `R`, `S`, or `I`. The original AST and genotype text stays in the table.

### Step 4 — Make exploratory gene flags

Run:

```powershell
& $PY outputs/scripts/03_flag_carbapenemase_genes.py
```

This counts annotations with names from the selected `blaKPC`, `blaNDM`, and selected OXA-48-like variants. This early flag includes calls even if NCBI labels them partial, so we inspect status before the main comparison.

### Step 5 — Review NCBI call status

Run:

```powershell
& $PY outputs/scripts/04_review_gene_call_status.py
```

This creates `selected_gene_call_details.csv`, where each selected gene annotation has its own row and its NCBI status, such as `COMPLETE` or `PARTIAL`.

### Step 6 — Compare the combined panel with R/S

Run:

```powershell
& $PY outputs/scripts/05_compare_panel_with_meropenem.py
```

For this main comparison, **panel positive** means at least one selected gene call is marked `COMPLETE`. The comparison uses the 663 resistant and 588 susceptible records. It sets the 55 intermediate records aside rather than combining them with susceptible records.

| Submitted result | Complete gene panel present | No complete panel call |
|---|---:|---:|
| Resistant (R) | 582 | 81 |
| Susceptible (S) | 46 | 542 |

The percentages are summaries of this database export. They are not the performance of a clinical test.

### Step 7 — Review mismatches

Run:

```powershell
& $PY outputs/scripts/06_review_discordant_isolates.py
```

This creates a list of the 81 resistant/panel-negative and 46 susceptible/panel-positive records, plus a frequency table of their other NCBI annotations. Some listed genes may be common background genes; a frequent annotation is not automatically the explanation for a mismatch. Review the original call and AST details before proposing a biological explanation.

### Step 8 — Compare each family separately

Run:

```powershell
& $PY outputs/scripts/08_compare_gene_families.py
```

This writes `gene_family_comparison_summary.csv`. In these records, KPC, NDM, and the selected OXA-48-like variants each have lower sensitivity when considered alone than the combined panel. One isolate may carry multiple families, so family comparisons overlap.

### Step 9 — Make the report and simple figure

Run:

```powershell
& $PY outputs/scripts/07_make_results_figure.py
```

This creates the SVG figure and refreshes `outputs/RESULTS_SUMMARY.md`. Run step 8 first if you want the family-specific table included in that report.

## 5. Make a polished figure in Jupyter (optional)

Open `outputs/plot_results_in_jupyter.ipynb` in Jupyter and run all cells from top to bottom. The notebook reads the summary created in step 6 and saves a 300-dpi PNG and SVG in `outputs/figures/`. It needs `pandas` and `matplotlib` in the selected Python 3 kernel.

The black PowerShell window that launched Jupyter is the server console. It stays open while the browser notebook is using it. To return to its prompt, save the notebook, click the console, press **Ctrl+C**, and confirm shutdown with `y`. If you want to type other commands while Jupyter runs, open a second PowerShell window and use the setup lines in section 3.

## 6. What the results mean—and do not mean

- The combined panel matched the submitted category in 1,124 of the 1,251 R/S records (about 89.8%).
- The 46 susceptible records with panel calls and the 81 resistant records without complete panel calls are worth reviewing; they do not automatically indicate incorrect AST or incorrect gene results.
- Other carbapenemases, porin changes, gene expression, AST methods, breakpoints, and data quality may matter. This project does not establish a cause for each mismatch.
- The 436 SNP clusters suggest related isolates. The report's Wilson confidence intervals treat records as independent and do not adjust for relatedness, so use them cautiously.
- NCBI records are public submissions and are not a random sample of all clinical infections. Database contents can change; keep the retrieval date and query with any report.

## 7. Small glossary

- **Isolate:** a sample of bacteria with a database identifier.
- **AST:** antimicrobial susceptibility testing, a laboratory test of how bacteria respond to an antibiotic.
- **AMR:** antimicrobial resistance.
- **AMRFinderPlus:** NCBI software that identifies known resistance genes and markers in submitted genome data.
- **Cohort:** the records selected for an analysis using stated rules.
- **Concordance/agreement:** how often the two categories match in the selected records.
- **Sensitivity/specificity:** here, mathematical summaries of the table using the gene panel as the positive/negative classification. They are not clinical validation results.

## 8. Current deliverables

- [Main results and limitations](RESULTS_SUMMARY.md)
- [Cohort and source notes](COHORT_FEASIBILITY.md)
- [Main SVG figure](figures/panel_vs_meropenem.svg)
- [Interactive Jupyter figure notebook](plot_results_in_jupyter.ipynb)
