# Cohort feasibility and source record

## Search and export

- Source: [NCBI Pathogen Detection Isolates Browser](https://www.ncbi.nlm.nih.gov/pathogens/isolates/)
- Retrieval date: 17 October 2026
- Search: `taxgroup_name:"Klebsiella pneumoniae" AND epi_type:clinical AND (AST_phenotypes:meropenem=R OR AST_phenotypes:meropenem=S OR AST_phenotypes:meropenem=I) AND AMR_genotypes:*`
- Export: metadata CSV, 1,306 rows and 18 columns
- Raw file: `work/data/raw/isolates.csv` (kept unchanged)

## Initial checks

Pandas successfully parsed the CSV, including its header that begins with `#AMR genotypes`.

- Organism group: 1,306 *Klebsiella pneumoniae* records
- Isolation type: 1,306 clinical records
- Unique isolate accessions: 1,306; duplicate accessions: 0
- Meropenem AST calls: 663 resistant, 588 susceptible, 55 intermediate
- AMR genotype field and AST phenotype field were populated for all 1,306 rows

The browser also showed 436 matched SNP clusters. Cluster counts are not isolate counts. These are cohort and metadata checks, not genotype–phenotype analysis results. The data are public and may change; retain this retrieval date and query when reporting the analysis.

## Missing metadata observed

Some descriptive fields are incomplete. In this export, missing counts include Assembly (271), Isolation source (193), Location (126), Strain (510), SNP cluster (465), Min-diff (1,194), and Min-same (481). Computed types, Food origin, and Serovar were missing for all rows. Missingness in optional descriptive fields does not remove an isolate from the initial cohort.

## Re-run the basic audit

From the project folder, run:

```powershell
python outputs/scripts/01_audit_metadata.py
```

This script checks row/column counts, duplicate isolate accessions, organism and isolation-type labels, meropenem phenotype counts, and missing fields. It does not calculate concordance or diagnostic performance.
