"""Expand selected gene annotations and summarize their NCBI call status."""

from pathlib import Path
import argparse
import re

import pandas as pd


DEFAULT_INPUT = Path("work/data/processed/isolates_analysis_ready.csv")
DEFAULT_OUTPUT = Path("work/data/processed/selected_gene_call_details.csv")
FAMILIES = {
    "KPC": re.compile(r"^blaKPC(?:-|$)", re.IGNORECASE),
    "NDM": re.compile(r"^blaNDM(?:-|$)", re.IGNORECASE),
    "OXA-48-like selected variants": re.compile(
        r"^blaOXA-(?:48|163|181|204|232)(?:-|$)", re.IGNORECASE
    ),
}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="List selected carbapenemase annotations and their NCBI status."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    data = pd.read_csv(args.input, encoding="utf-8-sig")
    required = {"Isolate", "meropenem_phenotype", "#AMR genotypes"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Required columns missing from input: {sorted(missing)}")

    records = []
    for row in data[["Isolate", "meropenem_phenotype", "#AMR genotypes"]].itertuples(index=False, name=None):
        isolate, phenotype, genotype_text = row
        if pd.isna(genotype_text):
            continue
        for call in str(genotype_text).split(","):
            gene, separator, status = call.partition("=")
            if not separator:
                continue
            for family, pattern in FAMILIES.items():
                if pattern.match(gene.strip()):
                    records.append({
                        "Isolate": isolate,
                        "meropenem_phenotype": phenotype,
                        "gene_family": family,
                        "gene_call": gene.strip(),
                        "NCBI_call_status": status.strip(),
                    })
                    break

    details = pd.DataFrame(records, columns=[
        "Isolate", "meropenem_phenotype", "gene_family", "gene_call", "NCBI_call_status"
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    details.to_csv(args.output, index=False, encoding="utf-8-sig")

    print(f"Created detail file: {args.output.resolve()}")
    print(f"Selected gene annotations listed: {len(details)}")
    print("Annotation counts by family, NCBI call status, and phenotype:")
    if details.empty:
        print("No selected annotations found.")
    else:
        summary = details.groupby(
            ["gene_family", "NCBI_call_status", "meropenem_phenotype"]
        ).size().unstack("meropenem_phenotype", fill_value=0)
        print(summary.to_string())
    print("These are annotation-call counts; one isolate can have multiple calls.")


if __name__ == "__main__":
    main()
