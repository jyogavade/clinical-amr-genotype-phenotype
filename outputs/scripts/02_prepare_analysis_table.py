"""Create a tidy one-row-per-isolate table from the unchanged NCBI export."""

from pathlib import Path
import argparse

import pandas as pd


DEFAULT_INPUT = Path("work/data/raw/isolates.csv")
DEFAULT_OUTPUT = Path("work/data/processed/isolates_analysis_ready.csv")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Copy NCBI metadata to a working table and extract meropenem AST."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    data = pd.read_csv(args.input, encoding="utf-8-sig")
    required = {"Isolate", "AST phenotypes", "#AMR genotypes"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Required columns missing from input: {sorted(missing)}")

    data.insert(
        0,
        "meropenem_phenotype",
        data["AST phenotypes"].str.extract(
            r"(?:^|,)meropenem=([^,]+)(?:,|$)", expand=False
        ),
    )
    data["meropenem_phenotype"] = data["meropenem_phenotype"].str.strip()

    if data["Isolate"].isna().any() or data["Isolate"].duplicated().any():
        raise ValueError("Expected exactly one row per isolate with a unique accession.")
    if data["meropenem_phenotype"].isna().any():
        raise ValueError("At least one row has no meropenem phenotype; inspect the source.")
    unexpected = sorted(set(data["meropenem_phenotype"]) - {"R", "S", "I"})
    if unexpected:
        raise ValueError(f"Unexpected meropenem labels found: {unexpected}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.output, index=False, encoding="utf-8-sig")

    print(f"Created: {args.output.resolve()}")
    print(f"Rows: {len(data)} (one row per isolate)")
    print("Meropenem phenotype counts:")
    print(data["meropenem_phenotype"].value_counts().sort_index().to_string())
    print("The original genotype and phenotype text columns are retained unchanged.")


if __name__ == "__main__":
    main()
