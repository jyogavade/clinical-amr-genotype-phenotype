"""Basic, repeatable audit of the NCBI isolate metadata CSV."""

from pathlib import Path
import argparse

import pandas as pd


DEFAULT_INPUT = Path("work/data/raw/isolates.csv")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check the NCBI export's rows, identifiers, labels, and missing values."
    )
    parser.add_argument(
        "csv_path",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"CSV file to audit (default: {DEFAULT_INPUT})",
    )
    args = parser.parse_args()

    data = pd.read_csv(args.csv_path, encoding="utf-8-sig")

    required_columns = {
        "#AMR genotypes",
        "AST phenotypes",
        "Isolate",
        "Isolation type",
        "Organism group",
    }
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Required columns missing from CSV: {sorted(missing_columns)}")

    meropenem_calls = data["AST phenotypes"].str.extract(
        r"(?:^|,)meropenem=([^,]+)(?:,|$)"
    )[0].rename("Meropenem phenotype")

    print(f"File: {args.csv_path.resolve()}")
    print(f"Rows: {len(data)}")
    print(f"Columns: {len(data.columns)}")
    print(f"Unique isolate accessions: {data['Isolate'].nunique(dropna=False)}")
    print(f"Duplicate isolate accessions: {data['Isolate'].duplicated().sum()}")
    print("Organism group counts:")
    print(data["Organism group"].value_counts(dropna=False).to_string())
    print("Isolation type counts:")
    print(data["Isolation type"].value_counts(dropna=False).to_string())
    print("Meropenem phenotype counts:")
    print(meropenem_calls.value_counts(dropna=False).sort_index().to_string())
    print("Columns with missing values:")
    missing = data.isna().sum()
    missing = missing[missing > 0].sort_values(ascending=False)
    print(missing.to_string() if not missing.empty else "None")


if __name__ == "__main__":
    main()
