"""Add simple gene-family presence flags for the preselected panel."""

from pathlib import Path
import argparse
import re

import pandas as pd


DEFAULT_INPUT = Path("work/data/processed/isolates_analysis_ready.csv")
DEFAULT_OUTPUT = Path("work/data/processed/isolates_with_gene_flags.csv")
GENE_PATTERNS = {
    "blaKPC_present": re.compile(r"^blaKPC(?:-|$)", re.IGNORECASE),
    "blaNDM_present": re.compile(r"^blaNDM(?:-|$)", re.IGNORECASE),
    # Specific NCBI annotations in the selected OXA-48-like family.
    "blaOXA48_like_present": re.compile(
        r"^blaOXA-(?:48|163|181|204|232)(?:-|$)", re.IGNORECASE
    ),
}


def gene_names(genotype_text: str) -> set[str]:
    """Return NCBI AMR gene labels, ignoring the call status after '='."""
    if pd.isna(genotype_text):
        return set()
    return {item.split("=", 1)[0].strip() for item in str(genotype_text).split(",") if item}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Flag KPC, NDM, and selected OXA-48-like gene calls."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    data = pd.read_csv(args.input, encoding="utf-8-sig")
    if "#AMR genotypes" not in data.columns:
        raise ValueError("Input is missing the '#AMR genotypes' column.")

    parsed = data["#AMR genotypes"].map(gene_names)
    for column, pattern in GENE_PATTERNS.items():
        data[column] = parsed.map(lambda names: int(any(pattern.match(n) for n in names)))
    flag_columns = list(GENE_PATTERNS)
    data["selected_carbapenemase_panel_present"] = data[flag_columns].any(axis=1).astype(int)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.output, index=False, encoding="utf-8-sig")

    print(f"Created: {args.output.resolve()}")
    print("Isolates with each selected gene-family call:")
    for column in flag_columns + ["selected_carbapenemase_panel_present"]:
        print(f"  {column}: {int(data[column].sum())}")
    print("\nCounts by meropenem result (descriptive only):")
    print(data.groupby("meropenem_phenotype")[flag_columns + ["selected_carbapenemase_panel_present"]].sum().to_string())
    print("\nFlags mean an NCBI annotation was present; they do not prove gene expression or clinical resistance.")


if __name__ == "__main__":
    main()
