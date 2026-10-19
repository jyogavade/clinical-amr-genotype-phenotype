"""Identify isolates that do not match the selected gene panel and summarize calls."""

from pathlib import Path
import argparse

import pandas as pd


DEFAULT_COHORT = Path("work/data/processed/isolates_analysis_ready.csv")
DEFAULT_CALLS = Path("work/data/processed/selected_gene_call_details.csv")
DEFAULT_DISCORDANT = Path("work/data/processed/discordant_isolates.csv")
DEFAULT_SUMMARY = Path("work/data/processed/discordant_gene_summary.csv")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="List R isolates without, and S isolates with, a complete selected-panel call."
    )
    parser.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    parser.add_argument("--calls", type=Path, default=DEFAULT_CALLS)
    parser.add_argument("--discordant-output", type=Path, default=DEFAULT_DISCORDANT)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY)
    args = parser.parse_args()

    cohort = pd.read_csv(args.cohort, encoding="utf-8-sig")
    calls = pd.read_csv(args.calls, encoding="utf-8-sig")
    complete_ids = set(calls.loc[calls["NCBI_call_status"].eq("COMPLETE"), "Isolate"])
    cohort = cohort[cohort["meropenem_phenotype"].isin(["R", "S"])].copy()
    cohort["selected_panel_complete"] = cohort["Isolate"].isin(complete_ids)
    cohort["discordance_group"] = ""
    cohort.loc[(cohort["meropenem_phenotype"] == "R") & ~cohort["selected_panel_complete"], "discordance_group"] = "Resistant, panel negative"
    cohort.loc[(cohort["meropenem_phenotype"] == "S") & cohort["selected_panel_complete"], "discordance_group"] = "Susceptible, panel positive"
    discordant = cohort[cohort["discordance_group"].ne("")].copy()

    keep = ["Isolate", "meropenem_phenotype", "discordance_group", "#AMR genotypes", "AST phenotypes"]
    args.discordant_output.parent.mkdir(parents=True, exist_ok=True)
    discordant[keep].to_csv(args.discordant_output, index=False, encoding="utf-8-sig")

    records = []
    for row in discordant[["Isolate", "discordance_group", "#AMR genotypes"]].itertuples(index=False, name=None):
        isolate, group, genotype_text = row
        if pd.isna(genotype_text):
            continue
        for item in str(genotype_text).split(","):
            gene, separator, status = item.partition("=")
            if separator:
                records.append({"Isolate": isolate, "discordance_group": group,
                                "gene_call": gene.strip(), "NCBI_call_status": status.strip()})
    detail = pd.DataFrame(records)
    if detail.empty:
        summary = pd.DataFrame(columns=["discordance_group", "gene_call", "isolates_with_call"])
    else:
        summary = (detail.drop_duplicates(["Isolate", "discordance_group", "gene_call"])
                   .groupby(["discordance_group", "gene_call"])["Isolate"]
                   .nunique().rename("isolates_with_call").reset_index()
                   .sort_values(["discordance_group", "isolates_with_call", "gene_call"],
                                ascending=[True, False, True]))
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_output, index=False, encoding="utf-8-sig")

    print(f"Discordant isolates saved: {args.discordant_output.resolve()}")
    print(f"Gene-call frequency table saved: {args.summary_output.resolve()}")
    print("Discordant isolate counts:")
    print(discordant["discordance_group"].value_counts().to_string())
    for group in ["Resistant, panel negative", "Susceptible, panel positive"]:
        print(f"\nMost frequent gene annotations in {group}:")
        part = summary[summary["discordance_group"].eq(group)].head(12)
        print(part[["gene_call", "isolates_with_call"]].to_string(index=False) if not part.empty else "No annotations")
    print("\nGene annotations are candidates to review, not proof of a resistance mechanism.")


if __name__ == "__main__":
    main()
