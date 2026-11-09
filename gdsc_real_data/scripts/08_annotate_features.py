"""Attach gene symbols to selected features and stability summaries."""

from __future__ import annotations

import zipfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw" / "rnaseq"
PROCESSED = ROOT / "processed"
REPEATED = ROOT / "results" / "repeated_cv"


def main() -> None:
    zip_path = sorted(RAW.glob("*_rnaseq.csv.zip"))[0]
    with zipfile.ZipFile(zip_path) as archive:
        with archive.open(archive.namelist()[0]) as stream:
            mapping = pd.read_csv(
                stream, usecols=["ensembl_gene_id", "gene_symbol"]
            ).drop_duplicates("ensembl_gene_id")

    features = pd.read_csv(PROCESSED / "selected_features.csv")
    features.merge(mapping, on="ensembl_gene_id", how="left", validate="one_to_one").to_csv(
        PROCESSED / "selected_features_annotated.csv", index=False
    )
    stability = pd.read_csv(REPEATED / "selection_stability.csv")
    stability.merge(mapping, on="ensembl_gene_id", how="left", validate="many_to_one").to_csv(
        REPEATED / "selection_stability_annotated.csv", index=False
    )


if __name__ == "__main__":
    main()
