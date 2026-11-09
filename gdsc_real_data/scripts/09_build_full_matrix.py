"""Build the 2,000-feature matrix used by the full GDSC experiment."""

from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw" / "rnaseq"
PROCESSED = ROOT / "processed"
MAX_FEATURES = 2000


def read_expression(model_id: str) -> tuple[pd.Series, pd.DataFrame]:
    path = RAW / f"{model_id}_rnaseq.csv.zip"
    with zipfile.ZipFile(path) as archive:
        with archive.open(f"{model_id}_rnaseq.csv") as stream:
            frame = pd.read_csv(
                stream,
                usecols=["ensembl_gene_id", "gene_symbol", "rsem_tpm"],
            )
    frame = frame.dropna(subset=["ensembl_gene_id", "rsem_tpm"])
    if frame["ensembl_gene_id"].duplicated().any():
        frame = frame.groupby("ensembl_gene_id", as_index=False).agg(
            gene_symbol=("gene_symbol", "first"),
            rsem_tpm=("rsem_tpm", "sum"),
        )
    mapping = frame[["ensembl_gene_id", "gene_symbol"]].drop_duplicates(
        "ensembl_gene_id"
    )
    values = frame.set_index("ensembl_gene_id")["rsem_tpm"].astype(float)
    values.name = model_id
    return values, mapping


def main() -> None:
    manifest = pd.read_csv(PROCESSED / "analysis_manifest.csv")
    expression: list[pd.Series] = []
    gene_mapping: pd.DataFrame | None = None
    for number, model_id in enumerate(manifest["model_id"], start=1):
        values, mapping = read_expression(model_id)
        expression.append(values)
        if gene_mapping is None:
            gene_mapping = mapping
        print(f"[{number:03d}/{len(manifest):03d}] read {model_id}", flush=True)

    tpm = pd.concat(expression, axis=1, join="inner").T
    tpm.index.name = "model_id"
    if tpm.isna().any().any():
        raise RuntimeError("Missing values remain after common-gene join")
    external_ids = manifest.loc[manifest["role"].eq("external"), "model_id"]
    external = tpm.loc[external_ids]
    expression_fraction = (external > 1).mean(axis=0)
    eligible = expression_fraction >= 0.20
    external_log = np.log2(external.loc[:, eligible] + 1)
    variance = external_log.var(axis=0, ddof=1).sort_values(ascending=False)
    selected = variance.head(MAX_FEATURES).index
    if len(selected) != MAX_FEATURES:
        raise RuntimeError(f"Only {len(selected)} eligible genes")

    matrix = np.log2(tpm.loc[:, selected] + 1)
    output = manifest.drop(columns=["rnaseq_url"]).merge(
        matrix.reset_index(), on="model_id", how="left", validate="one_to_one"
    )
    output.to_csv(
        PROCESSED / "gdsc_gemcitabine_top2000.csv.gz", index=False
    )

    feature_log = pd.DataFrame(
        {
            "ensembl_gene_id": selected,
            "external_variance_log2_tpm": variance.loc[selected].values,
            "external_fraction_tpm_gt_1": expression_fraction.loc[selected].values,
            "rank": np.arange(1, MAX_FEATURES + 1),
        }
    )
    if gene_mapping is None:
        raise RuntimeError("Gene mapping was not constructed")
    feature_log = feature_log.merge(
        gene_mapping, on="ensembl_gene_id", how="left", validate="one_to_one"
    )
    feature_log.to_csv(
        PROCESSED / "selected_features_top2000_annotated.csv", index=False
    )
    audit = pd.DataFrame(
        {
            "item": [
                "models",
                "external_models",
                "target_models",
                "common_ensembl_genes",
                "genes_passing_expression_filter",
                "selected_features",
            ],
            "value": [
                len(manifest),
                manifest["role"].eq("external").sum(),
                manifest["role"].eq("target").sum(),
                tpm.shape[1],
                int(eligible.sum()),
                len(selected),
            ],
        }
    )
    audit.to_csv(PROCESSED / "full_matrix_audit.csv", index=False)
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
