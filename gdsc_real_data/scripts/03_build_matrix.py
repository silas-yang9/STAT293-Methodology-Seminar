"""Build the external-selected GDSC gemcitabine expression matrix."""

from __future__ import annotations

import csv
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw" / "rnaseq"
PROCESSED = ROOT / "processed"
N_FEATURES = 1000


def read_tpm(model_id: str) -> pd.Series:
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
            gene_symbol=("gene_symbol", "first"), rsem_tpm=("rsem_tpm", "sum")
        )
    values = frame.set_index("ensembl_gene_id")["rsem_tpm"].astype(float)
    values.name = model_id
    return values


def main() -> None:
    manifest = pd.read_csv(PROCESSED / "model_manifest.csv")
    status = pd.read_csv(PROCESSED / "download_status.csv")
    available_ids = set(status.loc[status["status"].ne("failed"), "model_id"])
    excluded = manifest.loc[~manifest["model_id"].isin(available_ids)].merge(
        status[["model_id", "status", "message"]], on="model_id", how="left"
    )
    excluded.to_csv(PROCESSED / "excluded_missing_rnaseq.csv", index=False)
    manifest = manifest.loc[manifest["model_id"].isin(available_ids)].copy()
    if manifest["role"].eq("target").sum() != 30:
        raise RuntimeError("RNA-seq exclusion removed a target model")
    manifest.to_csv(PROCESSED / "analysis_manifest.csv", index=False)
    series: list[pd.Series] = []
    for number, model_id in enumerate(manifest["model_id"], start=1):
        series.append(read_tpm(model_id))
        print(f"[{number:03d}/{len(manifest):03d}] read {model_id}", flush=True)
    tpm = pd.concat(series, axis=1, join="inner").T
    tpm.index.name = "model_id"
    if tpm.isna().any().any():
        raise RuntimeError("Missing values remain after inner gene join")

    external_ids = manifest.loc[manifest["role"].eq("external"), "model_id"]
    external = tpm.loc[external_ids]
    expressed = (external > 1).mean(axis=0) >= 0.20
    log_external = np.log2(external.loc[:, expressed] + 1)
    variances = log_external.var(axis=0, ddof=1).sort_values(ascending=False)
    selected_genes = variances.head(N_FEATURES).index
    if len(selected_genes) != N_FEATURES:
        raise RuntimeError(f"Only {len(selected_genes)} eligible genes")

    matrix = np.log2(tpm.loc[:, selected_genes] + 1)
    output = manifest.drop(columns=["rnaseq_url"]).merge(
        matrix.reset_index(), on="model_id", how="left", validate="one_to_one"
    )
    output.to_csv(PROCESSED / "gdsc_gemcitabine_top1000.csv.gz", index=False)

    feature_log = pd.DataFrame(
        {
            "ensembl_gene_id": selected_genes,
            "external_variance_log2_tpm": variances.loc[selected_genes].values,
            "external_fraction_tpm_gt_1": (external.loc[:, selected_genes] > 1).mean(axis=0).values,
            "rank": np.arange(1, N_FEATURES + 1),
        }
    )
    feature_log.to_csv(PROCESSED / "selected_features.csv", index=False)
    audit = pd.DataFrame(
        {
            "item": [
                "models",
                "external_models",
                "target_models",
                "models_excluded_missing_rnaseq",
                "common_ensembl_genes",
                "genes_passing_expression_filter",
                "selected_features",
            ],
            "value": [
                len(manifest),
                manifest["role"].eq("external").sum(),
                manifest["role"].eq("target").sum(),
                len(excluded),
                tpm.shape[1],
                int(expressed.sum()),
                len(selected_genes),
            ],
        }
    )
    audit.to_csv(PROCESSED / "matrix_audit.csv", index=False)
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
