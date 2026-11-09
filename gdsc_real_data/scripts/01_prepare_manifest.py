"""Create the fixed GDSC2 gemcitabine analysis manifest.

The official source files remain unchanged under ``gdsc_real_data/raw/gdsc``.
The script writes only filtered, auditable outputs to ``processed``.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "raw" / "gdsc"
PROCESSED = ROOT / "processed"

TARGET = "Pancreatic Carcinoma"
EXTERNAL = [
    "Colorectal Carcinoma",
    "Gastric Carcinoma",
    "Ovarian Carcinoma",
    "Non-Small Cell Lung Carcinoma",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    response_path = SOURCE / "GDSC2_fitted_dose_response_27Oct23.xlsx"
    model_path = SOURCE / "model_list_20261013.csv"

    response = pd.read_excel(
        response_path,
        sheet_name="Sheet1",
        usecols=[
            "DATASET",
            "CELL_LINE_NAME",
            "SANGER_MODEL_ID",
            "CANCER_TYPE",
            "DRUG_ID",
            "DRUG_NAME",
            "LN_IC50",
            "AUC",
            "Z_SCORE",
        ],
    )
    keep_groups = [TARGET, *EXTERNAL]
    selected = response.loc[
        response["DRUG_NAME"].str.casefold().eq("gemcitabine")
        & response["CANCER_TYPE"].isin(keep_groups)
    ].copy()
    if selected["SANGER_MODEL_ID"].duplicated().any():
        duplicates = selected.loc[
            selected["SANGER_MODEL_ID"].duplicated(False), "SANGER_MODEL_ID"
        ].tolist()
        raise RuntimeError(f"Duplicate GDSC2 model rows: {duplicates}")

    models = pd.read_csv(
        model_path,
        usecols=[
            "model_id",
            "model_name",
            "tissue",
            "cancer_type",
            "cancer_type_detail",
            "model_type",
            "species",
        ],
    )
    manifest = selected.merge(
        models,
        left_on="SANGER_MODEL_ID",
        right_on="model_id",
        how="left",
        validate="one_to_one",
    )
    if manifest["model_id"].isna().any():
        missing = manifest.loc[manifest["model_id"].isna(), "SANGER_MODEL_ID"]
        raise RuntimeError(f"Missing model annotations: {missing.tolist()}")
    if not manifest["CANCER_TYPE"].eq(manifest["cancer_type"]).all():
        raise RuntimeError("GDSC2 and model-list cancer labels disagree")
    if manifest["LN_IC50"].isna().any():
        raise RuntimeError("Missing LN_IC50 in selected records")

    manifest["role"] = manifest["CANCER_TYPE"].map(
        lambda value: "target" if value == TARGET else "external"
    )
    manifest["rnaseq_url"] = manifest["model_id"].map(
        lambda model_id: (
            "https://cog.sanger.ac.uk/cmp/download/"
            f"{model_id}_rnaseq.csv.zip"
        )
    )
    manifest = manifest.rename(
        columns={
            "CELL_LINE_NAME": "gdsc_cell_line_name",
            "CANCER_TYPE": "cancer_group",
        }
    )
    columns = [
        "model_id",
        "model_name",
        "gdsc_cell_line_name",
        "cancer_group",
        "cancer_type_detail",
        "role",
        "DRUG_ID",
        "DRUG_NAME",
        "LN_IC50",
        "AUC",
        "Z_SCORE",
        "rnaseq_url",
    ]
    manifest = manifest[columns].sort_values(
        ["role", "cancer_group", "model_id"]
    )
    counts = manifest.groupby(["role", "cancer_group"]).size().rename("n")

    expected = {
        TARGET: 30,
        "Colorectal Carcinoma": 47,
        "Gastric Carcinoma": 28,
        "Ovarian Carcinoma": 40,
        "Non-Small Cell Lung Carcinoma": 88,
    }
    observed = manifest.groupby("cancer_group").size().to_dict()
    if observed != expected:
        raise RuntimeError(f"Unexpected group counts: {observed}")

    manifest.to_csv(PROCESSED / "model_manifest.csv", index=False)
    counts.reset_index().to_csv(PROCESSED / "group_counts.csv", index=False)
    provenance = pd.DataFrame(
        {
            "source_file": [response_path.name, model_path.name],
            "sha256": [sha256(response_path), sha256(model_path)],
        }
    )
    provenance.to_csv(PROCESSED / "source_checksums.csv", index=False)
    print(counts.to_string())
    print(f"Total models: {len(manifest)}")


if __name__ == "__main__":
    main()
