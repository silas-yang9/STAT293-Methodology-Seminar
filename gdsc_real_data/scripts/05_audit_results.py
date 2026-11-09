"""Audit the completed pilot without refitting any model."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "processed"
RESULTS = ROOT / "results"
FOLDS = RESULTS / "folds"
BOOTSTRAP_REPS = 20000
SEED = 2932027


def metrics(frame: pd.DataFrame) -> dict[str, float | int | str]:
    residual = frame["observed"] - frame["predicted"]
    denominator = ((frame["observed"] - frame["observed"].mean()) ** 2).sum()
    return {
        "method": frame["method"].iloc[0],
        "n": len(frame),
        "mse": float((residual**2).mean()),
        "rmse": float(np.sqrt((residual**2).mean())),
        "mae": float(residual.abs().mean()),
        "correlation": float(frame["observed"].corr(frame["predicted"])),
        "r_squared": float(1 - (residual**2).sum() / denominator),
    }


def main() -> None:
    data = pd.read_csv(PROCESSED / "gdsc_gemcitabine_top1000.csv.gz")
    folds = pd.read_csv(PROCESSED / "outer_fold_assignment.csv")
    target = data.loc[data["role"].eq("target"), ["model_id", "LN_IC50"]].merge(
        folds, on="model_id", validate="one_to_one"
    )
    baseline_rows = []
    for fold in sorted(target["outer_fold"].unique()):
        train_mean = target.loc[target["outer_fold"].ne(fold), "LN_IC50"].mean()
        test = target.loc[target["outer_fold"].eq(fold)]
        baseline_rows.append(
            pd.DataFrame(
                {
                    "model_id": test["model_id"],
                    "fold": fold,
                    "observed": test["LN_IC50"],
                    "method": "Target mean",
                    "predicted": train_mean,
                }
            )
        )
    baseline = pd.concat(baseline_rows, ignore_index=True)

    predictions = pd.read_csv(RESULTS / "predictions.csv")
    if predictions.duplicated(["model_id", "method"]).any():
        raise RuntimeError("Duplicate model-method predictions")
    counts = predictions.groupby("method")["model_id"].nunique().to_dict()
    if counts != {"Individual": 30, "Overall": 30, "Pretrained": 30}:
        raise RuntimeError(f"Incomplete predictions: {counts}")
    audited = pd.concat([predictions, baseline], ignore_index=True)
    audited.to_csv(RESULTS / "audited_predictions.csv", index=False)

    overall_metrics = pd.DataFrame(
        [metrics(frame) for _, frame in audited.groupby("method")]
    ).sort_values("mse")
    overall_metrics.to_csv(RESULTS / "audit_metrics.csv", index=False)

    fold_metrics = pd.DataFrame(
        [
            {"fold": fold, **metrics(frame)}
            for (fold, _), frame in audited.groupby(["fold", "method"])
        ]
    ).sort_values(["fold", "mse"])
    fold_metrics.to_csv(RESULTS / "fold_metrics.csv", index=False)

    wide = audited.pivot(index="model_id", columns="method", values=["observed", "predicted"])
    observed = wide["observed"].iloc[:, 0].to_numpy()
    squared_errors = {
        method: (observed - wide["predicted"][method].to_numpy()) ** 2
        for method in wide["predicted"].columns
    }
    rng = np.random.default_rng(SEED)
    indices = rng.integers(0, len(observed), size=(BOOTSTRAP_REPS, len(observed)))
    comparisons = [
        ("Individual", "Pretrained"),
        ("Overall", "Pretrained"),
        ("Target mean", "Pretrained"),
        ("Individual", "Overall"),
    ]
    bootstrap_rows = []
    for comparator, reference in comparisons:
        per_model_difference = squared_errors[comparator] - squared_errors[reference]
        draws = per_model_difference[indices].mean(axis=1)
        bootstrap_rows.append(
            {
                "contrast": f"MSE({comparator}) - MSE({reference})",
                "estimate": per_model_difference.mean(),
                "ci_2.5_percent": np.quantile(draws, 0.025),
                "ci_97.5_percent": np.quantile(draws, 0.975),
                "bootstrap_probability_positive": (draws > 0).mean(),
                "bootstrap_reps": BOOTSTRAP_REPS,
                "seed": SEED,
            }
        )
    pd.DataFrame(bootstrap_rows).to_csv(
        RESULTS / "paired_bootstrap_differences.csv", index=False
    )

    support_frames = [
        pd.read_csv(path) for path in sorted(FOLDS.glob("fold_*_support.csv"))
    ]
    supports = pd.concat(support_frames, ignore_index=True)
    support_counts = (
        supports.dropna(subset=["ensembl_gene_id"])
        .groupby(["fold", "method"])["ensembl_gene_id"]
        .nunique()
        .reindex(
            pd.MultiIndex.from_product(
                [range(1, 6), ["Individual", "Overall", "Pretrained"]],
                names=["fold", "method"],
            ),
            fill_value=0,
        )
        .rename("selected_features")
        .reset_index()
    )
    support_counts.to_csv(RESULTS / "support_counts.csv", index=False)

    outcome_summary = (
        data.groupby(["role", "cancer_group"])["LN_IC50"]
        .agg(["count", "mean", "std", "min", "median", "max"])
        .reset_index()
    )
    outcome_summary.to_csv(RESULTS / "outcome_summary_by_group.csv", index=False)
    print(overall_metrics.to_string(index=False))
    print("\nPaired bootstrap differences:")
    print(pd.DataFrame(bootstrap_rows).to_string(index=False))


if __name__ == "__main__":
    main()
