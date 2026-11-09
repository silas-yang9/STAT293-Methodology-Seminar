from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "conservative_alpha_experiment"
PROCESSED = ROOT / "processed"


def add_check(rows: list[dict[str, str]], name: str, condition: bool, detail: str) -> None:
    rows.append({"check": name, "status": "PASS" if condition else "FAIL", "detail": detail})


def main() -> None:
    predictions = pd.read_csv(RESULTS / "predictions.csv")
    diagnostics = pd.read_csv(RESULTS / "diagnostics.csv")
    alpha_errors = pd.read_csv(RESULTS / "alpha_errors.csv")
    summary = pd.read_csv(RESULTS / "method_summary.csv")
    manifest = pd.read_csv(PROCESSED / "analysis_manifest.csv")
    repeats = int(predictions["repeat_id"].nunique())
    expected_fits = repeats * 5
    methods = [
        "Pretrained global-alpha",
        "Pretrained target-min",
        "Pretrained target-1se",
        "Overall",
        "Individual",
        "Target mean",
    ]
    checks: list[dict[str, str]] = []

    add_check(checks, "target sample count", (manifest["role"] == "target").sum() == 30, "expected 30")
    add_check(checks, "external sample count", (manifest["role"] == "external").sum() == 199, "expected 199")
    add_check(
        checks,
        "prediction row count",
        len(predictions) == repeats * 30 * len(methods),
        f"observed {len(predictions)}, expected {repeats * 30 * len(methods)}",
    )
    add_check(
        checks,
        "unique prediction keys",
        not predictions.duplicated(["repeat_id", "model_id", "method"]).any(),
        "one row per repeat, target model, and method",
    )
    coverage = predictions.groupby(["repeat_id", "method"])["model_id"].nunique()
    add_check(checks, "complete target coverage", bool((coverage == 30).all()), "all repeat-method cells contain 30 targets")
    outcome_counts = predictions.groupby(["repeat_id", "model_id"])["observed"].nunique()
    add_check(checks, "outcome agreement", bool((outcome_counts == 1).all()), "all methods share each observed outcome")
    add_check(
        checks,
        "diagnostic row count",
        len(diagnostics) == expected_fits,
        f"observed {len(diagnostics)}, expected {expected_fits}",
    )
    add_check(
        checks,
        "alpha-error row count",
        len(alpha_errors) == expected_fits * 11,
        f"observed {len(alpha_errors)}, expected {expected_fits * 11}",
    )
    alpha_grid = {round(x / 10, 10) for x in range(11)}
    alpha_columns = ["alpha_global", "alpha_target_min", "alpha_target_1se"]
    alpha_ok = all(diagnostics[column].round(10).isin(alpha_grid).all() for column in alpha_columns)
    add_check(checks, "alpha grid", alpha_ok, "all alpha rules select from 0.0, 0.1, ..., 1.0")
    add_check(checks, "fixed feature count", bool((diagnostics["features"] == 2000).all()), "all fits use 2,000 external-selected genes")
    add_check(checks, "fixed lambda rule", bool((diagnostics["lambda_choice"] == "lambda.1se").all()), "both stages use lambda.1se")

    expected_1se: list[dict[str, float | int]] = []
    for (repeat_id, fold), frame in alpha_errors.groupby(["repeat_id", "fold"], sort=False):
        frame = frame.sort_values("alpha")
        best = frame.loc[frame["target_inner_mse"].idxmin()]
        threshold = float(best["target_inner_mse"] + best["target_inner_se"])
        selected = float(frame.loc[frame["target_inner_mse"] <= threshold + 1e-12, "alpha"].min())
        expected_1se.append({"repeat_id": repeat_id, "fold": fold, "expected_alpha_target_1se": selected})
    expected_frame = pd.DataFrame(expected_1se)
    rule_check = diagnostics.merge(expected_frame, on=["repeat_id", "fold"], how="left")
    add_check(
        checks,
        "one-SE alpha rule",
        bool(np.allclose(rule_check["alpha_target_1se"], rule_check["expected_alpha_target_1se"], atol=1e-12)),
        "saved alpha equals the smallest alpha within one SE of the pancreatic inner-CV minimum",
    )

    recomputed = (
        predictions.assign(sq_error=lambda x: (x["observed"] - x["predicted"]) ** 2)
        .groupby(["repeat_id", "method"], as_index=False)["sq_error"]
        .mean()
        .groupby("method", as_index=False)["sq_error"]
        .mean()
        .rename(columns={"sq_error": "recomputed_mse"})
    )
    reconciled = summary.merge(recomputed, on="method", how="left")
    add_check(
        checks,
        "MSE reconciliation",
        bool(((reconciled["mean_mse"] - reconciled["recomputed_mse"]).abs() < 1e-10).all()),
        "method_summary mean MSE matches independent recomputation",
    )

    checks_frame = pd.DataFrame(checks)
    checks_frame.to_csv(RESULTS / "audit_checks.csv", index=False)
    print(checks_frame.to_string(index=False))
    if checks_frame["status"].ne("PASS").any():
        raise SystemExit("Conservative alpha audit failed")


if __name__ == "__main__":
    main()

