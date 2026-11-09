from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "nested_tuning_experiment"
PROCESSED = ROOT / "processed"


def add_check(rows: list[dict[str, str]], name: str, condition: bool, detail: str) -> None:
    rows.append({"check": name, "status": "PASS" if condition else "FAIL", "detail": detail})


def main() -> None:
    predictions = pd.read_csv(RESULTS / "predictions.csv")
    diagnostics = pd.read_csv(RESULTS / "diagnostics.csv")
    candidates = pd.read_csv(RESULTS / "candidate_inner_metrics.csv")
    summary = pd.read_csv(RESULTS / "method_summary.csv")
    manifest = pd.read_csv(PROCESSED / "analysis_manifest.csv")
    repeats = int(predictions["repeat_id"].nunique())
    expected_fits = repeats * 5
    methods = ["Individual", "Overall", "Pretrained", "Target mean"]
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
        len(diagnostics) == expected_fits * 3,
        f"observed {len(diagnostics)}, expected {expected_fits * 3}",
    )
    add_check(
        checks,
        "candidate row count",
        len(candidates) == expected_fits * 4,
        f"observed {len(candidates)}, expected {expected_fits * 4}",
    )
    add_check(
        checks,
        "target-specific alpha grid",
        diagnostics.loc[diagnostics["method"] == "Pretrained", "selected_target_alpha"]
        .round(10)
        .isin([round(x / 10, 10) for x in range(11)])
        .all(),
        "all selected alpha values belong to 0.0, 0.1, ..., 1.0",
    )
    add_check(
        checks,
        "candidate choices",
        diagnostics["selected_candidate"].isin(["f1000_1se", "f2000_1se", "f1000_min", "f2000_min"]).all(),
        "all method choices belong to the frozen candidate grid",
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
        raise SystemExit("Nested tuning audit failed")


if __name__ == "__main__":
    main()

