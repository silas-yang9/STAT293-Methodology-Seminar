"""Independent integrity checks for the full experiment outputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "full_experiment"
PROCESSED = ROOT / "processed"


def check(name: str, condition: bool, detail: str) -> dict[str, str]:
    return {"check": name, "status": "PASS" if condition else "FAIL", "detail": detail}


def main() -> None:
    settings = pd.read_csv(RESULTS / "scenario_settings.csv")
    data = pd.read_csv(PROCESSED / "gdsc_gemcitabine_top2000.csv.gz")
    checks: list[dict[str, str]] = []
    checks.append(check("target sample count", (data["role"] == "target").sum() == 30, "expected 30"))
    checks.append(check("external sample count", (data["role"] == "external").sum() == 199, "expected 199"))
    alpha_grid = np.round(np.arange(0, 1.01, 0.1), 10)

    for setting in settings.itertuples(index=False):
        scenario_dir = RESULTS / setting.scenario
        predictions = pd.read_csv(scenario_dir / "predictions.csv")
        diagnostics = pd.read_csv(scenario_dir / "diagnostics.csv")
        summary = pd.read_csv(scenario_dir / "method_summary.csv")
        expected_predictions = setting.repeats * 30 * 4
        expected_fits = setting.repeats * 5
        checks.append(check(
            f"{setting.scenario}: prediction rows",
            len(predictions) == expected_predictions,
            f"observed {len(predictions)}, expected {expected_predictions}",
        ))
        checks.append(check(
            f"{setting.scenario}: unique predictions",
            not predictions.duplicated(["repeat_id", "model_id", "method"]).any(),
            "one row per repeat, target model, and method",
        ))
        counts = predictions.groupby(["repeat_id", "method"])["model_id"].nunique()
        checks.append(check(
            f"{setting.scenario}: target coverage",
            bool((counts == 30).all()),
            "all repeat-method cells contain 30 targets",
        ))
        observed_nunique = predictions.groupby(["repeat_id", "model_id"])["observed"].nunique()
        checks.append(check(
            f"{setting.scenario}: outcome agreement",
            bool((observed_nunique == 1).all()),
            "all methods share the same observed outcome",
        ))
        checks.append(check(
            f"{setting.scenario}: diagnostic rows",
            len(diagnostics) == expected_fits,
            f"observed {len(diagnostics)}, expected {expected_fits}",
        ))
        alpha_valid = diagnostics["selected_alpha"].round(10).isin(alpha_grid).all()
        checks.append(check(
            f"{setting.scenario}: alpha grid",
            bool(alpha_valid),
            "all selected alpha values belong to 0.0, 0.1, ..., 1.0",
        ))
        recomputed = []
        for (repeat_id, method), frame in predictions.groupby(["repeat_id", "method"]):
            recomputed.append(
                {
                    "repeat_id": repeat_id,
                    "method": method,
                    "mse": ((frame["observed"] - frame["predicted"]) ** 2).mean(),
                }
            )
        recomputed = pd.DataFrame(recomputed)
        recomputed_mean = recomputed.groupby("method")["mse"].mean().sort_index()
        reported_mean = summary.set_index("method")["mean_mse"].sort_index()
        checks.append(check(
            f"{setting.scenario}: MSE reconciliation",
            bool(np.allclose(recomputed_mean, reported_mean, rtol=0, atol=1e-12)),
            "method_summary mean MSE matches independent recomputation",
        ))

    audit = pd.DataFrame(checks)
    audit.to_csv(RESULTS / "audit_checks.csv", index=False)
    failed = audit[audit["status"] == "FAIL"]
    print(audit.to_string(index=False))
    if not failed.empty:
        raise RuntimeError(f"{len(failed)} full-experiment audit checks failed")


if __name__ == "__main__":
    main()
