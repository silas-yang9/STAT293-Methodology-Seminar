# GDSC gemcitabine real-data experiment

This directory contains the reproducible real-data analysis described in the
course final report. It is separate from the journal-submission workspace.
The complete cleaning and processing narrative is in `01_DataCleaning.Rmd`.
The official GDSC2 workbook and Cell Model Passports model list used to build
the manifest are retained under `raw/gdsc/`; per-model expression archives are
retained under `raw/rnaseq/`.

## Primary design

- Target group: Pancreatic Carcinoma (30 GDSC2 models).
- External groups: Colorectal Carcinoma, Gastric Carcinoma, Ovarian
  Carcinoma, and Non-Small Cell Lung Carcinoma.
- Outcome: GDSC2 `LN_IC50` for gemcitabine (Drug ID 1190).
- Predictors: `log2(rsem_tpm + 1)` from Cell Model Passports RNA-seq.
- Feature set: 1,000 genes with the largest variance among external samples
  after requiring TPM > 1 in at least 20% of external samples. Target outcomes
  and target expression are not used to select the feature set.
- Evaluation: 20 repeated five-fold evaluations defined only among pancreatic-cancer models.
- Comparators: Individual Lasso, Overall Lasso, and Pretrained Lasso.
- Alpha grid: 0, 0.1, ..., 1.0.
- Inner cross-validation: three folds, assigned within each cancer group.
- Lasso choice: `lambda.1se` for the overall and second-stage fits.

Sensitivity analyses use 500 and 2,000 genes, AUC instead of `LN_IC50`, and
`lambda.min` instead of `lambda.1se`. See `FULL_EXPERIMENT_RESULTS.md` for the
main numerical results and interpretation.

## Files

- `scripts/01_prepare_manifest.py`: creates the eligible model manifest from
  the local official GDSC2 workbook and Cell Model Passports model list.
- `scripts/02_download_rnaseq.py`: resumable RNA-seq downloader with ZIP
  integrity checking and a status file after every model.
- `scripts/03_build_matrix.py`: constructs the feature matrix and records the
  external-only feature selection.
- `scripts/04_run_pilot.R`: runs the fold-wise models and saves a checkpoint
  after each outer fold.
- `scripts/05_audit_results.py`: adds the target-training-mean baseline and
  paired error audit without refitting.
- `scripts/06_rebuild_supports.R`: reconstructs target support sets from the
  exact fitted coefficients used at `lambda.1se`.
- `scripts/07_run_repeated_pilot.R`: runs ten repeated five-fold target-domain
  evaluations with a separate checkpoint for every repeat and fold.
- `scripts/08_annotate_features.py`: adds gene symbols to the externally
  selected Ensembl features used in stability summaries.
- `scripts/09_build_full_matrix.py`: constructs the 2,000-gene matrix used by
  the primary and feature-count sensitivity analyses.
- `scripts/10_run_full_experiment.R`: runs the primary analysis and four
  sensitivity analyses with resumable checkpoints.
- `scripts/11_audit_full_experiment.py`: independently checks sample counts,
  target coverage, duplicate keys, alpha values, and MSE reconciliation.
- `scripts/12_make_full_figures.R`: creates the report figures from saved
  predictions, coefficients, and summaries.
- `scripts/13_run_nested_tuning_experiment.R`: evaluates a wider nested grid
  for feature count and second-stage regularization. This exploratory approach
  was unstable and is retained as a negative development result.
- `scripts/14_audit_nested_tuning.py`: audits the wider nested-grid experiment.
- `scripts/15_run_conservative_alpha_experiment.R`: evaluates global,
  target-minimum, and conservative target one-standard-error alpha rules with
  2,000 genes and `lambda.1se` in both stages.
- `scripts/16_audit_conservative_alpha.py`: independently verifies the
  conservative alpha rule and all prediction summaries.
- `scripts/17_make_conservative_figures.R`: creates the paired-MSE figure for
  the refined analysis.
- `raw/`: downloaded per-model RNA-seq ZIP files.
- `processed/`: analysis manifest, feature matrix, and feature-selection log.
- `results/full_experiment/`: tracked outer-fold assignments, consolidated
  predictions, coefficient summaries, audit checks, figures, and session
  information for the completed experiment. Per-fold checkpoint files are
  generated locally and intentionally excluded from Git.
- `results/nested_tuning_experiment/`: the retained wider-grid development
  experiment.
- `results/conservative_alpha_experiment/`: the 20-repeat refined analysis in
  which target one-standard-error alpha selection gives Pretrained Lasso the
  lowest average MSE.

The experiment is not a clinical validation study. Cell lines, rather than patients,
are the unit of analysis.

Models without a harmonized public RNA-seq ZIP are excluded before expression
processing. The exclusion list and server/schema-validation reason are retained in
`processed/excluded_missing_rnaseq.csv`; no target model may be removed by this
rule.
