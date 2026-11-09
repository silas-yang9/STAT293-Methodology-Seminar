# GDSC gemcitabine pilot results

## Status

The pilot completed successfully. All 50 repeat-by-fold model fits have
checkpoint files and the workflow resumes without refitting completed folds.

## Analysis set

- Target: 30 pancreatic-carcinoma cell lines.
- External: 199 cell lines: colorectal carcinoma (47), gastric carcinoma
  (28), non-small-cell lung carcinoma (84), and ovarian carcinoma (40).
- Four additional non-small-cell lung carcinoma models were excluded: two
  official endpoints returned HTTP 403, while two files used a legacy schema
  without the harmonized Ensembl/RSEM columns required by the common-feature
  pipeline. No target model was excluded.
- Outcome: GDSC2 gemcitabine `LN_IC50`.
- Features: 1,000 genes selected by expression prevalence and variance using
  only the 199 external expression profiles.

## Evaluation

- Ten repeated five-fold splits of the 30 pancreatic-carcinoma models.
- Every fold used 199 external observations and 24 pancreatic-carcinoma
  training observations; the other six pancreatic-carcinoma observations were
  used only for testing.
- Three-fold inner cross-validation within cancer groups.
- `alpha` grid: 0, 0.25, 0.50, 0.75, 1.00.
- `lambda.1se` used for the overall and second-stage lasso fits.

## Repeated-CV results

| Method | Mean MSE | SD across repeats | Mean MAE | Mean correlation | Mean R-squared |
|---|---:|---:|---:|---:|---:|
| Overall Lasso | 5.531 | 0.372 | 1.943 | 0.421 | 0.120 |
| Pretrained Lasso | 5.538 | 0.434 | 1.942 | 0.383 | 0.119 |
| Target-training mean | 6.759 | 0.275 | 2.157 | -0.346 | -0.075 |
| Individual Lasso | 6.771 | 0.319 | 2.158 | -0.286 | -0.077 |

Pretrained Lasso improved on Individual Lasso in all ten repeats. The mean
paired MSE difference, `MSE(Individual) - MSE(Pretrained)`, was 1.233, an 18.2%
reduction relative to the mean Individual-Lasso MSE. Overall and Pretrained
Lasso were effectively tied: their mean MSEs differed by 0.008, and each had
the lower repeat-level MSE in five of the ten repeats.

The selected `alpha` values across 50 outer fits were 0 in 13 fits (26%), 0.25
in 21 (42%), 0.50 in 10 (20%), 0.75 in three (6%), and 1 in three (6%). Thus,
the validation procedure usually selected a model closer to the overall end
of the Overall--Individual continuum.

## Sparsity and stability

The conservative target-only fit selected no gene in 44 of 50 fits, explaining
why its predictions were usually the target-training mean. Overall and
Pretrained Lasso each selected a median of 4.5 genes per fit and had identical
support sets in this pilot; their predictions differed through the pretrained
offset and coefficient fitting. The most frequently selected genes included
`SLFN11` and `WNT10A` (43/50 fits each), followed by `EVA1A` (29/50), `TBX18`
(27/50), and `SGK2` (24/50). These are predictive stability results, not
validated causal or clinical biomarkers.

## Interpretation and limitations

This pilot supports the narrow conclusion that external cancer-cell-line data
improved pancreatic-cell-line prediction relative to a conservative
target-only lasso. It does not show that Pretrained Lasso is better than simple
pooling: Overall Lasso and Pretrained Lasso were indistinguishable at the
current resolution.

The target sample contains only 30 cell lines, repeated splits reuse the same
models, the alpha grid is coarse, and there is no independent external target
cohort. Results are conditional on the selected external cancer types, the
external-only variance filter, GDSC2 assay conditions, and `lambda.1se`.
Cell-line drug sensitivity must not be described as patient treatment response.
