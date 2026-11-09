# Full GDSC gemcitabine experiment

Completed on 2026-11-04. This file records the analysis choices and the main numerical results used in the course report.

## Data

- Target: 30 pancreatic-cancer cell lines.
- External data: 199 cell lines: 47 colorectal, 28 gastric, 84 non-small-cell lung, and 40 ovarian.
- Four non-small-cell lung models were excluded before modelling because the public RNA-seq file was unavailable or failed schema validation. No target model was excluded.
- Outcome: GDSC2 gemcitabine `LN_IC50`; lower values mean greater sensitivity.
- Predictors: Cell Model Passports RNA-seq, transformed as `log2(TPM + 1)`.
- Feature screening used external samples only. Genes required TPM above 1 in at least 20% of external models and were ranked by external-sample variance.

## Evaluation

- Primary analysis: 1,000 genes, `lambda.1se`, 20 repeated five-fold target-domain evaluations.
- Four sensitivity analyses: 500 genes, 2,000 genes, AUC outcome, and `lambda.min`; each used 10 repeated five-fold evaluations.
- Inner tuning: three folds assigned within cancer groups; alpha grid 0, 0.1, ..., 1.
- The first stage used all external models and only the pancreatic training folds. Held-out pancreatic models were never used in fitting or tuning.
- Comparators: target-only Individual Lasso, pooled Overall Lasso, and grouped Pretrained Lasso. A target-training mean was retained as a reference.

## Primary results

| Method | Mean MSE | SD MSE | Mean MAE | Mean correlation | Mean R-squared |
|---|---:|---:|---:|---:|---:|
| Overall | 5.619 | 0.535 | 1.948 | 0.372 | 0.106 |
| Pretrained | 5.706 | 0.513 | 1.972 | 0.315 | 0.093 |
| Individual | 6.821 | 0.369 | 2.168 | -0.311 | -0.085 |
| Target mean | 6.825 | 0.353 | 2.167 | -0.367 | -0.085 |

Pretrained Lasso reduced mean MSE by 16.3% relative to Individual Lasso and won all 20 paired repetitions. Overall Lasso was slightly better than Pretrained Lasso on average and won 15 of 20 paired repetitions.

## Sensitivity results

| Analysis | Overall MSE | Pretrained MSE | Individual MSE |
|---|---:|---:|---:|
| 500 genes; LN_IC50; lambda.1se | 5.933 | 6.107 | 6.944 |
| 2,000 genes; LN_IC50; lambda.1se | 6.221 | 6.090 | 6.820 |
| 1,000 genes; AUC; lambda.1se | 0.01943 | 0.02011 | 0.02266 |
| 1,000 genes; LN_IC50; lambda.min | 4.382 | 5.114 | 6.948 |

Pretrained Lasso beat Individual Lasso in all ten paired repetitions of every sensitivity analysis. Its incremental advantage over Overall Lasso was not stable: only the 2,000-gene analysis favored Pretrained on average.

## Variable stability

Across 100 primary outer fits, the most frequently selected effective Pretrained coefficients were:

| Gene | Selection frequency | Median coefficient | Sign when selected |
|---|---:|---:|---|
| SLFN11 | 0.85 | -0.476 | Negative |
| WNT10A | 0.85 | 0.240 | Positive |
| EVA1A | 0.59 | 0.111 | Positive |
| TBX18 | 0.59 | -0.074 | Negative |
| SGK2 | 0.52 | 0.083 | Positive |
| KLK5 | 0.46 | 0.020 | Positive |

The negative SLFN11 coefficient is concordant with published cell-line evidence that higher SLFN11 expression is associated with greater gemcitabine sensitivity. This experiment does not establish any variable as a clinical biomarker.

## Interpretation

The most defensible conclusion is that external cancer groups improved pancreatic target-domain prediction relative to target-only fitting. The experiment did not show a stable advantage of Pretrained Lasso over simple pooling. Performance remained modest, and all claims are limited to in-vitro cell-line drug sensitivity.

## Refined conservative alpha selection

A subsequent refinement addressed the imbalance between 199 external models and only 24 pancreatic training models in each outer split. The original global alpha minimizes prevalidated error across all training observations and is therefore dominated by external groups. The refined rule uses only pancreatic training-fold error and selects the smallest alpha whose error is within one standard error of the pancreatic minimum. Choosing the smallest eligible alpha deliberately favors more sharing and reduces unstable target-specific adjustment.

The confirmation experiment fixed 2,000 external-selected genes and used `lambda.1se` in both lasso stages. It compared global alpha, pancreatic minimum-error alpha, and pancreatic one-standard-error alpha in 20 new repeated five-fold evaluations. All alpha choices were made inside each outer training set.

| Method | Mean MSE | SD MSE | Mean MAE | Mean correlation | Mean R-squared |
|---|---:|---:|---:|---:|---:|
| Pretrained, target one-SE alpha | 5.983 | 0.720 | 2.013 | 0.268 | 0.049 |
| Pretrained, target minimum alpha | 6.055 | 0.659 | 2.035 | 0.252 | 0.037 |
| Pretrained, global alpha | 6.062 | 0.512 | 2.028 | 0.197 | 0.036 |
| Overall | 6.136 | 0.412 | 2.040 | 0.143 | 0.024 |
| Individual | 6.752 | 0.344 | 2.156 | -0.277 | -0.074 |

The conservative Pretrained rule had the lowest mean MSE. It beat Overall in 14 of 20 repetitions and Individual in 18 of 20. Its mean MSE was 2.5% lower than Overall and 11.4% lower than Individual. This supports describing Pretrained Lasso as the best average method under the refined rule, but not as decisively or uniformly superior.

A broader exploratory grid that simultaneously selected 1,000 versus 2,000 genes and `lambda.1se` versus `lambda.min` was also tested for five development repetitions. That procedure over-selected `lambda.min` and was unstable; Overall remained best on average. It was therefore not adopted as the refined analysis. The failed development result is retained under `results/nested_tuning_experiment/` rather than being discarded.

