# Simulation study

This directory contains the R simulation experiments for the STAT 293 course
project. It does not contain the GDSC real-data analysis.

## Experiments

| Experiment | Purpose | Runner | Output directory |
|---|---|---|---|
| Experiment 1 (paper Figure 4) | Assess cross-validation error under different fold counts and lambda rules | `scripts/11_run_Figure4.R` | `results/figure4/` |
| Experiment 2 (paper Figure 6) | Assess the value of external pretraining as source-target sharing changes | `scripts/12_run_Figure6.R` | `results/figure6/` |
| Experiment 3 | Vary dimension, SNR, and sparsity one factor at a time | `scripts/13_run_ParameterSensitivity.R` | `results/parameter_sensitivity/` |

All experiments use 100 Monte Carlo repetitions in the final report. The
runners use fixed seeds and checkpoint their raw CSV output, so an interrupted
run can be resumed with the same command and output directory.

## Main files

- `scripts/00_requirements.R`: loads the project-local R packages.
- `scripts/01_funct_PretrainedLassoSimulations.R`: shared data-generation and
  fitting functions.
- `scripts/02_funct_ParameterSensitivity.R`: metrics, summaries, and plots for
  Experiment 3.
- `scripts/11_run_Figure4.R`, `scripts/12_run_Figure6.R`, and
  `scripts/13_run_ParameterSensitivity.R`: reproducible runners.
- `results/`: raw results, summaries, settings, figures, and session records.

## Run from the repository root

```powershell
Rscript simulation_study/scripts/11_run_Figure4.R --reps=100
Rscript simulation_study/scripts/12_run_Figure6.R --reps=100
Rscript simulation_study/scripts/13_run_ParameterSensitivity.R --reps=100
```

Experiment 3 uses the baseline `p=500`, `SNR=1.8`, and `s=10`, then separately
evaluates `p={250,500,1000}`, `SNR={0.8,1.8,3.0}`, and `s={5,10,20}` at
`rho=0.5`. It compares Individual, Overall, and Pretrained Lasso using test MSE,
MAE, TPR, FDP, support F1, and selected-model size.
