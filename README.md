# STAT 293 Methodology Seminar

This repository contains the STAT 293 course project on Pretrained Lasso for
grouped high-dimensional prediction. It includes an article-faithful method
summary, three simulation experiments, and a real-data application predicting
GDSC2 gemcitabine sensitivity from RNA-seq expression. The journal-submission
project is maintained separately and is not part of this repository.

## Start here

- Final report source: `report/FinalReport.Rmd`
- Knitted report: `output/pdf/FinalReport.pdf`
- Simulation code and results: `simulation_study/`
- GDSC real-data analysis: `gdsc_real_data/`
- Data-cleaning R Markdown: `gdsc_real_data/01_DataCleaning.Rmd`

## Repository structure

```text
.
|-- 00_requirements.R
|-- python_requirements.txt
|-- report/
|   `-- FinalReport.Rmd
|-- output/pdf/
|   `-- FinalReport.pdf
|-- simulation_study/
|   |-- scripts/
|   |-- report/
|   `-- results/
`-- gdsc_real_data/
    |-- 01_DataCleaning.Rmd
    |-- scripts/
    |-- raw/
    |-- processed/
    `-- results/
```

Local package libraries, rendering tools, temporary files, archived homework,
legacy TCGA attempts, and copyrighted reference PDFs are intentionally excluded
from the public repository.

## Software requirements

The analysis was developed with R 4.4.2 and Python 3.12.4. From the repository
root, install the recorded dependencies with:

```powershell
Rscript 00_requirements.R
python -m pip install -r python_requirements.txt
```

The R setup installs project packages into the ignored local `.Rlib/`
directory. The Pretrained Lasso implementation is installed from the authors'
public `erincr/ptLasso` repository.

## Reproduce the simulation experiments

Run the following commands from the repository root:

```powershell
Rscript simulation_study/scripts/11_run_Figure4.R --reps=100
Rscript simulation_study/scripts/12_run_Figure6.R --reps=100
Rscript simulation_study/scripts/13_run_ParameterSensitivity.R --reps=100
```

The three runners save raw results, summaries, figures, settings, and session
information under separate subdirectories of `simulation_study/results/`.
Resumable checkpoints are generated locally during long runs but are excluded
from Git because the consolidated raw and summary results are retained.

## Reproduce the GDSC real-data experiment

The data sources, cleaning decisions, variables, exclusion audit, and feature
selection procedure are documented in `gdsc_real_data/01_DataCleaning.Rmd`.
The main processing and analysis order is:

```powershell
python gdsc_real_data/scripts/01_prepare_manifest.py
python gdsc_real_data/scripts/02_download_rnaseq.py
python gdsc_real_data/scripts/09_build_full_matrix.py
Rscript gdsc_real_data/scripts/10_run_full_experiment.R
python gdsc_real_data/scripts/11_audit_full_experiment.py
Rscript gdsc_real_data/scripts/12_make_full_figures.R
Rscript gdsc_real_data/scripts/15_run_conservative_alpha_experiment.R
python gdsc_real_data/scripts/16_audit_conservative_alpha.py
Rscript gdsc_real_data/scripts/17_make_conservative_figures.R
```

The primary target group contains 30 pancreatic-cancer cell lines. The final
external dataset contains 199 colorectal, gastric, non-small-cell lung, and
ovarian cancer cell lines. The outcome is GDSC2 `LN_IC50`; predictors are
`log2(TPM + 1)` expression values selected without using target outcomes or
target expression measurements.

## Render the final report

After the saved experiment outputs are present, run:

```powershell
Rscript -e "rmarkdown::render('report/FinalReport.Rmd', output_file='FinalReport.pdf', output_dir='output/pdf')"
```

The report checks that required result files exist before rendering. The
real-data application is an in-vitro cell-line prediction study, not a clinical
validation analysis.
