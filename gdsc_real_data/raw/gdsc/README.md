# GDSC and Cell Model Passports source files

These files are unchanged public inputs for the STAT 293 real-data experiment.
They were downloaded and audited on October 31, 2026.

| File | Purpose | SHA-256 |
|---|---|---|
| `GDSC2_fitted_dose_response_27Oct23.xlsx` | GDSC2 fitted dose-response data, including gemcitabine `LN_IC50` and AUC | `0E99120A14A003DCBD3C27C020F1CFD35E38EE091067A2406F66C0D4DBC1C899` |
| `model_list_20261013.csv` | Cell Model Passports identifiers and cancer annotations | `5E2868D7C063DFD83B8965CBAB60F7531F5089968107999E316A2A38E07DA6D9` |

Official sources:

- GDSC and Cell Model Passports downloads: <https://cellmodelpassports.sanger.ac.uk/downloads/gdsc>
- Drug-sensitivity documentation: <https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/>
- Expression-data documentation: <https://depmap.sanger.ac.uk/documentation/datasets/expression/>

Per-model RNA-seq archives are stored in the adjacent `raw/rnaseq/` directory
and use the official endpoint pattern
`https://cog.sanger.ac.uk/cmp/download/<SANGER_MODEL_ID>_rnaseq.csv.zip`.
