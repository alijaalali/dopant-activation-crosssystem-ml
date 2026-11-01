# Cross-system electrical activation of implanted p-type dopants: literature-compiled dataset and analysis

Data and analysis scripts for the study of how implant temperature, anneal
temperature and anneal time relate to the measured electrical activation of
ion-implanted acceptors in Al/4H-SiC, Mg/GaN, B/diamond, Be/GaAs,
Zn,Be,Mg/InP and B/Si. The dataset is compiled entirely from published
literature (67 rows, 40 independent sources); no new atomistic, DFT or
molecular-dynamics calculations were performed.

## Layout

| Path | Content |
|---|---|
| `data/activation_dataset.csv` | Compiled dataset, one row per reported process condition, with DOI or arXiv identifier |
| `data/system_covariates.csv` | One value per host system, used to normalize anneal temperature across systems |
| `data/sources.md` | Sources per system and conventions used when compiling |
| `data/*.json` | Numerical output of the scripts |
| `figures/` | Figures produced by the scripts |

## Scripts

Run from the repository root, in this order (`python 01_sic_gan_baseline.py`, and so on). Every script reads from `data/` and writes to `data/` or `figures/` next to itself.

1. `01_sic_gan_baseline.py` - SiC multivariate model of activation on standardized anneal and implant temperature, with row-level bootstrap intervals and leave-one-source-out cross-validation; pooled SiC/GaN correlation with raw and system-normalized anneal temperature; Hall-scattering-factor reporting count and sensitivity bound.
2. `02_sic_bootstrap_robustness.py` - 17-model-family panel with permutation tests, source-level (cluster) bootstrap, thermal-budget covariate, and Hall-factor sensitivity re-fit.
3. `03_dataset_summary_figures.py` - per-system dataset summary, SiC fit on the full set, dose-covariate check, model panel, and the SiC figures.
4. `04_final_sic_model_panel.py` - SiC multivariate model, cluster bootstrap, leave-one-source-out error, and model panel on the final dataset (`data/sic_final_model_results.json`).
5. `05_six_system_summary.py` - per-system row, source and anneal-temperature summary and fittability check (`data/six_system_summary.json`).
6. `06_six_system_figure.py` - anneal-temperature range per system (`figures/sixsystem_comparison.*`).
7. `07_gan_models.py` - GaN models of activation on anneal temperature and on anneal temperature plus anneal time with source-level bootstrap (`data/gan_model_results.json`); reads `data/sic_final_model_results.json`.

## Requirements

Python 3.11 with the packages in `requirements.txt` (versions used: numpy 2.4, pandas 2.3, scipy 1.17, scikit-learn 1.9, matplotlib 3.11).

## Data notes

- `activation_pct` is the value stated in the source where one exists; rows marked *Derived* in `activation_pct_note` were computed from other reported quantities.
- Activation is not defined identically across sources (Hall ratios to implanted concentration or dose, resistivity-derived values, and net values equal to 100% minus a fitted compensation ratio); see `measurement_technique` and `activation_pct_note`.
- All values come from the public, peer-reviewed or preprint sources identified by the DOI or arXiv column; no proprietary data are used.
