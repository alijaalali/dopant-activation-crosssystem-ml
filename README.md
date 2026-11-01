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
| `common.py` | Data loading, standardization, exact source-level bootstrap, random-intercept (REML) fit |
| `data/activation_dataset.csv` | Compiled dataset, one row per reported process condition, with DOI or arXiv identifier |
| `data/system_covariates.csv` | One value per host system, used for the temperature-scale sensitivity analysis |
| `data/sources.md` | Sources per system and conventions used when compiling |
| `data/*_results.json` | Numerical output of the scripts |
| `figures/` | Figures produced by the scripts |

## Scripts

Run from the repository root, in this order. Each script reads from `data/` and writes to `data/` and `figures/`.

1. `01_sic_analysis.py` - SiC multivariate model of activation on standardized anneal and implant temperature; exact source-level bootstrap (all 256 ordered resamples of the 4 sources, with excluded resamples counted); row-level bootstrap for comparison; leave-one-source-out coefficients; random-intercept model; leave-one-source-out prediction error; thermal-budget covariate; per-source Hall-scattering-factor sensitivity (grid and Monte Carlo); 17-model specification panel with row-level permutation tests and leave-one-source-out sign checks (`data/sic_results.json`).
2. `02_gan_analysis.py` - GaN models of activation on anneal temperature, and on anneal temperature plus anneal time, with exactly enumerated source-level resampling and leave-one-source-out refits (`data/gan_results.json`).
3. `03_cross_system.py` - pooled SiC/GaN correlation under four temperature scales with a source-level bootstrap, ceiling sweep, per-system summary and figures (`data/cross_system_results.json`).

## Requirements

Python 3.11 with the packages in `requirements.txt` (versions used: numpy 2.4, pandas 2.3, scipy 1.17, scikit-learn 1.9, matplotlib 3.11).

## Data notes

- `activation_pct` is the value stated in the source where one exists; rows marked *Derived* in `activation_pct_note` were computed from other reported quantities.
- Activation is not defined identically across sources (Hall ratios to implanted concentration or dose, resistivity-derived values, and net values equal to 100% minus a fitted compensation ratio); see `measurement_technique` and `activation_pct_note`.
- All values come from the public, peer-reviewed or preprint sources identified by the DOI or arXiv column; no proprietary data are used.
- The 1175 C anneal of row 2 (Spera et al.) is the value printed in the source abstract and may be a typographical error for 1775 C; the row has no activation value and is not used in any fit.
