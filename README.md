# Cross-system electrical activation of implanted p-type dopants: literature-compiled dataset and analysis

Data and analysis scripts for the study of how implant temperature, anneal
temperature and anneal time relate to the measured electrical activation of
ion-implanted acceptors in Al/4H-SiC, Mg/GaN, B/diamond, Be/GaAs,
Zn,Be,Mg/InP and B/Si. The dataset is compiled entirely from published
literature (initial compilation: 67 rows, 40 independent sources; enlarged database: 137 rows, 51 sources); no new atomistic, DFT or
molecular-dynamics calculations were performed.

## Layout

| Path | Content |
|---|---|
| `common.py` | Data loading, standardization, exact source-level bootstrap, random-intercept (REML) fit |
| `data/activation_dataset.csv` | Compiled dataset, one row per reported process condition, with DOI or arXiv identifier, `value_basis` (stated, derived, net_of_compensation, approximate, none) and `text_basis` (full_text or abstract_or_partial) flags |
| `data/database.csv` | Enlarged database (initial compilation plus 70 verified rows from full-text papers) |
| `data/expansion_verified.csv` | The added rows with location in the paper |
| `data/reporting_audit.csv` | Reporting audit of the 13 full-text papers |
| `data/search_audit/` | Record of the structured literature search of 3 November 2026 (queries, hits, screening of records outside the database) |
| `explorer.html` | Static browser for the database |
| `data/system_covariates.csv` | One value per host system, used for the temperature-scale sensitivity analysis |
| `data/sources.md` | Sources per system and conventions used when compiling |
| `data/DATA_DICTIONARY.md` | Variables, units, missing-value convention, derivation and inclusion flags |
| `data/*_results.json` | Numerical output of the scripts |
| `figures/` | Figures produced by the scripts |

## Scripts

Run from the repository root, in this order. Each script reads from `data/` and writes to `data/` and `figures/`.

1. `01_sic_analysis.py` - SiC multivariate model of activation on standardized anneal and implant temperature; exact source-level bootstrap (all 256 ordered resamples of the 4 sources, with excluded resamples counted); row-level bootstrap for comparison; leave-one-source-out coefficients; random-intercept model; leave-one-source-out prediction error; thermal-budget covariate; per-source Hall-scattering-factor sensitivity (grid and Monte Carlo); 17-model specification panel with row-level permutation tests and leave-one-source-out sign checks; refits after removing sources with a different activation definition, the Saks et al. factorial alone, and a source fixed-effects model (`data/sic_results.json`).
2. `02_gan_analysis.py` - GaN models of activation on anneal temperature, and on anneal temperature plus anneal time, with exactly enumerated source-level resampling and leave-one-source-out refits (`data/gan_results.json`).
3. `03_cross_system.py` - pooled SiC/GaN correlation under four temperature scales with a source-level bootstrap, ceiling sweep, per-system summary, counts of rows by value basis and text basis, the Hall-factor reporting count by text basis, the pooled correlation excluding approximate rows, and figures (`data/cross_system_results.json`).

4. `04_build_database.py` - merges the initial compilation with `data/expansion_verified.csv` into `data/database.csv` (laboratory group, year, activation-definition class, assumed implant temperature flag).
5. `05_hierarchical_model.py` - Bayesian linear mixed model of logit activation with source and laboratory intercepts for Al-implanted SiC, with leave-one-source-out and prior sensitivity analyses (`data/hierarchical_results.json`, `figures/design_map.png`).
6. `06_kinetic_validation.py` - Arrhenius-equivalent-time model, and leave-one-source-out, leave-one-laboratory-out and temporal evaluation of five predictors (`data/kinetic_validation_results.json`).
7. `07_method_dependence_audit.py` - difference between two analysis methods applied to the same samples, and the reporting audit of the full-text papers (`data/method_audit_results.json`).
8. `08_build_explorer.py` - writes `explorer.html`, a static page for browsing and filtering the database.

9. `09_graphical_abstract.py` - draws the graphical abstract from the saved results (`figures/graphical_abstract.png`).
10. `10_sensitivity_analyses.py` - clipping bounds, variance priors, restricted subsets by definition of activation, measurement-error scenarios and the row counts at each selection step (`data/sensitivity_results.json`).
11. `11_dataset_flowchart.py` - draws the dataset-selection flowchart (`figures/dataset_flowchart.png`).

12. `12_revision_analyses.py` - effect definitions and baseline covariates of the hierarchical model, refits after removing each definition class, variance share by structure and prior, fold-level validation summary with covariate differences, and the paired comparison of the two fitting methods (`data/revision_v2_results.json`, `figures/method_pairs.pdf`).

Suggested order: 01, 04, 05, 10, 06, 12, 07, 02, 03, 11, 09, 08. Scripts 05, 06 and 10 use four Gibbs chains with seeds 1-4 and dispersed starting values; the bootstrap over held-out units uses seed 0.

`common.py` and `models.py` hold shared helpers (data loading, exact source-level bootstrap, REML and Gibbs mixed-model fits).

## Requirements

Figures are drawn in a Times-type font. If Times New Roman, Times or Nimbus Roman is not installed, point `FIGURE_FONT_DIR` to a folder containing TeX Gyre Termes `.otf` files; otherwise matplotlib falls back to DejaVu Serif.

Python 3.11.16 with the packages in `requirements.txt` (versions used: numpy 2.4.6, pandas 2.3.3, scipy 1.17.1, scikit-learn 1.9.1, matplotlib 3.11.1). The repository is tagged v1.1 for the analysis reported in the manuscript.

## Data notes

- `activation_pct` is the value stated in the source where one exists; rows marked *Derived* in `activation_pct_note` were computed from other reported quantities.
- `value_basis` marks each activation value as stated by the source, derived by us, net of compensation, approximate (qualitative wording, a bound, or a midpoint, entered as a number), or absent; `text_basis` marks whether the full text was read (17 rows) or only an abstract or partly accessible text (50 rows).
- Activation is not defined identically across sources (Hall ratios to implanted concentration or dose, resistivity-derived values, and net values equal to 100% minus a fitted compensation ratio); see `measurement_technique` and `activation_pct_note`.
- All values come from the public, peer-reviewed or preprint sources identified by the DOI or arXiv column; no proprietary data are used.
- The 1175 C anneal of row 2 (Spera et al.) is the value printed in the source abstract and may be a typographical error for 1775 C; the row has no activation value and is not used in any fit.
