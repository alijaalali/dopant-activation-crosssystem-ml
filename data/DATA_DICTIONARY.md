# Data dictionary for `database.csv`

One row is one reported process condition (or one sample) from one source. An empty cell means that the value is not reported by the source or is not applicable; no value is imputed. In `activation_dataset.csv` (the initial compilation) missing values are written `NA`.

| Column | Unit / values | Meaning and derivation | Source |
|---|---|---|---|
| `db_id` | integer | Row number in this file (1-137) | file |
| `origin` | `initial_compilation`, `expansion` | Initial 67-row compilation, or one of the 70 rows read from full-text papers and verified against the source table | file |
| `system` | SiC, GaN, Diamond, GaAs, InP, Si | Host-material system | source |
| `host` | text | Host material and polytype | source |
| `dopant` | text | Implanted acceptor species | source |
| `source_id` | text | Key of the independent source; the paired publications S1 and S2 report the same samples and share one key | file |
| `lab_group` | text | Research group or institution of the source; used as the laboratory grouping factor | assigned by us |
| `year` | year | Publication year | source |
| `doi_or_arxiv` | text | DOI or arXiv identifier of the source | source |
| `implant_temp_C` | deg C | Substrate temperature during implantation | source |
| `anneal_temp_C` | deg C | Anneal temperature | source |
| `anneal_time_min` | min | Anneal time | source |
| `anneal_method` | text | Furnace, microwave, ultra-high-pressure, rapid thermal, etc. | source |
| `cap` | text | Capping layer or ambient during the anneal | source |
| `co_implant` | text | Co-implanted species (N: none) | source |
| `measurement_technique` | text | Technique that produced the electrical value | source |
| `dose_cm2` | cm^-2 | Implanted dose | source |
| `implanted_conc_cm3` | cm^-3 | Implanted (box or peak) concentration as stated in the source | source |
| `hole_conc_cm3` | cm^-3 | Measured hole density | source |
| `sheet_hole_cm2` | cm^-2 | Sheet hole density (not populated in this release) | source |
| `mobility_cm2Vs` | cm^2 V^-1 s^-1 | Hole mobility | source |
| `EA_meV` | meV | Acceptor ionization energy reported or fitted by the source | source |
| `activation_pct` | percent | Electrical activation as defined by the source, or derived by us (see `value_basis` and `definition_class`) | source or derived |
| `value_basis` | `stated`, `derived`, `net_of_compensation`, `approximate`, `none` | How `activation_pct` was obtained: printed by the source; computed by us from other printed quantities; equal to 100% minus a fitted compensation ratio; a qualitative statement or bound entered as a number; absent | assigned by us |
| `text_basis` | `full_text`, `abstract_or_partial` | Whether the full text was read | assigned by us |
| `definition_class` | `stated_hall`, `fit_acceptor_fraction`, `net_of_compensation`, `derived_ratio`, `cv_activated_fraction`, `resistivity_adjusted`, `approximate`, `none` | Physical definition of the activation value (Table 1 of the manuscript) | assigned by us |
| `hall_factor_convention` | text | Whether and how the source states the Hall scattering factor | source |
| `implant_T_is_assumed` | True/False | Implant temperature was fixed from the methods of the paper, not given per sample | assigned by us |
| `two_step_anneal` | True/False | Two-step anneal; excluded from the hierarchical and validation analyses | assigned by us |
| `implant_T_is_range` | True/False | Source gives a range of implant temperatures; the midpoint is used | assigned by us |
| `compensation_ratio` | as reported | Compensation ratio reported or fitted by the source | source |
| `sheet_resistance_ohm_sq` | ohm per square | Sheet resistance | source |
| `ionization_correction` | text | Whether incomplete ionization was corrected for, and how | source |
| `overlap_note` | text | Samples that also appear in another source | assigned by us |
| `method_variant` | `p`, `p-mu` | Fit to the hole density alone, or to hole density and mobility (Nipoti et al. 2022) | source |
| `location_in_paper` | text | Table, figure or page of the value (expansion rows) | assigned by us |
| `notes` | text | Free text, including the verbatim evidence for the row | assigned by us |

## Inclusion flags used by the analyses

- Initial multivariate model (Section 3.2): `system` = SiC, rows of `activation_dataset.csv` with numeric anneal temperature, implant temperature and activation (15 rows, 4 sources).
- Hierarchical model and validation (Sections 3.3-3.5): `system` = SiC, `dopant` = Al, numeric `activation_pct`, `implant_temp_C`, `anneal_temp_C`, `anneal_time_min` and `implanted_conc_cm3`, `definition_class` not `approximate`, `two_step_anneal` False (32 rows, 6 sources).
- Exact counts at each step are in `sensitivity_results.json` (`dataset_flow`) and in Figure 1 of the manuscript.
