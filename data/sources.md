# Sources of the compiled activation dataset

Each row of `activation_dataset.csv` is identified by `source_id` and a DOI or arXiv identifier. The tables below list the sources per host system, the measurement technique, and the number of rows contributed.

## Conventions

- `activation_pct` is taken from the source where it is stated. Rows whose note begins with *Derived* were computed from other quantities reported in the source (for example hole concentration divided by implanted concentration or dose); the formula is given in `activation_pct_note`.
- Activation is not defined identically across sources: most rows are Hall-based ratios of hole concentration to implanted concentration or dose, some are resistivity-derived with 100% taken as the maximum, and two rows are net values equal to 100% minus a fitted compensation ratio. These distinctions are recorded in `measurement_technique` and `activation_pct_note`.
- `value_basis` is one of `stated`, `derived`, `net_of_compensation`, `approximate` (qualitative wording, a bound, or a midpoint, entered as a number) or `none`; `text_basis` is `full_text` for Negoro et al., Saks et al., Weiße et al. (AIP Advances) and Breckenridge et al., and `abstract_or_partial` otherwise.
- Rows with `NA` in `anneal_temp_C` or `activation_pct` are kept as qualitative context and are not used in fitted models.
- The implant temperature of the `S1_S2` rows is not entered in the CSV; the analysis scripts set it to 500 C, the value given in the methods of the source.
- `hall_factor_convention` records whether a source states the Hall scattering factor used; `NA` means that no statement was found in the text that was compiled.

## SiC

| source_id | rows | measurement technique | DOI / arXiv |
|---|---|---|---|
| S1_S2 | 3 | Hall | arXiv:2001.08021;10.3390/ma12213468 |
| S3 | 3 | Hall | 10.1063/1.1796518 |
| S4_MSF | 2 | Hall+admittance | 10.4028/www.scientific.net/msf.924.184 |
| Nipoti2011 | 1 | resistivity+Hall-implied | 10.1143/apex.4.111301 |
| Weng2010 | 4 | SCM(not_Hall) | 10.4028/www.scientific.net/msf.645-648.713 |
| Saks2004 | 8 | Hall(resistivity-derived) | 10.1063/1.1764934 |
| Bluet2000 | 2 | Hall(implied) | 10.1063/1.1305904 |
| Weisse2019_AIPAdv | 2 | Hall+admittance(two-defect_model) | 10.1063/1.5096440 |
| Watanabe2008_shuttle | 1 | Hall+CV | 10.1143/jjap.47.3027 |
| Nipoti2013_JMR | 1 | Hall | 10.1557/jmr.2012.207 |

## GaN

| source_id | rows | measurement technique | DOI / arXiv |
|---|---|---|---|
| G1 | 1 | Hall(T-dep_Arrhenius) | 10.1063/1.5116866 |
| G2_Hirukawa | 2 | Hall | 10.35848/1882-0786/abf4f3 |
| G3 | 1 | Hall | 10.1049/el.2013.3214 |
| G4 | 1 | Hall | 10.1016/j.jcrysgro.2011.12.016 |
| Sumida2025_Ncoimp | 1 | Hall | 10.35848/1882-0786/ade3ac |
| Sumida2025_Ccap | 1 | Hall | 10.35848/1882-0786/ae0866 |
| Breckenridge2021 | 5 | Hall | 10.1063/5.0038628 |
| Sakurai2020_ISPSD | 1 | Hall | 10.1109/ispsd46842.2020.9170174 |
| Niwa2017_APEX | 1 | Hall | 10.7567/apex.10.091002 |
| Rocco2020_JAP | 1 | Hall+TLM | 10.1063/5.0016358 |

## Diamond

| source_id | rows | measurement technique | DOI / arXiv |
|---|---|---|---|
| Seki2019_APL | 1 | Hall | 10.1063/1.5111882 |
| Vogel2004_DRM | 1 | Hall | 10.1016/j.diamond.2004.04.005 |
| Fontaine1996_APL | 1 | Hall | 10.1063/1.115879 |
| Imamura2025_JAP | 2 | Raman+RBS+resistivity | 10.1063/5.0297918 |
| Seki2023_pssa | 1 | Hall | 10.1002/pssa.202300522 |
| Seki2021_JAP | 1 | Hall | 10.1063/5.0048309 |
| Ueda2008_DRM | 1 | Hall(resistivity-derived) | 10.1016/j.diamond.2007.12.017 |
| Tsubouchi2008_JJAP | 1 | Hall | 10.1143/jjap.47.7047 |

## GaAs

| source_id | rows | measurement technique | DOI / arXiv |
|---|---|---|---|
| McLevige1977_JAP | 1 | Hall+resistivity | 10.1063/1.324218 |
| Chambon1985_APL | 1 | Hall(differential) | 10.1063/1.95722 |
| Maezawa1986_EDL | 2 | Hall(77K) | 10.1109/edl.1986.26276 |
| Anderson1979_APL | 1 | Hall | 10.1063/1.91027 |
| Opyd1990_JAP | 1 | Hall(vdP)+ECV | 10.1063/1.344531 |

## InP

| source_id | rows | measurement technique | DOI / arXiv |
|---|---|---|---|
| Donnelly1979_APL | 1 | Hall | 10.1063/1.90573 |
| Dinia2000_SSC | 1 | Hall(carrier_profile) | 10.1016/s0038-1098(99)00511-6 |
| Shen1991_NIMB | 2 | Hall(ECV-derived) | 10.1016/0168-583x(91)96283-q |

## Si

| source_id | rows | measurement technique | DOI / arXiv |
|---|---|---|---|
| Landi1988_APA | 2 | Hall(sheet_resistance) | 10.1007/bf00615499 |
| Powell1984_JAP | 2 | Sheet_resistance | 10.1063/1.333818 |
| Narayan1983_APL | 1 | Hall+RBS | 10.1063/1.94200 |
| Saavedra2004_JAP | 1 | Hall+4pt_probe | 10.1063/1.1769095 |
