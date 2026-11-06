# Literature search record (not a systematic review)

This is a transparent record of the literature search used to check coverage of the extraction database. We do not claim a systematic review (no protocol registration, no dual screening, no full-text screening of the search hits).

## Search details

- Search date: 2026-11-03 (all six queries run the same day)
- Method: structured, relevance-ranked semantic literature search over Semantic Scholar / OpenAlex / PubMed records
- Settings: first 100 results per query, no filters (no year, journal, study-type, open-access or domain filters)
- Fields returned per hit: title, authors, year, citation count, journal, abstract (or abstract excerpt), DOI (where available)
- Raw search outputs were saved; parsed into `search_hits.csv` (one row per query and hit)

### Queries (exact strings)

| ID | Query | Hits returned |
|---|---|---|
| Q1 | Al ion implantation 4H-SiC p-type electrical activation annealing Hall effect | 100 |
| Q2 | aluminum implanted 4H-SiC implantation temperature effect on activation and sheet resistance | 100 |
| Q3 | Mg ion implantation GaN p-type activation annealing hole concentration | 100 |
| Q4 | ion implantation acceptor activation fraction hole concentration versus anneal temperature wide bandgap semiconductor | 100 |
| Q5 (secondary) | boron implanted diamond p-type electrical activation annealing | 100 |
| Q6 (secondary) | high-temperature or microwave annealing of implanted SiC activation of acceptors | 99 |

Total hit rows: 599. Unique hits (deduplicated by lowercase DOI; 4 hits without a DOI deduplicated by title): 408. Of these, 285 were returned by one query only, 81 by two, 20 by three, 20 by four, and 2 by five or more query appearances (one paper appears twice within the same query result list). Near-duplicate records of the same work with different DOIs or years may still be counted separately (for example, a "comparative study of high-temperature aluminium post-implantation annealing in 6H- and 4H-SiC" appears under two DOIs).

## Recall check against the extraction database

Database: `code_release/data/database.csv`, 137 rows from 51 sources, each with a DOI (`doi_or_arxiv` column; matched by lowercase DOI). 40 sources are from the initial compilation and 11 from the expansion step. `data/expansion/existing_sources.txt` lists 40 source-to-DOI pairs; `reporting_audit.csv` was read for source IDs only.

| Query | Database sources retrieved (of 51) |
|---|---|
| Q1 | 14 |
| Q2 | 13 |
| Q3 | 10 |
| Q4 | 16 |
| Q5 | 13 |
| Q6 | 13 |
| Any query | 38 (74.5%) |

By origin: initial compilation 29 of 40; expansion 9 of 11. Of the 40 sources in `existing_sources.txt`, 28 were retrieved.

Database sources not retrieved by any query (13), by source ID and DOI (the Watanabe 2008 record, DOI 10.1143/jjap.47.2841, is matched by Q1 and Q6):

- McLevige1977_JAP, 10.1063/1.324218
- Chambon1985_APL, 10.1063/1.95722
- Maezawa1986_EDL, 10.1109/edl.1986.26276
- Anderson1979_APL, 10.1063/1.91027
- Opyd1990_JAP, 10.1063/1.344531
- Donnelly1979_APL, 10.1063/1.90573
- Dinia2000_SSC, 10.1016/s0038-1098(99)00511-6
- Shen1991_NIMB, 10.1016/0168-583x(91)96283-q
- Powell1984_JAP, 10.1063/1.333818
- Narayan1983_APL, 10.1063/1.94200
- Saavedra2004_JAP, 10.1063/1.1769095
- Itoh1998_APL, 10.1063/1.121965 (a hit "Comparison of the annealing behavior of high-dose nitrogen-, aluminum-, and boron-implanted 4H-SiC", 1998, has DOI 10.1063/1.121681; not verified)
- Flemish1995_JES, 10.1149/1.2048727

### What this recall does and does not show

- The database was built from earlier semantic literature searches and web searches plus manually downloaded open-access PDFs. Recall against it is therefore partly circular: it measures overlap with a set that was itself assembled partly by similar searches, not coverage of the true literature.
- Of the 13 unretrieved sources, 11 are GaAs (5), InP (3) or Si (3) papers from 1977 to 2004, systems the six queries (SiC, GaN and diamond focused) were not designed to target; the other 2 are SiC (Itoh1998_APL, Flemish1995_JES). The check therefore says little about recall for GaAs, InP and Si.
- Matching is strict on DOI. Alternate DOIs for the same work (preprint vs. journal, proceedings vs. journal) would count as misses, so recall may be slightly understated.
- The check does not show that papers outside the database and outside the search hits do not exist. A single semantic search with the top 100 results per query is a limited sample.

## Screening of search hits not in the database (title and abstract only)

370 unique hits are not in the database. Screening used only the title and abstract text returned by the tool; no full texts were fetched or read. Classification was made by a keyword rule on title and abstract (implanted acceptor species, host material, electrical measure, annealing), followed by manual review of the topic-eligible and borderline hits. It was not double-screened, and some abstracts are truncated or missing, so labels are provisional. `topic_ok` = ion-implanted acceptor activation measured electrically in SiC, GaN, diamond, GaAs, InP or Si. Per-hit labels are in `screening.csv`.

| topic_ok | Screening label | Count |
|---|---|---|
| yes | abstract suggests numeric activation vs process conditions (potentially eligible, not extracted) | 99 |
| yes | unclear | 45 |
| yes | review or theory/simulation only | 12 |
| no | not an implanted-acceptor activation study | 171 |
| no | n-type or donor | 37 |
| no | review or theory/simulation only | 6 |

Notes:
- "not an implanted-acceptor activation study" includes defect, positron, photoluminescence, TEM, diffusion, contact, MOS-interface and device papers without an electrical activation measurement in the abstract, and hosts outside the list (for example Ga2O3, AlN, Ge, h-BN).
- "Potentially eligible" means only that the abstract mentions activation or electrical results together with annealing conditions. It does not mean the paper reports usable numbers (implant temperature, anneal conditions, activation or hole concentration). That requires full-text review.
- Papers excluded for inaccessible full text and papers excluded for incomplete numerical reporting are NOT yet separated at this stage: no full texts were attempted for these hits. That distinction applies only to the database-extraction stage, using `reporting_audit.csv` and the full-text access record for sources actually read.

## Potentially eligible, not yet extracted (99 papers)

| DOI | Year | Title |
|---|---|---|
| 10.1080/00207216908938149 | 1969 | The measurement of electrical activity and Hall mobility of boron and phosphorus ion-implanted layers in silicon† |
| 10.1080/00337577008236299 | 1970 | Low energy boron and phosphorus implants in silicon (b) doping profiles |
| 10.1080/00337577108232558 | 1971 | The isothermal annealing of boron implanted silicon |
| 10.1007/bf00618767 | 1988 | Electrical activation of low-fluence boron implantation in silicon studied by PCV in combination with SIMS |
| 10.1063/1.101605 | 1989 | Doping of diamond by coimplantation of carbon and boron |
| 10.1063/1.350967 | 1992 | Hall mobility and carrier concentration versus temperature for type IIa natural insulating diamond doped with boron by ion implantation |
| 10.1016/0925-9635(92)90192-q | 1992 | Issues in diamond device fabrication |
| 10.1016/0921-5107(92)90223-v | 1992 | p-Type semiconducting structures in diamond implanted with boron ions |
| 10.1103/physrevb.47.2065 | 1993 | Carrier activation and mobility of boron-dopant atoms in ion-implanted diamond as a function of implantation conditions. |
| 10.1063/1.355098 | 1993 | Electrical activation of boron coimplanted with carbon in a silicon substrate |
| 10.1016/0925-9635(94)90237-2 | 1994 | Reality of doping by boron implantation of CVD polycrystalline diamond from a comparison of Raman and electrical measurements |
| 10.1063/1.358776 | 1995 | Al and B ion‐implantations in 6H‐ and 3C‐SiC |
| 10.1063/1.113713 | 1995 | Metastable acceptor centers in boron implanted silicon |
| 10.1557/proc-416-249 | 1995 | Temperature Sensor On Boron Ion Implanted Diamond |
| 10.1007/bf02666177 | 1996 | Al, Al/C and Al/Si implantations in 6H-SiC |
| 10.1002/1521-396x(199707)162:1<277::aid-pssa277>3.0.co;2-c | 1997 | Doping of SiC by Implantation of Boron and Aluminum |
| 10.1063/1.119946 | 1997 | Formation of buried p-type conducting layers in diamond |
| 10.1063/1.120173 | 1997 | Observation of boron bound excitons in boron-implanted and annealed natural IIa diamonds |
| 10.1063/1.365236 | 1997 | Phosphorus and boron implantation in 6H–SiC |
| 10.1063/1.121681 | 1998 | Comparison of the annealing behavior of high-dose nitrogen-, aluminum-, and boron-implanted 4H–SiC |
| 10.1007/s11664-998-0417-4 | 1998 | Dopant activation and surface morphology of ion implanted 4H- and 6H-silicon carbide |
| 10.1016/s0925-9635(98)00231-3 | 1998 | Formation of delta-doped, buried conducting layers in diamond, by high-energy, B-ion implantation |
| 10.1063/1.123429 | 1999 | Efficient p-type doping of 6H-SiC: Flash-lamp annealing after aluminum implantation |
| 10.1016/s0169-4332(99)00153-1 | 1999 | The effect of boron ion implantation and annealing on the microstructure and electrical characteristics of the diamond films |
| 10.4028/www.scientific.net/msf.353-356.571 | 2001 | High Electrical Activation of Aluminium and Nitrogen Implanted in 6H-SiC at Room Temperature by RF Annealing |
| 10.4028/www.scientific.net/msf.389-393.827 | 2002 | A Comparative Study of High-Temperature Aluminum Post-Implantation Annealing in 6H- and 4H-SiC, Non-Uniform Temperature Effects |
| 10.1063/1.1490145 | 2002 | Low-resistivity, p-type SiC layers produced by Al implantation and ion-beam-induced crystallization |
| 10.1063/1.1598631 | 2003 | Effect of ion implantation parameters on Al dopant redistribution in SiC after annealing: Defect recovery and electrical properties of p-type layers |
| 10.4028/www.scientific.net/msf.433-436.605 | 2003 | Electrical Activation of Ion-Implanted Nitrogen and Aluminum in 4H-SiC by Excimer Laser Annealing |
| 10.4028/www.scientific.net/msf.433-436.637 | 2003 | Gettering Effect with Al Implanted into 4H-SiC CVD Epitaxial Layers |
| 10.1002/pssa.200306704 | 2004 | Electrical, CL, EPR and RBS study of annealed SiC implanted with Al or Al and C |
| 10.4028/www.scientific.net/msf.457-460.913 | 2004 | Low Sheet Resistance of High-Dose Aluminium Implanted 4H-SiC using (11-20) Face |
| 10.4028/www.scientific.net/msf.483-485.621 | 2005 | Annealing of Aluminum Implanted 4H-SiC: Comparison of Furnace and Lamp Annealing |
| 10.1063/1.2030411 | 2005 | Carrier compensation near tail region in aluminum- or boron-implanted 4H-SiC(0001) |
| 10.4028/www.scientific.net/msf.483-485.617 | 2005 | Electrical Behavior of Implanted Aluminum and Boron near Tail Region in 4H-SiC after High-Temperature Annealing |
| 10.4028/www.scientific.net/msf.483-485.401 | 2005 | Electrical Properties of p-Type In-Situ Doped vs. Al-Implanted 4H-SiC |
| 10.1016/j.nimb.2005.04.079 | 2005 | Implantation and annealing of aluminum in 4H silicon carbide |
| 10.1109/isdrs.2005.1596103 | 2005 | The Effects of Implanting Al and Al and C at Different Temperatures in Different Concentrations into SiC |
| 10.1063/1.2196233 | 2006 | Effects of heating ramp rates on the characteristics of Al implanted 4H–SiC junctions |
| 10.4028/www.scientific.net/msf.527-529.831 | 2006 | Variations in the Effects of Implanting Al at Different Concentrations into SiC |
| 10.1063/1.2813022 | 2007 | Acceptor, compensation, and mobility profiles in multiple Al implanted 4H‐SiC |
| 10.1007/s11664-006-0032-1 | 2007 | Comparison of Solid-State Microwave Annealing with Conventional Furnace Annealing of Ion-Implanted SiC |
| 10.4028/www.scientific.net/msf.556-557.343 | 2007 | Electrical and Structural Properties of Al-Implanted and Annealed 4H-SiC |
| 10.1063/1.2715034 | 2007 | High-pressure and high-temperature annealing as an activation method for ion-implanted dopants in diamond |
| 10.1016/j.sse.2007.06.021 | 2007 | Ultra-low resistivity Al+ implanted 4H–SiC obtained by microwave annealing and a protective graphite cap |
| 10.4028/www.scientific.net/msf.600-603.603 | 2008 | Annealing Temperature Dependence of the Electrically Active Profiles and Surface Roughness in Multiple Al Implanted 4H-SiC |
| 10.4028/www.scientific.net/msf.600-603.585 | 2008 | Crystalline Recovery after Activation Annealing of Al Implanted 4H-SiC |
| 10.1109/edssc.2009.5394252 | 2009 | Annealing effects on structural, optical and electrical properties of Al implanted 4H-SiC |
| 10.4028/www.scientific.net/ssp.156-158.493 | 2009 | Correlation Study of Morphology, Electrical Activation and Contact formation of Ion Implanted 4H-SiC |
| 10.4028/www.scientific.net/msf.615-617.457 | 2009 | Electrical Characterization of Al Implanted 4H-SiC Layers by Four Point Probe and Scanning Capacitance Microscopy |
| 10.1016/j.diamond.2008.09.013 | 2009 | Electrical properties of a B doped layer in diamond formed by hot B implantation and high-temperature annealing |
| 10.4028/www.scientific.net/msf.645-648.705 | 2010 | Effects of Implantation Temperature on Sheet and Contact Resistance of Heavily Al Implanted 4H-SiC |
| 10.4028/www.scientific.net/msf.645-648.697 | 2010 | Electrical Activation of B+-Ions Implanted into 4H-SiC |
| 10.4028/www.scientific.net/msf.645-648.717 | 2010 | Influence of Heating and Cooling Rates of Post-Implantation Annealing Process on Al-Implanted 4H-SiC Epitaxial Samples |
| 10.1016/j.matlet.2010.08.048 | 2010 | Study of SiO2 encapsulation for aluminum and phosphorus implant activation in 4H-SiC |
| 10.1016/j.nimb.2011.11.025 | 2012 | Change of structural and electrical properties of diamond with high-dose ion implantation at elevated temperatures: Dependences on donor/acceptor impurity species |
| 10.4028/www.scientific.net/msf.717-720.817 | 2012 | High Dose Al+ Implanted and Microwave Annealed 4H-SiC |
| 10.1007/s11431-012-4827-4 | 2012 | Investigation of surface morphology and ion activation of aluminium implanted 4H-SiC |
| 10.4028/www.scientific.net/msf.740-742.767 | 2013 | Al+ Implanted 4H-SiC: Improved Electrical Activation and Ohmic Contacts |
| 10.4028/www.scientific.net/msf.740-742.581 | 2013 | Aluminum Implantation in 4H-SiC: Physical and Electrical Properties |
| 10.14723/tmrsj.38.431 | 2013 | Optical absorptions, electrical properties and related band structures of diamond with high-dose B, C, N or P implantation at elevated temperatures |
| 10.4028/www.scientific.net/msf.778-780.673 | 2014 | Development of a Novel Cap-Free Activation Annealing Technique of 4H-SiC by Si-Vapor Ambient Annealing Using TaC/Ta Composite Materials |
| 10.4028/www.scientific.net/msf.778-780.653 | 2014 | Microwave Annealing of Al+ Implanted 4H-SiC: Towards Device Fabrication |
| 10.4028/www.scientific.net/msf.821-823.416 | 2015 | About the Hole Transport Analysis in Heavy Doped p-Type 4H-SiC(Al) |
| 10.4028/www.scientific.net/msf.821-823.399 | 2015 | Al+ Ion Implanted On-Axis <0001> Semi-Insulating 4H-SiC |
| 10.4028/www.scientific.net/msf.821-823.407 | 2015 | Distribution of Secondary Defects and Electrical Activation after Annealing of Al-Implanted SiC |
| 10.1149/07512.0171ecst | 2016 | (Invited) 4H-SiC Ion Implanted Bipolar Junctions: Relevance of the 1950°C Temperature for Post Implantation Annealing |
| 10.4028/www.scientific.net/msf.858.523 | 2016 | 1950°C Annealing of Al+ Implanted 4H-SiC: Sheet Resistance Dependence on the Annealing Time |
| 10.1109/iwjt.2016.7486673 | 2016 | Ion implantation technology in SiC for high-voltage/high-temperature devices |
| 10.7567/jjap.56.070306 | 2017 | Electrical properties of n- and p-type 4H-SiC formed by ion implantation into high-purity semi-insulating substrates |
| 10.1002/pssa.201700230 | 2017 | Impact of Thermal Treatments in Crystalline Reconstruction and Electrical Properties of Diamond Ohmic Contacts Created by Boron Ion Implantation |
| 10.4028/www.scientific.net/msf.924.333 | 2018 | About the Electrical Activation of 1×1020 cm-3 Ion Implanted Al in 4H-SiC at Annealing Temperatures in the Range 1500 - 1950°C |
| 10.1149/2.0121909jss | 2019 | 1300°C Annealing of 1 × 1020 cm−3 Al+ Ion Implanted 3C-SiC/Si |
| 10.4028/www.scientific.net/msf.963.420 | 2019 | 1300°C Annealing of 1×1020 Al+ Ion Implanted 3C-SiC |
| 10.4028/www.scientific.net/msf.1004.683 | 2019 | 3× 1018 - 1 × 1019 cm -3 Al+ Ion Implanted 4H-SiC: Annealing Time Effect |
| 10.4028/www.scientific.net/msf.963.416 | 2019 | Activation Energy for the Post Implantation Annealing of 1019 cm-3 and 1020 cm-3 Ion Implanted Al in 4H SiC |
| 10.1088/1674-1056/ab3cc2 | 2019 | Defects and electrical properties in Al-implanted 4H-SiC after activation annealing |
| 10.1016/j.mssp.2019.01.019 | 2019 | Effect of high temperature annealing (T > 1650 °C) on the morphological and electrical properties of p-type implanted 4H-SiC layers |
| 10.3390/ma12203362 | 2019 | Laser Annealing of P and Al Implanted 4H-SiC Epitaxial Layers |
| 10.1002/pip.3186 | 2019 | Plasma‐immersion ion implantation: A path to lower the annealing temperature of implanted boron emitters and simplify PERT solar cell processing |
| 10.35848/1347-4065/ab699c | 2020 | Electrical properties and conduction mechanisms of heavily B+-ion-implanted type IIa diamond: effects of temperatures during the ion implantation and postannealing upon electrical conduction |
| 10.35848/1882-0786/abd308 | 2020 | Enhanced activation of Mg ion-implanted GaN at decreasing annealing temperature by prolonging duration |
| 10.35848/1347-4065/aba9a6 | 2020 | Initial activation behavior of boron at low temperatures with implantation doses below the amorphization threshold |
| n/a | 2021 | A comparative study of high temperature Aluminium post-implantation annealing in 6H and 4H-SiC, non-uniformity temperature effects |
| 10.35848/1882-0786/ac39b0 | 2021 | Effect of annealing time and pressure on electrical activation and surface morphology of Mg-implanted GaN annealed at 1300 °C in ultra-high-pressure nitrogen ambient |
| 10.1063/5.0049101 | 2021 | P-type conductivity and suppression of green luminescence in Mg/N co-implanted GaN by gyrotron microwave annealing |
| 10.4028/p-sg9dq0 | 2022 | Aluminum Activation in 4H-SiC Measured on Laterally Contacted MOS Capacitors with a Buried Current-Spreading Layer |
| 10.4028/p-98r1yh | 2022 | Graphite Assisted P and Al Implanted 4H-SiC Laser Annealing |
| 10.1557/s43580-022-00403-2 | 2022 | Ion implantation and activation of aluminum in bulk 3C-SiC and 3C-SiC on Si |
| 10.1063/5.0107921 | 2022 | Process engineering of GaN power devices via selective-area p-type doping with ion implantation and ultra-high-pressure annealing |
| 10.61640/ujeas.2023.1211 | 2023 | Growth of synthetic diamond films and their electrophysical properties |
| 10.1016/j.mssp.2024.108175 | 2024 | Exploring crystal recovery and dopant activation in coated laser annealing on ion implanted 4H–SiC epitaxial layers |
| 10.4028/p-q1ovxc | 2024 | Free-Standing 3C-SiC P-Type Doping by Al Ion Implantation |
| 10.1002/pssb.202400613 | 2025 | Depth Analyses of Acceptor and Compensating Donor Concentrations in Mg‐Ion‐Implanted GaN with Ultra‐High‐Pressure Annealing |
| 10.1002/pssr.202500346 | 2025 | P‐Type Electrical Activation in GaN by Mg Ion Implantation for PN Junctions |
| 10.1063/5.0272310 | 2025 | The origin of degradation in electrical properties by ion implantation doping of boron in diamond |
| 10.1063/5.0341213 | 2026 | Influence of implantation and annealing temperatures on electrical property of low-dose implanted/annealed SiC |
| 10.1016/j.mlblux.2026.100266 | 2026 | Influence of the Al-ion implantation temperature in 4H-SiC on the activation rate and the quality of the subsequent thermal oxide |
| 10.1016/j.mssp.2026.110502 | 2026 | Transient and long-term effects of Al implantation temperature on electrical and structural properties of 4H-SiC |

## Files

- `search_hits.csv`: one row per (query, hit), with rank, DOI, title, first author, authors field, year, journal, citation count and abstract
- `search_hits_unique.csv`: one row per DOI, with the queries that returned it
- `screening.csv`: label for each of the 370 unique hits not in the database
- Recall computed with python3 from `code_release/data/database.csv`
