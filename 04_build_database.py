import re

import numpy as np
import pandas as pd

from common import DATA

df = pd.read_csv(DATA / "activation_dataset.csv", na_values=["NA"])

LAB = {"S1_S2": "CNR-IMM-Catania", "Nipoti2011": "CNR-IMM-Bologna", "Nipoti2013_JMR": "CNR-IMM-Bologna", "Weng2010": "CNR-IMM-Catania",
       "Saks2004": "NRL-Cree",
       "S4_MSF": "Fraunhofer-IISB", "Weisse2019_AIPAdv": "Fraunhofer-IISB",
       "S3": "Kyoto", "Watanabe2008_shuttle": "Kyoto", "Niwa2017_APEX": "Kyoto",
       "G1": "Nagoya", "G2_Hirukawa": "Nagoya", "Sumida2025_Ncoimp": "Nagoya", "Sumida2025_Ccap": "Nagoya",
       "Sakurai2020_ISPSD": "Nagoya", "G3": "NRL", "G4": "NRL", "Breckenridge2021": "NCSU"}
YEAR = {"S1_S2": 2019, "S3": 2004, "S4_MSF": 2018, "G1": 2019, "G2_Hirukawa": 2021, "G3": 2014, "G4": 2012}


def year_of(sid):
    if sid in YEAR:
        return YEAR[sid]
    m = re.search(r"(19|20)\d\d", sid)
    return int(m.group(0)) if m else np.nan


def definition_class(row):
    vb = row.get("value_basis", "stated")
    sid = row["source_id"]
    if vb == "net_of_compensation":
        return "net_of_compensation"
    if sid == "Saks2004":
        return "resistivity_adjusted"
    if vb == "derived":
        return "derived_ratio"
    if vb == "approximate":
        return "approximate"
    return "stated_hall"


base = df.copy()
base["origin"] = "initial_compilation"
base["lab_group"] = base.source_id.map(LAB)
base["year"] = base.source_id.map(year_of)
base["definition_class"] = base.apply(definition_class, axis=1)
for c in ["dose_cm-2", "hole_conc_cm-3", "mobility_cm2Vs", "EA_meV"]:
    base[c] = pd.to_numeric(base[c], errors="coerce")
base = base.rename(columns={"dose_cm-2": "dose_cm2", "hole_conc_cm-3": "hole_conc_cm3"})
base["implanted_conc_cm3"] = np.nan
base["sheet_hole_cm2"] = np.nan
cols = ["origin", "system", "host", "dopant", "source_id", "lab_group", "year", "doi_or_arxiv", "implant_temp_C",
        "anneal_temp_C", "anneal_time_min", "anneal_method", "cap", "co_implant", "measurement_technique",
        "dose_cm2", "implanted_conc_cm3", "hole_conc_cm3", "sheet_hole_cm2", "mobility_cm2Vs", "EA_meV",
        "activation_pct", "value_basis", "text_basis", "definition_class", "hall_factor_convention"]
base.loc[base.source_id == "S1_S2", "implant_temp_C"] = 500
for c in ["implant_temp_C", "anneal_temp_C", "anneal_time_min", "activation_pct"]:
    base[c] = pd.to_numeric(base[c], errors="coerce")
base["two_step_anneal"] = False
base["implant_T_is_assumed"] = base.source_id == "S1_S2"
CONC = {"Saks2004": 1.8e20, "S3": 1.5e21, "S1_S2": 1.0e20, "Nipoti2011": 8e20, "Nipoti2013_JMR": 1e19}
for sid, c in CONC.items():
    base.loc[base.source_id == sid, "implanted_conc_cm3"] = c
base.loc[(base.source_id == "Weisse2019_AIPAdv") & (base.activation_pct == 61), "implanted_conc_cm3"] = 1e17
base.loc[(base.source_id == "Weisse2019_AIPAdv") & (base.activation_pct == 67), "implanted_conc_cm3"] = 1e18
cols += ["implant_T_is_assumed"]
out = base[[c for c in cols if c in base.columns]].copy()

for c in ["two_step_anneal", "implant_T_is_range", "compensation_ratio", "sheet_resistance_ohm_sq", "ionization_correction",
          "overlap_note", "method_variant", "location_in_paper", "notes"]:
    if c not in out.columns:
        out[c] = np.nan
try:
    ex = pd.read_csv(DATA / "expansion_verified.csv", na_values=["NA"])
except FileNotFoundError:
    ex = None
if ex is not None:
    ex["origin"] = "expansion"
    out = pd.concat([out, ex[[c for c in out.columns if c in ex.columns]]], ignore_index=True)

for c in ["two_step_anneal", "implant_T_is_assumed", "implant_T_is_range"]:
    out[c] = out[c].fillna(False).astype(bool)
out = out.reset_index(drop=True)
out.insert(0, "db_id", np.arange(1, len(out) + 1))
out.to_csv(DATA / "database.csv", index=False)
print(len(out), out.source_id.nunique(), out.origin.value_counts().to_dict())
