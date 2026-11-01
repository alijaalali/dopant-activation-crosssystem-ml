import json

import numpy as np
import pandas as pd

from common import DATA

db = pd.read_csv(DATA / "database.csv")
OUT = {}

m = db[(db.source_id == "Nipoti2022_MSF") & db.activation_pct.notna()].copy()
m["sample"] = m.location_in_paper.astype(str) + m.implanted_conc_cm3.astype(str) + m.anneal_temp_C.astype(str) + m.anneal_time_min.astype(str)
m["sample"] = m.notes.str.extract(r"(\d{3}[a-z]?-[A-Za-z0-9-]+)", expand=False).fillna(m["sample"])
wide = m.pivot_table(index=["implanted_conc_cm3", "anneal_temp_C", "anneal_time_min", "sample"], columns="method_variant",
                     values="activation_pct", aggfunc="first").dropna()
d = (wide["p"] - wide["p-mu"]).to_numpy()
OUT["same_samples_two_fitting_methods"] = {
    "n_samples": int(len(d)), "mean_difference_pp": float(d.mean()), "min_pp": float(d.min()), "max_pp": float(d.max()),
    "sd_pp": float(d.std(ddof=1)), "p_gt_pmu_in_samples": int((d > 0).sum())}
print(OUT["same_samples_two_fitting_methods"])

OUT["same_nominal_conditions_two_papers"] = [
    {"conc_cm3": 1e19, "anneal": "1600 C, 30 min", "ECS_Trans_92": 40, "MSF_1062_p": 60, "MSF_1062_p-mu": 31},
    {"conc_cm3": 3e18, "anneal": "1600 C, 30 min", "ECS_Trans_92": 81, "MSF_1062_p": 96, "MSF_1062_p-mu": 77}]

a = pd.read_csv(DATA / "reporting_audit.csv")
hall = a[~a.hall_factor_stated.astype(str).str.startswith("not applicable") & (a.n_conditions.astype(str) != "0")]
yes = lambda s: s.astype(str).str.lower().str.startswith("yes").sum()
OUT["audit"] = {
    "n_papers_read_in_full": int(len(a)),
    "reports_implant_temperature": int(yes(a.reports_implant_temp)),
    "reports_anneal_time": int(yes(a.reports_anneal_time)),
    "temperature_dependent_hall": int(yes(a.temperature_dependent_hall)),
    "ionization_correction": int(yes(a.ionization_correction)),
    "compensation_reported": int(yes(a.compensation_reported)),
    "activation_defined_in_text": int(yes(a.activation_defined_in_text)),
    "data_availability_statement": int(yes(a.data_availability_statement)),
    "hall_factor_value_or_model_stated": int(a.hall_factor_stated.astype(str).str.lower().str.match(r"^(value|1|r_exp|yes|partly)").sum()),
    "hall_factor_not_stated": int(a.hall_factor_stated.astype(str).str.lower().str.startswith("not stated").sum()),
}
print(OUT["audit"])

counts = db[db.system == "SiC"].groupby("source_id").size()
OUT["database"] = {"rows": int(len(db)), "sources": int(db.source_id.nunique()), "sic_rows": int((db.system == "SiC").sum()),
                   "sic_sources": int(db[db.system == "SiC"].source_id.nunique())}
json.dump(OUT, open(DATA / "method_audit_results.json", "w"), indent=2)
