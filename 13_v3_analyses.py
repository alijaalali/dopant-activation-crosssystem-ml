import json

import pandas as pd

from common import DATA
from models import FEATS, fit_hier, load_sic_al

al = load_sic_al(DATA / "database.csv")
db = pd.read_csv(DATA / "database.csv")
OUT = {}


def brief(r):
    return {"n_rows": r["n_rows"], "n_sources": r["n_sources"], "n_labs": r["n_labs"],
            "implant_T_pp": r["marginal_effect_pp_per_sd"]["implant_temp_C"], "anneal_T_pp": r["marginal_effect_pp_per_sd"]["anneal_temp_C"],
            "sd": r["sd"], "share": r["share"]}


assumed_depth_sources = ["S3"]
subsets = {
    "all_modelled_rows": al,
    "without_assumed_implanted_depth": al[~al.source_id.isin(assumed_depth_sources)],
    "stated_values_only": al[al.value_basis == "stated"],
    "stated_values_without_resistivity_referenced": al[(al.value_basis == "stated") & (al.definition_class != "resistivity_adjusted")],
    "without_assumed_depth_and_resistivity_referenced": al[(~al.source_id.isin(assumed_depth_sources)) & (al.definition_class != "resistivity_adjusted")],
}
res = {}
for name, frame in subsets.items():
    frame = frame.reset_index(drop=True)
    entry = {"rows": int(len(frame)), "sources": int(frame.source_id.nunique()), "labs": int(frame.lab_group.nunique())}
    for fac in ["source", "both"]:
        entry[fac] = brief(fit_hier(frame, factors=fac))
    res[name] = entry
    print(name, entry["rows"], entry["sources"], entry["source"]["implant_T_pp"]["mean"], entry["source"]["share"]["between_groups_total"]["mean"], flush=True)
OUT["subset_refits"] = res

lab = al.groupby("lab_group").source_id.nunique()
OUT["sources_per_laboratory"] = {k: int(v) for k, v in lab.items()}
OUT["rows_per_source"] = {k: int(v) for k, v in al.groupby("source_id").size().items()}
OUT["assumed_depth_rows"] = {"sources": assumed_depth_sources, "rows": int(al.source_id.isin(assumed_depth_sources).sum()),
                             "reason": "implanted concentration derived from the dose with an assumed 0.2 um box depth"}

sic = db[db.system == "SiC"]
OUT["sic_rows_by_dopant"] = {k: int(v) for k, v in sic.dopant.value_counts().items()}
OUT["sic_rows_by_dopant_source"] = {f"{a}|{b}": int(v) for (a, b), v in sic.groupby(["dopant", "source_id"]).size().items()}
OUT["classification_check"] = {"sic_rows": int(len(sic)), "al_rows": int((sic.dopant == "Al").sum()), "non_al_rows": int((sic.dopant != "Al").sum())}

kv = json.load(open(DATA / "kinetic_validation_results.json"))
OUT["temporal_test_sources"] = kv["temporal"]["test_sources"]
OUT["temporal_train_rows_sources"] = {k: v for k, v in kv["temporal"].items() if k in ("n_train", "n_test", "train_sources")}
json.dump(OUT, open(DATA / "revision_v3_results.json", "w"), indent=2)
print("DONE")
