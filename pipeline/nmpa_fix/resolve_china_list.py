"""Attach structures to the web-verified China innovative-drug list and append them to
pipeline/curation/verified_additions_Claude.csv (NMPA column only); also adds their alternative
names to pipeline/curation/name_aliases_Claude.csv."""
import csv, json, sys, time, urllib.parse, urllib.request
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "round1_agents" / "claude"))
import collect_drugs as cd  # noqa: E402

here = Path(__file__).parent
lst = pd.read_csv(here / "china_innovative_verified_Claude.csv", dtype=str)
parents, child_of, index = cd.load_chembl()


def pubchem(name):
    try:
        u = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + urllib.parse.quote(name)
             + "/property/IsomericSMILES,SMILES/JSON")
        p = json.load(urllib.request.urlopen(u, timeout=30))["PropertyTable"]["Properties"][0]
        return p.get("IsomericSMILES") or p.get("SMILES")
    except Exception:
        return None
    finally:
        time.sleep(0.25)


rows, unresolved = [], []
for r in lst.itertuples():
    if not r.NMPA_year.isdigit():
        continue
    names = [n.strip() for n in r.INN.split("/")]
    smi, via = None, ""
    for n in names:
        cid = index.lookup(n)
        if cid:
            smi, via = parents[cid]["smiles"], f"ChEMBL {cid}"
            break
    if not smi:
        for n in names:
            smi = pubchem(n)
            if smi:
                via = f"PubChem name '{n}'"
                break
    if not smi:
        unresolved.append(r.INN)
        continue
    nm = names[0]
    rows.append({"Name": nm[:1].upper() + nm[1:], "SMILES": smi, "NMPA": r.NMPA_year,
                 "source": f"NMPA {r.NMPA_year} ({r.Chinese_name}; {r.confidence}; {r.source}); structure: {via}"})

va = ROOT / "pipeline" / "curation" / "verified_additions_Claude.csv"
old = pd.read_csv(va, dtype=str, keep_default_na=False)
old = old[~old.source.str.match(r"NMPA \d{4} \(")]     # keep hand-made rows, refresh this list's rows
new = pd.DataFrame(rows).reindex(columns=old.columns).fillna("")
pd.concat([old, new]).to_csv(va, index=False)
alias = ROOT / "pipeline" / "curation" / "name_aliases_Claude.csv"
al = pd.read_csv(alias, dtype=str)
extra = [(a.strip(), b.strip()) for inn in lst.INN for a, b in
         [(inn.split("/")[0], x) for x in inn.split("/")[1:]]]
al = pd.concat([al, pd.DataFrame(extra, columns=al.columns)]).drop_duplicates()
al.to_csv(alias, index=False)
print(f"approved entries: {sum(lst.NMPA_year.str.isdigit())} | with structure: {len(rows)} | unresolved: {len(unresolved)}")
print("unresolved:", unresolved)
