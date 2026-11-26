"""Turn Wikidata zh-label hits into extra NMPA rows (source 'NM' for the merge script)."""
import pandas as pd
from pathlib import Path
HERE = Path(__file__).resolve().parent
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog("rdApp.*")
lfc, unch = rdMolStandardize.LargestFragmentChooser(), rdMolStandardize.Uncharger()
c = pd.read_csv(HERE / "nmpa_candidates_clean.csv")
c = c[c.status != "ok"]                       # only records DeepSeek could not resolve
h = pd.read_csv(HERE / "wikidata_zh_hits.csv", dtype=str)
import gzip, json
_ch = json.load(gzip.open(HERE.parents[1] / "data" / "raw" / "chembl_molecules.json.gz"))
_ik2smi = {c["inchikey"]: c["smiles"] for c in _ch if c.get("inchikey") and c.get("smiles")}
h["smiles"] = h.ismi.fillna(h.smi).fillna(h.ik.map(_ik2smi))
h = h[h.smiles.notna()]

def parent(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    p = unch.uncharge(lfc.choose(m))
    return p
recs = {}
for lab, g in h.groupby("lab"):
    ps = []
    for r in g.itertuples():
        p = parent(r.smiles)
        if p is not None: ps.append((r, p, Chem.MolToInchiKey(p)))
    if not ps: continue
    if len({ik[:14] for _, _, ik in ps}) > 1:   # ambiguous Chinese name -> skip
        continue
    r, p, ik = ps[0]
    recs[lab] = dict(label=r.itemLabel, smiles=Chem.MolToSmiles(p), ik=ik, qid=r.item.rsplit("/", 1)[-1],
                     mw=Descriptors.MolWt(p), nC=sum(a.GetSymbol() == "C" for a in p.GetAtoms()), ha=p.GetNumHeavyAtoms())
rows = []
for r in c.itertuples():
    hit = recs.get(r.zh_clean) or recs.get(r.zh_base)
    if not hit: continue
    if hit["mw"] > 2000 or hit["nC"] == 0 or hit["ha"] < 3: continue
    rows.append(dict(zh=r.zh, zh_match=r.zh_clean if r.zh_clean in recs else r.zh_base, year=int(r.year), **hit))
o = pd.DataFrame(rows)
agg = o.groupby("ik").agg(Name=("label", "first"), SMILES=("smiles", "first"), NMPA=("year", "min"),
                          zh=("zh", lambda s: "; ".join(sorted(set(s))[:5])), qid=("qid", "first")).reset_index()
agg["Name"] = agg.Name.str[:1].str.upper() + agg.Name.str[1:]
agg["FDA"] = agg["EMA"] = agg["PMDA"] = ""
agg["NMPA"] = agg.NMPA.astype(str)
agg[["Name", "SMILES", "FDA", "EMA", "NMPA", "PMDA", "zh", "qid"]].to_csv(HERE / "nmpa_rescued_Claude.csv", index=False)
print(f"unresolved NMPA records: {len(c)} | matched to a structure: {len(o)} | unique molecules: {len(agg)}")
print(agg[["Name", "NMPA", "zh"]].sample(min(25, len(agg)), random_state=0).to_string())
