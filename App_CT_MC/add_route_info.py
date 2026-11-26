"""Step 3: route of administration (oral vs non-oral) for the macrocycles only.

For each macrocycle in ApprovedDrug_CT_2026Nov_MC_only.csv, collect route evidence from:
  1. FDA      Drugs@FDA products (openFDA bulk file, field `route`)            data/raw/openfda_drugsfda.json.zip
  2. EU       EMA Article 57 database (route of every EU/EEA product)          data/raw/ema_article57.xlsx
  3. France   ANSM BDPM (voies d'administration)                             pipeline/eu_national/raw/
  4. China    dosage form in NMPA product names (片/胶囊 = oral, 注射 = injection) data/nmpa_candidates.csv,
              pipeline/nmpa_fix/*, pipeline/curation/verified_additions_Claude.csv
  5. Trials   ClinicalTrials.gov intervention/arm descriptions (cached in route_cache/ctgov_routes.json)
  6. ChEMBL   oral / parenteral / topical flags of approved drugs (cached in route_cache/chembl_flags.json)
Marketed products (1-4) decide the class when present; otherwise trials (5); otherwise ChEMBL flags (6).

Adds three columns to ApprovedDrug_CT_2026Nov_MC_only.csv and ApprovedDrug_CT_2026Nov_MC_MolDs_Onion.csv:
  Route_Class   Oral | Injection | Oral + Injection | Non-oral (other) | Unknown
                (Non-oral (other) = only topical / ophthalmic / inhaled / rectal / vaginal ... products)
  Routes        normalised routes found, e.g. "oral; intravenous; topical"
  Route_Source  which evidence decided it (and a per-source summary in ApprovedDrug_CT_2026Nov_MC_routes.csv)
"""
import gzip
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
os.chdir(HERE)
sys.path.insert(0, str(ROOT / "round1_agents" / "claude"))
sys.path.insert(0, str(ROOT / "pipeline" / "common"))
import collect_drugs as cd  # noqa: E402  (ChEMBL name index + cached downloads)
from name_mapper import map_names  # noqa: E402

CACHE = HERE / "route_cache"
CACHE.mkdir(exist_ok=True)
MC_FILE = HERE / "ApprovedDrug_CT_2026Nov_MC_only.csv"
MOLD_FILE = HERE / "ApprovedDrug_CT_2026Nov_MC_MolDs_Onion.csv"
OUT_DETAIL = HERE / "ApprovedDrug_CT_2026Nov_MC_routes.csv"
AG = ["FDA", "EMA", "EU_National", "NMPA", "PMDA", "ANVISA"]

# ----------------------------------------------------------------------------- route vocabulary
INJECTION = "injection"
ROUTE_RULES = [  # (regex on lower-case text, normalised route, category)
    (r"\boral|per ?os\b|orale|tablet|capsule|by mouth|orally", "oral", "oral"),
    (r"intraven|\biv\b|i\.v\.|infusion|intraveineuse", "intravenous", INJECTION),
    (r"intramusc|\bim\b|i\.m\.", "intramuscular", INJECTION),
    (r"subcutan|sous-cutan|\bsc\b|s\.c\.", "subcutaneous", INJECTION),
    (r"intravitreal|intravitr", "intravitreal", INJECTION),
    (r"intrathecal|intravesical|intra-?arterial|intralesional|intraperitoneal|intrapleural|"
     r"intracavern|intradermal|intra-?articular|epidural|perineural", "other injection", INJECTION),
    (r"injection|injectable|parenteral|injectable", "injection", INJECTION),
    (r"ophthalm|ocular|\beye\b|oculaire|ophtalmique|conjunctival", "ophthalmic", "other"),
    (r"(?<!sub)(?<!sous-)cutan|topical|(?<!intra)dermal|cream|ointment|gel\b|lotion", "topical", "other"),
    (r"inhal|nebuli|respiratory|pulmonary|nasal", "inhaled/nasal", "other"),
    (r"rectal|vaginal|urethral|auricular|otic|dental|oromucosal|buccal|sublingual|gingival", "other mucosal", "other"),
]


def norm_routes(text):
    out = set()
    t = str(text).lower()
    for rx, name, cat in ROUTE_RULES:
        if re.search(rx, t):
            out.add((name, cat))
    if any(n in ("intravenous", "intramuscular", "subcutaneous", "intravitreal", "other injection")
           for n, _ in out):
        out.discard(("injection", INJECTION))          # keep the specific injection route only
    return out


ZH_FORMS = [  # dosage form in a Chinese product name -> route
    (r"注射|输液|针", ("injection", INJECTION)),
    (r"片|胶囊|颗粒|口服|干混悬|混悬剂|散$|丸|糖浆|滴丸", ("oral", "oral")),
    (r"滴眼|眼膏|眼用", ("ophthalmic", "other")),
    (r"乳膏|软膏|凝胶|搽剂|贴|酊|洗剂|外用", ("topical", "other")),
    (r"栓|阴道|直肠", ("other mucosal", "other")),
    (r"吸入|喷雾|气雾", ("inhaled/nasal", "other")),
]


def zh_routes(name):
    out = set()
    for rx, r in ZH_FORMS:
        if re.search(rx, str(name)):
            out.add(r)
    return out


def name_keys(name):
    """Normalised lookup keys of a drug name (salt / hydrate / dosage-form words stripped)."""
    return set(cd.candidate_keys(name))


# ----------------------------------------------------------------------------- evidence collectors
def fda_routes(index):
    z = zipfile.ZipFile(cd.cached_download("openfda_drugsfda.json.zip", "https://download.open.fda.gov/drug/drugsfda/drug-drugsfda-0001-of-0001.json.zip"))
    data = json.load(z.open(z.namelist()[0]))["results"]
    res, by_key = defaultdict(set), defaultdict(set)
    cache = {}
    for app in data:
        for p in app.get("products") or []:
            routes = norm_routes(p.get("route", "")) or norm_routes(p.get("dosage_form", ""))
            for ai in p.get("active_ingredients") or []:
                n = ai.get("name", "")
                if n not in cache:
                    cache[n] = index.lookup(n)
                if cache[n]:
                    res[cache[n]] |= routes
                for k in name_keys(n):
                    by_key[k] |= routes
    return res, by_key


def article57_routes(index):
    x = pd.read_excel(cd.cached_download("ema_article57.xlsx",
                      "https://www.ema.europa.eu/en/documents/other/article-57-product-data_en.xlsx"),
                      header=None, skiprows=20, usecols=[0, 1, 2, 3], names=["product", "substance", "route", "country"],
                      dtype=str).dropna()
    x = x.drop_duplicates()[["substance", "route"]]     # one row per product and country
    # Article 57 has data-entry errors (e.g. a few generic caspofungin products coded "Oral Use"), so a
    # route counts only if >= 10% of that substance's EU/EEA product rows report it.
    n_rows = defaultdict(int)
    cnt = defaultdict(lambda: defaultdict(int))
    cache = {}
    for subst, route in x.itertuples(index=False):
        routes = norm_routes(route)
        for part in re.split(r"\s*\|\s*|,\s*|\s\+\s|/", subst):
            part = part.strip()
            if not (2 < len(part) < 80):
                continue
            if part not in cache:
                cache[part] = index.lookup(part)
            ids = ([("cid", cache[part])] if cache[part] else []) + [("key", k) for k in name_keys(part)]
            for i in ids:
                n_rows[i] += 1
                for r in routes:
                    cnt[i][r] += 1
    res, by_key = defaultdict(set), defaultdict(set)
    for i, c in cnt.items():
        keep = {r for r, n in c.items() if n >= 0.1 * n_rows[i]}
        (res if i[0] == "cid" else by_key)[i[1]] |= keep
    return res, by_key


def france_routes(index, parents):
    raw = ROOT / "pipeline" / "eu_national" / "raw"
    b = pd.read_csv(raw / "CIS_bdpm.txt", sep="\t", header=None, encoding="latin1", dtype=str, quoting=3,
                    on_bad_lines="skip")
    c = pd.read_csv(raw / "CIS_COMPO_bdpm.txt", sep="\t", header=None, encoding="latin1", dtype=str, quoting=3,
                    on_bad_lines="skip")
    c.columns = ["cis", "form", "code", "name", "dose", "ref", "nature", "link"]
    c = c.merge(b[[0, 3]].rename(columns={0: "cis", 3: "voies"}), on="cis")
    names = c.name.str.lower().unique()
    hit = map_names(names, "fr", index, parents, ["fr"],
                    ROOT / "pipeline" / "eu_national" / "wikidata_fr_hits.csv")
    res = defaultdict(set)
    for r in c.itertuples():
        cid = hit.get(str(r.name).lower())
        if cid:
            res[cid] |= norm_routes(r.voies) | norm_routes(r.form)
    return res


def china_routes():
    """Return {InChIKey: set(routes)} and {lower-case english name: set(routes)} from NMPA product names."""
    by_ik, by_name = defaultdict(set), defaultdict(set)      # by_ik is keyed by the InChIKey skeleton (14 chars)
    cand = pd.read_csv(ROOT / "data" / "nmpa_candidates.csv", dtype=str)
    pc = pd.read_csv(ROOT / "data" / "pubchem_results.csv", dtype=str)
    ik_of = dict(zip(pc["query"], pc.inchikey))
    for r in cand.itertuples():
        rs = zh_routes(r.zh)
        if not rs:
            continue
        ik = ik_of.get(r.key2)
        if isinstance(ik, str):
            by_ik[ik[:14]] |= rs
        for k in name_keys(r.key2):                         # English name of the NMPA record
            by_name[k] |= rs
    # records rescued through Chinese names (Wikidata) -> InChIKey of the matched item
    clean = pd.read_csv(ROOT / "pipeline" / "nmpa_fix" / "nmpa_candidates_clean.csv", dtype=str)
    wd = pd.read_csv(ROOT / "pipeline" / "nmpa_fix" / "wikidata_zh_hits.csv", dtype=str).dropna(subset=["ik"])
    lab2ik = wd.groupby("lab").ik.agg(lambda s: s.iloc[0] if s.nunique() == 1 else None).dropna().to_dict()
    for r in clean.itertuples():
        ik = lab2ik.get(r.zh_clean) or lab2ik.get(r.zh_base)
        if ik:
            by_ik[ik[:14]] |= zh_routes(r.zh)
    # web-verified China list and hand-made additions (Chinese product names / notes)
    china = pd.read_csv(ROOT / "pipeline" / "nmpa_fix" / "china_innovative_verified_Claude.csv", dtype=str)
    for r in china.itertuples():
        for n in r.INN.split("/"):
            for k in name_keys(n):
                by_name[k] |= zh_routes(r.Chinese_name)
    va = pd.read_csv(ROOT / "pipeline" / "curation" / "verified_additions_Claude.csv", dtype=str,
                     keep_default_na=False)
    for r in va.itertuples():
        for k in name_keys(r.Name):
            by_name[k] |= zh_routes(r.source)
    return by_ik, by_name


def chembl_flags(cids):
    f = CACHE / "chembl_flags.json"
    have = json.loads(f.read_text()) if f.exists() else {}
    todo = [c for c in cids if c not in have]
    for i in range(0, len(todo), 40):
        chunk = ";".join(todo[i:i + 40])
        d = cd.http_get(f"{cd.CHEMBL}/molecule/set/{chunk}.json").json()
        for m in d["molecules"]:
            have[m["molecule_chembl_id"]] = {k: m.get(k) for k in ("oral", "parenteral", "topical", "max_phase")}
        time.sleep(0.3)
    f.write_text(json.dumps(have, indent=0))
    return have


def ctgov_routes(names):
    """{query name: {'oral': n_studies, 'injection': n, 'other': n, 'ncts': [...]}} from intervention text."""
    f = CACHE / "ctgov_routes.json"
    have = json.loads(f.read_text()) if f.exists() else {}
    for q in names:
        if q in have:
            continue
        url = ("https://clinicaltrials.gov/api/v2/studies?pageSize=100&fields=NCTId,ArmsInterventionsModule"
               "&query.intr=" + urllib.parse.quote(q))
        try:
            d = json.load(urllib.request.urlopen(urllib.request.Request(
                url, headers={"User-Agent": "drug-table-crosscheck/1.0"}), timeout=60))
        except Exception as e:
            print("  ctgov failed", q, e)
            continue
        tally = {"oral": 0, INJECTION: 0, "other": 0, "routes": [], "ncts": []}
        key = re.sub(r"[^a-z0-9]", "", q.lower())
        for s in d.get("studies", []):
            mod = s["protocolSection"].get("armsInterventionsModule", {})
            texts = []
            for iv in mod.get("interventions", []):
                nm = " ".join([iv.get("name", "")] + (iv.get("otherNames") or []))
                if key and key in re.sub(r"[^a-z0-9]", "", nm.lower()):
                    texts.append(f"{iv.get('name', '')} {iv.get('description', '')}")
            if not texts:
                continue
            arms = " ".join(a.get("description", "") for a in mod.get("armGroups", []))
            rs = norm_routes(" ".join(texts)) or norm_routes(arms)
            if not rs:
                continue
            tally["n"] = tally.get("n", 0) + 1
            for c in {c for _, c in rs}:
                tally[c] += 1
            for n, _ in rs:
                tally.setdefault("route_n", {})
                tally["route_n"][n] = tally["route_n"].get(n, 0) + 1
            tally["ncts"].append(s["protocolSection"]["identificationModule"]["nctId"])
        tally["ncts"] = tally["ncts"][:5]
        have[q] = tally
        time.sleep(0.3)
    f.write_text(json.dumps(have, indent=0, ensure_ascii=False))
    return have


# ----------------------------------------------------------------------------- classification
def classify(routes):
    cats = {c for _, c in routes}
    if "oral" in cats and INJECTION in cats:
        return "Oral + Injection"
    if "oral" in cats:
        return "Oral"
    if INJECTION in cats:
        return "Injection"
    if "other" in cats:
        return "Non-oral (other)"
    return "Unknown"


def fmt(routes):
    order = ["oral", "intravenous", "intramuscular", "subcutaneous", "intravitreal", "other injection",
             "injection", "topical", "ophthalmic", "inhaled/nasal", "other mucosal"]
    return "; ".join(sorted({n for n, _ in routes}, key=order.index))


def main():
    mc = pd.read_csv(MC_FILE, dtype=str, keep_default_na=False)
    mc = mc.drop(columns=[c for c in ("Route_Class", "Routes", "Route_Source") if c in mc])
    det = pd.read_csv(ROOT / "results" / "final_small_molecule_drugs_detailed_Claude.csv", dtype=str,
                      keep_default_na=False).set_index("Name")
    parents, child_of, index = cd.load_chembl()
    ik2cid = {p["inchikey"]: cid for cid, p in parents.items() if p.get("inchikey")}

    print("collecting FDA / EU Article 57 / France / China product routes ...")
    (fda_c, fda_k), (a57_c, a57_k) = fda_routes(index), article57_routes(index)
    src = {"FDA": fda_c, "EU-A57": a57_c, "France": france_routes(index, parents)}
    by_key = {"FDA": fda_k, "EU-A57": a57_k}
    zh_ik, zh_name = china_routes()

    rows = []
    for r in mc.itertuples():
        names = det.loc[r.Name, "Names_in_sources"].split("; ") if r.Name in det.index else [r.Name]
        ik = det.loc[r.Name, "InChIKey"] if r.Name in det.index else ""
        cid = ik2cid.get(ik) or next((c for c in map(index.lookup, names) if c), None)
        ev = {k: set(v.get(cid, set())) if cid else set() for k, v in src.items()}
        keys = set().union(*(name_keys(n) for n in names + [r.Name]))
        for k, d in by_key.items():                     # also link products by drug name
            for key in keys:
                ev[k] |= d.get(key, set())
        ev["China"] = set(zh_ik.get(ik[:14], set())) if ik else set()
        for key in keys:
            ev["China"] |= zh_name.get(key, set())
        if r.Name.startswith("Carrimycin component"):
            ev["China"] |= {("oral", "oral")}               # 可利霉素片 (tablets)
        rows.append({"Name": r.Name, "cid": cid, "names": names, **ev})

    flags = chembl_flags(sorted({x["cid"] for x in rows if x["cid"]}))
    approved = mc[AG].apply(lambda c: c.str.fullmatch(r"\d{4}") | c.eq("Approved")).any(axis=1)
    queries = {}
    for x, appr in zip(rows, approved):
        product = set().union(*(x[k] for k in ("FDA", "EU-A57", "France", "China")))
        if not product:
            q = x["Name"].split(": ")[-1]
            q = "carrimycin" if x["Name"].startswith("Carrimycin") else q
            queries[x["Name"]] = q
    print(f"querying ClinicalTrials.gov for {len(queries)} macrocycles without marketed-product routes ...")
    ct = ctgov_routes(sorted(set(queries.values())))

    out = []
    for x, appr in zip(rows, approved):
        product = set().union(*(x[k] for k in ("FDA", "EU-A57", "France", "China")))
        fl = flags.get(x["cid"] or "", {})
        flag_routes = set()
        if fl.get("oral"):
            flag_routes.add(("oral", "oral"))
        if fl.get("parenteral"):
            flag_routes.add(("injection", INJECTION))
        if fl.get("topical"):
            flag_routes.add(("topical", "other"))
        t = ct.get(queries.get(x["Name"], ""), {})
        trial_routes = set()
        for c, n in (("oral", "oral"), (INJECTION, "injection"), ("other", "topical")):
            if t.get(c):
                trial_routes.add((n, c))
        cat_of = {rn: rc for _, rn, rc in ROUTE_RULES}
        n_st = t.get("n", 0)
        keep = {c for c in ("oral", INJECTION, "other") if n_st and t.get(c, 0) >= max(1, n_st / 3)}
        trial_named = {(n, cat_of.get(n, "other")) for n, k in (t.get("route_n") or {}).items()
                       if cat_of.get(n, "other") in keep and k >= max(1, n_st / 3)}
        if product:
            routes, why = product, "marketed products: " + ", ".join(k for k in ("FDA", "EU-A57", "France", "China")
                                                                       if x[k])
        elif trial_named:
            routes = trial_named
            why = "clinical trials (" + ", ".join(t.get("ncts", [])[:3]) + ")"
        elif flag_routes:
            routes, why = flag_routes, f"ChEMBL route flags ({x['cid']})"
        else:
            routes, why = set(), "no route evidence found"
        out.append({"Name": x["Name"], "Route_Class": classify(routes), "Routes": fmt(routes), "Route_Source": why,
                    "FDA_routes": fmt(x["FDA"]), "EU_A57_routes": fmt(x["EU-A57"]),
                    "France_routes": fmt(x["France"]), "China_routes": fmt(x["China"]),
                    "Trial_studies_with_route": t.get("n", ""), "Trial_studies_oral": t.get("oral", ""),
                    "Trial_studies_injection": t.get(INJECTION, ""), "Trial_studies_other": t.get("other", ""),
                    "ChEMBL_oral": fl.get("oral", ""), "ChEMBL_parenteral": fl.get("parenteral", ""),
                    "ChEMBL_topical": fl.get("topical", ""), "ChEMBL_ID": x["cid"] or "",
                    "Approved_somewhere": bool(appr)})
    o = pd.DataFrame(out)
    ov_file = HERE / "route_overrides_Claude.csv"
    if ov_file.exists():                          # web-checked corrections
        ov = pd.read_csv(ov_file, dtype=str, keep_default_na=False)
        for r in ov.itertuples():
            m = o.Name.eq(r.Name)
            o.loc[m, ["Route_Class", "Routes"]] = [r.Route_Class, r.Routes]
            o.loc[m, "Route_Source"] = "verified: " + r.source
    o.to_csv(OUT_DETAIL, index=False)

    cols = ["Route_Class", "Routes", "Route_Source"]
    add = o.set_index("Name")[cols]
    for f in (MC_FILE, MOLD_FILE):
        d = pd.read_csv(f, dtype=str, keep_default_na=False)
        d = d.drop(columns=[c for c in cols if c in d])
        pos = d.columns.get_loc("ANVISA") + 1
        for i, c in enumerate(cols):
            d.insert(pos + i, c, d.Name.map(add[c]))
        d.to_csv(f, index=False)
    print(o.Route_Class.value_counts().to_string())
    print(pd.crosstab(o.Route_Class, o.Approved_somewhere))


if __name__ == "__main__":
    main()
