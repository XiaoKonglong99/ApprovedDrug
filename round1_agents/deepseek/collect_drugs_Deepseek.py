#!/usr/bin/env python3
"""
Collect small-molecule drug approvals from FDA (US), EMA (EU), NMPA (China),
and PMDA (Japan), resolve SMILES, and annotate the highest clinical-trial
phase reached in each agency's jurisdiction.

Output columns:
    Name, SMILES, FDA, EMA, NMPA, PMDA

A cell is:
    * an approval year (e.g. 2019),
    * a clinical-trial phase (e.g. Phase 3) when the drug is not approved by
      that agency but has an interventional trial with a site in its region,
    * "NA" when no approval/trial record was found.

The script uses public sources and a local cache.  The first run is slow
(PubChem and ClinicalTrials.gov are queried); subsequent runs reuse the cache.

Dependencies:
    pandas, requests, rdkit, openpyxl, lxml, pg8000 (only for --refresh)

Usage:
    python collect_drugs_Deepseek.py --out drug_approvals_smiles_Deepseek.csv
    python collect_drugs_Deepseek.py --refresh          # re-download everything
    python collect_drugs_Deepseek.py --skip-trials      # approvals only
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import sys
import time
import urllib.parse
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]   # repo root (shared data/ folder); was .parent in round 1
DATA = ROOT / "data"
CACHE = ROOT / "data" / "cache_deepseek"
for p in (DATA, CACHE):
    p.mkdir(exist_ok=True)


def log(msg: str) -> None:
    print(f"[collect] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# Download helpers
# --------------------------------------------------------------------------- #
def download(url: str, path: Path, timeout: int = 180) -> Path:
    if path.exists() and path.stat().st_size > 0:
        return path
    log(f"downloading {url}")
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    path.write_bytes(r.content)
    return path


# --------------------------------------------------------------------------- #
# Source 1: DrugCentral public PostgreSQL (FDA / EMA / PMDA approvals)
# --------------------------------------------------------------------------- #
def fetch_drugcentral() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read approval dates and synonyms from DrugCentral's public database."""
    import pg8000  # pip install pg8000

    log("querying DrugCentral public PostgreSQL")
    con = pg8000.connect(
        host="unmtid-dbs.net", port=5433, database="drugcentral",
        user="drugman", password="dosage", timeout=120,
    )
    cur = con.cursor()
    cur.execute("""
        SELECT s.id, s.name, s.smiles, s.inchi, s.inchikey, s.cd_formula,
               s.cd_molweight, s.cas_reg_no,
               min(a.approval) FILTER (WHERE a.type='FDA')  AS fda_date,
               min(a.approval) FILTER (WHERE a.type='EMA')  AS ema_date,
               min(a.approval) FILTER (WHERE a.type='PMDA') AS pmda_date,
               min(a.approval) FILTER (WHERE a.type='China Food and Drug Administration (CFDA)') AS cfda_date,
               count(*) FILTER (WHERE a.type='FDA')  AS fda_n,
               count(*) FILTER (WHERE a.type='EMA')  AS ema_n,
               count(*) FILTER (WHERE a.type='PMDA') AS pmda_n,
               count(*) FILTER (WHERE a.type='China Food and Drug Administration (CFDA)') AS cfda_n
        FROM structures s
        JOIN approval a ON s.id = a.struct_id
        WHERE a.type IN ('FDA','EMA','PMDA','China Food and Drug Administration (CFDA)')
        GROUP BY s.id, s.name, s.smiles, s.inchi, s.inchikey, s.cd_formula, s.cd_molweight, s.cas_reg_no
        ORDER BY s.id
    """)
    rows = cur.fetchall()
    cols = [d[0] for d in cur.description]
    dc = pd.DataFrame(rows, columns=cols)

    cur.execute("""
        SELECT sy.id, sy.name, sy.preferred_name, s.name AS struct_name
        FROM synonyms sy
        JOIN structures s ON sy.id = s.id
        WHERE sy.id IN (
            SELECT struct_id FROM approval
            WHERE type IN ('FDA','EMA','PMDA','China Food and Drug Administration (CFDA)')
        )
    """)
    srows = cur.fetchall()
    scols = [d[0] for d in cur.description]
    syn = pd.DataFrame(srows, columns=scols)
    con.close()

    dc.to_csv(DATA / "dc_approvals.csv", index=False)
    syn.to_csv(DATA / "dc_synonyms.csv", index=False)
    return dc, syn


def load_drugcentral(refresh: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    ap, sp = DATA / "dc_approvals.csv", DATA / "dc_synonyms.csv"
    if refresh or not ap.exists() or not sp.exists():
        return fetch_drugcentral()
    return pd.read_csv(ap), pd.read_csv(sp)


# --------------------------------------------------------------------------- #
# Source 2: EMA download table
# --------------------------------------------------------------------------- #
def fetch_ema() -> pd.DataFrame:
    import openpyxl

    url = "https://www.ema.europa.eu/en/documents/report/medicines-output-medicines-report_en.xlsx"
    path = download(url, DATA / "ema_medicines.xlsx")
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb["Medicine"]
    out = []
    for r in ws.iter_rows(min_row=10, values_only=True):
        if not r or not r[1]:
            continue
        out.append([
            r[0], r[1], r[2], r[3], r[4], r[6], r[7], r[11], r[31], r[26], r[29], r[38],
        ])
    df = pd.DataFrame(out, columns=[
        "category", "name", "ema_no", "status", "opinion_status", "inn",
        "active_substance", "atc", "ma_date", "ec_date", "opinion_date", "url",
    ])
    df.to_csv(DATA / "ema_medicines.csv", index=False)
    return df


def load_ema(refresh: bool = False) -> pd.DataFrame:
    p = DATA / "ema_medicines.csv"
    if refresh or not p.exists():
        return fetch_ema()
    return pd.read_csv(p)


# --------------------------------------------------------------------------- #
# Source 3: NMPA public database snapshot (GitHub mirror)
# --------------------------------------------------------------------------- #
NMPA_DOM = ("https://raw.githubusercontent.com/xprincehf/nmpa-data/master/"
            "%E8%8D%AF%E5%93%81/2020/%E5%9B%BD%E4%BA%A7%E8%8D%AF%E5%93%81/"
            "%E5%9B%BD%E4%BA%A7%E8%8D%AF%E5%93%81.json.zip")
NMPA_IMP = ("https://raw.githubusercontent.com/xprincehf/nmpa-data/master/"
            "%E8%8D%AF%E5%93%81/2020/%E8%BF%9B%E5%8F%A3%E8%8D%AF%E5%93%81/"
            "%E8%BF%9B%E5%8F%A3%E8%8D%AF%E5%93%81.json")


def _approval_number_year(no: str) -> int | None:
    """Decode the year embedded in a Chinese drug approval number."""
    s = str(no or "").upper().replace("国药准字", "").strip()
    m = re.search(r"([A-Z]+)(\d{8})", s)
    if not m:
        return None
    d = m.group(2)
    first4 = int(d[:4])
    if 1980 <= first4 <= 2026:          # new format: H20123350
        return first4
    prov, yy = int(d[:2]), int(d[2:4])  # old format: H33020747
    if 11 <= prov <= 65 and (yy >= 80 or yy <= 26):
        return 1900 + yy if yy >= 80 else 2000 + yy
    return None


def fetch_nmpa() -> pd.DataFrame:
    dom = download(NMPA_DOM, DATA / "nmpa_domestic.json.zip", timeout=300)
    imp = download(NMPA_IMP, DATA / "nmpa_imported.json", timeout=300)

    def records():
        z = zipfile.ZipFile(dom)
        for n in z.namelist():
            for r in json.loads(z.read(n)):
                yield "domestic", r
        for r in json.loads(imp.read_bytes()):
            yield "imported", r

    agg: dict[str, dict] = {}
    for src, r in records():
        if (r.get("产品类别") or "").strip() != "化学药品":
            continue
        if src == "domestic":
            en = r.get("英文名称") or ""
            zh = r.get("产品名称") or ""
            dt = r.get("批准日期") or ""
            no = r.get("批准文号") or ""
        else:
            en = r.get("产品名称（英文）") or ""
            zh = r.get("产品名称（中文）") or ""
            dt = r.get("发证日期") or ""
            no = r.get("注册证号") or ""
        key = clean_name(en)
        if not key or not candidate_name(key):
            continue
        dy = int(dt[:4]) if re.match(r"^\d{4}", dt) else None
        ny = _approval_number_year(no)
        years = [y for y in (dy, ny) if y]
        if not years:
            continue
        y = min(years)
        a = agg.setdefault(key, {
            "key2": key, "en_raw": en, "zh": zh, "year": y,
            "date_year": dy, "num_year": ny, "n_records": 0, "sources": set(),
        })
        a["n_records"] += 1
        a["sources"].add(src)
        if y < a["year"]:
            a["year"], a["en_raw"], a["zh"] = y, en, zh
    rows = []
    for a in agg.values():
        rows.append({
            "key2": a["key2"], "en_raw": a["en_raw"], "zh": a["zh"],
            "year": a["year"], "date_year": a["date_year"],
            "num_year": a["num_year"], "n_records": a["n_records"],
            "sources": "|".join(sorted(a["sources"])),
        })
    df = pd.DataFrame(rows).sort_values(["year", "key2"]).reset_index(drop=True)
    df.to_csv(DATA / "nmpa_candidates.csv", index=False)
    return df


def load_nmpa(refresh: bool = False) -> pd.DataFrame:
    p = DATA / "nmpa_candidates.csv"
    if refresh or not p.exists():
        return fetch_nmpa()
    return pd.read_csv(p)


# --------------------------------------------------------------------------- #
# Source 4: FDA novel approvals 2025-2026
# --------------------------------------------------------------------------- #
def fetch_fda_recent() -> pd.DataFrame:
    frames = []
    for y in (2025, 2026):
        url = f"https://www.fda.gov/drugs/novel-drug-approvals-fda/novel-drug-approvals-{y}"
        t = pd.read_html(url)[0]
        t.columns = ["no", "brand", "active", "date", "use"]
        t["year"] = y
        frames.append(t)
    df = pd.concat(frames, ignore_index=True)
    df["date"] = pd.to_datetime(df["date"], format="%m/%d/%Y", errors="coerce")
    df = df[["year", "brand", "active", "date", "use"]]
    df.to_csv(DATA / "fda_recent.csv", index=False)
    return df


def load_fda_recent(refresh: bool = False) -> pd.DataFrame:
    p = DATA / "fda_recent.csv"
    if refresh or not p.exists():
        return fetch_fda_recent()
    return pd.read_csv(p)


# --------------------------------------------------------------------------- #
# Name normalisation
# --------------------------------------------------------------------------- #
DOSAGE = r"""(?:injection|injectable|infusion|tablets?|capsules?|granules?|granule|powder|suspension|solution|drops?|cream|ointment|gel|sprays?|syrup|emulsion|suppositor(?:y|ies)|patch(?:es)?|film(?:s)?|inhalation|aerosol|implant(?:s)?|lozenges?|mouthwash|enema|liniment|tincture|sachets?|lotion|dispersible|orally disintegrating|effervescent|enteric[- ]coated|sustained[- ]release|extended[- ]release|modified[- ]release|controlled[- ]release|delayed[- ]release|prolonged[- ]release|immediate[- ]release|liposome|microsphere|nanoparticle|depot|oral|eye|ophthalmic|nasal|topical|transdermal|rectal|vaginal|inhalational|nebul(?:e|iser|izer)|sr|cr|er|xr|od)"""
SALT = r"""(?:hydrochloride|hydrobromide|hydroiodide|sulfate|sulphate|phosphate|acetate|citrate|tartrate|maleate|fumarate|succinate|mesylate|tosylate|besylate|napsylate|ethylsuccinate|dihydrate|monohydrate|trihydrate|tetrahydrate|anhydrous|hydrate|hydrates|sesquihydrate|hemihydrate|calcium|sodium|potassium|magnesium|meglumine|tromethamine|lysine|arginine|choline)"""
INORGANIC_KEEP = {
    "sodium chloride", "potassium chloride", "magnesium sulfate", "magnesium sulphate",
    "calcium chloride", "sodium bicarbonate", "calcium carbonate", "ferrous sulfate",
    "ferrous sulphate", "zinc sulfate", "aluminum hydroxide", "magnesium hydroxide",
    "barium sulfate", "sodium phosphate", "potassium phosphate", "calcium gluconate",
    "sodium lactate", "ammonium chloride", "silver nitrate", "hydrogen peroxide",
}


def clean_name(raw) -> str:
    if not raw:
        return ""
    s = str(raw).strip().lower()
    if s in INORGANIC_KEEP:
        return s
    s = re.sub(r"\([^)]*\)", " ", s)
    s = s.replace("&", " and ")
    s = re.sub(r"\bfor (?:injection|infusion|oral|intravenous|iv|im|subcutaneous|external|topical|ophthalmic|nasal|inhalation|depot)\b.*$", " ", s)
    s = re.sub(r"\b" + DOSAGE + r"\b.*$", " ", s)
    s = re.sub(r"\b(?:sterile|micronized|micronised|usp|bp|ep|jp|ph\.?\s*eur\.?|grade)\b", " ", s)
    for _ in range(3):
        s2 = re.sub(r"[,;/\-\s]*\b" + SALT + r"\b\s*$", "", s)
        if s2 == s:
            break
        s = s2
    s = re.sub(r"\b(?:and|with)\s*$", " ", s)
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def clean_name2(raw) -> str:
    s = clean_name(raw)
    for _ in range(3):
        s2 = re.sub(r"[,;/\-\s]*\b(?:besilate|mesilate|dihydrochloride|dihydrobromide|hemifumarate|hemihydrate|trometamol)\b\s*$", "", s)
        if s2 == s:
            break
        s = s2
    return s


EXCLUDE_SUBSTR = [
    "compound", "vitamin", "amino acid", "glucose", "water", "oxygen",
    "sodium chloride", "potassium chloride", "calcium chloride", "sodium bicarbonate",
    "albumin", "heparin", "insulin", "vaccine", "toxoid", "immunoglobulin", "interferon",
    "interleukin", "plasma", "blood", "serum", "antibody", "monoclonal", "enzyme", "toxin",
    "allergen", "extract", "pollen", "bacteria", "virus", "yeast", "tissue", "protein",
    "dextran", "starch", "cellulose", "glycerin", "glycerol", "propylene glycol",
    "polyethylene glycol", "sorbitol", "mannitol", "lactose", "sucrose", "maltose",
    "cyclodextrin", "paraffin", "petrolatum", "lanolin", "silicone", "dimethicone",
    "ethanol", "alcohol", "iodine", "iodinated", "nutrient", "electrolyte", "enteral",
    "parenteral nutrition", "fat emulsion", "fish oil", "contrast medium", "radiopharmaceutical",
]


def candidate_name(s: str) -> bool:
    if not s or len(s) < 4:
        return False
    if any(x in s for x in EXCLUDE_SUBSTR):
        return False
    if " and " in s or " with " in s or "/" in s or "+" in s:
        return False
    return True


def title_name(s: str) -> str:
    s = str(s).strip()
    m = re.search(r"[A-Za-z]", s)
    if not m:
        return s
    i = m.start()
    return s[:i] + s[i].upper() + s[i + 1:]


# --------------------------------------------------------------------------- #
# Chemistry helpers
# --------------------------------------------------------------------------- #
def parent_mol(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
    except Exception:
        return None
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True)
    if not frags:
        return None
    return max(frags, key=lambda m: (
        sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C"), m.GetNumHeavyAtoms()
    ))


def struct_from_smiles(smi):
    m = parent_mol(smi)
    if m is None:
        return None
    mw = Descriptors.MolWt(m)
    c = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C")
    ha = m.GetNumHeavyAtoms()
    if c < 1 or ha < 3 or mw > 2000 or mw < 40:
        return None
    out = Chem.MolToSmiles(m, isomericSmiles=True)
    if "*" in out:
        return None
    try:
        return {
            "smiles": out,
            "inchikey": Chem.MolToInchiKey(m),
            "mw": round(mw, 2),
            "formula": rdMolDescriptors.CalcMolFormula(m),
        }
    except Exception:
        return None


def year_of(x):
    try:
        if pd.isna(x):
            return None
        s = str(x)
        if s[:4].isdigit():
            return int(s[:4])
    except Exception:
        pass
    return None


def min_year(a, b):
    vals = [x for x in (a, b) if x is not None]
    return min(vals) if vals else None


# --------------------------------------------------------------------------- #
# PubChem resolver
# --------------------------------------------------------------------------- #
PC_FIELDS = ["query", "cid", "smiles", "connectivity_smiles", "inchikey", "formula", "mw", "status"]


def pubchem_resolve(names: list[str]) -> dict[str, dict]:
    cache_path = DATA / "pubchem_results.csv"
    done = {}
    if cache_path.exists():
        for r in csv.DictReader(open(cache_path, encoding="utf-8")):
            done[r["query"]] = r
    todo = [n for n in dict.fromkeys(names) if n and n not in done]
    if todo:
        log(f"PubChem: resolving {len(todo)} names")
        lock = __import__("threading").Lock()
        next_t = [0.0]
        local = __import__("threading").local()

        def rate_wait():
            with lock:
                now = time.time()
                t = max(now, next_t[0])
                next_t[0] = t + 0.24
            if t > now:
                time.sleep(t - now)

        def get_session():
            if not hasattr(local, "s"):
                local.s = requests.Session()
                local.s.headers.update({"User-Agent": "drug-approval-collector/1.0 (research)"})
            return local.s

        def fetch(name):
            url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
                   + urllib.parse.quote(name, safe="")
                   + "/property/SMILES,ConnectivitySMILES,InChIKey,MolecularFormula,MolecularWeight/JSON")
            s = get_session()
            for attempt in range(5):
                rate_wait()
                try:
                    r = s.get(url, timeout=30)
                    if r.status_code == 200:
                        p = r.json()["PropertyTable"]["Properties"][0]
                        return {
                            "query": name, "cid": p.get("CID", ""),
                            "smiles": p.get("SMILES", ""),
                            "connectivity_smiles": p.get("ConnectivitySMILES", ""),
                            "inchikey": p.get("InChIKey", ""),
                            "formula": p.get("MolecularFormula", ""),
                            "mw": p.get("MolecularWeight", ""), "status": "ok",
                        }
                    if r.status_code == 404:
                        return {"query": name, "status": "not_found"}
                    if r.status_code in (429, 500, 502, 503, 504):
                        time.sleep(3 * (attempt + 1))
                        continue
                    return {"query": name, "status": f"http_{r.status_code}"}
                except Exception:
                    time.sleep(2 * (attempt + 1))
            return {"query": name, "status": "error"}

        new_file = not cache_path.exists()
        f = open(cache_path, "a", newline="", encoding="utf-8")
        w = csv.DictWriter(f, fieldnames=PC_FIELDS)
        if new_file:
            w.writeheader()
        with ThreadPoolExecutor(max_workers=4) as ex:
            futs = {ex.submit(fetch, n): n for n in todo}
            for i, fut in enumerate(as_completed(futs), 1):
                rec = fut.result()
                w.writerow({k: rec.get(k, "") for k in PC_FIELDS})
                f.flush()
                if i % 250 == 0:
                    log(f"PubChem: {i}/{len(todo)}")
        f.close()
        done = {}
        for r in csv.DictReader(open(cache_path, encoding="utf-8")):
            done[r["query"]] = r
    out = {}
    for q, r in done.items():
        if r.get("status") == "ok":
            s = struct_from_smiles(r.get("smiles"))
            if s:
                out[q.lower()] = s
    return out


# --------------------------------------------------------------------------- #
# Build the approval table
# --------------------------------------------------------------------------- #
def build_table(dc, syn, nmpa, ema, fda_recent, pc):
    groups: dict[str, dict] = {}
    name2ik: dict[str, str] = {}
    ok14 = set()
    p_ok = DATA / "dc_fda_undated_ok.csv"
    if p_ok.exists():
        ok = pd.read_csv(p_ok)
        ok14 = set(str(x)[:14] for x in ok.inchikey.dropna())

    def G(ik):
        if ik not in groups:
            groups[ik] = {
                "inchikey": ik, "name": None, "smiles": None, "mw": None,
                "formula": None, "fda": None, "ema": None, "nmpa": None,
                "pmda": None,
                "fda_approved": False, "ema_approved": False,
                "nmpa_approved": False, "pmda_approved": False,
                "aliases": set(), "nmpa_names": [],
            }
        return groups[ik]

    # 1. DrugCentral structures with SMILES
    dc_rows = []
    for r in dc.itertuples():
        if pd.isna(r.smiles):
            continue
        s = struct_from_smiles(r.smiles)
        if not s:
            continue
        s.update({"name": str(r.name).strip(), "fda": year_of(r.fda_date),
                  "ema": year_of(r.ema_date), "pmda": year_of(r.pmda_date),
                  "fda_approved": bool(r.fda_n > 0) and (pd.notna(r.fda_date) or str(r.inchikey)[:14] in ok14),
                  "ema_approved": bool(r.ema_n > 0), "pmda_approved": bool(r.pmda_n > 0)})
        dc_rows.append(s)
    for s in dc_rows:
        g = G(s["inchikey"])
        g["name"], g["smiles"], g["mw"], g["formula"] = s["name"], s["smiles"], s["mw"], s["formula"]
        g["fda"], g["ema"], g["pmda"] = s["fda"], s["ema"], s["pmda"]
        g["fda_approved"], g["ema_approved"], g["pmda_approved"] = s["fda_approved"], s["ema_approved"], s["pmda_approved"]
        g["aliases"].add(s["name"])
        name2ik[s["name"].lower()] = s["inchikey"]
        name2ik[clean_name2(s["name"])] = s["inchikey"]

    # 2. DrugCentral rows without SMILES, resolved by PubChem
    for r in dc[dc.smiles.isna()].itertuples():
        s = pc.get(str(r.name).lower()) or pc.get(clean_name2(str(r.name)))
        if not s:
            continue
        g = G(s["inchikey"])
        if g["name"] is None:
            g.update({"name": r.name, "smiles": s["smiles"], "mw": s["mw"], "formula": s["formula"]})
        g["fda"] = min_year(g["fda"], year_of(r.fda_date))
        g["ema"] = min_year(g["ema"], year_of(r.ema_date))
        g["pmda"] = min_year(g["pmda"], year_of(r.pmda_date))
        g["fda_approved"] = g["fda_approved"] or (bool(r.fda_n > 0) and (pd.notna(r.fda_date) or str(r.inchikey)[:14] in ok14))
        g["ema_approved"] = g["ema_approved"] or bool(r.ema_n > 0)
        g["pmda_approved"] = g["pmda_approved"] or bool(r.pmda_n > 0)
        g["aliases"].add(str(r.name))
        name2ik[str(r.name).lower()] = s["inchikey"]
        name2ik[clean_name2(str(r.name))] = s["inchikey"]

    # 3. Synonyms
    for r in syn.itertuples():
        ik = name2ik.get(str(r.struct_name).lower()) or name2ik.get(clean_name2(str(r.struct_name)))
        if ik:
            groups[ik]["aliases"].add(str(r.name))

    # 4. NMPA
    for r in nmpa.itertuples():
        key = str(r.key2).lower()
        s = pc.get(key)
        if s:
            g = G(s["inchikey"])
        else:
            ik = name2ik.get(key) or name2ik.get(clean_name2(key))
            g = groups.get(ik) if ik else None
        if g is None:
            continue
        if g["name"] is None and s:
            g.update({"name": r.key2, "smiles": s["smiles"], "mw": s["mw"], "formula": s["formula"]})
        y = int(r.year) if pd.notna(r.year) else None
        g["nmpa"] = min_year(g["nmpa"], y)
        g["nmpa_approved"] = True
        g["aliases"].add(str(r.key2))
        g["aliases"].add(str(r.en_raw))
        g["nmpa_names"].append(str(r.en_raw))

    # 5. EMA
    ema = ema[(ema.category == "Human") & ema.ma_date.notna()]
    ema = ema[ema.status.isin({"Authorised", "Withdrawn", "Suspended", "Revoked", "Expired", "Lapsed"})]
    for r in ema.itertuples():
        inn = str(r.inn).strip() if pd.notna(r.inn) else ""
        if not inn or ";" in inn:
            continue
        key = clean_name2(inn)
        s = pc.get(key) or pc.get(clean_name(inn))
        g = None
        if s:
            g = G(s["inchikey"])
        else:
            ik = name2ik.get(key) or name2ik.get(clean_name(inn))
            if ik:
                g = groups[ik]
        if g is None:
            continue
        try:
            y = int(str(r.ma_date).split("/")[2])
        except Exception:
            continue
        g["ema"] = min_year(g["ema"], y)
        g["ema_approved"] = True
        g["aliases"].add(inn)
        g["aliases"].add(str(r.name))
        if g["name"] is None:
            g.update({"name": key, "smiles": s["smiles"] if s else None,
                      "mw": s["mw"] if s else None, "formula": s["formula"] if s else None})

    # 6. FDA recent
    for r in fda_recent.itertuples():
        key = clean_name2(str(r.active))
        if not key:
            continue
        s = pc.get(key) or pc.get(clean_name(str(r.active)))
        g = G(s["inchikey"]) if s else (groups.get(name2ik.get(key)) if name2ik.get(key) else None)
        if g is None:
            continue
        if g["name"] is None and s:
            g.update({"name": str(r.active), "smiles": s["smiles"], "mw": s["mw"], "formula": s["formula"]})
        g["fda"] = min_year(g["fda"], year_of(r.date))
        g["fda_approved"] = True
        g["aliases"].add(str(r.active))
        g["aliases"].add(str(r.brand))

    # 7. Merge ester/ketal forms onto an existing parent
    ester_re = re.compile(r"\s+(?:acetate|benzoate|propionate|dipropionate|valerate|butyrate|isobutyrate|enanthate|heptanoate|decanoate|cypionate|undecanoate|phenylpropionate|furoate|palmitate|stearate|laurate|myristate|succinate|fumarate|maleate|tartrate|citrate|phosphate|sulfate|sulphate|mesylate|besylate|tosylate|napsylate|axetil|proxetil|pivoxil|cilexetil|medoxomil|etabonate|bunod|acetonide|pamoate|embonate|tannate)$", re.I)

    def ester_base(n):
        s = str(n).strip()
        s2 = ester_re.sub("", s)
        return clean_name2(s2) if s2 != s else None

    base_map = {}
    for ik, g in list(groups.items()):
        b = ester_base(g["name"]) if g["name"] else None
        if b:
            base_map.setdefault(b, []).append(ik)
    for b, iks in base_map.items():
        base_ik = next((ik for ik, g in groups.items() if g["name"] and str(g["name"]).strip().lower() == b), None)
        if not base_ik:
            continue
        base = groups[base_ik]
        for ik in iks:
            if ik == base_ik:
                continue
            o = groups[ik]
            for c in ("fda", "ema", "nmpa", "pmda"):
                base[c] = min_year(base[c], o[c])
            for c in ("fda_approved", "ema_approved", "nmpa_approved", "pmda_approved"):
                base[c] = base[c] or o[c]
            base["aliases"] |= o["aliases"]
            base["nmpa_names"] += o["nmpa_names"]
            del groups[ik]

    # 8. Rows and filters
    nutrient = {"alanine", "arginine", "aspartic acid", "asparagine", "cysteine", "cystine",
                "glutamic acid", "glutamine", "glycine", "histidine", "isoleucine", "leucine",
                "lysine", "methionine", "phenylalanine", "proline", "serine", "threonine",
                "tryptophan", "tyrosine", "valine", "ornithine", "citrulline", "taurine",
                "xylitol", "sorbitol", "mannitol", "glucose", "fructose", "dextrose", "sucrose",
                "maltose", "lactose", "glycerin", "glycerol", "propylene glycol", "ethanol",
                "glucosamine", "chondroitin"}
    bio_terms = ["gonadotrophin", "gonadotropin", "chorionic", "insulin", "heparin", "interferon",
                 "interleukin", "vaccine", "toxoid", "immunoglobulin", "antibody", "monoclonal",
                 "albumin", "globulin", "plasmin", "streptokinase", "urokinase", "alteplase",
                 "erythropoietin", "somatropin", "follitropin", "menotropin", "urofollitropin",
                 "prolactin", "corticotropin", "sermorelin", "calcitonin", "parathyroid",
                 "somatostatin", "fibrinogen", "thrombin", "collagenase", "hyaluronidase", "glucagon"]
    latin = ["suppositoria", "capsulae", "unguentum", "oeulentum", "tablettae", "ampullae",
             "injectio", "solutio", "suspensio", "emplastrum", "gelatum", "cremor", "pulvis",
             "mixtura", "tinctura", "syrupus", "guttae", "soft"]

    def bad_nmpa_only(name):
        k = name.lower().strip()
        if k in nutrient:
            return True
        if any(x in k for x in latin):
            return True
        return any(x in k for x in bio_terms)

    rows = []
    for ik, g in groups.items():
        if not g["smiles"]:
            continue
        if all(g[c] is None for c in ("fda", "ema", "nmpa", "pmda")) and not any(
                g[c] for c in ("fda_approved", "ema_approved", "nmpa_approved", "pmda_approved")):
            continue
        name = g["name"] or (g["nmpa_names"][0] if g["nmpa_names"] else ik)
        rows.append({"Name": title_name(name), "SMILES": g["smiles"],
                     "FDA": g["fda"], "EMA": g["ema"], "NMPA": g["nmpa"], "PMDA": g["pmda"],
                     "_fda_approved": g["fda_approved"], "_ema_approved": g["ema_approved"],
                     "_nmpa_approved": g["nmpa_approved"], "_pmda_approved": g["pmda_approved"],
                     "_ik": ik, "_aliases": g["aliases"]})
    out = pd.DataFrame(rows)

    # veterinary-only FDA products (no human EMA/NMPA/PMDA approval)
    vet_path = DATA / "dc_vet_only.csv"
    if vet_path.exists():
        vet = pd.read_csv(vet_path)
        vet14 = set(str(x)[:14] for x in vet.inchikey.dropna())
        out["_vet"] = out["_ik"].astype(str).str[:14].isin(vet14)
        vmask = (out["_vet"] & out["EMA"].isna() & out["NMPA"].isna() & out["PMDA"].isna()
                 & ~out["_ema_approved"] & ~out["_nmpa_approved"] & ~out["_pmda_approved"])
        out = out[~vmask].copy()

    mask = out[["FDA", "EMA", "PMDA"]].isna().all(axis=1) & out.Name.map(bad_nmpa_only)
    out = out[~mask].copy()

    # 9. Merge duplicate display names (salt/stereo variants resolved by name)
    out["_nk"] = out.Name.map(lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower()))
    merged = []
    for _, g in out.groupby("_nk", sort=False):
        if len(g) == 1:
            merged.append(g.iloc[0].to_dict())
            continue
        g = g.copy()
        g["_score"] = g[["FDA", "EMA", "NMPA", "PMDA"]].notna().sum(axis=1)
        g = g.sort_values("_score", ascending=False)
        base = g.iloc[0].to_dict()
        aliases = set()
        for a in g["_aliases"]:
            aliases.update(a)
        base["_aliases"] = aliases
        for c in ("FDA", "EMA", "NMPA", "PMDA"):
            vals = [x for x in g[c].tolist() if pd.notna(x)]
            base[c] = min(vals) if vals else None
        for c in ("_fda_approved", "_ema_approved", "_nmpa_approved", "_pmda_approved"):
            base[c] = bool(g[c].any()) if c in g.columns else False
        merged.append(base)
    out = pd.DataFrame(merged).sort_values("Name").reset_index(drop=True)
    return out


# --------------------------------------------------------------------------- #
# ClinicalTrials.gov phases
# --------------------------------------------------------------------------- #
EU = {"Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", "Czech Republic",
      "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Ireland",
      "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland",
      "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden", "Iceland",
      "Liechtenstein", "Norway"}
STATUS_BAD = {"WITHDRAWN", "UNKNOWN", "TERMINATED", "SUSPENDED"}
RANK_LABEL = {6: "Phase 4", 5: "Phase 3", 4: "Phase 2/3", 3: "Phase 2", 2: "Phase 1/2", 1: "Phase 1"}
STOP = {"placebo", "standard", "care", "therapy", "treatment", "drug", "saline", "water",
        "vehicle", "diet", "exercise", "surgery", "radiotherapy", "chemotherapy", "best",
        "supportive", "investigator", "open", "label"}


def norm_alias(s):
    s = str(s).lower()
    s = re.sub(r"\([^)]*\)", " ", s)
    s = re.sub(r"\b\d+(?:\.\d+)?\s*(?:mg|ml|g|mcg|ug|iu|%)\b", " ", s)
    s = re.sub(r"[^a-z0-9]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    t = clean_name(s)
    return t or s


def valid_alias(a):
    return a and len(a) >= 4 and a not in STOP and not a.isdigit()


def phase_rank(phases):
    p = {str(x).upper() for x in phases}
    if {"PHASE2", "PHASE3"} <= p:
        return 4
    if {"PHASE1", "PHASE2"} <= p:
        return 2
    best = 0
    for x in p:
        best = max(best, {"PHASE4": 6, "PHASE3": 5, "PHASE2": 3, "PHASE1": 1, "EARLY_PHASE1": 1}.get(x, 0))
    return best


def match_ngrams(n, alias_map):
    if n in alias_map:
        return {alias_map[n]}
    words = n.split()
    hits = set()
    for i in range(len(words)):
        for j in range(i + 1, min(i + 7, len(words)) + 1):
            d = alias_map.get(" ".join(words[i:j]))
            if d:
                hits.add(d)
    return hits


def ctgov_phases(df: pd.DataFrame) -> pd.DataFrame:
    need = df[df[["FDA", "EMA", "NMPA", "PMDA"]].isna().any(axis=1)].reset_index(drop=True)
    log(f"ClinicalTrials.gov: querying {len(need)} drugs with missing agency data")
    cdir = CACHE / "ct"
    cdir.mkdir(exist_ok=True)
    session = requests.Session()
    session.headers.update({"User-Agent": "drug-approval-collector/1.0"})
    results = {}
    batch_size = 8
    batches = [need.iloc[i:i + batch_size] for i in range(0, len(need), batch_size)]
    for bi, b in enumerate(batches, 1):
        key = hashlib.md5(("|".join(sorted(b["Name"].astype(str))) + f"|{batch_size}").encode()).hexdigest()[:16]
        cp = cdir / f"{key}.json"
        if cp.exists():
            res = json.loads(cp.read_text())
        else:
            names, alias_map = [], {}
            for _, r in b.iterrows():
                prim = norm_alias(r["Name"])
                if valid_alias(prim):
                    names.append(prim)
                for a in r["_aliases"]:
                    na = norm_alias(a)
                    if valid_alias(na) and na not in alias_map:
                        alias_map[na] = r["Name"]
            if not names:
                res = {}
            else:
                q = " OR ".join(f'"{n}"' for n in names[:20])
                params = {
                    "query.intr": q, "pageSize": 500, "countTotal": "true",
                    "fields": "NCTId,Phase,StudyType,OverallStatus,LocationCountry,InterventionName",
                }
                token, page = None, 0
                maxrank, counts = {}, {}
                while True:
                    if token:
                        params["pageToken"] = token
                    rr = None
                    for attempt in range(4):
                        try:
                            rr = session.get("https://clinicaltrials.gov/api/v2/studies", params=params, timeout=60)
                            if rr.status_code == 200:
                                break
                            if rr.status_code in (429, 503):
                                time.sleep(10)
                                continue
                        except Exception:
                            time.sleep(5)
                    if rr is None or rr.status_code != 200:
                        break
                    d = rr.json()
                    studies = d.get("studies", [])
                    token = d.get("nextPageToken")
                    page += 1
                    for st in studies:
                        ps = st.get("protocolSection", {})
                        dm = ps.get("designModule", {})
                        if dm.get("studyType") != "INTERVENTIONAL":
                            continue
                        if ps.get("statusModule", {}).get("overallStatus") in STATUS_BAD:
                            continue
                        rank = phase_rank(dm.get("phases") or [])
                        if rank == 0:
                            continue
                        locs = ps.get("contactsLocationsModule", {}).get("locations", [])
                        countries = {l.get("country") for l in locs if l.get("country")}
                        regions = []
                        if "United States" in countries:
                            regions.append("FDA")
                        if countries & EU:
                            regions.append("EMA")
                        if "China" in countries:
                            regions.append("NMPA")
                        if "Japan" in countries:
                            regions.append("PMDA")
                        if not regions:
                            continue
                        hit = set()
                        for it in ps.get("armsInterventionsModule", {}).get("interventions", []):
                            n = norm_alias(it.get("name", ""))
                            if n:
                                hit |= match_ngrams(n, alias_map)
                        for drug in hit:
                            counts[drug] = counts.get(drug, 0) + 1
                            for reg in regions:
                                maxrank.setdefault(drug, {})[reg] = max(maxrank.get(drug, {}).get(reg, 0), rank)
                    if not token or page >= 12:
                        break
                    time.sleep(1.25)
                res = {name: {"ranks": maxrank.get(name, {}), "n": counts.get(name, 0)} for name in b["Name"]}
            cp.write_text(json.dumps(res))
            time.sleep(1.25)
        results.update(res)
        if bi % 20 == 0:
            log(f"ClinicalTrials.gov: {bi}/{len(batches)} batches")

    out = df.copy()
    for c in ("FDA_phase", "EMA_phase", "NMPA_phase", "PMDA_phase"):
        out[c] = "NA"
    for i, r in out.iterrows():
        x = results.get(r["Name"], {})
        ranks = x.get("ranks", {})
        for c, reg in (("FDA_phase", "FDA"), ("EMA_phase", "EMA"), ("NMPA_phase", "NMPA"), ("PMDA_phase", "PMDA")):
            out.at[i, c] = RANK_LABEL.get(ranks.get(reg, 0), "NA")
    return out


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="drug_approvals_smiles_Deepseek.csv")
    ap.add_argument("--refresh", action="store_true", help="re-download all sources")
    ap.add_argument("--skip-trials", action="store_true", help="do not query ClinicalTrials.gov")
    args = ap.parse_args()

    dc, syn = load_drugcentral(args.refresh)
    ema = load_ema(args.refresh)
    nmpa = load_nmpa(args.refresh)
    fda_recent = load_fda_recent(args.refresh)

    # names needing a PubChem lookup
    names = []
    names += [r.name for r in dc[dc.smiles.isna()].itertuples()]
    names += [k for k in nmpa.key2.dropna().unique()]
    # EMA INNs not already covered by exact name/synonym
    covered = set()
    for x in dc.name.dropna():
        covered.add(clean_name2(x))
    for x in syn.name.dropna():
        covered.add(clean_name2(x))
    for x in syn.struct_name.dropna():
        covered.add(clean_name2(x))
    for x in nmpa.key2.dropna():
        covered.add(x)
    for r in ema[(ema.category == "Human") & ema.ma_date.notna()].itertuples():
        inn = str(r.inn).strip() if pd.notna(r.inn) else ""
        if inn and ";" not in inn and clean_name2(inn) not in covered:
            names.append(clean_name2(inn))
    # FDA recent names not covered
    for x in fda_recent.active.dropna():
        k = clean_name2(x)
        if k and k not in covered:
            names.append(k)
    pc = pubchem_resolve(list(dict.fromkeys(names)))

    out = build_table(dc, syn, nmpa, ema, fda_recent, pc)
    log(f"approval table: {len(out)} rows")

    if not args.skip_trials:
        out = ctgov_phases(out)

    out = out.sort_values("Name").reset_index(drop=True)

    def fmt(year, approved, phase):
        if pd.notna(year):
            try:
                return str(int(year))
            except Exception:
                pass
        if bool(approved):
            return "Approved"
        if isinstance(phase, str) and phase.startswith("Phase"):
            return phase
        return "NA"

    final = pd.DataFrame({
        "Name": out["Name"],
        "SMILES": out["SMILES"],
        "FDA": [fmt(out.at[i, "FDA"], out.at[i, "_fda_approved"], out.at[i, "FDA_phase"] if "FDA_phase" in out.columns else None) for i in out.index],
        "EMA": [fmt(out.at[i, "EMA"], out.at[i, "_ema_approved"], out.at[i, "EMA_phase"] if "EMA_phase" in out.columns else None) for i in out.index],
        "NMPA": [fmt(out.at[i, "NMPA"], out.at[i, "_nmpa_approved"], out.at[i, "NMPA_phase"] if "NMPA_phase" in out.columns else None) for i in out.index],
        "PMDA": [fmt(out.at[i, "PMDA"], out.at[i, "_pmda_approved"], out.at[i, "PMDA_phase"] if "PMDA_phase" in out.columns else None) for i in out.index],
    })
    final.to_csv(args.out, index=False)
    log(f"wrote {args.out} ({len(final)} rows)")

    counts = {c: (final[c] != "NA").sum() for c in ("FDA", "EMA", "NMPA", "PMDA")}
    log(f"non-NA cells: {counts}")


if __name__ == "__main__":
    main()
