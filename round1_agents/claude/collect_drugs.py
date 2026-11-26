#!/usr/bin/env python3
"""
Collect small-molecule drugs (approved or in clinical trials) with SMILES and their
status at four regulators: FDA (US), EMA (EU), NMPA (China), PMDA (Japan).

Output cell values per agency (highest evidence wins):
    <year>                      approved; year of first approval found for that agency
    Approved                    approved/marketed there, but no reliable year in the source
    Phase 4 / Phase 3 / ...     highest clinical-trial phase with a trial site in that region
    NA                          no record found

Sources (all public, downloaded and cached under data/raw/):
    Structures/universe : ChEMBL REST API (molecule_type = "Small molecule", any clinical phase)
    FDA approvals       : Drugs@FDA via openFDA bulk export (earliest original NDA/ANDA approval)
    EMA approvals       : KEGG BRITE br08329 (EU new active substances, 1995-, gives the year)
                          + EMA "Medicines output" report (all centralised products, 'Approved')
                          + EMA Article 57 database (national EU/EEA authorisations, 'Approved')
    PMDA approvals      : KEGG BRITE br08318 (Japan new drug approvals, 2000-) + br08301
                          (all drugs marketed in Japan, no date)
    NMPA approvals      : NMPA English site approval news (novel drugs, ~2021-), plus an
                          optional manual file data/manual/nmpa_approvals.csv (name,year)
    Clinical phases     : ClinicalTrials.gov API v2 (trial site countries -> region),
                          linked to molecules via ChEMBL drug_indication NCT refs + name matching

Usage:
    python3 collect_drugs.py                # uses cache when present
    python3 collect_drugs.py --refresh      # re-download everything
Requires: requests, pandas, openpyxl, lxml
"""
import argparse
import gzip
import io
import json
import re
import sys
import time
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]   # repo root (shared data/ folder); was .parent in round 1
RAW = ROOT / "data" / "raw"
MANUAL = ROOT / "data" / "manual"
OUT = Path(__file__).resolve().parent
UA = {"User-Agent": "Mozilla/5.0"}
REFRESH = False

AGENCIES = ["FDA", "EMA", "NMPA", "PMDA"]

EEA = {
    "Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", "Czech Republic",
    "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Ireland",
    "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland",
    "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden",
    "Iceland", "Liechtenstein", "Norway",
}
REGION_OF_COUNTRY = {"United States": "FDA", "Puerto Rico": "FDA", "China": "NMPA", "Japan": "PMDA"}
REGION_OF_COUNTRY.update({c: "EMA" for c in EEA})

PHASE_SCORE = {"EARLY_PHASE1": 0.5, "PHASE1": 1, "PHASE2": 2, "PHASE3": 3, "PHASE4": 4}
PHASE_LABEL = {0.5: "Early Phase 1", 1: "Phase 1", 1.5: "Phase 1/2", 2: "Phase 2",
               2.5: "Phase 2/3", 3: "Phase 3", 4: "Phase 4"}


# ----------------------------------------------------------------------------- helpers
def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def http_get(url, params=None, retries=5, timeout=120, **kw):
    for i in range(retries):
        try:
            r = requests.get(url, params=params, headers=UA, timeout=timeout, **kw)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                r.raise_for_status()
        except requests.RequestException as e:
            if i == retries - 1:
                raise
            log(f"  retry {i + 1} {url}: {e}")
        time.sleep(2 * (i + 1))
    raise RuntimeError(f"GET failed: {url}")


def cached_download(name, url):
    path = RAW / name
    if REFRESH or not path.exists():
        log(f"download {url}")
        path.write_bytes(http_get(url).content)
    return path


def cached_json(name, builder):
    path = RAW / name
    if not REFRESH and path.exists():
        with gzip.open(path, "rt") as f:
            return json.load(f)
    data = builder()
    with gzip.open(path, "wt") as f:
        json.dump(data, f)
    return data


# ----------------------------------------------------------------------------- names
SALT_WORDS = set("""
hydrochloride dihydrochloride trihydrochloride hcl hydrobromide hbr hydroiodide
sodium disodium trisodium potassium dipotassium calcium magnesium lithium zinc ammonium
mesylate mesilate dimesylate besylate besilate tosylate tosilate esylate edisylate
maleate dimaleate fumarate hemifumarate tartrate bitartrate hemitartrate citrate dihydrogen
succinate hemisuccinate sulfate sulphate bisulfate hemisulfate phosphate diphosphate
acetate diacetate bromide chloride iodide nitrate lactate gluconate malate oxalate
pamoate embonate palmitate stearate hyclate napsylate napadisylate xinafoate
trometamol tromethamine meglumine lysine arginine olamine diolamine erbumine choline
hydrate monohydrate dihydrate trihydrate tetrahydrate pentahydrate sesquihydrate
hemihydrate anhydrous monosodium monopotassium salt free base acid
""".split())
# 'acid' only stripped as last resort (handled in candidates)
SALT_WORDS.discard("acid")
# counter-ion acids written as "<drug> fumaric acid" etc.
SALT_ACIDS = set("""fumaric maleic succinic tartaric citric hydrochloric hydrobromic phosphoric
sulfuric sulphuric acetic methanesulfonic benzenesulfonic toluenesulfonic lactic malic oxalic
adipic benzoic edisilic""".split())
FORM_WORDS = set("""
tablet tablets capsule capsules injection injectable infusion solution suspension
oral iv im sc intravenous intramuscular subcutaneous topical cream gel ointment patch
film coated coatedtablets enteric extended release sustained modified delayed dispersible
chewable granules granule powder for concentrate syrup drops spray inhalation inhaler
nasal ophthalmic eye ear lotion suppository suppositories emulsion lyophilized
sterile vial dose doses dosage formulation liquid prefilled syringe pen co packaged
copackaged kit and
""".split())
DOSE_RE = re.compile(r"\b\d+(?:[.,]\d+)?\s*(?:mg|g|mcg|ug|µg|ml|iu|units?|%|mg/kg|mg/m2|mg/ml|mmol)\b", re.I)


def basic_norm(s):
    s = str(s).lower()
    s = s.replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    s = re.sub(r"\([^)]*\)|\[[^\]]*\]", " ", s)          # drop bracketed annotations
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return " ".join(s.split())


def key_of(words):
    return "".join(words)


def candidate_keys(name):
    """Normalised lookup keys for a name, from most to least specific."""
    base = basic_norm(DOSE_RE.sub(" ", str(name)))
    if not base:
        return []
    words = base.split()
    keys = [key_of(words)]
    # strip trailing dosage-form words, then trailing salt/hydrate words
    w = list(words)
    while len(w) > 1 and w[-1] in FORM_WORDS:
        w.pop()
    keys.append(key_of(w))
    w2 = list(w)
    while len(w2) > 1 and (w2[-1] in SALT_WORDS or
                           (w2[-1] == "acid" and len(w2) > 2 and w2[-2] in SALT_ACIDS)):
        w2 = w2[:-2] if w2[-1] == "acid" else w2[:-1]
    keys.append(key_of(w2))
    # leading cation: "sodium valproate" -> "valproate"
    w3 = list(w2)
    while len(w3) > 1 and w3[0] in SALT_WORDS:
        w3.pop(0)
    keys.append(key_of(w3))
    # "xxx acid" <-> "xxx"
    if len(w3) > 1 and w3[-1] == "acid":
        keys.append(key_of(w3[:-1]))
    out = []
    for k in keys:
        if k and len(k) >= 3 and k not in out:
            out.append(k)
    return out


class NameIndex:
    """Two tiers: preferred names (ChEMBL pref_name of parent/salt) beat synonyms."""

    def __init__(self):
        self.pref = defaultdict(set)
        self.idx = defaultdict(set)

    def add(self, name, cid, preferred=False):
        if not name:
            return
        k = key_of(basic_norm(name).split())
        if len(k) >= 3:
            self.idx[k].add(cid)
            if preferred:
                self.pref[k].add(cid)

    def lookup(self, name):
        for k in candidate_keys(name):
            for tier in (self.pref, self.idx):
                hit = tier.get(k)
                if hit and len(hit) == 1:
                    return next(iter(hit))
            if self.idx.get(k):  # ambiguous in both tiers -> do not guess
                return None
        return None


def split_ingredients(s):
    parts = re.split(r";|\+|/| plus | and | with |,|&", str(s), flags=re.I)
    return [p.strip() for p in parts if p.strip()]


# ----------------------------------------------------------------------------- ChEMBL
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"


def chembl_paged(endpoint, params, key):
    out, offset, limit = [], 0, 500
    while True:
        p = dict(params, limit=limit, offset=offset, format="json")
        d = http_get(f"{CHEMBL}/{endpoint}.json", params=p).json()
        out.extend(d[key])
        total = d["page_meta"]["total_count"]
        offset += limit
        log(f"  ChEMBL {endpoint}: {min(offset, total)}/{total}")
        if offset >= total:
            return out


def slim_mol(m):
    st = m.get("molecule_structures") or {}
    return {
        "id": m["molecule_chembl_id"],
        "name": m.get("pref_name"),
        "type": m.get("molecule_type"),
        "max_phase": m.get("max_phase"),
        "parent": (m.get("molecule_hierarchy") or {}).get("parent_chembl_id") or m["molecule_chembl_id"],
        "smiles": st.get("canonical_smiles"),
        "inchikey": st.get("standard_inchi_key"),
        "mw": (m.get("molecule_properties") or {}).get("full_mwt"),
        "syn": sorted({s.get("molecule_synonym") or s.get("synonyms") for s in (m.get("molecule_synonyms") or [])} - {None}),
    }


def load_chembl():
    def build():
        mols = chembl_paged("molecule", {"molecule_type": "Small molecule", "max_phase__isnull": "false"}, "molecules")
        mols = [slim_mol(m) for m in mols]
        have = {m["id"] for m in mols}
        missing = sorted({m["parent"] for m in mols} - have)
        for i in range(0, len(missing), 50):
            chunk = ";".join(missing[i:i + 50])
            d = http_get(f"{CHEMBL}/molecule/set/{chunk}.json").json()
            mols.extend(slim_mol(m) for m in d["molecules"])
        return mols

    mols = cached_json("chembl_molecules.json.gz", build)
    by_id = {m["id"]: m for m in mols}

    def eff_parent(m):
        # merge only salts/hydrates/mixtures ('.' in SMILES) into the ChEMBL parent; keep other
        # children (e.g. deuterated drugs, which ChEMBL's parent drops the isotopes of) separate
        if m["parent"] in by_id and (m["id"] == m["parent"] or "." in (m.get("smiles") or "") or not m.get("smiles")):
            return by_id[m["parent"]]
        return m

    parents = {}
    for m in mols:
        p = eff_parent(m)
        if p.get("type") != "Small molecule" or not p.get("smiles"):
            continue
        parents.setdefault(p["id"], p)
    index = NameIndex()
    child_of = {}
    for m in mols:
        pid = eff_parent(m)["id"]
        if pid not in parents:
            continue
        child_of[m["id"]] = pid
        index.add(m["name"], pid, preferred=True)
        for s in m["syn"]:
            index.add(s, pid)
    log(f"ChEMBL: {len(mols)} records -> {len(parents)} small-molecule parents, {len(index.idx)} name keys")
    return parents, child_of, index


def load_chembl_trial_links():
    def build():
        rows = chembl_paged("drug_indication", {}, "drug_indications")
        links = defaultdict(set)
        for r in rows:
            mid = r.get("molecule_chembl_id") or r.get("parent_molecule_chembl_id")
            for ref in r.get("indication_refs") or []:
                if ref.get("ref_type") == "ClinicalTrials":
                    for n in re.findall(r"NCT\d{8}", ref.get("ref_id") or ""):
                        links[mid].add(n)
        return {k: sorted(v) for k, v in links.items()}

    return cached_json("chembl_trial_links.json.gz", build)


# ----------------------------------------------------------------------------- FDA
def load_fda(index):
    # openFDA bulk export of Drugs@FDA (FDA's official endpoint for programmatic downloads)
    z = zipfile.ZipFile(cached_download(
        "openfda_drugsfda.json.zip",
        "https://download.open.fda.gov/drug/drugsfda/drug-drugsfda-0001-of-0001.json.zip"))
    apps = json.load(z.open(z.namelist()[0]))["results"]
    res, unmatched = {}, set()
    for app in apps:
        years = [int(s["submission_status_date"][:4]) for s in app.get("submissions") or []
                 if s.get("submission_type") == "ORIG" and s.get("submission_status") == "AP"
                 and str(s.get("submission_status_date", ""))[:4].isdigit()]
        if not years:
            continue
        y = min(years)
        y = y if y >= 1938 else "Approved"
        names = {i.get("name") for p in app.get("products") or [] for i in p.get("active_ingredients") or []}
        for ing in filter(None, names):
            for part in ([ing] if index.lookup(ing) else split_ingredients(ing)):
                cid = index.lookup(part)
                if cid:
                    res[cid] = best_year(res.get(cid), y)
                else:
                    unmatched.add(part)
    log(f"FDA: {len(res)} molecules matched ({len(unmatched)} unmatched ingredient names, mostly biologics/mixtures)")
    return res


def best_year(a, b):
    """Combine two approval values: earliest year wins; a year beats 'Approved'."""
    if a is None:
        return b
    if isinstance(a, int) and isinstance(b, int):
        return min(a, b)
    return a if isinstance(a, int) else b


# ----------------------------------------------------------------------------- EMA
def read_ema_excel(path, must_have):
    raw = pd.read_excel(path, header=None, dtype=str)
    hdr = next(i for i in range(60) if any(must_have in str(v) for v in raw.iloc[i].values))
    df = raw.iloc[hdr + 1:].copy()
    df.columns = [str(c).split("\n")[0].strip() for c in raw.iloc[hdr].values]
    return df


def load_ema(index, fda=None):
    # first-approval year of new active substances (EU centralised procedure, 1995-)
    res = kegg_new_approvals("br08329", index, fda)
    n_year = sum(isinstance(v, int) for v in res.values())
    kegg_ids, fallback = set(res), {}
    # any centrally authorised human medicine (incl. generics/combos) -> at least 'Approved'
    df = read_ema_excel(cached_download("ema_medicines.xlsx",
                        "https://www.ema.europa.eu/en/documents/report/medicines-output-medicines-report_en.xlsx"),
                        "Name of medicine")
    df = df[(df["Category"].str.strip().str.lower() == "human") & df["Marketing authorisation date"].notna()]
    for _, r in df.iterrows():
        names = set()
        for col in ["International non-proprietary name (INN) / common name", "Active substance"]:
            if pd.notna(r.get(col)):
                names.update(split_ingredients(r[col]))
        single = len(names) == 1 and str(r.get("Generic", "")).strip().lower() != "yes"
        ma = pd.to_datetime(r["Marketing authorisation date"], errors="coerce", dayfirst=True)
        for n in names:
            cid = index.lookup(n)
            if not cid:
                continue
            # fallback year for new substances missing from the KEGG table: single-ingredient,
            # non-generic product, and not a drug already approved by FDA >3 years earlier
            f = (fda or {}).get(cid)
            if single and pd.notna(ma) and cid not in kegg_ids and not (isinstance(f, int) and f <= ma.year - 3):
                fallback[cid] = min(fallback.get(cid, 9999), int(ma.year))
            res[cid] = best_year(res.get(cid), "Approved")
    for cid, y in fallback.items():
        res[cid] = best_year(res.get(cid), y)
    n_central = len(res)
    # national EU/EEA authorisations (no dates)
    a57 = read_ema_excel(cached_download("ema_article57.xlsx",
                         "https://www.ema.europa.eu/en/documents/other/article-57-product-data_en.xlsx"),
                         "Product authorisation country")
    subst = a57["Active substance"].dropna().unique()
    national = 0
    for s in subst:
        for n in split_ingredients(s.replace("|", ";")):
            cid = index.lookup(n)
            if cid and cid not in res:
                res[cid] = "Approved"
                national += 1
    log(f"EMA: {n_year} new active substances with year, {n_central - n_year} other central, +{national} national-only")
    return res


# ----------------------------------------------------------------------------- PMDA (via KEGG)
def kegg_names(cell):
    return [n.strip() for n in re.split(r";", re.sub(r"\([^)]*\)", "", str(cell))) if n.strip()]


def kegg_new_approvals(br, index, fda=None):
    """KEGG 'New drug approvals' table -> {chembl_id: first-approval year or 'Approved'}.
    Year only for rows flagged NME (new active substance); other rows are new
    formulations/indications of already-approved drugs."""
    html = cached_download(f"kegg_{br}.html", f"https://rest.kegg.jp/get/br:{br}").read_text()
    t = pd.read_html(io.StringIO(html))[0]
    t = t.iloc[1:].rename(columns={0: "Date", 3: "Ingredient"} if br == "br08329" else {0: "Date", 4: "Ingredient"})
    res, combo = {}, {}
    for _, r in t.iterrows():
        nme = any(str(v).strip() in ("NME", "EUA") for v in r.values)   # EUA = emergency/conditional approval
        y = int(str(r["Date"])[:4]) if nme else "Approved"
        ing = str(r["Ingredient"])
        cid = index.lookup(ing)
        if cid:
            res[cid] = best_year(res.get(cid), y)
            continue
        # combination row: the NME flag only applies to the new component(s)
        for part in split_ingredients(ing):
            cid = index.lookup(part)
            if cid:
                combo[cid] = best_year(combo.get(cid), y)
    for cid, y in combo.items():
        old_elsewhere = isinstance(y, int) and isinstance((fda or {}).get(cid), int) and fda[cid] <= y - 3
        if cid in res or old_elsewhere:
            res[cid] = best_year(res.get(cid), "Approved")
        else:
            res[cid] = y
    return res


def load_pmda(index, fda=None):
    res = kegg_new_approvals("br08318", index, fda)
    n_new = len(res)
    txt = cached_download("kegg_br08301.txt", "https://rest.kegg.jp/get/br:br08301").read_text()
    for m in re.finditer(r"^E\s+D\d{5}\s+(.+)$", txt, re.M):
        for n in kegg_names(m.group(1)):
            cid = index.lookup(n)
            if cid and cid not in res:
                res[cid] = "Approved"
    # NME year must beat 'Approved' from the marketed list; best_year already ensures that
    log(f"PMDA: {n_new} molecules from new-approval table (2000-), {len(res) - n_new} more marketed in Japan")
    return res


# ----------------------------------------------------------------------------- NMPA
MONTHS = "January|February|March|April|May|June|July|August|September|October|November|December"


def load_nmpa(index, fda=None):
    def build():
        base = "https://english.nmpa.gov.cn"
        items = {}
        for sec in ["news", "drugs", "regulatoryinformation"]:
            first = http_get(f"{base}/{sec}.html").text
            n = max([int(x) for x in re.findall(rf"{sec}_(\d+)\.html", first)] + [1])
            for i in range(1, n + 1):
                t = first if i == 1 else http_get(f"{base}/{sec}_{i}.html").text
                for href, title in re.findall(r'<a [^>]*href="([^"]*c_\d+\.htm)"[^>]*>([^<]+)', t):
                    items[href] = " ".join(title.split())
            log(f"  NMPA {sec}: {n} pages")
        out = []
        for href, title in items.items():
            if not re.search(r"approved", title, re.I) or re.search(
                    r"announcement|policy|interpretation|procedure|provisions|guideline|device|system|"
                    r"stent|software|kit|valve|machine|prosthesis|catheter|instrument", title, re.I):
                continue
            year = int(href[:4])
            try:
                body = http_get(f"{base}/{href}").text
                m = re.search(rf"(?:{MONTHS})\s+\d{{1,2}},\s*(20\d\d)", re.sub(r"<[^>]+>", " ", body))
                if m:
                    year = int(m.group(1))
            except Exception:
                pass
            out.append({"title": title, "year": year, "url": f"{base}/{href}"})
        return out

    news = cached_json("nmpa_news.json.gz", build)
    res = {}
    for it in news:
        prod = re.split(r"\s+(?:approved|is approved|gets approved)", it["title"], flags=re.I)[0]
        prod = re.sub(r"^(innovative( tcm)?|importation of|china'?s first|first)\s+", "", prod, flags=re.I)
        prod = re.sub(r"^[A-Za-z]+'s\s+", "", prod)
        parts = split_ingredients(prod)
        for part in parts:
            cid = index.lookup(part)
            if cid:
                f = (fda or {}).get(cid)
                old_drug = len(parts) > 1 and isinstance(f, int) and f <= it["year"] - 3
                res[cid] = best_year(res.get(cid), "Approved" if old_drug else it["year"])
    n_news = len(res)
    manual = MANUAL / "nmpa_approvals.csv"
    if manual.exists():
        m = pd.read_csv(manual, dtype=str)
        for _, r in m.iterrows():
            cid = index.lookup(r["name"])
            if cid:
                y = int(r["year"]) if str(r.get("year", "")).strip().isdigit() else "Approved"
                res[cid] = best_year(res.get(cid), y)
    log(f"NMPA: {n_news} molecules from NMPA approval news, {len(res) - n_news} more from manual file")
    return res


def apply_manual(agency, res, index):
    """Optional extra approvals: data/manual/<agency>_approvals.csv with columns name,year."""
    f = MANUAL / f"{agency.lower()}_approvals.csv"
    if agency == "NMPA" or not f.exists():
        return
    for _, r in pd.read_csv(f, dtype=str).iterrows():
        cid = index.lookup(r["name"])
        if cid:
            y = int(r["year"]) if str(r.get("year", "")).strip().isdigit() else "Approved"
            res[cid] = best_year(res.get(cid), y)


# ----------------------------------------------------------------------------- ClinicalTrials.gov
def load_trials():
    def build():
        url = "https://clinicaltrials.gov/api/v2/studies"
        phases = " OR ".join(f"AREA[Phase]{p}" for p in PHASE_SCORE)
        params = {
            "pageSize": 1000,
            "fields": "NCTId,Phase,OverallStatus,LocationCountry,InterventionType,InterventionName,InterventionOtherName",
            "filter.advanced": f"AREA[StudyType]INTERVENTIONAL AND ({phases})",
        }
        out, token, page = [], None, 0
        while True:
            p = dict(params, pageToken=token) if token else params
            d = http_get(url, params=p).json()
            for s in d.get("studies", []):
                ps = s.get("protocolSection", {})
                ph = (ps.get("designModule") or {}).get("phases") or []
                locs = (ps.get("contactsLocationsModule") or {}).get("locations") or []
                ivs = (ps.get("armsInterventionsModule") or {}).get("interventions") or []
                out.append({
                    "nct": ps["identificationModule"]["nctId"],
                    "phases": ph,
                    "status": (ps.get("statusModule") or {}).get("overallStatus"),
                    "countries": sorted({l.get("country") for l in locs if l.get("country")}),
                    "iv": [[i.get("type"), i.get("name"), i.get("otherNames") or []] for i in ivs],
                })
            page += 1
            if page % 20 == 0:
                log(f"  ClinicalTrials.gov: {len(out)} studies")
            token = d.get("nextPageToken")
            if not token:
                return out

    return cached_json("ctgov_studies.json.gz", build)


def trial_phase(phases):
    sc = sorted(PHASE_SCORE[p] for p in phases if p in PHASE_SCORE)
    if not sc:
        return None
    if len(sc) == 2:  # PHASE1/PHASE2 or PHASE2/PHASE3
        return (sc[0] + sc[1]) / 2 if sc[1] - sc[0] == 1 else sc[-1]
    return sc[-1]


def trial_phases_by_region(studies, chembl_links, child_of, index):
    nct_to_mols = defaultdict(set)
    for mid, ncts in chembl_links.items():
        pid = child_of.get(mid)
        if pid:
            for n in ncts:
                nct_to_mols[n].add(pid)
    ok_types = {"DRUG", "COMBINATION_PRODUCT", "DIETARY_SUPPLEMENT", "BIOLOGICAL"}
    cache = {}
    res = defaultdict(dict)  # cid -> agency -> phase score
    for s in studies:
        if s["status"] == "WITHDRAWN":
            continue
        ph = trial_phase(s["phases"])
        regions = {REGION_OF_COUNTRY[c] for c in s["countries"] if c in REGION_OF_COUNTRY}
        if ph is None or not regions:
            continue
        mols = set(nct_to_mols.get(s["nct"], ()))
        for typ, name, others in s["iv"]:
            if typ not in ok_types:
                continue
            for nm in [name] + list(others):
                for part in split_ingredients(nm) + [nm]:
                    if part not in cache:
                        cache[part] = index.lookup(part)
                    if cache[part]:
                        mols.add(cache[part])
        for cid in mols:
            for reg in regions:
                if ph > res[cid].get(reg, 0):
                    res[cid][reg] = ph
    log(f"Trials: {len(res)} molecules with at least one trial in US/EU/China/Japan")
    return res


# ----------------------------------------------------------------------------- assemble
def cell(approval, phase):
    if isinstance(approval, int):
        return str(approval)
    if approval == "Approved":
        return "Approved"
    if phase:
        return PHASE_LABEL.get(phase, f"Phase {phase}")
    return "NA"


def main():
    global REFRESH
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true", help="re-download all sources")
    args = ap.parse_args()
    REFRESH = args.refresh
    for d in (RAW, MANUAL, OUT):
        d.mkdir(parents=True, exist_ok=True)

    parents, child_of, index = load_chembl()
    approvals = {"FDA": load_fda(index)}
    approvals["EMA"] = load_ema(index, approvals["FDA"])
    approvals["NMPA"] = load_nmpa(index, approvals["FDA"])
    approvals["PMDA"] = load_pmda(index, approvals["FDA"])
    for a in AGENCIES:
        apply_manual(a, approvals[a], index)
    phases = trial_phases_by_region(load_trials(), load_chembl_trial_links(), child_of, index)

    rows = []
    for cid, p in parents.items():
        vals = {a: cell(approvals[a].get(cid), phases.get(cid, {}).get(a)) for a in AGENCIES}
        if all(v == "NA" for v in vals.values()):
            continue
        rows.append({"Name": (p["name"] or cid).title() if p["name"] else cid, "SMILES": p["smiles"], **vals,
                     "ChEMBL_ID": cid, "InChIKey": p["inchikey"], "MW": p["mw"],
                     "ChEMBL_max_phase": p["max_phase"]})
    df = pd.DataFrame(rows).sort_values("Name", key=lambda s: s.str.lower())
    df[["Name", "SMILES"] + AGENCIES].to_csv(OUT / "small_molecule_drugs_by_agency_Claude.csv", index=False)
    df.to_csv(OUT / "small_molecule_drugs_by_agency_detailed_Claude.csv", index=False)

    log(f"Wrote {len(df)} molecules -> {OUT / 'small_molecule_drugs_by_agency_Claude.csv'}")
    for a in AGENCIES:
        v = df[a]
        log(f"  {a}: approved w/ year {v.str.fullmatch(r'\d{4}').sum()}, approved no year "
            f"{(v == 'Approved').sum()}, clinical only {v.str.startswith(('Phase', 'Early')).sum()}, NA {(v == 'NA').sum()}")


if __name__ == "__main__":
    sys.exit(main())
