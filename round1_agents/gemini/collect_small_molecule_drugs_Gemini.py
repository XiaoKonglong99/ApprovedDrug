#!/usr/bin/env python3
"""
Collect small-molecule drugs (excluding proteins, antibodies, RNA/oligonucleotides,
vaccines, enzymes, cell/gene therapies, etc.) and their canonical SMILES structures,
along with regulatory approval years or highest clinical trial stages across:
  - FDA  (United States)
  - EMA  (European Union / EEA)
  - NMPA (China)
  - PMDA (Japan)

Data Sources (all public, open APIs / official bulk files):
  1. ChEMBL REST API (molecule_type="Small molecule", max_phase >= 1):
     Master small-molecule catalog, canonical SMILES, parent-salt hierarchy, synonyms.
  2. FDA Drugs@FDA (openFDA bulk JSON ZIP) + DrugCentral FDA_Approved.csv:
     Official US FDA approval years.
  3. EMA Medicines Report (official EMA JSON) + DrugCentral EMA_Approved.csv:
     Official EU EMA marketing authorisation years.
  4. PMDA List of Approved Drugs (official PMDA April 2004 - Feb 2026 PDF +
     PMDA Review Reports HTML) + DrugCentral PMDA_Approved.csv:
     Official Japan PMDA approval years.
  5. ClinicalTrials.gov API v2 (interventional DRUG studies, Phase 1-4):
     Regional clinical trial stages (Stage 1 / Stage 2 / Stage 3) for US (FDA),
     EU/EEA (EMA), China (NMPA), and Japan (PMDA), plus Phase 4 post-marketing
     signals for NMPA (China) where open bulk NMPA approval-year feeds are unavailable.

Output CSV columns:
  Name,SMILES,FDA,EMA,NMPA,PMDA

Usage:
  python3 collect_small_molecule_drugs_Gemini.py --output small_molecule_drugs_Gemini.csv
"""

import argparse
import concurrent.futures
import csv
import io
import json
import os
import re
import ssl
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict

# ---------------------------------------------------------------------------
# Configuration & Constants
# ---------------------------------------------------------------------------

SSL_CTX = ssl.create_default_context()
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) SmallMoleculeDrugCollector/1.0",
    "Accept": "application/json, text/plain, */*",
}

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
OPENFDA_BULK_URL = "https://download.open.fda.gov/drug/drugsfda/drug-drugsfda-0001-of-0001.json.zip"
EMA_JSON_URL = "https://www.ema.europa.eu/en/documents/report/medicines-output-medicines_json-report_en.json"
PMDA_HTML_URL = "https://www.pmda.go.jp/english/review-services/reviews/approved-information/drugs/0001.html"
PMDA_PDF_URL = "https://www.pmda.go.jp/files/000281190.pdf"
DRUGCENTRAL_FDA_URL = "https://drugcentral.org/static/FDA_Approved.csv"
DRUGCENTRAL_EMA_URL = "https://drugcentral.org/static/EMA_Approved.csv"
DRUGCENTRAL_PMDA_URL = "https://drugcentral.org/static/PMDA_Approved.csv"
CTGOV_API_URL = "https://clinicaltrials.gov/api/v2/studies"

EU_EEA_COUNTRIES = {
    "Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia",
    "Czech Republic", "Denmark", "Estonia", "Finland", "France", "Germany",
    "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Latvia",
    "Liechtenstein", "Lithuania", "Luxembourg", "Malta", "Netherlands",
    "Norway", "Poland", "Portugal", "Romania", "Slovakia", "Slovenia",
    "Spain", "Sweden",
}

# Common pharmaceutical salt / hydrate / ester / counterion suffixes to strip
# during drug-name normalization so salt forms map to their parent small molecule.
SALT_WORDS = {
    "hydrochloride", "hcl", "dihydrochloride", "trihydrochloride", "monohydrochloride",
    "hydrobromide", "hbr", "sodium", "disodium", "trisodium", "monosodium",
    "potassium", "dipotassium", "calcium", "magnesium", "zinc", "aluminum", "aluminium",
    "lithium", "ammonium", "meglumine", "tromethamine", "trometamol", "olamine",
    "diolamine", "trolamine", "ethanolamine", "diethanolamine", "triethanolamine",
    "mesylate", "mesilate", "dimesylate", "dimesilate", "besylate", "besilate",
    "tosylate", "tosilate", "ditosylate", "esylate", "esilate", "edisylate", "edisilate",
    "maleate", "dimaleate", "fumarate", "hemifumarate", "difumarate",
    "tartrate", "bitartrate", "hemitartrate", "l-tartrate",
    "citrate", "succinate", "hemisuccinate", "oxalate", "malate", "lactate",
    "acetate", "diacetate", "triacetate", "propionate", "dipropionate",
    "butyrate", "valerate", "pivalate", "benzoate", "salicylate",
    "sulfate", "sulphate", "bisulfate", "hemisulfate", "hydrogen",
    "phosphate", "diphosphate", "dihydrogen", "monohydrogen",
    "nitrate", "chloride", "bromide", "iodide", "fluoride", "carbonate", "bicarbonate",
    "gluconate", "glucoheptonate", "lactobionate", "pamoate", "embonate",
    "xinafoate", "napadisylate", "napadisilate", "napsylate", "napsilate",
    "isethionate", "isetionate", "estolate", "stearate", "palmitate",
    "decanoate", "undecanoate", "enanthate", "cypionate", "caproate",
    "furoate", "acetonide", "hexacetonide", "axetil", "proxetil", "pivoxil",
    "medoxomil", "cilexetil", "alafenamide", "disoproxil", "fosamil",
    "hydrate", "monohydrate", "dihydrate", "trihydrate", "tetrahydrate",
    "pentahydrate", "hexahydrate", "hemihydrate", "sesquihydrate", "anhydrous",
    "ethanolate", "solvate", "triflate", "trifluoroacetate", "adipate",
    "aspartate", "glutamate", "glycinate", "lysinate", "arginate",
    "camphorsulfonate", "camsylate", "camsilate", "cyclamate", "saccharinate",
    "hippurate", "mandelate", "nicotinate", "orotate", "pantothenate",
    "subsalicylate", "subcitrate", "subgallate", "tannate", "teoclate",
    "bromohydrate", "chlorhydrate", "base", "free", "acid", "salt",
    "injection", "tablets", "tablet", "capsules", "capsule", "oral", "solution",
    "suspension", "cream", "ointment", "gel", "patch", "infusion", "powder",
    "iv", "po", "sc", "im", "topical", "ophthalmic", "nasal", "inhalation",
}

# Patterns indicating non-small-molecule biologics in external agency feeds
BIOLOGIC_SUFFIX_RE = re.compile(
    r"(mab|cept|kinra|fusp|tug|bart|mig|ment|pafusp|parib_no|tide|relin|pressin|"
    r"tropin|poetin|grastim|plostim|feron|leukin|kin|ase|zyme|virsen|siran|"
    r"nersen|tarsen|mersen|gen|vec|cel|cogene|vaccine|toxoid|immunoglobulin|"
    r"antitoxin|insulin|glucagon|somatropin|interferon|pegfilgrastim)\b",
    re.I,
)


# ---------------------------------------------------------------------------
# HTTP & Caching Helpers
# ---------------------------------------------------------------------------

def http_get_bytes(url, retries=4, timeout=45):
    """Fetch URL bytes with retry/backoff."""
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=timeout, context=SSL_CTX) as resp:
                return resp.read()
        except Exception as e:
            last_err = e
            wait = 1.5 * (attempt + 1)
            time.sleep(wait)
    raise RuntimeError(f"Failed to fetch {url} after {retries} attempts: {last_err}")


def fetch_cached(url, cache_path, timeout=60):
    """Download URL to cache_path if not already cached, and return raw bytes."""
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 0:
        with open(cache_path, "rb") as f:
            return f.read()
    data = http_get_bytes(url, timeout=timeout)
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, "wb") as f:
        f.write(data)
    return data


# ---------------------------------------------------------------------------
# Drug Name Normalization
# ---------------------------------------------------------------------------

def normalize_name(raw: str) -> str:
    """Normalize a drug name to a canonical lowercase lookup key."""
    if not raw:
        return ""
    s = raw.strip().lower()
    # Remove parenthetical notes like "(genetical recombination)", "(USAN)", "(INN)"
    s = re.sub(r"\([^)]*\)", " ", s)
    s = re.sub(r"\[[^\]]*\]", " ", s)
    # Remove dosage numbers like "100 mg", "5%", "20mg/ml"
    s = re.sub(r"\b\d+(\.\d+)?\s*(mg|mcg|ug|g|ml|l|iu|%|units?|mg/ml|mg/m2)\b", " ", s)
    # Replace punctuation with spaces (keep hyphens inside drug codes like "azd-9291")
    s = re.sub(r"[,;:/\\+®™*]", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def strip_salts(norm: str) -> str:
    """Strip salt/hydrate/formulation words from a normalized drug name."""
    if not norm:
        return ""
    tokens = norm.split()
    kept = [t for t in tokens if t not in SALT_WORDS]
    if not kept:
        return norm
    return " ".join(kept)


def generate_name_keys(raw: str):
    """Yield normalized lookup keys for a drug name (exact normalized + salt-stripped + hyphen-free)."""
    norm = normalize_name(raw)
    if not norm or len(norm) < 3:
        return set()
    keys = {norm}
    stripped = strip_salts(norm)
    if stripped and len(stripped) >= 3:
        keys.add(stripped)
    # Also add compact form without hyphens/spaces for investigational codes (e.g. "azd 9291" -> "azd9291")
    compact = re.sub(r"[\s\-]+", "", stripped)
    if compact and len(compact) >= 4:
        keys.add(compact)
    return keys


def split_ingredients(text: str):
    """Split multi-ingredient strings (e.g. 'sofosbuvir; velpatasvir' or 'a / b')."""
    if not text:
        return []
    parts = re.split(r"\s*(?:;|/|(?:\band\b)|\+|,)\s*", text)
    return [p.strip() for p in parts if p and len(p.strip()) >= 3]


def extract_year(date_str: str):
    """Extract a plausible 4-digit approval year (1930-2026) from a date string."""
    if not date_str:
        return None
    s = str(date_str).strip()
    # Handle compact YYYYMMDD format (e.g. openFDA "19870526")
    m_compact = re.fullmatch(r"(19[3-9]\d|20[0-2]\d)\d{4}", s)
    if m_compact:
        yr = int(m_compact.group(1))
        if 1930 <= yr <= 2026:
            return yr
    m = re.search(r"\b(19[3-9]\d|20[0-2]\d)\b", s)
    if m:
        yr = int(m.group(1))
        if 1930 <= yr <= 2026:
            return yr
    return None


def update_earliest_year(mapping: dict, name: str, year: int):
    """Record the earliest approval year for all normalized keys of a drug name."""
    if not year or not name:
        return
    for key in generate_name_keys(name):
        prev = mapping.get(key)
        if prev is None or year < prev:
            mapping[key] = year


# ---------------------------------------------------------------------------
# Step 1: ChEMBL Master Small-Molecule Backbone
# ---------------------------------------------------------------------------

def fetch_chembl_small_molecules(cache_dir: str):
    """
    Fetch all Small molecule entries with max_phase >= 1 from ChEMBL REST API.
    Deduplicates salts into parent molecules using parent_chembl_id and canonical_smiles.
    Returns:
      molecules: dict[mol_id -> {name, smiles, max_phase, first_approval, synonyms}]
      name_to_mol: dict[normalized_key -> set[mol_id]]
    """
    cache_file = os.path.join(cache_dir, "chembl_small_molecules.json")
    if os.path.exists(cache_file) and os.path.getsize(cache_file) > 0:
        print(f"[ChEMBL] Loading cached small molecules from {cache_file}...")
        with open(cache_file, "r", encoding="utf-8") as f:
            raw_mols = json.load(f)
    else:
        print("[ChEMBL] Querying total count of small molecules (max_phase >= 1)...")
        first_url = f"{CHEMBL_BASE}?molecule_type=Small%20molecule&max_phase__gte=1&limit=1000&offset=0"
        first_page = json.loads(http_get_bytes(first_url).decode("utf-8"))
        total = first_page["page_meta"]["total_count"]
        print(f"[ChEMBL] Total records to fetch: {total}")

        raw_mols = list(first_page["molecules"])
        offsets = list(range(1000, total, 1000))

        def fetch_page(offset):
            url = f"{CHEMBL_BASE}?molecule_type=Small%20molecule&max_phase__gte=1&limit=1000&offset={offset}"
            data = json.loads(http_get_bytes(url, timeout=60).decode("utf-8"))
            return data["molecules"]

        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
            futures = {ex.submit(fetch_page, off): off for off in offsets}
            done = 1
            for fut in concurrent.futures.as_completed(futures):
                raw_mols.extend(fut.result())
                done += 1
                print(f"  [ChEMBL] Fetched page {done}/{len(offsets) + 1} ({len(raw_mols)} records)")

        os.makedirs(cache_dir, exist_ok=True)
        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(raw_mols, f)

    # Group by parent_chembl_id (so salt forms merge into one parent row)
    by_parent = {}
    for m in raw_mols:
        # Strict small-molecule guard
        if m.get("molecule_type") != "Small molecule":
            continue
        if m.get("biotherapeutic") is not None:
            continue
        structs = m.get("molecule_structures") or {}
        smiles = (structs.get("canonical_smiles") or "").strip()
        if not smiles:
            continue

        chembl_id = m.get("molecule_chembl_id")
        hier = m.get("molecule_hierarchy") or {}
        parent_id = hier.get("parent_chembl_id") or chembl_id
        is_parent = (chembl_id == parent_id)

        pref_name = (m.get("pref_name") or "").strip()
        max_phase = m.get("max_phase")
        try:
            max_phase_val = float(max_phase) if max_phase is not None else 0.0
        except (ValueError, TypeError):
            max_phase_val = 0.0
        first_app = m.get("first_approval")
        try:
            first_app_yr = int(first_app) if first_app is not None else None
        except (ValueError, TypeError):
            first_app_yr = None

        syns = set()
        official_names = set()
        if pref_name:
            syns.add(pref_name)
            official_names.add(pref_name)
        for s in (m.get("molecule_synonyms") or []):
            syn_txt = (s.get("molecule_synonym") or "").strip()
            syn_type = (s.get("syn_type") or "").upper()
            if not syn_txt or len(syn_txt) < 4:
                continue
            # Skip metabolite annotations like "Mirabegron (m9)" that collide with parent drugs
            if re.search(r"\(m\d+\)|\bmetabolite\b", syn_txt, re.I):
                continue
            syns.add(syn_txt)
            if syn_type in ("INN", "USAN", "BAN", "JAN", "FDA", "USP"):
                official_names.add(syn_txt)

        if parent_id not in by_parent:
            by_parent[parent_id] = {
                "chembl_id": parent_id,
                "name": pref_name,
                "smiles": smiles,
                "is_parent_smiles": is_parent,
                "is_parent_name": is_parent and bool(pref_name),
                "max_phase": max_phase_val,
                "first_approval": first_app_yr,
                "synonyms": syns,
                "official_names": official_names,
            }
        else:
            rec = by_parent[parent_id]
            rec["synonyms"].update(syns)
            rec["official_names"].update(official_names)
            if max_phase_val > rec["max_phase"]:
                rec["max_phase"] = max_phase_val
            if first_app_yr and (rec["first_approval"] is None or first_app_yr < rec["first_approval"]):
                rec["first_approval"] = first_app_yr
            # Prefer parent molecule's canonical SMILES and preferred name
            if is_parent and not rec["is_parent_smiles"]:
                rec["smiles"] = smiles
                rec["is_parent_smiles"] = True
            if is_parent and pref_name:
                rec["name"] = pref_name
                rec["is_parent_name"] = True
            elif not rec["name"] and pref_name:
                rec["name"] = pref_name

    # Second-pass deduplication by canonical SMILES
    by_smiles = {}
    for pid, rec in by_parent.items():
        smi = rec["smiles"]
        if not rec["name"]:
            # Pick shortest readable synonym if pref_name is missing
            readable = sorted(
                [s for s in rec["synonyms"] if len(s) >= 3],
                key=lambda x: (len(x), x),
            )
            rec["name"] = readable[0] if readable else pid

        if smi not in by_smiles:
            by_smiles[smi] = rec
        else:
            existing = by_smiles[smi]
            existing["synonyms"].update(rec["synonyms"])
            existing["official_names"].update(rec["official_names"])
            if rec["max_phase"] > existing["max_phase"]:
                existing["max_phase"] = rec["max_phase"]
            if rec["first_approval"] and (
                existing["first_approval"] is None
                or rec["first_approval"] < existing["first_approval"]
            ):
                existing["first_approval"] = rec["first_approval"]
            # Prefer non-CHEMBL pref_name
            if existing["name"].startswith("CHEMBL") and not rec["name"].startswith("CHEMBL"):
                existing["name"] = rec["name"]

    # Build lookup index: normalized_key -> set of smiles (molecule IDs)
    molecules = {}
    name_to_mols = defaultdict(set)
    for smi, rec in by_smiles.items():
        mol_id = rec["chembl_id"]
        molecules[mol_id] = rec
        all_names = set(rec["synonyms"])
        all_names.add(rec["name"])
        for nm in all_names:
            for key in generate_name_keys(nm):
                name_to_mols[key].add(mol_id)

    print(f"[ChEMBL] Deduplicated to {len(molecules)} unique small molecules with SMILES.")
    return molecules, name_to_mols


# ---------------------------------------------------------------------------
# Step 2: FDA Approvals (openFDA Drugs@FDA + DrugCentral)
# ---------------------------------------------------------------------------

def load_fda_approvals(cache_dir: str):
    """Parse openFDA Drugs@FDA bulk JSON ZIP and DrugCentral FDA_Approved.csv."""
    print("[FDA] Downloading/loading openFDA Drugs@FDA bulk dataset...")
    zip_bytes = fetch_cached(
        OPENFDA_BULK_URL,
        os.path.join(cache_dir, "openfda_drugsfda.json.zip"),
        timeout=90,
    )
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        json_name = zf.namelist()[0]
        data = json.loads(zf.read(json_name).decode("utf-8"))

    fda_years = {}
    for rec in data.get("results", []):
        app_num = (rec.get("application_number") or "").upper()
        # Exclude BLA (Biologics License Applications)
        if app_num.startswith("BLA"):
            continue

        # Find earliest approved submission date (submission_status == "AP")
        earliest_yr = None
        for sub in rec.get("submissions", []):
            if (sub.get("submission_status") or "").upper() == "AP":
                yr = extract_year(sub.get("submission_status_date", ""))
                if yr and (earliest_yr is None or yr < earliest_yr):
                    earliest_yr = yr

        if not earliest_yr:
            continue

        # Collect active ingredient names
        ingredients = set()
        for prod in rec.get("products", []):
            for ai in prod.get("active_ingredients", []):
                nm = (ai.get("name") or "").strip()
                if nm:
                    ingredients.add(nm)
        openfda = rec.get("openfda") or {}
        for nm in openfda.get("generic_name", []):
            for part in split_ingredients(nm):
                ingredients.add(part)
        for nm in openfda.get("substance_name", []):
            for part in split_ingredients(nm):
                ingredients.add(part)

        for ing in ingredients:
            update_earliest_year(fda_years, ing, earliest_yr)

    # Also load DrugCentral FDA_Approved.csv as fallback approval indicator (no year)
    fda_dc_approved = set()
    try:
        dc_bytes = fetch_cached(
            DRUGCENTRAL_FDA_URL,
            os.path.join(cache_dir, "dc_fda_approved.csv"),
        )
        for line in dc_bytes.decode("utf-8", errors="replace").splitlines():
            parts = line.split(",", 1)
            if len(parts) == 2:
                for key in generate_name_keys(parts[1]):
                    fda_dc_approved.add(key)
    except Exception as e:
        print(f"  [FDA] Warning: DrugCentral FDA fallback skipped: {e}")

    print(f"[FDA] Indexed {len(fda_years)} normalized ingredient keys with approval years.")
    return fda_years, fda_dc_approved


# ---------------------------------------------------------------------------
# Step 3: EMA Approvals (Official EMA Medicines JSON + DrugCentral)
# ---------------------------------------------------------------------------

def load_ema_approvals(cache_dir: str):
    """Parse official EMA Medicines JSON report and DrugCentral EMA_Approved.csv."""
    print("[EMA] Downloading/loading official EMA medicines JSON...")
    raw = fetch_cached(
        EMA_JSON_URL,
        os.path.join(cache_dir, "ema_medicines.json"),
        timeout=60,
    )
    ema_data = json.loads(raw.decode("utf-8"))
    records = ema_data.get("data", [])

    ema_years = {}
    for rec in records:
        if (rec.get("category") or "").lower() != "human":
            continue
        # Only count medicines that actually received marketing authorisation
        ma_date = rec.get("marketing_authorisation_date") or ""
        yr = extract_year(ma_date)
        if not yr:
            continue
        status = (rec.get("medicine_status") or "").lower()
        # Exclude refused / withdrawn-before-approval applications
        if "refused" in status or "pending" in status:
            continue

        names = set()
        for field in ("active_substance", "international_non_proprietary_name_common_name"):
            val = rec.get(field) or ""
            for part in split_ingredients(val):
                names.add(part)

        for nm in names:
            update_earliest_year(ema_years, nm, yr)

    ema_dc_approved = set()
    try:
        dc_bytes = fetch_cached(
            DRUGCENTRAL_EMA_URL,
            os.path.join(cache_dir, "dc_ema_approved.csv"),
        )
        for line in dc_bytes.decode("utf-8", errors="replace").splitlines():
            parts = line.split(",", 1)
            if len(parts) == 2:
                for key in generate_name_keys(parts[1]):
                    ema_dc_approved.add(key)
    except Exception as e:
        print(f"  [EMA] Warning: DrugCentral EMA fallback skipped: {e}")

    print(f"[EMA] Indexed {len(ema_years)} normalized ingredient keys with approval years.")
    return ema_years, ema_dc_approved


# ---------------------------------------------------------------------------
# Step 4: PMDA Approvals (Official PMDA PDF + HTML + DrugCentral)
# ---------------------------------------------------------------------------

def extract_pdf_text_macos(pdf_path: str) -> str:
    """Extract full text from a PDF using macOS built-in PDFKit via JXA (osascript)."""
    jxa = f"""
ObjC.import("PDFKit");
var doc = $.PDFDocument.alloc.initWithURL($.NSURL.fileURLWithPath({json.dumps(pdf_path)}));
var out = [];
var n = doc.pageCount;
for (var i = 0; i < n; i++) {{
    var p = doc.pageAtIndex(i);
    if (p && p.string) {{
        out.push(p.string.js);
    }}
}}
out.join("\\n---PAGE---\\n");
"""
    res = subprocess.run(
        ["osascript", "-l", "JavaScript", "-e", jxa],
        capture_output=True,
        text=True,
        timeout=120,
    )
    if res.returncode != 0:
        raise RuntimeError(f"PDFKit extraction failed: {res.stderr}")
    return res.stdout


def load_pmda_approvals(cache_dir: str):
    """
    Parse:
      1. Official PMDA Approved Drugs PDF (April 2004 - Feb 2026)
      2. Official PMDA Review Reports HTML table (0001.html)
      3. DrugCentral PMDA_Approved.csv (fallback for pre-2004 / additional drugs)
    """
    pmda_years = {}

    # 4a. PMDA Review Reports HTML (structured table: Brand | Non-proprietary | Approved In)
    print("[PMDA] Downloading/loading official PMDA Review Reports HTML...")
    html = fetch_cached(
        PMDA_HTML_URL,
        os.path.join(cache_dir, "pmda_0001.html"),
        timeout=45,
    ).decode("utf-8", errors="replace")

    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", html, re.S):
        cells = [
            re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", c)).strip()
            for c in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)
        ]
        if len(cells) >= 3:
            non_prop = cells[1]
            yr = extract_year(cells[2])
            if non_prop and yr:
                for part in split_ingredients(non_prop):
                    update_earliest_year(pmda_years, part, yr)

    # 4b. PMDA Official List of Approved Drugs PDF (April 2004 - Feb 2026)
    print("[PMDA] Downloading/loading official PMDA Approved Drugs PDF (2004-2026)...")
    pdf_path = os.path.join(cache_dir, "pmda_approved.pdf")
    txt_cache = os.path.join(cache_dir, "pmda_approved.txt")
    fetch_cached(PMDA_PDF_URL, pdf_path, timeout=90)

    pdf_text = ""
    if os.path.exists(txt_cache) and os.path.getsize(txt_cache) > 0:
        with open(txt_cache, "r", encoding="utf-8") as f:
            pdf_text = f.read()
    elif sys.platform == "darwin":
        try:
            pdf_text = extract_pdf_text_macos(pdf_path)
            with open(txt_cache, "w", encoding="utf-8") as f:
                f.write(pdf_text)
        except Exception as e:
            print(f"  [PMDA] Warning: PDFKit extraction failed ({e}), using HTML + DrugCentral only.")

    if pdf_text:
        # Parse entries across all PMDA PDF historical date formats:
        #   1. "May 19, 2025" / "Jan. 20, 2022" / "Sept. 25, 2018"
        #   2. "23-Apr-04" (older FY2004-FY2008 tables)
        date_re = re.compile(
            r"(?:"
            r"(?:January|February|March|April|May|June|July|August|September|October|November|December|"
            r"Jan\.|Feb\.|Mar\.|Apr\.|Jun\.|Jul\.|Aug\.|Sep\.|Sept\.|Oct\.|Nov\.|Dec\.)"
            r"\s+\d{1,2},\s+(20\d{2})"
            r"|"
            r"\b\d{1,2}-(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-(\d{2})\b"
            r")"
        )
        matches = list(date_re.finditer(pdf_text))
        for idx, m in enumerate(matches):
            if m.group(1):
                yr = int(m.group(1))
            else:
                yy = int(m.group(2))
                yr = 2000 + yy if yy <= 30 else 1900 + yy
            start = m.end()
            end = matches[idx + 1].start() if idx + 1 < len(matches) else len(pdf_text)
            chunk = pdf_text[start:end]

            # Extract text after "Approval" or "Change" marker up to Notes sentence
            ac_matches = list(re.finditer(r"\b(?:Approval|Change)\b", chunk))
            if not ac_matches:
                continue
            after_ac = chunk[ac_matches[-1].end():].strip()
            # Cut at the start of the Notes column
            note_cut = re.search(
                r"\b(?:A drug|Drugs|A new|An? [a-z]+ drug|Indication|Dosage|Orphan|\[Orphan|"
                r"Re-examination|This drug|New |A combination|A biologic|A vaccine|"
                r"Additional|Change in|Partial|The drug|Products? with|A product)",
                after_ac,
            )
            ing_block = after_ac[:note_cut.start()].strip() if note_cut else after_ac[:140].strip()
            ing_block = re.sub(r"\s+", " ", ing_block).strip()
            if ing_block and len(ing_block) <= 180:
                for part in split_ingredients(ing_block):
                    update_earliest_year(pmda_years, part, yr)

    # 4c. DrugCentral PMDA_Approved.csv (for older pre-2004 PMDA drugs)
    pmda_dc_approved = set()
    try:
        dc_bytes = fetch_cached(
            DRUGCENTRAL_PMDA_URL,
            os.path.join(cache_dir, "dc_pmda_approved.csv"),
        )
        for line in dc_bytes.decode("utf-8", errors="replace").splitlines():
            parts = line.split(",", 1)
            if len(parts) == 2:
                for key in generate_name_keys(parts[1]):
                    pmda_dc_approved.add(key)
    except Exception as e:
        print(f"  [PMDA] Warning: DrugCentral PMDA fallback skipped: {e}")

    print(f"[PMDA] Indexed {len(pmda_years)} normalized ingredient keys with approval years.")
    return pmda_years, pmda_dc_approved


# ---------------------------------------------------------------------------
# Step 5: ClinicalTrials.gov v2 (Clinical Stages by Region: US, EU, China, Japan)
# ---------------------------------------------------------------------------

PHASE_TO_NUM = {
    "EARLY_PHASE1": 1,
    "PHASE1": 1,
    "PHASE2": 2,
    "PHASE3": 3,
    "PHASE4": 4,
}


def fetch_ctgov_phase(phase_code: str, name_to_mols: dict, cache_dir: str):
    """
    Paginate ClinicalTrials.gov API v2 for all interventional DRUG studies of a given phase
    (PHASE1, PHASE2, PHASE3, PHASE4) and map matched small molecules to jurisdictions:
      FDA  <- United States
      EMA  <- EU/EEA member states
      NMPA <- China
      PMDA <- Japan
    Returns dict[agency -> dict[mol_id -> max_phase_num]]
    """
    cache_file = os.path.join(cache_dir, f"ctgov_{phase_code.lower()}_hits.json")
    if os.path.exists(cache_file) and os.path.getsize(cache_file) > 0:
        print(f"[CT.gov] Loading cached {phase_code} matches from {cache_file}...")
        with open(cache_file, "r", encoding="utf-8") as f:
            return json.load(f)

    phase_num = PHASE_TO_NUM[phase_code]
    adv_filter = (
        f"AREA[StudyType]INTERVENTIONAL AND AREA[InterventionType]DRUG "
        f"AND AREA[Phase]{phase_code}"
    )
    fields = (
        "protocolSection.designModule.phases,"
        "protocolSection.armsInterventionsModule.interventions.name,"
        "protocolSection.armsInterventionsModule.interventions.type,"
        "protocolSection.contactsLocationsModule.locations.country"
    )

    hits = {"FDA": {}, "EMA": {}, "NMPA": {}, "PMDA": {}}
    next_token = None
    page = 0
    total_studies = 0

    while True:
        params = {
            "filter.advanced": adv_filter,
            "fields": fields,
            "pageSize": 1000,
        }
        if next_token:
            params["pageToken"] = next_token

        url = f"{CTGOV_API_URL}?{urllib.parse.urlencode(params)}"
        data = json.loads(http_get_bytes(url, timeout=60).decode("utf-8"))
        studies = data.get("studies", [])
        page += 1
        total_studies += len(studies)

        for st in studies:
            ps = st.get("protocolSection") or {}
            # Determine countries
            locs = (ps.get("contactsLocationsModule") or {}).get("locations") or []
            agencies = set()
            for loc in locs:
                c = loc.get("country")
                if not c:
                    continue
                if c == "United States":
                    agencies.add("FDA")
                elif c == "China":
                    agencies.add("NMPA")
                elif c == "Japan":
                    agencies.add("PMDA")
                elif c in EU_EEA_COUNTRIES:
                    agencies.add("EMA")
            if not agencies:
                continue

            # Determine highest phase of this study (e.g. PHASE1 + PHASE2 -> 2)
            st_phases = (ps.get("designModule") or {}).get("phases") or []
            st_phase_num = phase_num
            for p in st_phases:
                pn = PHASE_TO_NUM.get(p, 0)
                if pn > st_phase_num:
                    st_phase_num = pn

            # Match DRUG interventions against ChEMBL small-molecule index
            interventions = (ps.get("armsInterventionsModule") or {}).get("interventions") or []
            matched_mols = set()
            for iv in interventions:
                if iv.get("type") != "DRUG":
                    continue
                iv_name = iv.get("name") or ""
                for part in split_ingredients(iv_name):
                    for key in generate_name_keys(part):
                        mols = name_to_mols.get(key)
                        if mols:
                            matched_mols.update(mols)

            for ag in agencies:
                ag_map = hits[ag]
                for mid in matched_mols:
                    prev = ag_map.get(mid, 0)
                    if st_phase_num > prev:
                        ag_map[mid] = st_phase_num

        if page % 10 == 0:
            print(f"  [CT.gov {phase_code}] Processed {page} pages ({total_studies} studies)...")

        next_token = data.get("nextPageToken")
        if not next_token:
            break

    print(f"[CT.gov {phase_code}] Done: {total_studies} studies across {page} pages.")
    os.makedirs(cache_dir, exist_ok=True)
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump(hits, f)
    return hits


def load_clinical_trials(name_to_mols: dict, cache_dir: str):
    """Fetch Phase 1-4 clinical trial regional mappings in parallel."""
    print("[CT.gov] Fetching Phase 1-4 interventional drug trials across US/EU/China/Japan...")
    combined = {"FDA": defaultdict(int), "EMA": defaultdict(int), "NMPA": defaultdict(int), "PMDA": defaultdict(int)}
    phases = ["PHASE1", "PHASE2", "PHASE3", "PHASE4"]

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        futures = {ex.submit(fetch_ctgov_phase, ph, name_to_mols, cache_dir): ph for ph in phases}
        for fut in concurrent.futures.as_completed(futures):
            ph_hits = fut.result()
            for ag in ("FDA", "EMA", "NMPA", "PMDA"):
                for mid, pnum in ph_hits.get(ag, {}).items():
                    if pnum > combined[ag][mid]:
                        combined[ag][mid] = pnum

    return combined


# ---------------------------------------------------------------------------
# Step 6: Merge All Sources & Write Final CSV
# ---------------------------------------------------------------------------

def lookup_earliest_year(mol_rec: dict, agency_year_map: dict):
    """Find earliest approval year across official non-proprietary names of a molecule."""
    best = None
    all_names = set(mol_rec.get("official_names") or ())
    all_names.add(mol_rec["name"])
    for nm in all_names:
        for key in generate_name_keys(nm):
            yr = agency_year_map.get(key)
            if yr and (best is None or yr < best):
                best = yr
    return best


def lookup_in_set(mol_rec: dict, agency_set: set) -> bool:
    """Check if any official non-proprietary name of a molecule appears in an approval set."""
    all_names = set(mol_rec.get("official_names") or ())
    all_names.add(mol_rec["name"])
    for nm in all_names:
        for key in generate_name_keys(nm):
            if key in agency_set:
                return True
    return False


def format_status(approval_year, is_dc_approved, ct_phase_num, fallback_year=None):
    """
    Determine cell value per user specification:
      - Approval year (e.g. '2019') if approved with known year
      - 'Approved (Year N/A)' if confirmed approved/post-marketing without exact year
      - 'Stage 1' / 'Stage 2' / 'Stage 3' if in clinical trials and not approved
      - 'NA' otherwise
    """
    if approval_year:
        return str(approval_year)
    if is_dc_approved:
        return str(fallback_year) if fallback_year else "Approved (Year N/A)"
    if ct_phase_num == 4:
        return str(fallback_year) if fallback_year else "Approved (Year N/A)"
    if ct_phase_num in (1, 2, 3):
        return f"Stage {ct_phase_num}"
    return "NA"


def build_dataset(output_csv: str, cache_dir: str):
    t0 = time.time()
    molecules, name_to_mols = fetch_chembl_small_molecules(cache_dir)

    # Fetch regulatory approval datasets + clinical trials
    fda_years, fda_dc = load_fda_approvals(cache_dir)
    ema_years, ema_dc = load_ema_approvals(cache_dir)
    pmda_years, pmda_dc = load_pmda_approvals(cache_dir)
    ct_stages = load_clinical_trials(name_to_mols, cache_dir)

    rows = []
    stats = defaultdict(int)

    for mol_id, rec in molecules.items():
        fda_yr = lookup_earliest_year(rec, fda_years)
        ema_yr = lookup_earliest_year(rec, ema_years)
        pmda_yr = lookup_earliest_year(rec, pmda_years)

        fda_is_dc = lookup_in_set(rec, fda_dc)
        ema_is_dc = lookup_in_set(rec, ema_dc)
        pmda_is_dc = lookup_in_set(rec, pmda_dc)

        fda_ct = ct_stages["FDA"].get(mol_id, 0)
        ema_ct = ct_stages["EMA"].get(mol_id, 0)
        nmpa_ct = ct_stages["NMPA"].get(mol_id, 0)
        pmda_ct = ct_stages["PMDA"].get(mol_id, 0)

        # For NMPA (China): open bulk NMPA approval-year APIs do not exist without CAPTCHA.
        # If ChEMBL marks max_phase == 4 with a first_approval year, and the drug is NOT
        # approved by FDA/EMA/PMDA, AND has Phase 3 or Phase 4 trials in China,
        # attribute ChEMBL's first_approval year to NMPA (covers China-first innovative drugs
        # such as icotinib, apatinib, anlotinib, almonertinib, furmonertinib, chidamide, etc.).
        nmpa_yr = None
        if (
            rec["max_phase"] >= 4
            and rec["first_approval"]
            and not fda_yr and not fda_is_dc
            and not ema_yr and not ema_is_dc
            and not pmda_yr and not pmda_is_dc
            and nmpa_ct >= 3
        ):
            nmpa_yr = rec["first_approval"]

        # For FDA/EMA/PMDA confirmed via DrugCentral (pre-2004 PMDA or older NDA)
        # we only use ChEMBL first_approval if it's the sole agency or fallback is unambiguous;
        # otherwise label "Approved (Year N/A)" to avoid misattributing another agency's year.
        fda_val = format_status(fda_yr, fda_is_dc, fda_ct)
        ema_val = format_status(ema_yr, ema_is_dc, ema_ct)
        nmpa_val = format_status(nmpa_yr, False, nmpa_ct)
        pmda_val = format_status(pmda_yr, pmda_is_dc, pmda_ct)

        # Skip molecules with no record in any of the 4 agencies
        if fda_val == "NA" and ema_val == "NA" and nmpa_val == "NA" and pmda_val == "NA":
            continue

        # Clean display name (Title Case for all-caps ChEMBL pref_names, keep codes like AZD9291)
        disp_name = rec["name"]
        if disp_name.isupper() and not re.search(r"\d", disp_name) and len(disp_name) > 4:
            disp_name = disp_name.capitalize()

        rows.append({
            "Name": disp_name,
            "SMILES": rec["smiles"],
            "FDA": fda_val,
            "EMA": ema_val,
            "NMPA": nmpa_val,
            "PMDA": pmda_val,
        })

        for ag, val in (("FDA", fda_val), ("EMA", ema_val), ("NMPA", nmpa_val), ("PMDA", pmda_val)):
            if re.fullmatch(r"\d{4}", val):
                stats[f"{ag}_approved_year"] += 1
            elif val.startswith("Approved"):
                stats[f"{ag}_approved_no_year"] += 1
            elif val.startswith("Stage"):
                stats[f"{ag}_{val.lower().replace(' ', '_')}"] += 1

    # Sort: molecules with more agencies covered first, then alphabetically by Name
    def sort_key(r):
        non_na = sum(1 for k in ("FDA", "EMA", "NMPA", "PMDA") if r[k] != "NA")
        approved_cnt = sum(
            1 for k in ("FDA", "EMA", "NMPA", "PMDA")
            if re.fullmatch(r"\d{4}", r[k]) or r[k].startswith("Approved")
        )
        return (-approved_cnt, -non_na, r["Name"].lower())

    rows.sort(key=sort_key)

    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["Name", "SMILES", "FDA", "EMA", "NMPA", "PMDA"])
        writer.writeheader()
        writer.writerows(rows)

    elapsed = time.time() - t0
    print("\n================ SUMMARY ================")
    print(f"Output CSV: {output_csv}")
    print(f"Total small-molecule drugs written: {len(rows)}")
    print(f"Elapsed time: {elapsed:.1f}s")
    for ag in ("FDA", "EMA", "NMPA", "PMDA"):
        print(
            f"  {ag:4s}: "
            f"Approved(year)={stats[f'{ag}_approved_year']}, "
            f"Approved(year N/A)={stats[f'{ag}_approved_no_year']}, "
            f"Stage 3={stats[f'{ag}_stage_3']}, "
            f"Stage 2={stats[f'{ag}_stage_2']}, "
            f"Stage 1={stats[f'{ag}_stage_1']}"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Collect small-molecule drugs with SMILES and FDA/EMA/NMPA/PMDA status."
    )
    parser.add_argument(
        "--output",
        default="small_molecule_drugs_Gemini.csv",
        help="Path to output CSV file (default: small_molecule_drugs_Gemini.csv)",
    )
    parser.add_argument(
        "--cache-dir",
        default="/tmp/drug_collector_cache",
        help="Directory for caching downloaded source datasets",
    )
    args = parser.parse_args()
    build_dataset(args.output, args.cache_dir)


if __name__ == "__main__":
    main()
