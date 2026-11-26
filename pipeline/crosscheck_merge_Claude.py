#!/usr/bin/env python3
"""
Final merge: cross-check the three AI agents' round-1 tables and add the extra sources.

Inputs
  round1_agents/claude/small_molecule_drugs_by_agency_detailed_Claude.csv  -> "CL"
  round1_agents/gemini/small_molecule_drugs_Gemini.csv                     -> "GE"
  round1_agents/deepseek/drug_approvals_smiles_Deepseek.csv                -> "DS"
  pipeline/nmpa_fix/nmpa_rescued_Claude.csv        -> "NM"  China records rescued via Chinese names
  pipeline/anvisa/anvisa_rows_Claude.csv           -> "BR"  ANVISA (Brazil)
  pipeline/eu_national/eu_national_rows_Claude.csv -> "EU"  EU national authorisations
  pipeline/curation/verified_additions_Claude.csv  -> "VA"  web-verified additions (win over agents)
  pipeline/curation/verified_overrides_Claude.csv  manual, web-verified cell corrections
  pipeline/curation/name_aliases_Claude.csv        same drug under different names
  data/dc_approvals.csv                            DrugCentral, only to fill molecules no agent had

Outputs
  [BEST]final_small_molecule_drugs_Claude.csv       Name,SMILES,FDA,EMA,EU_National,NMPA,PMDA,ANVISA
  results/final_small_molecule_drugs_detailed_Claude.csv  + InChIKey, each source's raw value, evidence
  results/crosscheck_conflicts_Claude.csv           every cell where the three agents disagree

Cell vocabulary
  2019            approved by that agency, year of first approval
  Approved        approved / marketed there, no reliable first-approval year
  Stage 1..4      not (known to be) approved; highest clinical-trial phase with a site in that region
                  (Stage 1/2, Stage 2/3 = combined phases; Stage 0/1 = early phase 1; Stage 4 = post-marketing)
  NA              no record found
  EMA = EU centralised procedure only; EU_National = national / MRP / DCP authorisation in any EU/EEA state.

Usage:  python3 pipeline/crosscheck_merge_Claude.py      (pandas, rdkit)
"""
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]          # repo root
PIPE = ROOT / "pipeline"
R1 = ROOT / "round1_agents"                        # the three agents' round-1 tables
OUT = ROOT / "results"
OVERRIDES = PIPE / "curation" / "verified_overrides_Claude.csv"
ALIASES = PIPE / "curation" / "name_aliases_Claude.csv"          # same drug, different names
# Extra sources added after the cross-check (same columns as the agent files; blank = no info):
EXTRA = {
    # NMPA database records (DeepSeek's snapshot) whose names only resolved through their
    # Chinese name (Wikidata zh labels): only the NMPA column is filled.
    "NM": PIPE / "nmpa_fix" / "nmpa_rescued_Claude.csv",
    # Web-verified additions (China-first innovative drugs, Carrimycin components ...).
    # A non-blank cell here overrides the agents.
    "VA": PIPE / "curation" / "verified_additions_Claude.csv",
    # ANVISA (Brazil): first registration year from ANVISA open data + Brazil trial phase
    # (built by crosscheck/anvisa/build_anvisa.py); only the ANVISA column is filled.
    "BR": PIPE / "anvisa" / "anvisa_rows_Claude.csv",
    # EU national / MRP / DCP authorisations: France ANSM + Spain AEMPS dates, EMA Article 57 presence
    # (built by crosscheck/eu_national/build_eu_national.py); only the EU_National column is filled.
    "EU": PIPE / "eu_national" / "eu_national_rows_Claude.csv",
}
BEST = ROOT / "[BEST]final_small_molecule_drugs_Claude.csv"
SRC_RANK = {"CL": 0, "GE": 1, "VA": 2, "BR": 3, "EU": 4, "DS": 5, "NM": 6}
FILES = {
    "CL": "claude/small_molecule_drugs_by_agency_detailed_Claude.csv",
    "GE": "gemini/small_molecule_drugs_Gemini.csv",
    "DS": "deepseek/drug_approvals_smiles_Deepseek.csv",
}
AG = ["FDA", "EMA", "EU_National", "NMPA", "PMDA", "ANVISA"]
SALT_WORDS = {
    "sodium", "disodium", "potassium", "calcium", "magnesium", "hydrochloride", "hcl",
    "sulfate", "mesylate", "tosylate", "citrate", "maleate", "fumarate", "tartrate",
    "acetate", "phosphate", "salt", "anhydrous", "hydrate", "monohydrate", "dihydrate",
    "pentahydrate", "hemiheptahydrate", "besylate", "bromide", "chloride", "succinate",
}

# --------------------------------------------------------------------------- structure keys
_lfc = rdMolStandardize.LargestFragmentChooser()
_unch = rdMolStandardize.Uncharger()


def parent_key(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return "", 0
    try:
        p = _unch.uncharge(_lfc.choose(m))
    except Exception:
        p = m
    return Chem.MolToInchiKey(p) or "", p.GetNumHeavyAtoms()


def norm_name(n):
    return re.sub(r"[\s\-_,'’()\[\]\.]", "", n.lower())


# --------------------------------------------------------------------------- cell parsing
def parse(v):
    """-> (kind, num): kind in Y (year) / A (approved, no year) / P (phase) / N (NA)."""
    if not isinstance(v, str) or v.strip() == "":
        return ("M", None)  # molecule absent from that agent's file
    v = v.strip()
    if re.fullmatch(r"\d{4}", v):
        return ("Y", int(v))
    if v.lower().startswith("approved"):
        return ("A", None)
    m = re.search(r"(Phase|Stage)\s*(\d)(?:/(\d))?", v)
    if m:
        hi = int(m.group(3) or m.group(2))
        lo = int(m.group(2))
        return ("P", (hi + lo) / 2 if m.group(3) else float(hi))
    if "early phase 1" in v.lower():
        return ("P", 0.5)
    return ("N", None)


def phase_label(x):
    if x is None:
        return "NA"
    if x == 0.5:
        return "Stage 0/1"
    if x != int(x):
        return f"Stage {int(x)}/{int(x) + 1}"
    return f"Stage {int(x)}"


# Per-agency preference when the agents give different years and there is no majority.
# Derived from the web-verified sample (see crosscheck/verification_log_Claude.csv).
YEAR_PRIORITY = {
    "FDA": ["CL", "GE", "DS"],
    "EMA": ["CL", "DS", "GE"],
    "PMDA": ["CL", "DS", "GE"],
    "NMPA": ["CL", "GE", "DS"],
    "ANVISA": ["BR"],
    "EU_National": ["EU"],
}


def decide(agency, vals):
    """vals: {src: raw string or None}. Returns (final_value, evidence_tag)."""
    va = vals.get("VA")
    if va:
        return va, "verified-addition"
    p = {s: parse(v) for s, v in vals.items() if s != "VA"}
    years = {s: n for s, (k, n) in p.items() if k == "Y"}
    appr = {s for s, (k, _) in p.items() if k in ("Y", "A")}
    phases = [n for s, (k, n) in p.items() if k == "P"]

    # Gemini turns a Phase-4 trial into "Approved (Year N/A)". When Claude (same trial source)
    # reports only "Phase 4" for that region, treat Gemini's claim as a Phase-4 trial, not an approval.
    if "GE" in appr and "GE" not in years and p.get("CL", ("M",))[0] == "P" and p["CL"][1] == 4.0:
        appr.discard("GE")
        phases.append(4.0)

    # A lone Gemini "Approved (Year N/A)" that no other agent supports was wrong in every
    # verified case (huperzine A, idebenone, levofenfluramine, dexlansoprazole-PMDA ...):
    # fall back to the other agents' trial data when they have the molecule.
    others_have = any(k in ("P", "N") for s, (k, _) in p.items() if s != "GE")
    if appr == {"GE"} and "GE" not in years and others_have:
        appr = set()

    # PMDA: Claude's KEGG-based years were right in 15/17 verified conflicts; Gemini/DeepSeek
    # usually give the year of a later product or formulation.
    if agency == "PMDA" and "CL" in appr:
        if "CL" in years:
            return str(years["CL"]), f"year-CL-KEGG|appr:{''.join(sorted(appr))}"
        other = [y for s, y in years.items() if s != "CL"]
        if len(other) >= 2 and len(set(other)) == 1:
            return str(other[0]), f"year-GE+DS-agree|appr:{''.join(sorted(appr))}"
        return "Approved", f"approved-no-year|appr:{''.join(sorted(appr))}"

    if appr:
        if years:
            c = Counter(years.values()).most_common()
            if len(c) == 1 or c[0][1] > c[1][1]:
                y = c[0][0]
                tag = "year-agree" if len(years) > 1 and len(c) == 1 else (
                    "year-majority" if len(years) > 1 else "year-single")
            else:
                src = next((s for s in YEAR_PRIORITY[agency] if s in years), sorted(years)[0])
                y, tag = years[src], f"year-conflict->{src}"
            if agency == "NMPA" and set(years) <= {"DS", "NM"}:
                # DeepSeek's NMPA years come from a database snapshot whose dates cluster in the
                # 1998-2003 nationwide re-registration (979 of 1,563 years; 527 in 2002 alone) and
                # were wrong in 5 of 7 verified cases -> keep the approval, drop the year.
                if y <= 2003:
                    return "Approved", "nmpa-db-reregistration-year|appr:" + "".join(sorted(appr))
                tag += "|nmpa-db-year"
            return str(y), tag + f"|appr:{''.join(sorted(appr))}"
        return "Approved", f"approved-no-year|appr:{''.join(sorted(appr))}"
    if phases:
        return phase_label(max(phases)), "trial-only"
    return "NA", "none"


# --------------------------------------------------------------------------- main
def main():
    frames = []
    for src, fn in FILES.items():
        d = pd.read_csv(R1 / fn, dtype=str, keep_default_na=False)
        for a in AG:
            if a not in d:
                d[a] = ""
        # The agents' EMA column mixes EMA (centralised, dated) and national authorisations (Claude:
        # EMA Article 57, no year). Split it: an undated "Approved" goes to EU_National; trial phases
        # (EU sites) apply to both columns; years stay EMA-only.
        kind = d.EMA.map(lambda v: parse(v)[0])
        d["EU_National"] = d.EMA.where(kind.isin(["A", "P"]), "")
        d.loc[kind.eq("A"), "EMA"] = ""
        d = d[["Name", "SMILES"] + AG].copy()
        d["src"] = src
        frames.append(d)
    for src, fn in EXTRA.items():
        if Path(fn).exists():
            d = pd.read_csv(fn, dtype=str, keep_default_na=False)
            for a in AG:
                if a not in d:
                    d[a] = ""
            d = d[["Name", "SMILES"] + AG].copy()
            d["src"] = src
            frames.append(d)
    rows = pd.concat(frames, ignore_index=True)
    keys = rows.SMILES.map(parent_key)
    rows["ik"] = keys.str[0]
    rows["ha"] = keys.str[1]
    rows["nname"] = rows.Name.map(norm_name)
    rows["rid"] = range(len(rows))

    # ---- link rows into molecules (union-find): same name, or same parent InChIKey when it is
    #      unambiguous in each file, or same parent InChIKey of a real organic parent (>=8 heavy atoms;
    #      this folds salt/ester-name duplicates such as "Sodium cromoglicate" into "Cromolyn" but keeps
    #      e.g. "Lithium carbonate" and "Calcium carbonate" apart).
    par = list(range(len(rows)))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    def u(a, b):
        par[f(a)] = f(b)

    for _, g in rows.groupby("nname"):
        ids = g.rid.tolist()
        for x in ids[1:]:
            u(ids[0], x)
    if ALIASES.exists():
        first = rows.groupby("nname").rid.first()
        for a, b in pd.read_csv(ALIASES, dtype=str).itertuples(index=False):
            ka, kb = norm_name(a), norm_name(b)
            if ka in first.index and kb in first.index:
                u(first[ka], first[kb])
    dup_in_file = rows[rows.ik != ""].groupby(["src", "ik"]).size()
    for ik, g in rows[rows.ik != ""].groupby("ik"):
        uniq = all(dup_in_file[(s, ik)] == 1 for s in g.src.unique())
        if (uniq and g.src.nunique() > 1) or g.ha.iloc[0] >= 8:
            ids = g.rid.tolist()
            for x in ids[1:]:
                u(ids[0], x)
    rows["ent"] = rows.rid.map(f)

    # DrugCentral approvals (the source DeepSeek used), only to fill molecules that came in from the
    # China/Brazil/verified additions and are absent from all three agent tables (e.g. cisplatin).
    dc = pd.read_csv(ROOT / "data" / "dc_approvals.csv", dtype=str, keep_default_na=False)
    dc_by = {}
    for r in dc.itertuples():
        for k in (norm_name(r.name), r.inchikey):
            if k:
                dc_by[k] = {"FDA": r.fda_date[:4], "EMA": r.ema_date[:4], "PMDA": r.pmda_date[:4]}

    ov = pd.read_csv(OVERRIDES, dtype=str, keep_default_na=False) if OVERRIDES.exists() else pd.DataFrame(
        columns=["name", "field", "value", "source"])
    ov["key"] = ov.name.map(norm_name)

    out, conflicts = [], []
    for ent, g in rows.groupby("ent"):
        # display name: prefer a name without salt words, from CL > GE > DS, shortest
        def name_rank(r):
            words = set(re.findall(r"[a-z]+", r.Name.lower()))
            return (r.src != "VA", bool(words & SALT_WORDS), SRC_RANK[r.src], len(r.Name))
        best = min(g.itertuples(), key=name_rank)
        # spell out a leading Greek letter (Wikidata labels such as "β-lactose" are title-cased to "Β-Lactose")
        name = re.sub(r"^[Αα](?=-)", "alpha", re.sub(r"^[Ββ](?=-)", "beta", re.sub(r"^[Γγ](?=-)", "gamma", best.Name)))
        # SMILES: ChEMBL-based (CL/GE agree 100%) first, else DeepSeek
        smi_row = min(g.itertuples(), key=lambda r: (SRC_RANK[r.src], name_rank(r)))
        rec = {"Name": name, "SMILES": smi_row.SMILES, "InChIKey": smi_row.ik,
               "Names_in_sources": "; ".join(sorted(set(g.Name))),
               **{f"In_{s}": int(s in set(g.src)) for s in SRC_RANK}}
        for a in AG:
            vals = {}
            for s in SRC_RANK:
                gs = g[g.src == s]
                cells = [v for v in gs[a] if parse(v)[0] != "M"]
                if not cells:
                    vals[s] = None
                    continue
                # several rows from one source folded together: keep its strongest claim
                cand = sorted(cells, key=lambda v: {"Y": 0, "A": 1, "P": 2, "N": 3}[parse(v)[0]]
                              * 10000 + (parse(v)[1] or 0) * (1 if parse(v)[0] == "Y" else -1))
                vals[s] = cand[0]
            val, tag = decide(a, vals)
            rec[a] = val
            rec[a + "_evidence"] = tag
            for s in SRC_RANK:
                rec[f"{a}_{s}"] = vals[s] if vals[s] is not None else ""
            kinds = {s: parse(v) for s, v in vals.items() if v is not None and s in ("CL", "GE", "DS")}
            appr_set = {s for s, (k, _) in kinds.items() if k in ("Y", "A")}
            yrs = {n for k, n in kinds.values() if k == "Y"}
            if len(kinds) > 1 and ((0 < len(appr_set) < len(kinds)) or len(yrs) > 1):
                conflicts.append({"Name": name, "Agency": a, **{s: vals[s] or "(absent)" for s in ("CL", "GE", "DS")},
                                  "Final": val, "Rule": tag})
        if not ({"CL", "GE", "DS"} & set(g.src)):
            hit = next((dc_by[k] for k in [norm_name(n) for n in g.Name] + list(g.ik) if k in dc_by), None)
            if hit:
                for a, y in hit.items():
                    if y and not re.fullmatch(r"\d{4}", rec[a]):
                        rec[a], rec[a + "_evidence"] = y, "drugcentral-fill"
        # ANVISA years before 2000 are often only the earliest record in ANVISA's open data (records from
        # before ANVISA's creation in 1999 are incomplete): when the drug was approved elsewhere 10+ years
        # earlier, keep the approval but drop the year (verified: metamizole, captopril, oxamniquine).
        ya = rec["ANVISA"]
        if re.fullmatch(r"\d{4}", ya) and int(ya) < 2000:
            other = [int(rec[a]) for a in ("FDA", "EMA", "EU_National", "PMDA") if re.fullmatch(r"\d{4}", rec[a])]
            if other and min(other) <= int(ya) - 10:
                rec["ANVISA"], rec["ANVISA_evidence"] = "Approved", f"anvisa-early-record-{ya}"
        # manual, web-verified overrides
        for o in ov[ov.key.isin({norm_name(n) for n in g.Name})].itertuples():
            if o.field in AG:
                rec[o.field] = o.value
                rec[o.field + "_evidence"] = "verified:" + o.source
            elif o.field == "SMILES":
                rec["SMILES"] = o.value
                rec["InChIKey"] = parent_key(o.value)[0]
                rec["SMILES_note"] = "verified:" + o.source
            elif o.field == "Name":
                rec["Name"] = o.value
            elif o.field == "DROP":
                rec["_drop"] = o.source
        out.append(rec)

    df = pd.DataFrame(out)
    dropped = df[df.get("_drop", pd.Series(dtype=str)).notna()] if "_drop" in df else df.iloc[0:0]
    if "_drop" in df:
        df = df[df["_drop"].isna()].drop(columns="_drop")
    df = df[~(df[AG] == "NA").all(axis=1)]
    df = df.sort_values("Name", key=lambda s: s.str.lower()).reset_index(drop=True)

    df[["Name", "SMILES"] + AG].to_csv(BEST, index=False)
    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "final_small_molecule_drugs_detailed_Claude.csv", index=False)
    pd.DataFrame(conflicts).to_csv(OUT / "crosscheck_conflicts_Claude.csv", index=False)

    print(f"molecules: {len(df)}  (dropped by override: {len(dropped)})  conflict cells: {len(conflicts)}")
    for a in AG:
        v = df[a]
        print(f"  {a:5s} year {v.str.fullmatch(r'\d{4}').sum():5d}  Approved {(v == 'Approved').sum():5d}  "
              f"trial-only {v.str.startswith('Stage').sum():5d}  NA {(v == 'NA').sum():5d}")


if __name__ == "__main__":
    sys.exit(main())
