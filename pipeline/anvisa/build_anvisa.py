"""ANVISA (Brazil) column: first registration year per active ingredient + Brazil trial phase.

Sources
  raw/medicamentos.csv  ANVISA open data, DADOS_ABERTOS_MEDICAMENTOS.csv (all registered medicines)
  data/raw/              Claude round-1 caches: ChEMBL molecules + synonyms, ClinicalTrials.gov studies
Output
  anvisa_ingredients_Claude.csv   every ANVISA ingredient name -> matched ChEMBL parent (audit)
  anvisa_rows_Claude.csv          Name,SMILES,ANVISA  (source "BR" for crosscheck_merge_Claude.py)
"""
import re, sys, unicodedata
from collections import defaultdict
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "round1_agents" / "claude"))
import collect_drugs as cd  # noqa: E402  (Claude round-1 script: cached ChEMBL index + trial download)

# "Biológico" is kept because Brazil registers synthetic peptides (e.g. semaglutide) there; proteins drop out
# later because only small-molecule ChEMBL parents can be matched.
KEEP_CAT = {"Novo", "Genérico", "Similar", "Específico", "Radiofármaco", "Gases Medici", "Biológico"}


def process_year(p):
    """Filing year encoded in 15/17-digit ANVISA process numbers (25000.xxxxxx/YY-DD, 25351.xxxxxx/YYYY-DD)."""
    p = re.sub(r"\D", "", str(p))
    if len(p) == 17:
        y = int(p[11:15])
    elif len(p) == 15:
        y = int(p[11:13])
    else:  # old 13-digit numbers do not reliably end in the year (e.g. CETOMED -> "40" for a 2000 grant)
        return None
    if y < 100:
        y += 2000 if y <= 26 else 1900
    return y if 1940 <= y <= 2026 else None
PT_SALT_PRE = """cloridrato dicloridrato tricloridrato bromidrato iodidrato maleato dimaleato fumarato
hemifumarato besilato mesilato dimesilato tosilato esilato edisilato succinato hemissuccinato
tartarato hemitartarato bitartarato citrato dicitrato sulfato hemissulfato bissulfato fosfato difosfato
lactato gluconato malato oxalato pamoato embonato nitrato xinafoato hiclato napsilato
iodeto brometo cloreto trometamol meglumina lisinato""".split()
PT_ESTER_PRE = """acetato diacetato valerato propionato dipropionato butirato cipionato enantato decanoato
undecilato palmitato estearato benzoato furoato hexacetonido aceponato""".split()
PT_TRAIL = r"(\s+(de\s+)?(sodio|sodica|sodico|dissodico|dissodica|trissodico|potassio|potassica|potassico|calcio|calcica|calcico|magnesio|magnesica|magnesico|zinco|amonio|litio|monoidratad[oa]|mono-hidratad[oa]|di-?hidratad[oa]|tri-?hidratad[oa]|hemi-?hidratad[oa]|sesqui-?hidratad[oa]|pentaidratad[oa]|hidratad[oa]|anidr[oa]|base|micronizad[oa]|racemic[oa]))+$"
EN_SALT = {"cloridrato": "hydrochloride", "dicloridrato": "dihydrochloride", "bromidrato": "hydrobromide",
           "maleato": "maleate", "fumarato": "fumarate", "besilato": "besylate", "mesilato": "mesylate",
           "succinato": "succinate", "tartarato": "tartrate", "citrato": "citrate", "sulfato": "sulfate",
           "fosfato": "phosphate", "acetato": "acetate", "valerato": "valerate", "propionato": "propionate",
           "dipropionato": "dipropionate", "cipionato": "cypionate", "enantato": "enanthate",
           "decanoato": "decanoate", "furoato": "furoate", "benzoato": "benzoate", "lactato": "lactate",
           "gluconato": "gluconate", "pamoato": "pamoate", "nitrato": "nitrate", "tosilato": "tosylate"}


def deacc(s):
    return "".join(c for c in unicodedata.normalize("NFKD", s) if not unicodedata.combining(c))


def translit(w):
    """Portuguese DCB spelling -> English INN spelling (heuristic)."""
    rules = [(r"ph", "ph"), (r"ina$", "ine"), (r"ino$", "ine"), (r"ona$", "one"), (r"ida$", "ide"),
             (r"eno$", "ene"), (r"ano$", "ane"), (r"ana$", "an"), (r"ila$", "il"), (r"ato$", "ate"),
             (r"ol[a]$", "ol"), (r"ico$", "ic"), (r"(?<=[bcdfgkmptvx])e$", ""), (r"^anlod", "amlod")]
    out = w
    for a, b in rules:
        out = re.sub(a, b, out)
    out = out.replace("qu", "qu").replace("ç", "c")
    return out


def candidates(name):
    n = deacc(str(name).lower()).strip()
    n = re.sub(r"\s*\(.*?\)", "", n)
    n = re.sub(r"\s+", " ", n).strip(" ,.;")
    cands = [n]
    m = re.match(r"^(\w+)\s+de\s+(.+)$", n)
    core, salt = n, None
    if m and m.group(1) in PT_SALT_PRE + PT_ESTER_PRE:
        salt, core = m.group(1), m.group(2)
    core = re.sub(PT_TRAIL, "", core).strip()
    cands.append(core)
    if core.startswith("acido "):
        acid = core[6:]
        cands += [acid + " acid", translit(acid) + " acid", acid]
    words = core.split()
    tw = " ".join(translit(w) for w in words)
    cands.append(tw)
    if salt and salt in EN_SALT:
        cands.insert(0, f"{tw} {EN_SALT[salt]}")
    return [c for i, c in enumerate(cands) if c and c not in cands[:i]]


PT_TRAIL_WORDS = set("""sodio sodica sodico dissodico dissodica trissodico potassio potassica potassico calcio
calcica calcico magnesio magnesica magnesico zinco amonio litio monoidratado monoidratada mono-hidratado
mono-hidratada di-hidratado di-hidratada diidratado diidratada tri-hidratado tri-hidratada triidratado
hemi-hidratado hemiidratado sesqui-hidratado pentaidratado pentaidratada hexaidratado hexaidratada
hemieptaidratada hemieptaidratado hidratado hidratada anidro anidra base micronizado micronizada
racemico racemica dietilamonio resinato de""".split())


def pt_variants(name):
    """Accented Portuguese variants (for Wikidata pt labels): full name and the bare active moiety."""
    toks = re.sub(r"\s*\(.*?\)", "", str(name).lower()).split()
    out = [" ".join(toks)]
    t = list(toks)
    if len(t) > 2 and deacc(t[0]) in PT_SALT_PRE + PT_ESTER_PRE and t[1] == "de":
        t = t[2:]
    while len(t) > 1 and deacc(t[-1]) in PT_TRAIL_WORDS:
        t.pop()
    out.append(" ".join(t))
    out += [o[:1].upper() + o[1:] for o in out]
    return [o for i, o in enumerate(out) if o and o not in out[:i]]


def wikidata_pt(names, parents):
    """name -> ChEMBL parent id (or ('WD', smiles, label) when Wikidata knows a structure ChEMBL lacks)."""
    import subprocess
    here = Path(__file__).parent
    hits_file = here / "wikidata_pt_hits.csv"
    var = {n: pt_variants(n) for n in names}
    if not hits_file.exists():
        allv = sorted({v for vs in var.values() for v in vs})
        (here / "pt_names.txt").write_text("\n".join(allv), encoding="utf-8")
        subprocess.run([sys.executable, str(ROOT / "pipeline" / "common" / "wikidata_labels.py"),
                        str(here / "pt_names.txt"), "pt,pt-br,es", str(hits_file)], check=True)
    h = pd.read_csv(hits_file, dtype=str)
    ik2cid = {}
    for cid, p in parents.items():
        if p.get("inchikey"):
            ik2cid[p["inchikey"]] = cid
    by_lab = {}
    for lab, g in h.groupby("lab"):
        cids = {ik2cid.get(ik) for ik in g.ik.dropna()} - {None}
        if len(cids) == 1:
            by_lab[lab] = next(iter(cids))
            continue
        smis = g.ismi.fillna(g.smi).dropna()
        if len(cids) == 0 and len(smis) and g["item"].nunique() == 1:
            by_lab[lab] = ("WD", smis.iloc[0], g.itemLabel.iloc[0])
    res = {}
    for n, vs in var.items():
        for v in vs:
            if v in by_lab:
                res[n] = by_lab[v]
                break
    return res


def main():
    raw = Path(__file__).parent / "raw" / "medicamentos.csv"
    if not raw.exists():
        import urllib.request
        raw.parent.mkdir(exist_ok=True)
        urllib.request.urlretrieve("https://dados.anvisa.gov.br/dados/DADOS_ABERTOS_MEDICAMENTOS.csv", raw)
    parents, child_of, index = cd.load_chembl()
    d = pd.read_csv(Path(__file__).parent / "raw" / "medicamentos.csv", sep=";", encoding="latin1", dtype=str)
    d = d[d.CATEGORIA_REGULATORIA.isin(KEEP_CAT)].copy()
    d["date"] = pd.to_datetime(d.DATA_FINALIZACAO_PROCESSO, format="%d/%m/%Y", errors="coerce")
    d["year"] = d.date.dt.year
    d.loc[~d.year.between(1950, 2026), "year"] = pd.NA   # 1901/1917/0001 are placeholders
    # The grant date is sometimes that of a later renewal/re-registration. When the process number
    # was filed more than 3 years earlier, use the filing year instead.
    py = d.NUMERO_PROCESSO.map(process_year)
    use_py = py.notna() & (d.year.isna() | (py < d.year - 3))
    d.loc[use_py, "year"] = py[use_py]
    ing = d.PRINCIPIO_ATIVO.fillna("")
    gen = d.CATEGORIA_REGULATORIA.eq("Genérico") & ing.eq("")
    ing[gen] = d.NOME_PRODUTO[gen]                       # generics are named by their DCB
    d["ing"] = ing
    d = d[d.ing.ne("")]

    rows = []
    for r in d.itertuples():
        parts = [p.strip() for p in re.split(r",|\s\+\s|;", r.ing) if p.strip()]
        for p in parts:
            rows.append((p, r.year, len(parts) > 1, r.CATEGORIA_REGULATORIA, r.SITUACAO_REGISTRO))
    x = pd.DataFrame(rows, columns=["ing", "year", "combo", "cat", "status"])
    names = x.ing.str.lower().unique()
    manual = pd.read_csv(Path(__file__).parent / "dcb_to_inn_Claude.csv", dtype=str)
    manual = dict(zip(manual.dcb.str.lower(), manual.inn))
    hit = {}
    for n in names:
        cid = index.lookup(manual[n]) if n in manual else None
        for c in ([] if cid else candidates(n)):
            cid = index.lookup(c)
            if cid:
                break
        hit[n] = cid
    wd = wikidata_pt([n for n, c in hit.items() if not c], parents)
    extra = {}                                           # Wikidata-only structures
    for n, v in wd.items():
        if isinstance(v, tuple):
            extra["WD:" + v[2]] = {"id": "WD:" + v[2], "name": v[2], "smiles": v[1], "inchikey": ""}
            hit[n] = "WD:" + v[2]
        else:
            hit[n] = v
    parents = {**parents, **extra}
    print(f"matched via ChEMBL names: {sum(1 for n in hit if hit[n] and n not in wd)}, via Wikidata pt: {len(wd)}")
    x["cid"] = x.ing.str.lower().map(hit)
    audit = x.groupby(x.ing.str.lower()).agg(n=("year", "size"), first_year=("year", "min"),
                                              cid=("cid", "first")).reset_index().rename(columns={"ing": "ingredient"})
    audit["name"] = audit.cid.map(lambda c: parents[c]["name"] if isinstance(c, str) else "")
    audit["via"] = audit.ingredient.map(lambda n: "wikidata-pt" if n in wd else ("chembl-name" if hit.get(n) else ""))
    audit.sort_values("n", ascending=False).to_csv(Path(__file__).parent / "anvisa_ingredients_Claude.csv", index=False)

    # first registration: single-ingredient products first; an ingredient seen only in combos uses those
    m = x[x.cid.notna()]
    single = m[~m.combo].groupby("cid").year.min()
    anyp = m.groupby("cid").year.min()
    first = single.combine_first(anyp)

    # Brazil trial phases
    cd.REGION_OF_COUNTRY.clear()
    cd.REGION_OF_COUNTRY["Brazil"] = "ANVISA"
    ph = cd.trial_phases_by_region(cd.load_trials(), cd.load_chembl_trial_links(), child_of, index)

    out = []
    for cid in set(first.index) | {c for c, v in ph.items() if "ANVISA" in v}:
        p = parents[cid]
        if cid in first.index:
            y = first[cid]
            val = str(int(y)) if pd.notna(y) else "Approved"
        else:
            val = cd.PHASE_LABEL[ph[cid]["ANVISA"]]
        out.append({"Name": (p["name"] or cid).title(), "SMILES": p["smiles"], "ANVISA": val, "ChEMBL_ID": cid})
    o = pd.DataFrame(out)
    # same scope rule as the final table: small molecules up to 2,000 Da (drops e.g. semaglutide)
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    mw = o.SMILES.map(lambda s: Descriptors.MolWt(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else 0)
    o = o[mw <= 2000]
    o.to_csv(Path(__file__).parent / "anvisa_rows_Claude.csv", index=False)
    tot = x.ing.str.lower().nunique()
    print(f"ANVISA product-ingredient records: {len(x)} | unique ingredient names {tot} | matched {sum(1 for v in hit.values() if v)}")
    print(f"matched records: {x.cid.notna().mean():.1%}")
    print(o.ANVISA.str.fullmatch(r"\d{4}").sum(), "with year;", (o.ANVISA == "Approved").sum(), "approved no year;",
          o.ANVISA.str.startswith(("Phase", "Early")).sum(), "trial-only")


if __name__ == "__main__":
    main()
