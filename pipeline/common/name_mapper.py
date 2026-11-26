"""Map non-English drug ingredient names (French, Spanish, Portuguese) to ChEMBL parents.

Order: (1) ChEMBL name/synonym index built by collect_drugs.py, tried on several spelling variants
(salt words stripped, 'acide/ácido X' -> 'X acid', local INN endings -> English); (2) exact Wikidata
labels in that language (wikidata_labels.py), linked to ChEMBL by InChIKey.
"""
import re
import subprocess
import sys
import unicodedata
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent

SALTS = {
    "fr": """chlorhydrate dichlorhydrate bromhydrate iodhydrate maleate fumarate hemifumarate besilate mesilate
    dimesilate tosilate succinate tartrate hemitartrate bitartrate citrate sulfate hemisulfate phosphate acetate
    valerate propionate dipropionate butyrate lactate gluconate nitrate pamoate embonate xinafoate hyclate
    napsilate oxalate malate cipionate enantate decanoate undecylate palmitate stearate benzoate furoate""",
    "es": """clorhidrato diclorhidrato bromhidrato maleato fumarato hemifumarato besilato mesilato dimesilato
    tosilato succinato tartrato hemitartrato bitartrato citrato sulfato hemisulfato fosfato acetato valerato
    propionato dipropionato butirato lactato gluconato nitrato pamoato embonato xinafoato hiclato oxalato malato
    cipionato enantato decanoato undecilato palmitato estearato benzoato furoato""",
}
TRAIL = {
    "fr": """sodique disodique trisodique potassique dipotassique calcique magnesique monohydrate dihydrate
    trihydrate hemihydrate sesquihydrate pentahydrate anhydre base micronise racemique de sodium potassium
    calcium magnesium""",
    "es": """sodico sodica disodico disodica trisodico potasico potasica dipotasico calcico calcica magnesico
    magnesica monohidrato dihidrato trihidrato hemihidrato sesquihidrato pentahidrato anhidro anhidra
    hidratado hidratada base micronizado racemico de sodio potasio calcio magnesio""",
}
EN_ENDINGS = {
    "fr": [(r"ine$", "in"), (r"ide$", "id"), (r"one$", "on"), (r"ol$", "ol"), (r"e$", "")],
    "es": [(r"ina$", "ine"), (r"ino$", "ine"), (r"ina$", "in"), (r"ino$", "in"), (r"ona$", "one"), (r"ida$", "ide"), (r"eno$", "ene"),
           (r"ano$", "ane"), (r"ilo$", "yl"), (r"ato$", "ate"), (r"ico$", "ic"), (r"(?<=[bcdfgkmptvx])o$", ""),
           (r"(?<=[bcdfgkmptvx])e$", ""), (r"a$", "")],
}


def deacc(s):
    return "".join(c for c in unicodedata.normalize("NFKD", str(s)) if not unicodedata.combining(c))


def candidates(name, lang):
    n = deacc(name).lower()
    n = re.sub(r"\(.*?\)", " ", n)
    n = re.sub(r"[^a-z0-9' \-]", " ", n)
    n = re.sub(r"\s+", " ", n).strip()
    salts, trail = set(SALTS[lang].split()), set(TRAIL[lang].split())
    out = [n]
    toks = n.replace("d'", "de ").split()
    if len(toks) > 2 and toks[0] in salts and toks[1] == "de":
        toks = toks[2:]
    while len(toks) > 1 and toks[-1] in trail:
        toks.pop()
    core = " ".join(toks)
    out.append(core)
    for acid in ("acide ", "acido "):
        if core.startswith(acid):
            out += [core[len(acid):] + " acid", core[len(acid):]]
    for a, b in EN_ENDINGS[lang]:
        out.append(" ".join(re.sub(a, b, w) for w in core.split()))
    seen = []
    for o in out:
        if o and o not in seen:
            seen.append(o)
    return seen


def map_names(names, lang, index, parents, wd_langs, hits_file, manual=None):
    """-> {name: ChEMBL parent id}"""
    manual = manual or {}
    hit = {}
    for n in names:
        cid = index.lookup(manual[n]) if n in manual else None
        for c in ([] if cid else candidates(n, lang)):
            cid = index.lookup(c)
            if cid:
                break
        hit[n] = cid
    todo = [n for n, c in hit.items() if not c]
    if todo:
        hits_file = Path(hits_file)
        if not hits_file.exists():
            labels = sorted({v for n in todo for v in (n.strip().lower(), n.strip().capitalize(), n.strip())})
            lf = hits_file.with_suffix(".names.txt")
            lf.write_text("\n".join(labels), encoding="utf-8")
            subprocess.run([sys.executable, str(HERE / "wikidata_labels.py"), str(lf), ",".join(wd_langs),
                            str(hits_file)], check=True)
        h = pd.read_csv(hits_file, dtype=str)
        ik2cid = {p["inchikey"]: cid for cid, p in parents.items() if p.get("inchikey")}
        by_lab = {}
        for lab, g in h.groupby("lab"):
            cids = {ik2cid.get(ik) for ik in g.ik.dropna()} - {None}
            if len(cids) == 1:
                by_lab[lab.lower()] = next(iter(cids))
        for n in todo:
            hit[n] = by_lab.get(n.strip().lower())
    return hit
