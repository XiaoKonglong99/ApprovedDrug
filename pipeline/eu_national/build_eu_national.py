"""EU national authorisations (outside the EMA centralised procedure).

Sources
  raw/CIS_bdpm.txt, raw/CIS_COMPO_bdpm.txt  France, ANSM Base de données publique des médicaments
                                             (date d'AMM + procedure: nationale / décentralisée /
                                             reconnaissance mutuelle / centralisée)
  raw/cima_es.json.gz                        Spain, AEMPS CIMA register (authorisation date, 'ema' flag)
  data/raw/ema_article57.xlsx (Claude round-1 cache) EMA Article 57 database: every medicine authorised in each
                                             EU/EEA country (no dates) -> presence only
Output
  eu_national_rows_Claude.csv   Name,SMILES,EU_National,FR_year,ES_year,A57_countries  (source "EU")
  eu_national_ingredients_Claude.csv   audit of the name mapping
"""
import gzip, json, re, sys
from collections import defaultdict
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "round1_agents" / "claude"))
sys.path.insert(0, str(ROOT / "pipeline" / "common"))
import collect_drugs as cd  # noqa: E402
from name_mapper import map_names  # noqa: E402

NATIONAL_FR = {"Procédure nationale", "Procédure décentralisée", "Procédure de reconnaissance mutuelle"}


def download_raw():
    """Fetch the national registers if they are not cached in raw/ yet."""
    import subprocess, urllib.request
    (HERE / "raw").mkdir(exist_ok=True)
    for f in ("CIS_bdpm.txt", "CIS_COMPO_bdpm.txt"):
        if not (HERE / "raw" / f).exists():
            urllib.request.urlretrieve(f"https://base-donnees-publique.medicaments.gouv.fr/download/file/{f}", HERE / "raw" / f)
    if not (HERE / "raw" / "cima_es.json.gz").exists():
        subprocess.run([sys.executable, str(HERE / "fetch_cima.py")], check=True)


def france():
    b = pd.read_csv(HERE / "raw" / "CIS_bdpm.txt", sep="\t", header=None, encoding="latin1", dtype=str,
                    quoting=3, on_bad_lines="skip")
    b = b[b[5].isin(NATIONAL_FR)]
    b["year"] = pd.to_datetime(b[7], format="%d/%m/%Y", errors="coerce").dt.year
    c = pd.read_csv(HERE / "raw" / "CIS_COMPO_bdpm.txt", sep="\t", header=None, encoding="latin1", dtype=str,
                    quoting=3, on_bad_lines="skip")
    c.columns = ["cis", "form", "code", "name", "dose", "ref", "nature", "link"]
    # prefer the 'fraction thérapeutique' (active moiety) over the salt ('substance active') of each link
    ft = c[c.nature == "FT"].drop_duplicates(["cis", "link"]).set_index(["cis", "link"]).name
    sa = c[c.nature == "SA"].drop_duplicates(["cis", "link"]).set_index(["cis", "link"]).name
    comp = ft.combine_first(sa).reset_index()
    rows = comp.merge(b[[0, "year"]].rename(columns={0: "cis"}), on="cis")
    return rows[["name", "year"]].dropna()


def spain():
    d = json.load(gzip.open(HERE / "raw" / "cima_es.json.gz"))
    rows = []
    for r in d:
        if r.get("ema"):
            continue
        aut = (r.get("estado") or {}).get("aut")
        vtm = (r.get("vtm") or {}).get("nombre")
        if aut is None or not vtm:
            continue
        y = pd.Timestamp(aut, unit="ms").year
        for p in re.split(r"\s\+\s|,\s| y ", vtm):
            rows.append((p.strip(), y))
    return pd.DataFrame(rows, columns=["name", "year"])


def article57():
    path = cd.cached_download("ema_article57.xlsx", "https://www.ema.europa.eu/en/documents/other/article-57-product-data_en.xlsx")
    x = pd.read_excel(path, header=None, skiprows=20, usecols=[1, 3],
                      names=["substance", "country"], dtype=str).dropna()
    x = x[x.country != "European Union"]                       # centrally authorised products
    rows = []
    for s, ctry in x.drop_duplicates().itertuples(index=False):
        for part in re.split(r"\s*\|\s*|,\s*|\s\+\s|/", s):
            if 2 < len(part) < 80:
                rows.append((part.strip(), ctry))
    return pd.DataFrame(rows, columns=["name", "country"]).drop_duplicates()


def main():
    download_raw()
    parents, child_of, index = cd.load_chembl()
    fr, es, a57 = france(), spain(), article57()
    maps = {}
    maps["fr"] = map_names(fr.name.str.lower().unique(), "fr", index, parents, ["fr"],
                           HERE / "wikidata_fr_hits.csv")
    maps["es"] = map_names(es.name.str.lower().unique(), "es", index, parents, ["es"],
                           HERE / "wikidata_es_hits.csv")
    a57_hit = {n: index.lookup(n) for n in a57.name.unique()}
    fr["cid"] = fr.name.str.lower().map(maps["fr"])
    es["cid"] = es.name.str.lower().map(maps["es"])
    a57["cid"] = a57.name.map(a57_hit)

    audit = pd.concat([
        fr.groupby(fr.name.str.lower()).agg(n=("year", "size"), cid=("cid", "first")).assign(src="FR"),
        es.groupby(es.name.str.lower()).agg(n=("year", "size"), cid=("cid", "first")).assign(src="ES"),
    ]).reset_index().rename(columns={"name": "ingredient"})
    audit.sort_values("n", ascending=False).to_csv(HERE / "eu_national_ingredients_Claude.csv", index=False)

    fr_y = fr.dropna(subset=["cid"]).groupby("cid").year.min()
    es_y = es.dropna(subset=["cid"]).groupby("cid").year.min()
    a57_c = a57.dropna(subset=["cid"]).groupby("cid").country.apply(lambda s: "; ".join(sorted(set(s))))
    out = []
    for cid in set(fr_y.index) | set(es_y.index) | set(a57_c.index):
        p = parents[cid]
        ys = [int(v[cid]) for v in (fr_y, es_y) if cid in v.index]
        out.append({"Name": (p["name"] or cid).title(), "SMILES": p["smiles"],
                    "EU_National": str(min(ys)) if ys else "Approved",
                    "FR_year": str(int(fr_y[cid])) if cid in fr_y.index else "",
                    "ES_year": str(int(es_y[cid])) if cid in es_y.index else "",
                    "A57_countries": a57_c.get(cid, ""), "ChEMBL_ID": cid})
    o = pd.DataFrame(out)
    o.to_csv(HERE / "eu_national_rows_Claude.csv", index=False)
    for k, df in (("France", fr), ("Spain", es)):
        print(f"{k}: {len(df)} product-ingredient records, matched {df.cid.notna().mean():.1%}, "
              f"{df.cid.nunique()} molecules")
    print(f"Article 57: {a57.name.nunique()} substance names, {a57.cid.nunique()} molecules")
    print(f"EU_National rows: {len(o)}; with year {(o.EU_National != 'Approved').sum()}")


if __name__ == "__main__":
    main()
