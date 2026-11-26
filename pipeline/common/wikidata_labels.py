"""Look up drug names in Wikidata by exact label/alias in a given language.

usage: python3 wikidata_labels.py names.txt LANGS out.csv      (LANGS e.g. "zh,zh-hans,zh-cn" or "pt,pt-br")
Output columns: lab, item, itemLabel, smi (P233), ismi (P2017), ik (P235)
"""
import csv, json, sys, time, urllib.parse, urllib.request

names_file, langs, out = sys.argv[1], sys.argv[2].split(","), sys.argv[3]
names = [l.strip() for l in open(names_file, encoding="utf-8") if l.strip() and '"' not in l and "\\" not in l]
Q = """SELECT ?lab ?item ?itemLabel ?smi ?ismi ?ik WHERE {
 VALUES ?lab { %s }
 { ?item rdfs:label ?lab } UNION { ?item skos:altLabel ?lab }
 OPTIONAL{?item wdt:P233 ?smi} OPTIONAL{?item wdt:P2017 ?ismi} OPTIONAL{?item wdt:P235 ?ik}
 SERVICE wikibase:label { bd:serviceParam wikibase:language "en". } }"""
rows, B = [], 80
for i in range(0, len(names), B):
    vals = " ".join(f'"{n}"@{t}' for n in names[i:i + B] for t in langs)
    data = urllib.parse.urlencode({"query": Q % vals}).encode()
    for attempt in range(5):
        try:
            req = urllib.request.Request("https://query.wikidata.org/sparql", data=data, headers={
                "User-Agent": "drug-table-crosscheck/1.0 (research script)",
                "Accept": "application/sparql-results+json",
                "Content-Type": "application/x-www-form-urlencoded"})
            res = json.load(urllib.request.urlopen(req, timeout=90))
            break
        except Exception as e:
            print(f"batch {i}: retry {attempt + 1} ({e})", flush=True)
            time.sleep(5 * (attempt + 1))
            res = None
    if res:
        for b in res["results"]["bindings"]:
            rows.append({k: v["value"] for k, v in b.items()})
    print(f"{min(i + B, len(names))}/{len(names)} names, {len(rows)} hits", flush=True)
    time.sleep(0.5)
with open(out, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["lab", "item", "itemLabel", "smi", "ismi", "ik"])
    w.writeheader()
    w.writerows(rows)
