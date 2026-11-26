"""Download every medicine in Spain's AEMPS CIMA register (REST API, 200 per page) -> raw/cima_es.json.gz"""
import gzip, json, time, urllib.request
from pathlib import Path
RAW = Path(__file__).resolve().parent / "raw"
RAW.mkdir(exist_ok=True)
out, page = [], 1
while True:
    for attempt in range(5):
        try:
            req = urllib.request.Request(f"https://cima.aemps.es/cima/rest/medicamentos?pagina={page}",
                                         headers={"Accept": "application/json", "User-Agent": "drug-table-crosscheck/1.0"})
            d = json.load(urllib.request.urlopen(req, timeout=60))
            break
        except Exception as e:
            print("retry", page, e, flush=True); time.sleep(3 * (attempt + 1))
    rs = d.get("resultados") or []
    for r in rs:
        out.append({k: r.get(k) for k in ("nregistro", "nombre", "estado", "ema", "generico", "vtm", "comerc")})
    if page % 20 == 0:
        print(page, len(out), "/", d.get("totalFilas"), flush=True)
    if len(rs) < d.get("tamanioPagina", 200):
        break
    page += 1
    time.sleep(0.2)
json.dump(out, gzip.open(RAW / "cima_es.json.gz", "wt"))
print("records:", len(out))
