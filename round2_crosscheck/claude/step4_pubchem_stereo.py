import pandas as pd, urllib.request, urllib.parse, json, time
m=pd.read_pickle("wide.pkl")
b=m[m.DS_ik.notna() & (m.CL_ik.notna()|m.GE_ik.notna())].copy()
b["ref"]=b.CL_ik.fillna(b.GE_ik); b["refname"]=b.CL_name.fillna(b.GE_name)
b=b[(b.ref!=b.DS_ik)]
out=[]
for _,r in b.iterrows():
    pk=""
    for q in [r.DS_name, r.refname]:
        try:
            u="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"+urllib.parse.quote(q)+"/property/InChIKey/JSON"
            pk=json.load(urllib.request.urlopen(u,timeout=20))["PropertyTable"]["Properties"][0]["InChIKey"]; break
        except Exception: pass
        time.sleep(0.25)
    time.sleep(0.25)
    def t(ik): return "FULL" if ik==pk else ("SKEL" if ik[:14]==pk[:14] else "no")
    out.append(dict(name=r.refname,pubchem=pk,CLGE=t(r.ref) if pk else "?",DS=t(r.DS_ik) if pk else "?",skel_same=r.ref[:14]==r.DS_ik[:14]))
o=pd.DataFrame(out); o.to_csv("smiles_conflicts_pubchem.csv",index=False)
print(len(o)); print(pd.crosstab(o.CLGE,o.DS))
