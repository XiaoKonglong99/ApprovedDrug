import pandas as pd, re, numpy as np
A=["FDA","EMA","NMPA","PMDA"]
D={k:pd.read_pickle(f"{k}.pkl").reset_index(drop=True) for k in ["CL","GE","DS"]}
for k,d in D.items(): d["rid"]=k+":"+d.index.astype(str)
# union-find over rows
parent={}
def f(x):
    while parent.setdefault(x,x)!=x: x=parent[x]
    return x
def u(a,b): parent[f(a)]=f(b)
allr=pd.concat([d.assign(src=k) for k,d in D.items()],ignore_index=True)
for r in allr.rid: f(r)
# 1 exact normalized name across files
for n,g in allr.groupby("nname"):
    ids=g.rid.tolist()
    for x in ids[1:]: u(ids[0],x)
# 2 full InChIKey, only if unique within each file (avoid carbonate-style collisions)
uniq={k:set(d.ik[~d.ik.duplicated(keep=False)]) for k,d in D.items()}
for ik,g in allr[allr.ik!=""].groupby("ik"):
    g=g[g.apply(lambda r: r.ik in uniq[r.src],axis=1)]
    ids=g.rid.tolist()
    if len(set(g.src))>1:
        for x in ids[1:]: u(ids[0],x)
allr["ent"]=allr.rid.map(f)
# per entity: refuse entities that contain >1 row from the same source (keep but flag)
cnt=allr.groupby(["ent","src"]).size().unstack(fill_value=0)
print("entities",len(cnt)); print("presence patterns:");
pat=(cnt>0).astype(int).astype(str).agg(lambda r: "".join(f"{c}" if v=="1" else "-" for c,v in zip(r.index,r)),axis=1)
print(pat.value_counts())
print("entities with same-source multi rows:", (cnt>1).any(axis=1).sum())
allr.to_pickle("allrows.pkl")
