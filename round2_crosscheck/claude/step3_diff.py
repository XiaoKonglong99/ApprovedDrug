import pandas as pd, re, itertools
A=["FDA","EMA","NMPA","PMDA"]
a=pd.read_pickle("allrows.pkl")
# collapse same-source multi rows: keep first
a=a.drop_duplicates(["ent","src"])
W={s:a[a.src==s].set_index("ent") for s in ["CL","GE","DS"]}
ents=sorted(set(a.ent))
def st(v):
    if re.fullmatch(r"\d{4}",v): return "Y"
    if v.startswith("Approved"): return "A"
    if v=="NA": return "-"
    return "P"
rows=[]
for e in ents:
    r={"ent":e}
    for s,w in W.items():
        if e in w.index:
            x=w.loc[e]; r[s+"_name"]=x.Name; r[s+"_ik"]=x.ik; r[s+"_smi"]=x.SMILES
            for c in A: r[s+"_"+c]=x[c]
    rows.append(r)
m=pd.DataFrame(rows)
m.to_pickle("wide.pkl")
# SMILES agreement
print("== SMILES agreement (pairwise, entities in both)")
for s1,s2 in itertools.combinations(["CL","GE","DS"],2):
    b=m[m[s1+"_ik"].notna()&m[s2+"_ik"].notna()]
    full=(b[s1+"_ik"]==b[s2+"_ik"]).mean(); con=(b[s1+"_ik"].str[:14]==b[s2+"_ik"].str[:14]).mean()
    print(s1,s2,len(b),f"identical {full:.3f}  same skeleton(ignore stereo) {con:.3f}  n_diff_skeleton {(b[s1+'_ik'].str[:14]!=b[s2+'_ik'].str[:14]).sum()}")
print("== approval status agreement per agency (approved = year or 'Approved')")
for c in A:
    for s1,s2 in itertools.combinations(["CL","GE","DS"],2):
        b=m[m[s1+"_ik"].notna()&m[s2+"_ik"].notna()]
        x=b[s1+"_"+c].map(st).isin(["Y","A"]); y=b[s2+"_"+c].map(st).isin(["Y","A"])
        yy=b[s1+"_"+c].map(st).eq("Y")&b[s2+"_"+c].map(st).eq("Y")
        d=(b[s1+"_"+c].where(yy).astype(float)-b[s2+"_"+c].where(yy).astype(float)).abs()
        print(f"{c} {s1}-{s2}: both appr {int((x&y).sum())}, only {s1} {int((x&~y).sum())}, only {s2} {int((~x&y).sum())} | both years {int(yy.sum())}, same year {int((d==0).sum())}, diff>1y {int((d>1).sum())}")
