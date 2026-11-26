import pandas as pd, re
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.MolStandardize import rdMolStandardize
from pathlib import Path
R1=Path(__file__).resolve().parents[2]/"round1_agents"
files={"DS":R1/"deepseek"/"drug_approvals_smiles_Deepseek.csv","GE":R1/"gemini"/"small_molecule_drugs_Gemini.csv","CL":R1/"claude"/"small_molecule_drugs_by_agency_detailed_Claude.csv"}
lfc=rdMolStandardize.LargestFragmentChooser()
unch=rdMolStandardize.Uncharger()
def keys(smi):
    m=Chem.MolFromSmiles(smi)
    if m is None: return pd.Series({"ik":"","ik14":"","ha":None,"mw":None,"valid":False})
    try:
        p=unch.uncharge(lfc.choose(m))
    except Exception: p=m
    ik=Chem.MolToInchiKey(p) or ""
    from rdkit.Chem import Descriptors
    return pd.Series({"ik":ik,"ik14":ik[:14],"ha":p.GetNumHeavyAtoms(),"mw":Descriptors.MolWt(p),"valid":True,"nfrag":len(Chem.GetMolFrags(m))})
def norm(n):
    n=n.lower().strip()
    n=re.sub(r"[\s\-_,'’()\[\]]","",n)
    return n
for k,f in files.items():
    d=pd.read_csv(f,dtype=str,keep_default_na=False)
    d=d[["Name","SMILES","FDA","EMA","NMPA","PMDA"]].copy()
    d["nname"]=d.Name.map(norm)
    d=pd.concat([d,d.SMILES.apply(keys)],axis=1)
    print(k,len(d),"invalid smiles",(~d.valid).sum(),"multi-frag",(d.nfrag>1).sum(),"dup ik",d.ik[d.ik!=""].duplicated().sum(),"dup ik14",d.ik14[d.ik14!=""].duplicated().sum())
    d.to_pickle(f"{k}.pkl")
