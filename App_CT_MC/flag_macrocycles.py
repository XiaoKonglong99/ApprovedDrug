"""Step 1: flag macrocycles (MRS >= 12) in the approved-drug CSV.

MRS uses MacrocycleWizard's own definition (A07_Macrocycle_Ring_Size):
number of bonds in the largest SSSR ring, after stripping substituents.
"""
import os
import sys
import pandas as pd
from rdkit import Chem, RDLogger

os.chdir(os.path.dirname(os.path.abspath(__file__)))   # run from anywhere
sys.path.insert(0, "MacrocycleWizard_CP")
from MacrocycleWizard.MolD.base import _get_macrocycle_ring_mol

RDLogger.DisableLog("rdApp.*")
SRC = "[BEST]ApprovedDrug_CT_2026Nov.csv"


def mrs(smiles):
    if pd.isna(smiles) or not str(smiles).strip():
        return pd.NA, "missing SMILES"
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return pd.NA, "unparsable SMILES"
    if mol.GetRingInfo().NumRings() == 0:
        return 0, ""
    try:
        return _get_macrocycle_ring_mol(mol, strip=True).GetNumBonds(), ""
    except Exception as e:
        return pd.NA, f"{type(e).__name__}: {e}"


df = pd.read_csv(SRC, encoding="utf_8_sig", dtype=str, keep_default_na=False)   # keep literal "NA" cells
res = df["SMILES"].map(mrs)
df["MRS"] = [r[0] for r in res]
df["is_MC"] = df["MRS"].apply(lambda x: bool(pd.notna(x) and x >= 12))
errs = [(n, r[1]) for n, r in zip(df["Name"], res) if r[1]]

df.to_csv("[BEST]ApprovedDrug_CT_2026Nov_isMC.csv", index=False)
mc = df[df["is_MC"]].reset_index(drop=True)
mc.to_csv("ApprovedDrug_CT_2026Nov_MC_only.csv", index=False)

print(f"Total rows: {len(df)} | macrocycles: {len(mc)} | MRS errors: {len(errs)}")
for n, e in errs:
    print("  ", n, "->", e)
