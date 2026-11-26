"""Step 2: all MolDs + Onion layers for the macrocycles, via MacrocycleWizard's CalculateAll.

Workaround (no package files modified): base._get_R_fragment_from_macrocycle calls
MolFragmentToSmiles(kekuleSmiles=True) on the substituent fragment, which raises
KekulizeException when the fragment cuts a fused aromatic system (e.g. Rifaximin).
Only on that exception, retry on a copy of the molecule kekulized as a whole first.
Molecules that already worked follow the original code path unchanged.
"""
import os
import sys
import pandas as pd
from rdkit import Chem

os.chdir(os.path.dirname(os.path.abspath(__file__)))   # run from anywhere
sys.path.insert(0, os.environ.get("MACROCYCLEWIZARD_DIR", "MacrocycleWizard_CP"))   # package location
import MacrocycleWizard.MolD.base as mw_base
from CalculateAll import calculate_molds

_orig_get_R = mw_base._get_R_fragment_from_macrocycle


def _get_R_with_fallback(mol):
    try:
        return _orig_get_R(mol)
    except Chem.rdchem.KekulizeException:
        kmol = Chem.Mol(mol)
        Chem.Kekulize(kmol, clearAromaticFlags=True)
        print(f"  [fallback] kekulized whole mol for R-fragments: {Chem.MolToSmiles(mol)}")
        return _orig_get_R(kmol)


mw_base._get_R_fragment_from_macrocycle = _get_R_with_fallback

mc = pd.read_csv("ApprovedDrug_CT_2026Nov_MC_only.csv", dtype=str, keep_default_na=False)
mc.drop(columns=["MRS", "is_MC"]).to_csv("_mc_input.tmp.csv", index=False)
calculate_molds("_mc_input.tmp.csv", "ApprovedDrug_CT_2026Nov_MC_MolDs_Onion.csv",
                smiles_col="SMILES", id_col="Name")
os.remove("_mc_input.tmp.csv")

# CalculateAll reads the table with pandas defaults, which turns the literal "NA" cells into blanks;
# put the original text columns back (same row order).
OUT = "ApprovedDrug_CT_2026Nov_MC_MolDs_Onion.csv"
res = pd.read_csv(OUT)
text_cols = [c for c in mc.columns if c not in ("MRS", "is_MC")]
res[text_cols] = mc[text_cols].values
res.to_csv(OUT, index=False)
