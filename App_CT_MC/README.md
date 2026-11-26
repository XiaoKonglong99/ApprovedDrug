# Macrocycles in the final table

Run from anywhere (needs a Python that can import MacrocycleWizard, which is not included in this repository):

```bash
cp "../[BEST]final_small_molecule_drugs_Claude.csv" "[BEST]ApprovedDrug_CT_2026Nov.csv"
python flag_macrocycles.py        # -> [BEST]ApprovedDrug_CT_2026Nov_isMC.csv, ApprovedDrug_CT_2026Nov_MC_only.csv
python calc_mc_molds_onion.py     # -> ApprovedDrug_CT_2026Nov_MC_MolDs_Onion.csv
```

- `MRS` = MacrocycleWizard's macrocycle ring size (bonds in the largest SSSR ring after stripping substituents);
  `is_MC` = MRS >= 12. 196 of 6,260 molecules are macrocycles.
- The MolDs file has all 151 descriptors (groups A–K) plus the Macrocycle Onion layer SMILES, 0 failures.
- Put MacrocycleWizard in `MacrocycleWizard_CP/` next to these scripts (it is git-ignored).
