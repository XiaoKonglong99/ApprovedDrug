# Macrocycles in the final table

Run from anywhere (needs a Python that can import MacrocycleWizard, which is not included in this repository):

```bash
cp "../[BEST]final_small_molecule_drugs_Claude.csv" "[BEST]ApprovedDrug_CT_2026Nov.csv"
python flag_macrocycles.py        # -> [BEST]ApprovedDrug_CT_2026Nov_isMC.csv, ApprovedDrug_CT_2026Nov_MC_only.csv
python calc_mc_molds_onion.py     # -> ApprovedDrug_CT_2026Nov_MC_MolDs_Onion.csv
python add_route_info.py          # adds Route_Class / Routes / Route_Source to both macrocycle files
```

- `MRS` = MacrocycleWizard's macrocycle ring size (bonds in the largest SSSR ring after stripping substituents);
  `is_MC` = MRS >= 12. 195 of 6,256 molecules are macrocycles.
- The MolDs file has all 151 descriptors (groups A–K) plus the Macrocycle Onion layer SMILES, 0 failures.
- **Route of administration** (macrocycles only), columns inserted after `ANVISA` in
  `ApprovedDrug_CT_2026Nov_MC_only.csv` and `ApprovedDrug_CT_2026Nov_MC_MolDs_Onion.csv`:
  - `Route_Class`: `Oral` (80), `Injection` (83), `Oral + Injection` (24), `Non-oral (other)` (8: only topical,
    ophthalmic, nasal, vaginal… products, e.g. pimecrolimus, natamycin, lypressin).
  - `Routes`: routes found, e.g. `oral; intravenous`. `Route_Source`: the evidence that decided it.
  - Approved drugs: routes of marketed products from Drugs@FDA, EMA Article 57 (a route counts only if
    >= 10% of the substance's EU products report it, because the database has coding errors), France ANSM,
    and Chinese product names (片/胶囊 = oral, 注射 = injection). Trial-only compounds: ClinicalTrials.gov
    intervention descriptions (a route counts if it appears in >= 1/3 of the studies that state one).
    Last resort: ChEMBL route flags.
  - 38 doubtful or unresolved compounds were checked online; corrections and sources are in
    `route_overrides_Claude.csv`. Per-source evidence for every compound: `ApprovedDrug_CT_2026Nov_MC_routes.csv`.
    Online lookups are cached in `route_cache/`.
- Put MacrocycleWizard in `MacrocycleWizard_CP/` next to these scripts (it is git-ignored), or point
  `MACROCYCLEWIZARD_DIR` at a MacrocycleWizard checkout.
