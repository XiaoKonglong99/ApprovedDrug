# Small-molecule drugs by regulator (FDA / EMA / NMPA / PMDA)

```bash
pip install requests pandas openpyxl lxml
python3 collect_drugs.py            # reuses downloads cached in data/raw/
python3 collect_drugs.py --refresh  # re-download everything (~15 min; ClinicalTrials.gov is the slow part)
```

## Outputs (`output/`)

| File | Content |
|---|---|
| `small_molecule_drugs_by_agency_Claude.csv` | `Name, SMILES, FDA, EMA, NMPA, PMDA` |
| `small_molecule_drugs_by_agency_detailed_Claude.csv` | same + `ChEMBL_ID, InChIKey, MW, ChEMBL_max_phase` |

Each agency cell holds the strongest evidence found:

| Value | Meaning |
|---|---|
| `2019` | approved; year of first approval at that agency |
| `Approved` | approved/marketed there, but the source has no reliable first-approval year |
| `Phase 4` … `Early Phase 1` | highest phase of a ClinicalTrials.gov trial with a site in that region (`Phase 1/2`, `Phase 2/3` for combined phases). `Phase 4` = post-marketing study, i.e. very likely marketed there |
| `NA` | no record found |

**Scope:** ChEMBL `molecule_type = "Small molecule"` with any clinical phase. Proteins, antibodies, oligonucleotides (RNA/DNA), cells, and gene therapies are excluded. ChEMBL classifies cyclic peptides such as ciclosporin and octreotide as proteins, so they are also excluded. Some larger peptides such as daptomycin and vancomycin are classed as small molecules and stay in. Salt and hydrate forms are merged into the parent molecule, and the SMILES is the parent (salt-free) structure.

## Sources and rules

| Agency | Year comes from | `Approved` (no year) comes from |
|---|---|---|
| FDA | Drugs@FDA (openFDA bulk export): earliest original NDA/ANDA approval | — |
| EMA | KEGG br08329, EU new active substances since 1995; fallback is the EMA centralised MA date for single-ingredient, non-generic products | EMA centralised medicines list + EMA Article 57 database (all national EU/EEA authorisations) |
| PMDA | KEGG br08318, Japan new active ingredients since 2000 (incl. emergency approvals) | KEGG br08301 (all drugs marketed in Japan) |
| NMPA | NMPA English-site approval news (novel drugs, ~2021–2025) | `data/manual/nmpa_approvals.csv`, if you add it |

- **Trials:** ClinicalTrials.gov phase 0.5–4 interventional studies (withdrawn excluded). Site countries map to regions: US/Puerto Rico → FDA, EU/EEA → EMA, mainland China → NMPA, Japan → PMDA. Trials are linked to molecules through ChEMBL's curated NCT references and exact name/synonym matching.
- **Combination products:** a new combination only gives a year to its *new* component. An old partner, such as ritonavir in a nirmatrelvir/ritonavir pack, gets `Approved` instead.

## Known limitations

1. **NMPA approvals are thin** (~36 molecules with a year). NMPA and CDE databases use anti-bot protection, so they can't be scraped. China status mostly comes from China trial sites: `Phase 4` usually means marketed in China. To add the official list, create `data/manual/nmpa_approvals.csv` (columns `name,year`, English INN) and re-run. The same works for `fda_`, `ema_` and `pmda_approvals.csv`.
2. **EMA and PMDA years start in 1995 and 2000.** Older drugs show `Approved`, not a year.
3. **FDA pre-1982 dates** in Drugs@FDA are sometimes the date of a later re-filed application, not the true first approval.
4. **Name matching** depends on ChEMBL synonyms; ambiguous names are skipped rather than guessed. A few ChEMBL synonym errors can create false trial matches. One example: "PV-10" is a rose bengal product, but ChEMBL lists it under α-PNP.
5. ClinicalTrials.gov under-covers China and Japan. Trials registered only in ChiCTR, CDE or jRCT are missing.
6. **KEGG** data is free for academic use.
