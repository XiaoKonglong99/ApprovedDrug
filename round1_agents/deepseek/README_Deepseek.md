# Small-molecule drug approvals and SMILES — FDA, EMA, NMPA, PMDA

**Main file:** `drug_approvals_smiles_Deepseek.csv`
**Generated:** 26 November 2026
**Rows:** 2,633
**Approval records:** FDA 1,900; EMA 618; NMPA 1,563; PMDA 393
**Non-NA cells after trial-phase fill:** FDA 2,025; EMA 1,345; NMPA 1,811; PMDA 755
**Columns:** `Name,SMILES,FDA,EMA,NMPA,PMDA`

## What the table contains

One row per small-molecule active substance for which at least one of the four
agencies has an approval record.

Each agency cell is one of:

| Value | Meaning |
|---|---|
| `2019` (a year) | Approved by that agency; year of approval (see caveats for NMPA). |
| `Phase 3` | Not approved by that agency, but an interventional trial with a site in its jurisdiction reached that phase. |
| `NA` | No approval or trial record found. |

A molecule approved by several agencies still has one row; the approval years
appear in the corresponding agency columns.

## Sources

* **FDA** — Drugs@FDA approval records via the
  [DrugCentral](https://drugcentral.org/download) public PostgreSQL database,
  supplemented with the FDA
  [Novel Drug Approvals 2025/2026](https://www.fda.gov/drugs/novel-drug-approvals-fda/novel-drug-approvals-2026)
  tables.
* **EMA** — the EMA
  [Download medicine data](https://www.ema.europa.eu/en/medicines/download-medicine-data)
  workbook (generated 26 Nov 2026; central marketing-authorisation dates).
* **NMPA** — the public NMPA drug database snapshot mirrored at
  [xprincehf/nmpa-data](https://github.com/xprincehf/nmpa-data)
  (domestic + imported chemical drugs; snapshot ends Aug 2021).
* **PMDA** — DrugCentral's PMDA approval records.
* **SMILES** — DrugCentral structures; missing structures and NMPA/EMA/FDA
  names were resolved with the
  [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) name
  service. SMILES are canonical isomeric SMILES of the parent (largest organic
  fragment); salts and counter-ions are removed.
* **Trial phases** — [ClinicalTrials.gov API v2](https://clinicaltrials.gov/data-api/api).
  For each drug/agency, the highest phase among interventional studies with a
  site in the United States (FDA), EU/EEA (EMA), China (NMPA) or Japan (PMDA)
  was used. Withdrawn, terminated, suspended and unknown-status studies were
  ignored.

## Inclusion/exclusion rules

* Small molecules only: parent molecular weight ≤ 2000 Da, at least one carbon
  and at least three heavy atoms. Proteins, antibodies, vaccines, RNA/DNA
  therapeutics, cells and gene therapies are excluded.
* Fixed-dose combinations are excluded (names containing `and`, `/`, `;`,
  `compound`, or EMA INNs with `;`).
* Salts are collapsed to the parent. Ester/ketal prodrugs are collapsed onto
  the parent **only when the parent itself is present** in the data (e.g.
  abiraterone acetate → abiraterone); otherwise they remain separate.
* Common nutrients/excipients and free amino acids are removed when they are
  NMPA-only records.

## Important caveats

1. **NMPA years are the weakest part of the table.** The public snapshot's
   `批准日期` is often a 2014/2015 bulk re-registration date, and many legacy
   approval numbers were re-issued in 1998/2002. The pipeline therefore uses
   the earliest of the date and the year encoded in the approval number; for
   legacy products this can still be a re-registration year rather than the
   original first approval. NMPA records after August 2021 are not covered by
   this snapshot, so recent Chinese approvals may show a trial phase or `NA`.
2. Approval is counted historically: a product later withdrawn/suspended is
   still shown with its approval year. Cells reading `Approved` mean an
   approval record exists but no reliable year was available (mostly legacy
   FDA products).
3. `NA` means “not found in the public sources used here”, not proof that the
   molecule has never been developed or approved.
4. ClinicalTrials.gov is incomplete for China and Japan (ChiCTR and jRCT are
   separate registries), so NMPA/PMDA phase cells may understate activity.
   For runtime reasons the pipeline scans the first 500 studies returned for
   each 8-drug query; for very common drugs a later page could contain a higher
   phase.
5. Trial phase is a proxy for regulatory stage; a Phase 3 trial in a region is
   not the same as a regulatory submission.
6. The table starts from approved drugs. Molecules that are only
   investigational everywhere are not included.

## Reproducing

```bash
# Python 3.10+, RDKit, pandas, requests, openpyxl, lxml
pip install pandas requests rdkit openpyxl lxml pg8000

# Use the cached source tables in data/ (fast; still queries CT.gov once)
python collect_drugs_Deepseek.py --out drug_approvals_smiles_Deepseek.csv

# Re-download DrugCentral/EMA/NMPA/FDA, then rebuild
python collect_drugs_Deepseek.py --refresh --out drug_approvals_smiles_Deepseek.csv

# Approvals only, no ClinicalTrials.gov queries
python collect_drugs_Deepseek.py --skip-trials --out approvals_only.csv
```

The script caches PubChem responses in `data/pubchem_results.csv` and
ClinicalTrials.gov batch responses in `cache/ct/`, so interrupted runs resume.
