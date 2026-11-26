# Cross-check of the three AI-generated drug-approval tables

**Prepared:** 26 November 2026  
**Inputs:** `drug_approvals_smiles_Deepseek.csv`, `small_molecule_drugs_by_agency_detailed_Claude.csv`, `small_molecule_drugs_Gemini.csv`  
**Final table:** `small_molecule_drug_approvals_FINAL.csv` (6,183 molecules)

---

## 1. Bottom line

The three files are not three versions of the same table. After name/salt normalisation they contain **6,419 distinct molecule names**, but only **1,828 (28.5%) appear in all three**. They disagree on coverage, on what “EMA” means, on NMPA dates, and on PMDA dates.

The reconciled table contains **6,183 unique molecules**:

| | Molecules |
|---|---:|
| Total rows | 6,183 |
| Approved by ≥1 of FDA / EMA / NMPA / PMDA | 3,512 |
| Clinical-stage only (no approval found) | 2,671 |
| Approved by all four | 817 |
| Approved by ≥2 agencies | 2,284 |

**Ranking: Claude best, Gemini second, DeepSeek worst.** The detailed reasons are in §11. In short: Claude is the most accurate and best documented, but its NMPA coverage is almost empty; Gemini is the broadest all-round table but its PMDA years are unreliable and it has no documentation; DeepSeek has the best NMPA status coverage but systematically misleading NMPA/PMDA years and includes biologics in a “small molecule” table.

---

## 2. What was checked

1. **All rows and all four agency cells** were parsed and compared; SMILES were canonicalised with RDKit and compared by InChIKey.
2. **Names were normalised** (case, punctuation, dosage forms, salt/hydrate suffixes) and mapped through ChEMBL and DrugCentral synonyms. Molecules were merged by full InChIKey; conflicting duplicates were resolved against PubChem.
3. **Independent reference layers were rebuilt from the official bulk data already in the workspace:**
   - FDA: openFDA Drugs@FDA bulk export (25,496 applications with an ORIG approval), earliest original NDA/ANDA approval.
   - EMA: EMA “Download medicine data” workbook (2,365 human products) for centralised approval dates; EMA Article 57 database (164,785 product rows, 11,163 active substances) for national EU/EEA authorisations.
   - PMDA: KEGG BRITE br08318 (975 new active ingredients, 2000–2026) and br08301 (2,534 drugs marketed in Japan).
   - NMPA: public NMPA snapshot (100,717 chemical-drug records, domestic + imported) and 244 NMPA English news items (2021–2025).
   - Cross-check: DrugCentral approval dates.
4. **Twenty representative disagreements were checked online** against EMA EPAR pages, Drugs@FDA, PMDA labels/review reports, PubChem and Chinese launch sources. The verified cases are in `verified_discrepancies_FINAL.csv`.

The reference layers are not perfect either: openFDA has occasional missing/duplicated applications, Article 57 has name variants, KEGG starts in 2000, and the NMPA snapshot ends in 2020/21. They are nevertheless independent of the agents’ own summaries and are sufficient to expose the systematic differences below.

---

## 3. File overview

| Agent | Rows | FDA non-NA | EMA non-NA | NMPA non-NA | PMDA non-NA | FDA years | EMA years | NMPA years | PMDA years |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DeepSeek | 2,633 | 2,025 | 1,345 | 1,811 | 755 | 1,752 | 618 | 1,563 | 375 |
| Claude | 5,607 | 4,511 | 3,830 | 1,519 | 1,918 | 1,796 | 423 | 32 | 402 |
| Gemini | 4,997 | 4,399 | 3,022 | 1,433 | 1,397 | 1,884 | 667 | 4 | 717 |
| **Final** | **6,183** | **4,886** | **4,073** | **2,563** | **2,142** | **2,035** | **770** | **386** | **688** |

Non-NA includes years, “Approved” without a year, and clinical phases. The final table is broader than any single input on all four agencies.

### Overlap

| Pair / set | Names in common |
|---|---:|
| DeepSeek ∩ Claude | 1,938 |
| DeepSeek ∩ Gemini | 1,865 |
| Claude ∩ Gemini | 4,826 |
| All three | 1,828 |
| Union | 6,419 |

---

## 4. FDA — mostly agree, but old-drug dates differ

Where all three give a year, agreement is high: DeepSeek–Claude 97.4% exact, DeepSeek–Gemini 96.5%, Claude–Gemini 98.7%.

The disagreement is in **which approval** is being dated:

| Drug | DeepSeek | Claude | Gemini | Verified | Who was right |
|---|---|---|---|---|---|
| Digoxin | 1975 | 1954 | 1954 | **1954** (Lanoxin NDA 009330, 16 Nov 1954) | Claude + Gemini |
| Drospirenone | 2019 | 2001 | 2001 | **2001** (Yasmin NDA 021098, 11 May 2001) | Claude + Gemini |
| Ranitidine | 1977 | 1983 | 1983 | **1983** (Zantac NDA 018703, 9 Jun 1983) | Claude + Gemini |
| Iobenguane I-123 | 2008 | 1994 | 1994 | **2008** (AdreView NDA 022290, 19 Sep 2008) | DeepSeek |

DeepSeek uses DrugCentral `fda_date`, which is often a later approval/supplement rather than the first approval. Claude and Gemini use openFDA/Drugs@FDA and are more reliable for old drugs.

Status recall/precision against openFDA + DrugCentral:

| Agent | Recall | Precision |
|---|---:|---:|
| DeepSeek | 97.8% | 93.1% |
| Claude | 94.4% | 98.5% |
| Gemini | 99.1% | 86.8% |

Claude has the fewest false positives, Gemini the fewest false negatives, DeepSeek the most false positives.

---

## 5. EMA — three different definitions

The EMA column is ambiguous, and the three agents interpreted it differently:

- **DeepSeek:** only EMA centralised approvals with a date. 618 years, **no** “Approved” values. It misses national EU approvals almost completely (recall 39.5%).
- **Claude:** centralised approvals **plus** national EU/EEA authorisations from Article 57. 2,087 approved, but only 423 with a year. Recall 94.6%, precision 86.6%.
- **Gemini:** 1,247 approved, 667 with a year, also largely EU-wide. Recall 69.3%, precision 92.0%.

Year accuracy where a reference exists:

| Agent | Exact | Within 1 year | >2 years off |
|---|---:|---:|---:|
| DeepSeek | 93.6% | 95.1% | 25 |
| Claude | 96.1% | 97.8% | 9 |
| Gemini | 97.0% | 97.1% | 17 |

Verified EMA disagreements:

| Drug | DeepSeek | Claude | Gemini | Verified | Who was right |
|---|---|---|---|---|---|
| Opicapone (Ongentys) | 2016 | 2022 | 2016 | **2016** (EC decision 24 Jun 2016) | DeepSeek + Gemini |
| Eliglustat (Cerdelga) | 2001 | 2015 | 2015 | **2015** (authorised 19 Jan 2015) | Claude + Gemini |
| Maribavir (Livtencity) | 2013 | 2022 | 2022 | **2022** (authorised 9 Nov 2022) | Claude + Gemini |

DeepSeek’s EMA years are accurate but its coverage is the weakest. Claude’s “Approved” values are mostly true EU national approvals, not EMA centralised approvals — a definitional issue that should be stated, not a pure error.

In the final table, EMA years are centralised approval dates; “Approved” means an EU national authorisation (Article 57) where no centralised date was found.

---

## 6. NMPA — the largest and most consequential divergence

This is where the three files are furthest apart:

| Agent | NMPA approved (year or “Approved”) | NMPA years | NMPA recall vs snapshot |
|---|---:|---:|---:|
| DeepSeek | 1,563 | 1,563 | 97.4% |
| Claude | 36 | 32 | 0.3% |
| Gemini | 644 | 4 | 40.9% |

DeepSeek has by far the best NMPA **status** coverage. But its **years** are mostly not first-approval years. Of its 1,563 NMPA years:

- **527 (33.7%) are exactly 2002**;
- **929 (59.4%) fall in 1998–2002**.

That is the fingerprint of China’s 1998–2002 re-registration/numbering exercise, not of first approval. The NMPA snapshot’s `批准日期` is often a 2014/2015 bulk re-registration date, and the approval-number year (`国药准字H2002…`) is the year the number was issued, which for legacy products is also a re-registration year.

Verified examples:

| Drug | DeepSeek | Claude | Gemini | Actual first China approval | Comment |
|---|---|---|---|---|---|
| Atorvastatin | 1999 | Phase 4 | Approved | **1999** | DeepSeek correct; Claude/Gemini only status/phase |
| Imatinib | 2013 | Phase 4 | Approved | **2002** | All three miss the year; DeepSeek’s 2013 is a re-registration |
| Sorafenib | 2020 | Phase 4 | Approved | **2006** | DeepSeek’s 2020 is a later domestic generic |
| Abacavir | 2015 | NA | NA | **registered 2015, not marketed** | DeepSeek catches the import registration; the other two miss it |

Claude’s NMPA column is accurate for the few novel drugs it covers (NMPA English news, ~2021–2025) but misses roughly 1,300 approved molecules. Gemini’s 640 “Approved (Year N/A)” values are a reasonable safe compromise, but still incomplete.

The final table therefore uses NMPA **status** broadly, but only reports a **year** when it is from NMPA news or when the snapshot year is within 10 years of the molecule’s global first approval. Everything else is “Approved” without a year. This is deliberately conservative: it loses some correct old years (e.g. atorvastatin 1999 is kept, but clopidogrel/entecavir-style cases may become “Approved”) in exchange for not publishing thousands of re-registration years as first approvals.

---

## 7. PMDA — Claude is clearly the best, Gemini the weakest

| Agent | PMDA approved | PMDA years | Year accuracy vs KEGG | Recall | Precision |
|---|---:|---:|---:|---:|---:|
| DeepSeek | 393 | 375 | 88.6% | 33.7% | 96.2% |
| Claude | 1,382 | 402 | **99.7%** | 92.9% | 85.2% |
| Gemini | 812 | 717 | 80.3% | 60.5% | 85.2% |

Verified PMDA disagreements:

| Drug | DeepSeek | Claude | Gemini | Verified first Japan approval | Who was right |
|---|---|---|---|---|---|
| Imatinib | Phase 4 | 2001 | 2012 | **2001** | Claude |
| Bortezomib | 2019 | 2006 | 2011 | **2006** | Claude |
| Everolimus | 2019 | 2007 | 2011 | **2007** | Claude |
| Ambrisentan | 2021 | 2010 | 2010 | **2010** | Claude + Gemini |
| Metformin | 1961 | Approved | 2010 | **1961** | DeepSeek |

DeepSeek again dates later approvals: DrugCentral’s `pmda_date` is often a generic/new-indication approval, not the first new-ingredient approval. Gemini’s PMDA years look like later product/generic approvals (imatinib 2012, metformin 2010). Claude’s KEGG-derived years are almost exact, but its “Approved” cells overstate the number of *approved* products (KEGG’s marketed-drug list is broader than new-ingredient approvals); this is a definitional issue, not a year error.

---

## 8. SMILES and scope

### Structure agreement

For the 1,828 molecules present in all three files:

- **1,741 (95.2%)** have identical InChIKeys;
- **62 (3.4%)** share the same skeleton but differ in salt/protonation/stereochemistry;
- **25 (1.4%)** have genuinely different connectivity.

Against ChEMBL structures (a useful but imperfect benchmark):

| Agent | Exact InChIKey | Same skeleton only | Different |
|---|---:|---:|---:|
| DeepSeek | 88.7% | 4.0% | 7.4% |
| Claude | 98.5% | 0.8% | 0.7% |
| Gemini | 98.3% | 1.1% | 0.6% |

ChEMBL is not always right: for alogliptin and aceneuramic acid, DeepSeek matched PubChem while ChEMBL’s synonym-mapped structure did not. But a 7.4% structural mismatch rate is still a major quality problem for DeepSeek.

PubChem spot checks:

- **Eravacycline:** DeepSeek `AKLMFDDQ…` is correct; Claude and Gemini `HLFSMU…` are wrong.
- **Dalfopristin:** Claude and Gemini `SUYRLXYY…` are correct; DeepSeek `KLRSZNO…` is wrong.
- **Apilimod:** Claude and Gemini each contain both the correct and an incorrect structure under the same name; the final table keeps the PubChem-verified one.

### Scope violations

All three include **somapacitan** (a growth-hormone fusion protein) with a fragment SMILES; the final table removes it.

DeepSeek additionally includes clear biologics/RNA products that the prompt explicitly excludes:

- **Interferon alfa-2b** — SMILES MW 377, actual protein ~19.3 kDa; the SMILES is a small synthetic molecule, not interferon.
- **Ropeginterferon alfa-2b** — same problem.
- **Cemdisiran** — a GalNAc-conjugated siRNA, not a small molecule.
- **Elapegademase**, **Pegunigalsidase alfa**, **Somatrem** — enzyme/protein therapeutics.

The final table removes 10 such rows and keeps peptide/small-protein antibiotics (daptomycin, teicoplanin, oritavancin, etc.) only where MW ≤ 2,000 and ChEMBL classifies them as small molecules.

---

## 9. Final table: how it was built

`small_molecule_drug_approvals_FINAL.csv` has exactly the requested columns:

`Name, SMILES, FDA, EMA, NMPA, PMDA`

Cell rules:

| Agency | Year | “Approved” | Phase |
|---|---|---|---|
| FDA | openFDA earliest original approval; else DrugCentral date; else two-agent consensus | agent-reported approval with no verifiable date | highest ClinicalTrials.gov phase in a US site |
| EMA | EMA centralised marketing-authorisation date; else DrugCentral date; else two-agent consensus | EU national authorisation (Article 57) or agent-reported approval | highest phase in an EU/EEA site |
| NMPA | NMPA English news, or NMPA register year within 10 years of global first approval | NMPA register exists but the year is a re-registration/generic year | highest phase in a mainland-China site |
| PMDA | KEGG new-ingredient year; else DrugCentral date; else consensus | KEGG marketed-drug list or agent-reported approval | highest phase in a Japan site |

Other construction details:

- Molecules are unique by full InChIKey (6,183 unique structures, 6,183 unique names).
- Case/salt/ester duplicates were merged; conflicting duplicate structures were resolved against PubChem.
- “NA” means “not found in the public sources used here”, not proof of non-approval.
- Clinical phases are the **maximum** reported by any agent. A phase is a trial-stage proxy, not a regulatory submission.
- NMPA years are deliberately sparse (386 of 1,841 approved); many are “Approved” because the public snapshot cannot distinguish first approval from re-registration.
- EMA “Approved” may be a national EU/EEA authorisation rather than an EMA centralised approval; the distinction is recorded in `evidence_table_FINAL.csv`.

### Final coverage

| Agency | Year | Approved (no year) | Phase | NA |
|---|---:|---:|---:|---:|
| FDA | 2,035 | 293 | 2,547 | 1,297 |
| EMA | 770 | 1,596 | 1,704 | 2,110 |
| NMPA | 386 | 1,455 | 717 | 3,620 |
| PMDA | 688 | 919 | 534 | 4,041 |

---

## 10. What the final table does **not** claim

- It is not a regulatory determination. A “Phase 3” cell means a trial with a site in that jurisdiction, not an approval.
- NMPA “Approved” does not distinguish an originator approval from a generic/import registration.
- EMA “Approved” does not distinguish centralised from national authorisation.
- PMDA “Approved” can mean “marketed in Japan”, including products approved before 2000.
- FDA dates before 1982 can be re-filed application dates in Drugs@FDA.
- ClinicalTrials.gov under-covers China and Japan (ChiCTR/jRCT are separate registries), so NMPA/PMDA phases are understated.

---

## 11. Who performed best and worst

### Ranking

1. **Claude — best overall.**
2. **Gemini — second.**
3. **DeepSeek — worst overall**, despite having the best NMPA status coverage.

### Why Claude is best

- **Best structure quality:** 98.5% exact InChIKey match to ChEMBL, only 0.7% genuinely different; apart from somapacitan its scope is clean.
- **Best PMDA years:** 99.7% exact against KEGG, 1 case >2 years off; correctly dates imatinib 2001, bortezomib 2006, everolimus 2007, ambrisentan 2010.
- **Best FDA precision:** 98.5%, only 27 false positives.
- **Best documentation:** explicit source hierarchy, inclusion/exclusion rules, known limitations and a reproducible script.
- **Weakness:** NMPA is almost empty — 36 approvals, 4 years, 0.3% recall. It is not a complete four-agency table on its own. It also made a few year errors (opicapone EMA 2022 vs 2016; iobenguane FDA 1994 vs 2008) and has one structural error (eravacycline).

### Why Gemini is second

- **Broadest single-file coverage:** 4,997 rows, 2,121 FDA approvals, 1,247 EMA approvals, 644 NMPA approvals, 812 PMDA approvals.
- **Good FDA/EMA years:** 95.0% and 97.0% exact.
- **Safe NMPA strategy:** mostly “Approved (Year N/A)” rather than inventing years.
- **Weaknesses:** PMDA years are the worst of the three (80.3% exact; imatinib 2012 vs 2001, metformin 2010 vs 1961); no README or documented methodology; includes somapacitan; several “Approved (Year N/A)” cells are national approvals rather than EMA centralised approvals.

### Why DeepSeek is worst

- **NMPA years are systematically misleading:** 59.4% of its 1,563 NMPA years are 1998–2002, and 33.7% are exactly 2002 — a re-registration/numbering artefact, not first approval. Imatinib 2013 (actual 2002), sorafenib 2020 (actual 2006), abacavir 2015 (registered, not marketed).
- **PMDA years are often later approvals:** everolimus 2019 (actual 2007), bortezomib 2019 (actual 2006), ambrisentan 2021 (actual 2010). Coverage is also thin: only 393 PMDA approvals and 618 EMA years.
- **Scope violations:** includes interferon alfa-2b and ropeginterferon alfa-2b with small-molecule SMILES, cemdisiran (siRNA), elapegademase, pegunigalsidase alfa and somatrem. This directly contradicts the “exclude proteins, RNA” instruction.
- **Structure quality:** 7.4% of its ChEMBL-matched structures are genuinely different; e.g. dalfopristin is wrong.
- **Strengths:** it is the only file with broad NMPA status coverage (1,563 approvals, 97.4% recall) and its FDA years are close to the other two (94.6% exact). Those strengths are why it is not last on coverage — but they do not compensate for the year and scope problems.

---

## 12. Files

| File | Content |
|---|---|
| `small_molecule_drug_approvals_FINAL.csv` | Final reconciled table: `Name, SMILES, FDA, EMA, NMPA, PMDA` (6,183 rows) |
| `evidence_table_FINAL.csv` | Audit table: each agent’s original cell, the final value, and the basis/source for each agency |
| `verified_discrepancies_FINAL.csv` | 20 representative disagreements checked online, with verdicts and URLs |
| `CROSSCHECK_report.md` | This report |

The scripts used for the cross-check are in `work/` (`build_master.py`, `match_refs.py`, `build_final3.py`, `finalize_v3.py`, `merge_dupes.py`), and the final build can be reproduced from the raw sources in `data/raw/` and `data/nmpa/`.
