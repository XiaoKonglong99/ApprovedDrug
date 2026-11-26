# Approved & clinical small-molecule drugs: a three-agent build, cross-check and merge

A table of **6,256 small-molecule drugs** with their **SMILES** and their status at six regulators:
**FDA** (US), **EMA** (EU centralised), **EU national** authorities, **NMPA** (China), **PMDA** (Japan) and **ANVISA** (Brazil).
It also records how the table was made. Three AI coding agents were given the same task, each built its own
table, each then cross-checked all three, and the results were merged, checked online and extended.

Main file: **[`[BEST]final_small_molecule_drugs_Claude.csv`](%5BBEST%5Dfinal_small_molecule_drugs_Claude.csv)**

```
Name,SMILES,FDA,EMA,EU_National,NMPA,PMDA,ANVISA
Imatinib,Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1,2001,2001,2016,2013,2001,2001
Rivoceranib,N#CC1(c2ccc(NC(=O)c3cccnc3NCc3ccncc3)cc2)CCCC1,Stage 3,Stage 3,Stage 3,2014,Stage 3,Stage 3
Dexketoprofen,C[C@H](C(=O)O)c1cccc(C(=O)c2ccccc2)c1,NA,Stage 4,1996,Approved,NA,NA
```

---

## 中文摘要（不严谨测评：让 AI Agents 爬取已批准药物的分子结构）

我们把**完全相同的指令**交给三个 AI Agent：Claude Opus 5.5 High、Gemini Neon Medium 和 DeepSeek V4.1 Flash Max。每个 Agent 在各自的 harness 里独立编写脚本，爬取美国 FDA、欧盟 EMA、中国 NMPA、日本 PMDA 已批准及临床在研的小分子药物和 SMILES 结构。随后让三者**交叉核验**彼此的结果并联网复查。

- **第一轮（独立爬取）**：三份表差异极大。Claude 5,607 个分子，Gemini 4,997 个，DeepSeek 2,633 个（DeepSeek 漏掉了全部"仅在研"分子）。
- **第二轮（交叉核验）**：三个 Agent **独立给出了完全一致的排名：Claude > Gemini > DeepSeek**。
- **第三轮（Claude 补充修正）**：
  - 补上 Carrimycin（可利霉素，混合物，列出全部三种活性组分）和 2014–2026 年间漏掉的大批中国 1 类新药。漏掉的原因有两个：药监局更名（SFDA→CFDA→NMPA），以及中英文/INN 命名不一致。
  - 新增巴西 ANVISA 和欧盟成员国国家审批（EU_National）两列。
  - 最后提取大环分子（MRS ≥ 12，195 个），并计算全部 MolDs 和 Onion 层级。

### 总结

1. **排名一致，但没有赢家能单独交付。** 三个 Agent 在第二轮独立得出同一排名，这本身说明差距是真实的。但每一家都有别人替代不了的数据：Claude 的结构和日本首批年份最准；Gemini 的 EMA 中央审批年份最准、格式最规范；DeepSeek 是唯一挖出 NMPA 历史批准记录的，也最会找老药和最新的 FDA 批准。任何单一 Agent 的表都不能直接作为科研级数据集。
2. **错误主要不是"爬不到"，而是"年份的含义不同"。** 同一个药，不同来源给的可能是首次批准、再注册（中国 2002 年批量换发批准文号）、新剂型或补充适应症、CHMP 意见日期，或数据库里最早的一条记录。三方分歧大部分来自这里，其次来自命名、盐型、前药和混合物的归并。
3. **数据源的结构性盲点，比 Agent 能力更决定结果。** 中国数据缺失源于药监局更名和中文名映射。欧盟有大量只走成员国审批、从不出现在 EMA 中央数据库里的老药。巴西 ANVISA 1999 年成立前的记录不完整。不专门处理这些盲点，换哪个 Agent 都会漏。
4. **"多 Agent 交叉 + 联网抽检 + 规则清洗"有效，但抽检有限。** 联网核实了约 280 个单元格/条目，其余约 3,700 个冲突单元格按抽检得出的规则裁决（例如：日本年份以 Claude 为准；中国 1998–2003 年的年份视为再注册年份；单独由 Gemini 声称的批准不采信）。
5. **成本差异很大。** Gemini 最快（23.5 分钟）；DeepSeek 最慢（约 39.5 分钟），且因反复重读上下文，总 token 约 6,000 万；Claude 约 35 分钟。
6. **偏差声明：** 第三轮的修正和这份 README 由 Claude 完成。排名以联网核实的打分为依据（见下文），但读者应考虑"自己评自己"的可能偏差。

---

## Contents

```
[BEST]final_small_molecule_drugs_Claude.csv   final table (6 regulators)
results/
  final_small_molecule_drugs_detailed_Claude.csv   + InChIKey, every source's raw value, evidence for each cell
  crosscheck_conflicts_Claude.csv                  every cell where the three agents disagree (3,681)
round1_agents/            round 1: each agent's own script, README and table, as delivered
  claude/  gemini/  deepseek/
round2_crosscheck/        round 2: each agent's cross-check of all three tables
  claude/    analysis scripts, PubChem SMILES check, online-verification log with per-agent scores
  deepseek/  DeepSeek's report, its final table and evidence table
  gemini/    (Gemini's round-2 files were written outside this folder; see its README)
pipeline/                 round 3: final merge + the extra sources
  crosscheck_merge_Claude.py   builds the final table
  curation/      web-verified overrides, additions and name aliases (hand-checked)
  nmpa_fix/      China: Chinese-name rescue + verified list of China-first innovative drugs
  anvisa/        Brazil: ANVISA open data + Brazil trial sites
  eu_national/   EU national authorities: France ANSM, Spain AEMPS, EMA Article 57
  common/        multilingual name mapper, Wikidata label lookup
App_CT_MC/                macrocycles: is_MC flag, MolDs + Onion layers (MacrocycleWizard, not included)
data/                     shared downloads and caches (large raw files are not published)
run_all.sh                rebuilds everything after round 1
```

---

## Round 1: each agent builds its own table

**Prompt (identical for all three agents):**

> Please write a script and help me collect small-molecule drug (exlculde protein, RNA...) and their SMILES structure approved or in clinical trail by FDA (US), EMA (EU), NMPA (CHina), PMDA (Japan). Give me both the code and the csv file like:
> `Name. SMILES FDA EMA NMPA PMDA`
> If a molecule was approved by more than one agency, make it one col and leave the unapproved stage, approve for approve year. For example, if Mole 1 approved by FDA in 2019, in stage 1 in EMA and stage 2 in NMPA but no record in PMDA, like:
> `Mole1. CCCCCCXX. 2019. Stage 1 stage 2 NA`

| | Claude Opus 5.5 High | Gemini Neon Medium | DeepSeek V4.1 Flash Max |
|---|---|---|---|
| Wall-clock time | ~35 min (≈15 min of it the script running) | 23.5 min | ~39.5 min (25.3 min model, 14.4 min tools) |
| Tokens | context ≈165k at the end | 194k new (2.43 M incl. cache reads) | 358k new (≈60.3 M incl. cache reads; 245 steps) |
| Molecules | **5,607** | 4,997 | 2,633 (no trial-only molecules) |
| Backbone | ChEMBL small molecules + ClinicalTrials.gov | ChEMBL + ClinicalTrials.gov | DrugCentral + PubChem |
| FDA | Drugs@FDA (openFDA) | Drugs@FDA + DrugCentral | DrugCentral + FDA novel approvals |
| EMA | KEGG EU new substances, EMA lists, Article 57 (national, no year) | EMA JSON + DrugCentral | EMA medicines workbook |
| NMPA | NMPA English news (≈2021–2025) → only 36 approvals | Phase-4 trial in China → "Approved" | NMPA database snapshot to Aug 2021 → 1,563 years |
| PMDA | KEGG Japan approvals (2000–) + marketed list | PMDA PDF list + DrugCentral | DrugCentral |
| Clinical stage | trial sites per region | trial sites per region (Stage 1–3) | trial sites, approved drugs only |

Scripts and tables: [`round1_agents/`](round1_agents). In Claude's and DeepSeek's scripts only the path constants were changed so they find the shared `data/` folder; nothing else was edited.

## Round 2: each agent cross-checks all three tables

**Prompt:** "I gave the task and prompt to the three AI Agents I have: *(round-1 prompt)*. And here are the results provided by three AI agents. Please cross check the results provided by them, and double check if there is anything difference. If there is any, double check it online and give me the final version. Finally, give me the comments on who performed the best who performed the worst."

**All three agents independently reached the same ranking: 1 Claude, 2 Gemini, 3 DeepSeek.**

### Overlap

- Claude's check links rows by name, by salt-stripped parent InChIKey and by synonyms. It finds **6,200** distinct molecules, of which **2,008 (32%)** are in all three tables.
- DeepSeek's check uses normalised names only. It finds **6,419** names, of which **1,828 (28.5%)** are in all three.
- Claude's and Gemini's SMILES agree for 100% of shared molecules, because both come from ChEMBL. DeepSeek's differ on 82 molecules:
  - 16 have a different skeleton. Checked against PubChem, DeepSeek was wrong on 6 (macimorelin, uridine triacetate, vaniprevir, dalfopristin, ensifentrine, enarodustat) and Claude/Gemini on 3 (sodium gluconate, voxilaprevir, ulipristal vs ulipristal acetate).
  - The other 66 are stereochemistry or tautomer variants.

### Online verification of conflicting cells (Claude's check)

85 cells where the agents disagreed were checked online. The sample was taken from conflicts, so these scores are lower than overall accuracy. Details: [`round2_crosscheck/claude/verification_log_Claude.csv`](round2_crosscheck/claude/verification_log_Claude.csv).

| Agency | Cells checked | Claude | Gemini | DeepSeek |
|---|---|---|---|---|
| FDA | 31 | 55% | 58% | 47% |
| EMA | 16 | 62% | 56% | 44% |
| PMDA | 17 | **91%** | 21% | 18% |
| NMPA | 21 | **50%** | 12% | 24% |
| **All** | 85 | **62%** | 39% | 35% |

### What each review found

**Claude (best).**
- Strengths: widest coverage, cleanest structures, best Japanese first-approval years (DeepSeek's review: 99.7% agreement with KEGG), documented limitations.
- Weaknesses:
  - Almost no China data (36 approvals).
  - Misses FDA approvals filed under salt names or withdrawn long ago (ixazomib, lorcaserin, imlunestrant, amodiaquine).
  - Some EMA years are CHMP-opinion dates rather than EC decisions (Gemini's review).
  - 4 unrequested extra columns.

**Gemini (second).**
- Strengths: fastest, followed the requested format exactly, best EMA centralised years.
- Two systematic errors:
  - It turns any Phase 4 trial into "Approved", giving false approvals (huperzine A, idebenone and melatonin at FDA; conivaptan in the EU; dexlansoprazole in Japan).
  - It copies a parent drug's approval year onto its esters and prodrugs (tenofovir alafenamide → 2001).
- Its Japanese years are usually those of a later product (imatinib 2012 instead of 2001). Almost no China data.

**DeepSeek (third, but irreplaceable).**
- Weaknesses:
  - It left out every trial-only molecule.
  - It let in biologics and siRNA (interferon, cemdisiran, somapacitan) and had the most SMILES errors.
  - Its FDA years for old drugs are often late, and two EMA "years" are orphan-designation dates (maribavir 2013, eliglustat 2001).
- Strengths: it was the only agent to bring in a large historical **NMPA** dataset, and it found the most old and newest FDA approvals.
- Catch: 979 of its 1,563 China years fall in 1998–2003 (527 in 2002 alone). That is China's nationwide re-registration, not first approval.

Gemini's own reconciled table had 6,062 molecules and DeepSeek's 6,183; both are summarised in the agents' reports. Claude's round-2 table, 6,088 molecules, is the base of this repository.

## Round 3: corrections and additions (Claude)

### China (NMPA): renamed agency, Chinese names, mixtures

Carrimycin (可利霉素, approved 2019-06-24) and many China-first drugs were missing or shown only as trials. The causes:
- The NMPA was named **SFDA** before 2013 and **CFDA** in 2013–2018. Claude's round-1 source, the NMPA English news site, only covers about 2021 onward.
- DeepSeek's NMPA records use Chinese English names (chidamide, apatinib, almonertinib) that never matched international INNs (tucidinostat, rivoceranib, aumolertinib). It also dropped every record whose English name had a typo or dosage-form words, about 2,700 records.
- **Carrimycin is a mixture**, so neither PubChem nor ChEMBL has a structure under that name.

Fixes:
1. **Carrimycin is listed as its three active components**: isovalerylspiramycin I and III (PubChem CID 101485376 / 101485378) and II, built as spiramycin II 4''-O-isovalerate; its formula and core stereochemistry were checked against spiramycin II. Values: NMPA 2019; FDA Stage 3 (NCT04672564, US site).
2. The unresolved NMPA records were cleaned (dosage-form and salt words stripped from the Chinese names) and matched to structures through **Wikidata Chinese labels**. This rescued 1,214 records and 399 molecules.
3. A **web-verified list of 145 China-first innovative small molecules (2001–2026)** was compiled from the NMPA 创新药品名录, CDE annual review reports and company announcements. Structures came from ChEMBL or PubChem; 139 of 145 resolved. See [`pipeline/nmpa_fix/china_innovative_verified_Claude.csv`](pipeline/nmpa_fix/china_innovative_verified_Claude.csv).
4. A hand-checked alias list joins INN and local names, e.g. chidamide = tucidinostat, anlotinib = catequentinib.
5. DeepSeek's NMPA years from 1998–2003 are kept as **Approved** without a year (re-registration window; wrong in 5 of 7 checked cases).

### Brazil (ANVISA)

ANVISA is a full ICH member and the reference regulator in Latin America. It is a credible source for neglected-tropical-disease drugs and for some domestically developed molecules.

- **Approvals**: ANVISA open data, `DADOS_ABERTOS_MEDICAMENTOS.csv`, with 229,909 registered medicines. Kept: new, generic, similar and specific products, radiopharmaceuticals, medical gases, and "biological" (Brazil files synthetic peptides there). Excluded: herbal, homeopathic and notified low-risk products.
- **Name mapping**: Portuguese DCB names (e.g. *cloridrato de metformina*, *dipirona*) are mapped to ChEMBL through the ChEMBL name index (with Portuguese salt and spelling rules), Wikidata pt/es labels and a short DCB→INN dictionary. 87% of product records are matched.
- **Year**: the earliest registration of any product with that ingredient. When a product's grant date is a later renewal, the filing year in the 15/17-digit process number is used instead.
- **Rule**: a pre-2000 year that is ≥10 years later than the drug's first approval at FDA, EMA, an EU national authority or PMDA becomes *Approved*. Records from before ANVISA's creation in 1999 are incomplete; for example metamizole shows 1976 but has been sold in Brazil since 1922.
- **Trial stage**: ClinicalTrials.gov studies with a site in Brazil.
- **Checked online**: 31 cells ([`pipeline/anvisa/anvisa_verification_Claude.csv`](pipeline/anvisa/anvisa_verification_Claude.csv)). Years from about 2000 onward were all confirmed (tafenoquine 2019, sofosbuvir 2015, osimertinib 2016, selpercatinib 2024, nifurtimox 2026). Corrections: ceftaroline 2014, metamizole and oxamniquine → Approved, meglumine antimoniate added.

### EU national authorities (EU_National)

EMA only runs the **centralised procedure** (since 1995; mandatory mainly for biologics, oncology and orphan drugs). Hundreds of older small molecules, regional antibiotics and niche drugs were authorised through **national, mutual-recognition or decentralised procedures**, for example by BfArM, ANSM, AIFA or MPA. They appear only in national registers.

The old EMA column is therefore split:
- **EMA**: EU centralised authorisation year only.
- **EU_National**: first national, MRP or DCP authorisation found, from:
  - **France, ANSM** public medicines database: authorisation date and procedure type for every product. National, decentralised and mutual-recognition products only.
  - **Spain, AEMPS** CIMA register (REST API): authorisation date. Products authorised through EMA are excluded.
  - **EMA Article 57** database: every medicine authorised in every EU/EEA country, no dates. A match gives *Approved*.

Germany (BfArM), Italy (AIFA) and Sweden (MPA) do not publish dated bulk lists, so they are covered only through Article 57, without years.

### Smaller fixes

- Removed somapacitan (a protein).
- Some structures matched only through Wikidata labels were wrong or out of scope, e.g. ANVISA's *urofolitropina* (urinary FSH, a protein) had picked up an unrelated cyclic nonapeptide, and pegaptanib (an RNA aptamer) only its PEG linker. These Wikidata items are blocked in [`pipeline/curation/wikidata_blocklist_Claude.csv`](pipeline/curation/wikidata_blocklist_Claude.csv), which removed urofollitropin, pegaptanib, peroxide (wildcard SMILES), hypromellose and carboxymethyl cellulose (polymers) and restored invert sugar as its own row.
- Restored platinum in lobaplatin and nedaplatin, whose agent SMILES had lost the metal. Added cisplatin and miriplatin, which no agent had because ChEMBL lists them without a small-molecule structure.
- Molecules that came only from the new sources are filled from DrugCentral (the source DeepSeek used) for FDA/EMA/PMDA.
- Every manual change, with its source, is in [`pipeline/curation/`](pipeline/curation): 66 cell overrides, 145 additions, 84 name aliases.

## The final table

**Columns:** `Name, SMILES, FDA, EMA, EU_National, NMPA, PMDA, ANVISA`

| Cell | Meaning |
|---|---|
| `2019` | approved by that regulator; year of first approval |
| `Approved` | approved or marketed there, no reliable first-approval year |
| `Stage 1` … `Stage 4` | not known to be approved there; highest phase of a trial with a site in that region (`Stage 1/2`, `Stage 2/3` = combined phases, `Stage 0/1` = early phase 1, `Stage 4` = post-marketing study, i.e. very likely marketed) |
| `NA` | no record found in the sources used (not proof of absence) |

| Column | Year | Approved (no year) | Trial stage only | NA |
|---|---|---|---|---|
| FDA | 2,024 | 132 | 2,650 | 1,450 |
| EMA (centralised) | 730 | 0 | 2,444 | 3,082 |
| EU_National | 1,239 | 844 | 1,811 | 2,362 |
| NMPA | 717 | 998 | 715 | 3,826 |
| PMDA | 542 | 933 | 567 | 4,214 |
| ANVISA | 1,176 | 366 | 482 | 4,232 |

- **6,256** molecules: 3,575 approved by at least one regulator, 2,681 in clinical trials only, 231 approved by all six.
- **Scope**: small molecules up to 2,000 Da. Salts and hydrates are folded into the parent, and the SMILES is the parent structure. Proteins, antibodies, oligonucleotides, cells and gene therapies are excluded; ChEMBL classes cyclic peptides such as ciclosporin as proteins, so most of those are excluded too. Esters and prodrugs that are separate drugs stay separate rows (e.g. tenofovir alafenamide).
- **Mixtures** are listed one row per active component when the components are defined (Carrimycin).

### How each cell is decided (`pipeline/crosscheck_merge_Claude.py`)

1. **Linking**:
   - Rows from all sources are joined by normalised name, by the hand-checked alias list, and by parent InChIKey.
   - InChIKey joins are used when the key is unique in each file, or when the parent has ≥8 heavy atoms. This folds "sodium cromoglicate" into cromolyn but keeps lithium carbonate and calcium carbonate apart.
2. **SMILES**: ChEMBL-based (Claude/Gemini) first, then verified additions, then other sources; PubChem-verified corrections applied.
3. **Approval and year, per regulator**:
   - Majority year; on ties: FDA → Claude > Gemini > DeepSeek, EMA → Claude > DeepSeek > Gemini.
   - **PMDA**: Claude's KEGG year wins.
   - **NMPA**: DeepSeek's 1998–2003 years → Approved.
   - A lone Gemini "Approved (Year N/A)" that no other agent supports is ignored, because it was wrong in every checked case.
4. **Not approved**: the highest trial phase any source reports for that region.
5. **Web-verified overrides and additions** are applied last.

### Macrocycles (`App_CT_MC/`)

Uses the MacrocycleWizard package (not included here; it has its own repository).

| File | Content |
|---|---|
| `[BEST]ApprovedDrug_CT_2026Nov.csv` | copy of the final table |
| `[BEST]ApprovedDrug_CT_2026Nov_isMC.csv` | + `MRS` (macrocycle ring size: bonds in the largest ring after stripping substituents) and `is_MC` (MRS ≥ 12) |
| `ApprovedDrug_CT_2026Nov_MC_only.csv` | the **195** macrocycles (MRS 12–36; e.g. the three Carrimycin components and utidelone, MRS 16; sirolimus, MRS 29) |
| `ApprovedDrug_CT_2026Nov_MC_MolDs_Onion.csv` | all **151 MolDs** (groups A–K) and the **Macrocycle Onion** layer SMILES (ring, peripheral, substituent, …) for every macrocycle; 0 failures |
| `ApprovedDrug_CT_2026Nov_MC_routes.csv` | **route of administration** per macrocycle with the evidence from each source; the class (`Route_Class`: Oral 80 / Injection 83 / Oral + Injection 24 / Non-oral (other) 8), `Routes` and `Route_Source` are also added to the two macrocycle files above |

`add_route_info.py` decides the route from marketed products for approved drugs (Drugs@FDA, EMA Article 57, France ANSM, Chinese product dosage forms), from ClinicalTrials.gov intervention descriptions for trial-only compounds, and from ChEMBL route flags as a last resort; 38 doubtful cases were checked online (`route_overrides_Claude.csv`).

`calc_mc_molds_onion.py` calls MacrocycleWizard's `CalculateAll` without modifying the package. It adds one fallback: when cutting off a substituent splits an aromatic system that RDKit cannot kekulize (e.g. rifaximin), it kekulizes the whole molecule first.

## Reproducing

```bash
pip install -r requirements.txt
./run_all.sh                      # China rescue -> China list -> ANVISA -> EU national -> merge -> macrocycles -> routes
```

- `run_all.sh` uses the caches in `data/` and `pipeline/*/raw/` when present. Otherwise it downloads ANVISA open data, ANSM BDPM, AEMPS CIMA, the EMA Article 57 workbook, and ChEMBL / ClinicalTrials.gov through Claude's round-1 loader; the last one takes about 15 minutes.
- The macrocycle step needs an environment that can import MacrocycleWizard. Set `WIZARD_PY=/path/to/python`.
- Re-running the pipeline after the folder reorganisation reproduced the final table exactly.
- Round-1 tables are inputs and are not rebuilt. Each agent's script is in `round1_agents/` (Claude: `python3 collect_drugs.py`; DeepSeek: `python3 collect_drugs_Deepseek.py`, needs `pg8000`; Gemini: `python3 collect_small_molecule_drugs_Gemini.py --output …`).
- Online lookups (Wikidata, PubChem) are cached in the repository (`wikidata_*_hits.csv`), so results are stable.

## Limitations

- **Coverage of the online checks**: about 280 items were checked online (85 conflicting cells, 145 China-first drugs, 31 ANVISA cells, 16 SMILES skeletons). The other conflicts (3,681 cells) are settled by rules derived from those checks.
- **"First approval" years are only as good as each source**:
  - FDA dates before 1982 in Drugs@FDA are sometimes re-filing dates.
  - EU_National years are the earliest product *still listed* in France or Spain, so older withdrawn products are missing.
  - ANVISA records before 1999 are incomplete.
  - NMPA years are reliable for the verified innovative-drug list and for database records after 2003; older China approvals (998) have no year because the public snapshot only holds re-registration dates.
- **EU_National** has years only from France and Spain. Germany, Italy, Sweden and other members contribute presence only, through Article 57.
- **Metal-containing salts**: for salts such as lithium carbonate or silver nitrate, the agents' salt stripping kept the counter-ion and dropped the active metal. This was fixed for the platinum drugs, cisplatin and meglumine antimoniate only.
- **Trial coverage**: ClinicalTrials.gov under-covers China, Japan and Brazil; trials registered only in ChiCTR, jRCT or ReBEC are missed.
- **Name mapping** (Chinese, Portuguese, French, Spanish → INN) matched 81–90% of product records. The rest are mostly combinations, herbals, vitamins, excipients or misspellings.
- **Peptides**: inclusion is inconsistent between ChEMBL-based and DrugCentral-based rows (≤2,000 Da peptides from DrugCentral are kept; ChEMBL "protein" peptides are not).
- **Bias**: the round-3 work and this README were done by Claude, one of the three agents. Rankings rest on the verification scores above.

## Data sources and terms

| Source | Use | Terms |
|---|---|---|
| ChEMBL (EMBL-EBI) | structures, names, trial links | CC BY-SA 3.0 |
| ClinicalTrials.gov API | trial phases and sites | public domain (US NLM) |
| Drugs@FDA / openFDA | FDA approvals | public domain |
| DrugCentral | approvals (DeepSeek, fill-in) | CC BY-SA 4.0 |
| EMA medicines data, Article 57 database | centralised and national EU authorisations | EMA legal notice (reuse with attribution) |
| KEGG BRITE br08301 / br08318 / br08329 | Japan and EU approvals (Claude) | free for academic use |
| NMPA database snapshot (github xprincehf/nmpa-data), NMPA 创新药品名录, CDE annual reports | China | public notices |
| ANVISA open data (dados.anvisa.gov.br) | Brazil registrations | ANVISA open-data terms |
| ANSM Base de données publique des médicaments | France national authorisations | BDPM terms of reuse |
| AEMPS CIMA | Spain authorisations | AEMPS reuse conditions |
| PubChem, Wikidata | structures, multilingual names | public domain / CC0 |

Raw downloads are not redistributed; the scripts fetch them. Derived tables contain ChEMBL- and DrugCentral-derived structures, so ShareAlike terms may apply to redistributed tables.
