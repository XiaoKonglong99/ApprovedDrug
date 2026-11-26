#!/usr/bin/env bash
# Rebuild everything after round 1 (the three agents' tables in round1_agents/ are inputs, not rebuilt here).
# Uses cached downloads in data/ and pipeline/*/raw/ when present; otherwise downloads them.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}                  # pandas + rdkit + requests + openpyxl
WIZARD_PY=${WIZARD_PY:-python3}    # an environment that can import MacrocycleWizard (needs IPython)
export MACROCYCLEWIZARD_DIR=${MACROCYCLEWIZARD_DIR:-$PWD/App_CT_MC/MacrocycleWizard_CP}   # MacrocycleWizard checkout

echo "== China: rescue NMPA records through their Chinese names"
$PY pipeline/nmpa_fix/clean_zh.py
[ -f pipeline/nmpa_fix/wikidata_zh_hits.csv ] || \
  $PY pipeline/common/wikidata_labels.py pipeline/nmpa_fix/zh_names.txt zh,zh-hans,zh-cn pipeline/nmpa_fix/wikidata_zh_hits.csv
$PY pipeline/nmpa_fix/build_rescue.py
echo "== China: structures for the web-verified innovative-drug list"
$PY pipeline/nmpa_fix/resolve_china_list.py
echo "== Brazil: ANVISA registrations + Brazil trials"
$PY pipeline/anvisa/build_anvisa.py
echo "== EU national: ANSM (France), AEMPS (Spain), EMA Article 57"
$PY pipeline/eu_national/build_eu_national.py
echo "== Merge -> [BEST]final_small_molecule_drugs_Claude.csv"
$PY pipeline/crosscheck_merge_Claude.py
echo "== Macrocycles (MRS >= 12) + MolDs / Onion layers"
cp "[BEST]final_small_molecule_drugs_Claude.csv" "App_CT_MC/[BEST]ApprovedDrug_CT_2026Nov.csv"
$WIZARD_PY App_CT_MC/flag_macrocycles.py
$WIZARD_PY App_CT_MC/calc_mc_molds_onion.py
echo "== Macrocycles: oral vs non-oral route"
$PY App_CT_MC/add_route_info.py
