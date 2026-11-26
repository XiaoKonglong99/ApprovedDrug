# Gemini — round 2 (cross-check)

Gemini wrote its round-2 files (`build_final.py`, `small_molecule_drugs_final_verified.csv`, 6,062 molecules)
to its own scratch directory outside this project, so they are not included here.

Its review ranked the agents **Claude > Gemini > DeepSeek**, the same order as the other two reviews. Its main findings:
- Claude: most reliable overall, best PMDA first-approval years, sparse NMPA coverage, a few EMA years taken
  from CHMP opinions instead of EC decisions.
- Gemini (itself): best EMA years and format adherence, but Phase 4 → "Approved" false positives and parent-year
  propagation onto esters/prodrugs.
- DeepSeek: missed all trial-only molecules, included biologics/siRNA, EMA orphan-designation dates used as
  approval years; its NMPA dataset was the most valuable single contribution.
