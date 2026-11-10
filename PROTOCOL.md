# Cyclophilin A: real-data baseline pilot

Frozen before calculating residue-level scores, 2026-11-10. The experimental paper and qualitative rescue results have already been read; this is retrospective, not a blind or prospective study.

## Scope

The longer-term question is whether an all-atom generator recovers the coordinated structural response to suppressor mutations. This first pilot is a feasibility/measurement study, not a test of that full claim. We will ask whether simple local baselines or a small frozen open-weight sequence model identify which residues' NMR environments change when I97V is added to S99T/C115S.

No novel biological mechanism, kinetic prediction, causal coupling, activity prediction, calibrated populations, clinical utility, or broad model generalization will be claimed. No model training or paid compute. BioEmu2 requires Linux/NVIDIA and is not supported on this M2 Air; do not download a checkpoint that cannot be run here.

## Data and target

- BMRB 27217 (S99T/C115S), 27218 (S99T/C115S/I97V), raw NMR-STAR files. Match atom identity, residue numbering, and temperature using deposit metadata, never file/list order alone.
- PDB 6BTA is the source/double-mutant structure for distance baselines. 5WC7, 3K0N, 3K0O are reference/audit only. Never use target-mutant geometry to calculate predictive scores.
- Match backbone H and N per residue; target CSP = sqrt((delta H)^2 + (delta N / 5)^2). This is an operational ranking target, not a motion amplitude or experimental uncertainty.
- Analyze 283 K (10 C) and 298 K (25 C) separately. The second temperature is a sensitivity check on the same mutation, not an independent protein cohort.
- Primary evaluation excludes changed site 97 and two residues either side (95-99). Report all matched residues as sensitivity, plus excluding five either side (92-102).

## Frozen baselines and metrics

1. Sequence locality: negative sequence distance to site 97.
2. Geometry locality: negative C-alpha distance to site 97 in 6BTA; use highest occupancy CA when alternatives exist. Missing sites are excluded consistently for all methods.
3. If feasible, ESM2 t6 8M UR50D, frozen last-layer contextual embedding Euclidean change between double and triple sequences, excluding special tokens. It is a sequence representation baseline, not an ensemble/dynamics model. No layer/model choice based on target correlations.

Primary metric: Spearman rank correlation with CSP. Secondary: precision among the top ceil(n/4) ranked residues, using top ceil(n/4) measured CSP residues as targets (ties deterministically broken by residue index). Report both values and chance top-k precision k/n. Residues are correlated; no IID significance claims or population inference. Record per-residue scores and all exclusions. Zero-change baseline has constant score and no meaningful Spearman ranking; do not fabricate one.

## Measurement and structure checks

Assert deposited sequences differ only at canonical site 97 (I to V). Validate NMR sequence against each PDB and author numbering. Inspect chemical-shift reference/TROSY-offset metadata; shifts with missing/ambiguous assignments are excluded and counted. Compare H and N separately to expose scaling sensitivity. Audit crystal temperatures, alternate conformers and site 97/99/113/115 availability. Alternate-location labels are not proof of globally correlated molecular states; crystal occupancy is not a solution-state population.

## Stop conditions and next step

If mapping/temperature/reference checks fail, stop score calculation and document the failure. If sequence baseline is poor, retain the negative result. Do not grow model size, tune a layer or adjust exclusions to manufacture success. A possible later model run must test coordinated all-atom geometry with separate experimental support; this pilot alone cannot establish that endpoint.
