"""Render measured outputs without fitting or choosing a best model."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent


def main():
    s=json.loads((ROOT/'results/summary.json').read_text())
    rows=list(csv.DictReader((ROOT/'results/residue_scores.csv').open()))
    names={'sequence_locality':'Sequence distance','geometry_locality':'3D distance',
           'esm2_embedding_change':'ESM2 8M embedding change'}
    fig,axes=plt.subplots(2,2,figsize=(10,6),layout='constrained')
    table=[]
    for col,temp in enumerate((283,298)):
        metrics=s['metrics'][f'{temp}K_exclude_2']['csp']
        keys=list(metrics)
        axes[0,col].bar([names[x] for x in keys],[metrics[x]['spearman'] for x in keys],color=['#78909c','#176b87','#dc7e23'][:len(keys)])
        axes[0,col].set(title=f'{temp-273}°C — {metrics[keys[0]]["n"]} residues',ylabel='Spearman correlation',ylim=(-.2,1))
        axes[0,col].tick_params(axis='x',labelsize=8,rotation=12)
        r=[x for x in rows if int(float(x['temperature_K']))==temp and abs(int(x['residue'])-97)>2]
        axes[1,col].scatter([-float(x['geometry_locality']) for x in r],[float(x['csp']) for x in r],s=12,alpha=.7,color='#176b87')
        axes[1,col].set(xlabel='Source-structure Cα distance from mutation (Å)',ylabel='Measured CSP (ppm)')
        for k,m in metrics.items():
            table.append(f'| {temp-273}°C | {names[k]} | {m["n"]} | {m["spearman"]:.3f} | {m["top_quartile_precision"]:.3f} |')
    fig.suptitle('CypA I97V addition: retrospective NMR perturbation localization\nOne mutant pair; temperatures are sensitivity checks, not independent proteins',fontsize=12)
    fig.savefig(ROOT/'results/baselines.png',dpi=180)
    plt.close(fig)
    report='''# Cyclophilin A pilot: measured baseline results

## What was actually tested

Can sequence distance, source-structure distance, or a small frozen sequence representation rank the residues whose measured backbone NMR chemical shifts change when I97V is added to S99T/C115S? This is a retrospective localization test on one mutant pair. It does not test whether a model predicts motion, activity, coordinated side-chain transitions, or equilibrium probabilities.

## Results

Site 97 and the two residues on each side are excluded in the primary analysis. NMR assignments and experimental conditions are matched before calculating shifts. Larger scores predict larger measured perturbations.

| Temperature | Method | Residues | Spearman | Top-quarter precision |
|---|---|---:|---:|---:|
'''+ '\n'.join(table)+'''

Chance top-quarter precision is 0.253 at 10°C and 0.250 at 25°C. Residues are correlated and these are two temperatures of the same mutant pair, so neither is an independent validation cohort. All-site, ±5-residue exclusion, H-only and N-only sensitivities are in `summary.json`, with every residue in `residue_scores.csv`.

![Baseline results](baselines.png)

## Data integrity and limits

- 158 sites have paired H/N measurements plus source CA coordinates at 283 K; 156 at 298 K. After the primary local exclusion, 154 and 152 remain. Missing residues are listed explicitly; no imputation.
- Deposits 27217 and 27218 differ only by I97V, both at pH 6.5, and both document TROSY offsets. Their detail strings differ only in whitespace; initial literal metadata equality stopped the run, and whitespace normalization resolved it before scoring. No numerical shift values were changed.
- This does not prove the absence of a between-sample referencing offset. Error values include zeros; they are not a validated measurement noise floor. Small CSPs should not be interpreted as significant biological changes.
- Source structure 6BTA and target 5WC7 were measured at 273 K; WT 3K0N and S99T 3K0O at 288 K. Do not combine crystal occupancies with solution NMR as matched-temperature probabilities.
- The distance predictor uses 6BTA only. Target-mutant coordinates are reserved for identity/data auditing, not predictive scores. The pretrained model may have seen these sequences/structures during training; this is not a held-out generalization claim.
- NMR chemical shifts report local magnetic environments. Their changes are not direct measurements of flexibility, exchange rates, or catalytic effects. The published NMR exchange states also must not be equated automatically with individual crystallographic rotamers.
- BMRB assignments are machine-readable. Relaxation profiles are publicly plotted with fitted parameters, but raw instrument data were not located. No relaxation model was refitted.
- The ESM2 score, if present, is a frozen last-layer embedding difference. It has no learned calibration to ppm and generates no atomic ensemble. BioEmu2 was not run: the released software requires Linux and NVIDIA hardware, which this Mac lacks.

## What this enables next

Treat 3D proximity as a baseline to beat before claiming useful nonlocal mutation predictions. A separate all-atom study would have to recover compatible joint packing arrangements and test evidence beyond backbone RMSD. Four variants and these NMR tables cannot establish a new general benchmark, causal allostery, or publishable biological discovery.

## Provenance and related work

- Otten et al. (2018), experimental rescue series: https://doi.org/10.1038/s41467-018-03562-9
- BMRB 27217 and 27218; PDB 3K0N, 3K0O, 6BTA and 5WC7. Exact downloaded bytes and SHA-256 hashes are in `summary.json`.
- Wapeesittipan et al. (2019), simulations and alternative timescale interpretation: https://doi.org/10.1038/s42004-019-0136-1
- Cagiada et al. (2024), conformational free-energy prediction including CypA: https://doi.org/10.1002/pro.5233

The mutation-rescue mechanism is prior work. We have not established algorithmic or biological novelty. No exact public BioEmu2 test of this rescue series was located; absence of a named match does not exclude internal evaluation or training exposure.
'''
    (ROOT/'results/REPORT.md').write_text(report)


if __name__=='__main__':main()
