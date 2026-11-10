import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
s=json.loads((ROOT/'results.json').read_text())
rows=[]
for temp in (283,298):
    for name,m in s['metrics'][f'{temp}K_exclude_2']['csp'].items():
        rows.append(f'| {temp-273}°C | {name} | {m["n"]} | {m["spearman"]:.4f} | {m["top_quartile_precision"]:.4f} |')
carbon=[]
for name,m in s['carbon_298K_secondary'].items():
    carbon.append(f'| {name} | {m["n"]} | {m["spearman"]:.4f} | {m["top_quartile_precision"]:.4f} |')
interval=[]
for temp,d in s['paired_block_bootstrap'].items():
    lo,hi=d['descriptive_percentile_95']
    interval.append(f'- {int(temp)-273}°C: Δrho {d["spearman_difference"]:+.4f}, descriptive interval [{lo:+.4f}, {hi:+.4f}].')
text='''# Small geometry follow-up: a modest point gain, not a robust win

The primary, fixed method ranks each residue by the distance from the mutation's side-chain heavy atoms to its backbone nitrogen. It is a refinement of a distance baseline, not a trained model or a new biological mechanism. Parameters, primary candidate, exclusions and all secondary scores were recorded in `protocol.json` and committed before these follow-up calculations. Prior outcomes were already observed; this is retrospective and is not independent validation.

## H/N result

The primary candidate improves Spearman at both temperatures (0.359→0.369 and 0.523→0.540). Top-quarter precision ties the baseline at 10°C and improves by one residue at 25°C (25/38→26/38). Gains are small. This does not establish a reliable advantage.

| Temperature | Score | n | Spearman | Top-quarter precision |
|---|---|---:|---:|---:|
'''+ '\n'.join(rows)+'''

The graph diffusion score is secondary, not a newly selected primary winner. It also gives small point gains in rank correlation; its 25°C precision does not improve. Minimum all-heavy-atom distance does not improve rank correlation at either temperature. All candidates and all-site, ±5-exclusion, H-only and N-only sensitivities are in `results.json`; no unfavorable rows were discarded.

## Uncertainty

'''+ '\n'.join(interval)+'''

Intervals use 2,000 paired circular sequence-block bootstrap draws, block length 10 and a fixed seed. Both cross zero. Blocks include gaps among observed residues, so this only approximates local dependence; it is not evidence of independent proteins or a population-level effect. The two temperatures share the same mutation and should not be pooled as independent validation.

## Additional observable, secondary only

We had not analyzed the 25°C CA/CB shifts in the first pilot. Using paired CA and CB observations at 136 eligible residues, the primary distance refinement gives 0.4105→0.4217 in rank correlation. This is an additional nucleus/observable check on the same samples, not a new cohort. The carbon target is sqrt(delta_CA²+delta_CB²); it has different physical meaning from the H/N target. No method was selected based on the carbon result.

| Score | n | Carbon Spearman | Top-quarter precision |
|---|---:|---:|---:|
'''+ '\n'.join(carbon)+'''

## Compute and validation

No new weights, model training, API calls, MD or paid compute. The first complete calculation including bootstrap took 10.15 seconds on the user's M2 Air, with 117,817,344 bytes maximum resident set size. Timings vary with imports/cache; this is a host observation, not a hardware benchmark. Source-only score generation reads 6BTA; target-mutant coordinates and NMR values are not score inputs. Scores are joined to measurements only for evaluation. Numeric replay and score-byte hashes are in `verification.json`; tests check rigid-transform invariance, network permutation invariance, isolated nodes and rejection of the post-mutation structure.

## What we learned

Using atom locations nearer the measured nuclei can help a little over CA-to-CA distance. We have not established allosteric prediction or utility beyond this mutant pair, and the point gain is within uncertainty. The appropriate next validation would use separately specified mutation contrasts or other proteins with matched public NMR data, without selecting them because this score works. Do not keep tuning on this pair and present it as fresh evidence.

Contact-based perturbation propagation and mutation-response networks predate this pilot: [pPerturb](https://pmc.ncbi.nlm.nih.gov/articles/PMC6977024/) and [fast mutation-response scanning](https://pmc.ncbi.nlm.nih.gov/articles/PMC8067912/). The experimental rescue mechanism is from [Otten et al.](https://doi.org/10.1038/s41467-018-03562-9). Our work is AI-assisted implementation and transparent baseline evaluation.

## Reproduce

From the parent pilot directory, using the pinned environment:

```sh
python -m unittest discover -s followup -p 'test_followup.py' -v
python followup/verify.py
python followup/report.py
```
'''
(ROOT/'REPORT.md').write_text(text)
fig,ax=plt.subplots(figsize=(7.2,3.8),layout='constrained')
for x,temp in enumerate((283,298)):
    m=s['metrics'][f'{temp}K_exclude_2']['csp']
    ax.bar(x-.17,m['geometry_locality']['spearman'],.32,color='#176b87',label='CA distance' if x==0 else None)
    ax.bar(x+.17,m['sidechain_amide_distance']['spearman'],.32,color='#dc7e23',label='Side-chain to backbone N' if x==0 else None)
ax.set(xticks=[0,1],xticklabels=['10°C (n=154)','25°C (n=152)'],ylim=(0,.7),ylabel='Spearman correlation',
       title='Small retrospective gain on one mutant pair\nBoth descriptive improvement intervals cross zero')
ax.legend(frameon=False)
fig.savefig(ROOT/'comparison.png',dpi=180)
plt.close(fig)
