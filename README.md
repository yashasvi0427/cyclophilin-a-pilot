# Cyclophilin A: real experimental data, simple baselines first

This is an AI-assisted personal research pilot using published experiments. It does not claim the published mutation-rescue mechanism as a new discovery. Start with [the measured report](results/REPORT.md) and [the frozen scope](PROTOCOL.md).

The [small structure-based follow-up](followup/REPORT.md) tests four fixed refinements. The primary side-chain-to-backbone-N distance gives a small point gain over CA distance (0.359→0.369 at 10°C; 0.523→0.540 at 25°C), but both descriptive intervals cross zero. It is not a robust or independently validated win. All candidates, including failures, are retained. No additional model download or fitting was used.

Question: when I97V is added to the CypA S99T/C115S rescue mutant, can cheap scores identify residues whose measured backbone NMR environments change? We compare sequence distance, C-alpha distance in the source structure and a frozen ESM2 8M last-layer embedding-change score. No fitting, tuning or model training.

The 3D distance baseline outperformed the tested ESM2 score on this single mutant pair at both temperatures. This result does not establish anything about coordinated motion, model-generated ensembles or enzyme activity. Two temperatures are not two independent proteins. BioEmu2 was not run because the released package requires Linux/NVIDIA and this pilot ran on an 8 GB M2 Air.

## Reproduce

Python 3.12 was used. Install `requirements.lock.txt` in an isolated environment. Keep the model cache outside the repository.

```sh
python -m unittest -v test_pilot
python verify.py
python make_report.py
```

`verify.py` uses the saved 165-by-320 embeddings and needs no model download. It checks raw-data/protocol hashes and demands byte-identical residue scores and identical metrics. To regenerate the embeddings, install torch/fair-esm and run:

```sh
python pilot.py --esm --model-cache /path/to/external/model-cache
```

This downloads about 30 MB of open weights, plus a small contact-regression file from Meta's official release. Model weights remain outside Git; their hash is recorded. The source files are cached with original bytes, URLs and hashes. Deleting a data file causes `pilot.py` to fetch its official source again; changed upstream bytes will fail the receipt check in `verify.py`.

## Files

- `pilot.py`: data alignment, identity checks, frozen scores, metrics and machine-readable receipts.
- `test_pilot.py`: known deposited values, temperatures, structure-identity failure, metrics and tie handling.
- `verify.py`: source and result replay.
- `results/summary.json`, `residue_scores.csv`: all measurements, exclusions, sensitivities and model metadata.
- `results/esm_embeddings.npz`: cached derived embeddings, not weights or experimental measurements.
- Numerical verification receipts are published; machine-specific baseline run logs stay in the private continuation.
- `DATA_SOURCES.md`: attribution and usage notices.

## Contribution and next-step limits

The user's contributions are the research direction, biological interests, scope decisions and portfolio goals. Codex assisted literature search, implementation, analysis and documentation. The user should review and explain the work before representing personal expertise in interviews. No new wet-lab experiment or validated biological mechanism was produced.

The original proposal concerns coordinated side-chain packing across the damage/rescue series. This first measurement-localization test does not answer that full question. It establishes a cheap baseline, a checked data path and a negative result for one representation score. A future all-atom comparison needs a separately frozen protocol and appropriate experimental observables; a four-variant case study cannot support a broad benchmark claim.

## Public review export

This is the public project version, exported without private conversation history, application files, resume, credentials or employer material. See [release inventory](PUBLIC_RELEASE.json). The original negative ESM2 result and the small uncertain geometry gain are both preserved.
