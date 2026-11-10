"""Recalculate scores from deposited data and saved embeddings; fail on drift."""
import hashlib
import json
from pathlib import Path

import numpy as np

from pilot import ROOT, run, sha

out=ROOT/'results'
before=json.loads((out/'summary.json').read_text())
csv_before=(out/'residue_scores.csv').read_bytes()
assert sha(ROOT/'PROTOCOL.md')=='7995a760c3a57f2a6029c52ec7c76d6bf677e9233aefbbb65fd2b25d54c2ce74'
for name,metadata in before['data'].items():
    assert sha(ROOT/'data'/name)==metadata['sha256'], name
e=np.load(out/'esm_embeddings.npz')
assert e['double'].shape==e['triple'].shape==(165,320)
assert np.isfinite(e['double']).all() and np.isfinite(e['triple']).all()
after=run(cached_esm=True)
assert before['metrics']==after['metrics']
assert before['matching']==after['matching']
assert before['structure_audit']==after['structure_audit']
assert csv_before==(out/'residue_scores.csv').read_bytes()
checks={'protocol_hash_verified':True,'all_source_hashes_verified':True,'score_csv_byte_identical':True,
        'metrics_identical':True,'embedding_shape':[165,320],
        'score_csv_sha256':sha(out/'residue_scores.csv')}
(out/'verification.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
