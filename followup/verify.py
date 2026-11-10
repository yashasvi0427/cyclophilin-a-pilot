import json
from run import ROOT,run,sha

before=json.loads((ROOT/'results.json').read_text())
hashes={name:sha(ROOT/name) for name in ('scores.csv','carbon_scores.csv')}
after=run()
before.pop('elapsed_seconds');after.pop('elapsed_seconds')
assert before==after,'Follow-up numerical outputs drifted'
assert hashes=={name:sha(ROOT/name) for name in hashes},'Score bytes drifted'
receipt={'all_non_timing_results_identical':True,'score_csv_byte_identical':True,
         'sha256':hashes,'protocol_sha256':after['protocol_sha256']}
(ROOT/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
