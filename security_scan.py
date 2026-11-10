"""Bounded credential-pattern scan, complementary to manual file review."""
import json
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parent
PATTERNS={
 'private_key':r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',
 'github_token':r'\b(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{40,})\b',
 'aws_access_key':r'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b',
 'api_token':r'\bsk-(?:proj-|ant-)?[A-Za-z0-9_-]{24,}\b',
 'employer_internal_url':r'https?://[^\s/]*(?:salesforce-internal|sfdcinternal|\.corp\.)[^\s]*',
}
findings=[]
files=[]
for path in ROOT.rglob('*'):
 if not path.is_file() or '__pycache__' in path.parts or path.suffix in ('.png','.npz'):
  continue
 files.append(str(path.relative_to(ROOT)))
 text=path.read_text(errors='replace')
 for name,pattern in PATTERNS.items():
  if re.search(pattern,text):findings.append({'path':str(path.relative_to(ROOT)),'type':name})
receipt={
 'scope':'new cyclophilin_a directory only; raw public STAR/CIF included',
 'text_files_scanned':len(files),'patterns':list(PATTERNS),'findings':findings,
 'manual_review':'No employer code, internal logs, patient data, resume PDF, personal application answers, credentials or model weights added. Public deposit authors/contact metadata retained. Private conversation history and application files are excluded from this public export.',
 'limitation':'Pattern scans and manual review do not prove absence of every secret or proprietary detail.'}
(ROOT/'results/security_review.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
if findings:raise SystemExit(1)
