"""Copy reproducible cover assets, or verify the committed asset manifest."""
from pathlib import Path
import json,hashlib,shutil,sys
HERE=Path(__file__).resolve().parent
DOCS=HERE.parents[1]/'docs/ti-2d'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    manifest=json.loads((DOCS/'manifest.json').read_text())
    for name,digest in manifest['sha256'].items():assert sha(DOCS/name)==digest,name
    print('Cover asset manifest verified; this does not rerun the full movie.')
if '--verify-only' in sys.argv:verify()
else:
    DOCS.mkdir(exist_ok=True)
    for p in HERE.iterdir():
        if p.suffix in ('.py','.cpp','.json','.txt','.md') and p.name!='publish.py':shutil.copy(p,DOCS/p.name)
    for p in (HERE/'output').iterdir():
        if p.is_file():shutil.copy(p,DOCS/p.name)
    names=[p for p in DOCS.iterdir() if p.is_file() and p.name not in ('index.html','manifest.json')]
    (DOCS/'manifest.json').write_text(json.dumps(dict(sha256={p.name:sha(p) for p in sorted(names)}),indent=2)+'\n');verify()
