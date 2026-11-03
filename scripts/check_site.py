"""Check the publication directory without network access or third-party packages."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
DOCS=ROOT/"docs"


class Page(HTMLParser):
    def __init__(self,path):
        super().__init__();self.links=[];self.ids=set();self.language=None
        self.feed(path.read_text())
    def handle_starttag(self,tag,attrs):
        for key,val in attrs:
            if key in ("href","src"):self.links.append(val)
            if key=="id":self.ids.add(val)
            if tag=="html" and key=="lang":self.language=val


def require(condition,message):
    if not condition:raise RuntimeError(message)


def main():
    allowed={".html",".css",".js",".json",".svg",".png",".csv",".py",".txt",".md"}
    for path in DOCS.rglob("*"):
        require(not path.is_symlink(),f"Publication cannot contain symlinks: {path}")
        if path.is_file():
            require(path.suffix in allowed,f"Unexpected publication file: {path}")
            require(not any(p.startswith('.') for p in path.relative_to(DOCS).parts),f"Hidden publication file: {path}")
    pages={p.resolve():Page(p) for p in DOCS.rglob("*.html")}
    for path,page in pages.items():
        require(page.language=="en",f"Expected English page: {path}")
        for href in page.links:
            url=urlsplit(href)
            if url.scheme or url.netloc:continue
            require(not url.path.startswith('/'),f"Use project-relative links for Pages: {href}")
            target=(path.parent/unquote(url.path)).resolve() if url.path else path
            require(DOCS.resolve() in [target,*target.parents],f"Link leaves publication directory: {href}")
            if target.is_dir():target=target/'index.html'
            require(target.is_file(),f"Missing target from {path}: {href}")
            if url.fragment:
                require(target in pages and url.fragment in pages[target].ids,f"Missing anchor from {path}: {href}")
    catalog=json.loads((DOCS/'catalog.json').read_text())
    require(len({p['id'] for p in catalog})==len(catalog),'Duplicate paper identifiers')
    for paper in catalog:
        for reference in paper.get('resources',[]):
            url=urlsplit(reference['href']);target=(DOCS/unquote(url.path)).resolve()
            require(target in pages and (not url.fragment or url.fragment in pages[target].ids),f"Missing catalog reference: {reference['href']}")
        for exhibit in paper.get('exhibits',[]):
            target=(DOCS/exhibit['href']).resolve()
            if target.is_dir():target=target/'index.html'
            require(target in pages,f"Catalog exhibit missing: {target}")
    reports=[('field-1965','validation.json','config.json'),
             ('field-1965','validation-wave.json','config-wave.json'),
             ('field-1965','validation-map.json','config-map.json'),
             ('fgh-1969','validation.json','config.json'),
             ('fgh-1969','validation-cold.json','config-cold.json'),
             ('fgh-1969','validation-two-phase.json','config-two-phase.json'),
             ('fgh-1969','validation-comparison.json','config-comparison.json'),
             ('takayanagi-nishimura-1960','validation.json','config.json'),
             ('iwasaki-2012','validation.json','config.json'),
             ('schwarz-1972','validation.json','config.json')]
    total=0
    for chapter,report_name,config_name in reports:
        directory=DOCS/chapter
        report=json.loads((directory/report_name).read_text())
        require(bool(report['checks']) and all(c['passed'] for c in report['checks']),f"Validation failed: {report_name}")
        total+=len(report['checks'])
        hashes=report['code_sha256']
        if isinstance(hashes,str):
            hashes={('reproduce-map.py' if 'map' in report_name else 'reproduce.py'):hashes}
        for name,sha in hashes.items():
            require(hashlib.sha256((directory/name).read_bytes()).hexdigest()==sha,f"Stale downloadable source: {name}")
        require(hashlib.sha256((directory/config_name).read_bytes()).hexdigest()==report['config_sha256'],f"Stale configuration: {config_name}")
    for chapter in {r[0] for r in reports}:
        directory=DOCS/chapter
        for path in (ROOT/'experiments'/chapter).iterdir():
            if path.name=='publish.py' or path.suffix not in ('.py','.json','.txt','.md'):continue
            require((directory/path.name).is_file(),f"Download missing: {chapter}/{path.name}")
            require(path.read_bytes()==(directory/path.name).read_bytes(),f"Download differs from source: {chapter}/{path.name}")
    print(f"Site check passed: {len(pages)} English pages; {total} scientific checks; links, catalog, and source hashes valid.")
    print("Only docs/ will be uploaded; PDFs and private source material are excluded.")


if __name__=='__main__':main()
