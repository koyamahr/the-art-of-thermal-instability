"""Check the publication directory without network access or third-party packages."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
from argparse import ArgumentParser
from xml.etree import ElementTree as ET
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
DOCS=ROOT/"docs"


class Page(HTMLParser):
    def __init__(self,path):
        super().__init__();self.links=[];self.ids=set();self.language=None
        self.canonicals=[];self.metadata={};self.title='';self.in_title=False;self.paper_ids=[]
        self.feed(path.read_text())
    def handle_starttag(self,tag,attrs):
        attributes=dict(attrs)
        if tag=='title':self.in_title=True
        if tag=='link' and attributes.get('rel')=='canonical':self.canonicals.append(attributes.get('href'))
        if tag=='meta':self.metadata[attributes.get('name',attributes.get('property'))]=attributes.get('content','')
        if tag=='article' and 'data-group' in attributes:self.paper_ids.append(attributes.get('id'))
        for key,val in attrs:
            if key in ("href","src","poster","data-src"):self.links.append(val)
            if key=="srcset":self.links.extend(item.strip().split()[0] for item in val.split(',') if item.strip())
            if key=="id":self.ids.add(val)
            if tag=="html" and key=="lang":self.language=val
    def handle_endtag(self,tag):
        if tag=='title':self.in_title=False
    def handle_data(self,data):
        if self.in_title:self.title+=data


def require(condition,message):
    if not condition:raise RuntimeError(message)


def main():
    allowed={".html",".css",".js",".json",".svg",".png",".csv",".py",".txt",".md",".cpp",".mp4",".npz",".xml"}
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
    registry=json.loads((DOCS/'equilibrium-landmarks.json').read_text())
    for name,sha in registry['reports_sha256'].items():
        require(hashlib.sha256((DOCS/name).read_bytes()).hexdigest()==sha,f'Stale equilibrium landmark registry: {name}')
    catalog=json.loads((DOCS/'catalog.json').read_text())
    require(len({p['id'] for p in catalog})==len(catalog),'Duplicate paper identifiers')
    if DOCS.resolve() != (ROOT/'docs').resolve():
        from build_site import SITE_URL, canonical_url
        sitemap=ET.parse(DOCS/'sitemap.xml')
        urls=[node.text for node in sitemap.findall('{http://www.sitemaps.org/schemas/sitemap/0.9}url/{http://www.sitemaps.org/schemas/sitemap/0.9}loc')]
        expected_urls={canonical_url(path.relative_to(DOCS.resolve()).as_posix()) for path in pages}
        require(len(urls)==len(set(urls)) and set(urls)==expected_urls,'Sitemap does not match generated pages')
        require(set(pages[(DOCS/'index.html').resolve()].paper_ids)=={p['id'] for p in catalog},'Homepage must contain every reference without JavaScript')
        for path,page in pages.items():
            canonical=canonical_url(path.relative_to(DOCS.resolve()).as_posix())
            require(page.canonicals==[canonical],f'Missing or duplicate canonical: {path}')
            require(page.title and page.metadata.get('description'),f'Missing title or description: {path}')
            require(page.metadata.get('og:title')==page.title,f'Incorrect Open Graph title: {path}')
            require(page.metadata.get('og:url')==canonical,f'Incorrect Open Graph URL: {path}')
            require(page.metadata.get('og:description')==page.metadata['description'],f'Incorrect Open Graph description: {path}')
        for paper in catalog:
            path=(DOCS/'papers'/paper['id']/'index.html').resolve()
            require(path in pages,f'Missing paper overview: {paper["id"]}')
            html=path.read_text()
            from html import escape
            require(escape(paper['scope']) in html and escape(paper['topic']) in html,f'Paper scope missing: {paper["id"]}')
    for paper in catalog:
        for reference in paper.get('resources',[]):
            url=urlsplit(reference['href']);target=(DOCS/unquote(url.path)).resolve()
            require(target in pages and (not url.fragment or url.fragment in pages[target].ids),f"Missing catalog reference: {reference['href']}")
        for exhibit in paper.get('exhibits',[]):
            target=(DOCS/exhibit['href']).resolve()
            if target.is_dir():target=target/'index.html'
            require(target in pages,f"Catalog exhibit missing: {target}")
    reports=[('ti-2d','validation.json','config.json'),('ki-2002','validation.json','config.json'),('ki-2002','validation-potential.json','config-potential.json'),
             ('fgh-1969','validation-saturation.json','config-saturation.json'),
             ('field-1965','validation.json','config.json'),
             ('field-1965','validation-wave.json','config-wave.json'),
             ('field-1965','validation-map.json','config-map.json'),
             ('fgh-1969','validation.json','config.json'),
             ('fgh-1969','validation-cold.json','config-cold.json'),
             ('fgh-1969','validation-two-phase.json','config-two-phase.json'),
             ('fgh-1969','validation-comparison.json','config-comparison.json'),
             ('takayanagi-nishimura-1960','validation.json','config.json'),
             ('iwasaki-2012','validation.json','config.json'),
             ('schwarz-1972','validation.json','config.json'),
             ('schwarz-1972','validation-nonlinear.json','config-nonlinear.json')]
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
            if path.name=='publish.py' or path.suffix not in ('.py','.json','.txt','.md','.cpp'):continue
            require((directory/path.name).is_file(),f"Download missing: {chapter}/{path.name}")
            require(path.read_bytes()==(directory/path.name).read_bytes(),f"Download differs from source: {chapter}/{path.name}")
    print(f"Site check passed: {len(pages)} English pages; {total} scientific checks; links, catalog, and source hashes valid.")
    print(f"Checked publication directory: {DOCS.name}/")


if __name__=='__main__':
    parser=ArgumentParser()
    parser.add_argument('--directory',type=Path,default=DOCS)
    DOCS=parser.parse_args().directory.resolve()
    main()
