"""Build a crawlable archive from the public docs/ inputs (standard library only)."""
from html import escape
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit
from xml.etree import ElementTree as ET
import json
import posixpath
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
SITE_URL = "https://koyamahr.github.io/the-art-of-thermal-instability/"
STATUS = {
    "validated": "Numerically validated",
    "deferred": "Historical match unresolved · deferred",
    "in_progress": "Equilibrium model · in progress",
    "prototype": "Computed prototype · review pending",
    "audited": "Source rates audited · truncated model",
    "listed": "Reference listed · scope proposed",
    "planned": "Planned",
}


def paper_path(paper):
    return f"papers/{paper['id']}/"


def exhibits(paper):
    return paper.get("exhibits", [{"href": paper["exhibit"], "label": "Explore exhibit →"}] if paper.get("exhibit") else [])


def validated_count(paper):
    return sum(bool(entry.get("validated", paper["status"] == "validated")) for entry in exhibits(paper))


def local_link(href, current):
    """Rebase a catalog URL from the docs root without changing its fragment."""
    url = urlsplit(href)
    if url.scheme or url.netloc:
        return href
    path = posixpath.relpath(url.path or ".", posixpath.dirname(current) or ".")
    if url.path.endswith("/"):
        path += "/"
    return urlunsplit(("", "", path, url.query, url.fragment))


def catalog_row(paper):
    links = "".join(f'<a class="source" href="{escape(e["href"], quote=True)}">{escape(e["label"])}</a>' for e in exhibits(paper) + paper.get("resources", []))
    summary = "Proposed scope" if paper["status"] == "listed" else "About this exhibit" if exhibits(paper) else "View planned exhibit"
    source_label = "Book & bibliography ↗" if paper.get("reference_type") == "book" else "Paper & bibliography ↗"
    return (f'<article class="paper" id="{paper["id"]}" data-group="{escape(paper["group"], quote=True)}" data-ready-count="{validated_count(paper)}">'
            f'<div class="year">{paper["year"]}</div><div><span class="state">{STATUS[paper["status"]]}</span>'
            f'<p class="author">{escape(paper["authors"])}</p><h3><a class="paper-title" href="{paper_path(paper)}">{escape(paper["title"])}</a></h3>'
            f'<p class="topic">{escape(paper["topic"])}</p><details><summary>{summary}</summary><p>{escape(paper["scope"])}</p></details>'
            f'{links}<a class="source" href="{escape(paper["source"], quote=True)}" target="_blank" rel="noopener noreferrer">{source_label}</a></div></article>')


def paper_page(paper, catalog):
    current = paper_path(paper) + "index.html"
    link = lambda href: escape(local_link(href, current), quote=True)
    title = f'{paper["title"]} — {paper["authors"]} ({paper["year"]}) | Thermal Instability Archive'
    description = f'{paper["authors"]} ({paper["year"]}). {paper["topic"]} {paper["scope"]}'[:280]
    entries = exhibits(paper)
    exhibit_section = ''
    if entries:
        exhibit_section = '<section id="exhibits"><h2>Figures, calculations and validation.</h2><ul>' + ''.join(
            f'<li><a href="{link(entry["href"])}">{escape(entry["label"])}</a></li>' for entry in entries
        ) + '</ul></section>'
    resource_section = ''
    if paper.get('resources'):
        resource_section = '<section><h2>Source notes and related resources.</h2><ul>' + ''.join(
            f'<li><a href="{link(entry["href"])}">{escape(entry["label"])}</a></li>' for entry in paper['resources']
        ) + '</ul></section>'
    same_path = [other for other in catalog if other['group'] == paper['group'] and other['id'] != paper['id']]
    related = ''
    if same_path:
        related = '<section><h2>Continue along this research path.</h2><ul>' + ''.join(
            f'<li>{other["year"]} · <a href="{link(paper_path(other))}">{escape(other["authors"])} — {escape(other["title"])}</a></li>' for other in same_path
        ) + '</ul></section>'
    review = f'<p><strong>Source review:</strong> {escape(paper["review_status"])}</p>' if paper.get('review_status') else ''
    doi = f'<p><a href="https://doi.org/{escape(paper["doi"], quote=True)}" target="_blank" rel="noopener noreferrer">Publication / DOI ↗</a></p>' if paper.get('doi') else ''
    research_path = f'<p><a href="{link("research-paths.html#" + paper["id"])}">Read this reference in its research path →</a></p>' if paper.get('reference_section') else ''
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title><meta name="description" content="{escape(description, quote=True)}">
<link rel="stylesheet" href="{link('styles.css')}"><link rel="stylesheet" href="{link('field-1965/exhibit.css')}"></head>
<body><header><a class="brand" href="{link('./')}">The Art of Thermal Instability</a><nav aria-label="Paper navigation"><a href="{link('./#atlas')}">Papers &amp; Figures</a><a href="#scope">Scope</a><a href="#source">Original source</a></nav></header>
<main class="exhibit"><section class="hero"><p class="eyebrow">Paper and book atlas / {paper['year']}</p>
<h1>{escape(paper['title'])}</h1><p class="intro">{escape(paper['authors'])} · {paper['year']}</p>
{f'<p>{escape(paper["journal"])}</p>' if paper.get('journal') else ''}
<p class="state">{STATUS[paper['status']]}</p><p>{escape(paper['topic'])}</p></section>
<section id="scope"><h2>Scope in this archive.</h2><p>{escape(paper['scope'])}</p>{review}</section>
{exhibit_section}{resource_section}
<section id="source"><h2>Original publication and bibliography.</h2><p><a href="{escape(paper['source'], quote=True)}" target="_blank" rel="noopener noreferrer">Read the original source and bibliography ↗</a></p>{doi}{research_path}</section>
{related}</main><footer><p>The Art of Thermal Instability</p><a href="{link('./#' + paper['id'])}">Back to this reference in the paper atlas</a><a href="{link('research-paths.html')}">Research paths</a></footer></body></html>
'''


class Metadata(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.title = ''
        self.description = ''
        self.in_title = False
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'title':
            self.in_title = True
        if tag == 'meta' and attrs.get('name') == 'description':
            self.description = attrs.get('content', '')

    def handle_endtag(self, tag):
        if tag == 'title':
            self.in_title = False

    def handle_data(self, data):
        if self.in_title:
            self.title += data


def canonical_url(relative):
    return SITE_URL + (relative[:-10] if relative.endswith('index.html') else relative)


def metadata(html, relative):
    page = Metadata(html)
    if not page.title or not page.description:
        raise ValueError(f'Missing title or description: {relative}')
    canonical = canonical_url(relative)
    tags = [f'<link rel="canonical" href="{canonical}">'] + [
        f'<meta property="{key}" content="{escape(value, quote=True)}">' for key, value in {
            'og:title': page.title,
            'og:description': page.description,
            'og:url': canonical,
            'og:type': 'website' if relative == 'index.html' else 'article',
        }.items()
    ]
    return html.replace('</head>', '\n' + '\n'.join(tags) + '\n</head>', 1)


def build(source=ROOT / 'docs', destination=ROOT / 'dist'):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if source == destination or source in destination.parents or destination in source.parents:
        raise ValueError('Build output must be separate from the source directory')
    catalog = sorted(json.loads((source / 'catalog.json').read_text()), key=lambda p: p['year'])
    ids = [p['id'] for p in catalog]
    if len(ids) != len(set(ids)) or any(not re.fullmatch(r'[a-z0-9]+(?:-[a-z0-9]+)*', key) for key in ids):
        raise ValueError('Duplicate or unsafe catalog identifier')
    if any(path.is_symlink() for path in source.rglob('*')):
        raise ValueError('The public source must not contain symlinks')
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination)

    index = (destination / 'index.html').read_text()
    rows = '\n'.join(catalog_row(paper) for paper in catalog)
    index, replacements = re.subn(r'<div id="papers">[\s\S]*?</div>', lambda _: f'<div id="papers">{rows}</div>', index, count=1)
    if replacements != 1:
        raise ValueError('Cannot locate the paper catalog container')
    ready = sum(validated_count(paper) for paper in catalog)
    index = index.replace('Loading the paper catalog.', f'{len(catalog)} references · oldest first · {ready} exhibits numerically validated')
    index = index.replace('<noscript><p>Enable JavaScript to browse the paper catalog.', '<noscript><p>All references are listed above. JavaScript enables topic filtering.')
    (destination / 'index.html').write_text(index)

    for paper in catalog:
        path = destination / paper_path(paper) / 'index.html'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(paper_page(paper, catalog))
        for entry in exhibits(paper):
            relative = urlsplit(entry['href']).path
            if relative.endswith('/'):
                relative += 'index.html'
            exhibit = destination / relative
            html = exhibit.read_text()
            context = f'<p class="paper-context"><a href="{escape(local_link(paper_path(paper), relative), quote=True)}">Paper overview: {escape(paper["authors"])} ({paper["year"]}) →</a></p>'
            html = re.sub(r'(<main\b[^>]*>)', lambda match: match[0] + context, html, count=1)
            exhibit.write_text(html)

    sitemap = ET.Element('urlset', xmlns='http://www.sitemaps.org/schemas/sitemap/0.9')
    pages = sorted(destination.rglob('*.html'))
    for path in pages:
        relative = path.relative_to(destination).as_posix()
        html = path.read_text()
        # Link each bibliography entry in the grouped reading guide to its overview.
        if relative == 'research-paths.html':
            for paper in catalog:
                if paper.get('reference_section'):
                    html = html.replace(f'<h3>{escape(paper["title"])}</h3>', f'<h3><a href="{paper_path(paper)}">{escape(paper["title"])}</a></h3>')
        path.write_text(metadata(html, relative))
        ET.SubElement(ET.SubElement(sitemap, 'url'), 'loc').text = canonical_url(relative)
    ET.ElementTree(sitemap).write(destination / 'sitemap.xml', encoding='utf-8', xml_declaration=True)
    print(f'Built {len(pages)} static English pages, including {len(catalog)} paper overviews, and sitemap.xml.')
    return len(pages)


if __name__ == '__main__':
    build()
