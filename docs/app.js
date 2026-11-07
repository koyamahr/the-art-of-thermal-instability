"use strict";
const list = document.querySelector("#papers");
const statusLine = document.querySelector("#catalog-status");
const filters = [...document.querySelectorAll("[data-filter]")];
const staticRows = [...list.querySelectorAll(".paper[data-group]")];
let papers = [];
function element(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text) node.textContent = text;
  return node;
}
function render(filter = "all") {
  if (staticRows.length) {
    const visible = staticRows.filter(row => filter === "all" || row.dataset.group === filter);
    staticRows.forEach(row => { row.hidden = !visible.includes(row); });
    const ready = visible.reduce((count, row) => count + Number(row.dataset.readyCount), 0);
    statusLine.textContent = `${visible.length} references · oldest first · ${ready} exhibit${ready === 1 ? "" : "s"} numerically validated`;
    return;
  }
  const visible = papers.filter(p => filter === "all" || p.group === filter);
  list.replaceChildren();
  for (const p of visible) {
    const row = element("article", "paper");
    row.id = p.id;
    row.append(element("div", "year", String(p.year)));
    const body = element("div");
    const heading = element("h3");
    const overview = element("a", "paper-title", p.title);
    overview.href = `papers/${p.id}/`;
    heading.append(overview);
    body.append(element("span", "state", p.status === "validated" ? "Numerically validated" : p.status === "deferred" ? "Historical match unresolved · deferred" : p.status === "in_progress" ? "Equilibrium model · in progress" : p.status === "prototype" ? "Computed prototype · review pending" : p.status === "audited" ? "Source rates audited · truncated model" : p.status === "listed" ? "Reference listed · scope proposed" : "Planned"), element("p", "author", p.authors), heading, element("p", "topic", p.topic));
    const details = element("details");
    details.append(element("summary", "", p.status === "listed" ? "Proposed scope" : p.exhibit || p.exhibits ? "About this exhibit" : "View planned exhibit"), element("p", "", p.scope));
    body.append(details);
    for (const entry of p.exhibits || (p.exhibit ? [{href:p.exhibit,label:"Explore exhibit →"}] : [])) {
      const link = element("a", "source", entry.label);
      link.href = entry.href;
      body.append(link);
    }
    for (const resource of p.resources || []) {
      const reference = element("a", "source", resource.label);
      reference.href = resource.href;
      body.append(reference);
    }
    const link = element("a", "source", p.reference_type === "book" ? "Book & bibliography ↗" : "Paper & bibliography ↗");
    link.href = p.source;
    link.target = "_blank";
    link.rel = "noopener noreferrer";
    body.append(link);
    row.append(body);
    list.append(row);
  }
  const ready = visible.reduce((count,p) => count + (p.exhibits ? p.exhibits.filter(e => e.validated ?? p.status === "validated").length : p.status === "validated" && p.exhibit ? 1 : 0), 0);
  statusLine.textContent = `${visible.length} references · oldest first · ${ready} exhibit${ready === 1 ? "" : "s"} numerically validated`;
}
filters.forEach(button => button.addEventListener("click", () => {
  filters.forEach(b => b.setAttribute("aria-pressed", String(b === button)));
  render(button.dataset.filter);
}));
if (staticRows.length) render();
else fetch("catalog.json").then(response => {
  if (!response.ok) throw new Error("catalog unavailable");
  return response.json();
}).then(data => { papers = [...data].sort((a, b) => a.year - b.year); render(); }).catch(() => {
  statusLine.textContent = "The paper catalog could not be loaded. Open this site through an HTTP server, or use the JSON link in the footer.";
  filters.forEach(button => { button.disabled = true; });
});
