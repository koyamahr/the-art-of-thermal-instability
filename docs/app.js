"use strict";
const list = document.querySelector("#papers");
const statusLine = document.querySelector("#catalog-status");
const filters = [...document.querySelectorAll("[data-filter]")];
let papers = [];
function element(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text) node.textContent = text;
  return node;
}
function render(filter = "all") {
  const visible = papers.filter(p => filter === "all" || p.group === filter);
  list.replaceChildren();
  for (const p of visible) {
    const row = element("article", "paper");
    row.id = p.id;
    row.append(element("div", "year", String(p.year)));
    const body = element("div");
    body.append(element("span", "state", p.status === "validated" ? "Validated locally" : "Planned"), element("p", "author", p.authors), element("h3", "", p.title), element("p", "topic", p.topic));
    const details = element("details");
    details.append(element("summary", "", p.exhibit ? "About this exhibit" : "View planned exhibit"), element("p", "", p.scope));
    body.append(details);
    for (const entry of p.exhibits || (p.exhibit ? [{href:p.exhibit,label:"Explore exhibit →"}] : [])) {
      const link = element("a", "source", entry.label);
      link.href = entry.href;
      body.append(link);
    }
    const link = element("a", "source", "Paper & bibliography ↗");
    link.href = p.source;
    link.target = "_blank";
    link.rel = "noopener noreferrer";
    body.append(link);
    row.append(body);
    list.append(row);
  }
  const ready = visible.reduce((count,p) => count + (p.status === "validated" ? (p.exhibits ? p.exhibits.length : p.exhibit ? 1 : 0) : 0), 0);
  statusLine.textContent = `${visible.length} papers · ${ready} exhibit${ready === 1 ? "" : "s"} validated locally`;
}
filters.forEach(button => button.addEventListener("click", () => {
  filters.forEach(b => b.setAttribute("aria-pressed", String(b === button)));
  render(button.dataset.filter);
}));
fetch("catalog.json").then(response => {
  if (!response.ok) throw new Error("catalog unavailable");
  return response.json();
}).then(data => { papers = data; render(); }).catch(() => {
  statusLine.textContent = "The paper catalog could not be loaded. Open this site through an HTTP server, or use the JSON link in the footer.";
  filters.forEach(button => { button.disabled = true; });
});
