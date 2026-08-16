"""Render a :class:`~bbl.report.build.DesignReport` as one self-contained HTML file.

No JavaScript, no external stylesheet, no image files: the figures are inline SVG, so the
report is a single artefact that can be emailed, committed, or opened from a cluster share
years later and still look the same. It prints sensibly too -- a notebook entry gets taped
into a book about as often as it gets read on screen.
"""

from __future__ import annotations

import html as _html
import json
import re

#: A sequence-looking string is set in a monospace face and allowed to wrap anywhere.
_SEQUENCE = re.compile(r"^5?'?-?[ACGT]{8,}-?3?'?$")

STYLESHEET = """
:root {
  --ink: #1c1c1c; --muted: #6b6b6b; --rule: #e0e0e0; --panel: #f7f7f5;
  --accent: #7a1f1f; --warn: #8a5a00; --warn-bg: #fdf6e6;
}
* { box-sizing: border-box; }
body {
  font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  color: var(--ink); margin: 0; padding: 2.5rem 1.5rem 4rem; background: #fff;
}
main { max-width: 62rem; margin: 0 auto; }
h1 { font-size: 1.6rem; margin: 0 0 .2rem; letter-spacing: -.01em; }
h2 {
  font-size: 1.1rem; margin: 2.6rem 0 .9rem; padding-bottom: .35rem;
  border-bottom: 2px solid var(--ink); text-transform: uppercase; letter-spacing: .06em;
}
h3 { font-size: 1rem; margin: 1.8rem 0 .5rem; }
p, li { margin: .45rem 0; }
.subtitle { color: var(--muted); margin: 0 0 1.4rem; font-size: 1.02rem; }
.meta { display: grid; grid-template-columns: repeat(auto-fit, minmax(13rem, 1fr));
        gap: .1rem 1.4rem; padding: .9rem 1.1rem; background: var(--panel);
        border: 1px solid var(--rule); border-radius: 6px; margin-bottom: 1.6rem; }
.meta div { font-size: .88rem; }
.meta span { display: block; color: var(--muted); text-transform: uppercase;
             letter-spacing: .05em; font-size: .7rem; }
table { border-collapse: collapse; width: 100%; margin: .5rem 0 .3rem; font-size: .86rem; }
caption { text-align: left; font-weight: 600; padding: .5rem 0 .3rem; font-size: .92rem; }
th, td { border: 1px solid var(--rule); padding: .35rem .55rem; text-align: left;
         vertical-align: top; }
th { background: var(--panel); font-weight: 600; }
tbody tr:last-child td { font-weight: 600; }
.note { font-size: .82rem; color: var(--muted); margin: .1rem 0 1.2rem; }
.seq { font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace; font-size: .8rem;
       word-break: break-all; }
.figures { display: grid; grid-template-columns: 1fr 1fr; gap: 1.2rem; align-items: start; }
.figure { margin: 0 0 1rem; }
.figure.wide { grid-column: 1 / -1; }
.figure svg { width: 100%; height: auto; display: block; }
figcaption { font-size: .8rem; color: var(--muted); margin-top: .3rem; }
.step { border-left: 3px solid var(--rule); padding: 0 0 .1rem 1.1rem; margin: 1.6rem 0; }
.step h3 { margin-top: 0; }
.record { background: var(--panel); border: 1px dashed #c9c9c9; border-radius: 5px;
          padding: .6rem .9rem; margin: .7rem 0 .2rem; font-size: .84rem; }
.record strong { display: block; font-size: .7rem; text-transform: uppercase;
                 letter-spacing: .05em; color: var(--muted); margin-bottom: .35rem; }
.record li { margin: .15rem 0; }
.record .line { display: inline-block; border-bottom: 1px solid #b5b5b5; min-width: 8rem;
                margin-left: .4rem; }
.warnings { background: var(--warn-bg); border: 1px solid #e8d5a8; border-radius: 6px;
            padding: .8rem 1.1rem .9rem; }
.warnings h2 { border: 0; margin: 0 0 .4rem; color: var(--warn); font-size: .95rem; }
.warnings ul { margin: 0; padding-left: 1.1rem; }
.warnings li { font-size: .87rem; }
pre.plan { background: var(--panel); border: 1px solid var(--rule); border-radius: 6px;
           padding: .9rem 1.1rem; font-size: .84rem; white-space: pre-wrap;
           font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace; }
footer { margin-top: 3rem; padding-top: .8rem; border-top: 1px solid var(--rule);
         font-size: .78rem; color: var(--muted); }
input[type=number], input[type=text] {
  font: inherit; font-size: .85rem; width: 100%; min-width: 4.5rem; padding: .12rem .3rem;
  border: 1px solid #c9c9c9; border-radius: 3px; background: #fffdf6;
}
input:focus { outline: 2px solid #9db8d2; outline-offset: -1px; background: #fff; }
table.reaction .vol { font-variant-numeric: tabular-nums; }
table.reaction .vol.pending { color: var(--muted); }
table.reaction tr.total td { font-weight: 600; }
table.reaction tbody tr:last-child td { font-weight: 600; }
table.worksheet tbody tr:last-child td { font-weight: 400; }
.note-inline { color: var(--muted); font-size: .78rem; }
.hint { font-size: .78rem; color: var(--accent); margin: -.1rem 0 .5rem; }
.overflow { color: #b3261e; font-weight: 600; }
button.addrow {
  font: inherit; font-size: .8rem; padding: .2rem .6rem; margin: 0 0 .6rem;
  border: 1px solid #c9c9c9; border-radius: 4px; background: var(--panel); cursor: pointer;
}
button.addrow:hover { background: #ececea; }
@media print {
  body { padding: 0; font-size: 11pt; }
  .step, .figure, table { break-inside: avoid; }
  h2 { break-after: avoid; }
  .hint, button.addrow { display: none; }
  input { border: 0; border-bottom: 1px solid #999; border-radius: 0; background: none; }
}
"""

#: Recomputes every reaction mix from the concentrations typed into it. Deliberately tiny and
#: dependency-free: the report has to keep working from a file:// URL on a machine with no
#: network, years from now.
SCRIPT = """
(function () {
  function render(table) {
    var spec = JSON.parse(table.getAttribute('data-reaction'));
    var used = 0, known = true;
    spec.rows.forEach(function (row, i) {
      var cell = table.querySelector('.vol[data-row="' + i + '"]');
      var volume = row.fixed;
      if (volume === null) {
        var input = table.querySelector('.conc[data-row="' + i + '"]');
        var conc = input ? parseFloat(input.value) : NaN;
        volume = (conc > 0) ? row.ng / conc : null;
      }
      if (volume === null || isNaN(volume)) { known = false; cell.textContent = '____'; }
      else { used += volume; cell.textContent = volume.toFixed(2); }
      cell.classList.toggle('pending', volume === null || isNaN(volume));
    });
    var water = table.querySelector('.vol.water');
    var left = spec.total - used;
    water.textContent = known ? left.toFixed(2) : '____';
    // Overfilling is the failure this catches: too much dilute DNA and there is no room for
    // the buffer, which is a real mistake people make at 8pm.
    water.classList.toggle('overflow', known && left < 0);
    water.title = (known && left < 0)
      ? 'The DNA alone exceeds the reaction volume — concentrate it, or scale the reaction up.'
      : '';
  }
  document.querySelectorAll('table.reaction').forEach(function (table) {
    table.addEventListener('input', function () { render(table); });
    render(table);
  });
  document.querySelectorAll('button.addrow').forEach(function (button) {
    button.addEventListener('click', function () {
      var body = button.previousElementSibling.querySelector('tbody');
      var row = body.rows[body.rows.length - 1].cloneNode(true);
      row.querySelectorAll('input').forEach(function (input) { input.value = ''; });
      var first = row.cells[0];
      if (first && !first.querySelector('input')) {
        first.textContent = String(body.rows.length + 1);
      }
      body.appendChild(row);
    });
  });
})();
"""


def _escape(text) -> str:
    return _html.escape(str(text), quote=False)


def _cell(value) -> str:
    text = str(value)
    classes = ' class="seq"' if _SEQUENCE.match(text.strip()) else ""
    return f"<td{classes}>{_escape(text)}</td>"


def _table(table) -> str:
    from .bench import Reaction

    if isinstance(table, Reaction):
        return _reaction(table)

    head = "".join(f"<th>{_escape(column)}</th>" for column in table.columns)
    rows = []
    for row in table.rows:
        cells = []
        for index, cell in enumerate(row):
            if index in table.editable:
                cells.append(f'<td><input type="text" value="{_html.escape(str(cell))}"></td>')
            else:
                cells.append(_cell(cell))
        rows.append("<tr>" + "".join(cells) + "</tr>")
    note = f'<p class="note">{_escape(table.note)}</p>' if table.note else ""
    add = (
        '<button type="button" class="addrow">+ add row</button>' if table.expandable else ""
    )
    css = ' class="worksheet"' if table.editable else ""
    return (
        f"<table{css}><caption>{_escape(table.title)}</caption>"
        f"<thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>{add}{note}"
    )


def _reaction(reaction) -> str:
    """A live reaction mix: type in the stock concentrations, the volumes follow.

    Every number that can be computed ahead of time already is; the only inputs are the
    concentrations, which do not exist until the DNA has been through a Nanodrop.
    """
    rows = []
    for index, row in enumerate(reaction.rows):
        amount = "" if row.ng is None else f"{row.ng:.1f}"
        fmol = "" if row.fmol is None else f"{row.fmol:.1f}"
        if row.is_dna:
            value = "" if row.ng_per_ul is None else f"{row.ng_per_ul:g}"
            stock = (
                f'<input type="number" step="any" min="0" class="conc" '
                f'data-row="{index}" value="{value}" placeholder="ng/µL">'
            )
        else:
            stock = "—"
        volume = row.volume_ul
        rows.append(
            f"<tr><td>{_escape(row.name)}"
            + (f' <span class="note-inline">{_escape(row.note)}</span>' if row.note else "")
            + f"</td><td>{amount}</td><td>{fmol}</td><td>{stock}</td>"
            + f'<td class="vol" data-row="{index}">{_ul(volume)}</td></tr>'
        )

    water = reaction.water_ul
    rows.insert(
        0,
        f"<tr><td>{_escape(reaction.water_name)}</td><td></td><td></td><td>—</td>"
        f'<td class="vol water">{_ul(water)}</td></tr>',
    )
    rows.append(
        f'<tr class="total"><td>Total</td><td></td><td></td><td></td>'
        f"<td>{reaction.total_ul:.2f}</td></tr>"
    )
    note = f'<p class="note">{_escape(reaction.note)}</p>' if reaction.note else ""
    spec = _html.escape(json.dumps(reaction.spec()), quote=True)
    return (
        f'<table class="reaction" data-reaction="{spec}">'
        f"<caption>{_escape(reaction.title)}</caption><thead><tr>"
        "<th>Component</th><th>Amount (ng)</th><th>fmol</th>"
        "<th>Stock (ng/µL)</th><th>Volume (µL)</th>"
        f"</tr></thead><tbody>{''.join(rows)}</tbody></table>"
        '<p class="hint">Type the measured stock concentrations — the volumes recalculate.</p>'
        f"{note}"
    )


def _ul(value) -> str:
    from .bench import BLANK

    return BLANK if value is None else f"{value:.2f}"


def _figure(figure, wide=False) -> str:
    css = "figure wide" if wide else "figure"
    return (
        f'<figure class="{css}">{figure.svg}'
        f"<figcaption>{_escape(figure.caption)}</figcaption></figure>"
    )


def _record_block(fields) -> str:
    items = "".join(
        f'<li>{_escape(field)}<span class="line"></span></li>' for field in fields
    )
    return f'<div class="record"><strong>Record at the bench</strong><ul>{items}</ul></div>'


def _step(step) -> str:
    body = "".join(f"<p>{_escape(line)}</p>" for line in step.body)
    tables = "".join(_table(table) for table in step.tables)
    record = _record_block(step.record) if step.record else ""
    return (
        f'<section class="step"><h3>{_escape(step.title)}</h3>'
        f"{body}{tables}{record}</section>"
    )


def render_html(report) -> str:
    """The whole report as one HTML document."""
    meta = "".join(
        f"<div><span>{_escape(label)}</span>{_escape(value)}</div>"
        for label, value in report.meta
    )
    parts = [
        f"<h1>{_escape(report.title)}</h1>",
        f'<p class="subtitle">{_escape(report.subtitle)}</p>',
        f'<div class="meta">{meta}</div>',
    ]

    if report.aim:
        parts.append("<h2>Aim</h2>" + f"<p>{_escape(report.aim)}</p>")

    if report.strategy:
        items = "".join(f"<li>{_escape(line)}</li>" for line in report.strategy)
        parts.append(f"<h2>Cloning strategy</h2><ul>{items}</ul>")

    if report.plan_protocol:
        parts.append(
            "<h2>Verified plan</h2>"
            '<p class="note">Generated from the simulated and verified plan, reproduced here '
            "without alteration. The procedure below expands it; if the two ever disagree, "
            "this box is the one that was checked against the sequence.</p>"
            f'<pre class="plan">{_escape(report.plan_protocol)}</pre>'
        )

    if report.figures:
        # The builder emits the circular maps first, as a before/after pair, then the linear
        # zooms; the pair sits in two columns and the zooms run the full width.
        narrow, wide = report.figures[:2], report.figures[2:]
        figures = "".join(_figure(f) for f in narrow) + "".join(
            _figure(f, wide=True) for f in wide
        )
        parts.append(f'<h2>Plasmid maps</h2><div class="figures">{figures}</div>')

    if report.materials:
        parts.append("<h2>Materials</h2>" + "".join(_table(t) for t in report.materials))

    if report.steps:
        parts.append("<h2>Procedure</h2>" + "".join(_step(step) for step in report.steps))

    body = "\n".join(parts)
    return (
        "<!doctype html>\n"
        '<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>{_escape(report.title)}</title>"
        f"<style>{STYLESHEET}</style></head><body><main>\n{body}\n"
        "<footer>Generated by <strong>bbl</strong> from a simulation-verified plan. "
        "Volumes are calculated, not measured — check them against your own stocks before "
        "pipetting.</footer>"
        f"</main><script>{SCRIPT}</script></body></html>\n"
    )
