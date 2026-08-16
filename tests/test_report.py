"""Reports: figure editing, route coverage, and the properties that make one trustworthy.

The interesting assertions are not "it produced HTML". They are that the report never
contradicts the plan it was built from, never invents a number nobody measured, and never
silently drops the plan's warnings -- a protocol that reads cleanly because the caveats fell
off is worse than no protocol.
"""

import re
import tempfile
from pathlib import Path

import pytest

from bbl import (
    design_deletion_primers,
    excise_features,
    load_plasmid,
    plan_insertion,
)
from bbl.plasmid_io import feature_label, feature_span
from bbl.report import build_report, write_report
from bbl.report.bench import BLANK
from bbl.report.maps import prune_features

pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")

BOXB_SITE = (1491, 1512)


@pytest.fixture(scope="module")
def parent(plasmid_dir):
    return load_plasmid(plasmid_dir / "pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8.dna")


@pytest.fixture(scope="module")
def excision(parent):
    return excise_features(parent, ["NFKBRE"])


@pytest.fixture(scope="module")
def backbone(plasmid_dir):
    return load_plasmid(plasmid_dir / "pCLM3_pcDNA3.1_CMV-mCherry.dna")


@pytest.fixture(scope="module")
def donor(plasmid_dir):
    return load_plasmid(plasmid_dir / "pCLM1_pcDNA3.1_miniCMV-mCherry-LambdaBoxBx8.dna")


def _mix(report, title_contains):
    """The reaction mix whose title contains ``title_contains``."""
    for step in report.steps:
        for table in step.tables:
            if title_contains in table.title:
                return table
    raise AssertionError(f"no table matching {title_contains!r}")


def _text(html: str) -> str:
    """The document with markup and figures stripped, for prose assertions."""
    stripped = re.sub(r"<svg.*?</svg>", " ", html, flags=re.S)
    stripped = re.sub(r"<(style|script).*?</\1>", " ", stripped, flags=re.S)
    return re.sub(r"<[^>]+>", " ", stripped)


# --------------------------------------------------------------------------- #
# what a figure shows
# --------------------------------------------------------------------------- #


def test_pruning_drops_an_umbrella_but_keeps_a_repeat_array(parent):
    """Two kinds of nesting, opposite answers.

    'Insert Sequence' wraps a whole cassette whose parts have different names, so the parts
    are what a reader wants. 'Lambda BoxB x8' wraps eight identically-named monomers, so the
    wrapper is what a reader wants.
    """
    labels = {feature_label(f) for f in prune_features(parent.features, min_bp=60)}
    assert "Insert Sequence" not in labels
    assert "Lambda BoxB x8" in labels
    assert "BoxB RNA aptamer" not in labels
    assert {"mCherry", "AmpR", "ori"} <= labels


def test_pruning_keeps_a_protected_target_and_folds_away_its_children(parent):
    """NFKBRE is 54 bp and contains five sub-sites; a deletion report must still name it."""
    target = [feature_span(f) for f in parent.features if feature_label(f) == "NFKBRE"]
    labels = [feature_label(f) for f in prune_features(parent.features, min_bp=60, protect=target)]
    assert "NFKBRE" in labels
    assert not any(label.startswith("NFKB_") for label in labels)


def test_pruning_never_keeps_a_primer_binding(parent):
    kept = prune_features(parent.features)
    assert all(f.type != "primer_bind" for f in kept)


def test_figures_are_inline_svg_with_unique_internal_ids(excision, parent, tmp_path):
    """Four figures share one document, so a clip path in one must not capture another."""
    html = write_report(excision, tmp_path / "r.html", parent=parent).read_text()
    figures = re.findall(r"<svg.*?</svg>", html, flags=re.S)
    assert len(figures) == 4
    # Within one figure an id is referenced many times over; what must not happen is the
    # same id appearing in two figures, where the second would resolve against the first.
    per_figure = [
        {
            identifier
            for identifier in re.findall(r"url\(#([^)]+)\)", figure)
            if not identifier.startswith("DejaVu")
        }
        for figure in figures
    ]
    pooled = [identifier for ids in per_figure for identifier in ids]
    assert len(pooled) == len(set(pooled))
    assert "<?xml" not in html  # a prolog inside <body> is not legal


def test_a_report_can_be_built_without_figures(excision, parent):
    """``figures=False`` is the only way to skip them now -- the procedure is unaffected."""
    report = build_report(excision, parent=parent, figures=False)
    assert report.figures == []
    assert report.steps


# --------------------------------------------------------------------------- #
# the report agrees with the plan
# --------------------------------------------------------------------------- #


def test_the_plans_own_protocol_is_reproduced_verbatim(excision, parent):
    report = build_report(excision, parent=parent, figures=False)
    assert report.plan_protocol == excision.protocol
    assert excision.protocol in _text(report.to_html()).replace("&#x27;", "'")


def test_every_plan_warning_reaches_the_report(excision, parent):
    report = build_report(excision, parent=parent, figures=False)
    for warning in excision.warnings:
        assert warning in report.warnings


def test_lengths_and_enzymes_are_not_paraphrased(excision, parent):
    text = _text(build_report(excision, parent=parent, figures=False).to_html())
    assert f"{len(excision.product):,} bp" in text
    for enzyme in excision.enzyme_pair:
        assert enzyme in text
    assert excision.junction_seq in text


def test_a_deletion_by_religation_is_not_dephosphorylated(excision, parent):
    """The vector closing on itself is the product here, so suppressing it would be wrong."""
    report = build_report(excision, parent=parent, figures=False)
    digest = _mix(report, "Digestion")
    assert not any("FastAP" in row.name for row in digest.rows)
    assert any("Do not dephosphorylate" in line for step in report.steps for line in step.body)


def test_a_nondirectional_vector_is_dephosphorylated(backbone, donor, monkeypatch):
    """The opposite case: compatible vector ends can self-close, so they are treated."""
    plan = plan_insertion(backbone, donor, insert_features=["Lambda BoxB x8"], at=BOXB_SITE)
    monkeypatch.setattr(plan, "directional", False)
    report = build_report(plan, parent=backbone, donor=donor, figures=False)
    assert any("FastAP" in row.name for row in _mix(report, "Digestion").rows)


# --------------------------------------------------------------------------- #
# nothing is invented
# --------------------------------------------------------------------------- #


def test_unmeasured_concentrations_are_left_blank(excision, parent):
    text = _text(build_report(excision, parent=parent, figures=False).to_html())
    assert BLANK in text


def test_a_supplied_concentration_is_used(excision, parent):
    report = build_report(
        excision, parent=parent, figures=False, concentrations={parent.name: 250}
    )
    dna = next(t for t in report.materials if t.title == "DNA")
    assert dna.rows[0][3] == "250"


def test_enzyme_stock_comes_from_the_lab_config(excision, parent):
    report = build_report(
        excision, parent=parent, figures=False, config={"enzyme_stock": ["MfeI"]}
    )
    enzymes = next(t for t in report.materials if t.title == "Enzymes")
    assert dict(enzymes.rows) == {"MfeI": "yes", "EcoRI": "order"}


def test_no_stock_list_says_so_rather_than_assuming_yes(excision, parent):
    report = build_report(excision, parent=parent, figures=False)
    enzymes = next(t for t in report.materials if t.title == "Enzymes")
    assert all(row[1] == BLANK for row in enzymes.rows)
    assert "check the freezer" in enzymes.note


# --------------------------------------------------------------------------- #
# each route produces a procedure that fits it
# --------------------------------------------------------------------------- #


def test_restriction_deletion_reads_as_a_digest(excision, parent):
    report = build_report(excision, parent=parent, figures=False)
    titles = " ".join(step.title for step in report.steps)
    assert "Digest" in titles and "Religate" in titles
    assert "PCR" not in titles


def test_pcr_deletion_orders_oligos_and_kills_the_template(plasmid_dir, no_sites_phl391):
    template = load_plasmid(no_sites_phl391)
    plan = design_deletion_primers(template, ["NFKBRE"])
    report = build_report(plan, parent=template, figures=False)
    oligos = next(t for t in report.materials if "Oligos" in t.title)
    assert [row[0] for row in oligos.rows] == [plan.forward.name, plan.reverse.name]
    assert plan.forward.sequence in oligos.rows[0][1]
    assert any("DpnI" in line for step in report.steps for line in step.body)


def test_pcr_report_reconciles_its_extension_time_with_the_plan(no_sites_phl391):
    """The plan quotes a generic 30 s/kb; a repeat-bearing amplicon gets 45. Say which."""
    template = load_plasmid(no_sites_phl391)
    plan = design_deletion_primers(template, ["NFKBRE"])
    report = build_report(plan, parent=template, figures=False)
    assert any("verified plan above" in line for step in report.steps for line in step.body)


def test_insertion_digests_the_donor_plasmid_not_the_band(backbone, donor):
    """You cut the whole donor; the 300 bp fragment only exists afterwards."""
    plan = plan_insertion(backbone, donor, insert_features=["Lambda BoxB x8"], at=BOXB_SITE)
    report = build_report(plan, parent=backbone, donor=donor, figures=False)
    donor_digest = _mix(report, f"Digestion — {donor.name}")
    assert any(f"{len(donor):,} bp" in row.name for row in donor_digest.rows), (
        "the donor digest should be costed against the whole donor plasmid"
    )


def test_gibson_insertion_uses_an_assembly_table(backbone, donor):
    plan = plan_insertion(
        backbone, donor, insert_features=["Lambda BoxB x8"], at=BOXB_SITE, method="gibson"
    )
    report = build_report(plan, parent=backbone, donor=donor, figures=False)
    titles = [t.title for step in report.steps for t in step.tables]
    assert any("Gibson assembly" in title for title in titles)
    assert not any("Ligation" in title for title in titles)


def test_a_synthetic_insert_is_listed_as_something_to_order(backbone):
    plan = plan_insertion(backbone, "ATG" + "GCTAGCTAGG" * 30, at=BOXB_SITE, method="gibson")
    report = build_report(plan, parent=backbone, figures=False)
    order = next(t for t in report.materials if "Synthetic" in t.title)
    assert plan.insert.order_sequence in order.rows[0][2]


# --------------------------------------------------------------------------- #
# repeat arrays change what you do at the bench
# --------------------------------------------------------------------------- #


def test_a_repeat_bearing_product_is_grown_at_thirty_degrees(excision, parent):
    text = _text(build_report(excision, parent=parent, figures=False).to_html())
    assert "30 °C" in text
    assert "recombine" in text


def test_a_product_without_repeats_is_grown_at_thirtyseven(plasmid_dir):
    """The empty MCS vector carries no array, so the caution should not appear."""
    vector = load_plasmid(plasmid_dir / "pHL162_pcDNA3.1_MCS.dna")
    plan = plan_insertion(vector, "ATGGCTAGCTAGGACCTGATCAAGGGTCCATG" * 10, at=1000,
                          method="gibson")
    text = _text(build_report(plan, parent=vector, figures=False).to_html())
    assert "37 °C" in text
    assert "30 °C" not in text


# --------------------------------------------------------------------------- #
# the document itself
# --------------------------------------------------------------------------- #


def test_the_report_is_one_self_contained_file(excision, parent, tmp_path):
    path = write_report(excision, tmp_path / "report.html", parent=parent)
    html = path.read_text()
    assert list(tmp_path.iterdir()) == [path]  # no sidecar images
    assert html.startswith("<!doctype html>")
    # The reaction calculators need script, but nothing may be *fetched*: the report has to
    # work from a file:// URL on a machine with no network.
    assert "<img" not in html and "src=" not in html
    assert "<link" not in html and "@import" not in html
    assert html.count("<script>") == 1


def test_product_and_parent_never_share_a_name(backbone, donor):
    """The restriction route keeps the vector's name; a before/after document cannot."""
    plan = plan_insertion(backbone, donor, insert_features=["Lambda BoxB x8"], at=BOXB_SITE)
    report = build_report(plan, parent=backbone, donor=donor, figures=False)
    product = dict(report.meta)["Product"]
    assert not product.startswith(f"{backbone.name} —")


def test_a_verified_comparison_is_stated_in_the_header(excision, parent):
    report = build_report(
        excision, parent=parent, figures=False, verified_against=("pCLM1", True)
    )
    assert dict(report.meta)["Verified"] == "identical to pCLM1"


def test_html_escapes_a_hostile_aim(excision, parent):
    report = build_report(
        excision, parent=parent, figures=False, aim="<script>alert('x')</script>"
    )
    html = report.to_html()
    assert "<script>alert" not in html
    assert "&lt;script&gt;" in html


# --------------------------------------------------------------------------- #
# the in-page calculator agrees with the Python one
# --------------------------------------------------------------------------- #


def _node_recompute(html, concentrations_by_row=None):
    """Run the report's own JavaScript arithmetic over its own embedded spec, under node.

    The reaction tables are calculated twice -- once in Python when the report is written,
    once in the browser whenever somebody types a concentration. Two implementations of one
    formula drift; this pins them together.
    """
    import json
    import shutil
    import subprocess
    import textwrap

    node = shutil.which("node")
    if node is None:
        pytest.skip("node not available")

    script = textwrap.dedent(
        """
        const html = require('fs').readFileSync(process.argv[2], 'utf8');
        const typed = JSON.parse(process.argv[3]);
        const specs = [...html.matchAll(/data-reaction="([^"]+)"/g)].map(m =>
          JSON.parse(m[1].replace(/&quot;/g,'"').replace(/&amp;/g,'&')
                         .replace(/&lt;/g,'<').replace(/&gt;/g,'>')));
        const out = specs.map(spec => {
          let used = 0, known = true;
          const vols = spec.rows.map((row, i) => {
            let v = row.fixed;
            if (v === null) {
              const c = (typed[i] !== undefined) ? typed[i] : row.conc;
              v = (c > 0) ? row.ng / c : null;
            }
            if (v === null || isNaN(v)) { known = false; return null; }
            used += v; return v;
          });
          return {volumes: vols, water: known ? spec.total - used : null, total: spec.total};
        });
        console.log(JSON.stringify(out));
        """
    )
    with tempfile.TemporaryDirectory() as work:
        page = Path(work) / "page.html"
        page.write_text(html)
        runner = Path(work) / "calc.js"
        runner.write_text(script)
        result = subprocess.run(
            [node, str(runner), str(page), json.dumps(concentrations_by_row or {})],
            capture_output=True,
            text=True,
            check=True,
        )
    return json.loads(result.stdout)


def test_the_pages_calculator_matches_the_python_one(excision, parent, tmp_path):
    stock = {parent.name: 1364}
    report = build_report(excision, parent=parent, figures=False, concentrations=stock)
    mix = _mix(report, "Digestion")
    computed = _node_recompute(report.to_html())[0]

    assert computed["volumes"] == [
        pytest.approx(row.volume_ul) for row in mix.rows
    ]
    assert computed["water"] == pytest.approx(mix.water_ul)
    assert sum(computed["volumes"]) + computed["water"] == pytest.approx(mix.total_ul)


def test_typing_a_concentration_fills_the_volumes_in_the_page(excision, parent, tmp_path):
    """Written with no concentration, the page still computes once one is typed in."""
    report = build_report(excision, parent=parent, figures=False)
    assert _mix(report, "Digestion").water_ul is None  # blank as written

    blank = _node_recompute(report.to_html())[0]
    assert blank["water"] is None

    # row 1 is the DNA; 1 ug at 1364 ng/uL is 0.73 uL
    filled = _node_recompute(report.to_html(), {1: 1364})[0]
    assert filled["volumes"][1] == pytest.approx(0.733, abs=0.001)
    assert filled["water"] == pytest.approx(15.267, abs=0.001)


def test_the_page_flags_a_reaction_the_dna_overflows(excision, parent):
    """A very dilute stock needs more volume than the reaction has. Water goes negative."""
    report = build_report(excision, parent=parent, figures=False)
    overflowed = _node_recompute(report.to_html(), {1: 20})[0]  # 1 ug at 20 ng/uL is 50 uL
    assert overflowed["water"] < 0
