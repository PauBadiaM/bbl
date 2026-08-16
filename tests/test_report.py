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
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

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


def test_the_linear_zooms_are_the_full_width_ones(excision, parent):
    """Which figures span the page is the figure's own property, not its list index."""
    report = build_report(excision, parent=parent)
    assert [f.wide for f in report.figures] == [False, False, True, True]


def test_a_linear_panel_grows_to_fit_its_labels():
    """The regression the panels were unreadable for.

    DnaFeaturesViewer stacks colliding labels onto successive levels but only grows the figure
    to match when it owns the figure. ``linear_svg`` hands it an axes, so the height is ours to
    set -- and when it was a fixed 1.5 inches, a crowded window crushed every level into it.
    Built from a synthetic record rather than the inventory so it runs anywhere.
    """
    from bbl.report.maps import _LINEAR_MIN_HEIGHT, linear_svg

    def record(count):
        seq = Seq("ATGC" * 250)
        features = [
            SeqFeature(
                FeatureLocation(60 * i + 10, 60 * i + 40),
                type="CDS",
                qualifiers={"label": [f"a long enough feature name {i}"]},
            )
            for i in range(count)
        ]
        return SeqRecord(seq, name="crowded", features=features)

    def height(svg):
        box = re.search(r'viewBox="[\d.]+ [\d.]+ ([\d.]+) ([\d.]+)"', svg)
        return float(box.group(2)) / float(box.group(1))  # height as a fraction of width

    sparse, crowded = height(linear_svg(record(2), (0, 1000))), height(
        linear_svg(record(14), (0, 1000))
    )
    assert crowded > sparse * 1.5, "a crowded window must get more vertical room, not less"
    # ...and the sparse one is not left as a sliver.
    assert sparse > 0
    assert _LINEAR_MIN_HEIGHT >= 1.5


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
# the strategy says where both pieces of DNA come from
# --------------------------------------------------------------------------- #


def _strategy(report) -> str:
    return " ".join(report.strategy)


def test_the_strategy_names_the_plasmid_the_insert_is_cut_out_of(backbone, donor):
    """The backbone has always been named. The insert is half the reaction and was not."""
    plan = plan_insertion(backbone, donor, insert_features=["Lambda BoxB x8"], at=BOXB_SITE)
    text = _strategy(build_report(plan, parent=backbone, donor=donor, figures=False))
    assert f"from {donor.name}" in text
    assert f"{len(donor):,} bp" in text
    assert " + ".join(plan.insert.enzymes) in text


def test_the_strategy_names_the_plasmid_the_insert_is_amplified_off(backbone, donor):
    plan = plan_insertion(
        backbone, donor, insert_features=["Lambda BoxB x8"], at=BOXB_SITE, method="gibson"
    )
    text = _strategy(build_report(plan, parent=backbone, donor=donor, figures=False))
    assert f"PCR-amplified from {donor.name}" in text
    assert all(p.name in text for p in plan.insert.primers)


def test_the_strategy_uses_the_full_construct_label_for_the_donor(backbone, donor):
    """A 16-character LOCUS is an accession; the strategy section is where the label belongs."""
    label = "pCLM1_pcDNA3.1_miniCMV-mCherry-LambdaBoxBx8"
    plan = plan_insertion(backbone, donor, insert_features=["Lambda BoxB x8"], at=BOXB_SITE)
    report = build_report(
        plan, parent=backbone, donor=donor, donor_name=label, figures=False
    )
    assert f"from {label}" in _strategy(report)


def test_the_donor_survives_on_the_plan_when_the_caller_forgets_it(backbone, donor):
    """``build_report(donor=...)`` is optional; the provenance should not depend on it."""
    plan = plan_insertion(backbone, donor, insert_features=["Lambda BoxB x8"], at=BOXB_SITE)
    assert plan.insert.donor_name == donor.name
    assert f"from {donor.name}" in _strategy(build_report(plan, parent=backbone, figures=False))


def test_a_synthetic_insert_says_so_rather_than_naming_a_donor(backbone):
    plan = plan_insertion(backbone, "ATG" + "GCTAGCTAGG" * 30, at=BOXB_SITE, method="gibson")
    text = _strategy(build_report(plan, parent=backbone, figures=False))
    assert "synthetic fragment" in text
    assert "from" not in text.split("Insert:")[1].split(".")[0]


def test_the_strategy_says_when_the_backbone_itself_has_to_be_amplified(backbone, donor):
    """Otherwise the first mention of a vector PCR is in step 1, after the reader has planned."""
    plan = plan_insertion(
        backbone, donor, insert_features=["Lambda BoxB x8"], at=BOXB_SITE,
        method="gibson", vector_prep="pcr",
    )
    text = _strategy(build_report(plan, parent=backbone, donor=donor, figures=False))
    assert "inverse PCR" in text
    assert all(p.name in text for p in plan.vector.primers)


# --------------------------------------------------------------------------- #
# diluting the template
# --------------------------------------------------------------------------- #


def _dilutions(report):
    from bbl.report.bench import Dilution

    return [t for s in report.steps for t in s.tables if isinstance(t, Dilution)]


def test_every_pcr_step_carries_a_dilution_box(no_sites_phl391):
    template = load_plasmid(no_sites_phl391)
    plan = design_deletion_primers(template, ["NFKBRE"])
    report = build_report(plan, parent=template, figures=False)
    boxes = _dilutions(report)
    assert len(boxes) == 1
    assert boxes[0].template == template.name


def test_a_gibson_route_dilutes_both_templates(backbone, donor):
    """Two PCRs, two different plasmids on the bench, two boxes."""
    plan = plan_insertion(
        backbone, donor, insert_features=["Lambda BoxB x8"], at=BOXB_SITE,
        method="gibson", vector_prep="pcr",
    )
    report = build_report(plan, parent=backbone, donor=donor, figures=False)
    assert {b.template for b in _dilutions(report)} == {backbone.name, donor.name}


def test_the_dilution_box_uses_a_supplied_concentration(no_sites_phl391):
    template = load_plasmid(no_sites_phl391)
    plan = design_deletion_primers(template, ["NFKBRE"])
    report = build_report(
        plan, parent=template, figures=False, concentrations={template.name: 1364}
    )
    box = _dilutions(report)[0]
    assert box.ng_per_ul == 1364
    assert box.plan()["status"] == "ok"


def test_a_digest_only_route_has_no_dilution_box(excision, parent):
    """Nothing is amplified, so there is no template to dilute and no box to ignore."""
    assert _dilutions(build_report(excision, parent=parent, figures=False)) == []


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
    # The reaction calculators need script and the header carries a logo, but nothing may be
    # *fetched*: the report has to work from a file:// URL on a machine with no network. So the
    # property is not "no images", it is that every src is inline -- a ``data:`` URI is part of
    # the document, an http:// or a relative path is a dependency on something outside it.
    assert all(src.startswith("data:") for src in re.findall(r'src="([^"]*)"', html))
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


def _node_dilute(html, stocks):
    """Run the page's own ``dilutionPlan`` over its embedded specs, at the given stocks.

    Same pin as :func:`_node_recompute`, for the second calculator: the dilution arithmetic
    exists twice, once in Python and once in the page, and two copies of a formula drift.
    The script is lifted out of the document rather than re-typed here, so the thing under
    test is the thing that ships.
    """
    import json
    import shutil
    import subprocess
    import textwrap

    node = shutil.which("node")
    if node is None:
        pytest.skip("node not available")

    harness = textwrap.dedent(
        """
        const html = require('fs').readFileSync(process.argv[2], 'utf8');
        const stocks = JSON.parse(process.argv[3]);
        const script = html.match(/<script>([\\s\\S]*)<\\/script>/)[1];
        // The page's IIFE runs against a DOM; pull out just the two pure functions.
        const body = script.slice(script.indexOf('function dilutionPlan'),
                                  script.indexOf("document.querySelectorAll('table.dilution')"));
        const make = new Function(body + '; return {dilutionPlan, dilutionRows};')();
        const specs = [...html.matchAll(/data-dilution="([^"]+)"/g)].map(m =>
          JSON.parse(m[1].replace(/&quot;/g,'"').replace(/&amp;/g,'&')
                         .replace(/&lt;/g,'<').replace(/&gt;/g,'>')));
        console.log(JSON.stringify(specs.map((spec, i) => {
          const stock = (stocks[i] !== undefined) ? stocks[i] : spec.conc;
          const plan = make.dilutionPlan(stock, spec);
          return {plan: plan, rows: make.dilutionRows(plan, spec.target)};
        })));
        """
    )
    with tempfile.TemporaryDirectory() as work:
        page = Path(work) / "page.html"
        page.write_text(html)
        runner = Path(work) / "dilute.js"
        runner.write_text(harness)
        result = subprocess.run(
            [node, str(runner), str(page), json.dumps(stocks)],
            capture_output=True,
            text=True,
            check=True,
        )
    return json.loads(result.stdout)


@pytest.mark.parametrize("stock", [None, 0.4, 1, 20, 250, 1364, 9999])
def test_the_pages_dilution_matches_the_python_one(no_sites_phl391, stock):
    from bbl.report.html import _dilution_rows

    template = load_plasmid(no_sites_phl391)
    plan = design_deletion_primers(template, ["NFKBRE"])
    concentrations = {} if stock is None else {template.name: stock}
    report = build_report(
        plan, parent=template, figures=False, concentrations=concentrations
    )
    box = _dilutions(report)[0]
    computed = _node_dilute(report.to_html(), {} if stock is None else {0: stock})[0]

    assert computed["plan"]["status"] == box.plan()["status"]
    for js, py in zip(computed["plan"]["steps"], box.plan()["steps"]):
        assert js["source"] == py["source"]
        for key in ("take_ul", "water_ul", "final_ul", "gives"):
            assert js[key] == pytest.approx(py[key])
    # ...and the markup they each produce is the same markup, not just the same numbers.
    assert computed["rows"] == _dilution_rows(box.plan(), box.target_ng_per_ul)
