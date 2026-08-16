"""The bake-off itself: both backends on real sequences from the inventory.

The assertion that matters is **zero disagreements** -- two backends computing the same
physical quantity must produce the same number. Coverage asymmetry (checks only one backend
can do) is reported, not asserted, because that asymmetry is the finding.
"""

from __future__ import annotations

import pytest

from bbl import load_plasmid
from bbl.qc import (
    NATIVE,
    PROTO,
    PROTO_AVAILABLE,
    PROTO_IMPORT_ERROR,
    compare,
    score,
)

needs_proto = pytest.mark.skipif(not PROTO_AVAILABLE, reason=f"proto absent: {PROTO_IMPORT_ERROR}")

# The two primers from the around-the-horn design on pHL391 (tests/test_pcr.py).
DEL_F = "TAGGCGTGTACGGTGGGAGGCCTATATAAG"
DEL_R = "TCGGTCAAGCCTTGCCTTGTTGTAGC"


def test_compare_runs_without_proto_and_reports_coverage():
    """Without Proto the report is still well-formed; everything is native-only."""
    report = compare("ACGT" * 50, label="uniform")
    assert report.rows
    if not PROTO_AVAILABLE:
        assert not report.coverage("both")
        assert report.disagreements == []
    assert isinstance(report.summary(), str)


@needs_proto
@pytest.mark.parametrize("primer,name", [(DEL_F, "del_F"), (DEL_R, "del_R")])
def test_backends_agree_on_the_real_deletion_primers(primer, name):
    report = compare(primer, label=name)
    assert report.disagreements == [], report.summary()
    assert report.coverage("both"), "no check was comparable, so nothing was tested"


@needs_proto
def test_backends_agree_on_a_real_plasmid_region(phl391):
    """A genuine ~300 bp slice of pHL391, not a synthetic string."""
    record = load_plasmid(phl391)
    region = str(record.seq)[200:500]
    report = compare(region, label="pHL391[200:500]")
    assert report.disagreements == [], report.summary()


@needs_proto
def test_backends_agree_across_the_whole_inventory(plasmid_dir):
    """Sweep every plasmid; a single disagreement anywhere is a bug in one backend."""
    files = sorted(plasmid_dir.glob("*.dna"))
    if not files:
        pytest.skip("no .dna files in the inventory")
    problems = []
    for path in files[:12]:  # a representative slice keeps the suite quick
        try:
            record = load_plasmid(path)
        except Exception:  # noqa: BLE001 -- unreadable files are not this test's concern
            continue
        region = str(record.seq)[:600]
        report = compare(region, label=path.name)
        if report.disagreements:
            problems.append(report.summary())
    assert not problems, "\n\n".join(problems)


@needs_proto
def test_the_two_backends_disagree_on_penalty_shape_even_when_measurements_match():
    """Documents an expected, informative difference rather than a defect.

    Native penalties are step functions (``complexity.py`` adds a fixed weight per reason);
    Proto's homopolymer penalty is log-scaled. Both flag the same run length, but only Proto
    ranks a run of 12 as worse than a run of 10 -- which is the argument for adopting it in
    any context where candidates get *ranked* rather than merely accepted or rejected.
    """
    # Flanks deliberately free of A so the run is exactly 12 and not 12 + a neighbour.
    seq = "GCGC" * 10 + "A" * 12 + "GCGC" * 10
    params = {"max_homopolymer": {"max_length": 9}}
    native_card = score(seq, backend=NATIVE, checks=["max_homopolymer"], params=params)
    proto_card = score(seq, backend=PROTO, checks=["max_homopolymer"], params=params)

    native_result = native_card.get("max_homopolymer")
    proto_result = proto_card.get("max_homopolymer")

    # same measurement
    assert native_result.measurement == proto_result.measurement == 12.0
    # different judgement: native saturates, Proto grades
    assert native_result.penalty == 1.0
    assert 0.0 < proto_result.penalty < 1.0


@needs_proto
def test_ambiguous_site_shows_up_as_coverage_not_disagreement():
    """Proto refusing GTMKAC must not be recorded as the backends disagreeing."""
    report = compare(
        "AAAAGTATACAAAA",
        label="accI-site",
        params={"specific_kmer": {"kmer": "GTMKAC"}},
    )
    row = next(r for r in report.rows if r.check == "specific_kmer")
    assert not row.comparable
    assert row.coverage == "native-only"
    assert row not in report.disagreements
