"""Proto backend: the CPU gate, the graceful-absence path, and measurement agreement.

Every test that needs Proto skips cleanly without it, because ``bbl`` must remain installable
and testable with no Proto at all -- rejecting the dependency is a permitted outcome of this
bake-off.

Run the Proto side with (see docs/DECISIONS.md for why the venv is pinned):

    PYTHONPATH=<proto-language>:<proto-tools> $SCRATCH/bbl_proto_v3/bin/python -m pytest tests/test_qc_proto.py
"""

from __future__ import annotations

import pytest

from bbl.qc import (
    CHECKS,
    PROTO,
    PROTO_AVAILABLE,
    PROTO_IMPORT_ERROR,
    UNAVAILABLE,
    UNSUPPORTED,
    cpu_clean,
    describe_constraint,
    registry_report,
    score,
)

needs_proto = pytest.mark.skipif(not PROTO_AVAILABLE, reason=f"proto absent: {PROTO_IMPORT_ERROR}")

# A fragment with a known composition: one EcoRI site, no long runs.
FRAGMENT = "ATGGCGAATTCAAAGGCTAGCTAGCTAGGGCCCTTTTACGTACGTGGCCAATT" * 2


# --------------------------------------------------------------------------------------
# Behaviour when Proto is NOT installed. These must hold in the default environment.
# --------------------------------------------------------------------------------------


def test_bbl_qc_imports_without_proto():
    """Importing the QC layer must never require Proto."""
    import bbl.qc  # noqa: F401  -- the import itself is the assertion


@pytest.mark.skipif(PROTO_AVAILABLE, reason="proto is installed here")
def test_absent_proto_reports_unavailable_rather_than_raising():
    card = score(FRAGMENT, backend=PROTO)
    assert card.results, "a scorecard should still be produced"
    for result in card.results:
        assert result.verdict in (UNAVAILABLE, UNSUPPORTED)
        assert result.measurement is None
        # "we could not measure" must not be reported as a pass
        assert not result.ok


# --------------------------------------------------------------------------------------
# The CPU gate. This is the guard that keeps the layer demo-safe.
# --------------------------------------------------------------------------------------


@needs_proto
def test_every_registered_check_is_cpu_clean():
    """No check may reach a GPU constraint or one that shells out to a tool."""
    offenders = []
    for check in CHECKS:
        if check.proto_key is None:
            continue
        ok, reason = cpu_clean(check.proto_key)
        if not ok:
            offenders.append((check.name, check.proto_key, reason))
    assert not offenders, f"checks reaching non-CPU constraints: {offenders}"


@needs_proto
def test_the_gate_actually_rejects_the_constraints_we_wanted():
    """Regression on two constraints we would otherwise have used.

    ``longest-orf-length`` (frame integrity) and ``seq-motif`` (motif survival) both look
    ideal for this project and both fail the gate -- ORF prediction shells out to orfipy and
    seq-motif is PWM scanning via MEME FIMO, not the exact-substring test the name suggests.
    If a future Proto version makes either pure-Python, this test fails and we should adopt it.
    """
    ok, reason = cpu_clean("longest-orf-length")
    assert not ok and "orfipy" in reason

    ok, reason = cpu_clean("seq-motif")
    assert not ok and "fimo" in reason.lower()


@needs_proto
def test_gpu_constraints_are_rejected():
    ok, reason = cpu_clean("borzoi-track-activity")
    assert not ok and "gpu" in reason.lower()


@needs_proto
def test_unknown_constraint_key_is_not_clean_and_does_not_raise():
    ok, reason = cpu_clean("no-such-constraint-exists")
    assert not ok
    assert reason


@needs_proto
def test_registry_report_covers_every_proto_backed_check():
    rows = {row["check"] for row in registry_report() if "check" in row}
    expected = {check.name for check in CHECKS if check.proto_key}
    assert expected <= rows


@needs_proto
def test_describe_exposes_the_cost_declarations_we_gate_on():
    info = describe_constraint("gc-content")
    assert info["uses_gpu"] is False
    assert info["tools_called"] == []
    assert "dna" in info["supported_sequence_types"]
    assert info["cpu_clean"] is True


# --------------------------------------------------------------------------------------
# Measurement agreement and Proto's limits.
# --------------------------------------------------------------------------------------


@needs_proto
@pytest.mark.parametrize(
    "check,expected",
    [
        ("gc_content", pytest.approx(50.0)),
        ("sequence_length", 8.0),
        ("max_homopolymer", 1.0),
    ],
)
def test_proto_measurements_match_hand_computed_values(check, expected):
    card = score("ATCGATCG", backend=PROTO, checks=[check])
    assert card.get(check).measurement == expected


@needs_proto
def test_proto_cannot_express_an_ambiguous_restriction_site():
    """The prediction this bake-off was built to test.

    Proto's DNA alphabet is ``frozenset("ACGT")``, so AccI's ``GTMKAC`` is not representable.
    It fails loudly (ValueError) rather than silently miscounting, which is the good failure
    mode -- but it means the native Biopython route must be kept for restriction work.
    """
    card = score(
        "AAAAGTATACAAAA",
        backend=PROTO,
        checks=["specific_kmer"],
        params={"specific_kmer": {"kmer": "GTMKAC"}},
    )
    result = card.get("specific_kmer")
    assert result.verdict == UNSUPPORTED
    assert "invalid" in (result.note or "").lower()
    assert result.measurement is None


@needs_proto
def test_proto_scoring_is_deterministic():
    first = score(FRAGMENT, backend=PROTO)
    second = score(FRAGMENT, backend=PROTO)
    assert [(r.check, r.measurement, r.penalty) for r in first.results] == [
        (r.check, r.measurement, r.penalty) for r in second.results
    ]


@needs_proto
def test_proto_only_checks_produce_a_measurement():
    """The capability that would justify the dependency."""
    card = score(FRAGMENT, backend=PROTO, checks=["kmer_frequency", "dinucleotide_composition"])
    for name in ("kmer_frequency", "dinucleotide_composition"):
        result = card.get(name)
        assert result.measurement is not None, f"{name} produced no measurement"
        assert result.verdict not in (UNAVAILABLE, UNSUPPORTED)


@needs_proto
def test_a_bad_threshold_mapping_is_reported_not_raised():
    """A config our params cannot satisfy must degrade, not explode."""
    card = score(
        FRAGMENT,
        backend=PROTO,
        checks=["gc_content"],
        params={"gc_content": {"min_gc": 90.0, "max_gc": 10.0}},  # min > max
    )
    result = card.get("gc_content")
    assert result.verdict == UNAVAILABLE
    assert "config rejected" in (result.note or "")
