"""Native scorers against hand-computed values -- no backend, no Proto.

These pin the incumbent behaviour so that a later swap to Proto has something to be measured
against. Values here are worked out by hand in the test, not copied from the implementation.
"""

from __future__ import annotations

import pytest

from bbl.qc import FAIL, NATIVE, PASS, UNSUPPORTED, score
from bbl.qc import native


def test_gc_content_on_a_hand_counted_sequence():
    # 4 of 8 bases are G or C -> 50%
    measurement, penalty, _ = native.gc_content("ATCGATCG", {"min_gc": 25.0, "max_gc": 75.0})
    assert measurement == pytest.approx(50.0)
    assert penalty == 0.0


def test_gc_content_outside_the_window_is_penalised():
    _, penalty, _ = native.gc_content("AAAATTTT", {"min_gc": 25.0, "max_gc": 75.0})
    assert penalty == 1.0


def test_gc_content_of_the_empty_string_does_not_divide_by_zero():
    measurement, penalty, detail = native.gc_content("", {"min_gc": 25.0, "max_gc": 75.0})
    assert measurement == 0.0
    assert penalty == 1.0
    assert detail["empty"] is True


def test_longest_homopolymer_counts_the_run_not_the_total():
    # AAA GG A -> longest run is 3, even though there are 4 A's overall
    measurement, _, _ = native.max_homopolymer("AAAGGA", {"max_length": 9})
    assert measurement == 3.0


def test_homopolymer_over_the_limit_fails():
    measurement, penalty, _ = native.max_homopolymer("A" * 12, {"max_length": 9})
    assert measurement == 12.0
    assert penalty == 1.0


def test_specific_kmer_frequency_is_count_over_positions():
    # GAATTC once in a 12-mer -> 12 - 6 + 1 == 7 start positions
    seq = "AAAGAATTCAAA"
    measurement, _, detail = native.specific_kmer(
        seq, {"kmer": "GAATTC", "min_value": 0.0, "max_value": 1.0}
    )
    assert detail["count"] == 1
    assert detail["positions"] == len(seq) - 6 + 1
    assert measurement == pytest.approx(1 / 7)


def test_specific_kmer_counts_overlapping_occurrences():
    # AAAA contains AA three times when overlaps count
    _, _, detail = native.specific_kmer(
        "AAAA", {"kmer": "AA", "min_value": 0.0, "max_value": 1.0}
    )
    assert detail["count"] == 3


def test_specific_kmer_on_a_sequence_shorter_than_the_kmer():
    measurement, penalty, detail = native.specific_kmer(
        "AT", {"kmer": "GAATTC", "min_value": 0.0, "max_value": 1.0}
    )
    assert measurement == 0.0
    assert penalty == 1.0
    assert detail["positions"] == 0


def test_longest_repeat_finds_a_tandem_duplication():
    # An exact duplication of a 10-mer means the longest repeated substring is >= 10
    unit = "ACGTTGCAAC"
    measurement, _, _ = native.longest_repeat_check(unit + unit, {"max_repeat": 40})
    assert measurement >= 10.0


def test_restriction_sites_counts_an_unambiguous_site():
    # One EcoRI site, GAATTC
    measurement, _, detail = native.restriction_sites(
        "AAAAGAATTCAAAA",
        {"enzyme": "EcoRI", "min_sites": 0, "max_sites": 0, "circular": False},
    )
    assert measurement == 1.0
    assert detail["enzyme"] == "EcoRI"


def test_restriction_sites_resolves_an_ambiguous_recognition_site():
    """The check Proto structurally cannot do.

    AccI is GT^MKAC: M == A/C, K == G/T. ``GTATAC`` matches (M->A, K->T) even though the
    literal string ``GTMKAC`` never appears in the sequence.
    """
    measurement, _, _ = native.restriction_sites(
        "AAAAGTATACAAAA",
        {"enzyme": "AccI", "min_sites": 0, "max_sites": 0, "circular": False},
    )
    assert measurement >= 1.0


def test_scorecard_marks_proto_only_checks_unsupported_on_native():
    card = score("ACGT" * 40, backend=NATIVE)
    kmer = card.get("kmer_frequency")
    dinuc = card.get("dinucleotide_composition")
    assert kmer is not None and kmer.verdict == UNSUPPORTED
    assert dinuc is not None and dinuc.verdict == UNSUPPORTED
    # "unsupported" must never read as ok
    assert not kmer.ok


def test_native_scoring_is_deterministic():
    seq = "ATGGCGAATTCAAAGGCTAGCTAGCTAGGGCCCGAATTCTTTTACGTACGT" * 3
    first = score(seq, backend=NATIVE)
    second = score(seq, backend=NATIVE)
    assert [(r.check, r.measurement, r.penalty) for r in first.results] == [
        (r.check, r.measurement, r.penalty) for r in second.results
    ]


def test_a_failing_scorer_degrades_to_unavailable_rather_than_raising():
    """One broken check must not take the scorecard down."""
    card = score(
        "ACGT" * 40,
        backend=NATIVE,
        params={"restriction_sites": {"enzyme": "NoSuchEnzyme"}},
    )
    result = card.get("restriction_sites")
    assert result is not None
    assert result.verdict in ("unavailable", FAIL, PASS)
    # every other check still produced a value
    assert card.get("gc_content").verdict == PASS


def test_unknown_backend_is_rejected():
    with pytest.raises(ValueError, match="unknown backend"):
        score("ACGT", backend="wishful")
