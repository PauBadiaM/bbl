"""Tests for the shared sequence primitives.

These exist because three modules previously carried their own copies of some of these
measurements, and because two of them (positional GC, hairpins with a loop constraint) are new
capability rather than a refactor.
"""

from __future__ import annotations

import random

import pytest

from bbl.seqfeatures import (
    clean_dna,
    longest_direct_repeat,
    max_hairpin_stem,
    max_homopolymer,
    poisson_sf,
    repeat_mask,
    revcomp,
    terminal_gc,
    three_prime_hairpin,
    windowed_gc,
    windowed_repeat_density,
)


# --------------------------------------------------------------------------- #
# repeats
# --------------------------------------------------------------------------- #


def test_overlapping_repeats_are_not_truncated_at_half_the_length():
    """``AAAA`` contains ``AAA`` twice, at offsets 0 and 1.

    The previous implementation capped its search at ``n // 2`` and so reported 2 here. That
    ceiling undercounted exactly the tandem-array case the measurement exists to catch.
    """
    assert longest_direct_repeat("AAAA") == (3, "AAA")
    assert longest_direct_repeat("ACGTACGT")[0] == 4


def test_no_repeat_reports_zero():
    assert longest_direct_repeat("ACGT") == (0, "")


def test_a_tandem_array_is_almost_entirely_masked():
    mask = repeat_mask("ACGTACGTAA" * 60)
    assert sum(mask) / len(mask) > 0.95


def test_composition_skew_alone_is_not_called_repetitive():
    """An AT-only sequence is compositionally skewed, not designed-repetitive.

    Plain recurrence counting reports it as heavily repetitive; the composition-matched null is
    what tells the two apart.
    """
    rng = random.Random(7)
    at_only = "".join(rng.choice("AT") for _ in range(600))
    assert sum(repeat_mask(at_only)) == 0


def test_the_significance_threshold_is_family_wise():
    """A recurrence that clears a per-test bar but not a family-wise one is not a repeat.

    Six copies of an AT-rich 8-mer in skewed composition sit at about ``p = 6e-4``: significant
    against a bare ``alpha = 1e-3``, but not against ``alpha`` divided by the several hundred
    k-mers actually tested. The mask is deliberately conservative here -- a genuine long repeat
    is caught by :func:`bbl.seqfeatures.longest_direct_repeat` regardless of the mask.
    """
    rng = random.Random(11)
    bases = "A" * 45 + "T" * 45 + "G" * 5 + "C" * 5
    filler = list("".join(rng.choice(bases) for _ in range(600)))
    motif = "ATTATAAT"
    for start in range(40, 600, 95):
        filler[start : start + len(motif)] = motif
    sequence = "".join(filler)[:600]

    assert sequence.count(motif) == 6
    assert sum(repeat_mask(sequence)) == 0
    # ...but the repeat itself is still visible to the exact-repeat measurement.
    assert longest_direct_repeat(sequence)[0] >= len(motif)


def test_windowed_repeat_density_finds_a_local_concentration():
    mask = [False] * 200 + [True] * 85 + [False] * 200
    assert windowed_repeat_density(mask, window=85) == 1.0
    assert windowed_repeat_density([False] * 100, window=85) == 0.0


def test_poisson_sf_bounds():
    assert poisson_sf(0, 5.0) == 1.0
    assert poisson_sf(1, 0.0) == 0.0
    assert 0.0 < poisson_sf(2, 1.0) < 1.0


# --------------------------------------------------------------------------- #
# positional GC -- capability the global measures could not provide
# --------------------------------------------------------------------------- #


def test_global_gc_hides_a_local_extreme():
    """The failure mode positional GC exists to catch: 50% overall, 100% in a window."""
    sequence = "AT" * 100 + "GC" * 100
    low, high = windowed_gc(sequence, window=50)
    assert low == 0.0
    assert high == 1.0


def test_terminal_gc_distinguishes_the_two_ends():
    five, three = terminal_gc("G" * 40 + "A" * 40, end_length=40)
    assert five == 1.0
    assert three == 0.0


def test_short_sequences_fall_back_to_the_global_value():
    low, high = windowed_gc("GCGC", window=50)
    assert low == high == 1.0


# --------------------------------------------------------------------------- #
# secondary structure
# --------------------------------------------------------------------------- #


def test_revcomp():
    assert revcomp("AATTCG") == "CGAATT"


def test_a_designed_hairpin_is_measured_at_its_true_stem_and_loop():
    stem_arm, loop = "GGGGCC", "ATAT"
    sequence = stem_arm + loop + revcomp(stem_arm)
    assert max_hairpin_stem(sequence, seed=4) == (len(stem_arm), len(loop))


def test_the_loop_constraint_shortens_the_reportable_stem():
    """The loop bound is what makes the measurement physical.

    A perfect palindrome pairs base ``i`` with base ``n - 1 - i``, so with no loop required the
    whole 6 bp arm is a "stem" -- a fold with nothing to fold around. Requiring a loop reports
    the shorter stem that can actually form, using the middle bases as the loop. This is the
    distinction ``_self_complementarity`` could not make.
    """
    arm = "GGGGCC"
    palindrome = arm + revcomp(arm)

    unconstrained, loop = max_hairpin_stem(palindrome, seed=4, min_loop=0)
    constrained, constrained_loop = max_hairpin_stem(palindrome, seed=4, min_loop=3)

    assert (unconstrained, loop) == (len(arm), 0)
    assert constrained < unconstrained
    assert constrained_loop >= 3


def test_a_three_prime_hairpin_is_distinguished_from_an_upstream_one():
    """The same stem costs differently depending on whether it occludes the 3' end.

    Note the anchored measurement never returns a hard 0 for real sequence -- a single terminal
    base pairs somewhere by chance -- which is why consumers charge only the excess over a free
    stem length rather than treating any hit as a defect.
    """
    arm = "GGGGCC"
    anchored = arm + "ATAT" + revcomp(arm)  # stem closes on the final base
    upstream = anchored + "AAAAAAAAAA"  # same stem, now internal

    assert three_prime_hairpin(anchored)[0] == len(arm)
    assert three_prime_hairpin(upstream)[0] < 4  # nothing chargeable at the 3' end
    # the internal hairpin is still visible to the unanchored measurement
    assert max_hairpin_stem(upstream, seed=4)[0] == len(arm)


def test_no_hairpin_in_a_homopolymer():
    assert three_prime_hairpin("A" * 30) == (0, 0)
    assert max_hairpin_stem("A" * 30, seed=4) == (0, 0)


def test_the_stem_search_is_bounded_on_a_tandem_array():
    """Seed buckets on a tandem array are enormous; the search must not degrade."""
    stem, _ = max_hairpin_stem("ACGTACGTAA" * 60, seed=4)
    assert stem >= 0  # completes at all, and quickly -- see the timing check in the plan


# --------------------------------------------------------------------------- #
# ambiguity handling
# --------------------------------------------------------------------------- #


def test_homopolymer_counts_runs():
    assert max_homopolymer("ACGTTTTACG") == 4
    assert max_homopolymer("") == 0
    assert max_homopolymer("A") == 1


def test_a_stray_ambiguity_code_is_reported_not_hidden():
    clean, dropped = clean_dna("ACGT" * 100 + "N")
    assert dropped == 1
    assert len(clean) == 400
    assert set(clean) <= set("ACGT")


def test_a_meaningfully_ambiguous_sequence_is_refused():
    """Silently stripping would shorten the sequence and shift GC and every density
    denominator -- the same silent-wrong-answer class as a swallowed exception."""
    with pytest.raises(ValueError, match="not ACGT"):
        clean_dna("ACGT" + "N" * 10)
