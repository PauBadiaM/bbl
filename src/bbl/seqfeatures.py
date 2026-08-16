"""Pure sequence measurements -- no thresholds, no verdicts, no scoring.

This module is deliberately the *bottom* of the stack: it measures, and callers decide what a
measurement means. :mod:`bbl.complexity` turns these into a synthesisability verdict,
:mod:`bbl.pcr` turns some of them into oligo penalties, and :mod:`bbl.qc.native` scores them
against Proto. Keeping the primitives here is what stops those three from drifting -- the
longest-homopolymer routine previously existed as three independent copies.

Two measurements here are new to ``bbl`` and worth calling out, because the incumbent code has
no equivalent and therefore no way to see the failure they catch:

* :func:`windowed_gc` / :func:`terminal_gc` -- **positional** GC. Every other GC number in the
  codebase is global, and a fragment at 50% global GC can hide a 90% GC 50-mer. High 5' GC is
  a real vendor rejection driver.
* :func:`max_hairpin_stem` -- secondary structure with a *loop* constraint, so a 3'-blocking
  hairpin can be told apart from a harmless inverted repeat at the far end of the oligo.

:func:`repeat_mask` is background-corrected: a k-mer counts as repetitive only if its
recurrence is significant under a composition-matched i.i.d. null. Plain recurrence counting
(what :func:`bbl.complexity.repeat_fraction` did) reports AT-rich low-complexity DNA as
"repetitive" purely from chance k-mer collisions.
"""

from __future__ import annotations

import math
from collections import Counter

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCATGCA")

#: Fraction of non-ACGT characters tolerated before :func:`clean_dna` refuses. Silently
#: dropping ambiguity codes shortens the sequence and shifts GC and every density denominator,
#: so a sequence that is meaningfully ambiguous must fail loudly instead.
AMBIGUITY_TOLERANCE = 0.02


def revcomp(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def complement(base: str) -> str:
    return base.translate(_COMPLEMENT)


def clean_dna(sequence, *, tolerance: float = AMBIGUITY_TOLERANCE) -> tuple[str, int]:
    """Upper-cased ACGT-only sequence, plus how many characters were dropped.

    Returns the count rather than discarding it so callers can report it. Raises when the
    dropped fraction exceeds ``tolerance``: past that point every downstream measurement is
    computed over a sequence that is not the one the caller passed, which is the silent-wrong
    -answer failure mode (cf. ``enumerate_cut_sites`` reporting zero cut sites when handed a
    plain ``str``).
    """
    raw = str(sequence).upper()
    kept = "".join(character for character in raw if character in "ACGT")
    dropped = len(raw) - len(kept)
    if raw and dropped > tolerance * len(raw):
        raise ValueError(
            f"{dropped} of {len(raw)} characters are not ACGT "
            f"({dropped / len(raw):.1%}, tolerance {tolerance:.1%}); "
            "resolve the ambiguity codes before measuring"
        )
    return kept, dropped


# --------------------------------------------------------------------------- composition ---


def gc_percent(sequence: str) -> float:
    """GC as a percentage. Counts literal G/C only -- see :mod:`bbl.qc.native` for why the
    Biopython route is kept separately."""
    sequence = sequence.upper()
    if not sequence:
        return 0.0
    return 100.0 * (sequence.count("G") + sequence.count("C")) / len(sequence)


def windowed_gc(sequence: str, window: int = 50) -> tuple[float, float]:
    """``(min, max)`` GC fraction over all sliding windows. Global GC hides local extremes."""
    sequence = sequence.upper()
    n = len(sequence)
    if n == 0:
        return 0.0, 0.0
    if n <= window:
        fraction = gc_percent(sequence) / 100.0
        return fraction, fraction
    count = sum(1 for character in sequence[:window] if character in "GC")
    low = high = count / window
    for i in range(window, n):
        count += (sequence[i] in "GC") - (sequence[i - window] in "GC")
        fraction = count / window
        low = min(low, fraction)
        high = max(high, fraction)
    return low, high


def terminal_gc(sequence: str, end_length: int = 40) -> tuple[float, float]:
    """``(5' GC, 3' GC)`` as fractions over the terminal ``end_length`` bases."""
    sequence = sequence.upper()
    if len(sequence) < end_length:
        fraction = gc_percent(sequence) / 100.0
        return fraction, fraction
    return (
        gc_percent(sequence[:end_length]) / 100.0,
        gc_percent(sequence[-end_length:]) / 100.0,
    )


def max_homopolymer(sequence: str) -> int:
    """Longest single-base run. The one primitive that had three independent copies."""
    sequence = sequence.upper()
    if not sequence:
        return 0
    best = run = 1
    for previous, current in zip(sequence, sequence[1:]):
        run = run + 1 if current == previous else 1
        best = max(best, run)
    return best


# -------------------------------------------------------------------------------- repeats ---


def longest_direct_repeat(sequence: str, min_length: int = 1) -> tuple[int, str]:
    """``(length, substring)`` of the longest substring occurring at least twice.

    Binary search on length: if a repeat of length ``L`` exists then one of length ``L - 1``
    does too, so the predicate is monotone and bisection is valid.

    The upper bound is ``n - 1``, not ``n // 2``. Occurrences may **overlap** -- ``AAAA``
    contains ``AAA`` twice, at offsets 0 and 1 -- so a ``n // 2`` ceiling undercounts exactly
    the tandem-array case this measurement exists to catch.
    """
    sequence = sequence.upper()
    n = len(sequence)

    def repeat_of_length(length: int) -> str | None:
        seen: set[str] = set()
        for i in range(n - length + 1):
            chunk = sequence[i : i + length]
            if chunk in seen:
                return chunk
            seen.add(chunk)
        return None

    low, high, best = max(1, min_length), max(0, n - 1), None
    while low <= high:
        middle = (low + high) // 2
        found = repeat_of_length(middle)
        if found:
            best, low = found, middle + 1
        else:
            high = middle - 1
    return (len(best), best) if best else (0, "")


def poisson_sf(count: int, mean: float) -> float:
    """``P(X >= count)`` for ``X ~ Poisson(mean)``."""
    if count <= 0:
        return 1.0
    if mean <= 0:
        return 0.0
    if count > mean + 40 * math.sqrt(mean + 1):
        return 0.0
    term = math.exp(-mean)
    cdf = term
    for i in range(1, count):
        term *= mean / i
        cdf += term
    return max(0.0, 1.0 - cdf)


def repeat_mask(sequence: str, k: int = 8, alpha: float = 1e-3) -> list[bool]:
    """Per-position mask of *significantly* repeated k-mers.

    A k-mer's positions are marked only when it occurs at least twice **and** that many
    occurrences are surprising under an i.i.d. null matched to the sequence's own base
    composition. Without the null, AT-rich or low-complexity DNA scores as repetitive from
    chance collisions alone.

    ``alpha`` is Bonferroni-corrected by the number of distinct k-mers examined. Uncorrected,
    a 3 kb fragment tests thousands of hypotheses at ``alpha`` and expects a handful of false
    repeat calls -- which would then be reported as designed repeats.
    """
    sequence = sequence.upper()
    n = len(sequence)
    mask = [False] * n
    if n < k:
        return mask

    starts = n - k + 1
    composition = Counter(sequence)
    total = sum(composition.get(base, 0) for base in "ACGT") or 1
    probability = {base: composition.get(base, 0) / total for base in "ACGT"}

    positions: dict[str, list[int]] = {}
    for i in range(starts):
        positions.setdefault(sequence[i : i + k], []).append(i)

    threshold = alpha / max(1, len(positions))
    for kmer, locations in positions.items():
        if len(locations) < 2:
            continue
        expected = starts
        for base in kmer:
            expected *= probability.get(base, 0.0)
        if poisson_sf(len(locations), expected) < threshold:
            for i in locations:
                for j in range(i, i + k):
                    mask[j] = True
    return mask


def windowed_repeat_density(mask: list[bool], window: int = 85) -> float:
    """Highest fraction of masked positions in any window -- a local repeat concentration."""
    n = len(mask)
    if n == 0:
        return 0.0
    if n <= window:
        return sum(mask) / n
    current = sum(mask[:window])
    best = current
    for i in range(window, n):
        current += mask[i] - mask[i - window]
        best = max(best, current)
    return best / window


# ---------------------------------------------------------------------- secondary structure ---

#: Stop extending once a stem this long is found. Every consumer treats "at least this long" as
#: the worst band, so the exact length past it changes no decision -- and the search is
#: superlinear on tandem arrays, where seed buckets are enormous.
STEM_CAP = 30

#: Skip seeds whose reverse-complement occurs more often than this. On an 8x tandem array a
#: single seed can pair with hundreds of partners, none of which change the answer once
#: :data:`STEM_CAP` is reachable.
MAX_SEED_PARTNERS = 64


def max_hairpin_stem(
    sequence: str,
    seed: int = 6,
    min_loop: int = 3,
    max_loop: int = 200,
    stem_cap: int = STEM_CAP,
    max_seed_partners: int = MAX_SEED_PARTNERS,
) -> tuple[int, int]:
    """``(stem, loop)`` of the strongest hairpin: seed on a reverse-complement match, extend.

    The ``min_loop`` constraint is the point. A self-complementary stretch is only a hairpin if
    the intervening bases can actually form a loop; without the bound this cannot distinguish a
    fold-back that blocks a primer's 3' end from an inverted repeat at the opposite end, which
    is why ``pcr._self_complementarity`` could not be used for ranking.

    Returns a lower bound when ``stem_cap`` or ``max_seed_partners`` truncates the search.
    """
    sequence = sequence.upper()
    n = len(sequence)
    if n < 2 * seed:
        return 0, 0

    index: dict[str, list[int]] = {}
    for i in range(n - seed + 1):
        index.setdefault(sequence[i : i + seed], []).append(i)

    best_stem = best_loop = 0
    for i in range(n - seed + 1):
        partners = index.get(revcomp(sequence[i : i + seed]), ())
        if len(partners) > max_seed_partners:
            partners = partners[:max_seed_partners]
        for j in partners:
            if j < i + seed:
                continue
            loop = j - (i + seed)
            if loop < min_loop or loop > max_loop:
                continue
            stem = seed
            left, right = i - 1, j + seed
            while left >= 0 and right < n and complement(sequence[left]) == sequence[right]:
                stem += 1
                left -= 1
                right += 1
            if stem > best_stem:
                best_stem, best_loop = stem, loop
                if best_stem >= stem_cap:
                    return best_stem, best_loop
    return best_stem, best_loop


def three_prime_hairpin(sequence: str, min_loop: int = 3) -> tuple[int, int]:
    """``(stem, loop)`` of the longest hairpin **anchored at the 3' terminal base**.

    A hairpin that occludes the 3' end stops extension outright; one further upstream mostly
    costs efficiency. Consumers weight the two differently, so they must be measurable
    separately -- and :func:`max_hairpin_stem` cannot answer this, since it reports the strongest
    hairpin anywhere without regard to where the arms sit.

    Geometry: the right arm is exactly the 3'-terminal ``stem`` bases, so a stem of length ``L``
    exists iff ``revcomp(sequence[-L:])`` occurs at some offset ``a`` with at least ``min_loop``
    bases between ``a + L`` and ``len(sequence) - L``. Scanning ``L`` downwards returns the
    longest stem, and ``rfind`` picks the tightest loop for it. Intended for oligos, where the
    quadratic term is negligible.
    """
    sequence = sequence.upper()
    n = len(sequence)
    for stem in range((n - min_loop) // 2, 0, -1):
        arm = revcomp(sequence[n - stem :])
        upstream = sequence[: n - stem - min_loop]
        offset = upstream.rfind(arm)
        if offset >= 0:
            return stem, (n - stem) - (offset + stem)
    return 0, 0
