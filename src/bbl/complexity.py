"""Synthesisability screening -- **PLACEHOLDER**.

The lab will supply the real scoring function. Everything here exists so the planner has a
stable interface to branch on; replace the body of :func:`complexity_score` (or the whole
module) and nothing downstream changes.

The heuristics below are deliberately crude but not vacuous: they catch the failure mode that
actually matters in this library, which is tandem-repeat arrays. Vendors refuse those, and the
8x BoxB array is exactly such a case.

To swap in the real function, keep the return type: the planner branches on
``ComplexityReport.synthesizable`` and reports ``reasons``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio.SeqUtils import gc_fraction

#: Rough commercial limits for a synthetic double-stranded fragment. Vendor-specific in
#: reality -- a fragment one vendor refuses another may accept.
DEFAULT_LIMITS = {
    "min_length": 125,
    "max_length": 3000,
    "min_gc": 25.0,
    "max_gc": 75.0,
    "max_homopolymer": 9,
    "max_repeat": 40,  # longest substring allowed to occur more than once
    "max_repeat_fraction": 0.35,
}


@dataclass
class ComplexityReport:
    """Verdict on whether a fragment can be bought as synthetic DNA."""

    sequence_length: int
    score: float  # 0 = trivial, 1 = certainly refused
    synthesizable: bool
    reasons: list[str] = field(default_factory=list)
    metrics: dict = field(default_factory=dict)
    vendor: str | None = None
    is_placeholder: bool = True

    def __str__(self) -> str:
        verdict = "synthesizable" if self.synthesizable else "NOT synthesizable"
        detail = f" -- {'; '.join(self.reasons)}" if self.reasons else ""
        return f"{self.sequence_length} bp, score {self.score:.2f}, {verdict}{detail}"


def longest_homopolymer(sequence: str) -> int:
    best = run = 1
    for previous, current in zip(sequence, sequence[1:]):
        run = run + 1 if current == previous else 1
        best = max(best, run)
    return best if sequence else 0


def longest_repeat(sequence: str, cap: int = 200) -> int:
    """Length of the longest substring occurring more than once.

    Binary search on length with a k-mer set per probe: O(n log n) sets, fine at fragment
    scale. This is the signal that catches tandem arrays.
    """
    sequence = sequence.upper()
    low, high, best = 1, min(len(sequence) // 2, cap), 0
    while low <= high:
        mid = (low + high) // 2
        seen, hit = set(), False
        for i in range(len(sequence) - mid + 1):
            chunk = sequence[i : i + mid]
            if chunk in seen:
                hit = True
                break
            seen.add(chunk)
        if hit:
            best, low = mid, mid + 1
        else:
            high = mid - 1
    return best


def repeat_fraction(sequence: str, window: int = 20) -> float:
    """Fraction of positions whose ``window``-mer appears more than once."""
    sequence = sequence.upper()
    if len(sequence) < window:
        return 0.0
    counts: dict[str, int] = {}
    for i in range(len(sequence) - window + 1):
        chunk = sequence[i : i + window]
        counts[chunk] = counts.get(chunk, 0) + 1
    repeated = sum(count for chunk, count in counts.items() if count > 1)
    return repeated / (len(sequence) - window + 1)


def complexity_score(sequence: str, vendor: str | None = None, limits=None) -> ComplexityReport:
    """PLACEHOLDER synthesisability check. Replace with the lab's scoring function.

    Feed it the fragment **as it would be ordered** -- insert plus any homology arms, i.e.
    ``plan_insertion(...).insert.order_sequence`` -- not the bare insert.
    """
    limits = {**DEFAULT_LIMITS, **(limits or {})}
    sequence = str(sequence).upper()
    length = len(sequence)

    gc = gc_fraction(sequence) * 100 if length else 0.0
    homopolymer = longest_homopolymer(sequence)
    repeat = longest_repeat(sequence)
    fraction = repeat_fraction(sequence)
    metrics = {
        "gc_percent": round(gc, 1),
        "longest_homopolymer": homopolymer,
        "longest_repeat": repeat,
        "repeat_fraction": round(fraction, 3),
    }

    reasons, penalty = [], 0.0
    if length < limits["min_length"]:
        reasons.append(f"{length} bp is below the {limits['min_length']} bp synthesis minimum")
        penalty += 0.5
    if length > limits["max_length"]:
        reasons.append(f"{length} bp exceeds the {limits['max_length']} bp maximum")
        penalty += 0.5
    if gc < limits["min_gc"] or gc > limits["max_gc"]:
        reasons.append(f"GC {gc:.0f}% is outside {limits['min_gc']}-{limits['max_gc']}%")
        penalty += 0.3
    if homopolymer > limits["max_homopolymer"]:
        reasons.append(f"homopolymer run of {homopolymer}")
        penalty += 0.3
    if repeat > limits["max_repeat"]:
        reasons.append(f"{repeat} bp repeated within the fragment")
        penalty += 0.5
    if fraction > limits["max_repeat_fraction"]:
        reasons.append(f"{fraction:.0%} of the fragment is repetitive")
        penalty += 0.4

    return ComplexityReport(
        sequence_length=length,
        score=min(1.0, penalty),
        synthesizable=not reasons,
        reasons=reasons,
        metrics=metrics,
        vendor=vendor,
    )
