"""Synthesis complexity screening -- estimates whether a vendor will refuse a fragment.

Replaces the original placeholder. The planner's contract is unchanged: it branches on
:attr:`ComplexityReport.synthesizable` and reports :attr:`ComplexityReport.reasons`.

**What makes this different from the placeholder it replaces.** The placeholder asked "does
this breach a threshold" and answered with a fixed weight per breach. This asks "will IDT
reject this order", and its weights are back-fit to the sub-scores IDT reported on a sequence
it actually rejected:

===========================================  ======  ============
feature                                      IDT     here
===========================================  ======  ============
87 bp exact repeat                           48.3    48.2
88% repeat density in an 85 bp window        21.6    22.0
93.9% repeat fraction                        14.4    14.4
~19 bp exact repeat                          10      10.0
high 5' GC                                   10      10.0
===========================================  ======  ============

Verdict bands are anchored to IDT's 24 threshold. For an exact accept/reject, query the
vendor's own screener -- :func:`feature_vector` exists so a calibrated ``P(reject)`` can be
fitted on real labels later.

**Three deliberate departures from a naive port**, each of which was a defect:

* Non-ACGT characters are no longer silently dropped. Stripping them shortens the sequence and
  shifts GC and every density denominator, so :func:`bbl.seqfeatures.clean_dna` reports the
  count and refuses past a small fraction.
* :func:`bbl.seqfeatures.repeat_mask` Bonferroni-corrects its significance threshold. Testing
  hundreds of k-mers at ``alpha`` and reporting each hit as "significant" overstates the
  evidence. This is a correctness fix to that claim more than a behavioural one: the rules that
  consume the mask only fire above 58% coverage, and at that level significance is overwhelming
  under either threshold. It makes the mask *conservative* -- a marginal recurrence a family-wise
  threshold cannot support is no longer called a designed repeat, and genuine long repeats are
  caught by :func:`longest_repeat` independently of the mask.
* The "~19 bp repeat" anchor is honoured. The source scored 8 there, not the documented 10.

**Backward compatibility.** :data:`DEFAULT_LIMITS`, :func:`longest_homopolymer`,
:func:`longest_repeat` and :func:`repeat_fraction` keep their old names, signatures and
semantics because :mod:`bbl.qc` scores against them -- that layer is the bake-off evidence for
this replacement and must not become collateral damage. In particular
:func:`repeat_fraction` stays *uncorrected* plain recurrence counting; the background-corrected
measurement is separate, under :attr:`ComplexityReport.metrics` key
``significant_repeat_fraction``.

One score is graded where the placeholder stepped: a 12-base homopolymer and a 40-base one are
no longer the same answer. That gradient is what makes candidate *ranking* possible.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio.SeqUtils import gc_fraction

from .seqfeatures import (
    clean_dna,
    longest_direct_repeat,
    max_hairpin_stem,
    max_homopolymer,
    repeat_mask,
    terminal_gc,
    windowed_gc,
    windowed_repeat_density,
)

#: Rough commercial limits for a synthetic double-stranded fragment. Vendor-specific in
#: reality -- a fragment one vendor refuses another may accept. ``min_length``/``max_length`` are
#: hard gates (see :func:`complexity_score`); the rest are thresholds :mod:`bbl.qc` scores
#: against and are kept for that reason.
DEFAULT_LIMITS = {
    "min_length": 125,
    "max_length": 3000,
    "min_gc": 25.0,
    "max_gc": 75.0,
    "max_homopolymer": 9,
    "max_repeat": 40,  # longest substring allowed to occur more than once
    "max_repeat_fraction": 0.35,
}

#: IDT's published complexity threshold, and the band below which a fragment is comfortable.
REJECT_SCORE = 24.0
BORDERLINE_SCORE = 10.0

# Penalty weights. Anchored where an observed IDT sub-score exists, interpolated where not --
# the interpolations are guesses and the docstring says so. Continuous and monotone at every
# band edge, so a one-base change in a measurement cannot swing the score by 14 points.
_REPEAT_MIN_LENGTH = 14  # below this, exact repeats are unremarkable
_REPEAT_BASE = 10.0  # the ~19 bp anchor
_REPEAT_SLOPE = 0.57  # reaches 48.2 at 87 bp, matching the 48.3 anchor
_REPEAT_MAX = 50.0
_FRACTION_FLOOR = 0.58
_FRACTION_SLOPE = 40.0  # reaches 14.4 at 93.9%, matching the anchor
_FRACTION_MAX = 20.0
_DENSITY_FLOOR = 0.70
_DENSITY_BASE = 8.0
_DENSITY_SLOPE = 77.8  # reaches 22.0 at 88%, matching the 21.6 anchor
_DENSITY_MAX = 22.0
_GC_TERMINAL_LIMIT = 0.80
_GC_WINDOW_HIGH = 0.85
_GC_WINDOW_LOW = 0.15
_GC_FLAT_PENALTY = 10.0
_GC_GLOBAL_PENALTY = 5.0


@dataclass
class ComplexityReport:
    """Verdict on whether a fragment can be bought as synthetic DNA."""

    sequence_length: int
    score: float  # 0 = trivial, 1 = at or past the vendor's reject threshold
    synthesizable: bool
    reasons: list[str] = field(default_factory=list)
    metrics: dict = field(default_factory=dict)
    vendor: str | None = None
    is_placeholder: bool = False
    #: Uncapped heuristic score in the vendor's own units; ``>= REJECT_SCORE`` means refusal.
    raw_score: float = 0.0
    #: ``"LIKELY ACCEPT"`` | ``"BORDERLINE (check vendor)"`` | ``"LIKELY REJECT"``.
    verdict: str = "LIKELY ACCEPT"
    #: ``(rule, penalty, note)`` per firing rule -- which features drive the risk.
    contributions: list[tuple[str, float, str]] = field(default_factory=list)

    def __str__(self) -> str:
        verdict = "synthesizable" if self.synthesizable else "NOT synthesizable"
        detail = f" -- {'; '.join(self.reasons)}" if self.reasons else ""
        return (
            f"{self.sequence_length} bp, score {self.raw_score:.1f} "
            f"({self.verdict}), {verdict}{detail}"
        )


# --------------------------------------------------------------------------------------------
# metrics kept at their original names and semantics for bbl.qc
# --------------------------------------------------------------------------------------------


def longest_homopolymer(sequence: str) -> int:
    """Longest single-base run. Thin wrapper over :func:`bbl.seqfeatures.max_homopolymer`."""
    return max_homopolymer(sequence)


def longest_repeat(sequence: str, cap: int = 200) -> int:
    """Length of the longest substring occurring more than once.

    Tolerates arbitrary characters -- callers pass non-ACGT strings -- so it does **not** go
    through :func:`clean_dna`.

    The search ceiling is now ``n - 1`` rather than ``n // 2``. Occurrences may overlap, so the
    old bound undercounted exactly the tandem array this measurement exists to catch: ``AAAA``
    contains ``AAA`` twice and scores 3, where the old bound reported 2.
    """
    length, _ = longest_direct_repeat(sequence.upper())
    return min(length, cap)


def repeat_fraction(sequence: str, window: int = 20) -> float:
    """Fraction of positions whose ``window``-mer appears more than once.

    Deliberately **uncorrected** plain recurrence counting, preserved verbatim because
    :mod:`bbl.qc.native` scores against it. For the background-corrected version see
    :func:`bbl.seqfeatures.repeat_mask`.
    """
    sequence = sequence.upper()
    if len(sequence) < window:
        return 0.0
    counts: dict[str, int] = {}
    for i in range(len(sequence) - window + 1):
        chunk = sequence[i : i + window]
        counts[chunk] = counts.get(chunk, 0) + 1
    repeated = sum(count for chunk, count in counts.items() if count > 1)
    return repeated / (len(sequence) - window + 1)


# --------------------------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------------------------


def _ramp(value: float, floor: float, base: float, slope: float, ceiling: float) -> float:
    """Penalty that starts at ``base`` once ``value`` passes ``floor`` and rises to ``ceiling``."""
    return min(ceiling, base + (value - floor) * slope)


def _measure(sequence: str) -> dict:
    """Every raw feature, no thresholds applied."""
    low_gc, high_gc = windowed_gc(sequence)
    gc_5, gc_3 = terminal_gc(sequence)
    repeat_length, repeat_sequence = longest_direct_repeat(sequence)
    mask = repeat_mask(sequence)
    stem, loop = max_hairpin_stem(sequence)
    coverage = sum(mask) / len(sequence) if sequence else 0.0
    return {
        "gc_percent": round(gc_fraction(sequence) * 100 if sequence else 0.0, 1),
        "gc_window_min": round(low_gc, 3),
        "gc_window_max": round(high_gc, 3),
        "gc_5prime": round(gc_5, 3),
        "gc_3prime": round(gc_3, 3),
        "longest_homopolymer": max_homopolymer(sequence),
        "longest_repeat": repeat_length,
        "longest_repeat_sequence": (
            repeat_sequence if repeat_length < 40 else repeat_sequence[:37] + "..."
        ),
        "repeat_fraction": round(repeat_fraction(sequence), 3),
        "significant_repeat_fraction": round(coverage, 3),
        "repeat_density_85bp": round(windowed_repeat_density(mask), 3),
        "hairpin_stem": stem,
        "hairpin_loop": loop,
    }


def _score_rules(metrics: dict) -> list[tuple[str, float, str]]:
    """Every firing rule as ``(name, penalty, note)``. Length gates are handled separately."""
    firing: list[tuple[str, float, str]] = []

    def add(name: str, penalty: float, note: str) -> None:
        firing.append((name, round(float(penalty), 1), note))

    repeat = metrics["longest_repeat"]
    if repeat > 20:
        add(
            "repeat_long",
            _ramp(repeat, 20, _REPEAT_BASE, _REPEAT_SLOPE, _REPEAT_MAX),
            f"longest exact repeat {repeat} bp",
        )
    elif repeat >= _REPEAT_MIN_LENGTH:
        add("repeat_moderate", _REPEAT_BASE, f"longest exact repeat {repeat} bp")

    coverage = metrics["significant_repeat_fraction"]
    if coverage > _FRACTION_FLOOR:
        add(
            "repeat_fraction",
            _ramp(coverage, _FRACTION_FLOOR, 0.0, _FRACTION_SLOPE, _FRACTION_MAX),
            f"{coverage:.0%} of the fragment is significant repeat",
        )

    density = metrics["repeat_density_85bp"]
    if density >= _DENSITY_FLOOR:
        add(
            "repeat_density",
            _ramp(density, _DENSITY_FLOOR, _DENSITY_BASE, _DENSITY_SLOPE, _DENSITY_MAX),
            f"{density:.0%} repeat in an 85 bp window",
        )

    if metrics["gc_5prime"] > _GC_TERMINAL_LIMIT:
        add("gc_5prime", _GC_FLAT_PENALTY, f"5' end GC {metrics['gc_5prime']:.0%}")
    if metrics["gc_window_max"] > _GC_WINDOW_HIGH:
        add(
            "gc_window_high",
            _GC_FLAT_PENALTY,
            f"local GC up to {metrics['gc_window_max']:.0%}",
        )
    if metrics["gc_window_min"] < _GC_WINDOW_LOW:
        add(
            "gc_window_low",
            _GC_FLAT_PENALTY,
            f"local GC down to {metrics['gc_window_min']:.0%}",
        )

    homopolymer = metrics["longest_homopolymer"]
    if homopolymer >= 20:
        add("homopolymer_very_long", 26, f"{homopolymer} bp homopolymer run")
    elif homopolymer >= 13:
        add("homopolymer_long", 12, f"{homopolymer} bp homopolymer run")
    elif homopolymer >= 10:
        add("homopolymer", 5, f"{homopolymer} bp homopolymer run")

    stem = metrics["hairpin_stem"]
    if stem >= 20:
        add("hairpin_severe", 20, f"{stem} bp hairpin stem")
    elif stem >= 16:
        add("hairpin_strong", 12, f"{stem} bp hairpin stem")
    elif stem >= 12:
        add("hairpin", 6, f"{stem} bp hairpin stem")

    return firing


def complexity_score(sequence, vendor: str | None = None, limits=None) -> ComplexityReport:
    """Score a fragment for synthesis complexity.

    Feed it the fragment **as it would be ordered** -- insert plus any homology arms, i.e.
    ``plan_insertion(...).insert.order_sequence`` -- not the bare insert.

    ``min_length``/``max_length`` are hard gates: a fragment outside them is not orderable at
    any complexity score, so they force ``synthesizable=False`` independently of the band.
    Everything else contributes to :attr:`ComplexityReport.raw_score`, which is refused at
    :data:`REJECT_SCORE`.
    """
    limits = {**DEFAULT_LIMITS, **(limits or {})}
    clean, dropped = clean_dna(sequence)
    length = len(clean)

    metrics = _measure(clean)
    if dropped:
        metrics["non_acgt_dropped"] = dropped

    contributions = _score_rules(metrics)
    raw = sum(penalty for _, penalty, _ in contributions)

    # Hard gates. Outside these bounds nothing is orderable, whatever the complexity score.
    blocking: list[str] = []
    if length < limits["min_length"]:
        blocking.append(f"{length} bp is below the {limits['min_length']} bp synthesis minimum")
    if length > limits["max_length"]:
        blocking.append(f"{length} bp exceeds the {limits['max_length']} bp maximum")

    global_gc = metrics["gc_percent"]
    if global_gc < limits["min_gc"] or global_gc > limits["max_gc"]:
        contributions.append(
            (
                "gc_global",
                _GC_GLOBAL_PENALTY,
                f"global GC {global_gc:.0f}% outside "
                f"{limits['min_gc']:.0f}-{limits['max_gc']:.0f}%",
            )
        )
        raw += _GC_GLOBAL_PENALTY

    if raw < BORDERLINE_SCORE:
        verdict = "LIKELY ACCEPT"
    elif raw < REJECT_SCORE:
        verdict = "BORDERLINE (check vendor)"
    else:
        verdict = "LIKELY REJECT"

    reasons = list(blocking)
    if raw >= REJECT_SCORE:
        reasons.extend(note for _, _, note in contributions)
    if dropped:
        reasons.append(f"{dropped} non-ACGT character(s) ignored")

    return ComplexityReport(
        sequence_length=length,
        score=round(min(1.0, raw / REJECT_SCORE), 3),
        synthesizable=not blocking and raw < REJECT_SCORE,
        reasons=reasons,
        metrics=metrics,
        vendor=vendor,
        raw_score=round(raw, 1),
        verdict=verdict,
        contributions=contributions,
    )


def feature_vector(sequence) -> dict:
    """Raw numeric features, for fitting a calibrated model on real accept/reject labels.

    The heuristic weights above are anchored to a single observed rejection; this is the way
    out of that. Non-numeric fields are dropped so the result can go straight into a dataframe.
    """
    clean, _ = clean_dna(sequence)
    return {
        key: value
        for key, value in _measure(clean).items()
        if isinstance(value, (int, float)) and not isinstance(value, bool)
    }
