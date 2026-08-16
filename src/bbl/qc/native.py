"""Hand-written scorers -- the incumbent side of the bake-off.

These deliberately wrap **what `bbl` already does today** rather than reimplementing Proto's
formulas. That is the whole point: if a native scorer mirrored Proto's arithmetic, the
comparison would agree by construction and tell us nothing.

So:

* ``gc_content`` uses Biopython's :func:`gc_fraction`, exactly as ``pcr.py`` and
  ``complexity.py`` do. Proto counts only literal ``G``/``C``; Biopython has its own
  ambiguity handling. On plain ACGT they must agree, and on ambiguous input they may not --
  that divergence is a finding, not a defect.
* ``max_homopolymer`` reuses :func:`bbl.complexity.longest_homopolymer`.
* ``longest_repeat`` / ``repeat_fraction`` reuse the ``complexity`` implementations, which are
  the checks that actually catch the 8x BoxB array.
* ``restriction_sites`` goes through :mod:`bbl.enzymes` (Biopython ``Restriction``), which
  understands ambiguity codes, isoschizomers and cut coordinates. Proto cannot express this at
  all -- see :mod:`bbl.qc.proto_backend`.

Penalties here are mostly **step functions**, because that is what the incumbent code does:
``complexity.py`` appends a reason and adds a fixed weight, it does not scale with severity.
Proto's are graded. Expect ``penalty_delta`` to be nonzero even when measurements match.
"""

from __future__ import annotations

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

from ..complexity import longest_homopolymer, longest_repeat, repeat_fraction
from ..enzymes import build_pool, enumerate_cut_sites

#: A scorer returns ``(measurement, penalty, detail)``.
#: ``measurement`` may be ``None`` when the quantity is not a scalar.


def _step(inside: bool) -> float:
    """Incumbent penalty shape: 0.0 when acceptable, 1.0 when not. No gradation."""
    return 0.0 if inside else 1.0


def gc_content(sequence: str, params: dict) -> tuple[float, float, dict]:
    """GC percent via Biopython, as the rest of ``bbl`` computes it."""
    if not sequence:
        return 0.0, 1.0, {"empty": True}
    value = gc_fraction(sequence) * 100.0
    low, high = params["min_gc"], params["max_gc"]
    return value, _step(low <= value <= high), {"min_gc": low, "max_gc": high}


def max_homopolymer(sequence: str, params: dict) -> tuple[float, float, dict]:
    """Longest single-base run, via ``complexity.longest_homopolymer``."""
    value = float(longest_homopolymer(sequence))
    limit = params["max_length"]
    return value, _step(value <= limit), {"max_length": limit}


def sequence_length(sequence: str, params: dict) -> tuple[float, float, dict]:
    value = float(len(sequence))
    low, high = params["min_length"], params["max_length"]
    return value, _step(low <= value <= high), {"min_length": low, "max_length": high}


def specific_kmer(sequence: str, params: dict) -> tuple[float, float, dict]:
    """Overlapping count of one literal k-mer, reported as a frequency.

    Frequency rather than count so it is directly comparable with Proto's
    ``specific-kmer-frequency`` in ``frequency`` mode, which divides by the number of
    k-mer start positions.
    """
    kmer = params["kmer"].upper()
    seq = sequence.upper()
    k = len(kmer)
    positions = len(seq) - k + 1
    if positions <= 0:
        return 0.0, 1.0, {"kmer": kmer, "count": 0, "positions": 0}
    count = sum(1 for i in range(positions) if seq[i : i + k] == kmer)
    frequency = count / positions
    low, high = params["min_value"], params["max_value"]
    return (
        frequency,
        _step(low <= frequency <= high),
        {"kmer": kmer, "count": count, "positions": positions},
    )


def longest_repeat_check(sequence: str, params: dict) -> tuple[float, float, dict]:
    """Longest substring occurring more than once -- the tandem-array signal."""
    value = float(longest_repeat(sequence.upper()))
    limit = params["max_repeat"]
    return value, _step(value <= limit), {"max_repeat": limit}


def repeat_fraction_check(sequence: str, params: dict) -> tuple[float, float, dict]:
    """Fraction of positions whose window-mer recurs."""
    value = float(repeat_fraction(sequence.upper(), window=params["window"]))
    limit = params["max_repeat_fraction"]
    return value, _step(value <= limit), {"max_repeat_fraction": limit, "window": params["window"]}


def restriction_sites(sequence: str, params: dict) -> tuple[float, float, dict]:
    """Occurrences of one restriction enzyme's site, via Biopython ``Restriction``.

    This is the check Proto structurally cannot do: ``AccI`` is ``GT^MKAC``, and Proto's
    sequence alphabet is ``frozenset("ACGT")``, so an ambiguous recognition site is not
    expressible as a k-mer. Reported as a count, not a frequency -- a count is what matters
    for "does this enzyme cut exactly once".
    """
    name = params["enzyme"]
    pool = build_pool([name])
    # Biopython's Restriction.search requires a Seq, not a str -- with a bare string it raises
    # and enumerate_cut_sites swallows the exception, silently reporting zero sites.
    # Linear by default: the caller passes a fragment, not a whole plasmid.
    sites = enumerate_cut_sites(
        Seq(sequence.upper()), pool, circular=params.get("circular", False)
    )
    count = float(len(sites.get(name, [])))
    low, high = params["min_sites"], params["max_sites"]
    return count, _step(low <= count <= high), {"enzyme": name, "sites": int(count)}
