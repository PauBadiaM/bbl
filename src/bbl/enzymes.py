"""Restriction-enzyme site enumeration and end-compatibility logic.

Coordinate convention
---------------------
``Bio.Restriction``'s ``search()`` returns the **1-based position of the first base of
the downstream fragment**. Everywhere in this package we store the 0-based index of that
base and call it the *cut index* ``t``: the top strand is broken immediately before ``t``.

``enzyme.ovhg`` is negative for a 5' overhang, positive for a 3' overhang, zero for blunt.
The bottom-strand break sits at ``t - ovhg``, so the single-stranded footprint of the
double-strand break is ``[min(t, b), max(t, b))``.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio.Restriction import AllEnzymes, CommOnly, RestrictionBatch

#: Ligation strategies, most to least preferred. See docs/DECISIONS.md D6 and D13.
COMPATIBLE_OVERHANG = "compatible_overhang"
SAME_ENZYME = "same_enzyme"
BLUNT = "blunt"

STRATEGY_RANK = {COMPATIBLE_OVERHANG: 0, SAME_ENZYME: 1, BLUNT: 2}

#: Preferred vendor when choosing between equischizomers (enzymes with an identical
#: site *and* an identical cut position), so we report e.g. MfeI rather than MunI.
_SUPPLIER_PREFERENCE = (
    "New England Biolabs",
    "Thermo Fisher Scientific",
    "Takara Bio Inc.",
    "Promega",
)


@dataclass(frozen=True)
class CutSite:
    """A single cut of one enzyme on one plasmid."""

    enzyme: object
    top: int  # 0-based cut index on the top strand
    bottom: int  # 0-based cut index on the bottom strand

    @property
    def name(self) -> str:
        return str(self.enzyme)

    @property
    def footprint(self) -> tuple[int, int]:
        """0-based half-open span of single-stranded DNA left by the break."""
        return (min(self.top, self.bottom), max(self.top, self.bottom))


def build_pool(spec="CommOnly") -> RestrictionBatch:
    """Build the candidate enzyme pool.

    ``spec`` is ``"CommOnly"`` (commercially available), ``"All"`` (all of REBASE), or an
    iterable of enzyme names -- the hook for a lab stock list.
    """
    if isinstance(spec, str):
        if spec == "CommOnly":
            return RestrictionBatch(list(CommOnly))
        if spec == "All":
            return RestrictionBatch(list(AllEnzymes))
        raise ValueError(f"Unknown enzyme pool {spec!r}; use 'CommOnly', 'All', or a name list")
    return RestrictionBatch(list(spec))


def enumerate_cut_sites(seq, pool, circular: bool = True) -> dict[str, list[CutSite]]:
    """Map enzyme name -> its cut sites on ``seq``.

    Enzymes with an undefined cut position (``is_unknown()``) are dropped: we cannot
    reason about what they would delete.
    """
    batch = pool if isinstance(pool, RestrictionBatch) else build_pool(pool)
    results: dict[str, list[CutSite]] = {}
    for enzyme in batch:
        if enzyme.is_unknown():
            continue
        try:
            positions = enzyme.search(seq, linear=not circular)
        except Exception:  # pragma: no cover - a few REBASE entries are unsearchable
            continue
        if not positions:
            continue
        results[str(enzyme)] = [
            CutSite(enzyme=enzyme, top=p - 1, bottom=(p - 1) - enzyme.ovhg) for p in positions
        ]
    return results


def ligation_strategy(upstream, downstream) -> str | None:
    """How the two ends would be joined, or ``None`` if they cannot be ligated.

    Note this delegates to ``compatible_end()`` rather than comparing overhang strings:
    KpnI and Acc65I both leave a ``GTAC`` overhang but with opposite polarity (3' vs 5')
    and are *not* compatible.
    """
    if str(upstream) == str(downstream):
        return SAME_ENZYME
    if upstream.is_blunt() and downstream.is_blunt():
        return BLUNT
    if downstream in upstream.compatible_end():
        return COMPATIBLE_OVERHANG
    return None


def canonical_name(enzyme) -> str:
    """Preferred name among equischizomers (identical site *and* cut position).

    Ranked by preferred vendor, then by how many suppliers stock it, then alphabetically.
    The supplier count matters: XhoI and PaeR7I are the same enzyme and both are sold by NEB,
    but XhoI has nine suppliers to PaeR7I's one and is the name everyone writes on a tube.
    """
    candidates = [enzyme, *enzyme.equischizomers()]

    def sort_key(e):
        suppliers = e.supplier_list()
        best = min(
            (_SUPPLIER_PREFERENCE.index(s) for s in suppliers if s in _SUPPLIER_PREFERENCE),
            default=len(_SUPPLIER_PREFERENCE),
        )
        return (best, -len(suppliers), str(e))

    return str(sorted(candidates, key=sort_key)[0])


def supplier_count(enzyme) -> int:
    """How many vendors stock this enzyme (under its canonical name).

    A proxy for "is this the enzyme a lab actually has". Several enzymes can cut at the same
    position -- XhoI, PaeR7I and PspXI all cut ``CTCGAG`` here -- and without this the choice
    falls to alphabetical order, which picks rare enzymes like PspXI over XhoI.
    """
    best = enzyme
    for candidate in [enzyme, *enzyme.equischizomers()]:
        if str(candidate) == canonical_name(enzyme):
            best = candidate
            break
    return len(best.supplier_list())


def describe(enzyme) -> str:
    """Compact human-readable description, e.g. ``MfeI (C^AATT_G, 5' AATT)``."""
    if enzyme.is_blunt():
        end = "blunt"
    else:
        end = f"{5 if enzyme.is_5overhang() else 3}' {enzyme.ovhgseq}"
    return f"{enzyme} ({enzyme.elucidate()}, {end})"
