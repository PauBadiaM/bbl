"""Delete annotated features from a plasmid by restriction digest and religation.

Given a plasmid and the features to remove, find flanking restriction sites whose ends can
actually be ligated, excise the DNA between the two cuts, and re-circularize.

The design rationale -- including several places where the obvious approach is
biologically wrong -- is recorded in ``docs/DECISIONS.md``.
"""

from __future__ import annotations

from dataclasses import dataclass, field


from .enzymes import (
    BLUNT,
    COMPATIBLE_OVERHANG,
    SAME_ENZYME,
    STRATEGY_RANK,
    CutSite,
    build_pool,
    canonical_name,
    describe,
    enumerate_cut_sites,
    ligation_strategy,
)
from .plasmid_io import circular_equal, delete_span, load_plasmid, slice_circular
from .targets import classify_features, essential_warnings, resolve_target

_JUNCTION_CONTEXT = 6


class NoExcisionFound(RuntimeError):
    """Raised when no ligatable pair of flanking sites exists.

    ``reasons`` explains what was rejected, so the agent layer can fall back to a
    PCR-based route (around-the-horn / Gibson) instead of guessing.
    """

    def __init__(self, message: str, reasons: dict[str, int] | None = None):
        super().__init__(message)
        self.reasons = reasons or {}


@dataclass
class ExcisionPlan:
    """One way to delete the requested features."""

    enzymes: tuple[str, ...]
    strategy: str
    deleted_span: tuple[int, int]
    deleted_bp: int
    collateral_bp: int
    removed_features: list[str]
    junction_seq: str
    sites_regenerated: bool
    product: object | None = None
    warnings: list[str] = field(default_factory=list)
    alternatives: list["ExcisionPlan"] = field(default_factory=list)

    @property
    def enzyme_pair(self) -> tuple[str, str]:
        """The two enzymes to order, as a pair (both entries equal for a single digest)."""
        return self.enzymes if len(self.enzymes) == 2 else (self.enzymes[0], self.enzymes[0])

    @property
    def protocol(self) -> str:
        lines = []
        if self.strategy == SAME_ENZYME:
            lines.append(f"1. Digest with {self.enzymes[0]} (cuts twice, flanking the target).")
            lines.append(
                "2. Gel-purify the large backbone fragment away from the "
                f"{self.deleted_bp + 4} bp excised fragment -- the site is regenerated on "
                "ligation, so the insert can go back in. Dephosphorylate (rSAP/CIP) to "
                "suppress re-insertion."
            )
        else:
            lines.append(f"1. Double-digest with {self.enzymes[0]} + {self.enzymes[1]}.")
            lines.append("2. Gel-purify the large backbone fragment.")
        lines.append(
            f"3. Ligate (T4 DNA ligase){' -- blunt ligation, expect low efficiency' if self.strategy == BLUNT else ''}"
            f"; the junction becomes {self.junction_seq}."
        )
        lines.append(
            f"4. Verify: product is {len(self.product) if self.product is not None else '?'} bp, "
            f"{self.deleted_bp} bp shorter than the parent"
            + (
                "; neither site is regenerated, so a diagnostic digest is negative."
                if not self.sites_regenerated
                else "; the site is regenerated."
            )
        )
        return "\n".join(lines)

    def summary(self) -> str:
        pair = " + ".join(self.enzymes)
        return (
            f"{pair} [{self.strategy}] deletes {self.deleted_bp} bp "
            f"({self.collateral_bp} bp collateral) at {self.deleted_span}, "
            f"removing {', '.join(self.removed_features) or 'nothing annotated'}"
        )


# ---------------------------------------------------------------------------
# target / protected resolution
# ---------------------------------------------------------------------------


def _overlaps_protected(deleted, protected):
    start, end = deleted
    return [
        label
        for (p_start, p_end, label) in protected
        if p_start < end and start < p_end
    ]


# ---------------------------------------------------------------------------
# candidate generation
# ---------------------------------------------------------------------------


def _candidate_pairs(cut_sites, span, protected, allow_single_enzyme):
    """Ranked ``(upstream, downstream, strategy)`` triples plus rejection tallies."""
    target_start, target_end = span
    reasons: dict[str, int] = {}

    def note(reason):
        reasons[reason] = reasons.get(reason, 0) + 1

    single, doubles = {}, {}
    for name, sites in cut_sites.items():
        if len(sites) == 1:
            single[name] = sites[0]
        elif len(sites) == 2:
            doubles[name] = sites
        else:
            note("enzyme cuts more than twice")

    upstream = [s for s in single.values() if s.top <= target_start]
    downstream = [s for s in single.values() if s.top >= target_end]
    if not upstream:
        note("no unique site upstream of the target")
    if not downstream:
        note("no unique site downstream of the target")

    candidates = []
    for u in upstream:
        for d in downstream:
            strategy = ligation_strategy(u.enzyme, d.enzyme)
            if strategy is None:
                note("incompatible ends")
                continue
            damaged = _overlaps_protected((u.top, d.top), protected)
            if damaged:
                note(f"would damage {damaged[0]}")
                continue
            candidates.append((u, d, strategy))

    if allow_single_enzyme:
        for name, sites in doubles.items():
            first, second = sorted(sites, key=lambda s: s.top)
            if not (first.top <= target_start and second.top >= target_end):
                note("enzyme cuts twice but not flanking the target")
                continue
            damaged = _overlaps_protected((first.top, second.top), protected)
            if damaged:
                note(f"would damage {damaged[0]}")
                continue
            candidates.append((first, second, SAME_ENZYME))

    target_bp = target_end - target_start
    candidates.sort(
        key=lambda c: (
            STRATEGY_RANK[c[2]],
            (c[1].top - c[0].top) - target_bp,  # collateral bp
            canonical_name(c[0].enzyme),
            canonical_name(c[1].enzyme),
        )
    )

    # Equischizomers collapse to the same canonical name (MunI -> MfeI), which would
    # otherwise surface the winning plan again as its own "alternative".
    deduped, seen = [], set()
    for upstream, downstream, strategy in candidates:
        key = (
            frozenset({canonical_name(upstream.enzyme), canonical_name(downstream.enzyme)}),
            upstream.top,
            downstream.top,
        )
        if key in seen:
            note("duplicate of an equischizomer")
            continue
        seen.add(key)
        deduped.append((upstream, downstream, strategy))
    return deduped, reasons


# ---------------------------------------------------------------------------
# simulation and verification
# ---------------------------------------------------------------------------


def _signature(record, protected):
    """A distinctive stretch of retained sequence, used to pick the right fragment."""
    if not protected:
        return None
    start, end, _ = max(protected, key=lambda p: p[1] - p[0])
    start, end = max(0, start), min(len(record), end)
    middle = (start + end) // 2
    return str(record.seq[max(0, middle - 30) : middle + 30]).upper()


def _simulate(record, upstream, downstream, strategy, signature):
    """Cut, keep the arc holding ``signature``, and re-circularize. Returns the product."""
    enzymes = (
        [upstream.enzyme]
        if strategy == SAME_ENZYME
        else [upstream.enzyme, downstream.enzyme]
    )
    fragments = record.cut(*enzymes)
    if len(fragments) != 2:
        raise ValueError(f"expected 2 fragments, got {len(fragments)}")

    if signature is None:
        keep = max(fragments, key=len)
    else:
        holding = [f for f in fragments if signature in f.seq.watson.upper()]
        if len(holding) != 1:
            raise ValueError("could not identify the fragment to retain")
        keep = holding[0]
    return keep.looped()


def _reannotate(record, upstream_cut, downstream_cut, simulated):
    """Rebuild the product in the *parent's* coordinate frame, keeping annotations.

    ``pydna``'s ``cut()`` is the oracle for whether the ends ligate and what the sequence is,
    but it discards features that straddle a cut and returns the product rotated to the
    downstream cut. ``delete_span()`` keeps the parent's origin and truncates rather than drops
    straddling features; the result is cross-checked against ``simulated``.
    """
    product = delete_span(record, upstream_cut, downstream_cut)
    if not circular_equal(product, simulated):
        raise ValueError("re-annotated product disagrees with the digest simulation")
    return product


_VERIFY_CONTEXT = 20


def _verify(record, product, deleted_bp, protected, upstream_cut, downstream_cut):
    """Return a list of failures; empty means the product is sound."""
    failures = []
    if len(product) != len(record) - deleted_bp:
        failures.append(f"length {len(product)} != expected {len(record) - deleted_bp}")

    doubled = (product.seq.watson * 2).upper()
    for start, end, label in protected:
        start, end = max(0, start), min(len(record), end)
        original = str(record.seq[start:end]).upper()
        if original and original not in doubled:
            failures.append(f"protected feature {label!r} is not intact in the product")

    # Confirm the deletion actually happened by checking that the parent's sequence
    # *across each cut* is no longer contiguous. Searching for the target sequence itself
    # would be wrong: a plasmid may legitimately carry duplicates of it elsewhere (e.g. one
    # repeat of a tandem aptamer array), which would make a correct product look broken.
    tripled = (str(record.seq) * 3).upper()
    length = len(record)
    for cut in (upstream_cut, downstream_cut):
        context = tripled[length + cut - _VERIFY_CONTEXT : length + cut + _VERIFY_CONTEXT]
        if len(context) == 2 * _VERIFY_CONTEXT and context in doubled:
            failures.append(f"parent sequence across the cut at {cut} is still contiguous")
    return failures


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------


def plan_excisions(
    plasmid,
    features_to_remove,
    protect=None,
    enzyme_pool="CommOnly",
    margin=0,
    allow_single_enzyme=True,
    max_alternatives=10,
):
    """Ranked, simulation-verified excision plans. Returns ``(plans, reasons)``.

    ``plans`` may be empty; ``reasons`` then tallies why candidates were rejected.
    """
    record, span, target_labels = resolve_target(load_plasmid(plasmid), features_to_remove)
    removed_labels, protected, umbrella = classify_features(record, span, protect, margin)

    base_warnings = []
    for label in umbrella:
        base_warnings.append(
            f"{label!r} spans the deleted region and will be shortened in the product"
        )
    base_warnings.extend(essential_warnings(removed_labels))

    cut_sites = enumerate_cut_sites(record.seq, build_pool(enzyme_pool), circular=True)
    candidates, reasons = _candidate_pairs(cut_sites, span, protected, allow_single_enzyme)

    signature = _signature(record, protected)
    plans: list[ExcisionPlan] = []
    budget = max_alternatives + 15

    for upstream, downstream, strategy in candidates[:budget]:
        deleted_span = (upstream.top, downstream.top)
        deleted_bp = downstream.top - upstream.top
        try:
            simulated = _simulate(record, upstream, downstream, strategy, signature)
            product = _reannotate(record, upstream.top, downstream.top, simulated)
        except Exception as exc:  # candidate is unbuildable -> downgrade, keep going
            reasons[f"simulation failed ({exc})"] = reasons.get(f"simulation failed ({exc})", 0) + 1
            continue

        failures = _verify(
            record, product, deleted_bp, protected, upstream.top, downstream.top
        )
        if failures:
            for failure in failures:
                reasons[failure] = reasons.get(failure, 0) + 1
            continue

        enzymes = (
            (canonical_name(upstream.enzyme),)
            if strategy == SAME_ENZYME
            else (canonical_name(upstream.enzyme), canonical_name(downstream.enzyme))
        )
        used = {upstream.enzyme, downstream.enzyme}
        regenerated = any(e.search(product.seq, linear=False) for e in used)

        # slice_circular works on the plain string: a circular Dseq slice with start == stop
        # returns the whole circle, which would corrupt the junction at position 0.
        junction = (
            slice_circular(record, upstream.top - _JUNCTION_CONTEXT, upstream.top)
            + slice_circular(record, downstream.top, downstream.top + _JUNCTION_CONTEXT)
        ).upper()

        warnings = list(base_warnings)
        if strategy == BLUNT:
            warnings.append("blunt ligation is inefficient and non-directional")
        if regenerated and strategy != SAME_ENZYME:
            warnings.append("a recognition site survives in the product")
        warnings.append(
            "methylation sensitivity (Dam/Dcm), star activity and double-digest buffer "
            "compatibility are not modelled -- check the supplier's tables"
        )

        product.name = f"{record.name}_delta"[:16]
        product.description = (
            f"{record.description} | {'+'.join(enzymes)} deletion of "
            f"{', '.join(target_labels)}"
        )

        plans.append(
            ExcisionPlan(
                enzymes=enzymes,
                strategy=strategy,
                deleted_span=deleted_span,
                deleted_bp=deleted_bp,
                collateral_bp=deleted_bp - (span[1] - span[0]),
                removed_features=removed_labels,
                junction_seq=junction,
                sites_regenerated=bool(regenerated),
                product=product,
                warnings=warnings,
            )
        )
        if len(plans) > max_alternatives:
            break

    return plans, reasons


def excise_features(
    plasmid,
    features_to_remove,
    protect=None,
    enzyme_pool="CommOnly",
    margin=0,
    allow_single_enzyme=True,
    max_alternatives=10,
):
    """Delete ``features_to_remove`` from ``plasmid`` using two restriction enzymes.

    Returns the best :class:`ExcisionPlan`: ``plan.product`` is the modified circular
    plasmid and ``plan.enzyme_pair`` the two enzymes used. Runners-up are in
    ``plan.alternatives``.

    Raises :class:`NoExcisionFound` if no ligatable flanking pair exists -- restriction
    excision simply is not always possible, and the caller should then fall back to a
    PCR-based deletion.
    """
    plans, reasons = plan_excisions(
        plasmid,
        features_to_remove,
        protect=protect,
        enzyme_pool=enzyme_pool,
        margin=margin,
        allow_single_enzyme=allow_single_enzyme,
        max_alternatives=max_alternatives,
    )
    if not plans:
        top = sorted(reasons.items(), key=lambda kv: -kv[1])[:5]
        detail = "; ".join(f"{reason} (x{count})" for reason, count in top) or "no candidates"
        raise NoExcisionFound(
            f"No ligatable pair of flanking restriction sites found. Rejections: {detail}. "
            "Consider a PCR-based deletion (around-the-horn or Gibson) instead.",
            reasons,
        )

    best, rest = plans[0], plans[1:]
    for plan in rest:
        plan.product = None  # keep the return value light
    best.alternatives = rest
    return best
