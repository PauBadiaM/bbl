"""Insert a sequence of interest into an opened backbone.

Two cases, matching how the work is actually done at the bench:

**Case 1 -- the insert comes from another plasmid.** First check whether the same restriction
enzymes can be reused: cut the donor to release the fragment and cut the vector to open it, so
the ends ligate directly. This is classic subcloning and needs no polymerase at all. When it is
not available, PCR the insert out of the donor with primers carrying homology tails and join by
Gibson.

**Case 2 -- you only have the sequence.** Order it with homology arms already attached (or add
them by PCR) and join by Gibson.

The rule is deliberately two-way: **restriction when the sites are already there, Gibson when
they are not.** A hybrid -- PCR-ing the insert with restriction sites in the primer tails, then
digesting the amplicon -- was built and removed (docs/DECISIONS.md D61). It reproduces a
historical construct byte-for-byte, but it costs an extra digest, needs the sites to be absent
from the insert, and needs clamp bases outside each site, for a plasmid that is functionally
identical to the Gibson product.

An important asymmetry in the restriction route: the **vector** must be cut exactly once by each
enzyme or the backbone is destroyed, but the **donor** may be cut anywhere else as well -- you
gel-purify the insert band. Only a cut *inside* the fragment is disqualifying. A symmetric rule
would reject many perfectly good subclonings.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord

from .enzymes import (
    build_pool,
    canonical_name,
    enumerate_cut_sites,
    ligation_strategy,
    supplier_count,
)
from .pcr import Primer, _annotate, design_inward_pair, design_outward_pair
from .plasmid_io import (
    features_in_span,
    load_plasmid,
    replace_span,
    slice_circular,
)
from .targets import classify_features, handle_warnings, resolve_target

RESTRICTION = "restriction"
GIBSON_PCR = "gibson_pcr"
GIBSON_SYNTHETIC = "gibson_synthetic"

_JUNCTION_CONTEXT = 6
_MAX_FLANK_CANDIDATES = 30
_SHORT_INSERT_BP = 200


class NoInsertionRoute(RuntimeError):
    """Raised when no viable way to build the insertion was found."""

    def __init__(self, message: str, reasons: dict[str, int] | None = None):
        super().__init__(message)
        self.reasons = reasons or {}


# ---------------------------------------------------------------------------
# what is being inserted
# ---------------------------------------------------------------------------


@dataclass
class InsertSource:
    """The sequence of interest and where it came from."""

    sequence: str
    name: str
    donor: object | None = None
    span: tuple[int, int] | None = None
    features: list = field(default_factory=list)

    @property
    def synthetic(self) -> bool:
        return self.donor is None

    @property
    def length(self) -> int:
        return len(self.sequence)


def resolve_insert(insert, features=None, name=None) -> InsertSource:
    """Normalise the many ways of naming an insert into an :class:`InsertSource`.

    ``insert`` may be a raw DNA string (case 2), or a donor plasmid -- path, ``SeqRecord`` or
    ``Dseqrecord`` -- together with ``features`` naming the sequence of interest (case 1).
    A donor given without ``features`` is treated as a linear fragment in its entirety.
    """
    if isinstance(insert, InsertSource):
        return insert

    if isinstance(insert, str) and not Path(insert).suffix:
        cleaned = insert.strip().upper().replace(" ", "").replace("\n", "")
        if cleaned and set(cleaned) <= set("ACGTRYSWKMBDHVN"):
            return InsertSource(sequence=cleaned, name=name or "synthetic_insert")
        raise ValueError("insert looks like neither a DNA string nor a file path")

    donor = load_plasmid(insert)
    if features is None:
        return InsertSource(
            sequence=str(donor.seq).upper(),
            name=name or donor.name,
            donor=donor,
            span=(0, len(donor)),
            features=features_in_span(donor, 0, len(donor)),
        )

    donor, span, labels = resolve_target(donor, features)
    start, end = span
    return InsertSource(
        sequence=str(donor.seq)[start:end].upper(),
        name=name or ", ".join(labels),
        donor=donor,
        span=span,
        features=features_in_span(donor, start, end),
    )


# ---------------------------------------------------------------------------
# results
# ---------------------------------------------------------------------------


@dataclass
class Fragment:
    """One physical piece of DNA you will have in a tube."""

    role: str  # "vector" | "insert"
    name: str
    length: int
    preparation: str
    enzymes: tuple[str, ...] | None = None
    primers: tuple[Primer, Primer] | None = None
    order_sequence: str | None = None  # for a synthesised fragment

    def __str__(self) -> str:
        return f"{self.role}: {self.name} ({self.length} bp) -- {self.preparation}"


@dataclass
class InsertionPlan:
    """A complete plan for putting an insert into a backbone."""

    product: object
    strategy: str
    vector: Fragment
    insert: Fragment
    site: tuple[int, int]
    inserted_bp: int
    replaced_bp: int
    junctions: tuple[str, str]
    directional: bool | None = None
    sites_regenerated: bool = False
    warnings: list[str] = field(default_factory=list)

    @property
    def enzymes(self) -> tuple[str, ...]:
        return self.vector.enzymes or ()

    @property
    def protocol(self) -> str:
        lines = []
        if self.strategy == RESTRICTION:
            up, down = self.vector.enzymes
            lines.append(
                f"1. Digest the vector with {up} + {down} and gel-purify the large backbone "
                f"({self.vector.length} bp)."
            )
            lines.append(
                f"2. Digest the donor with {' + '.join(self.insert.enzymes)} and gel-purify the "
                f"{self.insert.length} bp insert band. Extra cuts in the donor backbone are "
                "fine -- they only help separate the band."
            )
            if self.directional:
                lines.append(
                    "3. Ligate (T4). The two vector ends are mutually incompatible, so the "
                    "insert can only go in one way and the vector cannot self-close."
                )
            else:
                lines.append(
                    "3. Dephosphorylate the vector (rSAP/CIP) before ligating -- its ends are "
                    "compatible with each other, so it can self-close, and the insert can go "
                    "in either orientation. Screen for orientation."
                )
        else:
            lines.append(f"1. {self.vector.preparation}")
            lines.append(f"2. {self.insert.preparation}")
            lines.append(
                "3. Assemble vector + insert with Gibson/NEBuilder HiFi (50 C, 15-60 min), "
                "then transform."
            )
        lines.append(
            f"4. Screen across both junctions. The product is {len(self.product)} bp "
            f"({self.inserted_bp} bp in, {self.replaced_bp} bp out)."
        )
        return "\n".join(lines)

    def summary(self) -> str:
        return (
            f"{self.strategy}: insert {self.insert.name} ({self.inserted_bp} bp) at "
            f"{self.site}, replacing {self.replaced_bp} bp -> {len(self.product)} bp product"
        )


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _resolve_site(vector, at):
    """Return ``(vector, (start, end))`` for the place the insert goes."""
    if isinstance(at, int):
        return vector, (at, at)
    if isinstance(at, (tuple, list)) and len(at) == 2:
        return vector, (int(at[0]), int(at[1]))
    vector, span, _ = resolve_target(vector, [at] if isinstance(at, str) else at)
    return vector, span


def _overlaps(span, protected):
    start, end = span
    return [
        label for (p_start, p_end, label) in protected if p_start < end and start < p_end
    ]


def _unique_cuts(record, pool):
    return {
        name: sites[0]
        for name, sites in enumerate_cut_sites(record.seq, pool).items()
        if len(sites) == 1
    }


def _flanking(cuts, boundary, side, limit=_MAX_FLANK_CANDIDATES):
    """Cut sites on one side of ``boundary``, nearest first."""
    if side == "upstream":
        chosen = [c for c in cuts if c.top <= boundary]
        chosen.sort(key=lambda c: -c.top)
    else:
        chosen = [c for c in cuts if c.top >= boundary]
        chosen.sort(key=lambda c: c.top)
    return chosen[:limit]


def _verify_product(vector, product, protected, expected_length, wanted_sequence):
    """Return a list of failures; empty means the assembled product is sound.

    Mirrors :func:`bbl.excise._verify`. Insertion previously checked length, insert presence,
    homology-arm uniqueness and site regeneration -- but never re-checked that the features the
    caller asked to protect actually survived. ``protected`` steered the choice of cut sites and
    was then dropped, so a route that destroyed a protected feature by some path other than
    cutting inside it would have been returned without complaint. Insertion is the more
    error-prone of the two operations, so it should not verify to a weaker standard than
    excision.

    Searches the **doubled** product so a feature spanning the origin still reads as intact.
    """
    failures = []
    if len(product) != expected_length:
        failures.append(f"length {len(product)} != expected {expected_length}")

    doubled = (str(product.seq) * 2).upper()
    if wanted_sequence and wanted_sequence.upper() not in doubled:
        failures.append("the sequence of interest is missing from the product")

    for start, end, label in protected:
        start, end = max(0, start), min(len(vector), end)
        original = str(vector.seq[start:end]).upper()
        if original and original not in doubled:
            failures.append(f"protected feature {label!r} is not intact in the product")
    return failures


def _junctions(vector_text, left_cut, right_cut, insert_sequence):
    upstream = (
        slice_circular(vector_text, left_cut - _JUNCTION_CONTEXT, left_cut)
        + insert_sequence[:_JUNCTION_CONTEXT]
    )
    downstream = (
        insert_sequence[-_JUNCTION_CONTEXT:]
        + slice_circular(vector_text, right_cut, right_cut + _JUNCTION_CONTEXT)
    )
    return upstream.upper(), downstream.upper()


# ---------------------------------------------------------------------------
# route 1: reuse restriction enzymes
# ---------------------------------------------------------------------------


def _restriction_route(vector, site, source, protected, pool, note):
    """Classic subcloning: cut donor and vector so the ends ligate directly."""
    if source.donor is None:
        note("insert is synthetic, so there is nothing to digest")
        return None

    site_start, site_end = site
    soi_start, soi_end = source.span
    vector_cuts = _unique_cuts(vector, pool)
    donor_all = enumerate_cut_sites(source.donor.seq, pool)

    upstream_vector = _flanking(list(vector_cuts.values()), site_start, "upstream")
    downstream_vector = _flanking(list(vector_cuts.values()), site_end, "downstream")
    if not upstream_vector:
        note("no unique vector site upstream of the insertion point")
    if not downstream_vector:
        note("no unique vector site downstream of the insertion point")

    compatible_cache: dict[str, list[str]] = {}

    def partners(enzyme):
        key = str(enzyme)
        if key not in compatible_cache:
            compatible_cache[key] = [key] + [
                str(other) for other in enzyme.compatible_end() if str(other) != key
            ]
        return compatible_cache[key]

    def donor_options(enzyme, boundary, side, limit=6):
        """Usable donor cuts on one side, best first.

        Ordered by: the *same* enzyme before a merely compatible one, then fewer cuts in the
        donor overall (a promiscuous partner shreds the gel), then proximity to the insert.
        Taking simply the closest cut is wrong -- it reaches for frequent cutters such as
        BfaI (C^TAG, compatible with XbaI) that then cut inside the fragment.
        """
        found = []
        for index, name in enumerate(partners(enzyme)):
            cuts = donor_all.get(name, [])
            for cut in cuts:
                if side == "upstream" and cut.top <= boundary:
                    distance = boundary - cut.top
                elif side == "downstream" and cut.top >= boundary:
                    distance = cut.top - boundary
                else:
                    continue
                found.append(((0 if index == 0 else 1, len(cuts), distance), cut))
        found.sort(key=lambda item: item[0])
        return [cut for _, cut in found[:limit]]

    def cuts_inside(up_cut, down_cut):
        return any(
            up_cut.top < cut.top < down_cut.top
            for name in {str(up_cut.enzyme), str(down_cut.enzyme)}
            for cut in donor_all.get(name, [])
        )

    candidates = []
    for up_vector in upstream_vector:
        for down_vector in downstream_vector:
            if str(up_vector.enzyme) == str(down_vector.enzyme):
                continue  # one enzyme cannot cut the vector uniquely at two places
            damaged = _overlaps((up_vector.top, down_vector.top), protected)
            if damaged:
                note(f"opening the vector would damage {damaged[0]}")
                continue

            # The two enzymes must not cut *inside* the released fragment. Cuts elsewhere in
            # the donor are fine -- the band is gel-purified. Try each candidate partner
            # rather than committing to the nearest one and failing.
            up_donor = down_donor = None
            up_options = donor_options(up_vector.enzyme, soi_start, "upstream")
            down_options = donor_options(down_vector.enzyme, soi_end, "downstream")
            if not up_options or not down_options:
                note("donor lacks a compatible flanking site")
                continue
            for up_candidate in up_options:
                for down_candidate in down_options:
                    if up_candidate.top >= down_candidate.top:
                        continue
                    if cuts_inside(up_candidate, down_candidate):
                        continue
                    up_donor, down_donor = up_candidate, down_candidate
                    break
                if up_donor is not None:
                    break
            if up_donor is None:
                note("every candidate enzyme cuts inside the insert fragment")
                continue

            directional = ligation_strategy(up_vector.enzyme, down_vector.enzyme) is None
            same_enzymes = (
                str(up_donor.enzyme) == str(up_vector.enzyme)
                and str(down_donor.enzyme) == str(down_vector.enzyme)
            )
            candidates.append(
                {
                    "up_vector": up_vector,
                    "down_vector": down_vector,
                    "up_donor": up_donor,
                    "down_donor": down_donor,
                    "directional": directional,
                    "score": (
                        0 if directional else 1,
                        0 if same_enzymes else 1,
                        (down_vector.top - up_vector.top) - (site_end - site_start),
                        (down_donor.top - up_donor.top) - (soi_end - soi_start),
                        # several enzymes can cut at the same position; prefer the one a lab
                        # is likely to own rather than whichever sorts first alphabetically
                        -supplier_count(up_vector.enzyme) - supplier_count(down_vector.enzyme),
                        canonical_name(up_vector.enzyme),
                        canonical_name(down_vector.enzyme),
                    ),
                }
            )

    if not candidates:
        return None
    return min(candidates, key=lambda c: c["score"])


# ---------------------------------------------------------------------------
# route 2: Gibson
# ---------------------------------------------------------------------------


def _vector_digest_pair(vector, site, protected, pool, max_collateral):
    """A pair of unique cutters bracketing the site, for opening the vector by digest."""
    site_start, site_end = site
    cuts = _unique_cuts(vector, pool)
    best = None
    for up in _flanking(list(cuts.values()), site_start, "upstream"):
        for down in _flanking(list(cuts.values()), site_end, "downstream"):
            if str(up.enzyme) == str(down.enzyme):
                continue
            removed = down.top - up.top
            if removed - (site_end - site_start) > max_collateral:
                continue
            if _overlaps((up.top, down.top), protected):
                continue
            rank = (
                removed,
                -supplier_count(up.enzyme) - supplier_count(down.enzyme),
                canonical_name(up.enzyme),
            )
            if best is None or rank < best[0]:
                best = (rank, up, down)
    return (best[1], best[2]) if best else None


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------


def plan_insertion(
    backbone,
    insert,
    at,
    insert_features=None,
    protect=None,
    method: str = "auto",
    enzyme_pool="CommOnly",
    homology: int = 25,
    vector_prep: str = "auto",
    max_digest_collateral: int = 200,
    target_tm: float = 60.0,
    min_length: int = 18,
    max_length: int = 32,
    min_gc: float = 35.0,
    max_gc: float = 65.0,
) -> InsertionPlan:
    """Plan an insertion of ``insert`` into ``backbone`` at ``at``.

    ``at`` is a feature label (the feature is replaced), an integer position (pure insertion),
    or a ``(start, end)`` span to replace. ``insert`` is a DNA string or a donor plasmid plus
    ``insert_features``.

    ``method`` is ``"auto"`` (reuse restriction enzymes if possible, else Gibson),
    ``"restriction"``, or ``"gibson"``. ``vector_prep`` is ``"auto"``, ``"digest"`` or
    ``"pcr"`` and only applies to the Gibson route; ``max_digest_collateral`` caps how much
    extra vector sequence a digest may remove before inverse PCR is used instead.
    """
    if method not in ("auto", RESTRICTION, "gibson"):
        raise ValueError(f"method must be 'auto', {RESTRICTION!r} or 'gibson', got {method!r}")

    vector = load_plasmid(backbone)
    vector, site = _resolve_site(vector, at)
    site_start, site_end = site
    source = resolve_insert(insert, insert_features)
    pool = build_pool(enzyme_pool)

    _, protected, umbrella = classify_features(vector, site, protect)
    reasons: dict[str, int] = {}

    def note(reason):
        reasons[reason] = reasons.get(reason, 0) + 1

    warnings = [
        f"{label!r} spans the insertion site and will be interrupted" for label in umbrella
    ]
    vector_text = str(vector.seq).upper()

    # ---- route 1: reuse restriction enzymes -------------------------------------------
    if method in ("auto", RESTRICTION):
        chosen = _restriction_route(vector, site, source, protected, pool, note)
        if chosen is not None:
            up_v, down_v = chosen["up_vector"], chosen["down_vector"]
            up_d, down_d = chosen["up_donor"], chosen["down_donor"]
            donor_text = str(source.donor.seq).upper()
            fragment = donor_text[up_d.top : down_d.top]

            product = replace_span(
                vector,
                up_v.top,
                down_v.top,
                fragment,
                features=features_in_span(source.donor, up_d.top, down_d.top),
                truncation_note="truncated by insertion",
            )
            failures = _verify_product(
                vector,
                product,
                protected,
                len(vector) - (down_v.top - up_v.top) + len(fragment),
                source.sequence,
            )
            if failures:
                raise NoInsertionRoute("; ".join(failures))
            warnings.extend(handle_warnings(vector, product))

            enzymes_used = {up_v.enzyme, down_v.enzyme, up_d.enzyme, down_d.enzyme}
            regenerated = any(e.search(product.seq, linear=False) for e in enzymes_used)

            if not chosen["directional"]:
                warnings.append(
                    "the vector's two ends are compatible with each other: it can self-close "
                    "and the insert can invert -- dephosphorylate and screen orientation"
                )
            extra_donor = (down_d.top - up_d.top) - source.length
            if extra_donor:
                warnings.append(
                    f"{extra_donor} bp of donor sequence flanking the insert comes along"
                )
            if regenerated:
                warnings.append("a recognition site is regenerated at a junction")

            return InsertionPlan(
                product=product,
                strategy=RESTRICTION,
                vector=Fragment(
                    role="vector",
                    name=f"{vector.name} opened",
                    length=len(vector) - (down_v.top - up_v.top),
                    preparation=(
                        f"digest with {canonical_name(up_v.enzyme)} + "
                        f"{canonical_name(down_v.enzyme)}, gel-purify the backbone"
                    ),
                    enzymes=(canonical_name(up_v.enzyme), canonical_name(down_v.enzyme)),
                ),
                insert=Fragment(
                    role="insert",
                    name=source.name,
                    length=len(fragment),
                    preparation=(
                        f"digest {source.donor.name} with {canonical_name(up_d.enzyme)} + "
                        f"{canonical_name(down_d.enzyme)}, gel-purify the "
                        f"{len(fragment)} bp band"
                    ),
                    enzymes=(canonical_name(up_d.enzyme), canonical_name(down_d.enzyme)),
                ),
                site=(up_v.top, down_v.top),
                inserted_bp=len(fragment),
                replaced_bp=down_v.top - up_v.top,
                junctions=_junctions(vector_text, up_v.top, down_v.top, fragment),
                directional=chosen["directional"],
                sites_regenerated=bool(regenerated),
                warnings=warnings,
            )
        if method == RESTRICTION:
            detail = "; ".join(f"{r} (x{c})" for r, c in sorted(reasons.items(), key=lambda kv: -kv[1])[:4])
            raise NoInsertionRoute(
                f"the restriction enzymes cannot be reused for this insert. {detail}. "
                "Use method='gibson' to amplify the insert instead.",
                reasons,
            )

    # ---- route 2: Gibson ---------------------------------------------------------------
    digest_pair = None
    if vector_prep in ("auto", "digest"):
        digest_pair = _vector_digest_pair(
            vector, site, protected, pool, max_digest_collateral
        )
        if digest_pair is None and vector_prep == "digest":
            raise NoInsertionRoute(
                "no pair of unique cutters brackets the insertion site; use vector_prep='pcr'"
            )

    if digest_pair is not None:
        up_v, down_v = digest_pair
        open_start, open_end = up_v.top, down_v.top
        vector_fragment = Fragment(
            role="vector",
            name=f"{vector.name} opened",
            length=len(vector) - (open_end - open_start),
            preparation=(
                f"Digest the vector with {canonical_name(up_v.enzyme)} + "
                f"{canonical_name(down_v.enzyme)} and gel-purify the backbone."
            ),
            enzymes=(canonical_name(up_v.enzyme), canonical_name(down_v.enzyme)),
        )
    else:
        open_start, open_end = site_start, site_end
        forward_option, reverse_option = design_outward_pair(
            vector_text,
            (open_start, open_end),
            min_length=min_length,
            max_length=max_length,
            target_tm=target_tm,
            min_gc=min_gc,
            max_gc=max_gc,
        )
        vector_forward = _annotate(forward_option, vector_text, "vec_F", +1, "", False)
        vector_reverse = _annotate(reverse_option, vector_text, "vec_R", -1, "", False)
        amplicon = len(vector) - (open_end - open_start)
        vector_fragment = Fragment(
            role="vector",
            name=f"{vector.name} opened",
            length=amplicon,
            preparation=(
                f"Inverse-PCR the vector with {vector_forward.name} + {vector_reverse.name} "
                f"({amplicon} bp), then DpnI-digest the template."
            ),
            primers=(vector_forward, vector_reverse),
        )
        for primer in (vector_forward, vector_reverse):
            warnings.extend(f"{primer.name}: {n}" for n in primer.notes)

    arm_up = slice_circular(vector_text, open_start - homology, open_start)
    arm_down = slice_circular(vector_text, open_end, open_end + homology)

    if source.donor is None:
        strategy = GIBSON_SYNTHETIC
        order = arm_up + source.sequence + arm_down
        insert_fragment = Fragment(
            role="insert",
            name=source.name,
            length=len(order),
            preparation=(
                f"Order {source.name} as a synthetic fragment with the homology arms already "
                f"attached ({homology} bp each side, {len(order)} bp total)."
            ),
            order_sequence=order,
        )
        if len(order) < 125:
            warnings.append(
                f"{len(order)} bp is below the usual synthesis minimum; annealed oligos or "
                "carrying the insert in primer tails may be cheaper"
            )
    else:
        strategy = GIBSON_PCR
        donor_text = str(source.donor.seq).upper()
        forward_option, reverse_option = design_inward_pair(
            donor_text,
            source.span,
            min_length=min_length,
            max_length=max_length,
            target_tm=target_tm,
            min_gc=min_gc,
            max_gc=max_gc,
        )
        insert_forward = _annotate(forward_option, donor_text, "ins_F", +1, arm_up, False)
        insert_reverse = _annotate(
            reverse_option,
            donor_text,
            "ins_R",
            -1,
            str(Seq(arm_down).reverse_complement()),
            False,
        )
        insert_fragment = Fragment(
            role="insert",
            name=source.name,
            length=source.length + 2 * homology,
            preparation=(
                f"PCR {source.name} out of {source.donor.name} with {insert_forward.name} + "
                f"{insert_reverse.name}; their 5' tails add {homology} bp of homology to each "
                "vector end. DpnI-digest the template."
            ),
            primers=(insert_forward, insert_reverse),
        )
        for primer in (insert_forward, insert_reverse):
            warnings.extend(f"{primer.name}: {n}" for n in primer.notes)

    product = replace_span(
        vector,
        open_start,
        open_end,
        source.sequence,
        features=source.features,
        truncation_note="truncated by insertion",
    )
    failures = _verify_product(
        vector,
        product,
        protected,
        len(vector) - (open_end - open_start) + source.length,
        source.sequence,
    )
    if failures:
        raise NoInsertionRoute("; ".join(failures))
    warnings.extend(handle_warnings(vector, product))

    product_text = str(product.seq).upper()
    for label, arm in (("upstream", arm_up), ("downstream", arm_down)):
        if product_text.count(arm) > 1:
            warnings.append(
                f"the {label} homology arm occurs more than once in the product -- "
                "assembly may be ambiguous; move the junction or lengthen the arm"
            )
    if source.length < _SHORT_INSERT_BP:
        warnings.append(
            f"a {source.length} bp insert assembles inefficiently by Gibson; consider "
            "carrying it entirely in primer tails, or annealed oligos if shorter still"
        )
    warnings.append(
        "Tm is nearest-neighbour at default salt and dimer screening is heuristic -- check "
        "the oligos in Primer3 or IDT OligoAnalyzer before ordering."
    )

    product.name = f"{vector.name}_ins"[:16]
    product.description = f"{vector.description} | {source.name} inserted at {site_start}"

    return InsertionPlan(
        product=product,
        strategy=strategy,
        vector=vector_fragment,
        insert=insert_fragment,
        site=(open_start, open_end),
        inserted_bp=source.length,
        replaced_bp=open_end - open_start,
        junctions=_junctions(vector_text, open_start, open_end, source.sequence),
        directional=True,  # homology arms are unique to each end
        warnings=warnings,
    )


def insert_sequence(backbone, insert, at, **kwargs) -> InsertionPlan:
    """Convenience alias for :func:`plan_insertion`."""
    return plan_insertion(backbone, insert, at, **kwargs)
