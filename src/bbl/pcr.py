"""Around-the-horn (inverse) PCR deletion.

When no ligatable pair of restriction sites flanks the target, amplify the *backbone*
outward from the target instead: two primers sit back-to-back at the deletion boundaries and
point away from each other, so the target is simply never copied. The linear amplicon is then
re-circularized, either by blunt self-ligation (KLD) or by Gibson assembly, where 5' primer
tails give the amplicon a terminal direct repeat spanning the new junction.

Why this is worth having even though restriction is preferred: the deletion is **exact**. The
boundaries are set by where the primers sit, not by where an enzyme happens to cut, so there
is no collateral loss and no dependence on convenient unique sites.

The 5' end of each primer is pinned to a deletion boundary -- that is what makes the junction
exact -- so the only free parameter per primer is its length, which moves the 3' end.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

from .plasmid_io import delete_span, load_plasmid, slice_circular
from .targets import classify_features, essential_warnings, resolve_target

KLD = "KLD"
GIBSON = "gibson"

_UNIQUENESS_PROBE = 15


class NoPrimerDesign(RuntimeError):
    """Raised when no acceptable primer pair can be built for the requested deletion."""


@dataclass
class Primer:
    """One oligo. ``sequence`` is what you order, 5'->3'."""

    name: str
    sequence: str
    anneal: str  # the template-complementary part (``sequence`` minus any 5' tail)
    tail: str
    tm: float
    gc: float
    binding_span: tuple[int, int]
    strand: int  # +1 forward, -1 reverse
    phosphorylated: bool
    notes: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    def __str__(self) -> str:
        tag = " [5'-phos]" if self.phosphorylated else ""
        return (
            f"{self.name}{tag}: 5'-{self.sequence}-3'  "
            f"({self.length} nt, Tm {self.tm:.1f} C, GC {self.gc:.0f}%)"
        )


@dataclass
class DeletionPCR:
    """A complete around-the-horn deletion design."""

    forward: Primer
    reverse: Primer
    product: object
    amplicon_bp: int
    deleted_bp: int
    deleted_span: tuple[int, int]
    removed_features: list[str]
    method: str
    warnings: list[str] = field(default_factory=list)

    @property
    def primer_pair(self) -> tuple[str, str]:
        return (self.forward.sequence, self.reverse.sequence)

    @property
    def homology_bp(self) -> int:
        """Total end homology carried by the primer tails (0 for a blunt KLD design)."""
        return len(self.forward.tail) + len(self.reverse.tail)

    @property
    def tm_difference(self) -> float:
        return abs(self.forward.tm - self.reverse.tm)

    @property
    def protocol(self) -> str:
        extension = max(15, round(self.amplicon_bp / 1000 * 30))
        anneal_tm = min(self.forward.tm, self.reverse.tm)
        lines = [
            f"1. Amplify the backbone from ~1 ng of parent template with a high-fidelity "
            f"polymerase (Q5/Phusion): {self.forward.name} + {self.reverse.name}, "
            f"annealing ~{anneal_tm - 1:.0f} C, extension ~{extension} s "
            f"({self.amplicon_bp} bp product).",
            "2. DpnI-digest the reaction to destroy the parental plasmid (it is Dam-methylated "
            "from a dam+ host; the PCR product is not).",
        ]
        if self.method == KLD:
            lines.append(
                "3. Phosphorylate and blunt self-ligate: order both primers 5'-phosphorylated, "
                "or treat with T4 PNK. A KLD (kinase/ligase/DpnI) mix does steps 2-3 in one."
            )
        else:
            lines.append(
                f"3. Circularize by Gibson/NEBuilder HiFi assembly -- the 5' tails "
                f"({len(self.forward.tail)} bp on {self.forward.name}, "
                f"{len(self.reverse.tail)} bp on {self.reverse.name}) give the amplicon a "
                f"{self.homology_bp} bp terminal direct repeat spanning the new junction. No "
                "phosphorylation needed."
            )
        lines.append(
            f"4. Transform, then screen across the new junction. The product is "
            f"{len(self.product)} bp, exactly {self.deleted_bp} bp shorter than the parent, "
            "with no collateral sequence loss."
        )
        return "\n".join(lines)

    def summary(self) -> str:
        return (
            f"around-the-horn {self.method}: deletes {self.deleted_bp} bp at "
            f"{self.deleted_span} ({', '.join(self.removed_features) or 'unannotated'}), "
            f"{self.amplicon_bp} bp amplicon, dTm {self.tm_difference:.1f} C"
        )


# ---------------------------------------------------------------------------
# oligo quality heuristics
# ---------------------------------------------------------------------------
#
# These are deliberately simple screens, good enough to reject obviously bad oligos and to
# choose between lengths. They are not a replacement for Primer3 or IDT OligoAnalyzer; the
# emitted design says so.


def _longest_homopolymer(sequence: str) -> int:
    best = run = 1
    for previous, current in zip(sequence, sequence[1:]):
        run = run + 1 if current == previous else 1
        best = max(best, run)
    return best


def _self_complementarity(sequence: str, window: int = 5) -> int:
    """Longest stretch of the primer that is reverse-complementary to elsewhere in itself."""
    rc = str(Seq(sequence).reverse_complement())
    best = 0
    for size in range(window, min(len(sequence), 12) + 1):
        for start in range(len(sequence) - size + 1):
            if sequence[start : start + size] in rc:
                best = max(best, size)
    return best


def _three_prime_self_dimer(sequence: str, window: int = 5) -> bool:
    """True if the 3' end can fold back onto the primer -- the worst kind of dimer."""
    tail = sequence[-window:]
    return str(Seq(tail).reverse_complement()) in sequence[:-window]


def _binding_sites(template: str, probe: str) -> int:
    """Occurrences of ``probe`` on either strand of the circular ``template``.

    The template is extended by ``len(probe) - 1`` bases rather than doubled: doubling would
    report every single-copy site twice and make every primer look like a mispriming risk.
    """
    probe = probe.upper()
    template = template.upper()
    extended = template + template[: len(probe) - 1]
    reverse = str(Seq(probe).reverse_complement())
    return extended.count(probe) + extended.count(reverse)


def _penalty(sequence: str, tm: float, target_tm: float, min_gc: float, max_gc: float) -> float:
    score = abs(tm - target_tm)
    gc = gc_fraction(sequence) * 100
    if gc < min_gc:
        score += 5 + (min_gc - gc) * 0.2
    elif gc > max_gc:
        score += 5 + (gc - max_gc) * 0.2
    if sequence[-1].upper() not in "GC":
        score += 2.0  # a G/C clamp stabilises the 3' end
    if sum(base in "GC" for base in sequence[-5:].upper()) > 3:
        score += 2.0  # too G/C-rich a 3' end promotes mispriming
    if _longest_homopolymer(sequence) >= 4:
        score += 2.0
    score += 1.0 * max(0, _self_complementarity(sequence) - 5)
    if _three_prime_self_dimer(sequence):
        score += 3.0
    return score


def _candidates(template, boundary, strand, min_length, max_length, target_tm, min_gc, max_gc):
    """Every allowed length for one primer, scored. 5' end is pinned to ``boundary``."""
    options = []
    for length in range(min_length, max_length + 1):
        if strand > 0:
            anneal = slice_circular(template, boundary, boundary + length)
            span = (boundary, boundary + length)
            sequence = anneal
        else:
            window = slice_circular(template, boundary - length, boundary)
            span = (boundary - length, boundary)
            sequence = str(Seq(window).reverse_complement())
        tm = mt.Tm_NN(Seq(sequence))
        options.append(
            {
                "sequence": sequence.upper(),
                "span": span,
                "tm": tm,
                "gc": gc_fraction(sequence) * 100,
                "penalty": _penalty(sequence, tm, target_tm, min_gc, max_gc),
            }
        )
    return options


def _select_pair(
    template, forward_boundary, reverse_boundary, min_length, max_length,
    target_tm, min_gc, max_gc,
):
    """Best (forward, reverse) length combination, balancing per-primer quality and dTm."""
    forward = _candidates(
        template, forward_boundary, +1, min_length, max_length, target_tm, min_gc, max_gc
    )
    reverse = _candidates(
        template, reverse_boundary, -1, min_length, max_length, target_tm, min_gc, max_gc
    )
    _, best_forward, best_reverse = min(
        (
            (f["penalty"] + r["penalty"] + 2.0 * abs(f["tm"] - r["tm"]), f, r)
            for f in forward
            for r in reverse
        ),
        key=lambda triple: triple[0],
    )
    return best_forward, best_reverse


def design_outward_pair(template, span, **kwargs):
    """Primers pointing *away* from ``span`` -- amplifies everything except it.

    Used to delete a span, or to open a vector at a span for assembly.
    """
    start, end = span
    return _select_pair(template, end, start, **kwargs)


def design_inward_pair(template, span, **kwargs):
    """Primers pointing *into* ``span`` -- amplifies the span itself.

    Used to amplify an insert out of a donor plasmid.
    """
    start, end = span
    return _select_pair(template, start, end, **kwargs)


def _annotate(option, template, name, strand, tail, phosphorylated) -> Primer:
    sequence = option["sequence"]
    notes = []
    probe = sequence[-_UNIQUENESS_PROBE:] if len(sequence) >= _UNIQUENESS_PROBE else sequence
    sites = _binding_sites(template, probe)
    if sites > 1:
        notes.append(
            f"3' {len(probe)}-mer occurs {sites}x in the plasmid -- risk of mispriming"
        )
    if sequence[-1] not in "GC":
        notes.append("no G/C clamp at the 3' end")
    if _longest_homopolymer(sequence) >= 4:
        notes.append(f"homopolymer run of {_longest_homopolymer(sequence)}")
    if _three_prime_self_dimer(sequence):
        notes.append("3' end is self-complementary -- check for primer-dimer")
    return Primer(
        name=name,
        sequence=tail + sequence,
        anneal=sequence,
        tail=tail,
        tm=option["tm"],
        gc=option["gc"],
        binding_span=option["span"],
        strand=strand,
        phosphorylated=phosphorylated,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------


def design_deletion_primers(
    plasmid,
    features_to_remove,
    protect=None,
    method: str = KLD,
    target_tm: float = 60.0,
    tm_tolerance: float = 3.0,
    max_tm_difference: float = 3.0,
    min_length: int = 18,
    max_length: int = 32,
    min_gc: float = 35.0,
    max_gc: float = 65.0,
    overlap: int = 20,
    overlap_placement: str = "split",
    name_prefix: str = "del",
    warn_if_restriction_possible: bool = True,
) -> DeletionPCR:
    """Design outward-facing primers that amplify the backbone minus ``features_to_remove``.

    ``method`` is ``"KLD"`` (blunt self-ligation of 5'-phosphorylated ends, the default) or
    ``"gibson"``, which adds ``overlap`` bp of end homology via 5' primer tails.
    ``overlap_placement="split"`` (the default) puts half on each primer so the shared region
    straddles the new junction and neither oligo grows by the full overlap;
    ``"forward"`` puts it all on the forward primer.

    Returns a :class:`DeletionPCR` whose ``product`` is the re-circularized plasmid and whose
    ``forward``/``reverse`` are the oligos to order. Raises :class:`NoPrimerDesign` if the
    target is too close to the origin of a usable primer window.

    Restriction excision is preferred when it is available; with
    ``warn_if_restriction_possible`` this checks for that route and says so in ``warnings``.
    """
    if method not in (KLD, GIBSON):
        raise ValueError(f"method must be {KLD!r} or {GIBSON!r}, got {method!r}")

    record, span, target_labels = resolve_target(load_plasmid(plasmid), features_to_remove)
    target_start, target_end = span
    deleted_bp = target_end - target_start
    removed, protected, umbrella = classify_features(record, span, protect)

    if deleted_bp >= len(record) - 2 * max_length:
        raise NoPrimerDesign(
            f"target is {deleted_bp} bp of a {len(record)} bp plasmid; too little backbone "
            "left to prime against"
        )

    warnings: list[str] = []
    warnings.extend(essential_warnings(removed))
    for start, end, label in protected:
        if start < target_end and target_start < end:
            warnings.append(
                f"{label!r} only partially overlaps the target and will be truncated"
            )
    for label in umbrella:
        warnings.append(f"{label!r} spans the deleted region and will be shortened")

    template = str(record.seq).upper()
    forward_option, reverse_option = design_outward_pair(
        template,
        span,
        min_length=min_length,
        max_length=max_length,
        target_tm=target_tm,
        min_gc=min_gc,
        max_gc=max_gc,
    )

    forward_tail = reverse_tail = ""
    if method == GIBSON:
        if overlap_placement not in ("split", "forward"):
            raise ValueError("overlap_placement must be 'split' or 'forward'")
        if overlap_placement == "split":
            # Half the homology on each primer, so the shared region straddles the new
            # junction and neither oligo grows by the full overlap. An odd base goes to the
            # forward tail.
            upstream_bp = overlap - overlap // 2
            downstream_bp = overlap // 2
        else:
            upstream_bp, downstream_bp = overlap, 0
        forward_tail = slice_circular(template, target_start - upstream_bp, target_start)
        reverse_tail = str(
            Seq(slice_circular(template, target_end, target_end + downstream_bp))
            .reverse_complement()
        )

    forward = _annotate(
        forward_option, template, f"{name_prefix}_F", +1, forward_tail, method == KLD
    )
    reverse = _annotate(
        reverse_option, template, f"{name_prefix}_R", -1, reverse_tail, method == KLD
    )

    for primer in (forward, reverse):
        if abs(primer.tm - target_tm) > tm_tolerance:
            warnings.append(
                f"{primer.name} Tm {primer.tm:.1f} C is more than {tm_tolerance} C from the "
                f"{target_tm} C target"
            )
        warnings.extend(f"{primer.name}: {note}" for note in primer.notes)
    if abs(forward.tm - reverse.tm) > max_tm_difference:
        warnings.append(
            f"primer Tms differ by {abs(forward.tm - reverse.tm):.1f} C; consider a "
            "touchdown or gradient anneal"
        )

    product = delete_span(record, target_start, target_end)
    amplicon_bp = len(record) - deleted_bp + len(forward_tail) + len(reverse_tail)
    if amplicon_bp > 6000:
        warnings.append(
            f"{amplicon_bp} bp amplicon -- use a long-range high-fidelity polymerase and "
            "extend the extension time"
        )
    warnings.append(
        "Tm is nearest-neighbour (Bio.SeqUtils) at default salt; confirm with your "
        "polymerase vendor's calculator. Dimer/hairpin screening here is heuristic -- run the "
        "pair through Primer3 or IDT OligoAnalyzer before ordering."
    )

    if warn_if_restriction_possible:
        from .excise import NoExcisionFound, excise_features

        try:
            alternative = excise_features(record, features_to_remove, protect=protect)
        except (NoExcisionFound, ValueError):
            pass
        else:
            warnings.insert(
                0,
                f"a restriction route exists ({' + '.join(alternative.enzyme_pair)}, "
                f"{alternative.collateral_bp} bp collateral) and is usually preferred over PCR",
            )

    product.name = f"{record.name}_delta"[:16]
    product.description = (
        f"{record.description} | around-the-horn deletion of {', '.join(target_labels)}"
    )

    return DeletionPCR(
        forward=forward,
        reverse=reverse,
        product=product,
        amplicon_bp=amplicon_bp,
        deleted_bp=deleted_bp,
        deleted_span=(target_start, target_end),
        removed_features=removed,
        method=method,
        warnings=warnings,
    )
