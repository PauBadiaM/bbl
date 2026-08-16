"""Turn a verified plan into a bench-ready notebook entry.

The plans already know the biology: which enzymes, which primers, how many base pairs move.
What they do not carry is everything between "ligate the vector and the insert" and a pipette
-- reaction volumes, gel percentages, what to write down afterwards. This module fills that
gap and nothing else. **No design decision is taken here**: every enzyme name, coordinate and
base-pair count is copied from the plan, and the plan's own ``protocol`` text is reproduced
verbatim alongside the expanded version, so the two can be checked against each other.

The shape is taken from a real lab-notebook export (eCLM24, June 2026): aims, cloning
strategy, then one section per bench step with its reagent table, then blanks for the results.
A report is a form to work from, not a summary to file.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from datetime import date

from ..enzymes import SAME_ENZYME
from ..plasmid_io import feature_label, feature_span
from . import bench
from .bench import BLANK, Component, Table, bench_config

#: How many features of one label make a tandem array. Arrays recombine at 37 C and slow the
#: polymerase down, which changes two numbers at the bench.
REPEAT_ARRAY_MIN = 3


@dataclass
class Figure:
    svg: str
    caption: str


@dataclass
class Step:
    """One thing you do at the bench."""

    title: str
    body: list[str] = field(default_factory=list)
    tables: list[Table] = field(default_factory=list)
    #: Values to write down once the step is done -- rendered as blank fields.
    record: list[str] = field(default_factory=list)


@dataclass
class DesignReport:
    title: str
    subtitle: str
    meta: list[tuple[str, str]] = field(default_factory=list)
    aim: str | None = None
    strategy: list[str] = field(default_factory=list)
    figures: list[Figure] = field(default_factory=list)
    materials: list[Table] = field(default_factory=list)
    steps: list[Step] = field(default_factory=list)
    #: Carried for callers and tests, but **not rendered**: bench-critical cautions belong in
    #: the step they apply to, and a warnings appendix is where those go to be skipped.
    warnings: list[str] = field(default_factory=list)
    #: The plan's own protocol text, reproduced without alteration.
    plan_protocol: str = ""

    def to_html(self) -> str:
        from .html import render_html

        return render_html(self)


# ---------------------------------------------------------------------------
# small readings of a record
# ---------------------------------------------------------------------------


def _repeat_array(record) -> str | None:
    """Label of a tandem repeat array in ``record``, if it has one.

    Detected structurally -- a feature containing three or more identically-labelled
    sub-features -- rather than by matching names like "x8", so it also catches an array
    whose umbrella is named something else.
    """
    features = list(getattr(record, "features", []))
    for outer in features:
        o_start, o_end = feature_span(outer)
        labels = [
            feature_label(f)
            for f in features
            if f is not outer
            and o_start <= feature_span(f)[0]
            and feature_span(f)[1] <= o_end
            and feature_span(f) != (o_start, o_end)
        ]
        if len(labels) >= REPEAT_ARRAY_MIN and len(set(labels)) == 1:
            return feature_label(outer)
    return None


def _spans_of(record, labels) -> list[tuple[int, int]]:
    """Spans of the features named in ``labels`` -- the parts a map must not fold away."""
    wanted = {str(label).lower() for label in labels or ()}
    found = [
        feature_span(f) for f in getattr(record, "features", []) if feature_label(f).lower() in wanted
    ]
    return [max(found, key=lambda s: s[1] - s[0])] if found else []


def _renamed(record, name):
    """A shallow view of ``record`` carrying ``name``.

    The caller's record is left alone -- it belongs to the inventory and other plans may hold
    it -- but everything in this module can then read ``.name`` and get the full construct
    name instead of the truncated GenBank LOCUS.
    """
    if record is None or not name or getattr(record, "name", None) == name:
        return record
    view = copy.copy(record)
    view.name = name
    return view


def _component(record_or_name, length=None, concentrations=None):
    """A :class:`Component` with its measured concentration attached, if we were given one."""
    name = getattr(record_or_name, "name", record_or_name)
    if length is None and not isinstance(record_or_name, str):
        length = len(record_or_name)
    return Component(str(name), length, (concentrations or {}).get(str(name)))


def _short(name) -> str:
    """The plasmid ID alone, for places a full name would be repeated to no purpose.

    Full names belong in the header, the tables and the figures, where the reader is
    identifying the tube. A gel lane list is not one of those places.
    """
    return str(name).split("_")[0] if name else str(name)


def _sentence(text: str) -> str:
    """Capitalise the first letter and nothing else -- ``str.capitalize`` flattens 'XhoI'."""
    text = str(text).strip()
    if not text:
        return text
    text = text[0].upper() + text[1:]
    return text if text.endswith((".", "!", "?", ":")) else text + "."


def _primer_table(primers, title="Oligos to order") -> Table:
    rows = []
    for primer in primers:
        rows.append(
            [
                primer.name,
                f"5'-{primer.sequence}-3'",
                f"{primer.length}",
                f"{primer.tm:.1f}",
                f"{primer.gc:.0f}",
                "yes" if primer.phosphorylated else "no",
            ]
        )
    return Table(
        title=title,
        columns=["Name", "Sequence", "nt", "Tm (°C)", "GC (%)", "5'-phosphate"],
        rows=rows,
        note="Tm is nearest-neighbour at default salt; check the oligos in Primer3 or IDT "
        "OligoAnalyzer before ordering.",
    )


def _two_band_caveat(*lengths) -> list[str]:
    """One gel cannot resolve a 300 bp band and a 6 kb band. Say so rather than pretend."""
    small, large = min(lengths), max(lengths)
    if small and large / small < 5:
        return []
    return [
        f"The percentage above is set by the {small:,} bp band; the {large:,} bp band barely "
        "moves into that gel. Run the two on separate gels if you need both clean."
    ]


def _clone_table(colonies: int) -> Table:
    """Filled in after sequencing comes back, so every cell is an input."""
    return Table(
        title="Clones",
        columns=["Colony", "Sequencing result", "Concentration (ng/µL)", "Keep?"],
        rows=[[f"{i}", "", "", ""] for i in range(1, colonies + 1)],
        editable=[1, 2, 3],
        expandable=True,
    )


def _transformation_step(cfg, repeat, junctions) -> Step:
    transform = cfg["transformation"]
    body = bench.transformation_steps(cfg, bool(repeat))
    for name, sequence in junctions:
        body.append(f"Screen across the {name}: the sequence there becomes {sequence}.")
    if repeat:
        # This used to live in a warnings section. It is an instruction, not a caveat: a
        # shortened array reads as a clean sequence unless somebody counts the repeats.
        body.append(
            f"Check the length of the {repeat} array on the sequencing read rather than "
            "assuming it survived — a recombined clone aligns cleanly, just shorter."
        )
    return Step(
        title="Transformation, plating and colony picking",
        body=body,
        tables=[_clone_table(transform["colonies"])],
        record=["Colonies on the plate", "Colonies on the negative control"],
    )


# ---------------------------------------------------------------------------
# route-specific procedures
# ---------------------------------------------------------------------------


def _excision_steps(plan, parent, cfg, repetitive, stock=None) -> list[Step]:
    """Delete by digest and religation."""
    enzymes = list(dict.fromkeys(plan.enzymes))
    single = plan.strategy == SAME_ENZYME
    parent_component = _component(parent, concentrations=stock)
    backbone_bp = len(plan.product)
    excised_bp = plan.deleted_bp

    steps = [
        Step(
            title="1. Digest the parent",
            body=[
                f"Cut {parent.name} with "
                + (
                    f"{enzymes[0]}, which cuts twice and brackets the target."
                    if single
                    else f"{enzymes[0]} + {enzymes[1]} in a double digest."
                ),
                "Run an undigested control alongside — it is the only way to tell a partial "
                "digest from a bad gel.",
            ],
            tables=[
                bench.enzyme_conditions_table(enzymes, cfg),
                bench.digestion_table(parent_component, enzymes, cfg, dephosphorylate=single),
            ],
            record=(
                []
                if parent_component.ng_per_ul
                else [f"{_short(parent.name)} stock concentration (ng/µL)"]
            ),
        ),
        Step(
            title="2. Gel-purify the backbone",
            body=[
                f"Two fragments come off the digest: the {backbone_bp:,} bp backbone you want "
                f"and the {excised_bp:,} bp excised piece. Cut out the large band.",
                bench.cleanup_note(cfg),
            ],
            tables=[
                bench.gel_table(
                    backbone_bp,
                    [
                        "Ladder",
                        f"{_short(parent.name)} digested",
                        f"{_short(parent.name)} undigested",
                    ],
                    cfg,
                )
            ],
            record=["Gel slice mass (mg)", "Eluate concentration (ng/µL)"],
        ),
        Step(
            title="3. Religate",
            body=[
                f"This is an intramolecular ligation — the backbone closes on itself, so there "
                f"is no insert and no molar ratio to set. Keep the DNA dilute (2–5 ng/µL in the "
                f"reaction) so the ends find each other rather than another molecule.",
                f"Ligate with {cfg['ligation']['enzyme']} at {cfg['ligation']['incubation']}. "
                f"The junction becomes {plan.junction_seq}.",
            ],
            record=["Ligation volume (µL)", "DNA in the ligation (ng)"],
        ),
    ]
    if single:
        steps[2].body.append(
            f"The {enzymes[0]} site is regenerated on ligation, so the excised fragment can go "
            "back in. The vector was dephosphorylated in step 1 to suppress that; expect some "
            "re-insertion anyway and screen for it."
        )
    else:
        steps[2].body.append(
            f"{enzymes[0]} and {enzymes[1]} leave compatible ends, which is what allows the "
            "backbone to close. Do not dephosphorylate here — self-closure is the product."
        )
    return steps


def _pcr_deletion_steps(plan, parent, cfg, repetitive) -> list[Step]:
    """Delete by inverse PCR and re-circularization."""
    anneal = min(plan.forward.tm, plan.reverse.tm) - 2
    gibson = plan.method == "gibson"
    body = [
        f"The two primers sit back-to-back at the deletion boundaries pointing away "
        f"from each other, so the {plan.deleted_bp:,} bp target is never copied. "
        f"Amplicon: {plan.amplicon_bp:,} bp from ~1 ng of {parent.name}.",
        "Use a high-fidelity polymerase — every base in this product is PCR-derived "
        "and will have to be sequenced.",
    ]
    if repetitive:
        # The plan quotes a generic 30 s/kb. Leaving the two numbers to be noticed
        # independently is how a reader ends up trusting the wrong one.
        body.append(
            f"The extension below is longer than the {plan.amplicon_bp / 1000 * 30:.0f} s in "
            "the verified plan above: that figure is the generic 30 s/kb, and this amplicon "
            "carries a repeat array."
        )
    steps = [
        Step(
            title="1. Amplify the backbone outward",
            body=body,
            tables=[
                bench.pcr_table(cfg),
                bench.thermocycler_table(anneal, plan.amplicon_bp, cfg, repetitive=repetitive),
            ],
            record=["Yield (ng/µL)"],
        ),
        Step(
            title="2. Check and purify the amplicon",
            body=[
                f"One band at {plan.amplicon_bp:,} bp. A smear or a second band means "
                "mispriming — re-run with a higher annealing temperature before going on.",
                bench.cleanup_note(cfg),
            ],
            tables=[
                bench.gel_table(plan.amplicon_bp, ["Ladder", "PCR product"], cfg),
            ],
            record=["Gel slice mass (mg)", "Eluate concentration (ng/µL)"],
        ),
        Step(
            title="3. Remove the template",
            body=[
                "DpnI-digest the reaction: the parental plasmid is Dam-methylated from a dam+ "
                "host and is cut, the PCR product is not. Skipping this is the usual reason a "
                "'deletion' clone comes back as the parent.",
            ],
        ),
    ]
    if gibson:
        steps.append(
            Step(
                title="4. Circularize by Gibson",
                body=[
                    f"The 5' tails give the amplicon a {plan.homology_bp} bp terminal direct "
                    "repeat spanning the new junction, so it closes on itself. No "
                    "phosphorylation needed.",
                    f"Incubate {cfg['gibson']['incubation']}.",
                ],
            )
        )
    else:
        steps.append(
            Step(
                title="4. Phosphorylate and blunt self-ligate",
                body=[
                    "Order both primers 5'-phosphorylated, or treat the amplicon with T4 PNK. "
                    "A KLD (kinase / ligase / DpnI) mix does steps 3 and 4 in one tube.",
                    "Blunt ligation is inefficient — keep the DNA dilute and expect fewer "
                    "colonies than a sticky-end reaction.",
                ],
            )
        )
    return steps


def _insertion_steps(plan, vector_record, cfg, repetitive, donor=None, stock=None) -> list[Step]:
    """Subclone by restriction, or assemble by Gibson."""
    from ..insert import RESTRICTION

    steps: list[Step] = []
    vector = _component(plan.vector.name, plan.vector.length, stock)
    insert = _component(plan.insert.name, plan.insert.length, stock)
    #: What actually goes in the donor digest is the whole donor plasmid, not the band.
    donor_component = (
        _component(donor, concentrations=stock) if donor is not None
        else Component("donor plasmid")
    )

    if plan.strategy == RESTRICTION:
        up, down = plan.vector.enzymes
        steps.append(
            Step(
                title="1. Digest the vector",
                body=[
                    _sentence(plan.vector.preparation),
                    f"The backbone runs at {plan.vector.length:,} bp.",
                ],
                tables=[
                    bench.enzyme_conditions_table([up, down], cfg),
                    bench.digestion_table(
                        _component(vector_record, concentrations=stock),
                        [up, down],
                        cfg,
                        dephosphorylate=not plan.directional,
                    ),
                ],
                record=(
                    []
                    if (stock or {}).get(vector_record.name)
                    else [f"{_short(vector_record.name)} stock concentration (ng/µL)"]
                ),
            )
        )
        steps.append(
            Step(
                title="2. Digest the donor and cut out the insert",
                body=[
                    _sentence(plan.insert.preparation),
                    "Extra cuts in the donor backbone are fine — they only help separate the "
                    "band you want.",
                ],
                tables=[
                    bench.digestion_table(
                        donor_component, list(plan.insert.enzymes or ()), cfg
                    )
                ],
                record=(
                    []
                    if donor_component.ng_per_ul
                    else [f"{_short(donor_component.name)} stock concentration (ng/µL)"]
                ),
            )
        )
        steps.append(
            Step(
                title="3. Gel-purify both fragments",
                body=[bench.cleanup_note(cfg)]
                + _two_band_caveat(plan.vector.length, plan.insert.length),
                tables=[
                    bench.gel_table(
                        min(plan.vector.length, plan.insert.length),
                        ["Ladder", "Digested vector", "Digested donor", "Undigested controls"],
                        cfg,
                    )
                ],
                record=[
                    "Vector eluate (ng/µL)",
                    "Insert eluate (ng/µL)",
                ],
            )
        )
        steps.append(
            Step(
                title="4. Ligate",
                body=[
                    (
                        "The two vector ends are mutually incompatible, so the insert can only "
                        "go in one way and the vector cannot self-close."
                        if plan.directional
                        else "The vector's ends are compatible with each other: it can "
                        "self-close and the insert can go in either orientation. The vector "
                        "was dephosphorylated in step 1 — screen the orientation regardless."
                    )
                ],
                tables=[bench.ligation_table(vector, insert, cfg)],
                record=["Ligation volume (µL)"],
            )
        )
        return steps

    # Gibson: prepare each fragment, then assemble.
    prepare = Step(title="1. Prepare the vector", body=[_sentence(plan.vector.preparation)])
    if plan.vector.primers:
        anneal = min(p.tm for p in plan.vector.primers) - 2
        prepare.tables = [
            bench.pcr_table(cfg),
            bench.thermocycler_table(anneal, plan.vector.length, cfg, repetitive=repetitive),
        ]
    elif plan.vector.enzymes:
        prepare.tables = [
            bench.enzyme_conditions_table(list(plan.vector.enzymes), cfg),
            bench.digestion_table(
                _component(vector_record, concentrations=stock),
                list(plan.vector.enzymes),
                cfg,
            ),
        ]
    prepare.record = ["Vector prep concentration (ng/µL)"]
    steps.append(prepare)

    second = Step(title="2. Prepare the insert", body=[_sentence(plan.insert.preparation)])
    if plan.insert.primers:
        anneal = min(p.tm for p in plan.insert.primers) - 2
        second.tables = [
            bench.pcr_table(cfg),
            bench.thermocycler_table(anneal, plan.insert.length, cfg, repetitive=repetitive),
        ]
    if plan.insert.order_sequence:
        second.body.append(
            "Nothing to do at the bench until it arrives; resuspend to 10 ng/µL on delivery."
        )
    second.record = ["Insert concentration (ng/µL)"]
    steps.append(second)

    steps.append(
        Step(
            title="3. Gel-purify",
            body=[bench.cleanup_note(cfg)]
            + _two_band_caveat(plan.vector.length, plan.insert.length),
            tables=[
                bench.gel_table(
                    min(plan.vector.length, plan.insert.length),
                    ["Ladder", "Vector", "Insert"],
                    cfg,
                )
            ],
            record=["Vector eluate (ng/µL)", "Insert eluate (ng/µL)"],
        )
    )
    steps.append(
        Step(
            title="4. Assemble",
            body=[
                "The homology arms are unique to each end, so the assembly is directional and "
                "the vector cannot close without the insert.",
            ],
            tables=[bench.gibson_table(vector, [insert], cfg)],
        )
    )
    return steps


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------


def _strategy_lines(plan, parent, product) -> list[str]:
    """The two-or-three-line summary a notebook entry opens with."""
    from ..insert import RESTRICTION
    from ..pcr import DeletionPCR

    parent_name = getattr(parent, "name", "the parent")
    if hasattr(plan, "enzyme_pair") and hasattr(plan, "deleted_bp"):
        enzymes = list(dict.fromkeys(plan.enzymes))
        return [
            f"Backbone: {parent_name} ({len(parent):,} bp).",
            f"Cut with {' + '.join(enzymes)}"
            + (" (cuts twice, flanking the target)" if len(enzymes) == 1 else "")
            + ", gel-purify the backbone, religate.",
            f"Removes {', '.join(plan.removed_features) or 'unannotated sequence'} — "
            f"{plan.deleted_bp:,} bp, of which {plan.collateral_bp:,} bp is collateral.",
            f"Product: {len(product):,} bp.",
        ]
    if isinstance(plan, DeletionPCR):
        return [
            f"Backbone: {parent_name} ({len(parent):,} bp).",
            f"Inverse-PCR outward from the deletion boundaries with {plan.forward.name} + "
            f"{plan.reverse.name} ({plan.amplicon_bp:,} bp amplicon), DpnI, then "
            + ("Gibson self-assembly." if plan.method == "gibson" else "blunt self-ligation."),
            f"Removes {', '.join(plan.removed_features) or 'unannotated sequence'} — "
            f"{plan.deleted_bp:,} bp exactly, no collateral loss.",
            f"Product: {len(product):,} bp.",
        ]
    lines = [
        f"Backbone: {parent_name} ({len(parent):,} bp), opened at {plan.site[0]:,}–"
        f"{plan.site[1]:,}.",
        f"Insert: {plan.insert.name}, {plan.inserted_bp:,} bp.",
    ]
    if plan.strategy == RESTRICTION:
        lines.append(
            f"Subclone with {' + '.join(plan.vector.enzymes)} — "
            + ("directional." if plan.directional else "non-directional, screen orientation.")
        )
    else:
        lines.append("Assemble by Gibson; the homology arms make it directional.")
    lines.append(
        f"Product: {len(product):,} bp ({plan.inserted_bp:,} bp in, {plan.replaced_bp:,} bp out)."
    )
    return lines


def _figures(plan, parent, product) -> list[Figure]:
    """Circular maps of both plasmids, plus a linear zoom on the edit."""
    from .maps import circular_svg, linear_svg, zoom_window

    from ..pcr import DeletionPCR

    figures: list[Figure] = []
    deletion = hasattr(plan, "deleted_span")
    if deletion:
        parent_span = tuple(plan.deleted_span)
        product_span = None
        junction_at = parent_span[0]
        if isinstance(plan, DeletionPCR):
            marks = [(parent_span[0], f"{plan.reverse.name} 5'"), (parent_span[1], f"5' {plan.forward.name}")]
        else:
            enzymes = plan.enzyme_pair
            marks = [(parent_span[0], f"{enzymes[0]} ↓"), (parent_span[1], f"↓ {enzymes[1]}")]
        focus = _spans_of(parent, plan.removed_features)
    else:
        parent_span = tuple(plan.site)
        product_span = (plan.site[0], plan.site[0] + plan.inserted_bp)
        junction_at = plan.site[0]
        enzymes = plan.vector.enzymes or ()
        marks = (
            [(parent_span[0], f"{enzymes[0]} ↓"), (parent_span[1], f"↓ {enzymes[1]}")]
            if len(enzymes) == 2
            else [(parent_span[0], "opened here"), (parent_span[1], "opened here")]
        )
        focus = _spans_of(parent, [plan.insert.name])

    if parent is not None:
        figures.append(
            Figure(
                svg=circular_svg(
                    parent,
                    title=f"{parent.name} — parent",
                    removed=parent_span,
                    protect=focus,
                    salt="parent-circular",
                ),
                caption=(
                    f"{parent.name}, {len(parent):,} bp. The region that "
                    + ("goes away" if deletion else "is replaced")
                    + " is picked out in red."
                ),
            )
        )
    product_focus = _spans_of(product, [plan.insert.name]) if not deletion else []
    figures.append(
        Figure(
            svg=circular_svg(
                product,
                title=f"{product.name} — product",
                added=product_span,
                protect=product_focus,
                salt="product-circular",
            ),
            caption=f"The designed product, {len(product):,} bp."
            + ("" if deletion else " New sequence is picked out in green."),
        )
    )
    if parent is not None:
        window = zoom_window(parent_span, len(parent))
        figures.append(
            Figure(
                svg=linear_svg(
                    parent,
                    window,
                    marks=marks,
                    removed=parent_span,
                    protect=focus,
                    salt="parent-linear",
                ),
                caption=f"{parent.name} around the edit, {window[0]:,}–{window[1]:,} bp.",
            )
        )
    product_window = zoom_window(
        product_span or (junction_at, junction_at), len(product)
    )
    figures.append(
        Figure(
            svg=linear_svg(
                product,
                product_window,
                marks=[(junction_at, "new junction")] if deletion else [],
                added=product_span,
                protect=product_focus,
                salt="product-linear",
            ),
            caption=f"The product across the same region, in the parent's coordinate frame.",
        )
    )
    return figures


def build_report(
    plan,
    parent=None,
    product=None,
    donor=None,
    aim: str | None = None,
    name: str | None = None,
    parent_name: str | None = None,
    donor_name: str | None = None,
    author: str | None = None,
    config=None,
    concentrations=None,
    figures: bool = True,
    verified_against=None,
) -> DesignReport:
    """Assemble a :class:`DesignReport` from a verified plan.

    ``parent`` is the record the plan started from -- without it the before/after figures and
    the digest tables cannot be drawn, so pass it whenever it is known. ``donor`` is the
    plasmid an insert is cut out of, needed for that digest's reagent table.

    ``name``, ``parent_name`` and ``donor_name`` are the **full construct names**
    (``pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8``). A loaded record only keeps a
    16-character GenBank LOCUS name, which is an accession, not a description of what is in
    the tube -- and at the bench the difference between pCLM21 and pCLM24 is the part after
    the underscore. Callers that know the inventory label should pass it.

    ``concentrations`` maps a plasmid name to ng/µL; anything missing is left blank in the
    tables rather than guessed. ``verified_against`` is a ``(target_name, identical)`` pair
    from ``compare_product``.
    """
    from ..insert import RESTRICTION, InsertionPlan
    from ..pcr import DeletionPCR

    product = product if product is not None else plan.product
    cfg = bench_config(config)
    concentrations = concentrations or {}
    repeat = _repeat_array(product)

    if isinstance(plan, InsertionPlan):
        kind = "insertion"
        route = "restriction subcloning" if plan.strategy == RESTRICTION else "Gibson assembly"
    elif isinstance(plan, DeletionPCR):
        kind = "deletion"
        route = f"inverse PCR ({plan.method})"
    else:
        kind = "deletion"
        route = f"restriction digest ({plan.strategy.replace('_', ' ')})"

    parent_name = parent_name or getattr(parent, "name", None)
    donor_name = donor_name or getattr(donor, "name", None)
    product_name = name or getattr(product, "name", None) or "designed construct"
    if product_name == parent_name or product_name == getattr(parent, "name", None):
        # The restriction route keeps the vector's name, which would print the parent and the
        # product under one label in a document whose whole job is telling them apart.
        product_name = f"{parent_name}_designed"

    # Rebind to renamed views so every downstream ``.name`` is the full construct name.
    parent = _renamed(parent, parent_name)
    donor = _renamed(donor, donor_name)
    product = _renamed(product, product_name)

    meta = [
        ("Product", f"{product_name} — {len(product):,} bp"),
        ("Parent", f"{parent_name} — {len(parent):,} bp" if parent is not None else "—"),
        ("Route", route),
        ("Generated", date.today().isoformat()),
    ]
    if author:
        meta.insert(0, ("Author", author))
    if verified_against:
        target, identical = verified_against
        meta.append(
            (
                "Verified",
                f"identical to {target}" if identical else f"differs from {target}",
            )
        )

    report = DesignReport(
        title=f"{product_name} — cloning report",
        subtitle=(
            f"{kind.capitalize()} by {route}"
            + (f", from {parent_name}" if parent_name else "")
        ),
        meta=meta,
        aim=aim,
        strategy=_strategy_lines(plan, parent, product),
        plan_protocol=getattr(plan, "protocol", ""),
        warnings=list(getattr(plan, "warnings", [])),
    )

    # -- materials -----------------------------------------------------------
    dna_rows = []
    if parent is not None:
        dna_rows.append(
            [
                parent.name,
                f"{len(parent):,} bp",
                "template / backbone",
                str(concentrations.get(parent.name, BLANK)),
            ]
        )
    if donor is not None:
        dna_rows.append(
            [
                donor.name,
                f"{len(donor):,} bp",
                f"donor for {plan.insert.name}",
                str(concentrations.get(donor.name, BLANK)),
            ]
        )
    elif isinstance(plan, InsertionPlan) and plan.insert.order_sequence is None:
        dna_rows.append(
            [
                plan.insert.name,
                f"{plan.insert.length:,} bp",
                "insert",
                str(concentrations.get(plan.insert.name, BLANK)),
            ]
        )
    if dna_rows:
        report.materials.append(
            Table(
                title="DNA",
                columns=["Plasmid / fragment", "Length", "Role", "Stock (ng/µL)"],
                rows=dna_rows,
                note=(
                    "Every reaction volume below is calculated from these concentrations. "
                    "Where one is missing, type it into the table and the volumes follow."
                ),
            )
        )

    primers = []
    if isinstance(plan, DeletionPCR):
        primers = [plan.forward, plan.reverse]
    elif isinstance(plan, InsertionPlan):
        primers = list(plan.vector.primers or ()) + list(plan.insert.primers or ())
    if primers:
        report.materials.append(_primer_table(primers))

    if isinstance(plan, InsertionPlan) and plan.insert.order_sequence:
        report.materials.append(
            Table(
                title="Synthetic fragment to order",
                columns=["Name", "Length", "Sequence (homology arms attached)"],
                rows=[
                    [
                        plan.insert.name,
                        f"{len(plan.insert.order_sequence):,} bp",
                        plan.insert.order_sequence,
                    ]
                ],
                note="Run this past the vendor's complexity screen before ordering.",
            )
        )

    enzymes = []
    if hasattr(plan, "enzymes") and not isinstance(plan, InsertionPlan):
        enzymes = list(dict.fromkeys(plan.enzymes))
    elif isinstance(plan, InsertionPlan):
        enzymes = list(dict.fromkeys((plan.vector.enzymes or ()) + (plan.insert.enzymes or ())))
    if enzymes:
        stock = (config or {}).get("enzyme_stock")
        report.materials.append(
            Table(
                title="Enzymes",
                columns=["Enzyme", "In stock?"],
                rows=[
                    [name, BLANK if stock is None else ("yes" if name in stock else "order")]
                    for name in enzymes
                ],
                note=None
                if stock is not None
                else "No enzyme stock list configured — check the freezer before starting.",
            )
        )

    # -- procedure -----------------------------------------------------------
    if isinstance(plan, InsertionPlan):
        steps = _insertion_steps(
            plan, parent, cfg, bool(repeat), donor=donor, stock=concentrations
        )
    elif isinstance(plan, DeletionPCR):
        steps = _pcr_deletion_steps(plan, parent, cfg, bool(repeat))
    else:
        steps = _excision_steps(plan, parent, cfg, bool(repeat), stock=concentrations)

    junctions = []
    if hasattr(plan, "junction_seq"):
        junctions = [("new junction", plan.junction_seq)]
    elif isinstance(plan, InsertionPlan):
        junctions = [
            ("upstream junction", plan.junctions[0]),
            ("downstream junction", plan.junctions[1]),
        ]
    steps.append(_transformation_step(cfg, repeat, junctions))
    for index, step in enumerate(steps, start=1):
        if not step.title[0].isdigit():
            step.title = f"{index}. {step.title}"
    report.steps = steps

    # -- figures -------------------------------------------------------------
    if figures:
        report.figures = _figures(plan, parent, product)
    return report
