"""Where each piece of DNA comes from.

Two ladders, both terminating in "order it" rather than in failure.

**Insert** (a sequence the backbone does not already have)::

    do you know the sequence?
      no  -> ASK, then search (Addgene / literature) and offer candidates
      yes -> already in the inventory?
               yes -> PCR or restriction it out of a donor
               no  -> complexity check
                        pass -> order a synthetic fragment, assemble by Gibson
                        fail -> order the whole plasmid; do not attempt to clone it

**Backbone**::

    a construct in the inventory needing the fewest edits
      -> else a blank/base vector in the inventory (config-tagged)
      -> else order one

The split on a failed complexity check is deliberate. A repetitive block that already exists
somewhere in the library (the 8x BoxB array, present in seven plasmids) is *moved*, never
synthesised, so it never reaches the complexity gate. A genuinely new sequence that a vendor
refuses has nothing to clone from -- hence ordering the finished plasmid.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .complexity import REJECT_SCORE, ComplexityReport, complexity_score
from .config import is_available, is_base_vector, load_lab_config
from .inventory import Provenance, find_sequence

ASK_USER = "ask_user"
FROM_INVENTORY = "from_inventory"
SYNTHESIZE = "synthesize"
ORDER_PLASMID = "order_plasmid"

FROM_CONSTRUCT = "construct"
FROM_BASE_VECTOR = "base_vector"
ORDER_BACKBONE = "order_backbone"


@dataclass
class SourcingDecision:
    """How to obtain one insert."""

    route: str
    name: str
    sequence: str | None = None
    donors: list[Provenance] = field(default_factory=list)
    complexity: ComplexityReport | None = None
    rationale: list[str] = field(default_factory=list)
    question: str | None = None

    @property
    def needs_user_input(self) -> bool:
        return self.route == ASK_USER

    @property
    def blocks_cloning(self) -> bool:
        """True when the whole cloning route is abandoned in favour of buying the plasmid."""
        return self.route == ORDER_PLASMID

    def summary(self) -> str:
        if self.route == ASK_USER:
            return f"{self.name}: sequence unknown -- ask the user"
        if self.route == FROM_INVENTORY:
            donors = ", ".join(sorted({d.label for d in self.donors}))
            return f"{self.name}: already in the library ({donors})"
        if self.route == SYNTHESIZE:
            return f"{self.name}: order as a synthetic fragment ({len(self.sequence or '')} bp)"
        return f"{self.name}: cannot be synthesised -- order the whole plasmid"


@dataclass
class BackboneDecision:
    """Where the vector comes from."""

    route: str
    label: str | None = None
    entry: object | None = None
    rationale: list[str] = field(default_factory=list)
    question: str | None = None

    @property
    def must_order(self) -> bool:
        return self.route == ORDER_BACKBONE

    def summary(self) -> str:
        if self.route == FROM_CONSTRUCT:
            return f"backbone: {self.label} (existing construct)"
        if self.route == FROM_BASE_VECTOR:
            return f"backbone: {self.label} (empty base vector)"
        return "backbone: none suitable in the library -- order one"


def source_insert(
    sequence,
    entries,
    name: str = "insert",
    description: str | None = None,
    vendor: str | None = None,
    allow_split: bool = False,
    limits=None,
) -> SourcingDecision:
    """Decide how to obtain ``sequence``. Pass ``sequence=None`` when it is not yet known."""
    if not sequence:
        what = description or name
        return SourcingDecision(
            route=ASK_USER,
            name=name,
            question=(
                f"What sequence do you want for {what}? If you already have it (a new barcode, "
                "for example), paste it. If not, I can search Addgene and the literature and "
                "come back with candidates for you to choose from."
            ),
            rationale=["no sequence supplied, and it cannot be invented"],
        )

    sequence = str(sequence).upper()
    donors = find_sequence(entries, sequence)
    if donors:
        labels = sorted({d.label for d in donors})
        flipped = [d.label for d in donors if d.strand == -1]
        rationale = [
            f"exact match in {len(labels)} plasmid(s): {', '.join(labels)}",
            "moving existing DNA avoids the synthesis complexity limit entirely",
        ]
        if flipped:
            rationale.append(
                f"present on the reverse strand in {', '.join(sorted(set(flipped)))} -- "
                "check the orientation needed"
            )
        return SourcingDecision(
            route=FROM_INVENTORY,
            name=name,
            sequence=sequence,
            donors=donors,
            rationale=rationale,
        )

    report = complexity_score(sequence, vendor=vendor, limits=limits)
    if report.synthesizable:
        return SourcingDecision(
            route=SYNTHESIZE,
            name=name,
            sequence=sequence,
            complexity=report,
            rationale=[
                "not present anywhere in the library",
                f"passes the synthesisability screen ({report})",
            ],
        )

    rationale = [
        "not present anywhere in the library, so there is nothing to clone from",
        f"fails the synthesisability screen: {'; '.join(report.reasons)}",
        "no workaround -- order the finished plasmid from a vendor instead of cloning",
    ]
    if allow_split and report.sequence_length >= 2 * 125:
        rationale.append(
            "optional: splitting into two shorter fragments and doing a 3-fragment Gibson "
            "sometimes passes where one long fragment does not"
        )
    rationale.append(
        f"NOTE: the screen is heuristic (score {report.raw_score:.1f} against a "
        f"{REJECT_SCORE:.0f} reject threshold, {report.verdict}) and its weights are anchored "
        "to a single observed rejection -- confirm with the vendor's own tool before ordering"
    )
    return SourcingDecision(
        route=ORDER_PLASMID,
        name=name,
        sequence=sequence,
        complexity=report,
        rationale=rationale,
    )


def source_backbone(candidates, entries, config=None, spec_backbone: str | None = None):
    """Pick where the vector comes from.

    ``candidates`` are inventory labels already ranked by edit distance to the goal, best
    first; pass an empty list when nothing is close enough to be worth editing.
    """
    config = config or load_lab_config()
    by_label = {entry.label: entry for entry in entries}

    for label in candidates:
        if not is_available(label, config):
            continue
        return BackboneDecision(
            route=FROM_CONSTRUCT,
            label=label,
            entry=by_label.get(label),
            rationale=["closest existing construct; inherits features already present"],
        )

    for label in config.get("base_vectors") or []:
        if label in by_label and is_available(label, config):
            return BackboneDecision(
                route=FROM_BASE_VECTOR,
                label=label,
                entry=by_label[label],
                rationale=[
                    "no construct was close enough to edit",
                    f"{label} is tagged as an empty base vector; every feature must be added",
                ],
            )

    wanted = spec_backbone or "a suitable vector"
    return BackboneDecision(
        route=ORDER_BACKBONE,
        rationale=[
            "no usable construct and no empty base vector in the library",
            f"order {wanted}; planning cannot place cut sites or homology arms until its "
            "sequence is known",
        ],
        question=(
            f"No backbone in the library fits. Order {wanted} -- do you have a reference "
            "sequence for it, or should the plan stop at the order and resume when it arrives?"
        ),
    )
