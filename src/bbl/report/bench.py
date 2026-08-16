"""Reagent tables: the arithmetic between a cloning plan and a pipette.

A plan says "ligate the vector and the insert". A notebook entry says *how many microlitres of
each*, and that is a molarity calculation the experimentalist otherwise does by hand in a
spreadsheet every time. This module does it, and only it -- no biology decisions are made here.

Two rules govern the numbers:

**Never invent a concentration.** Molar ratios are fixed by the design; nanograms per
microlitre come from a Nanodrop reading that does not exist until the DNA is in a tube. When a
concentration is unknown the volume cell is left blank (:data:`BLANK`) and the requirement --
"100 ng", "15 fmol" -- is still stated, so the row is a form to fill in rather than a fiction.

**Defaults are the lab's, not ours.** Kit names, buffer volumes and molar excesses come from
:data:`BENCH_DEFAULTS`, overridable through the ``bench`` key of ``config/lab.json``. They are
seeded from a real notebook entry (eCLM24, June 2026), and the Gibson and ligation calculators
reproduce that entry's own volumes to the microlitre -- see ``tests/test_bench.py``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import fastdigest

#: Average molecular weight of one base pair of dsDNA, in g/mol. 650 is the usual convention;
#: the lab's ligation sheet uses 660. The two differ by 1.5%, which is below pipetting error,
#: but each calculator keeps the constant its source sheet used so numbers reconcile.
MW_PER_BP = 650
MW_PER_BP_LIGATION = 660

#: Printed where a number depends on a measurement nobody has taken yet.
BLANK = "____"

BENCH_DEFAULTS = {
    "supplier_table_url": (
        "https://www.thermofisher.com/us/en/home/brands/thermo-scientific/molecular-biology/"
        "thermo-scientific-restriction-modifying-enzymes/restriction-enzymes-thermo-scientific/"
        "fastdigest-thermo-scientific/reaction-conditions-for-fastdigest-enzymes.html"
    ),
    "digest": {
        "buffer": "10x FastDigest Buffer",
        "buffer_ul": 2,
        "enzyme_ul": 1,
        "micrograms": 1.0,
        "total_ul": 20,
        "phosphatase": "FastAP",
        "phosphatase_ul": 1,
        "phosphatase_incubation": "37 °C for 10 min, then 65 °C for 15 min to inactivate",
    },
    "gel": {
        "volume_ml": 50,
        "buffer": "1x TAE",
        "stain": "SYBR Safe",
        "stain_ul": 5,
        "volts": 140,
        "minutes": 30,
        "ladder_ul": 5,
    },
    "cleanup": {
        "kit": "Zymoclean Gel DNA Recovery",
        "adb_per_mg": 3,  # µL of ADB buffer per mg of excised gel
    },
    "gibson": {
        "master_mix": "Gibson Assembly MasterMix (2x)",
        "master_mix_fraction": 0.5,
        "total_ul": 10,
        "backbone_ng": 100,
        "insert_molar_excess": 3,
        "incubation": "50 °C for 30 min",
    },
    "ligation": {
        "backbone_fmol": 15,
        "insert_molar_excess": 5,
        "enzyme": "T4 DNA ligase",
        "enzyme_ul": 1,
        "buffer_ul": 2,
        "total_ul": 20,
        "incubation": "room temperature for 30 min",
    },
    "pcr": {
        "master_mix": "2X KAPA HiFi HotStart ReadyMix",
        "reaction_ul": 25,
        "primer_um": 10,
        "primer_ul": 0.75,
        "template_ul": 1,
        "template_ng_per_ul": 1,
        "cycles": 30,
        "seconds_per_kb": 30,
        "seconds_per_kb_repetitive": 45,
    },
    "transformation": {
        "cells_ul": 25,
        "dna_ul": 2,
        "soc_ul": 500,
        "outgrowth": "1 h shaking",
        "plated_ul": 50,
        "colonies": 3,
        "sequencing": "Plasmidsaurus whole-plasmid sequencing",
    },
    #: Repeat arrays (BoxB, MS2) recombine at 37 °C. 30 °C and a recA- strain keep them intact.
    "growth_temp_c": 37,
    "growth_temp_c_repetitive": 30,
}


def bench_config(config=None) -> dict:
    """Merge the ``bench`` block of a lab config over :data:`BENCH_DEFAULTS`, one level deep."""
    override = (config or {}).get("bench") or {}
    merged = {}
    for key, value in BENCH_DEFAULTS.items():
        supplied = override.get(key)
        if isinstance(value, dict) and isinstance(supplied, dict):
            merged[key] = {**value, **supplied}
        else:
            merged[key] = value if supplied is None else supplied
    return merged


# ---------------------------------------------------------------------------
# quantities
# ---------------------------------------------------------------------------


@dataclass
class Component:
    """One piece of DNA about to go into a reaction."""

    name: str
    #: ``None`` when the length is not known -- a donor plasmid nobody has told us about.
    length_bp: int | None = None
    #: ng/µL from the Nanodrop. ``None`` until somebody measures it.
    ng_per_ul: float | None = None

    @property
    def described(self) -> str:
        return f"{self.name} ({self.length_bp:,} bp)" if self.length_bp else self.name


@dataclass
class Table:
    """A rendered reagent table. ``rows`` are already formatted strings.

    ``editable`` lists column indices the reader fills in at the bench, rendered as text
    inputs; ``expandable`` adds a button for another row (you picked six colonies, not three).
    """

    title: str
    columns: list[str]
    rows: list[list[str]] = field(default_factory=list)
    note: str | None = None
    editable: list[int] = field(default_factory=list)
    expandable: bool = False


@dataclass
class ReactionRow:
    """One line of a reaction mix.

    A row is either **fixed volume** (buffer, enzyme, master mix -- ``fixed_ul``) or **DNA**
    (a mass in ``ng``, which becomes a volume once ``ng_per_ul`` is known). Water makes up the
    balance, so it is not a row.
    """

    name: str
    fixed_ul: float | None = None
    ng: float | None = None
    fmol: float | None = None
    ng_per_ul: float | None = None
    note: str = ""

    @property
    def is_dna(self) -> bool:
        return self.ng is not None

    @property
    def volume_ul(self) -> float | None:
        if self.fixed_ul is not None:
            return self.fixed_ul
        return _volume(self.ng, self.ng_per_ul) if self.is_dna else None


@dataclass
class Reaction:
    """A reaction mix whose DNA volumes follow from concentrations measured later.

    This is the one table type that has to stay live: the molar ratios are fixed by the
    design, but every volume depends on a Nanodrop reading taken after the report is written.
    Rendered with the stock column as an input box and the volumes recomputed in the page, so
    the reader types two numbers instead of redoing the arithmetic in a spreadsheet.
    """

    title: str
    rows: list[ReactionRow]
    total_ul: float
    note: str = ""
    water_name: str = "Water, nuclease-free"

    @property
    def water_ul(self) -> float | None:
        used = 0.0
        for row in self.rows:
            volume = row.volume_ul
            if volume is None:
                return None
            used += volume
        return self.total_ul - used

    def spec(self) -> dict:
        """The calculation, as data, for the page to recompute against."""
        return {
            "total": self.total_ul,
            "rows": [
                {"fixed": row.fixed_ul, "ng": row.ng, "conc": row.ng_per_ul}
                for row in self.rows
            ],
        }


def ng_to_fmol(ng: float, length_bp: int, mw_per_bp: int = MW_PER_BP) -> float:
    return ng * 1e6 / (length_bp * mw_per_bp)


def fmol_to_ng(fmol: float, length_bp: int, mw_per_bp: int = MW_PER_BP) -> float:
    return fmol * length_bp * mw_per_bp / 1e6


def _volume(ng: float, ng_per_ul: float | None) -> float | None:
    """Microlitres delivering ``ng``, or ``None`` when the stock is unmeasured."""
    if not ng_per_ul:
        return None
    return ng / ng_per_ul


def _ul(value: float | None) -> str:
    return BLANK if value is None else f"{value:.2f}"


def _ng(value: float | None) -> str:
    return BLANK if value is None else f"{value:.1f}"


# ---------------------------------------------------------------------------
# digestion
# ---------------------------------------------------------------------------


def enzyme_conditions_table(enzymes, cfg) -> Table:
    """Incubation and heat-inactivation per enzyme, filled from the supplier's own table.

    The numbers come from :mod:`bbl.report.fastdigest`, a dated scrape of Thermo's published
    conditions -- not from memory. An enzyme Thermo does not sell as FastDigest gets blank
    cells and the link, because the alternative is a plausible-looking invention.
    """
    rows = []
    for name in enzymes:
        supplier_name, row = fastdigest.lookup(name)
        if row is None:
            rows.append([name, BLANK, BLANK, BLANK, BLANK])
            continue
        temp, minutes, bp_end, inactivation, star_hours = row
        label = name if supplier_name == name else f"{name} (sold as {supplier_name})"
        rows.append(
            [label, f"{temp} °C, {minutes} min", inactivation, bp_end, f"{star_hours} h"]
        )
    return Table(
        title="Enzyme conditions — 1 µg plasmid DNA in 20 µL",
        columns=[
            "Enzyme",
            "Incubation",
            "Thermal inactivation",
            "bp from end",
            "Star-activity-free",
        ],
        rows=rows,
        note=(
            f"From Thermo's FastDigest conditions table, retrieved {fastdigest.RETRIEVED}: "
            f"{cfg['supplier_table_url']} — check it if a lot behaves differently. "
            "'bp from end' is how much sequence must flank the site for a terminal cut to work; "
            "'star-activity-free' is how long the digest can run before the enzyme starts "
            "cutting the wrong sites."
        ),
    )


def digest_conditions(enzymes) -> dict:
    """Incubation and inactivation for a digest, reconciled across all of its enzymes.

    The slowest enzyme sets the time and the least heat-labile sets the inactivation, because
    a double digest runs in one tube. ``star_limit_h`` is the shortest star-activity-free
    window in the set -- the reason a FastDigest reaction is five minutes and not an hour.
    """
    minutes, temps, inactivations, star = [], [], [], []
    unknown, no_heat_kill = [], []
    for name in enzymes:
        supplier_name, row = fastdigest.lookup(name)
        if row is None:
            unknown.append(name)
            continue
        temp, plasmid_min, _bp, inactivation, star_hours = row
        try:
            minutes.append(int(plasmid_min))
        except ValueError:
            pass
        temps.append(temp)
        star.append(float(star_hours))
        if inactivation.lower().startswith("no"):
            no_heat_kill.append(supplier_name)
        else:
            inactivations.append(inactivation)
    return {
        "minutes": max(minutes) if minutes else None,
        "temp_c": temps[0] if len(set(temps)) == 1 else None,
        "temps": sorted(set(temps)),
        "inactivations": sorted(set(inactivations)),
        "no_heat_kill": no_heat_kill,
        "star_limit_h": min(star) if star else None,
        "unknown": unknown,
    }


def digest_incubation_note(enzymes) -> str:
    """One sentence saying how long to run the digest and how to stop it."""
    found = digest_conditions(enzymes)
    if found["minutes"] is None:
        return "Incubate at 37 °C; check the supplier's table for the time and inactivation."

    temperature = f"{found['temp_c']} °C" if found["temp_c"] else (
        " then ".join(f"{t} °C" for t in found["temps"])
        + " (start with the lower-temperature enzyme)"
    )
    text = f"Incubate at {temperature} for {found['minutes']} min."
    if found["star_limit_h"] is not None and found["star_limit_h"] < 1:
        text += (
            f" Do not extend it: one of these enzymes shows star activity beyond "
            f"{found['star_limit_h']:g} h, so a 'just to be sure' overnight digest cuts "
            "sites you did not plan for."
        )
    if found["no_heat_kill"]:
        text += (
            f" {', '.join(found['no_heat_kill'])} cannot be heat-inactivated — column-purify "
            "or phenol/chloroform-extract before ligating, rather than relying on a heat step."
        )
        if found["inactivations"]:
            text += f" ({', '.join(found['inactivations'])} would handle the other(s).)"
    elif found["inactivations"]:
        text += f" Inactivate: {', '.join(found['inactivations'])}."
    if found["unknown"]:
        text += f" No FastDigest entry for {', '.join(found['unknown'])} — look it up."
    return text


def digestion_table(component: Component, enzymes, cfg, dephosphorylate=False) -> Reaction:
    """Restriction digest mix, back-calculated from the target mass of DNA."""
    digest = cfg["digest"]
    micrograms = digest["micrograms"]

    rows = [
        ReactionRow(digest["buffer"], fixed_ul=digest["buffer_ul"]),
        ReactionRow(
            component.described,
            ng=micrograms * 1000,
            ng_per_ul=component.ng_per_ul,
            note=f"{micrograms:g} µg",
        ),
    ]
    rows += [ReactionRow(f"{name} (FastDigest)", fixed_ul=digest["enzyme_ul"]) for name in enzymes]
    if dephosphorylate:
        rows.append(
            ReactionRow(
                digest["phosphatase"],
                fixed_ul=digest["phosphatase_ul"],
                note="dephosphorylation",
            )
        )

    return Reaction(
        title=f"Digestion — {component.name}",
        rows=rows,
        total_ul=digest["total_ul"],
        note=(
            "Assemble at room temperature in the order listed, flick to mix, spin down. "
            + digest_incubation_note(enzymes)
            + (
                f" Add {digest['phosphatase']} with the enzymes and incubate "
                f"{digest['phosphatase_incubation']}."
                if dephosphorylate
                else ""
            )
        ),
    )


# ---------------------------------------------------------------------------
# gel and clean-up
# ---------------------------------------------------------------------------


def gel_percent(fragment_bp: int) -> float:
    """Agarose percentage that resolves a band of this size."""
    if fragment_bp < 500:
        return 1.5
    if fragment_bp < 1000:
        return 1.2
    return 1.0


def gel_table(fragment_bp: int, lanes, cfg) -> Table:
    """Gel recipe sized to the band you are trying to see, plus the lane plan."""
    gel = cfg["gel"]
    percent = gel_percent(fragment_bp)
    grams = percent / 100 * gel["volume_ml"]
    rows = [
        ["Agarose", f"{grams:.2f} g", f"{percent:g}% gel — resolves the {fragment_bp:,} bp band"],
        [gel["buffer"], f"{gel['volume_ml']} mL", "melt in 30 s bursts until dissolved"],
        [gel["stain"], f"{gel['stain_ul']} µL", "add once cooled slightly, swirl"],
        ["Run", f"{gel['volts']} V, {gel['minutes']} min", "check polarity — running to red"],
    ]
    return Table(
        title=f"Gel — {percent:g}% agarose",
        columns=["Component", "Amount", "Notes"],
        rows=rows,
        note="Lanes: " + "; ".join(f"{i}. {lane}" for i, lane in enumerate(lanes, start=1)),
    )


def cleanup_note(cfg) -> str:
    clean = cfg["cleanup"]
    return (
        f"Excise the band and purify with the {clean['kit']} kit: weigh the gel slice and add "
        f"{clean['adb_per_mg']} µL of ADB per mg. Record the gel mass and the eluate "
        f"concentration below — later reactions are calculated from it."
    )


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------


def gibson_table(backbone: Component, inserts, cfg) -> Reaction:
    """Gibson/NEBuilder mix at a fixed molar excess of each insert over the backbone."""
    gibson = cfg["gibson"]
    total = gibson["total_ul"]
    excess = gibson["insert_molar_excess"]
    backbone_ng = gibson["backbone_ng"]
    backbone_fmol = ng_to_fmol(backbone_ng, backbone.length_bp)

    rows = [
        ReactionRow(
            gibson["master_mix"], fixed_ul=total * gibson["master_mix_fraction"]
        ),
        ReactionRow(
            f"{backbone.name} (backbone, {backbone.length_bp:,} bp)",
            ng=backbone_ng,
            fmol=backbone_fmol,
            ng_per_ul=backbone.ng_per_ul,
        ),
    ]
    for insert in inserts:
        fmol = backbone_fmol * excess
        rows.append(
            ReactionRow(
                f"{insert.name} (insert, {insert.length_bp:,} bp)",
                ng=fmol_to_ng(fmol, insert.length_bp),
                fmol=fmol,
                ng_per_ul=insert.ng_per_ul,
            )
        )
    return Reaction(
        title=f"Gibson assembly — {excess}:1 molar excess of insert",
        rows=rows,
        total_ul=total,
        note=(
            f"Incubate {gibson['incubation']}. Run a negative control with water in place of "
            "the insert — it tells you how much of the colony count is vector self-closure."
        ),
    )


def ligation_table(backbone: Component, insert: Component, cfg) -> Reaction:
    """Sticky-end ligation at a fixed backbone:insert molar ratio."""
    lig = cfg["ligation"]
    ratio = lig["insert_molar_excess"]
    backbone_fmol = lig["backbone_fmol"]

    rows = []
    for component, fmol, role in (
        (backbone, backbone_fmol, "backbone"),
        (insert, backbone_fmol * ratio, "insert"),
    ):
        rows.append(
            ReactionRow(
                f"{component.name} ({role}, {component.length_bp:,} bp)",
                ng=fmol_to_ng(fmol, component.length_bp, MW_PER_BP_LIGATION),
                fmol=fmol,
                ng_per_ul=component.ng_per_ul,
            )
        )
    rows.append(ReactionRow(f"10x {lig['enzyme']} buffer", fixed_ul=lig["buffer_ul"]))
    rows.append(ReactionRow(lig["enzyme"], fixed_ul=lig["enzyme_ul"]))
    return Reaction(
        title=f"Ligation — 1:{ratio} backbone:insert",
        rows=rows,
        total_ul=lig["total_ul"],
        note=(
            f"Leave at {lig['incubation']}. Include a no-insert control: it is the only way "
            "to read the colony count on the real plate."
        ),
    )


# ---------------------------------------------------------------------------
# PCR
# ---------------------------------------------------------------------------


def pcr_table(cfg) -> Table:
    """A single high-fidelity PCR reaction."""
    pcr = cfg["pcr"]
    volume = pcr["reaction_ul"]
    master_mix = volume / 2
    water = volume - master_mix - 2 * pcr["primer_ul"] - pcr["template_ul"]
    rows = [
        ["PCR-grade water", f"{water:.2f}"],
        [pcr["master_mix"], f"{master_mix:.2f}"],
        [f"{pcr['primer_um']} µM forward primer", f"{pcr['primer_ul']:.2f}"],
        [f"{pcr['primer_um']} µM reverse primer", f"{pcr['primer_ul']:.2f}"],
        [f"Template DNA ({pcr['template_ng_per_ul']} ng/µL)", f"{pcr['template_ul']:.2f}"],
        ["Total", f"{volume:.2f}"],
    ]
    return Table(
        title=f"PCR master mix — {volume:g} µL reaction",
        columns=["Reagent", "Volume (µL)"],
        rows=rows,
        note=(
            "Set up on ice: KAPA HiFi HotStart has enough proofreading activity to chew the "
            "primers at room temperature. Dilute the template to 1 ng/µL first."
        ),
    )


def thermocycler_table(anneal_c: float, amplicon_bp: int, cfg, repetitive=False) -> Table:
    """Cycling conditions, with the extension time set by amplicon length."""
    pcr = cfg["pcr"]
    per_kb = pcr["seconds_per_kb_repetitive"] if repetitive else pcr["seconds_per_kb"]
    extension = max(15, round(amplicon_bp / 1000 * per_kb))
    reason = (
        f"{per_kb} s/kb — the amplicon carries a repeat array, which the polymerase reads "
        "slowly; dropping to 30 s/kb risks a truncated product"
        if repetitive
        else f"{per_kb} s/kb"
    )
    rows = [
        ["Initial denaturation", "95 °C", "3 min", "1"],
        ["Denaturation", "98 °C", "20 s", ""],
        ["Annealing", f"{anneal_c:.0f} °C", "15 s", f"{pcr['cycles']}"],
        ["Extension", "72 °C", f"{extension} s", ""],
        ["Final extension", "72 °C", "2 min", "1"],
        ["Hold", "4 °C", "∞", ""],
    ]
    return Table(
        title="Thermocycler programme",
        columns=["Step", "Temperature", "Time", "Cycles"],
        rows=rows,
        note=f"Annealing is 2 °C below the lower primer Tm. Extension {reason}.",
    )


# ---------------------------------------------------------------------------
# transformation and screening
# ---------------------------------------------------------------------------


def growth_temperature(cfg, repetitive: bool) -> tuple[int, str]:
    """Outgrowth temperature, and why."""
    if repetitive:
        return (
            cfg["growth_temp_c_repetitive"],
            "tandem repeats recombine out at 37 °C, so use a recA- strain and keep the "
            "culture short",
        )
    return cfg["growth_temp_c"], "nothing repetitive in the product; standard growth"


def transformation_steps(cfg, repetitive: bool) -> list[str]:
    """Transformation, plating and colony picking as a numbered list."""
    transform = cfg["transformation"]
    temperature, reason = growth_temperature(cfg, repetitive)
    return [
        f"Thaw {transform['cells_ul']} µL of competent cells on ice and add "
        f"{transform['dna_ul']} µL of the assembly. Carry the negative control through "
        "alongside it.",
        f"Recover in {transform['soc_ul']} µL SOC, {transform['outgrowth']} at "
        f"{temperature} °C.",
        f"Plate {transform['plated_ul']} µL on selective agar and grow overnight at "
        f"{temperature} °C — {reason}.",
        f"Pick {transform['colonies']} colonies, miniprep, and send for "
        f"{transform['sequencing']}. Record the alignment and concentration below.",
    ]
