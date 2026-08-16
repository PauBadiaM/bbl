"""The check table: one row per quantity, declaring how each backend computes it.

Deliberately a flat table rather than a class hierarchy -- there are ~9 checks and the only
thing that varies is three small functions, so inheritance would add indirection and no
leverage.

Coverage is intentionally asymmetric, and that asymmetry *is* the result being measured:

``both``
    ``gc_content``, ``max_homopolymer``, ``sequence_length``, ``specific_kmer``.
    Real comparisons -- the measurements must agree.

``native-only``
    ``restriction_sites`` (Proto's alphabet is ``frozenset("ACGT")``, so ``AccI``'s ``GTMKAC``
    is not expressible), ``longest_repeat`` and ``repeat_fraction`` (no Proto equivalent, and
    these are what catch the 8x BoxB array that vendors refuse).

``proto-only``
    ``kmer_frequency`` and ``dinucleotide_composition`` -- capability ``bbl`` does not have
    today. These are the checks that would justify taking the dependency.

Thresholds default to :data:`bbl.complexity.DEFAULT_LIMITS` wherever an equivalent exists, so
the QC layer and the incumbent synthesisability screen cannot drift apart.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from ..complexity import DEFAULT_LIMITS
from . import native

#: Per-check tolerance when comparing measurements across backends. Floats computed by two
#: different code paths will not be bit-identical; counts and lengths should be exact.
TOLERANCES = {
    "gc_content": 1e-6,
    "max_homopolymer": 0.0,
    "sequence_length": 0.0,
    "specific_kmer": 1e-9,
    "restriction_sites": 0.0,
    "longest_repeat": 0.0,
    "repeat_fraction": 1e-9,
    "kmer_frequency": 1e-9,
    "dinucleotide_composition": 1e-9,
}


def tolerance_for(check: str) -> float:
    return TOLERANCES.get(check, 1e-9)


@dataclass(frozen=True)
class Check:
    """One quantity, and how to obtain it from each backend.

    ``native`` takes ``(sequence, params)`` and returns ``(measurement, penalty, detail)``.
    ``proto_config`` maps our params onto the Proto config model's kwargs.
    ``proto_measure`` pulls the measurement out of Proto's metadata, which for the k-mer
    checks has *dynamic* keys (e.g. ``"GAATTC_frequency"``), hence a callable.
    """

    name: str
    what: str
    defaults: dict
    native: Callable[[str, dict], tuple[float | None, float, dict]] | None = None
    proto_key: str | None = None
    proto_config: Callable[[dict], dict] | None = None
    proto_measure: Callable[[dict, dict], float | None] | None = None
    tags: tuple[str, ...] = field(default_factory=tuple)


def _max_kmer_frequency(metadata: dict, params: dict) -> float | None:
    """Highest observed k-mer frequency, reduced from Proto's per-k-mer dict.

    Proto reports a mapping over *observed* k-mers only; the maximum is the scalar that
    actually corresponds to "is any k-mer over-represented".
    """
    key = f"{params['k']}mer_frequencies"
    table = metadata.get(key)
    if not isinstance(table, dict) or not table:
        return None
    return float(max(table.values()))


CHECKS: list[Check] = [
    Check(
        name="gc_content",
        what="GC percent",
        defaults={"min_gc": DEFAULT_LIMITS["min_gc"], "max_gc": DEFAULT_LIMITS["max_gc"]},
        native=native.gc_content,
        proto_key="gc-content",
        proto_config=lambda p: {"min_gc": p["min_gc"], "max_gc": p["max_gc"]},
        proto_measure=lambda md, p: (
            float(md["gc_content"]) if md.get("gc_content") is not None else None
        ),
        tags=("composition",),
    ),
    Check(
        name="max_homopolymer",
        what="longest single-base run",
        defaults={"max_length": DEFAULT_LIMITS["max_homopolymer"]},
        native=native.max_homopolymer,
        proto_key="max-homopolymer",
        proto_config=lambda p: {"max_length": p["max_length"]},
        proto_measure=lambda md, p: (
            float(md["max_homopolymer_length"])
            if md.get("max_homopolymer_length") is not None
            else None
        ),
        tags=("composition", "synthesis"),
    ),
    Check(
        name="sequence_length",
        what="fragment length in bp",
        defaults={
            "min_length": DEFAULT_LIMITS["min_length"],
            "max_length": DEFAULT_LIMITS["max_length"],
        },
        native=native.sequence_length,
        proto_key="sequence-length",
        proto_config=lambda p: {"min_length": p["min_length"], "max_length": p["max_length"]},
        proto_measure=lambda md, p: float(md["length"]) if md.get("length") is not None else None,
        tags=("synthesis",),
    ),
    Check(
        name="specific_kmer",
        what="frequency of one literal k-mer",
        # GAATTC == EcoRI, chosen because it is unambiguous and therefore the one restriction
        # site Proto *can* express. Compare with the restriction_sites row below.
        defaults={"kmer": "GAATTC", "min_value": 0.0, "max_value": 1.0},
        native=native.specific_kmer,
        proto_key="specific-kmer-frequency",
        proto_config=lambda p: {
            "kmer": p["kmer"],
            "min_value": p["min_value"],
            "max_value": p["max_value"],
            "scoring_mode": "frequency",
        },
        proto_measure=lambda md, p: (
            float(md[f"{p['kmer'].upper()}_frequency"])
            if md.get(f"{p['kmer'].upper()}_frequency") is not None
            else None
        ),
        tags=("cloning",),
    ),
    # ---- native only -------------------------------------------------------------------
    Check(
        name="restriction_sites",
        what="cut-site count for one enzyme (ambiguity-aware)",
        defaults={"enzyme": "AccI", "min_sites": 0, "max_sites": 0, "circular": False},
        native=native.restriction_sites,
        proto_key=None,  # GT^MKAC is not expressible over frozenset("ACGT")
        tags=("cloning",),
    ),
    Check(
        name="longest_repeat",
        what="longest substring occurring more than once",
        defaults={"max_repeat": DEFAULT_LIMITS["max_repeat"]},
        native=native.longest_repeat_check,
        proto_key=None,
        tags=("synthesis",),
    ),
    Check(
        name="repeat_fraction",
        what="fraction of positions inside a repeated window",
        defaults={"max_repeat_fraction": DEFAULT_LIMITS["max_repeat_fraction"], "window": 20},
        native=native.repeat_fraction_check,
        proto_key=None,
        tags=("synthesis",),
    ),
    # ---- proto only --------------------------------------------------------------------
    Check(
        name="kmer_frequency",
        what="highest observed k-mer frequency",
        defaults={"k": 3, "min_value": 0.0, "max_value": 0.25},
        native=None,
        proto_key="kmer-frequency",
        proto_config=lambda p: {
            "k": p["k"],
            "min_value": p["min_value"],
            "max_value": p["max_value"],
            "scoring_mode": "frequency",
        },
        proto_measure=_max_kmer_frequency,
        tags=("composition",),
    ),
    Check(
        name="dinucleotide_composition",
        what="distance from a reference dinucleotide profile",
        # Uniform reference (1/16 each) is a neutral default; the lab's own profile can be
        # substituted through params without touching this table.
        defaults={
            "reference_frequencies": {
                a + b: 1.0 / 16.0 for a in "ACGT" for b in "ACGT"
            },
            "distance_metric": "total_variation",
            "scale": 1.0,
        },
        native=None,
        proto_key="dinucleotide-composition",
        proto_config=lambda p: {
            "reference_frequencies": p["reference_frequencies"],
            "distance_metric": p["distance_metric"],
            "scale": p["scale"],
        },
        proto_measure=lambda md, p: (
            float(md["dinucleotide_distance"])
            if md.get("dinucleotide_distance") is not None
            else None
        ),
        tags=("composition",),
    ),
]

CHECKS_BY_NAME = {check.name: check for check in CHECKS}


def resolve_params(check: Check, overrides: dict | None = None) -> dict:
    """Check defaults merged with per-call overrides."""
    params = dict(check.defaults)
    if overrides:
        params.update({k: v for k, v in overrides.items() if k in params or True})
    return params
