"""Lab-owned configuration: what is on the shelf, and what counts as a blank vector.

Kept as data rather than code because it changes with the freezer, not with the algorithm.
Resolution order: explicit path -> ``$BBL_LAB_CONFIG`` -> ``config/lab.json`` beside the repo
-> built-in defaults.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

DEFAULTS = {
    # Empty/backbone-only vectors: valid starting points when no construct is close enough.
    # A judgement about intent, so it cannot be inferred from sequence.
    "base_vectors": ["pHL162_pcDNA3.1_MCS"],
    # None means "assume anything commercially available is obtainable". A list restricts
    # planning to enzymes actually in the freezer.
    "enzyme_stock": None,
    # Plasmid labels present as files but not as tubes.
    "unavailable": [],
}

_SEARCH = ("config/lab.json", "../config/lab.json")


def load_lab_config(path=None) -> dict:
    """Load lab configuration, falling back to defaults for anything unspecified."""
    candidates = []
    if path:
        candidates.append(Path(path))
    if os.environ.get("BBL_LAB_CONFIG"):
        candidates.append(Path(os.environ["BBL_LAB_CONFIG"]))
    here = Path(__file__).resolve().parents[2]
    candidates += [here / name for name in _SEARCH]

    for candidate in candidates:
        if candidate.is_file():
            return {**DEFAULTS, **json.loads(candidate.read_text())}
    return dict(DEFAULTS)


def is_base_vector(label: str, config=None) -> bool:
    config = config or load_lab_config()
    return label in set(config.get("base_vectors") or [])


def is_available(label: str, config=None) -> bool:
    """A `.dna` file is a design archive, not proof of a tube in the freezer."""
    config = config or load_lab_config()
    return label not in set(config.get("unavailable") or [])
