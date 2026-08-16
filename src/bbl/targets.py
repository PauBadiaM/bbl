"""Resolving which sequence to delete, and which annotations must survive.

Shared by every deletion route (restriction excision, around-the-horn PCR) so that
"delete NFKBRE" means exactly the same thing regardless of how it is built.
"""

from __future__ import annotations

from .plasmid_io import NON_FUNCTIONAL_TYPES, feature_label, feature_span

#: Labels that look like an origin or a selection marker; deleting one is flagged.
ESSENTIAL_HINTS = ("ori", "ampr", "neor", "kanr", "puro", "blast", "bsd", "hygro")


def tight_arc(intervals: list[tuple[int, int]], length: int) -> tuple[int, int]:
    """Smallest circular arc covering every interval.

    Computed as the complement of the largest gap between consecutive intervals, which
    handles targets that straddle the origin. ``end`` may exceed ``length`` to signal a wrap.
    """
    ordered = sorted(intervals)
    best_gap, best_index = -1, 0
    for i, (start, _) in enumerate(ordered):
        prev_end = ordered[i - 1][1]
        gap = (start - prev_end) % length
        if gap > best_gap:
            best_gap, best_index = gap, i
    start = ordered[best_index][0]
    end = ordered[best_index - 1][1]
    if end <= start:
        end += length
    return start, end


def resolve_target(record, features_to_remove, _rotated: bool = False):
    """Return ``(record, (start, end), labels)`` with the target not wrapping the origin.

    ``features_to_remove`` is a list of feature-label substrings or ``(start, end)`` pairs.
    If the target arc crosses the origin the record is rotated so it does not, which keeps
    every downstream coordinate comparison a plain interval test.

    Rotating remaps annotations automatically, but *coordinate* targets have to be remapped
    by hand -- otherwise the same wrapping span is re-derived on the rotated record and the
    rotation repeats forever.
    """
    if not features_to_remove:
        raise ValueError("features_to_remove is empty; nothing to delete")

    coordinates_given = all(
        isinstance(f, (tuple, list)) and len(f) == 2 for f in features_to_remove
    )
    if coordinates_given:
        intervals = [(int(s), int(e)) for s, e in features_to_remove]
        labels = [f"{s}..{e}" for s, e in intervals]
    else:
        wanted = [str(f).lower() for f in features_to_remove]
        matches = [
            f
            for f in record.features
            if f.type not in NON_FUNCTIONAL_TYPES
            and any(w in feature_label(f).lower() for w in wanted)
        ]
        if not matches:
            available = sorted({feature_label(f) for f in record.features})
            raise ValueError(
                f"No features matched {list(features_to_remove)}. Available: {available}"
            )
        intervals = [feature_span(f) for f in matches]
        labels = [feature_label(f) for f in matches]

    length = len(record)
    start, end = tight_arc(intervals, length)
    if end > length:  # wraps the origin -> rotate so it does not
        if _rotated:
            raise ValueError(
                f"target {start}..{end} still wraps the origin after rotation; "
                "give explicit coordinates that do not span it"
            )
        rotated = record.shifted(start)
        if coordinates_given:
            remapped = [
                ((s - start) % length, (s - start) % length + (e - s)) for s, e in intervals
            ]
            return resolve_target(rotated, remapped, _rotated=True)
        return resolve_target(rotated, features_to_remove, _rotated=True)
    return record, (start, end), sorted(set(labels))


def classify_features(record, span, protect=None, margin=0):
    """Split annotations into ``(removed, protected, umbrella)``.

    ``removed``   -- fully inside the target, so intentionally deleted.
    ``protected`` -- ``(start, end, label)`` triples, expanded by ``margin``, that must survive.
    ``umbrella``  -- features that *contain* the target (e.g. "Insert Sequence"); these are
                     expected to shrink and must not block the design.

    ``primer_bind``-style annotations are ignored entirely: they are bookkeeping, and treating
    them as protected makes correct designs unreachable (see docs/DECISIONS.md D20).
    """
    target_start, target_end = span
    protect_names = [p.lower() for p in (protect or [])]

    removed, protected, umbrella = [], [], []
    for feature in record.features:
        if feature.type in NON_FUNCTIONAL_TYPES:
            continue
        label = feature_label(feature)
        start, end = feature_span(feature)
        explicitly_protected = any(p in label.lower() for p in protect_names)

        if start >= target_start and end <= target_end:
            if explicitly_protected:
                raise ValueError(
                    f"{label!r} lies inside the region being deleted but was listed in "
                    "protect=; the request is contradictory"
                )
            if label not in removed:
                removed.append(label)
        elif start <= target_start and end >= target_end:
            umbrella.append(label)
        else:
            protected.append((start - margin, end + margin, label))
    return removed, protected, umbrella


def essential_warnings(removed_labels) -> list[str]:
    """Warnings for deleting anything that looks like an origin or selection marker."""
    return [
        f"{label!r} looks essential (origin or selection marker) and is being deleted; "
        "the product may not propagate"
        for label in removed_labels
        if any(hint in label.lower() for hint in ESSENTIAL_HINTS)
    ]
