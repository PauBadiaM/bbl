"""Sequence-identity layer of the plasmid inventory.

Filenames and feature labels are unreliable identifiers (see docs/DECISIONS.md D41), so the
library is keyed by **sequence**. Two plasmids are the same molecule if their sequences match
under any rotation and either orientation, which is what :func:`fingerprint` normalises away.

Near-duplicates are found with canonical k-mer sets: ``containment`` answers "is A essentially
contained in B" (the usual relationship between a parent and a derivative that gained an
insert), while ``jaccard`` answers "are these the same construct".
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

from Bio.Seq import Seq

from .plasmid_io import load_plasmid
from .sources import as_source

DEFAULT_K = 31


def least_rotation(text: str) -> int:
    """Index of the lexicographically smallest rotation (Booth's algorithm, O(n)).

    A plasmid has no natural origin, so the smallest rotation gives every circular sequence a
    single canonical starting point that survives however the file happened to be saved.
    """
    doubled = text + text
    n = len(doubled)
    failure = [-1] * n
    k = 0
    for j in range(1, n):
        char = doubled[j]
        i = failure[j - k - 1]
        while i != -1 and char != doubled[k + i + 1]:
            if char < doubled[k + i + 1]:
                k = j - i - 1
            i = failure[i]
        if char != doubled[k + i + 1]:
            if char < doubled[k]:
                k = j
            failure[j - k] = -1
        else:
            failure[j - k] = i + 1
    return k % len(text)


def canonical_sequence(record) -> str:
    """Rotation- and orientation-independent form of a circular sequence."""
    text = str(getattr(record, "seq", record)).upper()
    if not text:
        return ""
    forward = text[least_rotation(text) :] + text[: least_rotation(text)]
    reverse = str(Seq(text).reverse_complement())
    reverse = reverse[least_rotation(reverse) :] + reverse[: least_rotation(reverse)]
    return min(forward, reverse)


def fingerprint(record) -> str:
    """Short stable identity for a circular plasmid."""
    return hashlib.sha256(canonical_sequence(record).encode()).hexdigest()[:16]


def canonical_kmers(record, k: int = DEFAULT_K) -> set[str]:
    """Strand-independent k-mer set of a circular sequence."""
    text = str(getattr(record, "seq", record)).upper()
    if len(text) < k:
        return set()
    circular = text + text[: k - 1]
    kmers = set()
    for i in range(len(text)):
        kmer = circular[i : i + k]
        rc = str(Seq(kmer).reverse_complement())
        kmers.add(min(kmer, rc))
    return kmers


@dataclass
class InventoryEntry:
    """One plasmid in the library, independent of where it is stored."""

    label: str
    name: str
    length: int
    fingerprint: str
    kmers: set = field(repr=False, default_factory=set)
    sequence: str = field(repr=False, default="")
    source: object = field(repr=False, default=None)
    source_id: str = ""

    def load(self):
        """The record. Fetched through the source, so remote stores work unchanged."""
        return self.source.load(self.source_id)

    @property
    def path(self) -> Path | None:
        """Filesystem location, when the source has one. ``None`` for remote sources."""
        candidate = Path(self.source_id) if self.source_id else None
        return candidate if candidate and candidate.exists() else None


@dataclass
class Provenance:
    """Where a sequence already exists in the library."""

    entry: InventoryEntry
    start: int
    strand: int  # +1 forward, -1 on the reverse strand

    @property
    def label(self) -> str:
        return self.entry.label


@dataclass
class Relationship:
    a: InventoryEntry
    b: InventoryEntry
    jaccard: float
    containment: float  # fraction of the smaller entry's k-mers present in the larger

    @property
    def kind(self) -> str:
        if self.jaccard >= 0.99:
            return "identical"
        if self.containment >= 0.98:
            return "subset"  # one is the other plus an insert
        if self.jaccard >= 0.80:
            return "variant"
        return "related"


def scan_inventory(source, pattern: str = "*.dna", k: int = DEFAULT_K):
    """Index a library. Returns ``(entries, failures)``.

    ``source`` is anything :func:`bbl.sources.as_source` accepts: a directory, a list of files,
    a mapping of records, or a :class:`~bbl.sources.PlasmidSource` such as a Benchling adapter.
    """
    source = as_source(source, pattern)
    entries, failures = [], []
    for plasmid_id in source.ids():
        try:
            record = source.load(plasmid_id)
        except Exception as exc:  # one unreadable plasmid should not stop the scan
            failures.append((plasmid_id, str(exc)))
            continue
        entries.append(
            InventoryEntry(
                label=source.label(plasmid_id),
                name=record.name,
                length=len(record),
                fingerprint=fingerprint(record),
                kmers=canonical_kmers(record, k),
                sequence=str(record.seq).upper(),
                source=source,
                source_id=plasmid_id,
            )
        )
    entries.sort(key=lambda e: e.label)
    return entries, failures


def find_sequence(entries, query: str) -> list[Provenance]:
    """Every exact occurrence of ``query`` in the library, either strand, wrapping the origin.

    This is the "do we already own this part?" lookup, and it must be sequence-based: the
    library annotates the same part under different labels and sometimes not at all.
    """
    query = str(query).upper()
    if not query:
        return []
    reverse = str(Seq(query).reverse_complement())
    hits = []
    for entry in entries:
        text = entry.sequence
        if not text or len(query) > len(text):
            continue
        extended = text + text[: len(query) - 1]
        for probe, strand in ((query, 1), (reverse, -1)):
            if strand == -1 and probe == query:
                continue  # palindrome: do not report it twice
            start = extended.find(probe)
            while start != -1:
                hits.append(Provenance(entry, start % len(text), strand))
                start = extended.find(probe, start + 1)
    return hits


def duplicate_groups(entries) -> dict[str, list[InventoryEntry]]:
    """Fingerprint -> entries, for fingerprints seen more than once."""
    groups: dict[str, list[InventoryEntry]] = {}
    for entry in entries:
        groups.setdefault(entry.fingerprint, []).append(entry)
    return {fp: group for fp, group in groups.items() if len(group) > 1}


def relationships(entries, min_containment: float = 0.5) -> list[Relationship]:
    """Pairs sharing enough sequence to be worth reporting, strongest first."""
    found = []
    for i, a in enumerate(entries):
        for b in entries[i + 1 :]:
            if not a.kmers or not b.kmers:
                continue
            shared = len(a.kmers & b.kmers)
            if not shared:
                continue
            union = len(a.kmers | b.kmers)
            smaller = min(len(a.kmers), len(b.kmers))
            containment = shared / smaller
            if containment < min_containment:
                continue
            found.append(Relationship(a, b, shared / union, containment))
    found.sort(key=lambda r: (-r.containment, -r.jaccard))
    return found


def lineage(entries, min_containment: float = 0.98) -> list[tuple[InventoryEntry, InventoryEntry]]:
    """Pairs related closely enough that one was plausibly built from the other.

    Returned smaller-first purely for stable ordering -- this is **not** a claim about which
    is the parent. Size cannot decide direction: both known derivations in this library are
    deletions, so the derivative is the *smaller* member in each case (pCLM1 from pHL391,
    pCLM3 from pCLM2). Direction needs external evidence -- file dates, lab records, or the
    experimentalist. See docs/DECISIONS.md D55.

    Every plasmid here shares a backbone, so containment has a floor near 0.5; only the top of
    the range is informative.
    """
    edges = []
    for rel in relationships(entries, min_containment):
        if rel.jaccard >= 0.99:
            continue  # same construct, not a derivation
        edges.append(tuple(sorted((rel.a, rel.b), key=lambda e: e.length)))
    return edges


def report(source, pattern: str = "*.dna", min_containment: float = 0.5) -> str:
    """Human-readable summary of what is actually distinct in the library."""
    resolved = as_source(source, pattern)
    entries, failures = scan_inventory(resolved, pattern)
    lines = [f"{len(entries)} plasmids loaded from {resolved.name}"]
    if failures:
        lines.append(f"{len(failures)} failed to load:")
        lines += [f"    {plasmid_id}: {why}" for plasmid_id, why in failures]

    distinct = {entry.fingerprint for entry in entries}
    lines.append(f"{len(distinct)} distinct sequences ({len(entries) - len(distinct)} redundant)")

    groups = duplicate_groups(entries)
    if groups:
        lines.append("\nEXACT DUPLICATES (identical under rotation/orientation)")
        for fp, group in groups.items():
            lines.append(f"  {fp}  {group[0].length} bp")
            for entry in group:
                lines.append(f"      {entry.label}")

    related = [r for r in relationships(entries, min_containment) if r.kind != "identical"]
    if related:
        lines.append("\nNEAR-DUPLICATES")
        for rel in related:
            lines.append(
                f"  [{rel.kind:8}] containment {rel.containment:.3f}  jaccard {rel.jaccard:.3f}"
            )
            lines.append(
                f"      {rel.a.label} ({rel.a.length} bp)\n      {rel.b.label} ({rel.b.length} bp)"
            )
    return "\n".join(lines)


#: Backwards-friendly short alias.
scan = scan_inventory
