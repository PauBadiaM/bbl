"""Loading, writing and comparing circular plasmid records.

SnapGene ``.dna`` files are read via Biopython's ``snapgene`` parser. That parser is
**read-only**, so products are written as GenBank, which SnapGene opens natively.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord

_SUFFIX_FORMAT = {
    ".dna": "snapgene",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".genbank": "genbank",
    ".fa": "fasta",
    ".fasta": "fasta",
}

#: Feature types that are annotation bookkeeping rather than functional elements.
#: These are never treated as "protected" when planning a deletion -- see docs/DECISIONS.md D20.
NON_FUNCTIONAL_TYPES = frozenset({"primer_bind", "primer", "source"})


def load_plasmid(source, circular: bool = True) -> Dseqrecord:
    """Load a plasmid as a circular :class:`Dseqrecord`.

    ``source`` may be a path to a ``.dna``/``.gb``/``.fasta`` file, a ``SeqRecord``,
    or an already-built ``Dseqrecord`` (returned as-is if already circular).
    """
    if isinstance(source, Dseqrecord):
        return source if source.circular == circular else Dseqrecord(source, circular=circular)
    if isinstance(source, SeqRecord):
        return Dseqrecord(source, circular=circular)

    path = Path(source)
    fmt = _SUFFIX_FORMAT.get(path.suffix.lower())
    if fmt is None:
        raise ValueError(
            f"Unrecognised plasmid file type {path.suffix!r}; expected one of "
            f"{sorted(_SUFFIX_FORMAT)}"
        )
    record = SeqIO.read(str(path), fmt)
    record.annotations.setdefault("molecule_type", "ds-DNA")
    record.annotations["topology"] = "circular" if circular else "linear"
    plasmid = Dseqrecord(record, circular=circular)
    # SnapGene files routinely carry junk in the name field (a GenScript URL, "GenBank",
    # "<unknown name>"), which then shows up in protocols. Fall back to the file stem.
    name = (plasmid.name or "").strip()
    if not name or name == "<unknown name>" or not name.replace("_", "").replace("-", "").isalnum():
        name = path.stem
    plasmid.name = name.split("_")[0][:16] if name.startswith(path.stem.split("_")[0]) else name[:16]
    return plasmid


def write_genbank(record, path) -> Path:
    """Write ``record`` as GenBank. Returns the path written."""
    path = Path(path)
    record.annotations.setdefault("molecule_type", "ds-DNA")
    path.write_text(record.format("gb"))
    return path


def feature_label(feature) -> str:
    """Best human-readable name for a feature."""
    for key in ("label", "note", "gene", "product"):
        if key in feature.qualifiers:
            value = feature.qualifiers[key]
            return (value[0] if isinstance(value, (list, tuple)) else value).strip()
    return feature.type


def feature_span(feature) -> tuple[int, int]:
    """0-based half-open ``(start, end)`` of a feature."""
    return int(feature.location.start), int(feature.location.end)


def find_features(record, names) -> list:
    """Features whose label matches any of ``names`` (case-insensitive substring)."""
    wanted = [n.lower() for n in names]
    return [
        f
        for f in record.features
        if f.type not in NON_FUNCTIONAL_TYPES
        and any(w in feature_label(f).lower() for w in wanted)
    ]


def slice_circular(sequence, start: int, end: int) -> str:
    """``sequence[start:end]`` with wrap-around; ``start`` may be negative."""
    text = str(getattr(sequence, "seq", sequence))
    length = len(text)
    return (text * 3)[length + start : length + start + (end - start)]


def replace_span(record, start: int, end: int, sequence="", features=None,
                 truncation_note: str = "truncated by edit"):
    """Return ``record`` with ``[start:end)`` replaced by ``sequence``, in the parent's frame.

    ``features`` are extra ``SeqFeature`` objects whose coordinates are relative to
    ``sequence``; they are shifted into place. Existing features entirely inside the replaced
    span are dropped, those straddling it are truncated (and noted) -- what SnapGene does.
    Keeping the parent's origin means product coordinates stay comparable to the parent's.

    Setting ``start == end`` is a pure insertion; ``sequence == ""`` is a pure deletion.
    """
    # Slice the plain string, never the Dseq: for a *circular* Dseq, ``seq[:0]`` returns the
    # whole circle rather than an empty string, which silently doubles the product whenever
    # the edit starts at position 0.
    text = str(record.seq)
    sequence = str(sequence)
    delta = len(sequence) - (end - start)

    # start and end coordinates break the tie differently, so that a feature butting up
    # against a pure insertion point stays put while one beginning there is pushed along.
    def clamp_start(position):
        if position < start:
            return position
        if position >= end:
            return position + delta
        return start

    def clamp_end(position):
        if position <= start:
            return position
        if position >= end:
            return position + delta
        return start

    new_record = SeqRecord(
        Seq(text[:start] + sequence + text[end:]),
        id=record.id,
        name=record.name,
        description=record.description,
        annotations=dict(record.annotations),
    )
    new_record.annotations["molecule_type"] = "ds-DNA"
    new_record.annotations["topology"] = "circular"

    for feature in record.features:
        parts, truncated = [], False
        for part in feature.location.parts:
            part_start, part_end = int(part.start), int(part.end)
            new_start, new_end = clamp_start(part_start), clamp_end(part_end)
            if new_end <= new_start:
                continue  # entirely inside the replaced span
            if (new_end - new_start) != (part_end - part_start):
                truncated = True
            parts.append(SimpleLocation(new_start, new_end, strand=part.strand))
        if not parts:
            continue
        if truncated and feature.type in NON_FUNCTIONAL_TYPES:
            continue  # a clipped primer no longer anneals; drop rather than mislead
        location = parts[0] if len(parts) == 1 else CompoundLocation(parts)
        qualifiers = {
            key: list(value) if isinstance(value, list) else value
            for key, value in feature.qualifiers.items()
        }
        if truncated:
            qualifiers.setdefault("note", []).append(truncation_note)
        new_record.features.append(
            SeqFeature(location=location, type=feature.type, qualifiers=qualifiers)
        )

    for feature in features or []:
        shifted = [
            SimpleLocation(int(part.start) + start, int(part.end) + start, strand=part.strand)
            for part in feature.location.parts
        ]
        new_record.features.append(
            SeqFeature(
                location=shifted[0] if len(shifted) == 1 else CompoundLocation(shifted),
                type=feature.type,
                qualifiers=dict(feature.qualifiers),
            )
        )
    return Dseqrecord(new_record, circular=True)


def delete_span(record, start: int, end: int):
    """Return ``record`` with ``[start:end)`` removed, in the parent's coordinate frame."""
    return replace_span(record, start, end, "", truncation_note="truncated by deletion")


def features_in_span(record, start: int, end: int) -> list:
    """Features fully inside ``[start:end)``, rebased so coordinates start at 0."""
    extracted = []
    for feature in record.features:
        if feature.type in NON_FUNCTIONAL_TYPES:
            continue
        f_start, f_end = int(feature.location.start), int(feature.location.end)
        if f_start < start or f_end > end:
            continue
        parts = [
            SimpleLocation(int(p.start) - start, int(p.end) - start, strand=p.strand)
            for p in feature.location.parts
        ]
        extracted.append(
            SeqFeature(
                location=parts[0] if len(parts) == 1 else CompoundLocation(parts),
                type=feature.type,
                qualifiers=dict(feature.qualifiers),
            )
        )
    return extracted


def rotation_offset(product, expected) -> int | None:
    """Offset at which ``expected`` occurs in the doubled ``product``, else ``None``.

    Circular records have an arbitrary origin, so plasmids that are biologically
    identical can differ by a rotation. Compares the forward strand only.
    """
    a = str(getattr(product, "seq", product)).upper()
    b = str(getattr(expected, "seq", expected)).upper()
    if len(a) != len(b):
        return None
    index = (a + a).find(b)
    return None if index < 0 else index


def circular_equal(product, expected) -> bool:
    """True if two circular sequences are identical up to rotation or reverse-complement."""
    if rotation_offset(product, expected) is not None:
        return True
    rc = str(Seq(str(getattr(expected, "seq", expected)).upper()).reverse_complement())
    return rotation_offset(product, rc) is not None
