"""Where plasmids come from.

Nothing in the package assumes a particular filesystem, a particular directory, or that the
plasmids are files at all. A source answers three questions -- what have you got, what is it
called, and give me the record -- and everything else is built on that.

Built in:

* :class:`DirectorySource` -- a folder of ``.dna``/``.gb`` files, wherever it lives.
* :class:`FileSource` -- an explicit list of files (uploads, a hand-picked subset).
* :class:`RecordSource` -- records already in memory, for callers that parsed their own.
* :class:`CallableSource` -- a list function plus a fetch function, for anything remote.

**Benchling and other remote stores** plug in as a :class:`CallableSource`::

    source = CallableSource(
        name="benchling:my-project",
        ids=lambda: [e["id"] for e in benchling.dna_sequences.list(folder_id=FOLDER)],
        load=lambda seq_id: benchling_to_record(benchling.dna_sequences.get(seq_id)),
        label=lambda seq_id: name_cache[seq_id],
    )

``load`` must return something :func:`bbl.plasmid_io.load_plasmid` accepts -- a ``SeqRecord``,
a ``Dseqrecord``, or a path. Fetches are lazy and cached per source, so a large remote library
costs one call per plasmid actually touched.
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol, runtime_checkable

from .plasmid_io import _SUFFIX_FORMAT, load_plasmid


@runtime_checkable
class PlasmidSource(Protocol):
    """Minimal contract every source satisfies."""

    name: str

    def ids(self) -> list[str]:
        """Stable identifiers for everything in this source."""

    def label(self, plasmid_id: str) -> str:
        """Human-readable name, used for matching and display."""

    def load(self, plasmid_id: str):
        """The plasmid itself, as anything ``load_plasmid`` accepts."""


class _Cached:
    """Mixin: fetch each plasmid at most once per source."""

    def __init__(self):
        self._cache: dict[str, object] = {}

    def load(self, plasmid_id: str):
        if plasmid_id not in self._cache:
            self._cache[plasmid_id] = load_plasmid(self._fetch(plasmid_id))
        return self._cache[plasmid_id]


class DirectorySource(_Cached):
    """Every plasmid file in a directory. The ordinary local or server case."""

    def __init__(self, directory, pattern: str = "*.dna", recursive: bool = False):
        super().__init__()
        self.directory = Path(directory).expanduser()
        self.pattern = pattern
        self.recursive = recursive
        self.name = str(self.directory)

    def ids(self) -> list[str]:
        glob = self.directory.rglob if self.recursive else self.directory.glob
        return sorted(str(p) for p in glob(self.pattern))

    def label(self, plasmid_id: str) -> str:
        return Path(plasmid_id).stem

    def _fetch(self, plasmid_id: str):
        return Path(plasmid_id)


class FileSource(_Cached):
    """An explicit list of files -- uploads, or a hand-picked subset of a bigger folder."""

    def __init__(self, paths, name: str = "files"):
        super().__init__()
        self.paths = [Path(p).expanduser() for p in paths]
        self.name = name

    def ids(self) -> list[str]:
        return [str(p) for p in self.paths]

    def label(self, plasmid_id: str) -> str:
        return Path(plasmid_id).stem

    def _fetch(self, plasmid_id: str):
        return Path(plasmid_id)


class RecordSource:
    """Records already in memory, keyed by name. No filesystem involved."""

    def __init__(self, records: dict, name: str = "memory"):
        self.records = dict(records)
        self.name = name

    def ids(self) -> list[str]:
        return list(self.records)

    def label(self, plasmid_id: str) -> str:
        return plasmid_id

    def load(self, plasmid_id: str):
        return load_plasmid(self.records[plasmid_id])


class CallableSource(_Cached):
    """A remote store, described by functions. This is the Benchling shape.

    Args:
        name: Identifies the source in reports and errors, e.g. ``"benchling:my-project"``.
        ids: Returns the identifiers available.
        load: Given an identifier, returns a record (or a path) for it.
        label: Given an identifier, returns its display name. Defaults to the identifier.
    """

    def __init__(self, name: str, ids, load, label=None):
        super().__init__()
        self.name = name
        self._ids = ids
        self._load = load
        self._label = label or (lambda plasmid_id: plasmid_id)

    def ids(self) -> list[str]:
        return list(self._ids())

    def label(self, plasmid_id: str) -> str:
        return self._label(plasmid_id)

    def _fetch(self, plasmid_id: str):
        return self._load(plasmid_id)


class MultiSource:
    """Several sources as one library -- e.g. a local folder plus a shared drive."""

    def __init__(self, sources, name: str = "combined"):
        self.sources = [as_source(s) for s in sources]
        self.name = name
        self._owner: dict[str, PlasmidSource] = {}
        for source in self.sources:
            for plasmid_id in source.ids():
                self._owner.setdefault(plasmid_id, source)

    def ids(self) -> list[str]:
        return list(self._owner)

    def label(self, plasmid_id: str) -> str:
        return self._owner[plasmid_id].label(plasmid_id)

    def load(self, plasmid_id: str):
        return self._owner[plasmid_id].load(plasmid_id)


def as_source(spec, pattern: str = "*.dna") -> PlasmidSource:
    """Coerce whatever the caller passed into a source.

    Accepts a source (returned as-is), a directory path, a list of file paths, or a mapping of
    name to record. This is what lets every public function keep taking a plain directory while
    the machinery underneath is storage-agnostic.
    """
    if isinstance(spec, PlasmidSource) and hasattr(spec, "ids"):
        return spec
    if isinstance(spec, dict):
        return RecordSource(spec)
    if isinstance(spec, (list, tuple, set)):
        items = list(spec)
        if items and all(isinstance(i, (str, Path)) and Path(i).suffix.lower() in _SUFFIX_FORMAT
                         for i in items):
            return FileSource(items)
        return MultiSource(items)
    path = Path(spec).expanduser()
    if path.is_dir():
        return DirectorySource(path, pattern)
    if path.exists():
        return FileSource([path])
    raise FileNotFoundError(
        f"{spec!r} is not a source, a directory, or a file. Pass a directory of .dna files, "
        "a list of files, or a PlasmidSource (see bbl.sources for the Benchling shape)."
    )
