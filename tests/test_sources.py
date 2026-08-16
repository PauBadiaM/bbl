"""Tests for storage independence.

Nothing in the package may assume a particular filesystem, a particular directory, or that
plasmids are files at all. These pin that: the same inventory, deletion and comparison work
from a directory, an explicit file list, in-memory records, and a synthetic remote store.
"""

import pytest

from bbl import load_plasmid
from bbl.inventory import scan_inventory
from bbl.llm import DesignSession
from bbl.sources import (
    CallableSource,
    DirectorySource,
    FileSource,
    MultiSource,
    PlasmidSource,
    RecordSource,
    as_source,
)

PCLM1 = "pCLM1_pcDNA3.1_miniCMV-mCherry-LambdaBoxBx8"
PHL391 = "pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8"


@pytest.fixture(scope="module")
def files(plasmid_dir):
    return sorted(plasmid_dir.glob("*.dna"))


# --------------------------------------------------------------------------- #
# every built-in source satisfies the protocol
# --------------------------------------------------------------------------- #


def test_all_sources_satisfy_the_protocol(plasmid_dir, files):
    records = {"one": load_plasmid(files[0])}
    sources = [
        DirectorySource(plasmid_dir),
        FileSource(files[:2]),
        RecordSource(records),
        CallableSource("remote", ids=lambda: ["x"], load=lambda i: records["one"]),
    ]
    for source in sources:
        assert isinstance(source, PlasmidSource)
        assert source.name
        for plasmid_id in source.ids()[:1]:
            assert source.label(plasmid_id)
            assert len(source.load(plasmid_id)) > 1000


def test_coercion_accepts_every_spec(plasmid_dir, files):
    assert isinstance(as_source(plasmid_dir), DirectorySource)
    assert isinstance(as_source(str(plasmid_dir)), DirectorySource)
    assert isinstance(as_source(files[:3]), FileSource)
    assert isinstance(as_source(files[0]), FileSource)          # a lone file
    assert isinstance(as_source({"a": load_plasmid(files[0])}), RecordSource)
    already = DirectorySource(plasmid_dir)
    assert as_source(already) is already


def test_coercion_failure_explains_the_options():
    with pytest.raises(FileNotFoundError, match="Benchling"):
        as_source("/nonexistent/path/xyz")


# --------------------------------------------------------------------------- #
# the same library, four ways
# --------------------------------------------------------------------------- #


def test_directory_and_file_list_agree(plasmid_dir, files):
    from_dir, _ = scan_inventory(plasmid_dir)
    from_files, _ = scan_inventory(files)
    assert [e.label for e in from_dir] == [e.label for e in from_files]
    assert [e.fingerprint for e in from_dir] == [e.fingerprint for e in from_files]


def test_in_memory_records_need_no_filesystem(files):
    records = {p.stem: load_plasmid(p) for p in files[:4]}
    entries, failures = scan_inventory(records)
    assert not failures and len(entries) == 4
    assert all(entry.path is None for entry in entries)   # nothing on disk
    assert all(entry.length > 1000 for entry in entries)


def test_a_synthetic_remote_store_behaves_like_a_directory(files):
    """Stands in for Benchling: identifiers are opaque, records arrive over a callable."""
    catalog = {f"seq_{i}": load_plasmid(p) for i, p in enumerate(files[:5])}
    names = {f"seq_{i}": p.stem for i, p in enumerate(files[:5])}
    calls = []

    remote = CallableSource(
        name="benchling:test-project",
        ids=lambda: list(catalog),
        load=lambda i: (calls.append(i), catalog[i])[1],
        label=lambda i: names[i],
    )
    entries, failures = scan_inventory(remote)
    assert not failures and len(entries) == 5
    assert {e.label for e in entries} == set(names.values())
    assert all(e.path is None for e in entries)
    assert all(e.source_id.startswith("seq_") for e in entries)


def test_remote_fetches_are_cached(files):
    calls = []
    record = load_plasmid(files[0])
    remote = CallableSource(
        "remote", ids=lambda: ["a"], load=lambda i: (calls.append(i), record)[1]
    )
    remote.load("a"), remote.load("a"), remote.load("a")
    assert len(calls) == 1, "a remote library must not refetch on every access"


def test_multisource_combines_libraries(plasmid_dir, files):
    combined = MultiSource([FileSource(files[:3]), FileSource(files[3:6])])
    assert len(combined.ids()) == 6
    entries, _ = scan_inventory(combined)
    assert len(entries) == 6


# --------------------------------------------------------------------------- #
# the design session is storage-agnostic end to end
# --------------------------------------------------------------------------- #


def test_full_design_from_in_memory_records(plasmid_dir):
    """The ground-truth deletion, with no directory anywhere in the call path."""
    records = {
        PHL391: load_plasmid(plasmid_dir / f"{PHL391}.dna"),
        PCLM1: load_plasmid(plasmid_dir / f"{PCLM1}.dna"),
    }
    session = DesignSession(RecordSource(records, name="in-memory"))
    assert len(session.entries) == 2

    result = session.plan_deletion("pHL391", ["NFKBRE"])
    assert set(result["enzymes"]) == {"MfeI", "EcoRI"}
    assert result["length_bp"] == 5704
    assert session.compare_product(result["product_id"], "pCLM1")["identical"]


def test_full_design_from_a_synthetic_remote_store(plasmid_dir):
    records = {
        "bnch_001": load_plasmid(plasmid_dir / f"{PHL391}.dna"),
        "bnch_002": load_plasmid(plasmid_dir / f"{PCLM1}.dna"),
    }
    labels = {"bnch_001": PHL391, "bnch_002": PCLM1}
    remote = CallableSource(
        "benchling:demo",
        ids=lambda: list(records),
        load=records.__getitem__,
        label=labels.__getitem__,
    )
    session = DesignSession(remote)
    result = session.plan_deletion("pHL391", ["NFKBRE"])
    assert result["length_bp"] == 5704
    assert session.compare_product(result["product_id"], "pCLM1")["identical"]
    assert "benchling:demo" in session.source.name


def test_lookup_error_names_the_source_not_a_directory(plasmid_dir):
    records = {PHL391: load_plasmid(plasmid_dir / f"{PHL391}.dna")}
    session = DesignSession(RecordSource(records, name="benchling:demo"))
    with pytest.raises(LookupError, match="benchling:demo"):
        session.resolve("pNOPE")


# --------------------------------------------------------------------------- #
# library resolution for the REPL
# --------------------------------------------------------------------------- #


def test_find_library_prefers_the_explicit_argument(plasmid_dir, monkeypatch):
    from bbl.llm.chat import find_library

    monkeypatch.delenv("BBL_PLASMID_DIR", raising=False)
    source = find_library(plasmid_dir)
    assert source is not None and len(source.ids()) == 28


def test_find_library_falls_back_to_the_env_var(plasmid_dir, monkeypatch):
    from bbl.llm.chat import find_library

    monkeypatch.setenv("BBL_PLASMID_DIR", str(plasmid_dir))
    assert len(find_library(None).ids()) == 28


def test_find_library_asks_when_it_cannot_tell(plasmid_dir, monkeypatch, tmp_path, capsys):
    from bbl.llm import chat

    monkeypatch.delenv("BBL_PLASMID_DIR", raising=False)
    monkeypatch.chdir(tmp_path)                      # no ./plasmid here
    answers = iter([str(tmp_path / "empty"), str(plasmid_dir)])
    (tmp_path / "empty").mkdir()
    monkeypatch.setattr(chat, "input", lambda _: next(answers), raising=False)

    source = chat.find_library(None)
    assert source is not None and len(source.ids()) == 28
    output = capsys.readouterr().out
    assert "Point me at your plasmids" in output
    assert "no plasmid files there" in output       # re-asked after the empty directory


def test_find_library_can_be_declined(monkeypatch, tmp_path):
    from bbl.llm import chat

    monkeypatch.delenv("BBL_PLASMID_DIR", raising=False)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(chat, "input", lambda _: "/quit", raising=False)
    assert chat.find_library(None) is None
