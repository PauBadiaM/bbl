"""Tests for the sequence-identity layer of the inventory.

Filenames and feature labels are unreliable, so the library is keyed by sequence. These tests
pin the invariants that make that safe, and check the index rediscovers the two cloning
relationships we validated independently in test_excise.py and test_insert.py.
"""

from Bio.Seq import Seq

from bbl.inventory import (
    canonical_kmers,
    canonical_sequence,
    duplicate_groups,
    fingerprint,
    least_rotation,
    lineage,
    relationships,
    scan,
)

PHL391 = "pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8"
PCLM1 = "pCLM1_pcDNA3.1_miniCMV-mCherry-LambdaBoxBx8"
PCLM2 = "pCLM2_pcDNA3.1_CMV-mCherry-LambdaBoxBx8"
PCLM3 = "pCLM3_pcDNA3.1_CMV-mCherry"

SEQ = "ATGGCCTTAAGGCATTACGATCCGGATTACAGGCATTTACGGCATAGCATTGACCA"


# --------------------------------------------------------------------------- #
# fingerprint invariants
# --------------------------------------------------------------------------- #


def test_least_rotation():
    assert least_rotation("bbaa") == 2
    assert least_rotation("aaaa") == 0
    assert least_rotation("abab") == 0
    rotated = "cab"
    assert rotated[least_rotation(rotated) :] + rotated[: least_rotation(rotated)] == "abc"


def test_fingerprint_is_rotation_invariant():
    """A plasmid has no natural origin, so every rotation is the same molecule."""
    prints = {fingerprint(SEQ[i:] + SEQ[:i]) for i in range(len(SEQ))}
    assert len(prints) == 1


def test_fingerprint_is_orientation_invariant():
    assert fingerprint(str(Seq(SEQ).reverse_complement())) == fingerprint(SEQ)


def test_fingerprint_distinguishes_real_differences():
    mutated = SEQ[:-1] + ("C" if SEQ[-1] != "C" else "G")
    assert fingerprint(mutated) != fingerprint(SEQ)
    assert fingerprint(SEQ + "A") != fingerprint(SEQ)


def test_canonical_sequence_preserves_length():
    assert len(canonical_sequence(SEQ)) == len(SEQ)


def test_canonical_kmers_are_strand_independent():
    assert canonical_kmers(SEQ, k=11) == canonical_kmers(
        str(Seq(SEQ).reverse_complement()), k=11
    )


def test_kmers_wrap_the_origin():
    """A k-mer spanning the origin must be found, or circular identity breaks."""
    kmers = canonical_kmers(SEQ, k=11)
    assert len(kmers) == len(SEQ)  # one k-mer per start position on a circle


# --------------------------------------------------------------------------- #
# the real inventory
# --------------------------------------------------------------------------- #


def test_inventory_loads_cleanly(plasmid_dir):
    entries, failures = scan(plasmid_dir)
    assert not failures, failures
    assert len(entries) >= 25
    assert all(entry.length > 1000 for entry in entries)


def test_no_two_files_hold_the_same_sequence(plasmid_dir):
    """Checked because two filenames carry sync-conflict artifacts; they are not duplicates."""
    entries, _ = scan(plasmid_dir)
    assert duplicate_groups(entries) == {}


def test_rediscovers_the_deletion_lineage(plasmid_dir):
    """pCLM1 was made from pHL391; it should fall out as a containment relationship."""
    entries, _ = scan(plasmid_dir)
    by_label = {e.label: e for e in entries}
    pairs = {
        frozenset((parent.label, child.label)) for parent, child in lineage(entries)
    }
    assert frozenset((PCLM1, PHL391)) in pairs


def test_rediscovers_the_insertion_lineage(plasmid_dir):
    """pCLM2 = pCLM3 + the BoxB array, so pCLM3 is contained in pCLM2."""
    entries, _ = scan(plasmid_dir)
    edges = {(parent.label, child.label) for parent, child in lineage(entries)}
    assert (PCLM3, PCLM2) in edges  # parent is the smaller one


def test_known_pairs_rank_at_the_top(plasmid_dir):
    """The two validated cloning steps should be among the strongest relationships."""
    entries, _ = scan(plasmid_dir)
    ranked = relationships(entries, min_containment=0.9)
    top = [frozenset((r.a.label, r.b.label)) for r in ranked[:5]]
    assert frozenset((PCLM1, PHL391)) in top
    assert frozenset((PCLM3, PCLM2)) in top


def test_backbone_sharing_sets_a_containment_floor(plasmid_dir):
    """Everything shares pcDNA3.1/pePB, so low containment is uninformative -- documents why
    the lineage threshold has to be high."""
    entries, _ = scan(plasmid_dir)
    assert len(relationships(entries, min_containment=0.5)) > len(
        relationships(entries, min_containment=0.98)
    )
