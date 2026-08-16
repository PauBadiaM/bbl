"""Tests for the sourcing ladders and the complexity placeholder."""

import pytest
from Bio.Seq import Seq

from bbl.complexity import complexity_score, longest_homopolymer, longest_repeat, repeat_fraction
from bbl.config import DEFAULTS, is_available, is_base_vector, load_lab_config
from bbl.inventory import find_sequence, scan_inventory
from bbl.sourcing import (
    ASK_USER,
    FROM_BASE_VECTOR,
    FROM_CONSTRUCT,
    FROM_INVENTORY,
    ORDER_BACKBONE,
    ORDER_PLASMID,
    SYNTHESIZE,
    source_backbone,
    source_insert,
)

BASE_VECTOR = "pHL162_pcDNA3.1_MCS"
PCLM1 = "pCLM1_pcDNA3.1_miniCMV-mCherry-LambdaBoxBx8"


@pytest.fixture(scope="module")
def entries(plasmid_dir):
    found, failures = scan_inventory(plasmid_dir)
    assert not failures
    return found


@pytest.fixture(scope="module")
def boxb(entries):
    """The 8x BoxB array -- repetitive enough that a vendor would refuse it."""
    pclm1 = next(e for e in entries if e.label == PCLM1)
    return pclm1.sequence[961:1255]


def plain_sequence(length, seed=7):
    """Deterministic non-repetitive filler that passes the screen."""
    import random

    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


# --------------------------------------------------------------------------- #
# complexity placeholder
# --------------------------------------------------------------------------- #


def test_metrics():
    assert longest_homopolymer("ACGTTTTACG") == 4
    assert longest_repeat("ABCDEFABCDEF") == 6
    assert longest_repeat("ACGTACGA") >= 3
    assert repeat_fraction("ACGT" * 20, window=8) == 1.0


def test_ordinary_fragment_passes():
    report = complexity_score(plain_sequence(600))
    assert report.synthesizable
    assert not report.reasons
    assert report.is_placeholder  # loudly flagged until the real function lands


def test_boxb_array_is_refused(boxb):
    """The real motivating case: tandem arrays are what vendors actually reject."""
    report = complexity_score(boxb)
    assert not report.synthesizable
    assert any("repeat" in reason or "repetitive" in reason for reason in report.reasons)
    assert report.metrics["longest_repeat"] > 40


def test_too_short_is_refused():
    report = complexity_score("ATGGCCGGCTAA")
    assert not report.synthesizable
    assert any("minimum" in reason for reason in report.reasons)


def test_limits_are_overridable():
    sequence = plain_sequence(600)
    assert complexity_score(sequence).synthesizable
    assert not complexity_score(sequence, limits={"max_length": 100}).synthesizable


# --------------------------------------------------------------------------- #
# provenance
# --------------------------------------------------------------------------- #


def test_boxb_is_found_across_the_library(entries, boxb):
    """Sequence search finds copies that a label search would miss."""
    labels = {hit.label.split("_")[0] for hit in find_sequence(entries, boxb)}
    assert {"pCLM1", "pCLM2", "pHL391"} <= labels
    assert len(labels) >= 5  # also pCLM22, pCLM24, pHL229


def test_reverse_strand_matches_are_found(entries, boxb):
    reverse = str(Seq(boxb).reverse_complement())
    hits = find_sequence(entries, reverse)
    assert hits and all(hit.strand == -1 for hit in hits)


def test_absent_sequence_is_not_found(entries):
    assert find_sequence(entries, plain_sequence(300, seed=99)) == []


# --------------------------------------------------------------------------- #
# insert ladder
# --------------------------------------------------------------------------- #


def test_unknown_sequence_asks_the_user(entries):
    decision = source_insert(None, entries, name="barcode", description="a new barcode")
    assert decision.route == ASK_USER
    assert decision.needs_user_input
    assert "Addgene" in decision.question
    assert "barcode" in decision.question


def test_sequence_in_the_library_is_cloned_not_synthesised(entries, boxb):
    """The repetitive array never reaches the complexity gate -- it is moved, not made."""
    decision = source_insert(boxb, entries, name="Lambda BoxB x8")
    assert decision.route == FROM_INVENTORY
    assert decision.complexity is None
    assert {d.label.split("_")[0] for d in decision.donors} >= {"pCLM1", "pHL391"}


def test_novel_simple_sequence_is_synthesised(entries):
    decision = source_insert(plain_sequence(400, seed=11), entries, name="new part")
    assert decision.route == SYNTHESIZE
    assert decision.complexity.synthesizable


def test_novel_repetitive_sequence_forces_ordering_the_plasmid(entries, boxb):
    """A repetitive sequence that is *not* in the library has no workaround."""
    novel = boxb[:150] + "TTGACCAGTCA" + boxb[150:]  # perturbed so it is not found
    assert find_sequence(entries, novel) == []
    decision = source_insert(novel, entries, name="novel array")
    assert decision.route == ORDER_PLASMID
    assert decision.blocks_cloning
    assert any("order the finished plasmid" in reason for reason in decision.rationale)


def test_split_hint_is_opt_in(entries, boxb):
    novel = boxb[:150] + "TTGACCAGTCA" + boxb[150:]
    without = source_insert(novel, entries, name="x")
    with_split = source_insert(novel, entries, name="x", allow_split=True)
    assert not any("3-fragment" in r for r in without.rationale)
    assert any("3-fragment" in r for r in with_split.rationale)


# --------------------------------------------------------------------------- #
# backbone ladder
# --------------------------------------------------------------------------- #


def test_prefers_an_existing_construct(entries):
    decision = source_backbone([PCLM1], entries)
    assert decision.route == FROM_CONSTRUCT
    assert decision.label == PCLM1
    assert not decision.must_order


def test_falls_back_to_the_base_vector(entries):
    decision = source_backbone([], entries)
    assert decision.route == FROM_BASE_VECTOR
    assert decision.label == BASE_VECTOR


def test_skips_unavailable_candidates(entries):
    config = {**DEFAULTS, "unavailable": [PCLM1]}
    decision = source_backbone([PCLM1], entries, config=config)
    assert decision.route == FROM_BASE_VECTOR


def test_orders_a_backbone_when_nothing_fits(entries):
    config = {**DEFAULTS, "base_vectors": []}
    decision = source_backbone([], entries, config=config, spec_backbone="pcDNA3.1")
    assert decision.route == ORDER_BACKBONE
    assert decision.must_order
    assert "pcDNA3.1" in decision.question
    assert any("sequence is known" in reason for reason in decision.rationale)


# --------------------------------------------------------------------------- #
# config
# --------------------------------------------------------------------------- #


def test_lab_config_tags_the_base_vector():
    config = load_lab_config()
    assert BASE_VECTOR in config["base_vectors"]
    assert is_base_vector(BASE_VECTOR, config)
    assert not is_base_vector(PCLM1, config)


def test_everything_available_by_default():
    assert is_available(PCLM1, load_lab_config())
