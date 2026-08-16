"""Synthetic-plasmid tests for cases the pHL391 -> pCLM1 regression cannot reach.

pHL391 exercises only a 5'-overhang, two-enzyme deletion. These construct plasmids that
force the other paths: a 3'-overhang single-enzyme digest, a blunt pair, the
essential-element warning, and the no-viable-route path.

Each test passes an explicit ``enzyme_pool``, which both makes the expected answer
deterministic and exercises the lab-stock-list hook.
"""

import random

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord

from bbl import BLUNT, SAME_ENZYME, excise_features, plan_excisions
from bbl.enzymes import build_pool, enumerate_cut_sites

FORBIDDEN = ("CTGCAG", "CCCGGG", "GAGCTC")  # sites we place deliberately
BACKBONE_BP = 200
TARGET_BP = 30


def filler(length, seed):
    """Deterministic non-repetitive sequence containing none of the placed sites."""
    rng = random.Random(seed)
    while True:
        sequence = "".join(rng.choice("ACGT") for _ in range(length))
        if not any(site in sequence for site in FORBIDDEN):
            return sequence


def build_plasmid(upstream_site, downstream_site, middle_label="junk"):
    """``ori | site | target | site | AmpR``, circular."""
    backbone_a = filler(BACKBONE_BP, seed=1)
    middle = filler(TARGET_BP, seed=2)
    backbone_b = filler(BACKBONE_BP, seed=3)
    sequence = backbone_a + upstream_site + middle + downstream_site + backbone_b

    start_mid = len(backbone_a) + len(upstream_site)
    start_b = start_mid + len(middle) + len(downstream_site)

    record = SeqRecord(
        Seq(sequence),
        id="synthetic",
        name="synthetic",
        description="synthetic test plasmid",
        annotations={"molecule_type": "ds-DNA", "topology": "circular"},
    )
    for start, end, label, kind in (
        (0, len(backbone_a), "ori", "rep_origin"),
        (start_mid, start_mid + len(middle), middle_label, "misc_feature"),
        (start_b, len(sequence), "AmpR", "CDS"),
    ):
        record.features.append(
            SeqFeature(
                SimpleLocation(start, end, strand=1),
                type=kind,
                qualifiers={"label": [label]},
            )
        )
    return Dseqrecord(record, circular=True)


def labels_of(record):
    return {
        f.qualifiers["label"][0]: (int(f.location.start), int(f.location.end))
        for f in record.features
        if "label" in f.qualifiers
    }


def test_single_enzyme_three_prime_overhang():
    """Two PstI sites (C_TGCA^G) flanking the target.

    This is the arithmetic the real case never checks: for a 3' overhang the bottom-strand
    break sits *upstream* of the top-strand break, so a footprint-based length calculation
    would be wrong. What the product actually loses is the top-strand span.
    """
    plasmid = build_plasmid("CTGCAG", "CTGCAG")
    sites = enumerate_cut_sites(plasmid.seq, build_pool(["PstI"]))
    assert len(sites["PstI"]) == 2
    first, second = sorted(sites["PstI"], key=lambda s: s.top)
    assert first.bottom < first.top  # 3' overhang: bottom break precedes top break

    plan = excise_features(plasmid, ["junk"], enzyme_pool=["PstI"])
    assert plan.strategy == SAME_ENZYME
    assert plan.enzyme_pair == ("PstI", "PstI")
    assert plan.deleted_bp == second.top - first.top
    assert len(plan.product) == len(plasmid) - plan.deleted_bp
    assert plan.sites_regenerated is True  # identical ends restore the site
    assert "junk" not in labels_of(plan.product)
    assert "dephosphorylate" in plan.protocol.lower()


def test_blunt_pair():
    """SmaI (CCC^GGG) and Ecl136II (GAG^CTC) both leave blunt ends."""
    plasmid = build_plasmid("CCCGGG", "GAGCTC")
    plan = excise_features(plasmid, ["junk"], enzyme_pool=["SmaI", "Ecl136II"])
    assert plan.strategy == BLUNT
    assert len(plan.product) == len(plasmid) - plan.deleted_bp
    assert any("blunt ligation is inefficient" in w for w in plan.warnings)


def test_incompatible_pair_is_rejected():
    """SacI leaves a 3' AGCT overhang; it cannot be ligated to blunt SmaI."""
    plasmid = build_plasmid("CCCGGG", "GAGCTC")
    plans, reasons = plan_excisions(plasmid, ["junk"], enzyme_pool=["SmaI", "SacI"])
    assert not plans
    assert "incompatible ends" in reasons


def test_essential_element_removal_warns():
    plasmid = build_plasmid("CTGCAG", "CTGCAG", middle_label="f1 ori")
    plans, _ = plan_excisions(plasmid, ["f1 ori"], enzyme_pool=["PstI"])
    assert plans, "expected a viable digest"
    assert any("essential" in w for w in plans[0].warnings)


def test_protected_neighbour_blocks_the_only_route():
    """A large margin makes the flanking cuts collide with ori/AmpR, so nothing is viable."""
    plasmid = build_plasmid("CTGCAG", "CTGCAG")
    plans, reasons = plan_excisions(plasmid, ["junk"], enzyme_pool=["PstI"], margin=50)
    assert not plans
    assert any("would damage" in reason for reason in reasons)


def test_coordinates_are_reported_in_the_parent_frame():
    plasmid = build_plasmid("CTGCAG", "CTGCAG")
    plan = excise_features(plasmid, ["junk"], enzyme_pool=["PstI"])
    labels = labels_of(plan.product)
    start_b = BACKBONE_BP + 6 + TARGET_BP + 6
    # ori is upstream of the cut, so unchanged...
    assert labels["ori"] == (0, BACKBONE_BP)
    # ...and AmpR shifts down by exactly the deleted length.
    assert labels["AmpR"] == (start_b - plan.deleted_bp, len(plasmid) - plan.deleted_bp)


def test_deleting_a_tandem_repeat_is_not_confused_by_duplicates():
    """A target whose sequence recurs elsewhere must still verify as deleted.

    Guards the junction-based check in ``_verify``: searching for the target sequence
    itself would wrongly report failure here.
    """
    repeat = filler(TARGET_BP, seed=2)  # same sequence the target uses
    plasmid = build_plasmid("CTGCAG", "CTGCAG")
    extended = Dseqrecord(
        str(plasmid.seq) + repeat, circular=True
    )  # a second, unannotated copy
    for feature in plasmid.features:
        extended.features.append(feature)
    assert str(extended.seq).count(repeat) == 2

    plan = excise_features(extended, ["junk"], enzyme_pool=["PstI"])
    assert len(plan.product) == len(extended) - plan.deleted_bp
    assert str(plan.product.seq).count(repeat) == 1  # only the unannotated copy survives
