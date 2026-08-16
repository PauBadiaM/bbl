"""Tests for around-the-horn (inverse) PCR deletion.

The case: ``test_pHL391...dna`` is pHL391 with the MfeI and EcoRI sites removed, so the
restriction route that built pCLM1 is unavailable and the target must be deleted by PCR.
"""

import pytest
from Bio.Restriction import EcoRI, MfeI
from Bio.Seq import Seq

from bbl import (
    GIBSON,
    KLD,
    NoExcisionFound,
    NoPrimerDesign,
    delete_span,
    design_deletion_primers,
    excise_features,
    load_plasmid,
)
from bbl.pcr import (
    _binding_sites,
    _candidates,
    _mispriming_sites,
    _pair_heterodimer_stem,
    _penalty,
)
from bbl.plasmid_io import feature_label

NFKBRE_START, NFKBRE_END = 160, 214
NFKBRE_BP = NFKBRE_END - NFKBRE_START


def test_the_test_plasmid_really_has_no_restriction_route(no_sites_phl391):
    """Precondition: this file exists precisely because restriction cannot be used."""
    record = load_plasmid(no_sites_phl391)
    assert MfeI.search(record.seq, linear=False) == []
    assert EcoRI.search(record.seq, linear=False) == []
    with pytest.raises(NoExcisionFound):
        excise_features(record, ["NFKBRE"])


def test_primers_flank_the_target_exactly(no_sites_phl391):
    record = load_plasmid(no_sites_phl391)
    template = str(record.seq).upper()
    design = design_deletion_primers(no_sites_phl391, ["NFKBRE"])

    # The 5' end of each primer is pinned to a deletion boundary; they point outward.
    assert design.forward.binding_span[0] == NFKBRE_END
    assert design.reverse.binding_span[1] == NFKBRE_START
    assert design.forward.strand == 1 and design.reverse.strand == -1

    forward_len = len(design.forward.anneal)
    reverse_len = len(design.reverse.anneal)
    assert design.forward.anneal == template[NFKBRE_END : NFKBRE_END + forward_len]
    assert design.reverse.anneal == str(
        Seq(template[NFKBRE_START - reverse_len : NFKBRE_START]).reverse_complement()
    )


def test_deletion_is_exact_with_no_collateral(no_sites_phl391):
    record = load_plasmid(no_sites_phl391)
    design = design_deletion_primers(no_sites_phl391, ["NFKBRE"])

    assert design.deleted_bp == NFKBRE_BP == 54
    assert design.deleted_span == (NFKBRE_START, NFKBRE_END)
    assert len(design.product) == len(record) - NFKBRE_BP
    # the product is precisely parent-minus-target, base for base
    expected = str(record.seq[:NFKBRE_START]) + str(record.seq[NFKBRE_END:])
    assert str(design.product.seq).upper() == expected.upper()


def test_nfkb_features_are_gone_and_the_rest_survive(no_sites_phl391):
    record = load_plasmid(no_sites_phl391)
    design = design_deletion_primers(no_sites_phl391, ["NFKBRE"])
    produced = {feature_label(f) for f in design.product.features}
    assert not {name for name in produced if name.upper().startswith("NFKB")}
    for keep in ("miniCMV", "mCherry", "AmpR", "ori", "bGH poly(A) signal"):
        assert keep in produced
    assert {feature_label(f) for f in record.features if feature_label(f).startswith("NFKB")}


def test_amplicon_and_protocol(no_sites_phl391):
    record = load_plasmid(no_sites_phl391)
    design = design_deletion_primers(no_sites_phl391, ["NFKBRE"])
    assert design.amplicon_bp == len(record) - NFKBRE_BP
    assert design.method == KLD
    assert design.forward.phosphorylated and design.reverse.phosphorylated
    assert "DpnI" in design.protocol
    assert "phosphorylate" in design.protocol.lower()


def test_primer_quality(no_sites_phl391):
    design = design_deletion_primers(no_sites_phl391, ["NFKBRE"], target_tm=60.0)
    for primer in (design.forward, design.reverse):
        assert 18 <= primer.length <= 32
        assert 35.0 <= primer.gc <= 65.0
        assert abs(primer.tm - 60.0) <= 3.0
    assert design.tm_difference <= 3.0
    # neither primer should have a mispriming note on this template
    assert not [note for p in (design.forward, design.reverse) for note in p.notes]


# --------------------------------------------------------------------------- #
# Gibson variant
# --------------------------------------------------------------------------- #


def test_gibson_overlap_is_split_evenly(no_sites_phl391):
    """20 bp of homology, 10 on each primer, straddling the new junction."""
    record = load_plasmid(no_sites_phl391)
    template = str(record.seq).upper()
    design = design_deletion_primers(no_sites_phl391, ["NFKBRE"], method=GIBSON, overlap=20)

    assert len(design.forward.tail) == len(design.reverse.tail) == 10
    assert design.homology_bp == 20
    assert design.forward.tail == template[NFKBRE_START - 10 : NFKBRE_START]
    assert design.reverse.tail == str(
        Seq(template[NFKBRE_END : NFKBRE_END + 10]).reverse_complement()
    )
    assert design.forward.sequence == design.forward.tail + design.forward.anneal
    assert design.reverse.sequence == design.reverse.tail + design.reverse.anneal
    assert not design.forward.phosphorylated  # Gibson needs no kinase
    assert design.amplicon_bp == len(record) - NFKBRE_BP + 20
    assert "NEBuilder" in design.protocol


def test_gibson_amplicon_has_a_terminal_direct_repeat(no_sites_phl391):
    """The functional requirement: both amplicon ends must share the same 20 bp."""
    template = str(load_plasmid(no_sites_phl391).seq).upper()
    design = design_deletion_primers(no_sites_phl391, ["NFKBRE"], method=GIBSON, overlap=20)

    amplicon = (
        design.forward.tail
        + template[NFKBRE_END:]
        + template[:NFKBRE_START]
        + str(Seq(design.reverse.tail).reverse_complement())
    )
    assert len(amplicon) == design.amplicon_bp
    assert amplicon[:20] == amplicon[-20:]
    # and that shared region is exactly the sequence spanning the new junction
    assert amplicon[:20] == template[NFKBRE_START - 10 : NFKBRE_START] + template[
        NFKBRE_END : NFKBRE_END + 10
    ]


def test_gibson_overlap_can_go_all_on_the_forward_primer(no_sites_phl391):
    template = str(load_plasmid(no_sites_phl391).seq).upper()
    design = design_deletion_primers(
        no_sites_phl391, ["NFKBRE"], method=GIBSON, overlap=20, overlap_placement="forward"
    )
    assert design.forward.tail == template[NFKBRE_START - 20 : NFKBRE_START]
    assert design.reverse.tail == ""
    assert design.homology_bp == 20


def test_odd_overlap_splits_with_the_extra_base_forward(no_sites_phl391):
    design = design_deletion_primers(no_sites_phl391, ["NFKBRE"], method=GIBSON, overlap=15)
    assert (len(design.forward.tail), len(design.reverse.tail)) == (8, 7)
    assert design.homology_bp == 15


def test_bad_overlap_placement_rejected(no_sites_phl391):
    with pytest.raises(ValueError, match="overlap_placement"):
        design_deletion_primers(
            no_sites_phl391, ["NFKBRE"], method=GIBSON, overlap_placement="middle"
        )


def test_kld_has_no_tails(no_sites_phl391):
    design = design_deletion_primers(no_sites_phl391, ["NFKBRE"], method=KLD)
    assert design.homology_bp == 0
    assert design.forward.tail == design.reverse.tail == ''


def test_gibson_and_kld_give_the_same_product(no_sites_phl391):
    kld = design_deletion_primers(no_sites_phl391, ["NFKBRE"], method=KLD)
    gibson = design_deletion_primers(no_sites_phl391, ["NFKBRE"], method=GIBSON)
    assert str(kld.product.seq) == str(gibson.product.seq)


# --------------------------------------------------------------------------- #
# relationship to the restriction route
# --------------------------------------------------------------------------- #


def test_warns_when_restriction_would_work(phl391):
    """On the original pHL391, MfeI+EcoRI is available and should be preferred."""
    design = design_deletion_primers(phl391, ["NFKBRE"])
    assert any("restriction route exists" in w for w in design.warnings)
    assert any("MfeI" in w and "EcoRI" in w for w in design.warnings)


def test_pcr_is_exact_where_restriction_is_not(phl391):
    """Same parent, same request: PCR removes 54 bp, restriction removes 66 bp."""
    pcr = design_deletion_primers(phl391, ["NFKBRE"], warn_if_restriction_possible=False)
    cut = excise_features(phl391, ["NFKBRE"])
    assert pcr.deleted_bp == 54
    assert cut.deleted_bp == 66 and cut.collateral_bp == 12
    assert len(pcr.product) == len(cut.product) + 12


def test_mispriming_in_a_repeat_is_flagged(phl391):
    """Deleting one BoxB repeat puts the primers inside the 8x array -- must warn."""
    record = load_plasmid(phl391)
    boxb = [
        (int(f.location.start), int(f.location.end))
        for f in record.features
        if feature_label(f) == "BoxB RNA aptamer"
    ]
    assert len(boxb) == 8
    design = design_deletion_primers(
        record, [boxb[3]], warn_if_restriction_possible=False
    )
    assert any("mispriming" in w for w in design.warnings)


# --------------------------------------------------------------------------- #
# helpers and error paths
# --------------------------------------------------------------------------- #


def test_specificity_changes_which_primer_is_chosen(phl391):
    """The point of the whole exercise: mispriming now *ranks*, it does not just annotate.

    ``_binding_sites`` used to run only in ``_annotate``, after ``_select_pair`` had already
    committed -- so a primer whose 3' probe occurs 8x in the plasmid could be returned as the
    best available with a note attached. Inside the 8x BoxB array at a 56 C target, that is
    exactly what the old objective did.
    """
    record = load_plasmid(phl391)
    template = str(record.seq).upper()
    boxb = [
        (int(f.location.start), int(f.location.end))
        for f in record.features
        if feature_label(f) == "BoxB RNA aptamer"
    ]
    options = _candidates(template, boxb[3][1], +1, 18, 32, 56.0, 35.0, 65.0)

    def specificity_blind(option):
        return _penalty(option["sequence"], option["tm"], 56.0, 35.0, 65.0, template=None)

    chosen = min(options, key=lambda option: option["penalty"])
    would_have_chosen = min(options, key=specificity_blind)

    blind_sites = _mispriming_sites(template, would_have_chosen["sequence"])
    chosen_sites = _mispriming_sites(template, chosen["sequence"])
    assert blind_sites == 8
    assert chosen_sites == 2
    assert chosen["sequence"] != would_have_chosen["sequence"]


def test_an_unavoidable_mispriming_is_still_reported(phl391):
    """Ranking must not silence a warning it cannot act on.

    Every candidate at a BoxB boundary sits inside the array, so no length escapes mispriming.
    The penalty picks the least bad one and the warning still fires.
    """
    record = load_plasmid(phl391)
    template = str(record.seq).upper()
    boxb = [
        (int(f.location.start), int(f.location.end))
        for f in record.features
        if feature_label(f) == "BoxB RNA aptamer"
    ]
    options = _candidates(template, boxb[3][1], +1, 18, 32, 60.0, 35.0, 65.0)
    assert all(_mispriming_sites(template, o["sequence"]) > 1 for o in options)

    design = design_deletion_primers(record, [boxb[3]], warn_if_restriction_possible=False)
    assert any("mispriming" in warning for warning in design.warnings)


def test_heterodimers_are_scored_on_the_annealing_region_only(no_sites_phl391):
    """A Gibson design's 5' tails are template-complementary by construction.

    Scoring the ordered oligo would flag every correct assembly primer as a dimer, so the pair
    term must see ``anneal``, not ``tail + anneal``. The KLD and Gibson designs share their
    annealing regions, so their heterodimer measurement must be identical.
    """
    kld = design_deletion_primers(no_sites_phl391, ["NFKBRE"], method=KLD)
    gibson = design_deletion_primers(no_sites_phl391, ["NFKBRE"], method=GIBSON, overlap=20)

    assert gibson.forward.tail and gibson.reverse.tail  # precondition: tails exist
    assert kld.forward.anneal == gibson.forward.anneal

    on_anneal = _pair_heterodimer_stem(gibson.forward.anneal, gibson.reverse.anneal)
    on_ordered = _pair_heterodimer_stem(gibson.forward.sequence, gibson.reverse.sequence)
    assert on_anneal == _pair_heterodimer_stem(kld.forward.anneal, kld.reverse.anneal)
    # the tails inflate the apparent dimer, which is why they are excluded
    assert on_ordered >= on_anneal


def test_binding_sites_does_not_double_count(no_sites_phl391):
    template = str(load_plasmid(no_sites_phl391).seq).upper()
    assert _binding_sites(template, template[1000:1020]) == 1
    # a probe spanning the origin is still found exactly once
    assert _binding_sites(template, template[-10:] + template[:10]) == 1


def test_bad_method_rejected(no_sites_phl391):
    with pytest.raises(ValueError, match="method must be"):
        design_deletion_primers(no_sites_phl391, ["NFKBRE"], method="golden-gate")


def test_unknown_feature_rejected(no_sites_phl391):
    with pytest.raises(ValueError, match="No features matched"):
        design_deletion_primers(no_sites_phl391, ["not_a_feature"])


def test_target_too_large_for_priming(no_sites_phl391):
    record = load_plasmid(no_sites_phl391)
    with pytest.raises(NoPrimerDesign, match="too little backbone"):
        design_deletion_primers(record, [(0, len(record) - 10)])


def test_deletion_starting_at_position_zero(no_sites_phl391):
    """Guards the circular-Dseq slicing trap: ``seq[:0]`` returns the whole circle.

    Before the fix this silently produced a double-length plasmid.
    """
    record = load_plasmid(no_sites_phl391)
    assert len(delete_span(record, 0, 40)) == len(record) - 40

    design = design_deletion_primers(record, [(0, 40)], warn_if_restriction_possible=False)
    assert len(design.product) == len(record) - 40
    assert str(design.product.seq).upper() == str(record.seq)[40:].upper()


def test_target_spanning_the_origin(no_sites_phl391):
    """A target that wraps the origin is rotated, and the primers still abut it."""
    record = load_plasmid(no_sites_phl391)
    length = len(record)
    design = design_deletion_primers(
        record, [(length - 30, length), (0, 20)], warn_if_restriction_possible=False
    )
    assert design.deleted_bp == 50
    assert len(design.product) == length - 50
