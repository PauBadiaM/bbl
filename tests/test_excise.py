"""Tests for restriction-based feature excision.

The headline test is a regression against a real cloning step: the experimentalist made
pCLM1 from pHL391 by cutting with MfeI + EcoRI and religating, to delete the NF-kB
response element and leave a promoter-only control.
"""

import pytest
from Bio.Restriction import Acc65I, EcoRI, KpnI, MfeI

from bbl import (
    COMPATIBLE_OVERHANG,
    NoExcisionFound,
    circular_equal,
    excise_features,
    ligation_strategy,
    load_plasmid,
    plan_excisions,
)
from bbl.enzymes import build_pool, enumerate_cut_sites
from bbl.plasmid_io import feature_label

# --------------------------------------------------------------------------- #
# ground truth
# --------------------------------------------------------------------------- #


def test_reproduces_pclm1_from_phl391(phl391, pclm1):
    """The function must independently pick MfeI+EcoRI and build pCLM1."""
    plan = excise_features(phl391, ["NFKBRE"])

    assert set(plan.enzyme_pair) == {"MfeI", "EcoRI"}
    assert plan.strategy == COMPATIBLE_OVERHANG

    expected = load_plasmid(pclm1)
    assert len(plan.product) == len(expected) == 5704
    assert circular_equal(plan.product, expected)

    # 5770 -> 5704; NFKBRE is 54 bp, so 12 bp of unannotated spacer goes too.
    assert plan.deleted_bp == 66
    assert plan.collateral_bp == 12

    # MfeI (C^AATTG) x EcoRI (G^AATTC) fuse to CAATTC: neither site survives.
    assert "CAATTC" in plan.junction_seq
    assert plan.sites_regenerated is False


def test_product_feature_set_matches_pclm1(phl391, pclm1):
    plan = excise_features(phl391, ["NFKBRE"])
    expected = load_plasmid(pclm1)

    def functional(record):
        return {
            feature_label(f)
            for f in record.features
            if f.type not in {"primer_bind", "source"}
        }

    produced, target = functional(plan.product), functional(expected)
    assert "NFKBRE" in produced ^ target or "NFKBRE" not in produced
    # every functional element of the real product is present in ours
    assert target - produced == set()
    # and no NF-kB element survived
    assert not {name for name in produced if name.upper().startswith("NFKB")}


def test_removed_features_include_nested_elements(phl391):
    """Asking for NFKBRE must also drop the sub-features nested inside it."""
    plan = excise_features(phl391, ["NFKBRE"])
    removed = {name.upper() for name in plan.removed_features}
    assert "NFKBRE" in removed
    assert {"NFKB_CANONICAL1", "NFKB_CANONICAL2", "NFKB_VARIANT"} <= removed


# --------------------------------------------------------------------------- #
# why "closest site" is the wrong rule
# --------------------------------------------------------------------------- #


def test_nearest_upstream_site_is_a_trap(phl391):
    """KpnI is closer to the target than MfeI but cannot ligate to EcoRI.

    This is the case that makes ligatable-first ranking mandatory rather than cosmetic.
    """
    record = load_plasmid(phl391)
    sites = enumerate_cut_sites(record.seq, build_pool(["MfeI", "KpnI", "EcoRI"]))
    mfei, kpni, ecori = sites["MfeI"][0], sites["KpnI"][0], sites["EcoRI"][0]

    target_start = 172  # NFKBRE
    assert kpni.top > mfei.top  # KpnI is nearer the feature
    assert kpni.top <= target_start

    assert ligation_strategy(kpni.enzyme, ecori.enzyme) is None
    assert ligation_strategy(mfei.enzyme, ecori.enzyme) == COMPATIBLE_OVERHANG


def test_same_overhang_sequence_is_not_sufficient():
    """KpnI and Acc65I both leave GTAC, but 3' vs 5' -- they must not be paired."""
    assert KpnI.ovhgseq == Acc65I.ovhgseq == "GTAC"
    assert KpnI.is_3overhang() and Acc65I.is_5overhang()
    assert ligation_strategy(KpnI, Acc65I) is None


def test_mfei_ecori_are_compatible():
    assert MfeI.ovhgseq == EcoRI.ovhgseq == "AATT"
    assert MfeI.ovhg == EcoRI.ovhg == -4
    assert ligation_strategy(MfeI, EcoRI) == COMPATIBLE_OVERHANG


# --------------------------------------------------------------------------- #
# guard rails
# --------------------------------------------------------------------------- #


def test_alternatives_are_all_ligatable(phl391):
    plan = excise_features(phl391, ["NFKBRE"])
    for alt in plan.alternatives:
        assert alt.strategy in {"compatible_overhang", "same_enzyme", "blunt"}
        assert alt.collateral_bp >= plan.collateral_bp or alt.strategy != plan.strategy


def test_protecting_a_target_feature_is_rejected(phl391):
    with pytest.raises(ValueError, match="contradictory"):
        excise_features(phl391, ["NFKBRE"], protect=["NFKB_Variant"])


def test_unknown_feature_name_lists_alternatives(phl391):
    with pytest.raises(ValueError, match="No features matched"):
        excise_features(phl391, ["not_a_real_feature"])


def test_no_solution_reports_a_reason_not_a_crash(phl391):
    """A feature hemmed in by protected elements has no restriction route."""
    with pytest.raises(NoExcisionFound) as excinfo:
        excise_features(phl391, ["Kozak sequence"], margin=10)
    message = str(excinfo.value)
    assert "PCR-based deletion" in message
    assert excinfo.value.reasons


def test_ampr_has_no_restriction_route_in_this_plasmid(phl391):
    """AmpR is hemmed in by ori and the AmpR promoter, so there is honestly no route.

    The essential-element *warning* is covered deterministically in test_synthetic.py; what
    matters here is that an impossible request fails loudly rather than returning something
    that would not propagate.
    """
    plans, reasons = plan_excisions(phl391, ["AmpR"])
    assert not plans
    assert reasons


def test_margin_shrinks_the_option_space(phl391):
    loose, _ = plan_excisions(phl391, ["NFKBRE"], margin=0)
    tight, _ = plan_excisions(phl391, ["NFKBRE"], margin=200)
    assert len(tight) <= len(loose)


# --------------------------------------------------------------------------- #
# inventory smoke test
# --------------------------------------------------------------------------- #


def test_every_inventory_plasmid_loads(plasmid_dir):
    files = sorted(plasmid_dir.glob("*.dna"))
    assert files, "no .dna files found"
    for path in files:
        record = load_plasmid(path)
        assert len(record) > 1000, path.name
        assert record.circular


# --------------------------------------------------------------------------- #
# readout-handle survival
# --------------------------------------------------------------------------- #


def test_readout_handles_are_checked_without_being_asked(plasmid_dir):
    """Losing a sequencing handle is invisible until sequencing comes back unreadable.

    The plasmid still propagates and still sequences clean -- it is just no longer readable at
    that end. ``protect`` catches this only if the caller remembers to name the handle, so the
    check runs unprompted instead.
    """
    from bbl.plasmid_io import delete_span
    from bbl.targets import READOUT_HANDLE_HINTS, handle_warnings

    record = load_plasmid(
        plasmid_dir
        / "pHL394_pLV-U6rev-Tornado-TruseqR2-oFH155-BCS-oFH99-10XCS1-LambdaBoxBx2-EF1-mScarlet.dna"
    )
    handles = [
        f
        for f in record.features
        if (label := feature_label(f))
        and any(hint in label.lower() for hint in READOUT_HANDLE_HINTS)
    ]
    assert len(handles) > 3  # precondition: this plasmid is annotated with handles

    # an untouched product must be silent
    assert handle_warnings(record, record) == []

    truseq = next(
        f for f in record.features if feature_label(f) == "Illumina Truseq Read 2 Primer"
    )
    damaged = delete_span(record, int(truseq.location.start), int(truseq.location.end))
    warnings = handle_warnings(record, damaged)
    assert len(warnings) == 1
    assert "Illumina Truseq Read 2 Primer" in warnings[0]


def test_the_hint_list_matches_annotated_labels_not_filenames(plasmid_dir):
    """Regression on a live bug in the first version of the registry.

    ``pHL394``'s *filename* says ``10XCS1``; its *feature* is labelled ``10X Capture
    Sequence 1``. A hint of ``"10xcs1"`` therefore matched nothing at all, which is the silent
    kind of failure -- a registry that looks populated and detects nothing.
    """
    from bbl.targets import READOUT_HANDLE_HINTS

    record = load_plasmid(
        plasmid_dir
        / "pHL394_pLV-U6rev-Tornado-TruseqR2-oFH155-BCS-oFH99-10XCS1-LambdaBoxBx2-EF1-mScarlet.dna"
    )
    labels = {feature_label(f) for f in record.features if feature_label(f)}
    matched = {
        label
        for label in labels
        if any(hint in label.lower() for hint in READOUT_HANDLE_HINTS)
    }
    assert "10X Capture Sequence 1" in matched
    assert "Illumina Truseq Read 2 Primer" in matched
    assert "oFH155_HC1F" in matched


def test_the_real_deletion_gains_no_spurious_handle_warning(phl391, pclm1):
    """pHL391 -> pCLM1 touches no readout handle, so the new check must stay quiet."""
    plan = excise_features(phl391, ["NFKBRE"])
    assert circular_equal(plan.product, load_plasmid(pclm1))
    assert not any("readout handle" in warning for warning in plan.warnings)
