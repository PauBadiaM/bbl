"""Tests for insertion.

Ground truth: **pCLM2 = pCLM3 with the Lambda BoxB x8 array put in.** pCLM3 carries a 21 bp
stub at [1491:1512] where pCLM2 has the 294 bp array, and the two plasmids are byte-identical
either side. The array is also identical in pCLM1 and pHL391, so pCLM1 serves as the donor.

That case happens to have a working restriction route (XhoI + XbaI), which is the "rare" branch
-- it is not rare here because both plasmids descend from the same pcDNA3.1 MCS.
"""

import pytest

from bbl import (
    GIBSON_PCR,
    GIBSON_SYNTHETIC,
    RESTRICTION,
    NoInsertionRoute,
    load_plasmid,
    plan_insertion,
    resolve_insert,
)
from bbl.plasmid_io import feature_label

SITE = (1491, 1512)  # the 21 bp stub in pCLM3 that pCLM2 replaces
BOXB_BP = 294


@pytest.fixture(scope="module")
def vector(plasmid_dir):
    return plasmid_dir / "pCLM3_pcDNA3.1_CMV-mCherry.dna"


@pytest.fixture(scope="module")
def donor(plasmid_dir):
    return plasmid_dir / "pCLM1_pcDNA3.1_miniCMV-mCherry-LambdaBoxBx8.dna"


@pytest.fixture(scope="module")
def expected(plasmid_dir):
    return plasmid_dir / "pCLM2_pcDNA3.1_CMV-mCherry-LambdaBoxBx8.dna"


@pytest.fixture(scope="module")
def boxb(donor):
    return str(load_plasmid(donor).seq)[961:1255].upper()


def test_ground_truth_relationship(vector, donor, expected, boxb):
    """Precondition: pCLM2 really is pCLM3 with the array swapped into the stub."""
    v, e = str(load_plasmid(vector).seq).upper(), str(load_plasmid(expected).seq).upper()
    assert len(boxb) == BOXB_BP
    assert v[: SITE[0]] == e[: SITE[0]]
    assert v[SITE[1] :] == e[SITE[0] + BOXB_BP :]


# --------------------------------------------------------------------------- #
# case 1: insert from another plasmid, reusing restriction enzymes
# --------------------------------------------------------------------------- #


def test_restriction_reuse_reproduces_pclm2(vector, donor, expected):
    plan = plan_insertion(vector, donor, insert_features=["Lambda BoxB x8"], at=SITE)

    assert plan.strategy == RESTRICTION
    assert set(plan.vector.enzymes) == {"XhoI", "XbaI"}
    assert set(plan.insert.enzymes) == {"XhoI", "XbaI"}

    target = load_plasmid(expected)
    assert len(plan.product) == len(target) == 6234
    assert str(plan.product.seq).upper() == str(target.seq).upper()


def test_restriction_route_is_directional(vector, donor):
    """XhoI and XbaI leave incompatible ends, so the insert cannot flip."""
    plan = plan_insertion(vector, donor, insert_features=["Lambda BoxB x8"], at=SITE)
    assert plan.directional is True
    assert "only go in one way" in plan.protocol
    # same enzyme on both sides of each junction, so the sites come back
    assert plan.sites_regenerated is True


def test_restriction_route_reports_donor_flank(vector, donor):
    """Cutting at the donor's sites brings 6 bp of donor sequence along with the array."""
    plan = plan_insertion(vector, donor, insert_features=["Lambda BoxB x8"], at=SITE)
    assert plan.inserted_bp == 300 == BOXB_BP + 6
    assert any("6 bp of donor sequence" in w for w in plan.warnings)


def test_product_carries_the_insert_annotations(vector, donor):
    plan = plan_insertion(vector, donor, insert_features=["Lambda BoxB x8"], at=SITE)
    labels = [feature_label(f) for f in plan.product.features]
    assert "Lambda BoxB x8" in labels
    assert labels.count("BoxB RNA aptamer") == 8
    for keep in ("CMV promoter", "mCherry", "AmpR", "ori"):
        assert keep in labels


def test_prefers_restriction_over_gibson(vector, donor):
    """method='auto' must pick the enzyme route when one exists."""
    assert plan_insertion(
        vector, donor, insert_features=["Lambda BoxB x8"], at=SITE
    ).strategy == RESTRICTION


# --------------------------------------------------------------------------- #
# case 1 fallback: PCR the insert and Gibson it in
# --------------------------------------------------------------------------- #


def test_gibson_pcr_inserts_exactly_the_feature(vector, donor):
    plan = plan_insertion(
        vector, donor, insert_features=["Lambda BoxB x8"], at=SITE, method="gibson"
    )
    assert plan.strategy == GIBSON_PCR
    assert plan.inserted_bp == BOXB_BP  # no donor flank tags along
    assert len(plan.product) == len(load_plasmid(vector)) - 27 + BOXB_BP


def test_gibson_primers_carry_vector_homology(vector, donor):
    plan = plan_insertion(
        vector, donor, insert_features=["Lambda BoxB x8"], at=SITE,
        method="gibson", homology=25,
    )
    forward, reverse = plan.insert.primers
    vector_text = str(load_plasmid(vector).seq).upper()
    open_start, open_end = plan.site

    assert forward.tail == vector_text[open_start - 25 : open_start]
    assert forward.sequence == forward.tail + forward.anneal
    assert len(reverse.tail) == 25
    # the annealing halves come from the donor, at the array's edges
    donor_text = str(load_plasmid(donor).seq).upper()
    assert donor_text[961:].startswith(forward.anneal)


def test_gibson_and_synthetic_agree(vector, donor, boxb):
    """Ordering the fragment gives the same plasmid as amplifying it."""
    by_pcr = plan_insertion(
        vector, donor, insert_features=["Lambda BoxB x8"], at=SITE, method="gibson"
    )
    by_order = plan_insertion(vector, boxb, at=SITE)
    assert str(by_pcr.product.seq) == str(by_order.product.seq)


# --------------------------------------------------------------------------- #
# case 2: you only have the sequence
# --------------------------------------------------------------------------- #


def test_synthetic_insert_gets_homology_arms(vector, boxb):
    plan = plan_insertion(vector, boxb, at=SITE, homology=25)
    assert plan.strategy == GIBSON_SYNTHETIC
    assert plan.insert.primers is None  # nothing to amplify
    order = plan.insert.order_sequence
    assert len(order) == BOXB_BP + 50
    assert order[25:-25] == boxb

    # the arms match the *opened vector* ends, which here come from an XhoI+XbaI digest
    # rather than the requested site
    vector_text = str(load_plasmid(vector).seq).upper()
    open_start, open_end = plan.site
    assert order.startswith(vector_text[open_start - 25 : open_start])
    assert order.endswith(vector_text[open_end : open_end + 25])


def test_short_synthetic_insert_is_flagged(vector):
    plan = plan_insertion(vector, "ATGGCCGGCTAA", at=1491, homology=25)
    assert plan.inserted_bp == 12
    assert any("cheaper" in w or "inefficient" in w for w in plan.warnings)


def test_pure_insertion_replaces_nothing(vector):
    """An integer site inserts without deleting."""
    parent = load_plasmid(vector)
    plan = plan_insertion(vector, "ATGGCCGGCTAA", at=1491, vector_prep="pcr")
    assert plan.replaced_bp == 0
    assert len(plan.product) == len(parent) + 12
    assert str(plan.product.seq).upper() == (
        str(parent.seq)[:1491] + "ATGGCCGGCTAA" + str(parent.seq)[1491:]
    ).upper()


# --------------------------------------------------------------------------- #
# insert resolution and error paths
# --------------------------------------------------------------------------- #


def test_resolve_insert_from_string(boxb):
    source = resolve_insert(boxb)
    assert source.synthetic and source.donor is None
    assert source.sequence == boxb


def test_resolve_insert_from_donor(donor, boxb):
    source = resolve_insert(donor, ["Lambda BoxB x8"])
    assert not source.synthetic
    assert source.span == (961, 1255)
    assert source.sequence == boxb
    assert len(source.features) == 9  # the array plus its 8 aptamers


def test_rejects_nonsense_insert(vector):
    with pytest.raises(ValueError, match="neither a DNA string nor a file path"):
        plan_insertion(vector, "not a sequence!", at=1491)


def test_bad_method_rejected(vector, boxb):
    with pytest.raises(ValueError, match="method must be"):
        plan_insertion(vector, boxb, at=SITE, method="goldengate")


def test_restriction_route_unavailable_for_synthetic_insert(vector, boxb):
    """You cannot digest a sequence you have not built yet."""
    with pytest.raises(NoInsertionRoute, match="cannot be reused"):
        plan_insertion(vector, boxb, at=SITE, method=RESTRICTION)


def test_vector_prep_digest_can_fail_cleanly(vector, boxb):
    """No digest may remove extra vector sequence, so no pair qualifies."""
    with pytest.raises(NoInsertionRoute, match="vector_prep='pcr'"):
        plan_insertion(
            vector, boxb, at=SITE, method="gibson",
            vector_prep="digest", max_digest_collateral=0,
        )


def test_auto_vector_prep_falls_back_to_pcr(vector, boxb):
    """With no digest allowed, the vector is opened by inverse PCR -- exactly at the site."""
    plan = plan_insertion(
        vector, boxb, at=SITE, method="gibson", max_digest_collateral=0
    )
    assert plan.site == SITE
    assert plan.vector.primers is not None
    assert plan.replaced_bp == SITE[1] - SITE[0]
    assert len(plan.product) == len(load_plasmid(vector)) - 21 + BOXB_BP


# --------------------------------------------------------------------------- #
# the pCLM2 build, via the route we would actually propose
#
# Historical route: pHL391 cut MfeI + NheI (removing NFKBRE1 + miniCMV, 132 bp), CMV amplified
# from pHL162 with primers adding matching restriction sites, ligate -> 6234 bp.
#
# We deliberately do NOT reproduce that byte-for-byte. Gibson gives the same construct with
# seamless junctions and 12 fewer bp -- the 12 being vestigial MfeI/NheI sites plus a spare
# KpnI the experimentalist parked for future cloning. See docs/DECISIONS.md D61.
# --------------------------------------------------------------------------- #

PHL391 = "pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8.dna"
PHL162 = "pHL162_pcDNA3.1_MCS.dna"
CMV_SPAN = (234, 818)  # CMV enhancer + CMV promoter in pHL162, annotated boundaries
PCLM2_SITE = (161, 293)
CMV_BP = CMV_SPAN[1] - CMV_SPAN[0]


@pytest.fixture(scope="module")
def phl391(plasmid_dir):
    return plasmid_dir / PHL391


@pytest.fixture(scope="module")
def phl162(plasmid_dir):
    return plasmid_dir / PHL162


def test_finds_the_backbone_enzymes_used_at_the_bench(phl391, phl162):
    """MfeI + NheI, chosen unprompted -- the pair the experimentalist actually used."""
    plan = plan_insertion(
        phl391, phl162, insert_features=[CMV_SPAN], at=PCLM2_SITE, method="gibson"
    )
    assert set(plan.vector.enzymes) == {"MfeI", "NheI"}
    assert plan.replaced_bp == 132  # NFKBRE1 + miniCMV


def test_gibson_builds_a_functionally_equivalent_pclm2(phl391, phl162, expected):
    parent = load_plasmid(phl391)
    donor = load_plasmid(phl162)
    plan = plan_insertion(
        phl391, phl162, insert_features=[CMV_SPAN], at=PCLM2_SITE, method="gibson"
    )
    assert plan.strategy == GIBSON_PCR
    assert plan.inserted_bp == CMV_BP == 584

    # exactly parent-minus-132 plus the CMV block, seamless
    parent_text, donor_text = str(parent.seq).upper(), str(donor.seq).upper()
    assert str(plan.product.seq).upper() == (
        parent_text[:161] + donor_text[CMV_SPAN[0] : CMV_SPAN[1]] + parent_text[293:]
    )
    # 12 bp shorter than the real construct: vestigial MfeI/NheI sites plus a spare KpnI
    assert len(plan.product) == len(load_plasmid(expected)) - 12


def test_product_keeps_the_cmv_and_loses_the_nfkb_element(phl391, phl162):
    plan = plan_insertion(
        phl391, phl162, insert_features=[CMV_SPAN], at=PCLM2_SITE, method="gibson"
    )
    labels = {feature_label(f) for f in plan.product.features}
    assert "CMV promoter" in labels and "CMV enhancer" in labels
    assert not {name for name in labels if name.upper().startswith("NFKB")}
    assert "miniCMV" not in labels
    for keep in ("mCherry", "Lambda BoxB x8", "AmpR", "ori"):
        assert keep in labels


def test_only_two_insertion_chemistries_remain(phl391, phl162):
    """Restriction when the sites are already there, Gibson when they are not."""
    with pytest.raises(ValueError, match="method must be"):
        plan_insertion(
            phl391, phl162, insert_features=[CMV_SPAN], at=PCLM2_SITE,
            method="pcr_restriction",
        )


# --------------------------------------------------------------------------- #
# product verification -- the standard excision already held itself to
# --------------------------------------------------------------------------- #


def test_a_route_that_destroys_a_protected_feature_is_refused(vector, donor, monkeypatch):
    """``protect`` is now a post-condition, not only a hint for choosing cut sites.

    Before this, ``protected`` steered ``_restriction_route`` and was then dropped: the
    post-assembly checks covered length, insert presence, arm uniqueness and site regeneration,
    but never re-checked that the protected features survived. ``excise._verify`` had always
    done so, and insertion is the more error-prone operation of the two.

    Simulated by corrupting the product after assembly, which is the failure the check exists to
    catch however it arises.
    """
    from bbl import insert as insert_module

    real_replace = insert_module.replace_span

    def corrupting_replace(record, start, end, *args, **kwargs):
        product = real_replace(record, start, end, *args, **kwargs)
        # excise a chunk of the surviving backbone, i.e. destroy something protected
        product.seq = product.seq[:200] + product.seq[400:]
        return product

    monkeypatch.setattr(insert_module, "replace_span", corrupting_replace)
    with pytest.raises(NoInsertionRoute) as raised:
        plan_insertion(vector, donor, insert_features=["Lambda BoxB x8"], at=SITE)
    message = str(raised.value)
    assert "length" in message or "not intact" in message


def test_a_sound_product_verifies_clean(vector, donor):
    """The same check must not fire on the known-good pCLM3 + pCLM1 -> pCLM2 route."""
    plan = plan_insertion(vector, donor, insert_features=["Lambda BoxB x8"], at=SITE)
    assert plan.inserted_bp >= BOXB_BP  # the restriction route carries 6 bp of donor flank
    assert not any("not intact" in warning for warning in plan.warnings)
