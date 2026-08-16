"""Reagent arithmetic, checked against a real notebook entry.

The numbers below are lifted from the eCLM24 cloning entry (June 2026) -- a spreadsheet an
experimentalist filled in by hand and then worked from at the bench. If our calculators
reproduce those volumes to the microlitre, they are doing the same arithmetic the lab does.
That is a stronger test than asserting against numbers we made up ourselves.
"""

import pytest

from bbl.report.bench import (
    BLANK,
    MW_PER_BP,
    MW_PER_BP_LIGATION,
    Component,
    bench_config,
    digestion_table,
    dilution_plan,
    fmol_to_ng,
    gel_percent,
    gibson_table,
    growth_temperature,
    ligation_table,
    ng_to_fmol,
    pcr_table,
    template_dilution,
    thermocycler_table,
)


@pytest.fixture
def cfg():
    return bench_config()


def _cell(table, row_contains, column):
    for row in table.rows:
        if row_contains in row[0]:
            return row[column]
    raise AssertionError(f"no row matching {row_contains!r} in {table.title!r}")


def _row(reaction, name_contains):
    """One line of a reaction mix, by substring of its component name."""
    for row in reaction.rows:
        if name_contains in row.name:
            return row
    raise AssertionError(f"no row matching {name_contains!r} in {reaction.title!r}")


# --------------------------------------------------------------------------- #
# molarity
# --------------------------------------------------------------------------- #


def test_ng_and_fmol_round_trip():
    assert ng_to_fmol(fmol_to_ng(42, 1500), 1500) == pytest.approx(42)


def test_hundred_ng_of_the_pHL391_backbone_is_twentyseven_fmol():
    assert ng_to_fmol(100, 5770) == pytest.approx(26.66, abs=0.01)


# --------------------------------------------------------------------------- #
# Gibson: reproduces the lab's own sheet
# --------------------------------------------------------------------------- #


def test_gibson_reproduces_the_notebook_volumes(cfg):
    """eCLM24, pCLM22 intermediate: 5770 bp backbone at 43.6 ng/uL, 137 bp insert at 10.

    The entry's own figures are 2.29 uL of backbone, 0.71 uL of insert, 1.99 uL of water
    into a 10 uL reaction with 5 uL of master mix.
    """
    mix = gibson_table(Component("pHL391", 5770, 43.6), [Component("gCLM_4", 137, 10)], cfg)
    assert _row(mix, "pHL391").volume_ul == pytest.approx(2.29, abs=0.005)
    assert _row(mix, "gCLM_4").volume_ul == pytest.approx(0.71, abs=0.005)
    assert mix.water_ul == pytest.approx(1.99, abs=0.005)
    assert _row(mix, "MasterMix").volume_ul == pytest.approx(5.0)
    assert mix.total_ul == 10


def test_gibson_insert_mass_follows_the_molar_excess(cfg):
    """7.12 ng of a 137 bp insert is three times the backbone's molarity, not its mass."""
    mix = gibson_table(Component("pHL391", 5770, 43.6), [Component("gCLM_4", 137, 10)], cfg)
    assert _row(mix, "gCLM_4").ng == pytest.approx(7.12, abs=0.01)
    assert _row(mix, "gCLM_4").fmol == pytest.approx(3 * _row(mix, "pHL391").fmol)


def test_gibson_leaves_volumes_blank_when_the_stock_is_unmeasured(cfg):
    mix = gibson_table(Component("backbone", 5770), [Component("insert", 300)], cfg)
    assert _row(mix, "backbone").volume_ul is None
    assert mix.water_ul is None
    # the requirement is still stated -- the row is a form, not a shrug
    assert _row(mix, "backbone").ng == 100


# --------------------------------------------------------------------------- #
# ligation: also from the notebook
# --------------------------------------------------------------------------- #


def test_ligation_reproduces_the_notebook_backbone_volume(cfg):
    """eCLM24, final pCLM24: 7462 bp backbone at 41.1 ng/uL, 15 fmol -> 1.797 uL."""
    mix = ligation_table(
        Component("pCLM24_intermediate", 7462, 41.1), Component("pHL391", 300, 3.3), cfg
    )
    assert _row(mix, "pCLM24_intermediate").fmol == 15
    assert _row(mix, "pCLM24_intermediate").volume_ul == pytest.approx(1.797, abs=0.005)


def test_ligation_insert_is_exactly_the_stated_ratio(cfg):
    """The notebook's own sheet says 1:5 and then uses 80 fmol against a 15 fmol backbone.

    That is 1:5.33 -- a rounding-up in the spreadsheet, not a different rule. We take the
    stated ratio literally (75 fmol, 4.50 uL); the arithmetic underneath is identical, as the
    matching backbone volume above shows, and 80 fmol reproduces the sheet's 4.80 uL exactly.
    """
    mix = ligation_table(
        Component("pCLM24_intermediate", 7462, 41.1), Component("pHL391", 300, 3.3), cfg
    )
    assert _row(mix, "pHL391").fmol == 75
    assert _row(mix, "pHL391").volume_ul == pytest.approx(4.50, abs=0.005)
    assert fmol_to_ng(80, 300, MW_PER_BP_LIGATION) / 3.3 == pytest.approx(4.80, abs=0.005)


def test_the_two_molecular_weights_stay_where_they_belong():
    """650 for Gibson, 660 for ligation -- each matching the sheet it came from."""
    assert MW_PER_BP == 650
    assert fmol_to_ng(15, 7462, MW_PER_BP_LIGATION) == pytest.approx(73.87, abs=0.01)


# --------------------------------------------------------------------------- #
# digestion
# --------------------------------------------------------------------------- #


def test_digest_volumes_sum_to_the_reaction(cfg):
    mix = digestion_table(Component("pHL391", 5770, 100), ["MfeI", "EcoRI"], cfg)
    assert _row(mix, "pHL391").volume_ul == pytest.approx(10.0)  # 1 ug at 100 ng/uL
    assert mix.water_ul == pytest.approx(20.0 - 10.0 - 2 - 1 - 1)
    assert sum(r.volume_ul for r in mix.rows) + mix.water_ul == pytest.approx(mix.total_ul)


def test_digest_adds_phosphatase_only_when_asked(cfg):
    plain = digestion_table(Component("p", 100, 50), ["MfeI"], cfg)
    treated = digestion_table(Component("p", 100, 50), ["MfeI"], cfg, dephosphorylate=True)
    assert not any("FastAP" in row.name for row in plain.rows)
    assert any("FastAP" in row.name for row in treated.rows)
    # the extra microlitre comes out of the water, not out of the total
    assert treated.water_ul == pytest.approx(plain.water_ul - 1)
    assert treated.total_ul == plain.total_ul


def test_digest_of_an_unmeasured_stock_is_blank_not_guessed(cfg):
    mix = digestion_table(Component("donor plasmid"), ["XhoI"], cfg)
    assert _row(mix, "donor plasmid").volume_ul is None
    assert mix.water_ul is None


# --------------------------------------------------------------------------- #
# gel, PCR, growth
# --------------------------------------------------------------------------- #


def test_gel_percentage_tracks_the_band_you_are_chasing():
    assert gel_percent(300) == 1.5  # the notebook switched to 1.5% for a 300 bp fragment
    assert gel_percent(700) == 1.2
    assert gel_percent(5704) == 1.0


def test_pcr_reaction_sums_to_its_stated_volume(cfg):
    table = pcr_table(cfg)
    parts = [float(row[1]) for row in table.rows[:-1]]
    assert sum(parts) == pytest.approx(float(table.rows[-1][1]))


def test_extension_time_lengthens_for_a_repeat_array(cfg):
    plain = thermocycler_table(60, 5698, cfg)
    repeats = thermocycler_table(60, 5698, cfg, repetitive=True)
    assert _cell(plain, "Extension", 2) == "171 s"
    assert _cell(repeats, "Extension", 2) == "256 s"
    assert "repeat array" in repeats.note


def test_annealing_temperature_is_printed_as_given(cfg):
    assert _cell(thermocycler_table(61.6, 1000, cfg), "Annealing", 1) == "62 °C"


# --------------------------------------------------------------------------- #
# diluting the template
# --------------------------------------------------------------------------- #


def _plan(stock, cfg, target=1):
    pcr = cfg["pcr"]
    return dilution_plan(stock, target, pcr["dilution_ladder_ul"], pcr["min_pipette_ul"])


def test_a_dilute_stock_needs_one_step(cfg):
    plan = _plan(20, cfg)
    assert plan["status"] == "ok"
    assert len(plan["steps"]) == 1
    step = plan["steps"][0]
    assert (step["take_ul"], step["water_ul"], step["final_ul"]) == (1.0, 19.0, 20.0)


def test_the_smallest_workable_volume_wins(cfg):
    """Nothing is made up to 100 uL when 10 uL would do -- that is wasted stock and water."""
    assert _plan(10, cfg)["steps"][0]["final_ul"] == 10.0


def test_a_concentrated_miniprep_gets_a_serial_dilution(cfg):
    """0.04 uL of a 500 ng/uL stock is a number, not an action. Two steps, both pipettable."""
    plan = _plan(500, cfg)
    assert [round(s["gives"], 3) for s in plan["steps"]] == [5.0, 1.0]
    assert [s["source"] for s in plan["steps"]] == ["stock", "step 1"]
    assert plan["steps"][0]["take_ul"] == 1.0


@pytest.mark.parametrize("stock", [1.1, 2, 5, 20, 87.3, 250, 500, 1364, 2000, 9999])
def test_no_step_ever_asks_for_an_unpipettable_volume(cfg, stock):
    """The property that makes this box worth having at all."""
    plan = _plan(stock, cfg)
    assert plan["status"] == "ok"
    for step in plan["steps"]:
        assert step["take_ul"] >= cfg["pcr"]["min_pipette_ul"]
        assert 0 <= step["water_ul"] <= step["final_ul"]
        assert step["take_ul"] + step["water_ul"] == pytest.approx(step["final_ul"])


@pytest.mark.parametrize("stock", [1.1, 20, 500, 1364, 9999])
def test_the_last_step_lands_on_the_target(cfg, stock):
    """Every step is a real dilution: what goes in over what it is made up to."""
    source = stock
    for step in _plan(stock, cfg)["steps"]:
        assert source * step["take_ul"] / step["final_ul"] == pytest.approx(step["gives"])
        source = step["gives"]
    assert source == pytest.approx(cfg["pcr"]["template_ng_per_ul"])


def test_an_unmeasured_stock_is_a_form_not_a_guess(cfg):
    assert _plan(None, cfg) == {"status": "unknown", "steps": []}
    assert _plan(0, cfg) == {"status": "unknown", "steps": []}


def test_a_stock_already_at_the_target_is_used_neat(cfg):
    assert _plan(1, cfg)["status"] == "neat"
    assert _plan(0.4, cfg)["status"] == "neat"


def test_the_box_picks_up_a_measured_concentration(cfg):
    box = template_dilution(Component("pHL391", 5770, 1364), cfg)
    assert box.ng_per_ul == 1364
    assert box.plan()["status"] == "ok"
    assert template_dilution(Component("pHL391", 5770), cfg).plan()["status"] == "unknown"


def test_the_lab_can_move_the_target_and_the_ladder():
    cfg = bench_config({"bench": {"pcr": {"template_ng_per_ul": 0.1,
                                          "dilution_ladder_ul": [50]}}})
    plan = template_dilution(Component("p", 5000, 50), cfg).plan()
    assert plan["steps"][-1]["gives"] == pytest.approx(0.1)
    assert plan["steps"][-1]["final_ul"] == 50


def test_repeat_arrays_are_grown_cool():
    cfg = bench_config()
    assert growth_temperature(cfg, repetitive=True)[0] == 30
    assert growth_temperature(cfg, repetitive=False)[0] == 37


# --------------------------------------------------------------------------- #
# configuration
# --------------------------------------------------------------------------- #


def test_lab_config_overrides_one_key_without_losing_the_rest():
    cfg = bench_config({"bench": {"gibson": {"total_ul": 20}}})
    assert cfg["gibson"]["total_ul"] == 20
    assert cfg["gibson"]["insert_molar_excess"] == 3  # untouched default survives
    assert cfg["digest"]["total_ul"] == 20


def test_a_larger_gibson_reaction_scales_its_master_mix():
    cfg = bench_config({"bench": {"gibson": {"total_ul": 20}}})
    mix = gibson_table(Component("b", 5000, 50), [Component("i", 500, 10)], cfg)
    assert _row(mix, "MasterMix").volume_ul == pytest.approx(10.0)
