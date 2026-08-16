"""Tests for the LLM boundary.

All of these run offline with no credentials: the session layer is deliberately free of any
model-SDK dependency, so the whole model-facing contract is testable without a network call.
The agent layer is driven by a fake connected client rather than a live one.
"""

import asyncio
import json
import os
from pathlib import Path

import pytest

from bbl.llm import DesignSession, Spec, build_system_blocks, inventory_digest

PHL391 = "pHL391"
PCLM1 = "pCLM1"
BOXB_SITE = "1491:1512"


@pytest.fixture(scope="module")
def session(plasmid_dir):
    return DesignSession(plasmid_dir)


# --------------------------------------------------------------------------- #
# the handles-not-payloads invariant
# --------------------------------------------------------------------------- #


def test_no_tool_ever_returns_a_sequence(session):
    """The model must not be able to get DNA into context -- it could then paraphrase it."""
    results = [
        session.search_inventory("BoxB"),
        session.inspect_plasmid(PHL391),
        session.plan_deletion(PHL391, ["NFKBRE"]),
        session.plan_insertion("pCLM3", at=BOXB_SITE, donor=PCLM1,
                               donor_features=["Lambda BoxB x8"]),
    ]
    for payload in results:
        blob = json.dumps(payload, default=str)
        # any run of 40+ ACGT is a sequence leak; primer strings are far shorter
        longest = max(
            (len(chunk) for chunk in _acgt_runs(blob)), default=0
        )
        assert longest < 40, f"{longest} bp of sequence leaked into a tool result"


def _acgt_runs(text):
    run = []
    for char in text.upper():
        if char in "ACGT":
            run.append(char)
        else:
            if run:
                yield "".join(run)
            run = []
    if run:
        yield "".join(run)


def test_products_are_addressed_by_handle(session):
    result = session.plan_deletion(PHL391, ["NFKBRE"])
    handle = result["product_id"]
    assert handle.startswith("prod_")
    assert handle in session.products
    assert len(session.products[handle].record) == result["length_bp"]


# --------------------------------------------------------------------------- #
# plasmid name resolution
# --------------------------------------------------------------------------- #


def test_plasmid_id_beats_prefix(session):
    """'pCLM2' must not be ambiguous against pCLM21/22/23/24 -- it is the name people type."""
    assert session.resolve("pCLM2").label.startswith("pCLM2_")
    assert session.resolve("pCLM21").label.startswith("pCLM21_")


def test_unknown_plasmid_raises_clearly(session):
    with pytest.raises(LookupError, match="no plasmid named"):
        session.resolve("pNOPE999")


def test_genuinely_ambiguous_name_lists_candidates(session):
    with pytest.raises(LookupError, match="matches"):
        session.resolve("pCLM")


# --------------------------------------------------------------------------- #
# the two ground-truth cases, through the tool layer
# --------------------------------------------------------------------------- #


def test_deletion_ground_truth(session):
    result = session.plan_deletion(PHL391, ["NFKBRE"])
    assert result["feasible"] and result["verified"]
    assert result["method"] == "restriction"
    assert set(result["enzymes"]) == {"MfeI", "EcoRI"}
    assert (result["length_bp"], result["deleted_bp"], result["collateral_bp"]) == (5704, 66, 12)
    assert result["pcr_derived_bp"] == 0 and result["new_oligos"] == 0
    assert session.compare_product(result["product_id"], PCLM1)["identical"]


def test_insertion_ground_truth(session):
    result = session.plan_insertion(
        "pCLM3", at=BOXB_SITE, donor=PCLM1, donor_features=["Lambda BoxB x8"]
    )
    assert result["feasible"] and result["verified"]
    assert set(result["enzymes"]) == {"XhoI", "XbaI"}
    assert result["length_bp"] == 6234 and result["directional"] is True
    assert session.compare_product(result["product_id"], "pCLM2")["identical"]


# --------------------------------------------------------------------------- #
# refusals are results, not errors
# --------------------------------------------------------------------------- #


def test_auto_falls_back_to_pcr_without_being_told(session):
    result = session.plan_deletion("test_pHL391", ["NFKBRE"])
    assert result["feasible"] and result["method"] == "pcr"
    assert result["deleted_bp"] == 54 and result["collateral_bp"] == 0
    assert result["new_oligos"] == 2


def test_refusal_carries_reasons_the_model_can_act_on(session):
    result = session.plan_deletion("test_pHL391", ["NFKBRE"], method="restriction")
    assert result["feasible"] is False
    assert result["reasons"]["restriction"]  # a tally, not a bare string
    assert "error" not in result  # a refusal is information, not a tool failure


def test_bad_feature_name_is_an_error_not_a_refusal(session):
    result = session.plan_deletion(PHL391, ["not_a_feature"])
    assert result["feasible"] is False
    assert "No features matched" in result["error"]


# --------------------------------------------------------------------------- #
# sourcing and the ask-the-user gate
# --------------------------------------------------------------------------- #


def test_unknown_sequence_returns_a_question_not_a_guess(session):
    result = session.source_sequence("barcode")
    assert result["route"] == "ask_user"
    assert "Addgene" in result["question"]


def test_repetitive_part_already_in_the_library_is_moved(session):
    boxb = session.resolve(PCLM1).sequence[961:1255]
    result = session.source_sequence("Lambda BoxB x8", boxb)
    assert result["route"] == "from_inventory"
    assert result["complexity"] is None  # never reached the synthesis gate
    assert len(result["donors"]) >= 5


# --------------------------------------------------------------------------- #
# export
# --------------------------------------------------------------------------- #


def test_export_writes_a_readable_genbank(session, tmp_path):
    from bbl import load_plasmid

    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    out = session.export_product(handle, str(tmp_path / "designed.gb"))
    assert out["length_bp"] == 5704
    assert load_plasmid(out["path"]).circular


def test_export_rejects_an_unknown_handle(session):
    assert "error" in session.export_product("prod_999", "/tmp/x.gb")


def test_save_protocol_writes_the_plan_s_own_text(session, tmp_path):
    """The model cannot write this file itself -- the protocol is appended after its message,
    so it never has the text. This tool is how "save the protocol" gets honoured."""
    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    out = session.save_protocol(handle, str(tmp_path / "protocol.md"))
    written = Path(out["path"]).read_text()
    assert written == session.protocol_for(handle)   # verbatim, not paraphrased
    assert "MfeI" in written and "5704 bp" in written


def test_save_protocol_is_gated_and_rejects_unknown_handles(plasmid_dir, tmp_path):
    declined = DesignSession(plasmid_dir, confirm=lambda prompt: False)
    handle = declined.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    assert declined.save_protocol(handle, str(tmp_path / "no.md"))["declined"] is True
    assert not (tmp_path / "no.md").exists()
    assert "error" in declined.save_protocol("prod_999", str(tmp_path / "x.md"))


# --------------------------------------------------------------------------- #
# reports
# --------------------------------------------------------------------------- #


def test_report_is_written_from_the_stored_plan(session, tmp_path):
    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    out = session.generate_report(handle, str(tmp_path / "report.html"))
    assert out["sections"] >= 4
    html = (tmp_path / "report.html").read_text()
    assert session.protocol_for(handle) in html.replace("&#x27;", "'")


def test_a_report_cites_a_comparison_made_earlier(session, tmp_path):
    """The verdict is remembered by the session, so the model cannot misreport it later."""
    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    session.compare_product(handle, PCLM1)
    session.generate_report(handle, str(tmp_path / "r.html"))
    assert "identical to pCLM1" in (tmp_path / "r.html").read_text()


def test_a_report_draws_the_parent_it_was_designed_from(session, tmp_path):
    """Only worth asserting because the plan itself does not carry the parent record."""
    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    session.generate_report(handle, str(tmp_path / "r.html"))
    html = (tmp_path / "r.html").read_text()
    assert "pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8 — 5,770 bp" in html
    assert html.count("<svg") == 4  # both maps and both zooms; figures are not optional (D83)


def test_a_report_names_plasmids_in_full(session, tmp_path):
    """pCLM21 and pCLM24 differ only after the underscore; the short name is an accession."""
    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    session.generate_report(handle, str(tmp_path / "r.html"))
    html = (tmp_path / "r.html").read_text()
    assert "NFKBRE1-miniCMV-mCherry-LambdaBoxBx8" in html


def test_the_session_says_which_concentrations_it_needs(session):
    """The prompt to ask the user is built from the design, not guessed at by the model."""
    handle = session.plan_insertion(
        "pCLM3", at=BOXB_SITE, donor=PCLM1, donor_features=["Lambda BoxB x8"]
    )["product_id"]
    needed = session.dna_needing_concentration(handle)
    assert [item["role"] for item in needed] == ["parent / backbone", "donor"]
    assert all(item["length_bp"] > 0 for item in needed)


def test_supplied_concentrations_fill_the_volumes(session, tmp_path):
    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    label = session.dna_needing_concentration(handle)[0]["plasmid"]
    out = session.generate_report(
        handle, str(tmp_path / "r.html"), concentrations={label: 100}
    )
    assert out["concentrations_supplied"] == [label]
    assert out["concentrations_missing"] == []
    # 1 ug at 100 ng/uL is 10 uL of DNA, leaving 6 uL of water in a 20 uL digest
    assert ">10.00<" in (tmp_path / "r.html").read_text()


def test_a_junk_concentration_is_dropped_rather_than_coerced(session, tmp_path):
    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    label = session.dna_needing_concentration(handle)[0]["plasmid"]
    out = session.generate_report(
        handle, str(tmp_path / "r.html"), concentrations={label: "unknown"}
    )
    assert out["concentrations_supplied"] == []
    assert out["concentrations_missing"] == [label]


def test_report_writing_is_gated_like_any_other_write(plasmid_dir, tmp_path):
    declining = DesignSession(plasmid_dir, confirm=lambda prompt: False)
    handle = declining.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    target = tmp_path / "nope.html"
    assert declining.generate_report(handle, str(target))["declined"] is True
    assert not target.exists()


def test_report_rejects_an_unknown_handle(session):
    assert "error" in session.generate_report("prod_999", "/tmp/x.html")


def test_report_generation_is_declared_outward_facing():
    from bbl.llm.tools import OUTWARD_FACING

    assert "generate_report" in OUTWARD_FACING


def test_a_written_report_is_remembered_for_the_harness_to_announce(session, tmp_path):
    """The model asks for the file; the harness says what landed. Same split as D64."""
    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    before = len(session.reports)
    session.generate_report(handle, str(tmp_path / "r.html"))
    assert len(session.reports) == before + 1
    assert session.reports[-1]["steps"], "the digest is built from the step titles"
    assert session.reports[-1]["path"].endswith("r.html")


def test_a_report_digest_never_repeats_the_protocol(session, tmp_path):
    """It prints beside the D64 protocol block, so overlapping would double it on screen."""
    from bbl.llm.chat import report_digest

    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    session.generate_report(handle, str(tmp_path / "r.html"))
    digest = report_digest(session.reports[-1])
    assert "r.html" in digest
    assert session.protocol_for(handle) not in digest


# --------------------------------------------------------------------------- #
# prompt assembly
# --------------------------------------------------------------------------- #


def test_digest_is_built_from_the_inventory(session):
    digest = inventory_digest(session.entries, session.config)
    assert f"{len(session.entries)} plasmids" in digest
    assert "pHL162_pcDNA3.1_MCS" in digest
    assert "empty base vector" in digest  # config-tagged, not hard-coded


def test_system_blocks_put_the_breakpoint_after_the_digest(session):
    blocks = build_system_blocks(session.entries, session.config)
    assert len(blocks) == 2
    assert "cache_control" not in blocks[0]      # volatile-free but not the breakpoint
    assert blocks[1]["cache_control"] == {"type": "ephemeral"}


def test_system_prompt_states_the_invariants(session):
    from bbl.llm import SYSTEM

    assert "Never write, quote, or reconstruct a DNA sequence" in SYSTEM
    assert "Do not write the bench protocol" in SYSTEM
    assert "is information, not an error" in SYSTEM


# --------------------------------------------------------------------------- #
# structured outputs
# --------------------------------------------------------------------------- #


def test_spec_records_ambiguity_rather_than_guessing():
    spec = Spec(
        backbone="pcDNA3.1",
        parts=["miniCMV", "mCherry", "LambdaBoxBx8", "polyA"],
        marker="AmpR",
        purpose="a control for my sensor",
        ambiguities=["unlisted features: must-be-absent or don't-care?"],
    )
    assert spec.must_not_contain == []      # nothing inferred
    assert spec.ambiguities                 # the question is surfaced instead


def test_spec_schema_tells_the_model_not_to_infer_exclusions():
    described = Spec.model_json_schema()["properties"]["must_not_contain"]["description"]
    assert "Do not infer" in described


# --------------------------------------------------------------------------- #
# tool wiring (needs the SDK, but no credentials)
# --------------------------------------------------------------------------- #

TOOL_NAMES = {
    "search_inventory", "inspect_plasmid", "plan_deletion", "plan_insertion",
    "source_sequence", "compare_product", "export_product", "save_protocol",
    "report_inputs", "generate_report",
}


def test_tools_build_with_usable_schemas(session):
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.tools import OUTWARD_FACING, build_tools

    tools = {t.name: t for t in build_tools(session)}
    assert set(tools) == TOOL_NAMES
    assert OUTWARD_FACING <= set(tools)
    for name, tool in tools.items():
        assert len(tool.description) > 200, f"{name} is under-described"
        assert tool.input_schema, name


def test_tool_descriptions_say_when_to_call(session):
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.tools import build_tools

    for tool in build_tools(session):
        description = tool.description.lower()
        assert "call this" in description or "use this" in description, tool.name


def test_tool_names_reach_the_allowlist_qualified(session):
    """A tool that is defined but not allowlisted is invisible; derive one from the other."""
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.tools import qualified_names, unqualify

    names = qualified_names(session)
    assert set(names) == {f"mcp__bbl__{n}" for n in TOOL_NAMES}
    assert {unqualify(n) for n in names} == TOOL_NAMES
    assert unqualify("Write") == "Write"  # built-ins pass through untouched


def test_handlers_return_mcp_content_not_raw_payloads(session):
    """The SDK expects {"content": [...]}; a bare dict would be silently dropped."""
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.tools import build_tools

    tools = {t.name: t for t in build_tools(session)}
    result = asyncio.run(tools["inspect_plasmid"].handler({"plasmid": PHL391}))
    assert result["content"][0]["type"] == "text"
    assert json.loads(result["content"][0]["text"])["length_bp"] == 5770


def test_the_report_handler_round_trips_a_concentrations_dict(session, tmp_path):
    """The only tool taking a dict argument -- a schema that flattened it would be silent."""
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.tools import build_tools

    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    label = session.dna_needing_concentration(handle)[0]["plasmid"]
    tools = {t.name: t for t in build_tools(session)}
    result = asyncio.run(
        tools["generate_report"].handler(
            {
                "product_id": handle,
                "path": str(tmp_path / "r.html"),
                "concentrations": {label: 100},
            }
        )
    )
    payload = json.loads(result["content"][0]["text"])
    assert payload["concentrations_supplied"] == [label]
    assert payload["concentrations_missing"] == []


# --------------------------------------------------------------------------- #
# the harness renders the protocol -- exercised offline with a fake message stream
# --------------------------------------------------------------------------- #


class _Block:
    def __init__(self, **kw):
        self.__dict__.update(kw)


class _Msg:
    def __init__(self, content):
        self.content = content


class _Result:
    """Stands in for the SDK's ResultMessage."""

    content: list = []

    def __init__(self, session_id="sess-1", num_turns=3, is_error=False, errors=None):
        self.session_id = session_id
        self.num_turns = num_turns
        self.is_error = is_error
        self.errors = errors
        self.total_cost_usd = 0.12
        self.duration_ms = 4200


def test_harness_appends_the_protocol_verbatim(session):
    """The model names a handle; the harness -- not the model -- renders the protocol."""
    from bbl.llm.agent import new_protocols

    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]

    # a product that already existed before the turn is not re-rendered
    assert new_protocols(session, set(session.products), f"Use {handle}.") == ""

    # one created during the turn, and named in the reply, is -- verbatim from the plan
    before = set(session.products)
    fresh = session.plan_deletion("test_pHL391", ["NFKBRE"])["product_id"]
    rendered = new_protocols(session, before, f"Recommend {fresh}.")
    assert f"--- protocol for {fresh} ---" in rendered
    assert rendered.endswith(session.protocol_for(fresh))
    assert "5698 bp" in rendered  # a number the model never got the chance to paraphrase

    # created but never named: still rendered, because a good answer can describe the route
    # without typing the handle, and dropping the protocol from that turn is the worse failure
    assert f"--- protocol for {fresh} ---" in new_protocols(session, before, "here is the route")

    # but when the model *does* name one of several, only that one is rendered
    before = set(session.products)
    a = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    b = session.plan_deletion(PHL391, ["NFKBRE"], method="pcr")["product_id"]
    only_a = new_protocols(session, before, f"I recommend {a} over the alternative.")
    assert f"protocol for {a}" in only_a and f"protocol for {b}" not in only_a


def test_renderer_streams_tool_calls_and_collects_the_text(capsys):
    from bbl.llm.render import Renderer, Style

    renderer = Renderer(Style(enabled=False))
    for message in [
        _Msg([_Block(type="thinking", thinking="weighing the two routes")]),
        _Msg([_Block(type="tool_use", name="mcp__bbl__plan_deletion",
                     input={"plasmid": "pHL391", "features": ["NFKBRE"]})]),
        _Msg([_Block(type="tool_result", tool_use_id="t1", content="x" * 412)]),
        _Msg([_Block(type="text", text="Delete NFKBRE with MfeI + EcoRI.")]),
        _Result(),
    ]:
        renderer.handle(message)

    out = capsys.readouterr().out
    assert "· weighing the two routes" in out
    assert "→ plan_deletion: pHL391" in out          # qualified name is stripped for the user
    assert "✓ 412 chars" in out                      # results are sized, never inlined
    assert "3 turns" in out and "$0.12" in out
    assert renderer.text == "Delete NFKBRE with MfeI + EcoRI."
    assert renderer.tool_calls[0][0] == "plan_deletion"


def test_renderer_never_inlines_a_tool_result(capsys):
    """Tool payloads carry sequence; printing them would put DNA in front of the model's user."""
    from bbl.llm.render import Renderer, Style

    renderer = Renderer(Style(enabled=False))
    renderer.handle(_Msg([_Block(type="tool_result", tool_use_id="t", content="ACGT" * 50)]))
    assert "ACGTACGT" not in capsys.readouterr().out


def test_style_is_plain_when_not_a_terminal():
    from bbl.llm.render import Style

    assert Style(enabled=False).dim == "" and Style(enabled=True).dim


# --------------------------------------------------------------------------- #
# the confirmation gate
# --------------------------------------------------------------------------- #


def _gate_decision(gate, tool, args):
    return asyncio.run(gate(tool, args, None))


def test_an_outward_facing_tool_is_never_allowlisted(session, tmp_path):
    """Allowlisting a whole tool auto-approves it *before* can_use_tool runs, so an entry for
    export_product would silently kill the confirmation gate."""
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import build_options
    from bbl.llm.tools import OUTWARD_FACING

    allowed = build_options(session, env={}, cwd=tmp_path).allowed_tools
    for name in OUTWARD_FACING:
        assert f"mcp__bbl__{name}" not in allowed
    assert "mcp__bbl__search_inventory" in allowed  # read-only tools still run unattended


def test_write_is_never_allowlisted(session, tmp_path):
    """Path-scoped Write(...) entries are silently ignored in allowed_tools (measured against
    CLI 2.1.x), so the working-directory scope is enforced in the gate instead. A bare entry
    here would allow writing anywhere and bypass it."""
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import build_options

    allowed = build_options(session, env={}, cwd=tmp_path).allowed_tools
    assert not any(entry.startswith("Write") for entry in allowed)


def test_the_agent_sees_only_bbl_s_own_mcp_servers(session, tmp_path):
    """Without this the agent picks up the user's account connectors and reasons about them."""
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import build_options

    options = build_options(session, env={}, cwd=tmp_path)
    assert options.strict_mcp_config is True
    assert set(options.mcp_servers) == {"bbl"}


def test_a_write_inside_the_working_directory_is_not_gated(tmp_path):
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import confirmation_gate

    gate = confirmation_gate(lambda q: pytest.fail(f"should not ask: {q}"), tmp_path)
    for path in (str(tmp_path / "protocol.md"), "protocol.md", "sub/dir/out.gb"):
        assert _gate_decision(gate, "Write", {"file_path": path}).behavior == "allow"


def test_a_write_outside_the_working_directory_is_gated(tmp_path):
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import confirmation_gate

    asked = []
    gate = confirmation_gate(lambda q: asked.append(q) or False, tmp_path)
    for path in ("/etc/passwd", "../escape.txt", str(tmp_path / ".." / "sneaky.txt")):
        assert _gate_decision(gate, "Write", {"file_path": path}).behavior == "deny"
    assert len(asked) == 3
    assert all("outside the working directory" in question for question in asked)


def test_a_confirmed_write_outside_still_goes_through(tmp_path):
    """The gate asks; it does not forbid. Saying yes is the user's call."""
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import confirmation_gate

    gate = confirmation_gate(lambda q: True, tmp_path)
    assert _gate_decision(gate, "Write", {"file_path": "/tmp/x.txt"}).behavior == "allow"


def test_the_agent_does_not_inherit_the_user_s_settings(session, tmp_path):
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import build_options

    # the SDK's default of None loads user settings *and* memories
    assert build_options(session, env={}, cwd=tmp_path).setting_sources == []


def test_gate_stops_an_export_and_tells_the_model_not_to_retry(tmp_path):
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import DECLINED, confirmation_gate

    asked = []
    gate = confirmation_gate(lambda q: asked.append(q) or False, tmp_path)
    args = {"product_id": "prod_1", "path": "out.gb"}
    decision = _gate_decision(gate, "mcp__bbl__export_product", args)
    assert decision.behavior == "deny" and DECLINED in decision.message
    assert str(tmp_path / "out.gb") in asked[0]  # the resolved path, not the bare argument
    assert "prod_1" in asked[0]


def test_gate_stops_a_report_the_model_asked_for(tmp_path):
    """The report is the real enforcement point now: it is what the model reaches for."""
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import DECLINED, confirmation_gate

    asked = []
    gate = confirmation_gate(lambda q: asked.append(q) or False, tmp_path)
    args = {"product_id": "prod_1", "path": "r.html"}
    decision = _gate_decision(gate, "mcp__bbl__generate_report", args)
    assert decision.behavior == "deny" and DECLINED in decision.message
    assert str(tmp_path / "r.html") in asked[0]


def test_gate_does_not_stand_between_the_model_and_report_inputs(tmp_path):
    """Asking which concentrations are needed writes nothing, so it must not stop to confirm."""
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import confirmation_gate

    gate = confirmation_gate(lambda q: False, tmp_path)
    decision = _gate_decision(gate, "mcp__bbl__report_inputs", {"product_id": "prod_1"})
    assert decision.behavior == "allow"


def test_gate_lets_a_confirmed_export_through(tmp_path):
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import confirmation_gate

    gate = confirmation_gate(lambda q: True, tmp_path)
    assert _gate_decision(gate, "mcp__bbl__export_product", {"path": "y.gb"}).behavior == "allow"


def test_gate_does_not_stand_between_the_model_and_a_read_only_tool(tmp_path):
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm.agent import confirmation_gate

    gate = confirmation_gate(lambda q: pytest.fail("read-only tools must not prompt"), tmp_path)
    for tool in ("mcp__bbl__search_inventory", "mcp__bbl__plan_deletion", "Read"):
        assert _gate_decision(gate, tool, {}).behavior == "allow"


def test_declining_an_export_tells_the_model_not_to_retry(plasmid_dir, tmp_path):
    """The session layer keeps its own gate, so it stays usable outside the agent."""
    declined = DesignSession(plasmid_dir, confirm=lambda prompt: False)
    handle = declined.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    result = declined.export_product(handle, str(tmp_path / "nope.gb"))
    assert result["declined"] is True
    assert "do not retry" in result["note"]
    assert not (tmp_path / "nope.gb").exists()


def test_confirming_an_export_writes_the_file(plasmid_dir, tmp_path):
    prompts = []
    approved = DesignSession(plasmid_dir, confirm=lambda p: prompts.append(p) or True)
    handle = approved.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    out = approved.export_product(handle, str(tmp_path / "yes.gb"))
    assert out["length_bp"] == 5704
    assert (tmp_path / "yes.gb").exists()
    assert "5704 bp" in prompts[0]


def test_products_persist_across_turns(session):
    """Handles designed early stay addressable, so routes can be compared later."""
    first = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    second = session.plan_deletion(PHL391, ["NFKBRE"], method="pcr")["product_id"]
    assert first != second
    assert {first, second} <= set(session.products)
    assert session.compare_product(first, PCLM1)["identical"]


# --------------------------------------------------------------------------- #
# one turn, over a fake connected client
# --------------------------------------------------------------------------- #


class _FakeClient:
    """Stands in for ClaudeSDKClient: records the prompt, replays a message stream."""

    def __init__(self, messages, raise_after=None):
        self.messages = messages
        self.raise_after = raise_after
        self.prompts = []

    async def query(self, prompt):
        self.prompts.append(prompt)

    async def receive_response(self):
        for message in self.messages:
            yield message
        if self.raise_after is not None:
            raise self.raise_after


def _run_turn(client):
    from bbl.llm.agent import run_turn
    from bbl.llm.render import Renderer, Style

    return asyncio.run(run_turn(client, "delete NFKBRE", Renderer(Style(enabled=False))))


def test_a_turn_reports_its_text_and_session_id():
    client = _FakeClient([_Msg([_Block(type="text", text="Use MfeI + EcoRI.")]), _Result()])
    turn = _run_turn(client)
    assert turn.text == "Use MfeI + EcoRI."
    assert turn.session_id == "sess-1" and not turn.is_error
    assert client.prompts == ["delete NFKBRE"]


def test_a_trailing_process_error_does_not_discard_a_delivered_result():
    """The CLI exits non-zero *after* reporting a cap breach. The result is still real."""
    client = _FakeClient(
        [_Msg([_Block(type="text", text="partial answer")]), _Result(is_error=True)],
        raise_after=RuntimeError("exited with code 1"),
    )
    turn = _run_turn(client)
    assert turn.is_error and turn.text == "partial answer"
    assert turn.session_id == "sess-1"
    assert "exited with code 1" in turn.error


def test_an_error_before_any_result_still_raises():
    client = _FakeClient([], raise_after=RuntimeError("connection refused"))
    with pytest.raises(RuntimeError, match="connection refused"):
        _run_turn(client)


def test_an_interrupted_turn_is_still_resumable():
    """The session id must come off any message, not just the result -- a turn cut short
    before its ResultMessage is exactly the one worth resuming."""
    from bbl.llm.render import Renderer, Style

    renderer = Renderer(Style(enabled=False))
    renderer.handle(_Msg([_Block(type="text", text="thinking out loud")]))
    assert renderer.session_id is None

    init = _Msg([])
    init.session_id = "sess-early"
    renderer.handle(init)
    assert renderer.session_id == "sess-early"


def test_resuming_remembers_the_last_conversation(tmp_path):
    from bbl.llm.chat import last_session_id, save_session_id

    assert last_session_id(tmp_path) is None
    save_session_id("sess-42", tmp_path)
    assert last_session_id(tmp_path) == "sess-42"


# --------------------------------------------------------------------------- #
# auth -- the subscription first, an API key as fallback, else a loud failure
# --------------------------------------------------------------------------- #


@pytest.fixture
def clean_env(monkeypatch, tmp_path):
    """No inherited credentials, and a config dir under our control."""
    for var in ("ANTHROPIC_API_KEY", "ANTHROPIC_AUTH_TOKEN", "CLAUDE_CODE_OAUTH_TOKEN",
                "CLAUDE_CODE_USE_BEDROCK", "CLAUDE_CODE_USE_VERTEX"):
        monkeypatch.delenv(var, raising=False)
    config = tmp_path / "claude-config"
    config.mkdir()
    monkeypatch.setenv("CLAUDE_CONFIG_DIR", str(config))
    return monkeypatch


def _write_credentials(payload):
    path = Path(os.environ["CLAUDE_CONFIG_DIR"]) / ".credentials.json"
    path.write_text(json.dumps(payload))
    return path


def test_a_subscription_login_beats_an_api_key(clean_env):
    """A key present alongside a subscription must not shadow it -- the plan is prepaid."""
    from bbl.llm import auth

    _write_credentials({"claudeAiOauth": {"accessToken": "x"}})
    clean_env.setenv("ANTHROPIC_API_KEY", "sk-ant-should-not-be-used")
    assert auth.resolve_auth_mode("auto") == "session"


def test_the_oauth_token_variable_counts_as_a_subscription(clean_env):
    from bbl.llm import auth

    clean_env.setenv("CLAUDE_CODE_OAUTH_TOKEN", "portable-tok")
    assert auth.session_auth_available()
    assert auth.resolve_auth_mode("auto") == "session"


def test_an_api_key_is_the_fallback_when_there_is_no_subscription(clean_env):
    from bbl.llm import auth

    clean_env.setenv("ANTHROPIC_API_KEY", "sk-ant-test")
    assert not auth.session_auth_available()
    assert auth.resolve_auth_mode("auto") == "api"


def test_an_api_only_credentials_file_is_not_a_subscription(clean_env):
    """The claudeAiOauth object is the signal, not the file's existence."""
    from bbl.llm import auth

    _write_credentials({"someOtherKey": {}})
    assert not auth.session_auth_available()
    with pytest.raises(auth.AuthError, match="no Claude credentials found"):
        auth.resolve_auth_mode("auto")


def test_an_unreadable_credentials_file_is_not_a_subscription(clean_env):
    from bbl.llm import auth

    (Path(os.environ["CLAUDE_CONFIG_DIR"]) / ".credentials.json").write_text("{not json")
    assert not auth.session_auth_available()


def test_an_empty_variable_does_not_count(clean_env):
    from bbl.llm import auth

    clean_env.setenv("ANTHROPIC_API_KEY", "")
    clean_env.setenv("CLAUDE_CODE_OAUTH_TOKEN", "")
    assert not auth.api_auth_available() and not auth.session_auth_available()


def test_forcing_a_mode_requires_that_mode_s_credential(clean_env):
    from bbl.llm import auth

    clean_env.setenv("ANTHROPIC_API_KEY", "sk-ant-test")
    assert auth.resolve_auth_mode("api") == "api"
    with pytest.raises(auth.AuthError, match="no Claude subscription login found"):
        auth.resolve_auth_mode("session")

    _write_credentials({"claudeAiOauth": {"accessToken": "x"}})
    clean_env.delenv("ANTHROPIC_API_KEY")
    assert auth.resolve_auth_mode("session") == "session"
    with pytest.raises(auth.AuthError, match="no API credential found"):
        auth.resolve_auth_mode("api")


def test_no_credentials_at_all_says_how_to_fix_it(clean_env):
    from bbl.llm import auth

    with pytest.raises(auth.AuthError) as excinfo:
        auth.resolve_auth_mode("auto")
    message = str(excinfo.value)
    assert "`claude`" in message and "ANTHROPIC_API_KEY" in message


def test_the_bundled_cli_counts_even_with_nothing_on_path(clean_env, monkeypatch):
    """The SDK prefers its own bundled binary over PATH, so `pip install bbl[claude]` alone is
    a complete install. Gating on shutil.which would reject a setup that works."""
    pytest.importorskip("claude_agent_sdk")
    from bbl.llm import auth

    monkeypatch.setattr(auth.shutil, "which", lambda name: None)
    import claude_agent_sdk

    bundled = Path(claude_agent_sdk.__file__).parent / "_bundled" / "claude"
    assert auth.claude_cli_found() is bundled.is_file()


def test_a_missing_transport_says_which_half_is_missing(clean_env, monkeypatch):
    from bbl.llm import auth

    monkeypatch.setattr(auth, "claude_sdk_available", lambda: False)
    with pytest.raises(auth.AuthError, match="bbl\\[claude\\]"):
        auth.check_agent_cli()

    monkeypatch.setattr(auth, "claude_sdk_available", lambda: True)
    monkeypatch.setattr(auth, "claude_cli_found", lambda: False)
    with pytest.raises(auth.AuthError, match="no `claude` CLI available"):
        auth.check_agent_cli()


def test_a_cloud_provider_flag_counts_as_an_api_credential(clean_env):
    from bbl.llm import auth

    clean_env.setenv("CLAUDE_CODE_USE_BEDROCK", "1")
    assert auth.api_auth_available()
    assert auth.resolve_auth_mode("auto") == "api"


# --------------------------------------------------------------------------- #
# the environment handed to the subprocess
# --------------------------------------------------------------------------- #


def test_the_unused_credential_is_blanked_not_dropped(clean_env):
    """The SDK merges our mapping over os.environ, so omitting a key leaves it live."""
    from bbl.llm import auth

    clean_env.setenv("ANTHROPIC_API_KEY", "sk-ant-should-not-be-used")
    clean_env.setenv("CLAUDE_CODE_OAUTH_TOKEN", "tok")

    session_env = auth.agent_env("session")
    assert session_env["ANTHROPIC_API_KEY"] == ""      # present and empty, not absent
    assert "ANTHROPIC_API_KEY" in session_env
    assert session_env["CLAUDE_CODE_OAUTH_TOKEN"] == "tok"

    api_env = auth.agent_env("api")
    assert api_env["CLAUDE_CODE_OAUTH_TOKEN"] == ""
    assert "CLAUDE_CODE_OAUTH_TOKEN" in api_env
    assert api_env["ANTHROPIC_API_KEY"] == "sk-ant-should-not-be-used"


def test_the_agent_does_not_inherit_the_user_s_claude_mds(clean_env):
    """Project discovery walks *up* from cwd; bbl runs wherever the plasmids are."""
    from bbl.llm import auth

    assert auth.agent_env("session")["CLAUDE_CODE_DISABLE_CLAUDE_MDS"] == "1"


def test_the_environment_is_otherwise_carried_through(clean_env):
    """bbl is not a sandbox -- the agent needs the real PATH, proxy and TLS settings."""
    from bbl.llm import auth

    clean_env.setenv("HTTPS_PROXY", "http://proxy.example:8080")
    env = auth.agent_env("session")
    assert env["HTTPS_PROXY"] == "http://proxy.example:8080"
    assert env["PATH"]
