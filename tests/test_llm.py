"""Tests for the LLM boundary.

All of these run offline with no API key: the session layer is deliberately free of any
Anthropic dependency, so the whole model-facing contract is testable without a network call.
"""

import json

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
# tool wiring (needs the SDK, but no API key)
# --------------------------------------------------------------------------- #


def test_tools_build_with_usable_schemas(session):
    pytest.importorskip("anthropic")
    from bbl.llm.tools import OUTWARD_FACING, build_tools

    tools = {t.name: t for t in build_tools(session)}
    assert set(tools) == {
        "search_inventory", "inspect_plasmid", "plan_deletion", "plan_insertion",
        "source_sequence", "compare_product", "export_product",
    }
    assert OUTWARD_FACING <= set(tools)
    for name, tool in tools.items():
        schema = tool.to_dict()
        assert len(schema["description"]) > 200, f"{name} is under-described"
        assert schema["input_schema"]["properties"], name


def test_tool_descriptions_say_when_to_call(session):
    pytest.importorskip("anthropic")
    from bbl.llm.tools import build_tools

    for tool in build_tools(session):
        description = tool.to_dict()["description"].lower()
        assert "call this" in description or "use this" in description, tool.name


# --------------------------------------------------------------------------- #
# interactive machinery -- exercised offline with a stub client
# --------------------------------------------------------------------------- #


class _Block:
    def __init__(self, **kw):
        self.__dict__.update(kw)


class _Msg:
    def __init__(self, content):
        self.content = content


class _Runner:
    """Stands in for client.beta.messages.tool_runner."""

    def __init__(self, messages):
        self._messages = messages

    def __iter__(self):
        yield from self._messages

    def generate_tool_call_response(self):
        return None


def _stub_client(messages):
    from types import SimpleNamespace

    return SimpleNamespace(
        beta=SimpleNamespace(
            messages=SimpleNamespace(tool_runner=lambda **kw: _Runner(messages))
        )
    )


def test_harness_appends_the_protocol_verbatim(session):
    """The model names a handle; the harness -- not the model -- renders the protocol."""
    from bbl.llm.tools import design

    handle = session.plan_deletion(PHL391, ["NFKBRE"])["product_id"]
    before = dict(session.products)
    session.products.clear()
    session.products.update(before)  # keep it, but make `design` see it as pre-existing

    reply, _ = design(
        session,
        [{"role": "user", "content": "delete it"}],
        client=_stub_client([_Msg([_Block(type="text", text=f"Use {handle}.")])]),
    )
    # pre-existing products are not re-rendered, only ones created during the turn
    assert "protocol for" not in reply

    fresh = session.plan_deletion("test_pHL391", ["NFKBRE"])["product_id"]
    session.products.pop(fresh)  # simulate: created *during* the next design() call
    def _make_it(**kw):
        session.products[fresh] = before[handle]
        return _Runner([_Msg([_Block(type="text", text=f"Recommend {fresh}.")])])

    from types import SimpleNamespace

    reply, _ = design(
        session,
        [{"role": "user", "content": "again"}],
        client=SimpleNamespace(
            beta=SimpleNamespace(messages=SimpleNamespace(tool_runner=_make_it))
        ),
    )
    assert f"--- protocol for {fresh} ---" in reply
    assert "Double-digest" in reply and "MfeI" in reply


def test_tool_calls_are_surfaced_to_the_ui(session):
    from bbl.llm.tools import design

    seen = []
    design(
        session,
        [{"role": "user", "content": "x"}],
        client=_stub_client(
            [
                _Msg([_Block(type="tool_use", name="plan_deletion",
                             input={"plasmid": "pHL391", "features": ["NFKBRE"]})]),
                _Msg([_Block(type="text", text="done")]),
            ]
        ),
        on_tool=lambda name, args: seen.append((name, args)),
    )
    assert seen == [("plan_deletion", {"plasmid": "pHL391", "features": ["NFKBRE"]})]


def test_declining_an_export_tells_the_model_not_to_retry(plasmid_dir, tmp_path):
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
# credentials -- active session first, API key as fallback, else offer to log in
# --------------------------------------------------------------------------- #


class _FakeResponse:
    status_code = 401
    headers: dict = {}
    request = None


@pytest.fixture
def clean_env(monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    monkeypatch.delenv("ANTHROPIC_AUTH_TOKEN", raising=False)
    return monkeypatch


def test_active_session_wins_over_an_api_key(clean_env):
    """A key present alongside a live session must not shadow the session."""
    from bbl.llm import credentials

    clean_env.setenv("ANTHROPIC_API_KEY", "sk-ant-should-not-be-used")
    clean_env.setattr(credentials, "ant_available", lambda: True)
    clean_env.setattr(credentials, "session_token", lambda profile=None: "sess-tok")

    cred = credentials.resolve()
    assert cred.kind == credentials.SESSION
    assert cred.client.auth_token == "sess-tok" and cred.client.api_key is None
    assert cred.client.default_headers["anthropic-beta"] == credentials.OAUTH_BETA


def test_exported_session_token_is_used_directly(clean_env):
    from bbl.llm import credentials

    clean_env.setenv("ANTHROPIC_AUTH_TOKEN", "exported-tok")
    cred = credentials.resolve()
    assert cred.kind == credentials.SESSION and cred.client.auth_token == "exported-tok"


def test_api_key_is_the_fallback_when_no_session(clean_env):
    from bbl.llm import credentials

    clean_env.setenv("ANTHROPIC_API_KEY", "sk-ant-test")
    clean_env.setattr(credentials, "session_token", lambda profile=None: None)

    cred = credentials.resolve()
    assert cred.kind == credentials.API_KEY
    assert cred.client.api_key == "sk-ant-test"
    assert "no session" in cred.description
    assert not cred.refreshable


def test_offers_login_when_nothing_is_available(clean_env):
    from bbl.llm import credentials

    tokens = iter([None, "after-login"])
    asked, logged_in = [], []
    clean_env.setattr(credentials, "ant_available", lambda: True)
    clean_env.setattr(credentials, "session_token", lambda profile=None: next(tokens))
    clean_env.setattr(credentials, "login", lambda profile=None: logged_in.append(1) or True)

    cred = credentials.resolve(ask=lambda prompt: asked.append(prompt) or True)
    assert cred.kind == credentials.SESSION and cred.client.auth_token == "after-login"
    assert "Log in now?" in asked[0] and logged_in


def test_declining_login_reports_both_options(clean_env):
    from bbl.llm import credentials

    clean_env.setattr(credentials, "ant_available", lambda: True)
    clean_env.setattr(credentials, "session_token", lambda profile=None: None)

    with pytest.raises(credentials.NoCredentials) as excinfo:
        credentials.resolve(ask=lambda prompt: False)
    message = str(excinfo.value)
    assert "ant auth login" in message and "ANTHROPIC_API_KEY" in message


def test_missing_cli_says_how_to_install(clean_env):
    from bbl.llm import credentials

    clean_env.setattr(credentials, "ant_available", lambda: False)
    with pytest.raises(credentials.NoCredentials, match="anthropic-cli"):
        credentials.resolve(ask=lambda prompt: True)


def test_expired_session_token_is_reminted_once(clean_env):
    """A long conversation must not die when the token ages out mid-session."""
    import anthropic

    from bbl.llm import credentials

    clean_env.setattr(credentials, "session_token", lambda profile=None: "fresh")
    cred = credentials.Credential(
        anthropic.Anthropic(auth_token="stale"), credentials.SESSION, "active Claude session"
    )
    calls = []

    def flaky(client):
        calls.append(client.auth_token)
        if len(calls) == 1:
            raise anthropic.AuthenticationError("expired", response=_FakeResponse(), body=None)
        return "ok"

    assert credentials.with_refresh(flaky, cred) == "ok"
    assert calls == ["stale", "fresh"]
    assert cred.client.auth_token == "fresh"      # the caller keeps the refreshed client
    assert "refreshed" in cred.description


def test_an_api_key_failure_is_not_retried(clean_env):
    """Retrying makes no sense for a credential that cannot expire -- it would hide the cause."""
    import anthropic

    from bbl.llm import credentials

    clean_env.setenv("ANTHROPIC_API_KEY", "sk-ant-test")
    cred = credentials.Credential(anthropic.Anthropic(), credentials.API_KEY, "key")
    attempts = []

    def always_fails(client):
        attempts.append(1)
        raise anthropic.AuthenticationError("bad key", response=_FakeResponse(), body=None)

    with pytest.raises(anthropic.AuthenticationError):
        credentials.with_refresh(always_fails, cred)
    assert len(attempts) == 1


def test_describe_reports_the_source_without_a_client(clean_env):
    from bbl.llm import credentials

    clean_env.setattr(credentials, "session_token", lambda profile=None: None)
    assert credentials.describe() == "none"
    clean_env.setenv("ANTHROPIC_API_KEY", "sk-ant-test")
    assert "no session" in credentials.describe()
    clean_env.setattr(credentials, "session_token", lambda profile=None: "tok")
    assert credentials.describe() == "active Claude session"
