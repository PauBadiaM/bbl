"""The design tools, as an in-process SDK MCP server.

The wrappers are thin: they serialise a :class:`DesignSession` method to JSON and carry a
description written for the model. Descriptions say *when* to call, not just what the tool does
-- that is what actually drives tool choice, so they live in module-level constants where they
can be reviewed and tested as prose rather than buried in decorators.

The server runs in bbl's own process, so the handlers close over the live
:class:`DesignSession` and product handles stay valid across turns. Importing this module
requires ``claude-agent-sdk``; ``session.py`` does not.

Every design handler is pure-CPU biopython work that neither blocks nor touches the terminal;
the confirmation gate for the outward-facing ones lives in ``agent.py``'s ``can_use_tool``
callback. ``ask_user`` is the one deliberate exception -- it waits on the person at the
keyboard, for as long as they take (docs/DECISIONS.md D84) -- and it does that on a worker
thread so the turn keeps streaming meanwhile.
"""

from __future__ import annotations

import asyncio
import json

from .session import DesignSession

#: Tools that change something outside the process. Gated before execution.
OUTWARD_FACING = {"export_product", "save_protocol", "generate_report"}

#: The MCP server name. Tool names reach ``allowed_tools`` as ``mcp__bbl__<tool>``.
SERVER = "bbl"


SEARCH_INVENTORY = """\
Search the plasmid library by name, feature label, or exact DNA sequence.

Call this first for any request that names a plasmid or a part, and whenever you need to know
whether a sequence already exists somewhere in the library -- a part that is already on a
plasmid can be moved rather than synthesised, which sidesteps the synthesis complexity limit
entirely. A DNA string of 12 bp or more is searched as a sequence on both strands.

query: A plasmid name ("pHL391"), a feature label ("BoxB"), or a DNA sequence.
"""

INSPECT_PLASMID = """\
List a plasmid's length, functional features with coordinates, and unique cutters.

Call this before planning an edit, to see what is present and where. Feature spans from here
are what you pass to plan_deletion and plan_insertion. It is also the cheapest way to check an
assumption about a construct before committing to a route that depends on it.

plasmid: Plasmid name, e.g. "pHL391". A prefix is enough if it is unambiguous.
"""

PLAN_DELETION = """\
Design a way to delete features from a plasmid, and verify the product.

Call this whenever the user wants something removed from a construct -- an unwanted response
element, a promoter being swapped out, a spacer. Also call it when a candidate backbone is right
except for extra features: deleting them is usually far cheaper than building the construct from
parts.

Use "auto" (the default) unless the user asked for a specific chemistry: it tries restriction
first and falls back to inverse PCR, which is the preference order you should follow anyway. A
result with feasible=false carries `reasons` explaining what was rejected -- read it and
consider the other method rather than calling again with the same arguments.

plasmid: Plasmid to edit, e.g. "pHL391".
features: Feature labels to remove, e.g. ["NFKBRE"]. Nested features inside the target are
    removed automatically; do not list them separately.
method: "auto", "restriction", or "pcr".
"""

PLAN_INSERTION = """\
Design a way to put a sequence into a backbone, and verify the product.

Call this when a construct is missing a part the user wants, and to replace one part with
another -- a replacement is a single insertion with a non-empty `at` span, not a deletion
followed by an insertion, and doing it in one step saves a cloning round.

Give the insert either as `donor` plus `donor_features` (taken from another plasmid --
preferred, since it can often be subcloned with no PCR) or as `sequence` (ordered
synthetically). "auto" tries reusing restriction sites before falling back to Gibson.

backbone: Plasmid to insert into, e.g. "pCLM3".
at: Where it goes -- a feature label to replace ("miniCMV"), a position ("1491"), or a span to
    replace ("1491:1512").
sequence: The insert as DNA, when it is not already in the library.
donor: Plasmid to take the insert from, e.g. "pCLM1".
donor_features: Which feature(s) of the donor to take, e.g. ["Lambda BoxB x8"].
method: "auto", "restriction", or "gibson".
"""

SOURCE_SEQUENCE = """\
Decide how to obtain a sequence: from the library, by synthesis, or by ordering.

Call this before planning an insertion whose insert may not be in the library. Leave `sequence`
empty when the user has not given it -- the tool returns a question to put to them rather than
guessing. Provenance is checked before synthesisability, so a repetitive part that already
exists somewhere is moved rather than refused.

name: What the part is called, for the user's benefit, e.g. "barcode".
sequence: The DNA, if known. Leave empty if not.
"""

COMPARE_PRODUCT = """\
Check a designed product against an existing plasmid, ignoring rotation and strand.

Use this when the user has a reference construct the design should match, or to show how a
proposed route differs from something already on the shelf. It is also how you confirm that a
multi-step route actually lands on the target, rather than asserting that it does.

product_id: A handle returned by a plan tool, e.g. "prod_1".
target: Plasmid to compare against, e.g. "pCLM1".
"""

EXPORT_PRODUCT = """\
Write a designed product to a GenBank file that SnapGene can open.

Only call this once the user has picked a route and asked for the file -- it writes to disk, so
it is not a way to inspect a product. Use compare_product for that. The file carries the full
annotation set with coordinates in the parent plasmid's frame, so it can be diffed against the
parent directly.

product_id: A handle returned by a plan tool, e.g. "prod_1".
path: Where to write, e.g. "pCLM1_designed.gb".
"""


SAVE_PROTOCOL = """\
Write a product's bench protocol to a text file.

Call this when the user asks for the protocol saved, written down, or handed over as a file --
it is the only way you can produce that file. You never have the protocol text yourself: it is
rendered from the verified plan and appended below your message after you have finished
writing, so there is nothing for you to copy into a generic file-writing tool. This tool takes
the protocol straight from the plan, so what lands on disk is exact.

This is the plain-text one: the steps, nothing else, for pasting into notes or a message. When
the user is actually going to build the construct, generate_report is the better answer -- it
carries the same protocol plus the reaction volumes, gels and maps they will need at the bench.

product_id: A handle returned by a plan tool, e.g. "prod_1".
path: Where to write, e.g. "protocol.md".
"""


ASK_USER = """\
Put a question to the user and wait for their answer.

This is the only way you can reach the user before your message ends, and it does not return
until they have replied -- there is no time limit, so a question you actually need answered is
never lost. Call this when the answer changes what gets built and you cannot responsibly guess:
which backbone to start from, whether a feature the user did not mention must be removed,
whether a second cloning round is acceptable, what a stock read on the Nanodrop.

Ask everything you need in one call. Do not use it for anything you can decide yourself, or to
report progress, or to seek approval for writing a file -- that is confirmed separately. A
question you can leave to the end of your message belongs there instead, in prose; the user
answers that at the prompt.

questions: The questions to ask, in order. Each is an object:
    question: The question itself, one sentence.
    header: Two or three words naming the choice, e.g. "Backbone". Optional.
    options: The choices, each {"label": ..., "description": ...}. Leave it out for an open
        question. The user can always answer in their own words instead of picking, and can
        skip, so do not add "other" or "no preference" options yourself.
An empty answer means that question was skipped: say what you are assuming and carry on.
"""


#: Returned when nothing is attached to the terminal -- a library caller driving the tools
#: directly, or a test. The model is told the answer is unavailable rather than left to read an
#: empty result as agreement.
NO_ONE_TO_ASK = "there is no interactive user in this session; decide without an answer"


REPORT_INPUTS = """\
List the DNA stocks whose concentration the report's reaction volumes need.

Call this **before** generate_report. Every volume in the digest and ligation tables is a mass
divided by a ng/uL reading, so ask the user for the ones listed here -- `ask_user` gets you the
numbers before you write the report, rather than a turn later -- and pass them to
generate_report. If they do not know them yet, say so and generate the report anyway -- those
tables become input boxes they can fill in at the bench.

product_id: A handle returned by a plan tool, e.g. "prod_1".
"""


GENERATE_REPORT = """\
Write a bench-ready cloning report for a designed product: procedure, reagent tables and
plasmid maps, as one self-contained HTML file.

Call this when the user is going to actually build the construct, or asks for a protocol, a
write-up, a notebook entry or something to send to a colleague. It writes to disk, so only call
it once a route has been chosen. Everything in it is generated from the verified plan -- you do
not need to, and should not, draft the procedure yourself. Prefer this over save_protocol
whenever the destination is the bench rather than a chat message.

Check report_inputs first and ask the user for the stock concentrations; passing them turns the
reaction tables from a form into a protocol. Anything you leave out stays an input box in the
report, which is fine -- never invent a concentration.

product_id: A handle returned by a plan tool, e.g. "prod_1".
path: Where to write, e.g. "pCLM1_report.html".
aim: One sentence on why this construct is being made, in the user's own terms. Leave empty if
    they have not said. Do not restate the cloning strategy here.
name: What to call the construct in the report, e.g. "pCLM1_designed". Leave empty to use the
    name the plan gave it.
concentrations: Measured stocks in ng/uL, keyed by the plasmid names that report_inputs
    returned, e.g. {"pHL391_pcDNA3.1_NFKBRE1-...": 1364}. Omit any the user has not measured.
"""


def build_tools(session: DesignSession, ask=None) -> list:
    """Tool objects bound to one session.

    ``ask`` is the blocking callable that puts a question on the terminal and waits -- see
    :func:`bbl.llm.chat.ask_questions`. Left as ``None``, ``ask_user`` reports that there is
    nobody to ask, which is what a library caller or a test wants.
    """
    from claude_agent_sdk import tool

    def ok(payload) -> dict:
        return {"content": [{"type": "text", "text": json.dumps(payload, default=str)}]}

    @tool("search_inventory", SEARCH_INVENTORY, {"query": str})
    async def search_inventory(args):
        return ok(session.search_inventory(args["query"]))

    @tool("inspect_plasmid", INSPECT_PLASMID, {"plasmid": str})
    async def inspect_plasmid(args):
        return ok(session.inspect_plasmid(args["plasmid"]))

    @tool("plan_deletion", PLAN_DELETION, {"plasmid": str, "features": list, "method": str})
    async def plan_deletion(args):
        return ok(
            session.plan_deletion(
                args["plasmid"], args["features"], args.get("method") or "auto"
            )
        )

    @tool(
        "plan_insertion",
        PLAN_INSERTION,
        {
            "backbone": str,
            "at": str,
            "sequence": str,
            "donor": str,
            "donor_features": list,
            "method": str,
        },
    )
    async def plan_insertion(args):
        return ok(
            session.plan_insertion(
                args["backbone"],
                args["at"],
                sequence=args.get("sequence") or None,
                donor=args.get("donor") or None,
                donor_features=args.get("donor_features") or None,
                method=args.get("method") or "auto",
            )
        )

    @tool("source_sequence", SOURCE_SEQUENCE, {"name": str, "sequence": str})
    async def source_sequence(args):
        return ok(session.source_sequence(args["name"], args.get("sequence") or None))

    @tool("compare_product", COMPARE_PRODUCT, {"product_id": str, "target": str})
    async def compare_product(args):
        return ok(session.compare_product(args["product_id"], args["target"]))

    @tool("export_product", EXPORT_PRODUCT, {"product_id": str, "path": str})
    async def export_product(args):
        return ok(session.export_product(args["product_id"], args["path"]))

    @tool("save_protocol", SAVE_PROTOCOL, {"product_id": str, "path": str})
    async def save_protocol(args):
        return ok(session.save_protocol(args["product_id"], args["path"]))

    @tool("ask_user", ASK_USER, {"questions": list})
    async def ask_user(args):
        questions = args.get("questions") or []
        if ask is None:
            return ok({"answers": [], "note": NO_ONE_TO_ASK})
        # On a worker thread, for the same reason the confirmation gate uses one: the SDK is
        # streaming this turn on our event loop and the user may take minutes to answer. The
        # CLI puts no practical timeout on an in-process MCP call (it is ~28 h, and the idle
        # timeout is off for "sdk" transports), so the wait itself is safe.
        return ok({"answers": await asyncio.to_thread(ask, questions)})

    @tool("report_inputs", REPORT_INPUTS, {"product_id": str})
    async def report_inputs(args):
        return ok(session.dna_needing_concentration(args["product_id"]))

    @tool(
        "generate_report",
        GENERATE_REPORT,
        {
            "product_id": str,
            "path": str,
            "aim": str,
            "name": str,
            "concentrations": dict,
        },
    )
    async def generate_report(args):
        return ok(
            session.generate_report(
                args["product_id"],
                args["path"],
                aim=args.get("aim") or None,
                name=args.get("name") or None,
                concentrations=args.get("concentrations") or None,
            )
        )

    return [
        search_inventory,
        inspect_plasmid,
        plan_deletion,
        plan_insertion,
        source_sequence,
        compare_product,
        export_product,
        save_protocol,
        ask_user,
        report_inputs,
        generate_report,
    ]


def build_server(session: DesignSession, ask=None):
    """The in-process MCP server carrying this session's tools."""
    from claude_agent_sdk import create_sdk_mcp_server

    return create_sdk_mcp_server(
        name=SERVER, version="0.1.0", tools=build_tools(session, ask=ask)
    )


def qualified_names(session: DesignSession) -> list[str]:
    """Tool names as ``allowed_tools`` sees them: ``mcp__bbl__search_inventory``.

    Derived from the tool objects rather than a hand-kept list, so a tool cannot be added
    without becoming callable -- or removed and left dangling in the allowlist.
    """
    return [f"mcp__{SERVER}__{t.name}" for t in build_tools(session)]


def unqualify(name: str) -> str:
    """``mcp__bbl__export_product`` -> ``export_product``; other tool names pass through."""
    prefix = f"mcp__{SERVER}__"
    return name[len(prefix) :] if name.startswith(prefix) else name
