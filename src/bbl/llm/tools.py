"""Tool definitions and the conversation loop.

The wrappers are thin: they serialise a :class:`DesignSession` method to JSON and carry a
docstring written for the model. Descriptions say *when* to call, not just what the tool does --
that is what actually drives tool choice.

Importing this module requires the ``anthropic`` SDK; ``session.py`` does not.
"""

from __future__ import annotations

import json

from .prompts import build_system_blocks
from .session import DesignSession

MODEL = "claude-opus-5"
MAX_TOKENS = 16000

#: Tools that change something outside the process. Gated before execution.
OUTWARD_FACING = {"export_product", "generate_report"}


def build_tools(session: DesignSession) -> list:
    """Tool objects bound to one session."""
    from anthropic import beta_tool

    def dump(payload) -> str:
        return json.dumps(payload, default=str)

    @beta_tool
    def search_inventory(query: str) -> str:
        """Search the plasmid library by name, feature label, or exact DNA sequence.

        Call this first for any request that names a plasmid or a part, and whenever you need
        to know whether a sequence already exists somewhere in the library -- a part that is
        already on a plasmid can be moved rather than synthesised, which sidesteps the
        synthesis complexity limit entirely. A DNA string of 12 bp or more is searched as a
        sequence on both strands.

        Args:
            query: A plasmid name ("pHL391"), a feature label ("BoxB"), or a DNA sequence.
        """
        return dump(session.search_inventory(query))

    @beta_tool
    def inspect_plasmid(plasmid: str) -> str:
        """List a plasmid's length, functional features with coordinates, and unique cutters.

        Call this before planning an edit, to see what is present and where. Feature spans
        from here are what you pass to plan_deletion and plan_insertion.

        Args:
            plasmid: Plasmid name, e.g. "pHL391". A prefix is enough if it is unambiguous.
        """
        return dump(session.inspect_plasmid(plasmid))

    @beta_tool
    def plan_deletion(plasmid: str, features: list[str], method: str = "auto") -> str:
        """Design a way to delete features from a plasmid, and verify the product.

        Call this whenever the user wants something removed from a construct -- an unwanted
        response element, a promoter being swapped out, a spacer. Also call it when a candidate
        backbone is right except for extra features: deleting them is usually far cheaper than
        building the construct from parts.

        Use "auto" (the default) unless the user asked for a specific chemistry: it tries
        restriction first and falls back to inverse PCR, which is the preference order you
        should follow anyway. A result with feasible=false carries `reasons` explaining what
        was rejected -- read it and consider the other method rather than calling again with
        the same arguments.

        Args:
            plasmid: Plasmid to edit, e.g. "pHL391".
            features: Feature labels to remove, e.g. ["NFKBRE"]. Nested features inside the
                target are removed automatically; do not list them separately.
            method: "auto", "restriction", or "pcr".
        """
        return dump(session.plan_deletion(plasmid, features, method))

    @beta_tool
    def plan_insertion(
        backbone: str,
        at: str,
        sequence: str = "",
        donor: str = "",
        donor_features: list[str] | None = None,
        method: str = "auto",
    ) -> str:
        """Design a way to put a sequence into a backbone, and verify the product.

        Call this when a construct is missing a part the user wants, and to replace one part
        with another -- a replacement is a single insertion with a non-empty `at` span, not a
        deletion followed by an insertion, and doing it in one step saves a cloning round.

        Give the insert either as `donor` plus `donor_features` (taken from another plasmid --
        preferred, since it can often be subcloned with no PCR) or as `sequence` (ordered
        synthetically). "auto" tries reusing restriction sites before falling back to Gibson.

        Args:
            backbone: Plasmid to insert into, e.g. "pCLM3".
            at: Where it goes -- a feature label to replace ("miniCMV"), a position ("1491"),
                or a span to replace ("1491:1512").
            sequence: The insert as DNA, when it is not already in the library.
            donor: Plasmid to take the insert from, e.g. "pCLM1".
            donor_features: Which feature(s) of the donor to take, e.g. ["Lambda BoxB x8"].
            method: "auto", "restriction", or "gibson".
        """
        return dump(
            session.plan_insertion(
                backbone,
                at,
                sequence=sequence or None,
                donor=donor or None,
                donor_features=donor_features,
                method=method,
            )
        )

    @beta_tool
    def source_sequence(name: str, sequence: str = "") -> str:
        """Decide how to obtain a sequence: from the library, by synthesis, or by ordering.

        Call this before planning an insertion whose insert may not be in the library. Leave
        `sequence` empty when the user has not given it -- the tool returns a question to put
        to them rather than guessing. Provenance is checked before synthesisability, so a
        repetitive part that already exists somewhere is moved rather than refused.

        Args:
            name: What the part is called, for the user's benefit, e.g. "barcode".
            sequence: The DNA, if known. Leave empty if not.
        """
        return dump(session.source_sequence(name, sequence or None))

    @beta_tool
    def compare_product(product_id: str, target: str) -> str:
        """Check a designed product against an existing plasmid, ignoring rotation and strand.

        Use this when the user has a reference construct the design should match, or to show
        how a proposed route differs from something already on the shelf.

        Args:
            product_id: A handle returned by a plan tool, e.g. "prod_1".
            target: Plasmid to compare against, e.g. "pCLM1".
        """
        return dump(session.compare_product(product_id, target))

    @beta_tool
    def export_product(product_id: str, path: str) -> str:
        """Write a designed product to a GenBank file that SnapGene can open.

        Only call this once the user has picked a route and asked for the file -- it writes to
        disk, so it is not a way to inspect a product. Use compare_product for that. The file
        carries the full annotation set with coordinates in the parent plasmid's frame, so it
        can be diffed against the parent directly.

        Args:
            product_id: A handle returned by a plan tool, e.g. "prod_1".
            path: Where to write, e.g. "pCLM1_designed.gb".
        """
        return dump(session.export_product(product_id, path))

    @beta_tool
    def report_inputs(product_id: str) -> str:
        """List the DNA stocks whose concentration the report's reaction volumes need.

        Call this **before** generate_report. Every volume in the digest and ligation tables
        is a mass divided by a ng/µL reading, so ask the user for the ones listed here and
        pass them to generate_report. If they do not know them yet, say so and generate the
        report anyway — the tables become input boxes they can fill in at the bench.

        Args:
            product_id: A handle returned by a plan tool, e.g. "prod_1".
        """
        return dump(session.dna_needing_concentration(product_id))

    @beta_tool
    def generate_report(
        product_id: str,
        path: str,
        aim: str = "",
        name: str = "",
        concentrations: dict | None = None,
    ) -> str:
        """Write a bench-ready cloning report for a designed product: procedure, reagent
        tables and plasmid maps, as one self-contained HTML file.

        Call this when the user is going to actually build the construct, or asks for a
        protocol, a write-up, a notebook entry or something to send to a colleague. It writes
        to disk, so only call it once a route has been chosen. Everything in it is generated
        from the verified plan — you do not need to, and should not, draft the procedure
        yourself.

        Check report_inputs first and ask the user for the stock concentrations; passing them
        turns the reaction tables from a form into a protocol. Anything you leave out stays an
        input box in the report, which is fine — never invent a concentration.

        Args:
            product_id: A handle returned by a plan tool, e.g. "prod_1".
            path: Where to write, e.g. "pCLM1_report.html".
            aim: One sentence on why this construct is being made, in the user's own terms.
                Leave empty if they have not said. Do not restate the cloning strategy here.
            name: What to call the construct in the report, e.g. "pCLM1_designed". Leave empty
                to use the name the plan gave it.
            concentrations: Measured stocks in ng/µL, keyed by the plasmid names that
                report_inputs returned, e.g. {"pHL391_pcDNA3.1_NFKBRE1-...": 1364}. Omit any
                the user has not measured.
        """
        return dump(
            session.generate_report(
                product_id,
                path,
                aim=aim or None,
                name=name or None,
                concentrations=concentrations,
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
        report_inputs,
        generate_report,
    ]


def design(session, messages, client=None, effort: str = "high", on_tool=None):
    """Run one turn of the design conversation.

    ``on_tool(name, args)`` is called as each tool call is issued, so a UI can show work in
    progress. Gating happens inside the tool wrappers via ``session.confirm``.

    Returns ``(reply_text, messages)`` with the protocol for any product the model settled on
    appended verbatim -- the harness renders it, not the model, so enzyme names and base-pair
    counts cannot be paraphrased into something wrong.
    """
    import anthropic

    client = client or anthropic.Anthropic()
    before = set(session.products)

    runner = client.beta.messages.tool_runner(
        model=MODEL,
        max_tokens=MAX_TOKENS,
        thinking={"type": "adaptive"},
        output_config={"effort": effort},
        system=build_system_blocks(session.entries, session.config),
        tools=build_tools(session),
        messages=messages,
    )

    final = None
    for message in runner:
        final = message
        if on_tool is not None:
            for block in message.content:
                if block.type == "tool_use":
                    on_tool(block.name, block.input)
        messages.append({"role": "assistant", "content": message.content})
        response = runner.generate_tool_call_response()
        if response is not None:
            messages.append(response)

    text = "\n".join(b.text for b in (final.content if final else []) if b.type == "text")

    fresh = [h for h in session.products if h not in before]
    for handle in fresh:
        if handle in text:
            text += f"\n\n--- protocol for {handle} ---\n{session.protocol_for(handle)}"
    return text, messages
