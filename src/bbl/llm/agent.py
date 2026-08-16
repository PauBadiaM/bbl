"""The agent: options, the confirmation gate, and one turn of the design conversation.

Conversation state lives in the ``claude`` subprocess, not in a Python list. One
:class:`~claude_agent_sdk.ClaudeSDKClient` is held open for the whole REPL, so a turn that
fails partway cannot leave a half-written transcript behind -- there is no transcript here to
corrupt. That is the bug this module exists to not have.

See docs/DECISIONS.md D62 for why bbl talks to the CLI rather than the Anthropic SDK directly.
"""

from __future__ import annotations

import asyncio
import os
from dataclasses import dataclass
from pathlib import Path

from .prompts import build_system_prompt
from .session import DesignSession
from .tools import OUTWARD_FACING, build_server, qualified_names, unqualify

DEFAULT_MODEL = os.environ.get("BBL_MODEL") or "claude-opus-5"

#: What the model is told when the user refuses a write. Matches the wording
#: ``session.export_product`` uses for the same refusal, so the model sees one behaviour
#: whichever layer catches it (docs/DECISIONS.md D69).
DECLINED = "the user declined; do not retry, ask what they want instead"


@dataclass
class Turn:
    """What one exchange produced."""

    text: str
    session_id: str | None = None
    is_error: bool = False
    error: str | None = None


def build_options(
    session: DesignSession,
    *,
    env: dict[str, str],
    cwd: Path,
    model: str = DEFAULT_MODEL,
    effort: str = "high",
    can_use_tool=None,
    resume: str | None = None,
):
    """Options for a bbl design session.

    What ``allowed_tools`` *omits* is the load-bearing part. An entry that allows a whole tool
    auto-approves it before ``can_use_tool`` is ever consulted, so the outward-facing tools and
    ``Write`` are deliberately left off: they fall through to the gate. The read-only design
    tools and ``Read`` are allowlisted outright, because stopping to confirm an inventory
    search would make the session unusable.

    ``Write`` is *not* scoped with a path rule here. Measured against CLI 2.1.x and
    ``claude-agent-sdk`` 0.2.139: a path-scoped entry in ``allowed_tools``
    (``Write(//abs/path/**)``, the form acumen uses, and every other spelling tried) is
    silently ignored -- only a whole-tool pattern such as ``Write(*)`` takes effect there.
    Path scoping works through a settings file's ``permissions.allow``, which is how acumen
    does it. Rather than carry a temp settings file, bbl leaves ``Write`` off the allowlist and
    enforces the working-directory scope in :func:`confirmation_gate`, where it is an ordinary
    path comparison we can test.

    Three settings keep the session self-contained, and each closes a different door.
    ``setting_sources=[]`` is set explicitly rather than left to default, because the SDK's
    default of ``None`` loads the user's own settings *and* memories. ``strict_mcp_config``
    drops every MCP server the CLI would otherwise pick up from project, user or plugin
    config -- without it the agent sees whatever connectors the user has on their account and
    starts reasoning about why they are unauthorized, mid-cloning-design.
    ``CLAUDE_CODE_DISABLE_CLAUDE_MDS`` is the third, set in :func:`bbl.llm.auth.agent_env`,
    because project discovery walks *up* from cwd.
    """
    from claude_agent_sdk import ClaudeAgentOptions

    unattended = [n for n in qualified_names(session) if unqualify(n) not in OUTWARD_FACING]
    return ClaudeAgentOptions(
        cwd=str(cwd),
        env=env,
        model=model,
        system_prompt=build_system_prompt(session.entries, session.config),
        mcp_servers={"bbl": build_server(session)},
        strict_mcp_config=True,
        allowed_tools=[*unattended, "Read"],
        setting_sources=[],
        permission_mode="default",
        can_use_tool=can_use_tool,
        resume=resume,
        effort=effort,
    )


#: Reasoning effort levels the CLI accepts, in order.
EFFORTS = ("low", "medium", "high", "xhigh", "max")


def _inside(target: Path, root: Path) -> bool:
    """Whether ``target`` lands within ``root``, symlinks and ``..`` resolved.

    ``strict=False`` matters: the file being written does not exist yet, so the check is on
    where it *would* land. Anything unresolvable is treated as outside and therefore asked
    about -- the safe direction for a question the user can always answer yes to.
    """
    try:
        return target.resolve(strict=False).is_relative_to(root.resolve())
    except (OSError, ValueError):
        return False


def confirmation_gate(confirm, cwd: Path):
    """A ``can_use_tool`` callback that stops before anything leaves the process.

    Two things reach this gate: the design tools in :data:`~bbl.llm.tools.OUTWARD_FACING`, and
    every ``Write``. Exports always ask -- a GenBank file is the deliverable, and which one the
    user meant is worth a keystroke. A ``Write`` inside the working directory is allowed
    silently, because saving a protocol next to the plasmids it was designed from is the
    expected thing; one that escapes the directory asks. See :func:`build_options` for why the
    scope is checked here rather than expressed as an ``allowed_tools`` path rule.

    ``confirm`` is a blocking ``input()``-style callable, so it is run on a worker thread: the
    SDK is streaming the model's output on this event loop, and blocking it would stall the
    session while the user reads the prompt.
    """
    from claude_agent_sdk import PermissionResultAllow, PermissionResultDeny

    async def gate(tool_name: str, tool_input: dict, context):
        bare = unqualify(tool_name)
        if bare in OUTWARD_FACING:
            # The resolved path, not the bare argument: the model gives a relative one, and
            # which directory it lands in is the part worth confirming.
            question = f"save {tool_input.get('product_id', 'the product')} to "
            question += f"{Path(cwd, tool_input.get('path', '?'))}?"
        elif bare == "Write":
            target = Path(cwd, tool_input.get("file_path", ""))
            if _inside(target, cwd):
                return PermissionResultAllow()
            question = f"write outside the working directory: {target}?"
        else:
            return PermissionResultAllow()

        if confirm is None or await asyncio.to_thread(confirm, question):
            return PermissionResultAllow()
        return PermissionResultDeny(message=DECLINED)

    return gate


def new_protocols(session: DesignSession, before: set[str], text: str) -> str:
    """Protocols for products the model settled on this turn, rendered by the harness.

    The model names a handle; the bench text is generated from the verified plan and appended
    verbatim. It is never written by the model, so enzyme names and base-pair counts cannot be
    paraphrased into something that does not work at the bench (docs/DECISIONS.md D64).

    Transport-independent by design -- this is the invariant, not the plumbing around it.
    """
    made = [handle for handle in session.products if handle not in before]
    named = [handle for handle in made if handle in text]
    # Prefer the ones the model actually pointed at -- when it plans three routes and
    # recommends one, only that one's protocol belongs under the message. But a good answer can
    # describe the route without ever typing "prod_1", and requiring the literal handle then
    # silently drops the protocol from the very turn that earned it. So: named if any, else all.
    return "".join(
        f"\n\n--- protocol for {handle} ---\n{session.protocol_for(handle)}"
        for handle in (named or made)
    )


async def run_turn(client, prompt: str, renderer) -> Turn:
    """Send one message and render the reply as it streams back.

    The trailing-exception case is load-bearing. The CLI exits non-zero on purpose after
    reporting an ``is_error`` result -- a turn cap, a quota breach -- and the SDK surfaces that
    as an exception *after* the structured result has already arrived. Raising would throw away
    the turns, cost and session id of an exchange that really happened, and would lose the
    actionable message, which is often only on stderr.
    """
    result = None
    await client.query(prompt)
    try:
        async for message in client.receive_response():
            renderer.handle(message)
            if getattr(message, "num_turns", None) is not None:
                result = message
    except Exception as exc:
        if result is None:
            raise
        return Turn(
            text=renderer.text,
            session_id=getattr(result, "session_id", None),
            is_error=True,
            error=str(exc),
        )

    return Turn(
        text=renderer.text,
        session_id=getattr(result, "session_id", None),
        is_error=bool(getattr(result, "is_error", False)),
        error="; ".join(getattr(result, "errors", None) or []) or None,
    )
