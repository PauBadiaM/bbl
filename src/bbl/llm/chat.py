"""The interactive design session.

The library can live anywhere -- a local folder, a mounted share, a scratch directory on a
cluster. If none is given, the session asks. Remote stores (Benchling and friends) plug in as a
``bbl.sources.PlasmidSource``; see that module.

State persists across turns in two places. Products designed early stay addressable by handle in
the :class:`~bbl.llm.session.DesignSession`, so you can plan several routes, compare them, and
export the one you settle on. The conversation itself lives in the ``claude`` subprocess and is
written to disk by it, so ``bbl --resume`` picks the design back up where it stopped.
"""

from __future__ import annotations

import asyncio
import contextlib
import json
import os
import signal
import sys
import warnings
from pathlib import Path

from .agent import DEFAULT_MODEL, build_options, confirmation_gate, new_protocols, run_turn
from .render import Renderer, Style
from .session import DesignSession

BANNER = """\
bbl -- {n} plasmids from {dir}

  Commands:  /products    list designed products
             /protocol N  show the protocol for prod_N
             /report N    write a bench-ready HTML report for prod_N
             /library P   switch to a different plasmid library
             /reset       start the conversation over (products are kept)
             /quit
"""

#: Sent automatically as the first turn so the agent opens the conversation rather than
#: waiting at an empty prompt. What it says is specified in the system prompt, not here.
OPENING = "Begin the session."

#: Where the last conversation's id is remembered, so `bbl --resume` needs no argument.
STATE_DIR = ".bbl"
STATE_FILE = "session.json"

STYLE = Style()


async def _ask_concentrations(session, handle) -> dict:
    """Ask for the Nanodrop readings the reaction volumes need. Blank means 'not yet'.

    Worth interrupting for: two numbers turn every volume in the report from a blank into a
    figure. Skipping is free -- the tables stay live in the browser either way.

    Read on a worker thread for the same reason ``converse`` reads the prompt there: this runs
    on the event loop the SDK client is attached to, and a blocking ``input()`` would stall it
    for as long as the user takes to find the number.
    """
    needed = session.dna_needing_concentration(handle)
    if not needed:
        return {}
    _dim("stock concentrations for the reaction volumes (Enter to skip)")
    supplied = {}
    for item in needed:
        try:
            answer = await asyncio.to_thread(
                input, f"  {item['plasmid']} ({item['role']}) ng/µL: "
            )
        except (EOFError, KeyboardInterrupt):
            print()
            return supplied
        if answer.strip():
            supplied[item["plasmid"]] = answer.strip()
    return supplied


def _ask(prompt: str) -> bool:
    try:
        return input(f"{STYLE.bold}  ⚠ {prompt} [y/N] {STYLE.reset}").strip().lower() in (
            "y",
            "yes",
        )
    except (EOFError, KeyboardInterrupt):
        print()
        return False


def save_session_id(session_id: str, cwd: Path | None = None) -> None:
    """Remember the conversation so it can be resumed from this directory."""
    directory = Path(cwd or Path.cwd()) / STATE_DIR
    try:
        directory.mkdir(parents=True, exist_ok=True)
        (directory / STATE_FILE).write_text(json.dumps({"session_id": session_id}))
    except OSError:
        pass  # a read-only working directory is not a reason to lose the turn


def last_session_id(cwd: Path | None = None) -> str | None:
    """The most recent conversation started in this directory, if any."""
    try:
        data = json.loads((Path(cwd or Path.cwd()) / STATE_DIR / STATE_FILE).read_text())
    except (OSError, ValueError):
        return None
    return data.get("session_id") if isinstance(data, dict) else None


def find_library(given=None):
    """Resolve where the plasmids live, asking the user if we cannot tell.

    Order: what was passed, then ``$BBL_PLASMID_DIR``, then a ``plasmid/`` folder in the
    working directory, then ask. Nothing is assumed about the machine this runs on.
    """
    from ..sources import as_source

    candidates = [given, os.environ.get("BBL_PLASMID_DIR"), Path.cwd() / "plasmid"]
    for candidate in candidates:
        if not candidate:
            continue
        try:
            source = as_source(candidate)
        except (FileNotFoundError, OSError):
            continue
        if source.ids():
            return source

    print(
        "Point me at your plasmids: a directory of .dna/.gb files, or a list of files "
        "separated by spaces.\n"
    )
    while True:
        try:
            answer = input(f"{STYLE.bold}  library › {STYLE.reset}").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return None
        if not answer:
            continue
        if answer in ("/quit", "/q", "quit", "exit"):
            return None
        parts = [p.strip("\"'") for p in answer.split()]
        try:
            source = as_source(parts[0] if len(parts) == 1 else parts)
        except (FileNotFoundError, OSError) as exc:
            print(f"{STYLE.dim}  {exc}{STYLE.reset}")
            continue
        found = source.ids()
        if not found:
            print(f"{STYLE.dim}  no plasmid files there -- is the extension .dna or .gb?{STYLE.reset}")
            continue
        print(f"{STYLE.dim}  found {len(found)} plasmid(s){STYLE.reset}\n")
        return source


def _dim(text: str) -> None:
    print(f"{STYLE.dim}  {text}{STYLE.reset}")


class _Interrupted(Exception):
    """Ctrl-C during a turn: abandon the turn, keep the session."""


async def _interruptible_turn(client, prompt: str, renderer):
    """Run one turn with Ctrl-C bound to "abort this turn", not "kill bbl".

    Python's default SIGINT raises ``KeyboardInterrupt`` in the main thread, which inside an
    ``await`` surfaces as a ``CancelledError`` that no ordinary ``except`` catches and takes
    the whole process down -- a long design turn would be unabortable except by losing the
    conversation. Installing a loop-level handler for the duration of the turn replaces that
    with a plain task cancellation we can act on: tell the agent to stop, then return to the
    prompt with the session still connected.

    The handler is removed afterwards, so Ctrl-C at an idle prompt still exits, which is what
    it should do when there is nothing to abort.
    """
    loop = asyncio.get_running_loop()
    task = asyncio.ensure_future(run_turn(client, prompt, renderer))

    try:
        loop.add_signal_handler(signal.SIGINT, task.cancel)
    except (NotImplementedError, RuntimeError):  # not a Unix main-thread loop
        return await task

    try:
        return await task
    except asyncio.CancelledError:
        # Best-effort: the CLI shares our process group, so a terminal Ctrl-C may already have
        # reached it. Either way the turn is over as far as this session is concerned.
        with contextlib.suppress(Exception):
            await client.interrupt()
        raise _Interrupted from None
    finally:
        with contextlib.suppress(NotImplementedError, RuntimeError):
            loop.remove_signal_handler(signal.SIGINT)


def _resolve_handle(argument: str, session: DesignSession) -> str:
    """``"2"`` | ``"prod_2"`` | ``""`` -> a product handle. Empty means the most recent one."""
    handle = argument or (list(session.products) or [""])[-1]
    return handle if handle.startswith("prod_") else f"prod_{handle}"


def report_digest(record: dict) -> str:
    """What a written report looks like in the terminal: the path, then what is in it.

    The step titles, not the protocol body. The protocol is already appended after any turn
    that designs a product (D64), and printing it a second time here would double it on the
    common path where a report is written in the same turn the product was made.
    """
    lines = [
        f"  wrote {record['path']} — {record['sections']} steps, {record['figures']} figures"
    ]
    lines += [f"    {title}" for title in record["steps"]]
    if record["concentrations_missing"]:
        missing = ", ".join(record["concentrations_missing"])
        lines.append(f"    blank until measured: {missing}")
    return "\n".join(lines)


async def _handle_command(line: str, session: DesignSession) -> str | None:
    """Run a slash command. Returns a control word for the caller, or None if handled here."""
    command, _, argument = line[1:].partition(" ")
    argument = argument.strip()

    if command in ("quit", "exit", "q"):
        return "quit"
    if command == "reset":
        return "reset"
    if command == "library":
        if not argument:
            print(f"  {session.source.name}  ({len(session.entries)} plasmids)")
            return None
        return f"library {argument}"
    if command == "products":
        if not session.products:
            _dim("none yet")
        for handle, product in session.products.items():
            print(f"  {handle}  {len(product.record)} bp  ({product.origin})")
        return None
    if command == "protocol":
        handle = _resolve_handle(argument, session)
        if handle in session.products:
            print(session.protocol_for(handle))
        else:
            _dim(f"no such product: {handle}")
        return None
    if command == "report":
        wanted, _, destination = argument.partition(" ")
        handle = _resolve_handle(wanted.strip(), session)
        if handle not in session.products:
            _dim(f"no such product: {handle}")
            return None
        path = destination.strip() or f"{handle}_report.html"
        concentrations = await _ask_concentrations(session, handle)
        result = session.generate_report(handle, path, concentrations=concentrations)
        # No confirmation gate on this path: the user typed the command and the destination,
        # which is the thing the gate exists to establish. The gate covers the model.
        if result.get("declined"):
            _dim("not written")
        else:
            print(report_digest(session.reports[-1]))
        return None

    _dim("unknown command")
    return None


async def converse(
    session: DesignSession,
    *,
    env: dict,
    cwd: Path,
    model: str,
    effort: str,
    resume: str | None,
    debug: bool = False,
    banner: bool = True,
) -> int:
    """Hold one connected agent open and drive it from the prompt."""
    from claude_agent_sdk import CanUseToolShadowedWarning, ClaudeSDKClient, CLIConnectionError

    # The SDK warns that allowlisted tools bypass can_use_tool. Here that is the design: the
    # read-only design tools are meant to run unattended, and only the outward-facing ones are
    # left off the allowlist so they reach the gate. See agent.build_options.
    warnings.filterwarnings("ignore", category=CanUseToolShadowedWarning)

    options = build_options(
        session,
        env=env,
        cwd=cwd,
        model=model,
        effort=effort,
        can_use_tool=confirmation_gate(_ask, cwd),
        resume=resume,
    )

    async with ClaudeSDKClient(options=options) as client:
        if banner:
            print(BANNER.format(n=len(session.entries), dir=session.source.name))
        pending: str | None = None if resume else OPENING
        completed = 0

        while True:
            if pending is None:
                try:
                    line = await asyncio.to_thread(input, f"{STYLE.bold}› {STYLE.reset}")
                except (EOFError, KeyboardInterrupt):
                    print()
                    return 0
                line = line.strip()
                if not line:
                    continue
                if line.startswith("/"):
                    control = await _handle_command(line, session)
                    if control == "quit":
                        return 0
                    if control is not None:
                        raise _Restart(control)
                    continue
                pending = line

            before = set(session.products)
            reports_before = len(session.reports)
            renderer = Renderer(STYLE)
            try:
                turn = await _interruptible_turn(client, pending, renderer)
            except _Interrupted:
                # A terminal Ctrl-C reaches the whole process group, and the SDK spawns the
                # CLI in ours, so the agent process is probably gone. Reconnect and resume
                # rather than leaving a session that looks alive and fails on the next turn.
                _dim("interrupted")
                raise _Reconnect(renderer.session_id, bool(completed)) from None
            except CLIConnectionError:
                _dim("lost the agent connection")
                raise _Reconnect(renderer.session_id, bool(completed)) from None
            except Exception as exc:
                if debug:
                    raise
                print(
                    f"{STYLE.dim}  error: {type(exc).__name__}: {exc}{STYLE.reset}",
                    file=sys.stderr,
                )
                pending = None
                continue
            pending = None
            completed += 1

            if turn.session_id:
                save_session_id(turn.session_id, cwd)
            if turn.is_error and turn.error:
                _dim(f"error: {turn.error}")

            # The harness renders the protocol, never the model (docs/DECISIONS.md D64).
            protocols = new_protocols(session, before, turn.text)
            if protocols:
                print(protocols)

            # Same principle for reports: the model asked for the file, the harness says what
            # landed. Identical output to the /report path, which shares report_digest.
            for record in session.reports[reports_before:]:
                print(report_digest(record))


class _Restart(Exception):
    """Raised to tear the agent down and rebuild it -- /reset and /library both need this."""

    def __init__(self, control: str):
        super().__init__(control)
        self.control = control


class _Reconnect(Exception):
    """The agent process is gone. Rebuild it and resume the same conversation."""

    def __init__(self, session_id: str | None, progressed: bool = False):
        super().__init__(session_id)
        self.session_id = session_id
        #: Whether any turn completed before the connection was lost. A reconnect that makes
        #: no progress is how a broken `claude` install would otherwise spin forever.
        self.progressed = progressed


#: Consecutive reconnects with no completed turn in between before bbl gives up.
MAX_DEAD_RECONNECTS = 3


async def run(
    library=None,
    *,
    env: dict,
    model: str = DEFAULT_MODEL,
    effort: str = "high",
    resume: str | None = None,
    debug: bool = False,
) -> int:
    """Resolve the library, then converse until the user quits.

    ``/reset`` and ``/library`` restart the agent rather than mutating it: the conversation
    lives in the subprocess and the tools are bound to one library, so a clean rebuild is the
    only way to make either of them mean what they say.
    """
    from claude_agent_sdk import CLIConnectionError, CLINotFoundError

    cwd = Path.cwd()
    source = find_library(library)
    if source is None:
        return 1

    session = DesignSession(source)
    banner, dead_reconnects = True, 0
    while True:
        try:
            return await converse(
                session,
                env=env,
                cwd=cwd,
                model=model,
                effort=effort,
                resume=resume,
                debug=debug,
                banner=banner,
            )
        except (CLIConnectionError, CLINotFoundError) as exc:
            # Raised while starting the agent, before any turn -- retrying would just fail the
            # same way, so this is fatal rather than a reconnect.
            print(f"bbl: could not start the agent: {exc}", file=sys.stderr)
            return 1
        except _Reconnect as lost:
            dead_reconnects = 0 if lost.progressed else dead_reconnects + 1
            if dead_reconnects > MAX_DEAD_RECONNECTS:
                print("bbl: cannot keep the agent connected; giving up", file=sys.stderr)
                return 1
            resume = lost.session_id or last_session_id(cwd)
            if resume:
                save_session_id(resume, cwd)
            _dim("reconnecting, resuming where we left off" if resume else "reconnecting")
            banner = False
            continue
        except _Restart as restart:
            resume, banner = None, True
            if restart.control == "reset":
                # A new conversation, but the same session object: products designed before the
                # reset stay addressable by handle, which is the point of keeping them.
                _dim(f"conversation cleared; {len(session.products)} product(s) kept")
                continue
            replacement = find_library(restart.control.removeprefix("library "))
            if replacement is not None:
                source = replacement
                session = DesignSession(source)
                _dim(f"switched to {source.name} ({len(session.entries)} plasmids)")
