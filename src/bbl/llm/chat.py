"""Interactive design session.

    python -m bbl.llm                            # asks where your plasmids are
    python -m bbl.llm --plasmids ~/constructs    # a directory
    python -m bbl.llm --plasmids a.dna b.dna     # specific files
    python -m bbl.llm --effort xhigh             # harder problems

The library can live anywhere -- a local folder, a mounted share, a scratch directory on a
cluster. If none is given, the session asks. Remote stores (Benchling and friends) plug in as
a ``bbl.sources.PlasmidSource``; see that module.

State persists across turns: products designed early stay addressable by handle, so you can
plan several routes, compare them, and export the one you settle on. Tool calls are printed as
they happen, and anything that writes to disk stops for confirmation.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from .session import DesignSession

BANNER = """\
bbl design session -- {n} plasmids from {dir}

  Describe the construct you want, or ask about what's in the library.
  Commands:  /products    list designed products
             /protocol N  show the protocol for prod_N
             /report N    write a bench-ready HTML report for prod_N
             /library P   switch to a different plasmid library
             /reset       clear the conversation (products are kept)
             /quit
"""

DIM, BOLD, RESET = "\033[2m", "\033[1m", "\033[0m"


def _fmt_args(args: dict) -> str:
    parts = []
    for key, value in args.items():
        if value in (None, "", [], {}):
            continue
        text = str(value)
        parts.append(f"{key}={text[:48] + '…' if len(text) > 48 else text}")
    return ", ".join(parts)


def _ask_concentrations(session, handle) -> dict:
    """Ask for the Nanodrop readings the reaction volumes need. Blank means 'not yet'.

    Worth interrupting for: two numbers turn every volume in the report from a blank into a
    figure. Skipping is free -- the tables stay live in the browser either way.
    """
    needed = session.dna_needing_concentration(handle)
    if not needed:
        return {}
    print(f"{DIM}  stock concentrations for the reaction volumes (Enter to skip){RESET}")
    supplied = {}
    for item in needed:
        try:
            answer = input(f"  {item['plasmid']} ({item['role']}) ng/µL: ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return supplied
        if answer:
            supplied[item["plasmid"]] = answer
    return supplied


def _ask(prompt: str) -> bool:
    try:
        return input(f"{BOLD}  ⚠ {prompt} [y/N] {RESET}").strip().lower() in ("y", "yes")
    except (EOFError, KeyboardInterrupt):
        print()
        return False


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
            answer = input(f"{BOLD}  library › {RESET}").strip()
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
            print(f"{DIM}  {exc}{RESET}")
            continue
        found = source.ids()
        if not found:
            print(f"{DIM}  no plasmid files there -- is the extension .dna or .gb?{RESET}")
            continue
        print(f"{DIM}  found {len(found)} plasmid(s){RESET}\n")
        return source


def run(library=None, effort: str = "high", profile: str | None = None) -> int:
    from .credentials import NoCredentials, resolve, verify, with_refresh
    from .tools import design  # imported here so --help works without the SDK

    # An active Claude session first, an API key as fallback, and an offer to log in if
    # neither. Resolved by building a client rather than guessing at env vars -- a working
    # `ant auth login` profile sets nothing in the environment.
    try:
        credential = resolve(profile, ask=_ask)
        using = verify(credential)
    except NoCredentials as exc:
        print(exc, file=sys.stderr)
        return 1
    print(f"{DIM}auth: {using}{RESET}")

    source = find_library(library)
    if source is None:
        return 1

    session = DesignSession(source, confirm=_ask)
    messages: list[dict] = []
    print(BANNER.format(n=len(session.entries), dir=session.source.name))

    while True:
        try:
            line = input(f"{BOLD}› {RESET}").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return 0
        if not line:
            continue

        if line.startswith("/"):
            command, _, argument = line[1:].partition(" ")
            if command in ("quit", "exit", "q"):
                return 0
            if command == "reset":
                messages = []
                print(f"{DIM}  conversation cleared; {len(session.products)} product(s) kept{RESET}")
            elif command == "products":
                if not session.products:
                    print(f"{DIM}  none yet{RESET}")
                for handle, product in session.products.items():
                    print(f"  {handle}  {len(product.record)} bp  ({product.origin})")
            elif command in ("protocol", "report"):
                argument, _, destination = argument.strip().partition(" ")
                handle = argument.strip() or (list(session.products) or [""])[-1]
                handle = handle if handle.startswith("prod_") else f"prod_{handle}"
                if handle not in session.products:
                    print(f"{DIM}  no such product: {handle}{RESET}")
                elif command == "protocol":
                    print(session.protocol_for(handle))
                else:
                    path = destination.strip() or f"{handle}_report.html"
                    result = session.generate_report(
                        handle, path, concentrations=_ask_concentrations(session, handle)
                    )
                    if result.get("declined"):
                        print(f"{DIM}  not written{RESET}")
                    else:
                        print(
                            f"{DIM}  wrote {result['path']} — {result['sections']} steps, "
                            f"{result['figures']} figures{RESET}"
                        )
                        if not result["maps_included"]:
                            print(f'{DIM}  (no plasmid maps: pip install "bbl[report]"){RESET}')
            elif command == "library":
                if not argument.strip():
                    print(f"  {session.source.name}  ({len(session.entries)} plasmids)")
                else:
                    replacement = find_library(argument.strip())
                    if replacement is not None:
                        session = DesignSession(replacement, confirm=_ask)
                        messages = []
                        print(
                            f"{DIM}  switched to {session.source.name} "
                            f"({len(session.entries)} plasmids); conversation cleared{RESET}"
                        )
            else:
                print(f"{DIM}  unknown command{RESET}")
            continue

        messages.append({"role": "user", "content": line})
        try:
            show = lambda name, args: print(  # noqa: E731
                f"{DIM}  → {name}({_fmt_args(args)}){RESET}", flush=True
            )
            # Session tokens expire; re-mint once and retry rather than dying mid-conversation.
            reply, messages = with_refresh(
                lambda client: design(
                    session, messages, client=client, effort=effort, on_tool=show
                ),
                credential,
                profile,
            )
        except KeyboardInterrupt:
            print(f"\n{DIM}  interrupted{RESET}")
            messages.pop()
            continue
        except Exception as exc:  # keep the session alive on a transient API failure
            print(f"{DIM}  error: {type(exc).__name__}: {exc}{RESET}", file=sys.stderr)
            messages.pop()
            continue
        print(f"\n{reply}\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="bbl.llm", description=__doc__)
    parser.add_argument(
        "--plasmids",
        nargs="*",
        default=None,
        help="directory of .dna/.gb files, or a list of files. If omitted, uses "
        "$BBL_PLASMID_DIR, then ./plasmid, then asks.",
    )
    parser.add_argument(
        "--profile",
        default=None,
        help="ant auth profile to mint a session token from (default: the active profile)",
    )
    parser.add_argument(
        "--effort",
        default="high",
        choices=["low", "medium", "high", "xhigh", "max"],
        help="reasoning effort (default: high)",
    )
    args = parser.parse_args(argv)
    given = args.plasmids
    if given is not None:
        given = given[0] if len(given) == 1 else given
    return run(given, args.effort, args.profile)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
