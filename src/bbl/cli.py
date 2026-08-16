"""``bbl`` -- start an interactive plasmid design session.

    bbl                                  # asks where your plasmids are
    bbl --plasmids ~/constructs          # a directory
    bbl --plasmids a.dna b.dna           # specific files
    bbl --resume                         # pick up the last design in this directory
    bbl --effort xhigh                   # harder problems

One command, one agent. It opens by asking what plasmid you are trying to build, works the
route against your library, and ends with a protocol for getting there.
"""

from __future__ import annotations

import argparse
import asyncio
import sys


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bbl",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--plasmids",
        nargs="*",
        default=None,
        metavar="P",
        help="directory of .dna/.gb files, or a list of files. If omitted, uses "
        "$BBL_PLASMID_DIR, then ./plasmid, then asks.",
    )
    parser.add_argument(
        "--auth",
        choices=("auto", "session", "api"),
        default="auto",
        help="which credential to bill: 'session' (Claude subscription), 'api' (Anthropic API), "
        "or 'auto' (default: session if you're logged in, else the API)",
    )
    parser.add_argument(
        "--effort",
        default="high",
        choices=("low", "medium", "high", "xhigh", "max"),
        help="reasoning effort (default: high)",
    )
    parser.add_argument(
        "--model",
        default=None,
        help="model to design with (default: claude-opus-5, or $BBL_MODEL)",
    )
    parser.add_argument(
        "--resume",
        nargs="?",
        const="",
        default=None,
        metavar="ID",
        help="continue a previous design. With no argument, the last one in this directory.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="re-raise errors with a full traceback instead of a one-line summary",
    )
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    # Imported here so --help works without the SDK installed.
    from .llm.agent import DEFAULT_MODEL
    from .llm.auth import AuthError, agent_env, check_agent_cli, describe, resolve_auth_mode
    from .llm.chat import STYLE, last_session_id, run

    # Fail before the library scan and the banner: an unauthenticated session would otherwise
    # start, look healthy, and only break on the first turn.
    try:
        check_agent_cli()
        auth_mode = resolve_auth_mode(args.auth)
    except AuthError as exc:
        print(f"bbl: {exc}", file=sys.stderr)
        return 1
    print(f"{STYLE.dim}auth: {describe(auth_mode)}{STYLE.reset}")

    resume = args.resume
    if resume == "":  # bare --resume
        resume = last_session_id()
        if resume is None:
            print("bbl: no previous design in this directory to resume", file=sys.stderr)
            return 1

    plasmids = args.plasmids
    if plasmids is not None:
        plasmids = plasmids[0] if len(plasmids) == 1 else plasmids

    try:
        return asyncio.run(
            run(
                plasmids,
                env=agent_env(auth_mode),
                model=args.model or DEFAULT_MODEL,
                effort=args.effort,
                resume=resume,
                debug=args.debug,
            )
        )
    except KeyboardInterrupt:
        print()
        return 130


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
