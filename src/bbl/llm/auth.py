"""Claude credential resolution for an agent session.

Ported from ``acumen`` (``scverse/acumen``, ``src/acumen/env.py``), minus the Codex provider
and minus the sandboxing: bbl runs the agent in the user's own working directory with their own
environment, so there is no throwaway ``HOME`` to seed credentials into.

Two credential classes, each with a detector. A *session* credential bills the Claude
subscription -- the ``claudeAiOauth`` object that ``claude`` login writes to
``~/.claude/.credentials.json``, or a portable ``CLAUDE_CODE_OAUTH_TOKEN`` from
``claude setup-token``. An *api* credential bills per token: ``ANTHROPIC_API_KEY``,
``ANTHROPIC_AUTH_TOKEN``, or a Bedrock/Vertex routing flag whose own cloud credentials apply.

Only the Agent SDK can use the subscription login, which is why bbl talks to Claude through the
``claude`` CLI rather than the raw ``anthropic`` SDK. See docs/DECISIONS.md D62.
"""

from __future__ import annotations

import importlib.util
import json
import os
import platform
import shutil
from pathlib import Path
from typing import Literal


class AuthError(RuntimeError):
    """No usable Claude credential, or the requested one is not reachable."""


#: A resolved authentication mode. ``"session"`` uses the Claude subscription login;
#: ``"api"`` uses a metered API credential.
AuthMode = Literal["session", "api"]

#: Variables that carry *metered* (API/cloud) authentication. Deliberately excludes
#: ``CLAUDE_CODE_OAUTH_TOKEN`` -- that token bills the plan, not the API. A run in ``"api"``
#: mode keeps these and drops the OAuth token; a ``"session"`` run does the reverse.
API_AUTH_ENV_VARS = (
    "ANTHROPIC_API_KEY",
    "ANTHROPIC_AUTH_TOKEN",
    "CLAUDE_CODE_USE_BEDROCK",
    "CLAUDE_CODE_USE_VERTEX",
)

#: The subscription (OAuth) token variable -- a portable ``claude setup-token`` credential.
SESSION_AUTH_ENV_VAR = "CLAUDE_CODE_OAUTH_TOKEN"


def claude_sdk_available() -> bool:
    """Whether the optional Claude Agent SDK is installed."""
    return importlib.util.find_spec("claude_agent_sdk") is not None


def claude_cli_dir() -> Path | None:
    """Return the directory holding a ``claude`` CLI on ``PATH``, if there is one.

    Only the ``PATH`` copy: this exists so :func:`agent_env` can keep that directory reachable
    after we shape the environment. It is *not* the answer to "can we run the agent" -- see
    :func:`claude_cli_found`.
    """
    found = shutil.which("claude")
    return Path(found).parent if found else None


def claude_cli_found() -> bool:
    """Whether the SDK will find a CLI to run, by the same rule the SDK itself uses.

    The SDK looks for its own bundled binary *first* and only then falls back to ``PATH``
    (``_internal/transport/subprocess_cli.py::_find_cli``). Recent wheels ship that binary, so
    ``pip install 'bbl[claude]'`` alone is a complete install and gating on ``shutil.which``
    would reject a setup that works perfectly.
    """
    if shutil.which("claude"):
        return True
    try:
        import claude_agent_sdk

        name = "claude.exe" if platform.system() == "Windows" else "claude"
        bundled = Path(claude_agent_sdk.__file__).parent / "_bundled" / name
        return bundled.is_file()
    except Exception:
        return False


def check_agent_cli() -> None:
    """Raise :class:`AuthError` unless the transport can actually start.

    The SDK is a Python wrapper that drives the ``claude`` CLI as a subprocess. Usually the
    wheel carries that binary; when it does not, it has to come from ``PATH``. The two failures
    have different fixes, so they are reported separately.
    """
    if not claude_sdk_available():
        raise AuthError(
            "the Claude Agent SDK is not installed -- run `pip install 'bbl[claude]'` before "
            "starting a design session"
        )
    if not claude_cli_found():
        raise AuthError(
            "no `claude` CLI available -- the Agent SDK runs it as a subprocess, and this "
            "install of the SDK does not bundle one. Install Claude Code from "
            "https://claude.com/claude-code so `claude` is on PATH, then run it once to log in."
        )


def _credentials_path() -> Path:
    """The user's Claude OAuth credentials file."""
    base = os.environ.get("CLAUDE_CONFIG_DIR") or Path.home() / ".claude"
    return Path(base) / ".credentials.json"


def _has_oauth_credentials() -> bool:
    """Whether the user's credentials file holds a Claude subscription (OAuth) login.

    ``claude`` login writes ``.credentials.json`` as ``{"claudeAiOauth": {...}}`` -- the
    presence of that object is what distinguishes a subscription login from a bare API-key
    setup, and it is the signal that "session usage" is available. A missing, unreadable, or
    API-only credentials file is not a subscription.
    """
    try:
        data = json.loads(_credentials_path().read_text())
    except (OSError, ValueError):
        return False
    return isinstance(data, dict) and isinstance(data.get("claudeAiOauth"), dict)


def session_auth_available() -> bool:
    """Whether the session could bill the Claude subscription ("session usage").

    True when a portable ``CLAUDE_CODE_OAUTH_TOKEN`` is set or a subscription OAuth login is on
    disk. An empty variable does not count.
    """
    if os.environ.get(SESSION_AUTH_ENV_VAR):
        return True
    return _has_oauth_credentials()


def api_auth_available() -> bool:
    """Whether the session could bill the Anthropic API (or a cloud provider) per token.

    True when any metered auth variable is set: a direct Anthropic key/token, or a
    Bedrock/Vertex routing flag whose own cloud credentials then apply. The subscription OAuth
    token is deliberately excluded -- it bills the plan, not the API.
    """
    return any(os.environ.get(var) for var in API_AUTH_ENV_VARS)


def resolve_auth_mode(requested: str = "auto") -> AuthMode:
    """Resolve a requested auth choice to the mode the session will use, or fail loudly.

    A preflight guard: it both validates that the chosen credential is actually reachable and
    reports which mode is live, so the choice is never silent. Without it an unauthenticated
    session would start, print a banner, and only fail on the first turn.

    Parameters
    ----------
    requested
        ``"auto"`` (prefer the subscription, else the API), ``"session"`` (force the
        subscription), or ``"api"`` (force the API).
    """
    if requested == "session":
        if not session_auth_available():
            raise AuthError(
                "no Claude subscription login found for --auth session. Log in with `claude` so "
                "~/.claude/.credentials.json exists, or set CLAUDE_CODE_OAUTH_TOKEN (from "
                "`claude setup-token`)."
            )
        return "session"
    if requested == "api":
        if not api_auth_available():
            raise AuthError(
                "no API credential found for --auth api. Set ANTHROPIC_API_KEY (or "
                "ANTHROPIC_AUTH_TOKEN), or enable a provider with CLAUDE_CODE_USE_BEDROCK / "
                "CLAUDE_CODE_USE_VERTEX."
            )
        return "api"

    # "auto": prefer the subscription (it's what the user is paying a flat rate for), fall
    # back to the metered API, and only then give up.
    if session_auth_available():
        return "session"
    if api_auth_available():
        return "api"
    raise AuthError(
        "no Claude credentials found -- an agent session cannot authenticate. Log in with "
        "`claude` so ~/.claude/.credentials.json exists (subscription), or set "
        "ANTHROPIC_API_KEY (or ANTHROPIC_AUTH_TOKEN / CLAUDE_CODE_OAUTH_TOKEN), or enable a "
        "provider with CLAUDE_CODE_USE_BEDROCK / CLAUDE_CODE_USE_VERTEX."
    )


def describe(mode: AuthMode) -> str:
    """How the resolved mode is reported to the user."""
    return "Claude subscription (session)" if mode == "session" else "Claude API key"


def agent_env(mode: AuthMode) -> dict[str, str]:
    """The environment for the ``claude`` subprocess, with exactly one auth path live.

    Everything the user has is carried through -- bbl is not a sandbox, and the agent needs the
    real PATH to reach ``claude`` and the real proxy/TLS settings to reach the network. The one
    thing shaped is the credential: the mode we are *not* using is neutralized, so billing is
    deterministic rather than a matter of which variable the CLI happens to check first.

    The cleared variables are set to ``""``, not deleted. The SDK builds the subprocess
    environment as ``{**os.environ, **options.env}`` -- it merges this mapping *over* the
    inherited environment -- so a credential we merely omit falls back through from
    ``os.environ`` unchanged and the session silently authenticates with the wrong one. An
    explicit empty value overrides the inherited one, which the CLI reads as unset.
    """
    env = dict(os.environ)
    for var in API_AUTH_ENV_VARS if mode == "session" else [SESSION_AUTH_ENV_VAR]:
        env[var] = ""

    cli_dir = claude_cli_dir()
    if cli_dir is not None:
        existing = env.get("PATH", "")
        env["PATH"] = os.pathsep.join([str(cli_dir), existing] if existing else [str(cli_dir)])

    # Project settings discovery walks *up* from cwd, auto-loading every CLAUDE.md it passes.
    # bbl runs wherever the user's plasmids are; that is not a place to inherit instructions
    # from. Paired with setting_sources=[] in agent.py, which covers the settings half.
    env["CLAUDE_CODE_DISABLE_CLAUDE_MDS"] = "1"
    return env
