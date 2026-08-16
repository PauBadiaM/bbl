"""Credential resolution: an active Claude session first, an API key as fallback.

Order, first match wins:

1. **Session token** -- ``ANTHROPIC_AUTH_TOKEN`` if already exported, otherwise one minted from
   an active ``ant auth login`` profile. Short-lived and refreshable, so nothing durable sits
   on disk. This is the right shape on a shared cluster.
2. **API key** -- ``ANTHROPIC_API_KEY``. A fallback: it lives on NFS, leaks into job logs and
   ``env`` dumps, and never expires.
3. **Ask the user to log in**, and run ``ant auth login`` for them if they say yes.

Two things the SDK (0.122) does not do for you, both handled here: it reads only the two
environment variables -- a bare ``Anthropic()`` does *not* pick up an ``ant`` profile -- and it
does not attach the ``oauth-2025-04-20`` beta header that ``/v1/messages`` requires for a
bearer token.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass

OAUTH_BETA = "oauth-2025-04-20"
SESSION, API_KEY = "session", "api_key"
_ANT = "ant"

INSTALL_HINT = (
    "Install the CLI first:\n"
    "    go install github.com/anthropics/anthropic-cli/cmd/ant@latest\n"
    "  or grab a release from github.com/anthropics/anthropic-cli/releases"
)


class NoCredentials(RuntimeError):
    """No usable credential, and the user did not log in."""


@dataclass
class Credential:
    """A live client plus which kind of credential is behind it."""

    client: object
    kind: str
    description: str

    @property
    def refreshable(self) -> bool:
        """Only a session token can be re-minted; an API key that fails is simply wrong."""
        return self.kind == SESSION


def ant_available() -> bool:
    return shutil.which(_ANT) is not None


def session_token(profile: str | None = None) -> str | None:
    """Mint a short-lived access token from an active ``ant`` session, or ``None``.

    ``--access-token`` is required: with no flag the CLI prints the whole credentials JSON,
    which is not a bearer token and fails with an opaque protocol error if used as one.
    Returns ``None`` when the CLI is absent or no session is active.
    """
    if not ant_available():
        return None
    command = [_ANT, "auth", "print-credentials", "--access-token"]
    if profile:
        command += ["--profile", profile]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    token = result.stdout.strip()
    return token if result.returncode == 0 and token else None


def session_active(profile: str | None = None) -> bool:
    """True when a Claude session can currently produce a token."""
    return session_token(profile) is not None


def login(profile: str | None = None) -> bool:
    """Run ``ant auth login``, inheriting the terminal so the browser/code flow works.

    Uses ``--no-browser`` when there is no display -- the usual case over SSH on a cluster,
    where it prints a URL to open elsewhere and takes the code back in the terminal.
    """
    if not ant_available():
        print(INSTALL_HINT)
        return False
    command = [_ANT, "auth", "login"]
    if profile:
        command += ["--profile", profile]
    if not os.environ.get("DISPLAY") and not os.environ.get("WAYLAND_DISPLAY"):
        command.append("--no-browser")
    try:
        return subprocess.run(command).returncode == 0
    except (OSError, subprocess.SubprocessError) as exc:
        print(f"could not run `ant auth login`: {exc}")
        return False


def _client(token: str | None = None):
    import anthropic

    if token is None:
        return anthropic.Anthropic()
    return anthropic.Anthropic(auth_token=token, default_headers={"anthropic-beta": OAUTH_BETA})


def resolve(profile: str | None = None, ask=None) -> Credential:
    """Find a credential. ``ask(prompt) -> bool`` enables the interactive login offer."""
    exported = os.environ.get("ANTHROPIC_AUTH_TOKEN")
    if exported:
        return Credential(_client(exported), SESSION, "session token (ANTHROPIC_AUTH_TOKEN)")

    token = session_token(profile)
    if token:
        name = f" ({profile})" if profile else ""
        return Credential(_client(token), SESSION, f"active Claude session{name}")

    if os.environ.get("ANTHROPIC_API_KEY"):
        return Credential(_client(), API_KEY, "ANTHROPIC_API_KEY (static key, no session)")

    if ask is not None and ant_available():
        if ask("No active Claude session. Log in now?") and login(profile):
            token = session_token(profile)
            if token:
                return Credential(_client(token), SESSION, "active Claude session (new)")

    message = [
        "No active Claude session and no API key.",
        "  Preferred:  ant auth login          (short-lived session token)",
        "  Fallback:   export ANTHROPIC_API_KEY=sk-ant-...",
    ]
    if not ant_available():
        message.append("  " + INSTALL_HINT.replace("\n", "\n  "))
    raise NoCredentials("\n".join(message))


def verify(credential: Credential) -> str:
    """Cheap round trip proving the credential works, before the first real turn."""
    import anthropic

    try:
        credential.client.models.list(limit=1)
    except anthropic.AuthenticationError as exc:
        raise NoCredentials(f"{credential.description} was rejected: {exc}") from exc
    return credential.description


def with_refresh(call, credential: Credential, profile: str | None = None):
    """Run ``call(client)``, re-minting an expired session token once and retrying.

    ``credential`` is mutated in place so the caller keeps using the refreshed client. An API
    key is never retried: it cannot expire, so a failure means something else is wrong and
    retrying would only hide it.
    """
    import anthropic

    try:
        return call(credential.client)
    except anthropic.AuthenticationError:
        if not credential.refreshable:
            raise
        token = session_token(profile)
        if not token:
            raise
        credential.client = _client(token)
        credential.description = "active Claude session (refreshed)"
        return call(credential.client)


def describe(profile: str | None = None) -> str:
    """Which source would be used, without building a client or revealing a secret."""
    if os.environ.get("ANTHROPIC_AUTH_TOKEN"):
        return "session token (ANTHROPIC_AUTH_TOKEN)"
    if session_active(profile):
        return "active Claude session"
    if os.environ.get("ANTHROPIC_API_KEY"):
        return "ANTHROPIC_API_KEY (static key, no session)"
    return "none"
