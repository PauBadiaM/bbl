"""Terminal rendering of an agent turn, as it happens.

A design turn can run for minutes across a dozen tool calls. Printing nothing until it finishes
makes the session look hung, so every message is rendered as it arrives:

    · thinking: the NFKBRE deletion is the cheaper route, but ...
    → inspect_plasmid: pHL391
    ✓ inspect_plasmid (412 chars)
    <the assistant's text>
    ● 6 turns · $0.31 · 42.1s

Tool *results* are reported by size, never inlined: they are JSON payloads hundreds of lines
long, and the model's prose is the part worth reading. Trimmed from ``acumen``'s ``logs.py``.
"""

from __future__ import annotations

import shutil
import sys

#: Arguments worth showing next to a tool call, most specific first. A design tool's first
#: argument is nearly always the thing being worked on.
_SUMMARY_KEYS = (
    "plasmid",
    "backbone",
    "product_id",
    "query",
    "name",
    "path",
    "file_path",
    "command",
    "url",
    "pattern",
)

_THINKING_CAP = 240
_ARG_CAP = 160


class Style:
    """ANSI codes, or empty strings when the output is not a terminal."""

    def __init__(self, enabled: bool | None = None):
        if enabled is None:
            enabled = sys.stdout.isatty()
        self.dim = "\033[2m" if enabled else ""
        self.bold = "\033[1m" if enabled else ""
        self.reset = "\033[0m" if enabled else ""


def _clip(text: str, cap: int) -> str:
    text = " ".join(str(text).split())
    return text if len(text) <= cap else text[: cap - 1] + "…"


def summarize_args(args) -> str:
    """The one argument worth showing beside a tool call.

    Scalars only. The fallback used to print whatever the first argument was, which for
    ``ask_user`` meant a line of Python-repr'd question objects immediately above the same
    questions rendered properly -- noise where the tool speaks for itself.
    """
    if not isinstance(args, dict) or not args:
        return ""
    for key in _SUMMARY_KEYS:
        if isinstance(args.get(key), (str, int, float)) and args[key]:
            return _clip(args[key], _ARG_CAP)
    first = next((v for v in args.values() if v and isinstance(v, (str, int, float))), None)
    return _clip(first, _ARG_CAP) if first is not None else ""


class Renderer:
    """Prints SDK messages as they stream, and collects the assistant's text.

    Deliberately duck-typed on block ``.type``/attribute presence rather than ``isinstance``
    against SDK classes, so it can be driven by a fake message stream in a test with the SDK
    absent.
    """

    def __init__(self, style: Style | None = None, out=None):
        self.style = style or Style()
        self.out = out or sys.stdout
        self.text_parts: list[str] = []
        self.tool_calls: list[tuple[str, dict]] = []
        self.session_id: str | None = None
        self.result = None

    @property
    def text(self) -> str:
        return "\n".join(self.text_parts).strip()

    def _line(self, text: str, dim: bool = True) -> None:
        style = self.style
        prefix, suffix = (style.dim, style.reset) if dim else ("", "")
        print(f"{prefix}{text}{suffix}", file=self.out, flush=True)

    def handle(self, message) -> None:
        """Render one SDK message."""
        from .tools import unqualify

        # Every message carries the id, and the init event arrives before any content. Taking
        # it here rather than only from the result means an interrupted turn is still
        # resumable -- which is exactly the turn you most want back.
        self.session_id = getattr(message, "session_id", None) or self.session_id

        for block in getattr(message, "content", None) or []:
            kind = getattr(block, "type", None)
            if kind == "thinking" or hasattr(block, "thinking"):
                thought = getattr(block, "thinking", "")
                if thought:
                    self._line(f"  · {_clip(thought, _THINKING_CAP)}")
            elif kind == "text" or hasattr(block, "text"):
                text = getattr(block, "text", "")
                if text:
                    self.text_parts.append(text)
                    print(f"\n{text}\n", file=self.out, flush=True)
            elif kind == "tool_use" or (hasattr(block, "name") and hasattr(block, "input")):
                name = unqualify(getattr(block, "name", "?"))
                args = getattr(block, "input", None) or {}
                self.tool_calls.append((name, args))
                detail = summarize_args(args)
                self._line(f"  → {name}{': ' + detail if detail else ''}")
            elif kind == "tool_result" or hasattr(block, "tool_use_id"):
                content = getattr(block, "content", None)
                size = len(str(content)) if content is not None else 0
                mark = "✗" if getattr(block, "is_error", False) else "✓"
                self._line(f"  {mark} {size} chars")

        if getattr(message, "total_cost_usd", None) is not None or hasattr(message, "num_turns"):
            self.result = message
            self._line(f"  ● {self.summary(message)}")

    def summary(self, result) -> str:
        parts = []
        turns = getattr(result, "num_turns", None)
        if turns:
            parts.append(f"{turns} turn{'s' if turns != 1 else ''}")
        cost = getattr(result, "total_cost_usd", None)
        if cost:
            parts.append(f"${cost:.2f}")
        duration = getattr(result, "duration_ms", None)
        if duration:
            parts.append(f"{duration / 1000:.1f}s")
        return " · ".join(parts) or "done"

    def rule(self, label: str) -> None:
        """A full-width separator, used to set the harness-rendered protocol apart."""
        width = min(shutil.get_terminal_size((80, 24)).columns, 88)
        pad = max(width - len(label) - 4, 0)
        self._line(f"── {label} {'─' * pad}")
