"""LLM interface: tool wrappers, structured schemas, and the design loop.

``session`` and ``schemas`` import cleanly without the Anthropic SDK, so the whole boundary is
testable offline. ``tools`` needs it.
"""

from .prompts import SYSTEM, build_system_blocks, inventory_digest
from .schemas import PartMatch, PartResolution, Spec
from .session import DesignSession, Product


def chat(*args, **kwargs):
    """Start an interactive design session (lazy import: needs the Anthropic SDK)."""
    from .chat import run

    return run(*args, **kwargs)

__all__ = [
    "DesignSession",
    "PartMatch",
    "PartResolution",
    "Product",
    "SYSTEM",
    "Spec",
    "build_system_blocks",
    "chat",
    "inventory_digest",
]
