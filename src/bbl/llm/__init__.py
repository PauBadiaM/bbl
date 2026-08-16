"""LLM interface: the design tools, structured schemas, and the interactive session.

``session``, ``schemas`` and ``prompts`` import cleanly without ``claude-agent-sdk``, so the
whole model-facing contract is testable offline. ``tools``, ``agent`` and ``chat`` need it.
"""

from .prompts import SYSTEM, build_system_blocks, build_system_prompt, inventory_digest
from .schemas import PartMatch, PartResolution, Spec
from .session import DesignSession, Product

__all__ = [
    "DesignSession",
    "PartMatch",
    "PartResolution",
    "Product",
    "SYSTEM",
    "Spec",
    "build_system_blocks",
    "build_system_prompt",
    "inventory_digest",
]
