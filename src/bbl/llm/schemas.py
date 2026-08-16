"""Structured outputs for the jobs that are single calls, not agent loops.

Parsing a request into a spec and mapping loose part names onto the catalog are extraction
tasks. A validated object beats a tool loop for both.
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class Spec(BaseModel):
    """A plasmid request, parsed from however the user phrased it."""

    backbone: str | None = Field(
        None, description="Vector family named or implied, e.g. 'pcDNA3.1'. Null if unstated."
    )
    parts: list[str] = Field(
        default_factory=list,
        description="Functional elements in the order they should appear, as the user named them.",
    )
    marker: str | None = Field(None, description="Selection marker, e.g. 'AmpR'. Null if unstated.")
    must_not_contain: list[str] = Field(
        default_factory=list,
        description="Elements the user explicitly said to exclude. Do not infer -- leave empty "
        "and record the question in `ambiguities` instead.",
    )
    purpose: str | None = Field(
        None, description="What the construct is for, if stated ('a control for my sensor')."
    )
    ambiguities: list[str] = Field(
        default_factory=list,
        description="Questions that would change the design if answered differently. Record "
        "anything you would otherwise have to guess -- especially whether features the user "
        "did not list must be absent or are simply don't-care.",
    )


class PartMatch(BaseModel):
    """One loose part name resolved against the inventory's vocabulary."""

    requested: str
    resolved_label: str | None = Field(
        None, description="The inventory feature label, or null if nothing matches."
    )
    plasmids: list[str] = Field(default_factory=list)
    confidence: float = Field(ge=0.0, le=1.0)
    note: str | None = Field(
        None, description="Why this match, or what makes it uncertain. Brief."
    )


class PartResolution(BaseModel):
    matches: list[PartMatch]
    unresolved: list[str] = Field(
        default_factory=list,
        description="Requested parts with no inventory match -- these need a sequence or a search.",
    )
