"""Bench-ready cloning reports: procedure, reagent tables and plasmid maps in one HTML file.

    from bbl.report import write_report

    plan = excise_features("pHL391.dna", ["NFKBRE"])
    write_report(plan, "pCLM1_report.html", parent=load_plasmid("pHL391.dna"),
                 aim="A sensor control with the NF-kB response element removed.")

The report is what a plan looks like once it has to survive contact with a pipette: reaction
volumes, gel percentages, what to write down afterwards. It is one self-contained HTML file --
figures inline, nothing fetched from the network.
"""

from __future__ import annotations

from pathlib import Path

from .bench import BENCH_DEFAULTS, Component, Table, bench_config
from .build import DesignReport, Figure, Step, build_report
from .maps import circular_svg, linear_svg

__all__ = [
    "BENCH_DEFAULTS",
    "Component",
    "DesignReport",
    "Figure",
    "Step",
    "Table",
    "bench_config",
    "build_report",
    "circular_svg",
    "linear_svg",
    "write_report",
]


def write_report(plan, path, **kwargs) -> Path:
    """Build a report for ``plan`` and write it to ``path``. Returns the path written.

    Keyword arguments are passed through to :func:`~bbl.report.build.build_report`.
    """
    path = Path(path)
    path.write_text(build_report(plan, **kwargs).to_html())
    return path
