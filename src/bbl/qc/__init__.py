"""Sequence QC with swappable scoring backends.

Built to answer one question with evidence rather than opinion: **is Proto worth depending
on?** So the layer never assumes an answer -- it runs the same checks through a hand-written
backend and a Proto backend and reports where they agree, where they differ, and where only
one of them can answer at all.

    from bbl.qc import score, compare

    score(fragment)                  # native backend (always available)
    score(fragment, backend="proto") # Proto, or UNAVAILABLE if not installed
    print(compare(fragment))         # side-by-side table + disagreements

Deliberately **not** wired into ``excise``/``pcr``/``insert`` yet: bake-off first, adoption
second. See ``docs/DECISIONS.md`` and the plan for the reasoning.
"""

from __future__ import annotations

from .checks import CHECKS, CHECKS_BY_NAME, Check, resolve_params, tolerance_for
from .proto_backend import (
    PROTO_AVAILABLE,
    PROTO_IMPORT_ERROR,
    cpu_clean,
    proto_version,
)
from .proto_backend import describe as describe_constraint
from .proto_backend import score_one as _proto_score_one
from .types import (
    FAIL,
    NATIVE,
    PASS,
    PROTO,
    UNAVAILABLE,
    UNSUPPORTED,
    WARN,
    CheckResult,
    DiffReport,
    DiffRow,
    Scorecard,
)

__all__ = [
    "CHECKS",
    "CHECKS_BY_NAME",
    "Check",
    "CheckResult",
    "DiffReport",
    "DiffRow",
    "FAIL",
    "NATIVE",
    "PASS",
    "PROTO",
    "PROTO_AVAILABLE",
    "PROTO_IMPORT_ERROR",
    "Scorecard",
    "UNAVAILABLE",
    "UNSUPPORTED",
    "WARN",
    "compare",
    "cpu_clean",
    "describe_constraint",
    "proto_version",
    "registry_report",
    "score",
    "tolerance_for",
]


def _coerce(sequence) -> str:
    """Accept a str, ``Seq``, or ``Dseqrecord`` and return plain uppercase DNA."""
    return str(getattr(sequence, "seq", sequence)).upper()


def _native_score_one(check: Check, sequence: str, params: dict) -> CheckResult:
    if check.native is None:
        return CheckResult(
            check=check.name,
            backend=NATIVE,
            measurement=None,
            penalty=0.0,
            verdict=UNSUPPORTED,
            note="no native implementation; capability Proto adds",
        )
    try:
        measurement, penalty, detail = check.native(sequence, params)
    except Exception as exc:  # noqa: BLE001 -- one bad check must not abort the scorecard
        return CheckResult(
            check=check.name,
            backend=NATIVE,
            measurement=None,
            penalty=0.0,
            verdict=UNAVAILABLE,
            note=f"{type(exc).__name__}: {exc}",
        )
    return CheckResult(
        check=check.name,
        backend=NATIVE,
        measurement=measurement,
        penalty=penalty,
        verdict=PASS if penalty <= 0.0 else FAIL,
        detail=detail,
    )


def score(
    sequence,
    backend: str = NATIVE,
    label: str = "sequence",
    params: dict | None = None,
    checks: list[str] | None = None,
) -> Scorecard:
    """Score one sequence with one backend.

    Args:
        sequence: str, ``Seq`` or ``Dseqrecord``.
        backend: ``"native"`` or ``"proto"``.
        label: name to show in reports.
        params: ``{check_name: {param: value}}`` overrides on top of each check's defaults.
        checks: restrict to these check names; default is all of them.
    """
    if backend not in (NATIVE, PROTO):
        raise ValueError(f"unknown backend {backend!r}; expected {NATIVE!r} or {PROTO!r}")
    seq = _coerce(sequence)
    selected = [CHECKS_BY_NAME[name] for name in checks] if checks else CHECKS
    overrides = params or {}

    results = []
    for check in selected:
        resolved = resolve_params(check, overrides.get(check.name))
        if backend == NATIVE:
            results.append(_native_score_one(check, seq, resolved))
        else:
            results.append(_proto_score_one(check, seq, resolved))
    return Scorecard(backend=backend, label=label, length=len(seq), results=results)


def compare(
    sequence,
    label: str = "sequence",
    params: dict | None = None,
    checks: list[str] | None = None,
) -> DiffReport:
    """Score with both backends and diff them.

    A nonzero ``disagreements`` list means the two backends computed *different values for the
    same physical quantity* -- that is a bug in one of them. Checks that only one backend can
    do are reported as coverage, not as disagreement.
    """
    seq = _coerce(sequence)
    native_card = score(seq, backend=NATIVE, label=label, params=params, checks=checks)
    proto_card = score(seq, backend=PROTO, label=label, params=params, checks=checks)
    native_by_name = native_card.by_name()
    proto_by_name = proto_card.by_name()
    rows = [
        DiffRow(check=name, native=native_by_name[name], proto=proto_by_name[name])
        for name in native_by_name
    ]
    return DiffReport(label=label, length=len(seq), rows=rows)


def registry_report() -> list[dict]:
    """Proto registry metadata for every check that names a constraint.

    The artifact for the adopt/reject decision: shows which constraints are CPU-clean and,
    for the ones that are not, exactly what they would drag in.
    """
    if not PROTO_AVAILABLE:
        return [
            {"key": check.proto_key, "cpu_clean": False, "cpu_clean_reason": PROTO_IMPORT_ERROR}
            for check in CHECKS
            if check.proto_key
        ]
    rows = []
    for check in CHECKS:
        if not check.proto_key:
            continue
        try:
            rows.append({"check": check.name, **describe_constraint(check.proto_key)})
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    "check": check.name,
                    "key": check.proto_key,
                    "cpu_clean": False,
                    "cpu_clean_reason": f"{type(exc).__name__}: {exc}",
                }
            )
    return rows
