"""Result types for the QC layer.

The layer exists to answer one question: *should we depend on Proto?* Everything here is
shaped by that. In particular a result carries **three** separate things, because collapsing
them is what makes backend comparisons meaningless:

``measurement``
    The physical quantity -- GC percent, the longest homopolymer run, a count. Two backends
    computing the same quantity **must** agree; a disagreement is a bug in one of them.

``penalty``
    Graded badness in ``[0, 1]``, ``0.0`` == satisfied. This is a *judgement*, not a
    measurement: backends may agree on the number and still disagree on how bad it is.
    Proto's composition constraints are often log- or ratio-scaled where a hand-written check
    is a step function, so penalty differences are expected and informative rather than wrong.

``verdict``
    Ours, never the backend's, applied from our own thresholds so that swapping backends
    cannot silently change a pass into a fail.
"""

from __future__ import annotations

from dataclasses import dataclass, field

#: Backend identifiers.
NATIVE = "native"
PROTO = "proto"

#: Verdicts. ``UNAVAILABLE`` means the backend could not answer at all -- Proto not
#: installed, or the check has no implementation on that backend. It is deliberately
#: distinct from ``FAIL``: "we don't know" must never read as "it's bad".
PASS = "pass"
WARN = "warn"
FAIL = "fail"
UNAVAILABLE = "unavailable"
UNSUPPORTED = "unsupported"


@dataclass(frozen=True)
class CheckResult:
    """One check, evaluated by one backend, on one sequence."""

    check: str
    backend: str
    measurement: float | None
    penalty: float
    verdict: str
    detail: dict = field(default_factory=dict)
    note: str | None = None

    @property
    def ok(self) -> bool:
        """True only for an affirmative pass. ``UNAVAILABLE`` is not ``ok``."""
        return self.verdict == PASS

    def __str__(self) -> str:
        if self.verdict in (UNAVAILABLE, UNSUPPORTED):
            return f"{self.check} [{self.backend}]: {self.verdict}" + (
                f" -- {self.note}" if self.note else ""
            )
        value = "n/a" if self.measurement is None else f"{self.measurement:g}"
        return (
            f"{self.check} [{self.backend}]: {value} "
            f"(penalty {self.penalty:.3f}, {self.verdict})"
        )


@dataclass(frozen=True)
class Scorecard:
    """Every check for one backend on one sequence."""

    backend: str
    label: str
    length: int
    results: list[CheckResult] = field(default_factory=list)

    def by_name(self) -> dict[str, CheckResult]:
        return {r.check: r for r in self.results}

    def get(self, name: str) -> CheckResult | None:
        return self.by_name().get(name)

    @property
    def failures(self) -> list[CheckResult]:
        return [r for r in self.results if r.verdict == FAIL]

    @property
    def answered(self) -> list[CheckResult]:
        """Results where the backend actually produced a number."""
        return [r for r in self.results if r.verdict not in (UNAVAILABLE, UNSUPPORTED)]

    def __str__(self) -> str:
        lines = [f"{self.label} ({self.length} bp) -- backend={self.backend}"]
        lines += [f"  {r}" for r in self.results]
        return "\n".join(lines)


@dataclass(frozen=True)
class DiffRow:
    """One check compared across both backends."""

    check: str
    native: CheckResult
    proto: CheckResult

    @property
    def comparable(self) -> bool:
        """Both backends produced a measurement, so a delta means something."""
        return (
            self.native.measurement is not None
            and self.proto.measurement is not None
            and self.native.verdict not in (UNAVAILABLE, UNSUPPORTED)
            and self.proto.verdict not in (UNAVAILABLE, UNSUPPORTED)
        )

    @property
    def delta(self) -> float | None:
        """Absolute difference in measurement, or ``None`` if not comparable."""
        if not self.comparable:
            return None
        assert self.native.measurement is not None and self.proto.measurement is not None
        return abs(self.native.measurement - self.proto.measurement)

    @property
    def penalty_delta(self) -> float | None:
        if self.native.verdict in (UNAVAILABLE, UNSUPPORTED):
            return None
        if self.proto.verdict in (UNAVAILABLE, UNSUPPORTED):
            return None
        return abs(self.native.penalty - self.proto.penalty)

    @property
    def verdicts_agree(self) -> bool | None:
        if not self.comparable:
            return None
        return self.native.verdict == self.proto.verdict

    @property
    def coverage(self) -> str:
        """Which backends can do this check at all -- the adopt/reject signal."""
        native_ok = self.native.verdict not in (UNAVAILABLE, UNSUPPORTED)
        proto_ok = self.proto.verdict not in (UNAVAILABLE, UNSUPPORTED)
        if native_ok and proto_ok:
            return "both"
        if native_ok:
            return "native-only"
        if proto_ok:
            return "proto-only"
        return "neither"


@dataclass(frozen=True)
class DiffReport:
    """Backend comparison for one sequence.

    ``disagreements`` is the thing to look at: checks where both backends produced a
    measurement but the numbers differ by more than the check's tolerance. Those are bugs.
    Checks covered by only one backend are *not* disagreements -- they are capability gaps,
    reported separately.
    """

    label: str
    length: int
    rows: list[DiffRow] = field(default_factory=list)

    @property
    def disagreements(self) -> list[DiffRow]:
        from .checks import tolerance_for

        return [
            row
            for row in self.rows
            if row.comparable and (row.delta or 0.0) > tolerance_for(row.check)
        ]

    def coverage(self, kind: str) -> list[DiffRow]:
        return [row for row in self.rows if row.coverage == kind]

    def summary(self) -> str:
        header = f"{self.label} ({self.length} bp)"
        width = max((len(r.check) for r in self.rows), default=5) + 1
        lines = [
            header,
            "-" * max(len(header), 68),
            f"{'check'.ljust(width)} {'native':>12} {'proto':>12} {'delta':>10}  coverage",
        ]
        for row in self.rows:
            native = (
                "n/a" if row.native.measurement is None else f"{row.native.measurement:.4g}"
            )
            proto = "n/a" if row.proto.measurement is None else f"{row.proto.measurement:.4g}"
            delta = "-" if row.delta is None else f"{row.delta:.4g}"
            lines.append(
                f"{row.check.ljust(width)} {native:>12} {proto:>12} {delta:>10}  {row.coverage}"
            )
        both = len(self.coverage("both"))
        lines.append("")
        lines.append(
            f"comparable: {both}   native-only: {len(self.coverage('native-only'))}   "
            f"proto-only: {len(self.coverage('proto-only'))}   "
            f"disagreements: {len(self.disagreements)}"
        )
        for row in self.disagreements:
            lines.append(
                f"  ! {row.check}: native={row.native.measurement!r} "
                f"proto={row.proto.measurement!r} (delta {row.delta:.6g})"
            )
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summary()
