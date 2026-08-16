"""Proto-backed scoring, in score-only mode.

Proto's constraint functions are plain callables:

    constraint_function(input_sequences: list[tuple[Sequence, ...]], config) -> list[ConstraintOutput]

``ConstraintOutput`` has exactly two fields, ``score`` (a penalty in ``[0, 1]``, ``0.0`` ==
satisfied) and ``metadata`` (the raw measurement). Nothing here constructs a ``Program``,
``Optimizer`` or ``Generator``: no optimisation, no GPU, no ``HF_TOKEN``, and no RNG, so
results are reproducible and diffable against the native backend.

Two guards matter:

**The CPU gate.** Proto self-declares cost on every constraint (``uses_gpu``,
``tools_called``). We refuse to run any constraint that declares either, so this layer cannot
silently start downloading model weights or shelling out to a binary mid-demo. Verified
casualties of that gate, both of which we would otherwise have wanted:

* ``longest-orf-length`` -- ``tools_called=["orfipy-prediction"]``
* ``seq-motif`` -- ``tools_called=["meme-fimo-scan"]`` (and it is PWM scanning via FIMO, not
  the exact-substring test one might assume from the name)

**The import guard.** Proto is an optional extra with a large dependency tree
(``proto_tools`` from git, rdkit, biotite, jupyter, scikit-learn). If it is absent or its
import fails for any reason, every Proto result is ``UNAVAILABLE`` and the native backend
still scores. ``bbl`` must never hard-depend on Proto -- rejecting it is a permitted outcome.
"""

from __future__ import annotations

from .types import FAIL, PASS, PROTO, UNAVAILABLE, UNSUPPORTED, CheckResult

try:  # pragma: no cover -- exercised by whichever environment runs the tests
    from proto_language.constraint import ConstraintRegistry
    from proto_language.core import Sequence

    PROTO_AVAILABLE = True
    PROTO_IMPORT_ERROR: str | None = None
except Exception as exc:  # noqa: BLE001 -- a heavy dep tree fails in many ways, not just ImportError
    ConstraintRegistry = None  # type: ignore[assignment]
    Sequence = None  # type: ignore[assignment]
    PROTO_AVAILABLE = False
    PROTO_IMPORT_ERROR = f"{type(exc).__name__}: {exc}"


def proto_version() -> str | None:
    """Installed Proto version, or ``None`` when unavailable."""
    if not PROTO_AVAILABLE:
        return None
    try:
        from importlib.metadata import version

        return version("proto-language")
    except Exception:  # noqa: BLE001
        return "unknown"


def cpu_clean(key: str) -> tuple[bool, str | None]:
    """Whether a Proto constraint is safe to run here: no GPU, no external tools.

    Returns ``(ok, reason_if_not)``. Raises nothing -- an unknown key is simply not clean.
    """
    if not PROTO_AVAILABLE:
        return False, "proto not installed"
    try:
        spec = ConstraintRegistry.get(key)
    except Exception as exc:  # noqa: BLE001
        return False, f"not in registry ({type(exc).__name__})"
    if getattr(spec, "uses_gpu", False):
        return False, "declares uses_gpu=True"
    tools = list(getattr(spec, "tools_called", []) or [])
    if tools:
        return False, f"calls external tools: {', '.join(tools)}"
    return True, None


def describe(key: str) -> dict:
    """Registry metadata for one constraint -- used by tests and the report header."""
    spec = ConstraintRegistry.get(key)
    ok, reason = cpu_clean(key)
    return {
        "key": key,
        "uses_gpu": getattr(spec, "uses_gpu", None),
        "tools_called": list(getattr(spec, "tools_called", []) or []),
        "category": getattr(spec, "category", None),
        "supported_sequence_types": list(getattr(spec, "supported_sequence_types", []) or []),
        "config_model": getattr(spec, "config_model", None).__name__
        if getattr(spec, "config_model", None)
        else None,
        "cpu_clean": ok,
        "cpu_clean_reason": reason,
    }


def _unavailable(check_name: str, note: str, verdict: str = UNAVAILABLE) -> CheckResult:
    return CheckResult(
        check=check_name,
        backend=PROTO,
        measurement=None,
        penalty=0.0,
        verdict=verdict,
        detail={},
        note=note,
    )


def score_one(check, sequence: str, params: dict) -> CheckResult:
    """Evaluate one check through Proto.

    Never raises: every failure mode becomes an ``UNAVAILABLE``/``UNSUPPORTED`` result so a
    single unsupported check cannot abort a whole scorecard.
    """
    if check.proto_key is None:
        return _unavailable(
            check.name, "no Proto constraint expresses this check", verdict=UNSUPPORTED
        )
    if not PROTO_AVAILABLE:
        return _unavailable(check.name, f"proto not importable ({PROTO_IMPORT_ERROR})")

    ok, reason = cpu_clean(check.proto_key)
    if not ok:
        return _unavailable(
            check.name, f"{check.proto_key}: {reason}", verdict=UNSUPPORTED
        )

    spec = ConstraintRegistry.get(check.proto_key)
    function = getattr(spec, "function", None)
    config_model = getattr(spec, "config_model", None)
    if function is None or config_model is None:
        return _unavailable(check.name, f"{check.proto_key}: registry entry incomplete")

    # Build the config from our thresholds. A pydantic ValidationError here means our
    # threshold mapping is wrong, which is our bug and worth surfacing loudly.
    try:
        config = config_model(**check.proto_config(params))
    except Exception as exc:  # noqa: BLE001
        return _unavailable(check.name, f"config rejected: {type(exc).__name__}: {exc}")

    # Proto's DNA alphabet is frozenset("ACGT"); ambiguity codes are not representable, and
    # both Sequence construction and the constraint itself can reject them.
    try:
        seq = Sequence(sequence.upper(), "dna")
        outputs = function([(seq,)], config)
    except Exception as exc:  # noqa: BLE001
        return _unavailable(
            check.name,
            f"{check.proto_key} rejected the input: {type(exc).__name__}: {exc}",
            verdict=UNSUPPORTED,
        )

    if not outputs:
        return _unavailable(check.name, f"{check.proto_key} returned no output")

    output = outputs[0]
    metadata = dict(getattr(output, "metadata", {}) or {})
    penalty = float(getattr(output, "score", 0.0))
    measurement = check.proto_measure(metadata, params) if check.proto_measure else None

    return CheckResult(
        check=check.name,
        backend=PROTO,
        measurement=measurement,
        penalty=penalty,
        # Proto's score is a penalty, not a verdict. We apply our own threshold so that
        # swapping backends cannot silently turn a pass into a fail.
        verdict=PASS if penalty <= 0.0 else FAIL,
        detail=metadata,
    )
