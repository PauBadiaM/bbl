# Proto bake-off: `bbl.qc`

Evidence for one decision: **should `bbl` depend on Proto** (Arc Institute's
`proto-language`) for sequence QC? The layer was built so that Proto-backed and hand-written
scorers are interchangeable, run on the same inputs, and produce a numeric diff.

Nothing is wired into `excise` / `pcr` / `insert`. Bake-off first, adoption second.

## How to run it

Native only (the default environment, no Proto):

```bash
python -m pytest tests/test_qc_native.py tests/test_qc_proto.py tests/test_qc_compare.py -q
# 18 passed, 19 skipped
```

With Proto (see the environment notes below for why it is a pinned venv + `PYTHONPATH`):

```bash
PYTHONPATH=<proto-language-src>:<proto-tools-src> \
  $SCRATCH/bbl_proto_v3/bin/python -m pytest tests/test_qc{_native,_proto,_compare}.py -q
# 36 passed, 1 skipped
```

```python
from bbl.qc import score, compare, registry_report
print(compare(fragment))          # side-by-side table + disagreements
print(score(fragment))            # native
print(score(fragment, backend="proto"))
```

## Result

**Where both backends can measure, they agree exactly.** Sweep over all 28 inventory
plasmids (first 600 bp each): **112 comparable checks, 0 disagreements**. The only nonzero
delta anywhere is `gc_content` at ~7e-15, i.e. float noise, well inside tolerance.

Coverage is asymmetric, and the asymmetry is the actual finding:

| Check | Coverage | Note |
|---|---|---|
| `gc_content` | both | agrees with Biopython `gc_fraction` |
| `max_homopolymer` | both | agrees with `complexity.longest_homopolymer` |
| `sequence_length` | both | |
| `specific_kmer` | both | literal k-mers only |
| `restriction_sites` | **native only** | Proto cannot express ambiguous sites |
| `longest_repeat` | **native only** | no Proto equivalent; catches the 8x BoxB array |
| `repeat_fraction` | **native only** | no Proto equivalent |
| `kmer_frequency` | **proto only** | capability `bbl` lacks |
| `dinucleotide_composition` | **proto only** | capability `bbl` lacks |

### What Proto genuinely adds

1. **Graded penalties instead of step functions.** Both backends measure a 12-base run
   identically, but native saturates at 1.0 while Proto returns ~0.53 (log-scaled). Anywhere
   candidates get *ranked* rather than merely accepted — primer length selection, choosing
   among sourcing routes — that gradient is real signal the incumbent throws away.
   Regression: `test_the_two_backends_disagree_on_penalty_shape_even_when_measurements_match`.
2. **`kmer_frequency` and `dinucleotide_composition`**, which `bbl` does not have and which
   map onto vendor synthesis-complexity criteria — the `complexity.py` placeholder's job.
3. **Self-describing configs and cost declarations** (`uses_gpu`, `tools_called`), which is
   what makes the CPU gate below possible at all.

### What Proto cannot do

1. **Ambiguous recognition sites.** Proto's DNA alphabet is `frozenset("ACGT")`, so AccI's
   `GT^MKAC` is not representable. It raises `ValueError` rather than silently miscounting —
   the good failure mode — but the Biopython `Restriction` route in `enzymes.py` must be kept.
   `GTATAC` is a real AccI site that no literal k-mer search will find.
2. **Frame integrity.** `longest-orf-length` declares `tools_called=["orfipy-prediction"]`.
3. **Motif survival.** `seq-motif` declares `tools_called=["meme-fimo-scan"]` and is PWM
   scanning via FIMO — *not* the exact-substring test the name suggests. So it is the wrong
   tool for "did the `oFH155` / `TruseqR2` / `10XCS1` readout handle survive surgery"; a plain
   substring check is correct there and Proto offers nothing better.

Both (2) and (3) are refused by the **CPU gate**: the layer resolves every constraint through
`ConstraintRegistry` and refuses any that declares `uses_gpu` or `tools_called`, so it cannot
silently start downloading weights or shelling out mid-demo. All six adopted constraints are
`gpu=False tools=[]`. `test_the_gate_actually_rejects_the_constraints_we_wanted` fails loudly
if a future Proto release makes either pure-Python — at which point we should adopt it.

## Environment: Proto does not install cleanly here

Non-trivial, and it belongs in any cost estimate.

**Root cause: this cluster is glibc 2.17** (CentOS 7); the newest supported wheel tag is
`manylinux_2_17`. Proto's dependency floors are unpinned (`numpy>=1.20`, `scipy>=1.10`,
`pandas>=1.5`), so pip resolves to the latest releases — `numpy 2.5.2`, `scipy 1.18.0`,
`pandas 3.0.5` — which publish only `manylinux_2_28`+ wheels. Pip falls back to sdists,
and the scipy source build dies on `Dependency "OpenBLAS" not found`. Same failure class as
`DECISIONS.md` D3 (`pydivsufsort` sdist-only).

The wall moves one dependency at a time: scipy → Pillow (`libjpeg` headers missing) → rdkit
(no release above `2024.3.2` for this platform). Working pin set:

```bash
python -m venv --system-site-packages $SCRATCH/bbl_proto_v3      # inherit conda's numpy/pydantic/biopython
$SCRATCH/bbl_proto_v3/bin/pip install docstring_parser filelock pyyaml requests
$SCRATCH/bbl_proto_v3/bin/pip install "pandas==2.2.3"
$SCRATCH/bbl_proto_v3/bin/pip install "numpy==2.2.6" "scipy==1.14.1" "biotite==1.0.1" \
                                      "Pillow==10.4.0" "rdkit==2024.3.2"
```

Plus **both** sources on `PYTHONPATH`: `proto-language` and `proto-tools` (a git submodule
that is absent from the sdist, and whose own dependency list adds numba, umap-learn, gemmi
and modal — none of which the composition path actually imports).

**`proto_language.constraint.__init__` imports every constraint eagerly**, including the
protein-structure ones. So there is no lightweight subset install: a 20-line GC-content
function is gated behind biotite and rdkit. This is the single biggest argument against the
dependency, and it is a packaging property, not a scientific one.

Because of all this, **no `[proto]` extra was added to `pyproject.toml`** — a plain
`pip install bbl[proto]` would fail on this cluster and an extra that does not work is worse
than none. `bbl.qc` treats Proto as strictly optional: absent or broken, every Proto result is
`UNAVAILABLE`, the native backend still scores, and `test_bbl_qc_imports_without_proto`
guarantees the import never requires it.

## Recommendation

**Adopt narrowly, do not depend broadly.** Take Proto's composition constraints for the
`complexity.py` placeholder, where graded penalties and the two proto-only checks are genuine
gains and the inputs are short synthetic fragments over plain ACGT. Keep the native Biopython
route for anything restriction- or annotation-aware. Do not put Proto on `bbl`'s install path
while the eager-import + glibc-2.17 situation stands — keep it an optional backend behind this
layer, which is exactly what the layer is shaped to allow.

## Design notes

- `CheckResult` separates `measurement` (physical quantity — backends *must* agree; a delta is
  a bug), `penalty` (graded judgement — differences are expected and informative), and
  `verdict` (ours, from our thresholds, never the backend's, so swapping backends cannot
  silently turn a pass into a fail).
- `UNAVAILABLE` / `UNSUPPORTED` are distinct from `FAIL`, and `.ok` is false for all three:
  "we could not measure this" must never read as "this is fine".
- Thresholds default to `complexity.DEFAULT_LIMITS` so the QC layer and the synthesisability
  screen cannot drift apart.
- Native scorers wrap **what `bbl` already does** (`gc_fraction`, `complexity.*`,
  `enzymes.enumerate_cut_sites`) rather than reimplementing Proto's formulas — a mirror would
  agree by construction and measure nothing.
- Scoring is deterministic on both sides (no optimiser, no RNG), so diffs are reproducible.

### Bug found while building this

`enzymes.enumerate_cut_sites` calls `enzyme.search(seq, ...)`, which requires a Biopython
`Seq`. Passed a plain `str` it raises, and the `except Exception: continue` inside
`enumerate_cut_sites` swallows it — silently reporting **zero cut sites**. `bbl.qc.native`
now wraps input in `Seq`. Callers elsewhere pass records so they are unaffected, but the
silent-zero path is worth hardening at source.
