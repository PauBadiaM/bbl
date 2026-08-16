# bbl

reAGENT hackathon — agentic design of optimized plasmids from an inventory of library
plasmids (SnapGene `.dna` files with sequences and annotations).

## First primitive: restriction-based feature deletion

`excise_features()` takes a plasmid and the annotated features to remove, finds two flanking
restriction sites whose ends can actually be ligated, excises the DNA between the cuts, and
re-circularizes.

```python
from bbl import excise_features, write_genbank

plan = excise_features("plasmid/pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8.dna",
                       ["NFKBRE"])

plan.enzyme_pair     # ('MfeI', 'EcoRI')
plan.product         # circular Dseqrecord, 5704 bp, features remapped
plan.deleted_bp      # 66   (the 54 bp feature + 12 bp of unannotated spacer)
plan.junction_seq    # 'ACCGACAATTCT'  -> the MfeI x EcoRI hybrid CAATTC
plan.sites_regenerated   # False: neither site survives the ligation
print(plan.protocol)
write_genbank(plan.product, "pCLM1_designed.gb")
```

Validated against a real cloning step: the experimentalist made **pCLM1** from **pHL391** by
cutting with MfeI + EcoRI and religating. The function proposes that pair independently and
reproduces pCLM1 exactly — same length, same coordinate frame, no feature-coordinate
mismatches.

### What it checks that a naive implementation would not

- **Ends must be ligatable.** Two different enzymes usually leave incompatible overhangs.
  MfeI and EcoRI both leave 5'-`AATT`, which is why this pair works.
- **Nearest is not best.** KpnI cuts closer to `NFKBRE` than MfeI does, but `GTAC` cannot be
  joined to `AATT`. Ranking is ligatable-first, then minimal collateral deletion.
- **Each enzyme must cut exactly once**, or the backbone is fragmented.
- **Cut coordinates, not site positions** — a site can sit outside a feature while its
  staggered cut lands inside.
- **No annotated non-target feature may be damaged**; unannotated spacer may go.
- **It says no when the answer is no** (`NoExcisionFound`), so the agent layer can fall back to
  a PCR-based deletion instead of proposing something that cannot be built.

## Fallback: around-the-horn PCR deletion

When no ligatable pair of sites flanks the target, amplify the backbone outward instead. Two
primers sit back-to-back at the deletion boundaries pointing away from each other, so the
target is never copied; the amplicon is then re-circularized.

```python
from bbl import design_deletion_primers

d = design_deletion_primers("plasmid/test_pHL391_....dna", ["NFKBRE"])   # method="KLD"
d.forward      # del_F [5'-phos]: 5'-TAGGCGTGTACGGTGGGAGGCCTATATAAG-3' (30 nt, Tm 62.1 C)
d.reverse      # del_R [5'-phos]: 5'-TCGGTCAAGCCTTGCCTTGTTGTAGC-3'     (26 nt, Tm 61.6 C)
d.deleted_bp   # 54 -- exactly the feature, no collateral
d.product      # 5698 bp circular plasmid
print(d.protocol)
```

`method="gibson"` instead adds 20 bp of end homology, split 10/10 across the two primer 5'
tails so the shared region straddles the new junction (no phosphorylation needed). Both methods
yield the same product; `overlap` and `overlap_placement` are configurable.

**Restriction is preferred when it is available** — cheaper, and no polymerase errors to
sequence out. `design_deletion_primers` checks for a restriction route and warns when one
exists. The PCR route's advantage is that it is **exact**: on pHL391 it deletes 54 bp where
MfeI + EcoRI removes 66 bp. A decision tree between the two comes later.

It also flags what will actually bite you: primer Tm mismatch, missing G/C clamp, homopolymer
runs, 3'-end self-dimers, amplicons needing long-range polymerase, and mispriming against
repeats — designing across pHL391's 8× BoxB array correctly warns.

## Insertion / subcloning

```python
from bbl import plan_insertion

# case 1: the insert lives in another plasmid
p = plan_insertion(backbone="plasmid/pCLM3_....dna",
                   insert="plasmid/pCLM1_....dna", insert_features=["Lambda BoxB x8"],
                   at=(1491, 1512))
p.strategy        # 'restriction' -- the enzymes could be reused
p.vector.enzymes  # ('XhoI', 'XbaI'), directional
p.product         # 6234 bp == pCLM2, exactly

# case 2: you only have the sequence
p = plan_insertion("plasmid/pCLM3_....dna", "GGGCCCTGAAGAAG...", at=(1491, 1512))
p.insert.order_sequence   # the fragment to order, homology arms already attached
```

`at` is a feature label (replaced), an integer (pure insertion) or a `(start, end)` span.
The rule is two-way — **restriction when the sites are already there, Gibson when they are
not.** `method="auto"` tries restriction reuse first and falls back to Gibson; `method="gibson"`
forces PCR/synthesis. Validated against **pCLM2 = pCLM3 + Lambda BoxB x8**, which it reproduces
base-for-base.

Note the two routes give *different, both-correct* products: restriction brings 6 bp of donor
flank along (6234 bp, what was actually built), Gibson inserts exactly the annotated 294 bp
feature (6228 bp). The plan reports which and why.

`docs/DECISIONS.md` records the reasoning, including several places where the obvious approach
is biologically wrong.

## Inventory index

The library is keyed by **sequence**, not by filename or feature label — both are unreliable
here (one file is named `...Tornado...` but annotates the part as `5'/3' ribozyme`).

```bash
python -m bbl.inventory /path/to/plasmid     # duplicates + lineage report
```

On the current 28-plasmid library: 28 distinct sequences, no duplicates. The strongest
containment relationships it finds are pCLM1 ⊂ pHL391 and pCLM3 ⊂ pCLM2 — precisely the two
cloning steps validated elsewhere in the tests.

## Sourcing

Where each piece of DNA comes from. Both ladders end in "order it", never in failure.

```python
from bbl import scan_inventory, source_insert, source_backbone
entries, _ = scan_inventory("plasmid/")

source_insert(None, entries, name="barcode")   # -> ASK_USER  (+ a question to put to the user)
source_insert(boxb_array, entries)             # -> FROM_INVENTORY (found in 7 plasmids)
source_insert(new_400bp, entries)              # -> SYNTHESIZE (passes the screen)
source_insert(new_repeat_array, entries)       # -> ORDER_PLASMID (refused; nothing to clone from)

source_backbone(ranked_labels, entries)        # construct -> base vector -> order one
```

Provenance is checked **before** complexity, which is what reconciles the two fallbacks: the
BoxB array scores 0.90 and would be refused by any vendor, but it already exists in seven
library plasmids, so it is moved rather than made and never reaches the gate.

`src/bbl/complexity.py` is a **placeholder** for the lab's own scoring function — every decision
made with it carries that caveat in its rationale. Lab state (base vectors, enzyme stock,
unavailable plasmids) lives in `config/lab.json`.

## Cloning reports

A plan says "ligate the vector and the insert". A report says how many microlitres of each.

```python
from bbl import excise_features, load_plasmid
from bbl.report import write_report

parent = load_plasmid("plasmid/pHL391_....dna")
write_report(excise_features(parent, ["NFKBRE"]), "pCLM1_report.html",
             parent=parent, aim="A promoter-only control for the NF-kB sensor.")
```

One self-contained HTML file — no sidecar images, nothing fetched from the network — laid out
like the lab's own notebook entries: aim, cloning strategy, plasmid maps, materials, then one
section per bench step with its reagent table, then fields to fill in as you go.

- **Plasmid maps** via [DnaFeaturesViewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer):
  circular parent and product side by side, plus a linear zoom on the edit, because a 66 bp
  change to a 5.7 kb plasmid is invisible at whole-plasmid scale. Sequence that leaves is red,
  sequence that arrives is green, cut sites and the new junction are marked. Labels are placed
  around the circle at their own feature's angle rather than stacked above it.
- **Live reaction tables.** Type the measured stock concentrations into the digest, Gibson or
  ligation table and the volumes recompute in the page — the water row turns red if the DNA
  alone overflows the reaction. No dependencies; it works from a `file://` URL offline.
- **Enzyme conditions are filled in, from the supplier.** Incubation, thermal inactivation, bp
  from the end and the star-activity-free window come from a dated scrape of Thermo's
  FastDigest table (176 enzymes), matched through isoschizomers — the planner says MfeI, the
  table says MunI. The digest's incubation time and inactivation are derived from it: on
  pHL391 that catches MfeI being un-heat-killable and EcoRI going star beyond 30 min.
- **Reagent maths is calculated, not templated.** The molar-ratio arithmetic is checked
  against a real notebook entry's own spreadsheet to the microlitre (`tests/test_bench.py`).
- **Nothing is invented.** A volume that depends on a Nanodrop reading nobody has taken shows
  the requirement it has to satisfy — "100 ng", "15 fmol" — beside an empty box. Pass
  `concentrations={...}`, or let the agent ask: `report_inputs` says which stocks it needs.
- **The verified plan is quoted, not regenerated.** `plan.protocol` appears verbatim above the
  expanded procedure; where the report overrides it, it says so and why.
- Two numbers move when the product carries a tandem repeat array — outgrowth drops to 30 °C
  and PCR extension goes to 45 s/kb. The array is detected structurally, not by matching "x8".

Lab-specific kit names, buffer volumes and molar excesses live under the `bench` key of
`config/lab.json`; the defaults are seeded from the eCLM24 entry in `plasmid/`.

Maps need the optional extras — on this cluster, wheels only:

```bash
pip install --only-binary=:all: -e ".[report]"
```

Without them the report still builds, minus the figures, and says what to install.

`examples/` holds one generated report per route — open
`examples/sensor_control_report.html` in a browser to see the output. Regenerate with
`python examples/make_reports.py`; that script drives the same `DesignSession` the model-facing
tools call, so it runs without an API key.

## Interactive design sessions

```bash
ant auth login                             # short-lived session token (preferred)
python -m bbl.llm                          # asks where your plasmids are
python -m bbl.llm --plasmids ~/constructs  # a directory, anywhere
python -m bbl.llm --plasmids a.dna b.dna   # specific files
```

Credentials: an **active Claude session** first (`ANTHROPIC_AUTH_TOKEN`, or a token minted from
your `ant auth login` profile), an **API key** as fallback, and otherwise an offer to log in.
The session deliberately outranks the key so a stale exported key can't silently shadow it —
and the startup line prints which source was used. Session tokens are re-minted automatically
if they expire mid-conversation.

```
› I need a control for my sensor: pcDNA3.1, miniCMV-mCherry-LambdaBoxBx8-polyA, AmpR.
  → search_inventory(query=NFKBRE)
  → plan_deletion(plasmid=pHL391, features=['NFKBRE'])
  → compare_product(product_id=prod_1, target=pCLM1)

  [recommendation, then the verified protocol appended verbatim]

› export prod_1 as pCLM1_designed.gb
  ⚠ write 5704 bp to pCLM1_designed.gb [y/N]

› write me a report for it
  → report_inputs(product_id=prod_1)
  I need the pHL391 miniprep concentration to fill in the digest volumes — what did it read?

› 1364 ng/uL
  → generate_report(product_id=prod_1, path=pCLM1_report.html, concentrations={...})
  ⚠ write a cloning report to pCLM1_report.html [y/N]
```

`/report N [path]` does the same thing without going through the model. A comparison run
earlier in the session is remembered, so the report's header can state that the design is
identical to pCLM1 without the model having to carry the verdict back in.

Nine tools over the primitives; state persists across turns, so routes can be planned,
compared and exported by handle. Two properties make it trustworthy rather than plausible:

- **No tool ever returns a sequence.** Products are handles (`prod_1`); the model physically
  cannot construct or paraphrase DNA. Enforced by a test that scans every tool result.
- **The harness renders the protocol, not the model** — enzyme names and bp counts come
  straight from the verified plan.

`session.py` has no Anthropic dependency, so the whole model-facing contract is tested offline
with no API key.

## Where plasmids come from

Nothing assumes a filesystem, a directory, or that plasmids are files. A source answers three
questions — what have you got, what is it called, give me the record:

```python
from bbl import DirectorySource, FileSource, RecordSource, CallableSource

DirectorySource("~/constructs")          # a folder, anywhere
FileSource(["a.dna", "b.dna"])           # uploads, or a hand-picked subset
RecordSource({"pHL391": record})         # already in memory, no disk
CallableSource(                          # anything remote — this is the Benchling shape
    name="benchling:my-project",
    ids=lambda: [s["id"] for s in client.dna_sequences.list(folder_id=FOLDER)],
    load=lambda i: to_record(client.dna_sequences.get(i)),
    label=lambda i: names[i],
)
```

Every public function accepts a directory, a file list, a mapping, or a source —
`as_source()` coerces. Remote fetches are cached, so a large library costs one call per plasmid
actually used. The pHL391 → pCLM1 ground truth is tested against a synthetic remote store with
opaque identifiers and no directory in the call path.

## Install

```bash
pip install -e ".[test]"
pytest -q
```

`pydivsufsort` is pinned to `0.0.18`; newer versions are sdist-only and fail to build on this
cluster. Tests locate the inventory via `BBL_PLASMID_DIR`.

## Layout

```
src/bbl/plasmid_io.py   load/write .dna and GenBank, replace_span, circular comparison
src/bbl/inventory.py    sequence-keyed library index: fingerprints, duplicates, lineage
src/bbl/sourcing.py     insert + backbone sourcing ladders
src/bbl/complexity.py   synthesisability screen (PLACEHOLDER)
src/bbl/config.py       lab state: base vectors, enzyme stock, availability
src/bbl/sources.py      storage abstraction: directory / files / memory / remote
src/bbl/llm/            LLM interface: session boundary, tools, prompts, REPL
src/bbl/report/         bench-ready reports: maps.py figures, bench.py reagent maths,
                        fastdigest.py supplier conditions, build.py plan -> notebook entry,
                        html.py single-file output
src/bbl/targets.py      target resolution + protected-feature classification (shared)
src/bbl/enzymes.py      site enumeration, cut coordinates, end compatibility
src/bbl/excise.py       excise_features() / plan_excisions()   -- deletion by restriction
src/bbl/pcr.py          design_deletion_primers()              -- deletion by inverse PCR
src/bbl/insert.py       plan_insertion()                       -- insertion / subcloning
tests/test_excise.py    pHL391 -> pCLM1 regression
tests/test_synthetic.py 3' overhangs, single-enzyme, blunt, no-route paths
tests/test_pcr.py       primer design on the site-free pHL391 variant
tests/test_insert.py    pCLM3 + BoxB -> pCLM2 regression
tests/test_inventory.py fingerprint invariants, dedup, lineage recovery
tests/test_sourcing.py  both ladders, provenance, complexity placeholder
tests/test_llm.py       tool contract, gates, protocol appending (offline)
tests/test_bench.py     reagent maths against the lab's own notebook spreadsheet
tests/test_report.py    figure pruning, route coverage, nothing-invented invariants
tests/test_sources.py   storage independence, incl. a synthetic Benchling store
docs/DECISIONS.md       design rationale and course corrections
```
