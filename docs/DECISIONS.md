# Decision log

Reasoning behind the design of `bbl`, especially the places where the obvious or originally
intended approach turned out to be **biologically or logically wrong**. Entries marked **↺**
are course corrections — read those first if you are picking this up cold.

## First milestone: `excise_features`

Delete annotated features from a plasmid by cutting at two flanking restriction sites and
religating. Validated against a real cloning step:

| | file |
|---|---|
| parent | `plasmid/pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8.dna` (5770 bp) |
| product | `plasmid/pCLM1_pcDNA3.1_miniCMV-mCherry-LambdaBoxBx8.dna` (5704 bp) |
| what was done | cut with MfeI + EcoRI, ligate — deletes the NF-κB response element to leave a promoter-only control |

`excise_features(pHL391, ["NFKBRE"])` independently proposes **MfeI + EcoRI** and produces a
5704 bp plasmid that matches pCLM1 exactly, in the same coordinate frame, with zero
feature-coordinate mismatches.

---

## Scope and setup

**D1. Repo.** Work inside `github.com/PauBadiaM/bbl`, which was empty at clone time
(`e7b907d`), so the layout is ours to define: `src/bbl/`, `tests/`, `docs/`.

**D2. Execution environment. ↺** Initially declined to run Python at all, because the Sherlock
policy forbids it on login nodes. This is an **interactive session (10 cores, 20 GB)**, so local
execution and `pip install` are fine. Don't re-litigate each session.

**D3. Dependency pin. ↺** `pip install pydna` **fails** here: pip resolves `pydivsufsort`
0.0.20, which is sdist-only and whose `build.sh` exits 127. Pin **`pydivsufsort==0.0.18`**,
which ships a `cp312` manylinux wheel. Recorded in `pyproject.toml`.

**D4. Output format. ↺** We wanted to emit `.dna` to match the inventory, but Biopython's
SnapGene support is **read-only**. Products are written as GenBank via `write_genbank()`;
SnapGene opens those natively.

**D5. Test fixture. ↺** The original plan was a synthetic plasmid ("assume we already have a
starting plasmid"). Superseded by real ground truth (pHL391 → pCLM1). Synthetic plasmids were
kept only for paths the real case cannot reach (see D22).

---

## Biology the first-pass spec got wrong

**D6. "Two restriction enzymes" implies *compatible ends*. ↺** The original output contract
(plasmid + two enzymes) omitted the constraint that makes it physically work: two different
enzymes generally leave different overhangs and will not ligate. `ligation_strategy()` returns
`compatible_overhang` / `same_enzyme` / `blunt`, or `None` to reject the pair.
**MfeI (`C^AATTG`) and EcoRI (`G^AATTC`) both leave a 5'-`AATT` overhang** — that is why the
real experiment worked, and it is the load-bearing fact of the whole first test.

**D7. "Pick the closest sites up/downstream" is the wrong primary criterion. ↺** Empirically
decisive here. Around `NFKBRE` (`[172:226]`) the cuts are **MfeI@161, Acc65I@167, KpnI@171,
EcoRI@227**. The *nearest* upstream site is KpnI, so a nearest-site rule returns **KpnI +
EcoRI** — `GTAC` against `AATT`, which cannot be ligated. MfeI is the nearest *AATT-compatible*
site. Ranking is therefore **ligatable first, then minimal collateral**.

**D8. Compatibility is not "same overhang string". ↺** KpnI and Acc65I both leave `GTAC`, but
3' vs 5' — incompatible. Hence `compatible_end()` rather than comparing `ovhgseq`.

**D9. Enzyme uniqueness is a hard filter** — absent from the original spec. Each enzyme must cut
the plasmid exactly once (or exactly twice, both flanking, for the single-enzyme strategy).
Anything else fragments the backbone.

**D10. "Enzymes that cut off extra bases", restated. ↺** The spec framed this as a property of
the enzyme. The real constraint is the **cut coordinate**: a recognition site can sit entirely
outside a feature while the staggered cut lands inside it. Gate on cuts, never on sites.
Refinement found during implementation: the gate is the **top-strand span `[u_cut, d_cut)`**,
which is provably what the product loses — not a two-strand footprint. `CutSite.footprint` is
still reported (it matters for blunting/Klenow work) but does not gate.
Note this also means KpnI above is rejected *purely* on incompatibility, not on footprint.

**D11. Circularity makes "remove the sequence in between" ambiguous. ↺** Two cuts on a circle
give two arcs and the spec never said which is discarded. We keep the arc containing the
protected features and re-circularize; that also guards against deleting the origin. Targets
that straddle the origin are handled by `_tight_arc()` (complement of the largest gap) plus a
rotation via `shifted()`.

**D12. What is deleted is cut-to-cut, not feature-to-feature.** Part of each recognition site
goes, and the half-sites fuse. MfeI × EcoRI → `CAATTC`, regenerating **neither** site.
Independent evidence: pHL391 carries an `MfeI` annotation and pCLM1 does not. Reported as
`junction_seq` / `sites_regenerated`.

**D13. Collateral deletion is asymmetric.** Removing extra *unannotated* spacer is unavoidable —
a strict no-collateral rule would make the real MfeI/EcoRI answer unreachable (it takes 12 bp
beyond the 54 bp feature). Removing an *annotated, non-targeted* feature is a failure.

**D14. Single enzyme cutting twice is allowed but ranked below a compatible pair.** Cheapest
digest, guaranteed-compatible ends, but religation regenerates the site and the excised fragment
can re-insert — so the emitted protocol calls for backbone gel-purification and rSAP/CIP.

**D15. Enzyme pool is `CommOnly` for now.** `build_pool()` also accepts an explicit name list,
which is the hook for the lab stock list. Stock-list *rules* deliberately deferred.

**D16. Restriction excision is not a general method.** It works only when convenient unique
sites happen to flank the target. `NoExcisionFound` carries a `reasons` tally so the agent layer
can fall back to a PCR route (around-the-horn / Gibson) instead of guessing. Confirmed reachable:
`AmpR` in pHL391 has no route, because ori and the AmpR promoter hem it in.

---

## Implementation choices that guard against our own errors

**D17. Simulate, then verify — do not trust coordinate arithmetic.** pydna performs the digest
and religation (`.looped()` raises on incompatible ends, a free check on the ranking), then
`_verify()` re-checks length, protected-feature integrity, and junction disruption. Failures
downgrade the candidate rather than returning a wrong answer. This is deliberate insurance
against off-by-one and `ovhg`-sign mistakes in the most error-prone part of the code.

**D18. Verify by junction disruption, not by target-sequence absence. ↺** The first
implementation asserted the target sequence was gone from the product. That is **unsound when
the sequence recurs elsewhere** — and plasmids genuinely have repeats; pHL391 carries an 8×
BoxB aptamer array, so deleting one repeat would have looked like a failure. `_verify()` now
checks that the parent's sequence *across each cut* is no longer contiguous. Regression test:
`test_deleting_a_tandem_repeat_is_not_confused_by_duplicates`.

**D19. Rebuild the product in the parent's coordinate frame. ↺** pydna's `cut()` **drops**
features that straddle a cut and returns the arc rotated to the downstream cut. SnapGene instead
*truncates* such features. `_reannotate()` reconstructs the product as `parent[:u] + parent[d:]`
and remaps features, cross-checking the result against the pydna simulation. This is what makes
the output match pCLM1's coordinates exactly rather than merely up to rotation.

**D20. `primer_bind` features are not protected, and umbrella features do not block. ↺** This
would have broken the ground-truth test outright. pHL391 annotates primers such as
`oCLM14_mCherry_Fragment_FOR` `[157:173]` that overlap the deletion, and an `Insert Sequence`
`[160:1327]` that *contains* the target. Treating either as protected rejects MfeI + EcoRI and
the correct answer becomes unreachable. Rules: skip `NON_FUNCTIONAL_TYPES`
(`primer_bind`/`primer`/`source`) entirely; a feature that *contains* the target is an umbrella
that gets truncated with a warning; only a feature *partially* overlapping the deleted span
blocks. Clipped primers are dropped from the product rather than left mis-annotated.

**D21. Deduplicate equischizomers.** MunI is an equischizomer of MfeI, so without dedup the
winning plan reappears as its own "alternative". `canonical_name()` picks a preferred vendor's
name (NEB first), and candidates are deduped on canonical pair + cut positions.

**D22. Circular comparison must be rotation- and reverse-complement-invariant.** Records have an
arbitrary origin, so string equality would fail spuriously — `circular_equal()` /
`rotation_offset()`.

**D23. Synthetic plasmids cover what the real case cannot.** pHL391 only exercises a
5'-overhang two-enzyme deletion. `tests/test_synthetic.py` forces the 3'-overhang single-enzyme
path (PstI, where the bottom-strand break precedes the top-strand break), blunt ligation,
incompatible-pair rejection, the essential-element warning, and the no-route path. Those tests
pass an explicit `enzyme_pool` so the expected answer is deterministic.

---

## Known limitations of the restriction route (reported as warnings, not modelled)

Dam/Dcm methylation blocking sites; star activity; partial digests; double-digest buffer and
temperature compatibility (Biopython carries no buffer data); frameshift risk if a cut lands in
a retained CDS.

---

# Second milestone: `design_deletion_primers` (around-the-horn PCR)

The fallback for when no ligatable pair of restriction sites exists. Two primers sit
back-to-back at the deletion boundaries pointing *outward*, so the backbone is amplified and
the target is simply never copied; the linear amplicon is then re-circularized.

Test case: `plasmid/test_pHL391_...dna` is pHL391 with the MfeI **and** EcoRI sites removed
(5752 bp vs 5770 bp, zero cuts for either enzyme), so `excise_features` correctly raises
`NoExcisionFound` and the PCR route is the only option.

**D24. The 5' end of each primer is pinned, not the 3' end.** After ligation the two primer
5' ends become the new junction, so pinning them is what makes the deletion exact. The only
free parameter per primer is length, which moves the 3' end — conveniently also the end whose
G/C clamp and uniqueness we want to optimise.

**D25. Exactness is the reason this route is worth having.** The boundaries come from where the
primers sit, not from where an enzyme happens to cut. On the original pHL391 the PCR route
deletes exactly 54 bp (the annotated feature) where MfeI + EcoRI removes 66 bp — 12 bp of
collateral. Restriction is still preferred when available (cheaper, no polymerase errors, no
sequencing of the whole amplicon), so `design_deletion_primers` checks for a restriction route
and warns when one exists. The full decision tree is deferred.

**D26. Two re-circularization methods.** `KLD` (default) blunt self-ligates 5'-phosphorylated
ends — order phosphorylated primers or treat with T4 PNK. `gibson` instead puts a 5' homology
tail on the forward primer, duplicating the `overlap` bases immediately upstream of the target
so the amplicon's two ends are homologous; no phosphorylation needed. Both give an identical
product, which is asserted in the tests.

**D27. DpnI is not optional.** The parental plasmid is Dam-methylated from a `dam+` host and the
PCR product is not, so template carry-over would otherwise dominate the transformants. In the
emitted protocol always.

**D28. Doubling a circular template double-counts every primer site. ↺** The first uniqueness
check searched `template * 2` to catch matches spanning the origin — but that reports every
single-copy site **twice**, so every primer looked like a mispriming risk. Fixed by extending
the template by only `len(probe) - 1` bases. Regression test:
`test_binding_sites_does_not_double_count`.

**D29. Slicing a circular `Dseq` with `seq[:0]` returns the whole circle. ↺** pydna's circular
slice convention means `start == stop` is a full rotation, not an empty string. Any deletion
whose upstream boundary was position 0 silently produced a **double-length** plasmid. Both
`delete_span()` and the excision junction now slice plain strings via `slice_circular()`.
This was latent in the restriction route too, not just the PCR one. Regression test:
`test_deletion_starting_at_position_zero`.

**D30. Rotating for an origin-spanning target must remap coordinate targets. ↺** `resolve_target`
rotates the record when the target wraps the origin. Feature-name targets survive this because
`shifted()` remaps annotations, but *coordinate* targets do not — the same wrapping span was
re-derived on the rotated record and it recursed until the stack blew. Coordinates are now
remapped explicitly, with a one-shot guard.

**D31. Primer QC screens are heuristic and say so.** Tm is nearest-neighbour at default salt
(`Bio.SeqUtils.MeltingTemp`); G/C clamp, homopolymer runs, 3'-end self-dimer and hairpin use
simple string checks. Good enough to rank lengths and reject bad oligos, not a substitute for
Primer3 or IDT OligoAnalyzer — the emitted design says as much. Mispriming *is* checked
properly against the actual template, which matters: designing across pHL391's 8× BoxB array
correctly warns (`test_mispriming_in_a_repeat_is_flagged`).

# Third milestone: `plan_insertion` (insertion / subcloning)

Two cases, as posed: (1) the insert comes from another plasmid — reuse the restriction enzymes
if possible, else PCR it out and Gibson it in; (2) you only have the sequence — attach homology
arms and Gibson it in.

Ground truth: **pCLM2 = pCLM3 with the Lambda BoxB x8 array swapped into a 21 bp stub at
[1491:1512]**; the two are byte-identical either side, and the array is identical in pCLM1 /
pCLM2 / pHL391 so pCLM1 serves as donor. `plan_insertion` independently proposes
**XhoI + XbaI** and reproduces pCLM2 exactly (6234 bp).

**D32. Vector and donor have asymmetric constraints. ↺** The natural instinct is to require
both plasmids to be cut exactly once by each enzyme. Wrong: the **vector** must be, or the
backbone is destroyed, but the **donor** may be cut anywhere else — you gel-purify the insert
band. Only a cut *inside* the released fragment disqualifies it. The symmetric rule would have
rejected the real XhoI + XbaI subcloning.

**D33. Restriction reuse is not as rare as expected.** The brief predicted it would be "very
rare". Within a curated library it is not: sibling plasmids descend from a shared backbone
(pcDNA3.1 here) and therefore share an MCS. Worth checking first every time.

**D34. Directionality is decided by the *vector's* two ends.** If they are mutually
incompatible the insert can only go in one way and the vector cannot self-close. If they are
compatible, the protocol must add rSAP/CIP and orientation screening. XhoI × XbaI are
incompatible, so this case is directional.

**D35. Same-enzyme junctions regenerate the site**, unlike the MfeI × EcoRI deletion junction
(D12). Reported as `sites_regenerated`, which also tells you a diagnostic digest will be
positive.

**D36. Restriction and Gibson give genuinely different products here, and both are right.**
Restriction brings 6 bp of donor flank along (cuts sit outside the annotated array), giving
6234 bp = pCLM2. Gibson inserts exactly the 294 bp feature, giving 6228 bp. The restriction
result is what the experimentalist actually made; the Gibson result is cleaner. Neither is a
bug — the plan reports `inserted_bp` and warns about the donor flank.

**D37. Choosing the nearest compatible donor site is wrong. ↺** The first implementation took
the closest cut among all compatible partners. Compatible sets include promiscuous cutters —
BfaI (`C^TAG`) is XbaI-compatible and cuts constantly — which were then found to cut inside the
fragment, and *all 900* candidate pairs were rejected. Donor partners are now ranked
same-enzyme first, then by how few times they cut the donor, then by proximity, and each
candidate is validated against the internal-cut constraint rather than committed to blindly.

**D38. Alphabetical tie-breaking picks the wrong enzyme. ↺** Several enzymes cut the same
position: XhoI, PaeR7I, SlaI, Sfr274I (equischizomers) and PspXI (a degenerate 8-bp site).
Ranking by vendor alone left an alphabetical tie that chose `PaeR7I`, then `PspXI` — both
technically correct and both names nobody uses. `supplier_count()` now breaks the tie by how
many vendors stock the enzyme (XhoI: 9, PaeR7I: 1), applied to both naming and candidate
scoring. MfeI still beats MunI because the vendor-preference key dominates.

**D39. SnapGene name fields contain junk. ↺** Records were carrying
`http://www.genscript.com/ORIGDB|GenBank` as their name, which then appeared in emitted
protocols ("digest http://www.genscript.com/..."). `load_plasmid` now falls back to the file
stem when the parsed name is not a plain identifier.

**D40. Insertion is modelled as replacement of a span.** `at` is a feature label (replace it),
an integer (pure insertion), or a `(start, end)` span. This makes `replace_span()` the single
primitive underneath deletion *and* insertion, and it means a pure insertion needs asymmetric
tie-breaking on feature coordinates: a feature ending exactly at the insertion point must stay
put while one beginning there must shift.

## Shared between routes

`targets.py` holds target resolution and feature classification so that "delete NFKBRE" means
the same thing regardless of how it is built; `plasmid_io.delete_span()` performs the sequence
surgery and feature remapping for both. The restriction route additionally cross-checks its
`delete_span` result against the pydna digest simulation (D17).


---

# Fourth milestone: `inventory.py` (sequence-identity layer)

**D41. The library is keyed by sequence, not by name.** Filenames and feature labels are both
unreliable in this inventory, demonstrated repeatedly:

- filename says `LambdaBoxBx8`, the feature label says `Lambda BoxB x8`; filename says
  `NFKBRE1`, the feature is `NFKBRE`
- `pHL511_..._Tornado-...` has **no feature named Tornado** -- it is annotated as
  `5' ribozyme` / `3' ribozyme`. A label-driven lookup fails outright here.
- labels are non-unique: 8x `BoxB RNA aptamer`, 2x `NFKB_Canonical1`, and `Insert Sequence`
  appears four times in pCLM3 at unrelated spans
- the same part carries the same label in three files but is only *provably* one part by
  sequence (the BoxB array is byte-identical in pCLM1/pCLM2/pHL391)

Layering, in decreasing order of trust: **sequence = truth, feature labels = annotation,
filename = declared intent**. Names resolve *what the user means* and nothing more; every
proposal is confirmed against sequence. Filenames remain valuable as the spec-mode vocabulary
(`pID_backbone_part-part-part` is exactly how the target gets requested).

**D42. Identity is rotation- and orientation-normalised.** `fingerprint()` takes the
lexicographically least rotation (Booth's algorithm, O(n)) of the sequence and its reverse
complement, whichever is smaller, then hashes. Without this the same construct saved with a
different origin looks like a different plasmid.

**D43. The duplicate hypothesis was wrong. ↺** Two filenames carry sync-conflict artifacts
(`pFH2.22 ... (1)`, `pFH2.48 ... (fig_s conflicted copy 2022-07-18)`), which suggested the
library held redundant copies. It does not: **28 files, 28 distinct sequences, zero exact
duplicates.** Those are simply the only copies of those constructs under untidy names. Nothing
to prune -- but the check is kept as a test, because the assumption was reasonable and will be
made again.

**D44. Containment recovers the cloning lineage, and it validates itself.** Canonical k-mer
containment ranks pCLM1 subset-of pHL391 (0.995) and pCLM3 subset-of pCLM2 (0.992) as the #1
and #3 relationships in the whole library -- the two cloning steps independently validated in
`test_excise.py` and `test_insert.py`. That is strong evidence containment is the right
backbone-ranking signal for the planner. It also surfaced pairs we had not noticed:
pHL335/pHL394 (0.983), pCLM21/pCLM22, pCLM23/pCLM24, pHL405/pHL406.

**D45. Containment has a floor set by the shared backbone.** Every plasmid here descends from
pcDNA3.1 or pePB, so unrelated pairs still score ~0.5 and the low range is noise. Only
containment above ~0.95 is informative; `lineage()` defaults to 0.98. The `Relationship.kind`
thresholds are provisional -- the raw `containment`/`jaccard` numbers are what should be
consumed.


---

# Fifth milestone: sourcing ladders (`sourcing.py`, `complexity.py`, `config.py`)

Both ladders terminate in "order it" rather than in failure.

**D46. A failed complexity check means order the whole plasmid -- but only for genuinely new
sequence.** The lab loop has two apparently conflicting fallbacks: an unsynthesisable fragment
should send us to intermediate cloning from the database (earlier), or to ordering the finished
plasmid (later). They reconcile on **whether the sequence already exists in the library**:

- *in the library* -> it is moved, not made, so it never reaches the complexity gate at all
- *genuinely new and refused* -> nothing to clone from, no workaround, order the plasmid

The ladder therefore checks provenance **before** complexity, not after. Concretely: the 8x
BoxB array scores 0.90 and is refused by the screen, yet is never a problem, because it exists
in seven library plasmids.

**D47. Provenance lookup must be sequence-based, and it pays immediately.** Searching for the
294 bp BoxB array by sequence finds it in **seven** plasmids -- pCLM1, pCLM2, pCLM22, pCLM24,
pHL229, pHL391, test_pHL391. We only knew of three from labels. `find_sequence()` searches both
strands and wraps the origin, and reports reverse-strand hits so orientation can be checked.

**D48. Complexity screening is an explicit placeholder** (`complexity.py`), to be replaced by
the lab's scoring function. `ComplexityReport` is the contract the planner branches on:
`synthesizable` (bool), `score`, `reasons`, `metrics`, `vendor`, and `is_placeholder`, which is
surfaced in the rationale so a decision made on placeholder logic is never mistaken for a real
vendor verdict. The heuristics are crude but not vacuous -- length, GC, homopolymer, longest
internal repeat and repeat fraction -- chosen so the one failure mode that matters here (tandem
arrays) is caught.

**D49. Feed the screen the fragment as ordered, not the bare insert.** Homology arms are part of
what the vendor synthesises, so `plan_insertion(...).insert.order_sequence` is the right input.

**D50. Backbone selection is a three-rung ladder, not two.** Closest existing *construct* ->
config-tagged empty *base vector* -> *order one*. Rung 2 exists because pHL162_pcDNA3.1_MCS is
an empty pcDNA3.1 and is always a valid, if expensive, answer for that vector family.

**D51. "Order a backbone" may block planning entirely.** Cut sites and homology arms cannot be
placed against a vector name. If a reference sequence can be obtained, planning continues with
the backbone marked as procured; otherwise the output is honestly two-phase: order, then re-plan
on arrival. `BackboneDecision.question` asks which.

**D52. Lab state is config, not code** (`config/lab.json`): `base_vectors`, `enzyme_stock`,
`unavailable`. Whether a plasmid is a "blank" is a judgement about intent that cannot be
inferred from sequence, and a `.dna` file is a design archive rather than proof of a tube --
hence an explicit `unavailable` list.

**D53. Four human-in-the-loop gates.** Everything between them is deterministic: (1) spec
ambiguity -- are unlisted features must-be-absent or don't-care; (2) a needed sequence is not
known -- ask, or search Addgene and offer candidates; (3) confirming an external candidate;
(4) choosing among ranked routes.


---

# Sixth milestone: the pCLM2 case (second human-verified ground truth)

Reported route: pHL391 cut **MfeI + NheI** (removing NFKBRE1 + miniCMV, 132 bp); CMV amplified
from **pHL162** with primers adding matching restriction sites; ligate. Product pCLM2, 6234 bp.
Every step was reconstructed and verified against the files.

**D54. I mislabelled a hypothetical route as history. ↺** I had described "pCLM2 = pCLM3 + the
BoxB array" as one of two validated *cloning steps*. The sequence relationship is exact and the
test remains a valid test of the insertion machinery -- but that is **not how pCLM2 was made**.
pCLM3 (file dated 8 Dec 2025) postdates pCLM2 (20 Nov 2025) and is almost certainly the
derivative, made by deleting BoxB. Sequence containment reveals *possible* routes, never
history. Any claim about what was actually done needs external evidence.

**D55. `lineage()` had the direction backwards, in both known cases. ↺** It returned
``(parent, child) = (smaller, larger)``, assuming derivation adds sequence. Both real
derivations here are **deletions**, so the derivative is the smaller member: pCLM1 from pHL391,
pCLM3 from pCLM2. Fixed by removing the directional claim entirely -- pairs are returned
smaller-first for stable ordering only, and direction requires file dates or lab records.

**D56. A fourth assembly method was missing: PCR with restriction-site primer tails.** The
insert is amplified with 5' tails carrying the recognition sites, the amplicon is digested, and
the fragment is ligated into the matching vector ends -- no Gibson. Neither existing route
reproduces pCLM2: pure subcloning gives 6372 bp (it drags 150 bp of extra pHL162 flank along,
because the donor sites sit well outside the CMV block) and Gibson gives 6222 bp. Implemented
as ``method="pcr_restriction"``, which reproduces pCLM2 **exactly**.

**D57. Primer tails need a clamp.** Most enzymes cut poorly at the very end of a fragment, so
the tail is ``clamp + site + insert``. The clamp is removed by the digest and never appears in
the product, so only its length matters; ``DEFAULT_CLAMP`` is 6 bp.

**D58. Build the fragment by simulating the digest on the amplicon string.** Rather than deriving
which bases of each site survive on which side, ``_pcr_restriction_fragment`` assembles the
amplicon and runs the enzymes over it linearly. This also gives the "does an enzyme cut inside
the insert" check for free -- a site inside means more than one cut, and the route is refused.

**D59. Experimentalists add spare sites, and it shows up in the sequence.** The real insert is
596 bp where the method alone predicts 590: the forward primer carried an extra **KpnI**
(``GGTACC``) site after the MfeI site, presumably for future cloning. With
``extra_upstream="GGTACC"`` the product is byte-identical to pCLM2; without it, it is 6 bp
shorter and functionally the same. Exposed as ``extra_upstream`` / ``extra_downstream`` rather
than guessed at.

**D60. The tool independently chose the right backbone enzymes.** Given pHL391 and the
insertion site, it proposes MfeI + NheI unprompted -- the pair actually used. And the insert
boundaries used at the bench, pHL162[234:818], are exactly the annotated ``CMV enhancer`` +
``CMV promoter`` span, so feature-level addressing would have found them.


---

# D61. The RE-tail insertion variant was built, then removed. ↺

`method="pcr_restriction"` -- PCR the insert with recognition sites in the primer tails, digest
the amplicon, ligate into matching vector ends -- reproduced pCLM2 byte-for-byte and was then
**removed at the user's direction**. Retaining restriction sites at the junctions is good
preparation for future cloning, but for the step in front of you it costs an extra digest,
requires the sites to be absent from the insert, and requires clamp bases outside each site --
expertise and failure modes spent on a plasmid that is functionally identical to the Gibson
product.

The rule is now two-way and stateable in one line: **restriction when the sites are already
there, Gibson when they are not.** No hybrid. This also collapses the duplication between the
two tail-based branches that was going to drift.

What this costs: the pCLM2 regression no longer asserts byte-identity. It asserts the route we
would actually propose -- Gibson, 6222 bp, seamless junctions -- and pins the 12 bp difference
against the historical construct, which is two vestigial MfeI/NheI sites plus a spare KpnI the
experimentalist parked for later. D56-D59 record how the removed variant worked, should a future
case genuinely need retained sites.

**Correction to D56, twice overstated.** I described the RE-tail variant as "a fourth assembly
method was missing" and later as a route that "wasn't in the toolkit". Both wrong. PCR
amplification of an insert with 5' primer tails already existed and was tested (`gibson_pcr`);
the addition was a different tail string plus one digest, ~60 lines on shared machinery. And
"neither existing method matched" referred to byte-identity with the historical file, not to
capability -- `gibson_pcr` would have produced a perfectly good pCLM2 all along.


---

# Seventh milestone: the LLM interface (`src/bbl/llm/`)

Interactive design sessions over the existing primitives. `python -m bbl.llm`.

**D62. Tool Runner, not the Claude Agent SDK.** They sound alike and are different packages.
The **Agent SDK** (`claude-agent-sdk`) is Claude Code as a library — built-in Read/Write/Edit/
Bash, for filesystem agents. The **Tool Runner** (`client.beta.messages.tool_runner`, in the
regular `anthropic` SDK) drives the loop over tools *you* define. We want the model calling
`excise_features`, not editing files, so Tool Runner is the fit. Managed Agents was also
rejected: it hosts the sandbox on Anthropic's infrastructure, and our tools need the `.dna`
files and biopython that live on OAK.

**D63. Handles, not payloads — the invariant made structural.** No tool returns a sequence.
Products live in a session-side registry keyed by `prod_1`, and the model passes handles. This
keeps ~1400 tokens of ACGT out of context per call, and — more importantly — means the model
*cannot* construct, paraphrase, or corrupt DNA even if it tried. "Only `replace_span` builds
sequence" stops being a rule and becomes a property of the interface. Tested by scanning every
tool result for runs of 40+ ACGT.

**D64. The harness renders the protocol, not the model.** `.protocol` already contains exact
enzyme names, bp counts and junctions. A paraphrase could turn MfeI into MunI or 6234 into
6324 — errors no reader would catch. The model writes the judgment; `design()` appends the
protocol verbatim for any product created during that turn, and the system prompt tells the
model not to restate it.

**D65. Refusals are successful tool results.** `{"feasible": false, "reasons": {...}}`, never
`is_error`. Marking a domain refusal as an error makes the model retry the identical call.
`plan_deletion(method="auto")` encodes the preference order in code — restriction, then PCR —
so the model relays a ranking rather than choosing chemistry.

**D66. Three of the four LLM jobs are single calls, not loops.** Parsing a request into a
`Spec` and resolving loose part names are extraction tasks — `client.messages.parse()` with a
Pydantic model beats a tool loop. `Spec.ambiguities` is the load-bearing field: the model
*declares* what it would otherwise have to guess (notably whether unlisted features must be
absent), and the harness turns that into the clarifying question. The schema says "Do not
infer" in the `must_not_contain` description.

**D67. Plasmid ID is its own resolution tier. ↺** The first resolver matched exact → prefix →
substring, which made `pCLM2` ambiguous against pCLM21/22/23/24 — exactly the name a user or
model types. Fixed by matching the label up to the first underscore as a tier above prefix.

**D68. Under-described tools, caught by a test. ↺** `plan_deletion` and `plan_insertion`
originally said what they did but never *when to call them*, and `export_product` was two
sentences. Since description quality is the main driver of tool selection, the fix was longer
descriptions, not a lower threshold. `test_tools_build_with_usable_schemas` enforces a floor.

**D69. Gating is a session callback, not loop surgery.** `DesignSession.confirm` is consulted
inside `export_product` — the only outward-facing tool. A decline returns
`{"declined": true, "note": "do not retry, ask what they want instead"}` so the model asks
rather than looping. Wiring it to `input()` in the REPL and to a lambda in tests keeps the
whole gate testable offline.

**D70. The session layer has no Anthropic dependency.** `session.py` and `schemas.py` import
cleanly without the SDK; only `tools.py` and `chat.py` need it. The entire model-facing
contract — every tool's JSON, both ground-truth cases, the gates, protocol appending — is
tested with no API key and no network.


---

# Eighth milestone: storage independence (`src/bbl/sources.py`)

**D71. Nothing may assume a filesystem, a directory, or that plasmids are files. ↺** The
inventory globbed a directory and `InventoryEntry` carried a `Path` that eight call sites
reopened. That worked only because the repo happened to sit next to `/oak/.../plasmid` -- an
accident of this machine, not a property of the design. Real users will point at a local
folder, a mounted share, a scratch directory on whatever cluster they are on, or eventually
Benchling.

A source answers three questions -- what have you got (`ids`), what is it called (`label`),
give me the record (`load`) -- and everything is built on that:

| source | for |
|---|---|
| `DirectorySource` | a folder of `.dna`/`.gb`, anywhere |
| `FileSource` | an explicit list -- uploads, a hand-picked subset |
| `RecordSource` | records already in memory; no filesystem at all |
| `CallableSource` | anything remote: an ids function plus a load function |
| `MultiSource` | several libraries as one |

`as_source()` coerces a directory, a file list, a mapping or a source, so every public function
still takes a plain path while the machinery underneath is storage-agnostic.

**D72. Benchling is a `CallableSource`, not a special case.** It needs an ids function and a
load function returning anything `load_plasmid` accepts. No API client is vendored and no stub
pretends to work -- the shape is documented in the module and exercised by
`test_full_design_from_a_synthetic_remote_store`, which runs the pHL391 to pCLM1 ground truth
against opaque `bnch_*` identifiers with no directory in the call path.

**D73. Remote fetches are cached per source.** A large remote library costs one call per
plasmid actually touched, not one per access -- the inventory scan alone reads every record
once, and the session then re-reads several of them.

**D74. `InventoryEntry.path` is now derived, and `None` for remote sources.** `label` became a
real field rather than `path.stem`, and callers use `entry.load()`. Errors name the source
(`benchling:demo`) instead of a directory that may not exist.

**D75. The session asks where the library is.** Resolution order: the `--plasmids` argument,
then `$BBL_PLASMID_DIR`, then `./plasmid` in the working directory, then it prompts -- and
re-prompts if the answer holds no plasmid files. `/library` shows or switches libraries
mid-session. The old default was a repo-relative path, which is exactly the assumption D71
removes.


---

# D76. Active Claude session first, API key as fallback, otherwise offer to log in

Precedence, first match wins:

1. **Session token** -- `ANTHROPIC_AUTH_TOKEN` if exported, else one minted from an active
   `ant auth login` profile. Short-lived and refreshable; nothing durable on disk.
2. **API key** -- `ANTHROPIC_API_KEY`. Works, but sits on NFS, leaks into job logs and `env`
   dumps, and never expires.
3. **Ask** -- offer to run `ant auth login` right there, then retry.

**The session deliberately outranks the key. ↺** The first version had it the other way round,
inheriting the SDK's own order. That reproduces the shadowing trap: a stale exported key
silently wins over a perfectly good session, and nothing says so. Here a live session is always
preferred, and the startup line prints which source was used so the choice is visible.

**The credential pre-check was wrong. ↺** `chat.py` gated on `ANTHROPIC_API_KEY` or
`ANTHROPIC_AUTH_TOKEN` being set -- so a working `ant auth login` profile, which sets neither,
would have been rejected with "no credentials found". Resolution now goes through `resolve()`,
plus a free `models.list(limit=1)` that proves the credential works before the first real turn
rather than failing mid-conversation.

**SDK 0.122 does not resolve an `ant` profile.** Its constructor reads only the two environment
variables; a bare `Anthropic()` leaves both `api_key` and `auth_token` as `None`. The token is
therefore minted explicitly via `ant auth print-credentials --access-token`. The flag is
required -- with no flag the CLI prints the whole credentials JSON, which is not a bearer token
and fails with an opaque protocol error.

**A bearer token needs the `oauth-2025-04-20` beta header** on `/v1/messages`, and the SDK does
not attach it. `_client()` adds it whenever a token rather than a key is in play.

**Only a session is refreshable.** `Credential.refreshable` drives `with_refresh()`: an expired
session token is re-minted once and the call retried, so an hour-long conversation does not die
mid-turn. An API-key failure is deliberately *not* retried -- a credential that cannot expire
failing means something else is wrong, and retrying would hide it.

**Login is interactive and headless-aware.** `login()` runs `ant auth login` with the terminal
inherited, adding `--no-browser` when there is no display -- the usual case over SSH on a
cluster, where it prints a URL to open elsewhere and takes the code back in the terminal. If
the CLI is missing, the error says how to install it.


---

# D77. Reports are a form to work from, not a summary to file

The report layer (`src/bbl/report/`) exists because a plan and a bench protocol are not the
same document. `plan.protocol` says "ligate the vector and the insert"; the person at the
bench needs to know how many microlitres of each, what percentage gel resolves the band, and
what to write down afterwards. The shape is taken from a real lab-notebook export in the
inventory folder (`eCLM24_Cloning_PiggyBac`, June 2026): aims, cloning strategy, then one
section per step with its reagent table, then blanks for the results.

**No design decision is taken in the report layer.** Every enzyme, coordinate and base-pair
count is copied from the plan, and the plan's own `protocol` is reproduced verbatim in a
"Verified plan" box above the expanded procedure. If the two ever disagree, the box is the one
that was checked against the sequence, and the document says so. This is the same split as
D-the-harness-renders-the-protocol in the LLM layer: the part that was verified is quoted, not
regenerated.

**Never invent a concentration. ↺** The first version filled the volume columns using a
nominal 100 ng/µL so the tables would look complete. That is a fabricated measurement sitting
in a column of real ones, and at the bench it is indistinguishable from a Nanodrop reading
somebody actually took. Unknown stocks now print `____` and the *requirement* — "100 ng",
"15 fmol" — is still stated, so the row is a form to fill in. `concentrations={}` supplies
real numbers when they exist and the volumes compute.

**The calculators are checked against the lab's own spreadsheet.** `tests/test_bench.py`
reproduces the eCLM24 entry's Gibson volumes (2.29 / 0.71 / 1.99 µL into 10 µL) and its
ligation backbone volume (1.797 µL) to the microlitre. Two findings came out of that:

- the Gibson sheet uses 650 g/mol per bp and the ligation sheet 660. The 1.5% difference is
  below pipetting error, but each calculator keeps the constant its source used, so the
  numbers reconcile against the sheet they came from rather than *almost* reconciling.
- the ligation sheet says 1:5 and then uses 80 fmol against a 15 fmol backbone, which is
  1:5.33. We take the stated ratio literally (75 fmol). The test records the discrepancy
  rather than quietly matching it.

**Two numbers change because of repeats, and both are derived, not typed.** A tandem array is
detected structurally — a feature containing three or more identically-labelled sub-features,
so it catches an array whose umbrella is not called "x8". When one is present the outgrowth
drops to 30 °C and the PCR extension goes to 45 s/kb. The second of those disagrees with the
generic 30 s/kb in `plan.protocol`, so the report says which figure it is overriding and why;
leaving the reader to notice two different extension times on one page is how they end up
trusting the wrong one.

**What the supplier's table knows, we do not.** FastDigest incubation times and inactivation
temperatures are buffer- and format-specific. The enzyme-conditions table names the enzymes,
leaves the cells blank, and links the Thermo table — the same refusal as the methylation and
star-activity warnings in `excise`.

---

# D78. Drawing a SnapGene record literally buries the figure

`DnaFeaturesViewer` draws what it is given, and a SnapGene record is not a clean feature list:
47 features on pHL391, of which 23 are primer bindings, one is an umbrella spanning a fifth of
the plasmid (`Insert Sequence`, 1167 bp over the whole cassette), five are 10 bp NF-κB
sub-sites inside the 54 bp part being deleted, and eight are identical monomers inside
`Lambda BoxB x8`. Plotted as-is, the label stack is taller than the plasmid.

`prune_features` is the editorial pass. Nesting resolves **in opposite directions depending on
whether the children agree**:

- children with **different** labels → the parent is an umbrella and the children are the
  informative layer, so the umbrella goes (`Insert Sequence`);
- children with **one** label → the parent is an array and it is the informative layer, so the
  children go (`Lambda BoxB x8` keeps its name, the eight monomers disappear).

**A protected span overrides all of it. ↺** The umbrella rule alone deletes `NFKBRE` — it
contains five differently-named sub-sites, so by the rule it is an umbrella. It is also the
feature the report is about. The target's span is passed in as `protect`, survives every pass,
and its children then fall out as nested-inside-a-retained-feature. Inferring the target from
the edit coordinates instead does not work: the smallest feature containing the cut-to-cut
span is `Insert Sequence`, not `NFKBRE`.

**Only features wholly inside the edit are coloured as removed.** A feature merely straddling
the span survives the edit, truncated; colouring it red would misreport it as deleted.

**Figures are inline SVG, not files.** A report is one artefact — nothing to keep together,
nothing to lose. Matplotlib derives clip-path ids from `svg.hashsalt`, so each figure is
rendered under its own salt; without that, four figures in one document can capture each
other's clip paths. A test asserts the referenced ids are disjoint across figures.

**The plotting stack is an optional extra**, so `import bbl` never needs matplotlib. A missing
extra degrades to a report with no maps and a warning saying what to install — it does not
raise. On this cluster it must be installed with `--only-binary=:all:`: glibc is 2.17, current
Pillow and matplotlib ship `manylinux_2_28` wheels, and without the flag pip falls back to
building them from source and fails on a missing libjpeg. Same class of problem as the
`pydivsufsort` pin in D3.


---

# D79. Second pass on the report, from bench feedback

Seven changes after the first draft was read by the person who would use it.

**Fill in the supplier's table rather than linking it. ↺** D77 argued that FastDigest
incubation and inactivation conditions were the supplier's to state and ours to leave blank.
The feedback was blunt: look it up. That was right, and the refusal was the wrong instinct
applied to the wrong thing — these are *published constants*, not judgement calls, and the
cost of the blank was a person opening a browser mid-protocol. `report/fastdigest.py` is now a
dated scrape of Thermo's own table, 176 enzymes, and every rendered table carries the
retrieval date and the link. What survives from D77 is the boundary: an enzyme not in the
table still renders blank rather than plausible.

Doing this surfaced two facts the hand-written defaults had wrong, both on the validated
pHL391 → pCLM1 digest:

- **MfeI cannot be heat-inactivated** (Thermo sells it as MunI; "No — chloroform extraction").
  The protocol previously said "then heat-inactivate" for every digest. It now says to
  column-purify, and says which enzyme is the reason.
- **EcoRI shows star activity beyond 0.5 h**, and the default incubation was "37 °C for 1 h",
  copied from the notebook entry. A FastDigest reaction is *five minutes*; the old default
  spent an hour past the point where EcoRI starts cutting the wrong sites. Incubation is now
  derived per digest — the slowest enzyme sets the time, the least heat-labile sets the
  inactivation.

The name mapping is derived, not maintained: our planner says MfeI and ClaI, FastDigest sells
MunI and Bsu15I, and Biopython's `isoschizomers()` already knows they are the same activity.

**Reaction tables are live. ↺** Every volume is a mass over a ng/µL reading that does not
exist when the report is written, so the first draft printed `____` and left arithmetic for
the reader. Now the stock column is an input and the volumes recompute in the page — 30 lines
of dependency-free JavaScript, because the report has to keep working from a `file://` URL
with no network years from now. Digest, Gibson and ligation collapsed into one `Reaction`
model to make that a single calculator rather than three. The water row goes red when the DNA
alone overflows the reaction, which is a real 8pm mistake.

The `_ = "no scripts"_` property from D77 is therefore gone, replaced by a narrower one the
test now enforces: no `src=`, no `<link>`, no `@import` — nothing is ever *fetched*.

**Ask for the concentrations rather than waiting to be given them.** `report_inputs` lists
the stocks a given design needs and the system prompt tells the model to ask before writing
the report. Supplied numbers fill the tables; junk ("unknown", "") is dropped rather than
coerced, because a zero silently divides.

**Full construct names.** A loaded record keeps a 16-character GenBank LOCUS name, which is an
accession — and pCLM21 through pCLM24 differ only after the underscore. The report now takes
the inventory label and renames a shallow *view* of the record, so the caller's copy is
untouched. Short IDs survive only where a full name would be noise, e.g. gel lane lists.

**Circular labels follow the circle.** DnaFeaturesViewer stacks them in a block above the
plot; on pHL391 that block was taller than the plasmid and impossible to match to a feature.
Labels are now placed by hand at each feature's own angle, reading outward, with crowded
neighbours pushed to an outer ring — so the map uses the bottom of the circle as well as the
top, the way SnapGene's does. Capped at 14 labels plus whatever is protected.

**Verification section and warnings appendix removed.** Both were requested, and both were
defensible: an appendix of caveats is where cautions go to be skipped. What was actionable in
them moved into the step it applies to — junction screening into transformation, the
repeat-array check into colony picking. `DesignReport.warnings` still exists for callers and
tests; it is simply not rendered.

**The clone table is editable and grows a row on demand.** You picked six colonies, not three.
