"""Design session: the boundary between the model and the biology.

Every method here returns plain JSON-serialisable data, and **never returns a sequence**.
Products live in a session-side registry keyed by a short handle (``prod_1``); the model
passes handles around and can therefore not construct, paraphrase, or corrupt DNA. That makes
the "only ``replace_span`` builds sequence" invariant structural rather than aspirational.

This module has no dependency on the Anthropic SDK, so all of it is testable offline.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..config import load_lab_config
from ..enzymes import build_pool, enumerate_cut_sites
from ..excise import NoExcisionFound, excise_features
from ..insert import NoInsertionRoute, plan_insertion
from ..inventory import InventoryEntry, find_sequence, scan_inventory
from ..pcr import NoPrimerDesign, design_deletion_primers
from ..plasmid_io import circular_equal, feature_label, feature_span, write_genbank
from ..sourcing import source_insert
from ..sources import as_source

MAX_HITS = 12


@dataclass
class Product:
    """One designed plasmid, plus the plan and the records that produced it."""

    handle: str
    record: object
    plan: object
    origin: str
    #: The records the plan started from. Kept because a report needs to draw the parent and
    #: cost out the donor digest, neither of which the plan itself carries.
    parent: object = None
    donor: object = None
    #: ``(target_label, identical)`` for every ``compare_product`` run against this product.
    comparisons: list = field(default_factory=list)

    @property
    def protocol(self) -> str:
        return getattr(self.plan, "protocol", "(no protocol)")


@dataclass
class DesignSession:
    """Holds the inventory and every product designed during a conversation."""

    #: Anything ``bbl.sources.as_source`` accepts: a directory, a list of files, a mapping of
    #: records, or a PlasmidSource such as a Benchling adapter.
    library: object
    entries: list[InventoryEntry] = field(default_factory=list)
    products: dict[str, Product] = field(default_factory=dict)
    config: dict = field(default_factory=load_lab_config)
    #: Called before any outward-facing action. Return False to decline. None = auto-approve.
    confirm: object = None
    _counter: int = 0

    def __post_init__(self):
        self.source = as_source(self.library)
        if not self.entries:
            self.entries, _ = scan_inventory(self.source)

    # -- helpers ---------------------------------------------------------------

    def resolve(self, name: str) -> InventoryEntry:
        """Find a plasmid by short name (``pHL391``) or full label.

        The plasmid ID -- the label up to the first underscore -- is its own match tier, above
        prefix matching. Without it ``pCLM2`` is ambiguous against pCLM21/22/23/24, which is
        exactly the name a user or model will type.
        """
        wanted = str(name).strip().lower()
        exact = [e for e in self.entries if e.label.lower() == wanted]
        ident = [e for e in self.entries if e.label.split("_")[0].lower() == wanted]
        prefix = [e for e in self.entries if e.label.lower().startswith(wanted)]
        contains = [e for e in self.entries if wanted in e.label.lower()]
        for group in (exact, ident, prefix, contains):
            if len(group) == 1:
                return group[0]
            if len(group) > 1:
                raise LookupError(
                    f"{name!r} matches {len(group)} plasmids: "
                    f"{', '.join(sorted(e.label for e in group)[:5])}"
                )
        raise LookupError(f"no plasmid named {name!r} in {self.source.name}")

    def _register(self, record, plan, origin: str, parent=None, donor=None) -> str:
        self._counter += 1
        handle = f"prod_{self._counter}"
        self.products[handle] = Product(handle, record, plan, origin, parent, donor)
        return handle

    def protocol_for(self, handle: str) -> str:
        """Verbatim protocol text. For the harness to append -- not a model-facing tool."""
        return self.products[handle].protocol

    # -- tools -----------------------------------------------------------------

    def search_inventory(self, query: str) -> dict:
        """Find plasmids by name, by feature label, or by exact sequence."""
        query = str(query).strip()
        lowered = query.lower()
        by_name, by_feature = [], []

        for entry in self.entries:
            if lowered in entry.label.lower():
                by_name.append({"plasmid": entry.label, "length_bp": entry.length})

        if set(query.upper()) <= set("ACGTRYSWKMBDHVN") and len(query) >= 12:
            hits = find_sequence(self.entries, query)
            return {
                "query": f"<{len(query)} bp sequence>",
                "sequence_hits": [
                    {"plasmid": h.label, "start": h.start, "strand": h.strand} for h in hits
                ][:MAX_HITS],
                "hit_count": len(hits),
            }

        for entry in self.entries:
            record = entry.load()
            for feature in record.features:
                if feature.type == "primer_bind":
                    continue
                label = feature_label(feature)
                if lowered in label.lower():
                    start, end = feature_span(feature)
                    by_feature.append(
                        {
                            "plasmid": entry.label,
                            "feature": label,
                            "span": [start, end],
                            "length_bp": end - start,
                        }
                    )
                    break

        return {
            "query": query,
            "name_matches": by_name[:MAX_HITS],
            "feature_matches": by_feature[:MAX_HITS],
            "base_vectors": self.config.get("base_vectors", []),
        }

    def inspect_plasmid(self, plasmid: str) -> dict:
        """Length, functional features, and unique cutters of one plasmid."""
        entry = self.resolve(plasmid)
        record = entry.load()
        features = [
            {"label": feature_label(f), "type": f.type, "span": list(feature_span(f))}
            for f in record.features
            if f.type != "primer_bind"
        ]
        cuts = enumerate_cut_sites(record.seq, build_pool("CommOnly"))
        unique = sorted(
            ({"enzyme": name, "cut": sites[0].top} for name, sites in cuts.items() if len(sites) == 1),
            key=lambda c: c["cut"],
        )
        return {
            "plasmid": entry.label,
            "length_bp": entry.length,
            "features": sorted(features, key=lambda f: f["span"][0]),
            "unique_cutter_count": len(unique),
            "unique_cutters": unique[:40],
        }

    def plan_deletion(self, plasmid: str, features: list[str], method: str = "auto") -> dict:
        """Remove features. Tries restriction first unless ``method`` says otherwise."""
        entry = self.resolve(plasmid)
        record = entry.load()
        attempts = ["restriction", "pcr"] if method == "auto" else [method]
        reasons: dict = {}

        for attempt in attempts:
            try:
                if attempt == "restriction":
                    plan = excise_features(record, features)
                    handle = self._register(
                        plan.product, plan, f"deletion from {entry.label}", parent=record
                    )
                    return {
                        "feasible": True,
                        "method": "restriction",
                        "product_id": handle,
                        "enzymes": list(plan.enzyme_pair),
                        "strategy": plan.strategy,
                        "length_bp": len(plan.product),
                        "deleted_bp": plan.deleted_bp,
                        "collateral_bp": plan.collateral_bp,
                        "pcr_derived_bp": 0,
                        "new_oligos": 0,
                        "removed_features": plan.removed_features,
                        "junction": plan.junction_seq,
                        "sites_regenerated": plan.sites_regenerated,
                        "warnings": plan.warnings,
                        "verified": True,
                    }
                design = design_deletion_primers(
                    record, features, warn_if_restriction_possible=False
                )
                handle = self._register(
                    design.product, design, f"deletion from {entry.label}", parent=record
                )
                return {
                    "feasible": True,
                    "method": "pcr",
                    "product_id": handle,
                    "primers": [str(design.forward), str(design.reverse)],
                    "length_bp": len(design.product),
                    "deleted_bp": design.deleted_bp,
                    "collateral_bp": 0,
                    "pcr_derived_bp": design.amplicon_bp,
                    "new_oligos": 2,
                    "removed_features": design.removed_features,
                    "warnings": design.warnings,
                    "verified": True,
                }
            except (NoExcisionFound, NoPrimerDesign) as exc:
                reasons[attempt] = getattr(exc, "reasons", None) or {"reason": str(exc)}
            except ValueError as exc:
                return {"feasible": False, "error": str(exc)}

        return {
            "feasible": False,
            "reasons": reasons,
            "suggestion": "no cloning route; consider ordering the construct",
        }

    def plan_insertion(
        self,
        backbone: str,
        at: str,
        sequence: str | None = None,
        donor: str | None = None,
        donor_features: list[str] | None = None,
        method: str = "auto",
    ) -> dict:
        """Insert a sequence (or a feature taken from a donor plasmid) into a backbone."""
        entry = self.resolve(backbone)
        record = entry.load()
        site = _parse_site(at)
        insert = sequence
        features = None
        donor_record = None
        if donor is not None:
            donor_record = self.resolve(donor).load()
            insert = donor_record
            features = donor_features
        if insert is None:
            return {"feasible": False, "error": "give either sequence= or donor="}

        try:
            plan = plan_insertion(
                record, insert, at=site, insert_features=features, method=method
            )
        except NoInsertionRoute as exc:
            return {"feasible": False, "reasons": exc.reasons or {"reason": str(exc)}}
        except (ValueError, LookupError) as exc:
            return {"feasible": False, "error": str(exc)}

        handle = self._register(
            plan.product,
            plan,
            f"insertion into {entry.label}",
            parent=record,
            donor=donor_record,
        )
        return {
            "feasible": True,
            "method": plan.strategy,
            "product_id": handle,
            "length_bp": len(plan.product),
            "inserted_bp": plan.inserted_bp,
            "replaced_bp": plan.replaced_bp,
            "site": list(plan.site),
            "directional": plan.directional,
            "vector_prep": plan.vector.preparation,
            "insert_prep": plan.insert.preparation,
            "enzymes": list(plan.vector.enzymes) if plan.vector.enzymes else None,
            "warnings": plan.warnings,
            "verified": True,
        }

    def source_sequence(self, name: str, sequence: str | None = None) -> dict:
        """Decide how to obtain a sequence: from the library, synthesis, or order the plasmid."""
        decision = source_insert(sequence, self.entries, name=name)
        return {
            "route": decision.route,
            "summary": decision.summary(),
            "donors": sorted({d.label for d in decision.donors}),
            "rationale": decision.rationale,
            "question": decision.question,
            "complexity": (
                {
                    "synthesizable": decision.complexity.synthesizable,
                    "score": decision.complexity.score,
                    "reasons": decision.complexity.reasons,
                    "is_placeholder": decision.complexity.is_placeholder,
                }
                if decision.complexity
                else None
            ),
        }

    def compare_product(self, product_id: str, target: str) -> dict:
        """Check a designed product against an existing plasmid, ignoring rotation."""
        if product_id not in self.products:
            return {"error": f"unknown product {product_id!r}"}
        entry = self.resolve(target)
        expected = entry.load()
        product = self.products[product_id].record
        identical = circular_equal(product, expected)
        # Remembered so a report written later can state what the design was checked against,
        # without the model having to carry the verdict back in through an argument.
        self.products[product_id].comparisons.append((entry.label, identical))
        return {
            "product_id": product_id,
            "target": target,
            "identical": identical,
            "product_bp": len(product),
            "target_bp": len(expected),
            "difference_bp": len(product) - len(expected),
        }

    def export_product(self, product_id: str, path: str) -> dict:
        """Write a product to GenBank. Outward-facing -- goes through ``confirm`` first."""
        if product_id not in self.products:
            return {"error": f"unknown product {product_id!r}"}
        product = self.products[product_id]
        if self.confirm is not None and not self.confirm(
            f"write {len(product.record)} bp to {path}"
        ):
            return {
                "declined": True,
                "note": "the user declined; do not retry, ask what they want instead",
            }
        written = write_genbank(product.record, path)
        return {
            "product_id": product_id,
            "path": str(written),
            "length_bp": len(product.record),
            "features": len(product.record.features),
        }

    def dna_needing_concentration(self, product_id: str) -> list[dict]:
        """Which stocks the report's volumes depend on, so the user can be asked for them.

        Every reaction volume is a mass divided by a concentration, and the concentration
        comes off a Nanodrop after the design exists. Asking for two numbers up front turns
        the tables from a form into a protocol.
        """
        if product_id not in self.products:
            return []
        product = self.products[product_id]
        needed = []
        for record, role in ((product.parent, "parent / backbone"), (product.donor, "donor")):
            if record is None:
                continue
            label = self._label_for(record)
            needed.append({"plasmid": label, "role": role, "length_bp": len(record)})
        return needed

    def _label_for(self, record) -> str:
        """The full inventory label for a record, falling back to its short name."""
        for entry in self.entries:
            if entry.length == len(record) and entry.label.startswith(str(record.name)):
                return entry.label
        return str(getattr(record, "name", "unknown"))

    def generate_report(
        self,
        product_id: str,
        path: str,
        aim: str | None = None,
        name: str | None = None,
        concentrations: dict | None = None,
    ) -> dict:
        """Write a bench-ready cloning report. Outward-facing -- goes through ``confirm``.

        Everything in the report comes from the stored plan and the records it was built from.
        ``aim`` is the user's sentence about why the construct exists; ``concentrations`` maps
        a plasmid name to ng/uL and fills in the reaction volumes. Anything not supplied stays
        an input box in the report rather than a guess.
        """
        from ..report import build_report, figures_available

        if product_id not in self.products:
            return {"error": f"unknown product {product_id!r}"}
        product = self.products[product_id]
        if self.confirm is not None and not self.confirm(f"write a cloning report to {path}"):
            return {
                "declined": True,
                "note": "the user declined; do not retry, ask what they want instead",
            }

        written = Path(path)
        parent_name = self._label_for(product.parent) if product.parent is not None else None
        donor_name = self._label_for(product.donor) if product.donor is not None else None
        report = build_report(
            product.plan,
            parent=product.parent,
            product=product.record,
            donor=product.donor,
            aim=aim,
            name=name,
            parent_name=parent_name,
            donor_name=donor_name,
            concentrations=_numeric(concentrations),
            config=self.config,
            verified_against=product.comparisons[-1] if product.comparisons else None,
        )
        written.write_text(report.to_html())
        return {
            "product_id": product_id,
            "path": str(written),
            "sections": len(report.steps),
            "figures": len(report.figures),
            "maps_included": figures_available(),
            "concentrations_supplied": sorted(_numeric(concentrations)),
            "concentrations_missing": [
                item["plasmid"]
                for item in self.dna_needing_concentration(product_id)
                if item["plasmid"] not in _numeric(concentrations)
            ],
        }


def _numeric(concentrations) -> dict:
    """Keep the entries that are actually numbers. A model may pass "unknown" or "".

    Dropping a junk value is right: it becomes an input box in the report, which is honest,
    where coercing it to zero would silently divide by nothing.
    """
    clean = {}
    for key, value in (concentrations or {}).items():
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if number > 0:
            clean[str(key)] = number
    return clean


def _parse_site(at: str):
    """``"NFKBRE"`` | ``"1491"`` | ``"1491:1512"`` -> the form plan_insertion expects."""
    text = str(at).strip()
    if ":" in text:
        start, end = text.split(":", 1)
        return (int(start), int(end))
    if text.lstrip("-").isdigit():
        return int(text)
    return text
