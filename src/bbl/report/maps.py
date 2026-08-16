"""Plasmid figures, drawn with DnaFeaturesViewer.

Two kinds of view, because no single one shows what a cloning step did:

* a **circular** whole-plasmid map, for orientation -- parent and product side by side;
* a **linear zoom** on the edit, where a 66 bp change to a 5.7 kb plasmid is actually visible.

Most of the work here is not drawing but *deciding what to draw*. A SnapGene record is not a
clean feature list: it carries primer bindings, umbrella annotations spanning half the plasmid
("Insert Sequence"), and sub-features nested inside the part you care about (five NF-kB sites
inside ``NFKBRE``, eight monomers inside ``Lambda BoxB x8``). Drawn literally, the labels bury
the figure. :func:`prune_features` is the editorial pass -- see ``docs/DECISIONS.md``.

``matplotlib`` and ``dna_features_viewer`` are optional extras -- ``import bbl`` must keep
working without them -- so both are imported lazily and their absence raises a message saying
what to install.

Figures come back as inline-ready SVG strings rather than files: a report is one
self-contained HTML document, so there is nothing to keep together and nothing to lose.
"""

from __future__ import annotations

import io
import re

from ..plasmid_io import NON_FUNCTIONAL_TYPES, feature_label, feature_span

#: On a whole-plasmid map, features below this are sub-pixel and their labels crowd out the
#: ones that matter. The edit and anything named in ``protect`` are kept regardless.
MIN_CIRCULAR_BP = 60

#: Feature colours, keyed by GenBank feature type.
FEATURE_COLORS = {
    "CDS": "#bcd7ee",
    "gene": "#bcd7ee",
    "promoter": "#c6e5c0",
    "enhancer": "#c6e5c0",
    "regulatory": "#c6e5c0",
    "terminator": "#f3cdcd",
    "polyA_signal": "#f3cdcd",
    "rep_origin": "#ecdfbe",
    "oriT": "#ecdfbe",
    "protein_bind": "#ddd0ec",
    "misc_RNA": "#ddd0ec",
    "ncRNA": "#ddd0ec",
}
DEFAULT_COLOR = "#e2e2e2"

#: Sequence that leaves in this step, sequence that arrives, and the cut/junction rules.
REMOVED_COLOR = "#eda3a3"
ADDED_COLOR = "#9fd39a"
MARK_COLOR = "#c0392b"
_TEXT = "#333333"


class MissingFigureDependency(RuntimeError):
    """Raised when the plotting extras are not installed."""


def _plotting():
    """Import the optional plotting stack, or explain how to get it."""
    try:
        import matplotlib

        matplotlib.use("Agg")  # no display on a cluster node; must precede pyplot
        import matplotlib.pyplot as plt
        from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
    except ImportError as exc:  # pragma: no cover - covered by the skip in the tests
        raise MissingFigureDependency(
            "plasmid figures need matplotlib and dna_features_viewer: "
            'pip install --only-binary=:all: "bbl[report]"'
        ) from exc
    return plt, BiopythonTranslator, CircularGraphicRecord


def figures_available() -> bool:
    """True when a report can carry plasmid maps."""
    try:
        _plotting()
    except MissingFigureDependency:
        return False
    return True


# ---------------------------------------------------------------------------
# editorial pass: which features to draw
# ---------------------------------------------------------------------------


def _contains(outer, inner) -> bool:
    (o_start, o_end), (i_start, i_end) = feature_span(outer), feature_span(inner)
    return o_start <= i_start and i_end <= o_end and (o_end - o_start) > (i_end - i_start)


def prune_features(features, min_bp=0, protect=(), limit=None):
    """Trim a SnapGene feature list down to what a reader can actually take in.

    Four passes, in order:

    1. drop bookkeeping types (primer bindings, ``source``);
    2. drop features shorter than ``min_bp``;
    3. resolve nesting one level: a feature containing two or more others is either an
       **umbrella** (its children have different labels -- "Insert Sequence" over a whole
       cassette, so the children are the informative ones and the umbrella goes) or an
       **array** (its children all share one label -- eight "BoxB RNA aptamer" monomers, so
       the umbrella "Lambda BoxB x8" is the informative one and the children go);
    4. drop whatever is still nested inside a retained feature.

    ``protect`` is a set of ``(start, end)`` spans that survive every pass -- the part being
    cut out or put in must be labelled even when the rules above would fold it away.
    ``limit`` caps how many survive, largest first, protected spans on top of the budget.
    """
    protect = {tuple(span) for span in protect}

    def is_protected(feature) -> bool:
        return feature_span(feature) in protect

    kept = [
        f
        for f in features
        if f.type not in NON_FUNCTIONAL_TYPES
        and (is_protected(f) or (feature_span(f)[1] - feature_span(f)[0]) >= min_bp)
    ]

    dropped = set()
    for outer in kept:
        children = [f for f in kept if _contains(outer, f) and id(f) not in dropped]
        if len(children) < 2:
            continue
        if len({feature_label(f) for f in children}) > 1:
            if not is_protected(outer):
                dropped.add(id(outer))  # umbrella: keep the parts, not the wrapper
        else:
            dropped.update(id(f) for f in children if not is_protected(f))  # repeat array

    kept = [f for f in kept if id(f) not in dropped]
    kept = [
        f
        for f in kept
        if is_protected(f) or not any(_contains(other, f) for other in kept)
    ]

    if limit is not None and len(kept) > limit:
        budget = sorted(
            (f for f in kept if not is_protected(f)),
            key=lambda f: feature_span(f)[0] - feature_span(f)[1],  # longest first
        )[:limit]
        keep = {id(f) for f in budget} | {id(f) for f in kept if is_protected(f)}
        kept = [f for f in kept if id(f) in keep]
    return kept


def _build_translator(base, removed=None, added=None, min_bp=0, protect=(), limit=None):
    """A ``BiopythonTranslator`` that speaks this library's conventions.

    ``removed`` and ``added`` are spans: a feature lying wholly inside one is coloured as
    leaving or arriving. A feature merely *straddling* the span keeps its normal colour --
    it survives the edit, truncated, and colouring it as deleted would misreport that.
    """

    def within(feature, span) -> bool:
        if span is None:
            return False
        start, end = feature_span(feature)
        return span[0] <= start and end <= span[1]

    class _Translator(base):
        def compute_filtered_features(self, features):
            return prune_features(features, min_bp=min_bp, protect=protect, limit=limit)

        def compute_feature_color(self, feature):
            if within(feature, removed):
                return REMOVED_COLOR
            if within(feature, added):
                return ADDED_COLOR
            return FEATURE_COLORS.get(feature.type, DEFAULT_COLOR)

        def compute_feature_label(self, feature):
            label = feature_label(feature)
            return label if len(label) <= 26 else label[:25] + "…"

        def compute_feature_linewidth(self, feature):
            return 1.3 if (within(feature, removed) or within(feature, added)) else 0.6

        def compute_feature_box_color(self, feature):
            return "#ffffff"

        def compute_feature_box_linewidth(self, feature):
            return 0

        def compute_feature_label_link_color(self, feature):
            return "#909090"

    return _Translator()


def _fig_to_svg(figure, salt: str) -> str:
    """Serialise ``figure`` as an SVG fragment safe to paste into an HTML document.

    Several figures share one page, so their internal ids (clip paths, gradients) must not
    collide -- matplotlib derives those from ``svg.hashsalt``, so a distinct salt per figure
    keeps them apart. The XML prolog and DOCTYPE are stripped: legal in a standalone file,
    not inside ``<body>``. Width and height are dropped so the page's CSS sets the size and
    the ``viewBox`` keeps the aspect ratio.
    """
    import matplotlib
    import matplotlib.pyplot as plt

    buffer = io.StringIO()
    with matplotlib.rc_context({"svg.hashsalt": salt}):
        figure.savefig(buffer, format="svg", bbox_inches="tight", transparent=True)
    plt.close(figure)
    svg = buffer.getvalue()
    svg = re.sub(r"<\?xml[^>]*\?>\s*", "", svg)
    svg = re.sub(r"<!DOCTYPE[^>]*>\s*", "", svg, flags=re.IGNORECASE)
    svg = re.sub(r"<!--.*?-->", "", svg, count=1, flags=re.DOTALL)
    svg = re.sub(r'(<svg[^>]*?)\swidth="[^"]*"', r"\1", svg, count=1)
    svg = re.sub(r'(<svg[^>]*?)\sheight="[^"]*"', r"\1", svg, count=1)
    return svg.strip()


# ---------------------------------------------------------------------------
# public figures
# ---------------------------------------------------------------------------


#: Most labels a circular map can carry before it stops being readable. The largest features
#: win; anything in ``protect`` is kept on top of the budget.
MAX_CIRCULAR_LABELS = 14

#: Angular room one label needs at the innermost ring, in degrees. Labels closer together than
#: this are pushed out to the next ring rather than overprinted.
_LABEL_ARC_DEGREES = 11.0
_LABEL_RINGS = (1.10, 1.34, 1.58)


def _wrap_name(name, width=22) -> list[str]:
    """Break a full construct name into lines that fit inside the circle.

    Construct names are structured -- ``pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8``
    -- so breaking on the separators keeps each line meaningful, where a plain character wrap
    would cut through the middle of a part name.
    """
    words, lines, current = re.split(r"(?<=[_-])", str(name)), [], ""
    for word in words:
        if current and len(current) + len(word) > width:
            lines.append(current.rstrip("_"))
            current = word
        else:
            current += word
    if current:
        lines.append(current.rstrip("_"))
    return lines or [str(name)]


def _place_radial_labels(axes, graphic, labels, colour=_TEXT):
    """Write feature names around the circle, reading outward like a plasmid map.

    DnaFeaturesViewer stacks circular labels in a block above the plot, which on a real
    SnapGene record is taller than the plasmid and impossible to match to its feature. Here
    each label sits at its own feature's angle with a short leader, so the map reads the way
    SnapGene's does: names follow the circle, and they use the bottom as well as the top.

    Crowded neighbours move to an outer ring instead of overprinting each other.
    """
    import numpy as np

    centre_y, radius = -graphic.radius, graphic.radius
    placed: list[tuple[float, int]] = []

    ordered = sorted(labels, key=lambda item: (item[0].start + item[0].end) / 2)
    for feature, text in ordered:
        if not text:
            continue
        angle = graphic.position_to_angle((feature.start + feature.end) / 2)
        radians = np.radians(angle)

        ring = 0
        while ring < len(_LABEL_RINGS) - 1 and any(
            other_ring == ring and abs((angle - other_angle + 180) % 360 - 180)
            < _LABEL_ARC_DEGREES
            for other_angle, other_ring in placed
        ):
            ring += 1
        placed.append((angle, ring))

        distance = _LABEL_RINGS[ring] * radius
        rightwards = np.cos(radians) >= 0
        axes.annotate(
            text,
            xy=(1.02 * radius * np.cos(radians), centre_y + 1.02 * radius * np.sin(radians)),
            xytext=(distance * np.cos(radians), centre_y + distance * np.sin(radians)),
            rotation=angle if rightwards else angle + 180,
            rotation_mode="anchor",
            ha="left" if rightwards else "right",
            va="center",
            fontsize=8.5,
            color=colour,
            annotation_clip=False,
            arrowprops={"arrowstyle": "-", "lw": 0.5, "color": "#b0b0b0",
                        "shrinkA": 1, "shrinkB": 0},
        )
    return _LABEL_RINGS[max((ring for _, ring in placed), default=0)]


def circular_svg(
    record,
    title=None,
    removed=None,
    added=None,
    protect=(),
    width=6.6,
    salt="circular",
) -> str:
    """Whole-plasmid circular map with radial labels, captioned with its name and length."""
    plt, translator_base, circular_class = _plotting()

    translator = _build_translator(
        translator_base,
        removed=removed,
        added=added,
        min_bp=MIN_CIRCULAR_BP,
        protect=protect,
        limit=MAX_CIRCULAR_LABELS,
    )
    graphic = translator.translate_record(record, record_class=circular_class)

    # Capture the labels, then blank them so DnaFeaturesViewer draws the arcs only and we own
    # the text placement entirely.
    labels = [(feature, feature.label) for feature in graphic.features]
    for feature in graphic.features:
        feature.label = None

    with plt.rc_context({"font.size": 9}):
        axes, _ = graphic.plot(figure_width=width)
        reach = _place_radial_labels(axes, graphic, labels)

        caption = _wrap_name(title or record.name) + [f"{len(record.seq):,} bp"]
        axes.text(
            0, -graphic.radius, "\n".join(caption), ha="center", va="center",
            fontsize=10.5 if len(caption) <= 2 else 9.5, linespacing=1.6, color=_TEXT,
        )
        # The labels stick out past whatever DnaFeaturesViewer sized the axes for.
        margin = reach * graphic.radius + 0.55
        axes.set_xlim(-margin, margin)
        axes.set_ylim(-graphic.radius - margin, -graphic.radius + margin)
        axes.set_aspect("equal")
        axes.axis("off")
        axes.figure.set_size_inches(width, width)
        return _fig_to_svg(axes.figure, salt)


def linear_svg(
    record,
    window,
    marks=(),
    removed=None,
    added=None,
    protect=(),
    width=9.4,
    salt="linear",
) -> str:
    """Zoomed linear view of ``window`` = ``(start, end)``.

    ``marks`` are ``(position, label)`` pairs drawn as dashed rules -- the two cut sites on the
    parent, the new junction on the product. Coordinates are in the record's own frame, which
    for a product built by ``replace_span`` is still the parent's frame, so the two panels of a
    before/after pair line up.
    """
    plt, translator_base, _ = _plotting()

    start, end = max(0, int(window[0])), min(len(record.seq), int(window[1]))
    translator = _build_translator(
        translator_base, removed=removed, added=added, min_bp=0, protect=protect
    )
    graphic = translator.translate_record(record).crop((start, end))
    graphic.labels_spacing = 6

    with plt.rc_context({"font.size": 9}):
        figure, axes = plt.subplots(figsize=(width, 1.5))
        graphic.plot(ax=axes, with_ruler=True, draw_line=True)

        for span, colour in ((removed, REMOVED_COLOR), (added, ADDED_COLOR)):
            if span and span[1] > span[0]:
                axes.axvspan(span[0], span[1], color=colour, alpha=0.28, zorder=0, lw=0)
        for position, label in marks:
            if not start <= position <= end:
                continue
            axes.axvline(position, color=MARK_COLOR, ls="--", lw=1.1, zorder=3)
            axes.annotate(
                label,
                xy=(position, axes.get_ylim()[1]),
                xytext=(0, 1),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8.5,
                color=MARK_COLOR,
            )
        return _fig_to_svg(figure, salt)


def zoom_window(span, length, context=300):
    """A ``(start, end)`` viewing window around ``span``, clamped to the record."""
    start, end = int(span[0]), int(span[1])
    if end < start:
        start, end = end, start
    pad = max(context, (end - start) // 4)
    return max(0, start - pad), min(int(length), end + pad)
