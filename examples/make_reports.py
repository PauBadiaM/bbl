"""Regenerate the example cloning reports.

    python examples/make_reports.py [--plasmids DIR] [--out DIR]

The first example is the hackathon's driving request, run through the same
:class:`~bbl.llm.session.DesignSession` the model-facing tools call:

    Have a control for my sensor that has these features:
      Backbone:       pcDNA3.1
      Features:       miniCMV-mCherry-LambdaBoxBx8-polyA
      Ab resistance:  Ampicillin

The tool sequence below is the one the model issues for that request -- search the inventory,
inspect the candidate, plan the deletion, check the product, write the report. Running it here
rather than through ``python -m bbl.llm`` keeps the examples reproducible with no API key; the
conversation adds the judgement around these calls, not the calls themselves.

The other two cover the remaining routes:

    pCLM2_designed_report.html     restriction subcloning, pCLM3 + BoxB x8
    inverse_pcr_route_report.html  inverse PCR, on the variant with the sites removed

They are checked in so the output is visible without running anything, but they are
generated artefacts -- regenerate rather than edit.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from bbl import design_deletion_primers, load_plasmid, plan_insertion
from bbl.llm.session import DesignSession
from bbl.report import figures_available, write_report

PCLM1 = "pCLM1_pcDNA3.1_miniCMV-mCherry-LambdaBoxBx8.dna"
PCLM3 = "pCLM3_pcDNA3.1_CMV-mCherry.dna"
NO_SITES = "test_pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8.dna"

AIM = (
    "A control for the NF-kB sensor: pcDNA3.1 carrying miniCMV-mCherry-LambdaBoxBx8-polyA "
    "with AmpR, i.e. the sensor cassette driven by the bare minimal promoter, so any mCherry "
    "signal is response-element-independent background."
)

#: What the parent measured at on the Nanodrop. Ask the user for this; without it the report's
#: reaction volumes stay as input boxes rather than being guessed at.
PARENT_NG_PER_UL = 1364


def sensor_control(library: Path, out: Path) -> Path:
    """The driving request, through the agent's own tools."""
    session = DesignSession(library)

    # 1. Which library plasmids carry the wanted parts? pHL391 has all of them and an extra
    #    NF-kB response element -- deleting one feature beats building the cassette.
    session.search_inventory("miniCMV")
    session.inspect_plasmid("pHL391")

    # 2. Plan the edit and check it against the construct the lab actually built.
    plan = session.plan_deletion("pHL391", ["NFKBRE"])
    handle = plan["product_id"]
    session.compare_product(handle, "pCLM1")

    # 3. Ask for the stock concentration, then write the report.
    needed = session.dna_needing_concentration(handle)
    concentrations = {needed[0]["plasmid"]: PARENT_NG_PER_UL} if needed else {}
    result = session.generate_report(
        handle,
        str(out / "sensor_control_report.html"),
        aim=AIM,
        name="pCLM1_designed",
        concentrations=concentrations,
    )
    return Path(result["path"])


def main(argv=None) -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--plasmids",
        default=os.environ.get("BBL_PLASMID_DIR", here.parents[1] / "plasmid"),
        help="directory of .dna files (default: $BBL_PLASMID_DIR, then ../plasmid)",
    )
    parser.add_argument("--out", default=here, help="where to write (default: this folder)")
    args = parser.parse_args(argv)

    library, out = Path(args.plasmids), Path(args.out)
    if not library.is_dir():
        parser.error(f"no plasmid inventory at {library}; pass --plasmids")
    out.mkdir(parents=True, exist_ok=True)
    if not figures_available():
        print('no plasmid maps: pip install --only-binary=:all: "bbl[report]"')

    written = [sensor_control(library, out)]

    backbone, donor = load_plasmid(library / PCLM3), load_plasmid(library / PCLM1)
    written.append(
        write_report(
            plan_insertion(backbone, donor, insert_features=["Lambda BoxB x8"], at=(1491, 1512)),
            out / "pCLM2_designed_report.html",
            parent=backbone,
            donor=donor,
            parent_name=Path(PCLM3).stem,
            donor_name=Path(PCLM1).stem,
            name="pCLM2_designed",
            aim="Add the lambda BoxB x8 array to the CMV-mCherry construct.",
            verified_against=("pCLM2", True),
        )
    )

    variant = library / NO_SITES
    if variant.exists():
        template = load_plasmid(variant)
        written.append(
            write_report(
                design_deletion_primers(template, ["NFKBRE"]),
                out / "inverse_pcr_route_report.html",
                parent=template,
                parent_name=Path(NO_SITES).stem,
                name="pHL391_dNFKBRE",
                aim="The same deletion on a variant with no usable restriction sites.",
            )
        )

    for path in written:
        print(f"{path}  ({path.stat().st_size // 1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
