import os
from pathlib import Path

import pytest

_DEFAULT = Path(__file__).resolve().parents[2] / "plasmid"

PHL391 = "pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8.dna"
PCLM1 = "pCLM1_pcDNA3.1_miniCMV-mCherry-LambdaBoxBx8.dna"
# Same construct with the MfeI and EcoRI sites removed, so no restriction route exists.
TEST_PHL391 = "test_pHL391_pcDNA3.1_NFKBRE1-miniCMV-mCherry-LambdaBoxBx8.dna"


@pytest.fixture(scope="session")
def plasmid_dir() -> Path:
    """Directory holding the SnapGene inventory. Override with ``BBL_PLASMID_DIR``."""
    path = Path(os.environ.get("BBL_PLASMID_DIR", _DEFAULT))
    if not path.is_dir():
        pytest.skip(f"plasmid inventory not found at {path}; set BBL_PLASMID_DIR")
    return path


@pytest.fixture(scope="session")
def phl391(plasmid_dir) -> Path:
    return plasmid_dir / PHL391


@pytest.fixture(scope="session")
def pclm1(plasmid_dir) -> Path:
    return plasmid_dir / PCLM1


@pytest.fixture(scope="session")
def no_sites_phl391(plasmid_dir) -> Path:
    """pHL391 variant with MfeI/EcoRI removed -- the PCR-only case."""
    path = plasmid_dir / TEST_PHL391
    if not path.exists():
        pytest.skip(f"{TEST_PHL391} not present")
    return path
