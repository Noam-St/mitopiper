from pathlib import Path

import pytest

from mitopiper.config import Workspace
from mitopiper.references import ReferenceRegistry, ReferenceSet

DATA = Path(__file__).parent / "data"


@pytest.fixture
def data_dir():
    return DATA


@pytest.fixture
def workspace(tmp_path):
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    ws.activate("default")
    return ws


@pytest.fixture
def registered(workspace):
    reg = ReferenceRegistry(workspace)
    reg.add("Danio_rerio", ReferenceSet(refpath="/refs/dr/chrM.fa", gtfpath="/refs/dr/chrM.gff3"))
    reg.add("Homo_sapiens", ReferenceSet(refpath="/refs/hs/chrM.fa", gtfpath="/refs/hs/chrM.gff3"))
    return workspace
