import sys
from pathlib import Path

import pytest

from mitopiper.config import Workspace
from mitopiper.references import ReferenceRegistry, ReferenceSet

DATA = Path(__file__).parent / "data"
REPO = Path(__file__).resolve().parent.parent
#: A ``mitopiper`` executable for job scripts run by the tests: this checkout, with this Python.
MITOPIPER_SHIM = f'#!/bin/sh\nPYTHONPATH="{REPO}" exec "{sys.executable}" -m mitopiper "$@"\n'


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


@pytest.fixture(autouse=True)
def no_site_lsf_settings(monkeypatch):
    """The tests must not depend on the machine: on an LSF cluster that reserves memory per slot, bparams made the LSF
    headers differ (test_headers failed there). Tests of the per-slot memory set memory_per_slot themselves."""
    from mitopiper import schedulers

    monkeypatch.setattr(schedulers, "lsf_reserves_per_slot", lambda: False)
