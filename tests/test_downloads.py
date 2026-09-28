import io
import os
import subprocess

import pytest

from mitopiper import monitor as mon
from mitopiper.cli import main
from mitopiper.config import Workspace
from mitopiper.local_samples import make_tsv
from mitopiper.references import ReferenceRegistry, ReferenceSet
from mitopiper.rerun import RerunOptions, plan, rerun

from test_arrays import FAKE_SBATCH


def workspace(tmp_path, mode, scheduler_queue=""):
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("default")
    wf.set_general("download_mode", mode)
    wf.set_general("download_queue", scheduler_queue)
    wf.set_general("tool_check", False)
    wf.pipeline["count"] = {"key": "wc", "description": "count", "switch": True, "path": "", "outdir": "counts",
                            "inputdir": "fastq", "input_ext": ".fastq.gz", "output_ext": ".txt", "scope": "sample",
                            "structure": "key -c input output_pipe"}
    wf.run_order = ["wget", "count"]
    ws.write_workflow(wf, "t")
    fq = tmp_path / "fq"
    fq.mkdir()
    for s in ("a", "b"):
        (fq / f"{s}.fastq").write_text("@r\nACGT\n+\nIIII\n")
    make_tsv(fq, "Homo sapiens", tmp_path / "exp.tsv")
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath="/x.fa"))
    return ws, fq


def build(ws, tmp_path, scheduler):
    assert main(["-w", str(ws.root), "build", str(tmp_path / "exp.tsv"), "-s", scheduler, "-p", "t"]) == 0


def states(ws):
    return {j.name: (j.state, j.kind) for j in mon.collect(ws, ["exp"], mon.ErrorScanner(), lambda s, i: None)}


def test_login_mode_downloads_before_the_samples(tmp_path):
    ws, _ = workspace(tmp_path, "login")
    build(ws, tmp_path, "local")
    scripts = ws.output_dir / "exp" / "scripts"
    assert "gzip -c" in (scripts / "dl_exp.sh").read_text()
    assert "pm_step download:a_0_HS" in (scripts / "dl_exp.sh").read_text()
    sample = (scripts / "s_a_0.sh").read_text()
    assert "gzip -c" not in sample and "for pm_dep in dl_exp" in sample
    master = (ws.root / "master.sh").read_text()
    assert "run_here " in master and master.index("run_here /") < master.index("for script in")
    subprocess.run(["bash", str(ws.root / "master.sh")], check=True, capture_output=True)
    assert states(ws) == {"dl_exp": (mon.DONE, "download"), "s_a_0": (mon.DONE, "sample"), "s_b_1": (mon.DONE, "sample")}


def test_queue_mode_submits_one_download_job_on_the_transfer_queue(tmp_path, monkeypatch):
    bin_dir, state = tmp_path / "bin", tmp_path / "state"
    bin_dir.mkdir()
    state.mkdir()
    (bin_dir / "sbatch").write_text(FAKE_SBATCH)
    (bin_dir / "sbatch").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_STATE", str(state))
    ws, _ = workspace(tmp_path, "queue", "transfer")
    build(ws, tmp_path, "slurm")
    dl = (ws.output_dir / "exp" / "scripts" / "dl_exp.sh").read_text()
    assert "#SBATCH --partition=transfer" in dl and "#SBATCH --cpus-per-task=1" in dl
    subprocess.run(["bash", str(ws.root / "master.sh")], check=True, capture_output=True)
    calls = (state / "calls").read_text().splitlines()
    assert calls[0].endswith("dl_exp.sh") and "afterok:101" in calls[1]      # samples wait for the download job
    assert all(s == mon.DONE for s, _ in states(ws).values())


def test_failed_download_skips_its_samples_and_a_rerun_recovers(tmp_path):
    ws, fq = workspace(tmp_path, "login")
    build(ws, tmp_path, "local")
    (fq / "b.fastq").unlink()                                   # the download of b fails
    subprocess.run(["bash", str(ws.root / "master.sh")], check=True, capture_output=True)
    got = states(ws)
    assert got["dl_exp"][0] == mon.FAILED and got["s_a_0"][0] == mon.SKIPPED and got["s_b_1"][0] == mon.SKIPPED
    assert {i.name for i in plan(ws, ["exp"], query=lambda s, i: None)} == {"dl_exp", "s_a_0", "s_b_1"}
    (fq / "b.fastq").write_text("@r\nACGT\n+\nIIII\n")
    assert rerun(ws, ["exp"], RerunOptions(), submit=True, stream=io.StringIO()) == 0
    assert all(s == mon.DONE for s, _ in states(ws).values())
    assert (ws.output_dir / "exp" / "counts" / "b_1_HS.txt").exists()


def test_invalid_download_mode(tmp_path):
    ws, _ = workspace(tmp_path, "sometimes")
    assert main(["-w", str(ws.root), "build", str(tmp_path / "exp.tsv"), "-p", "t"]) == 2
