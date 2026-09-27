"""Regression tests for the issues found in code review."""
import io
import json
import shutil
import subprocess

import pytest

from pipemaker import monitor as mon
from pipemaker.cli import main
from pipemaker.config import Workspace
from pipemaker.local_samples import make_tsv
from pipemaker.references import ReferenceRegistry, ReferenceSet, index_commands
from pipemaker.rerun import RerunOptions, plan

SHELLS = [s for s in ("bash", "sh", "tcsh", "csh") if shutil.which(s)]


def workspace_with(tmp_path, shell="bash", fail=True, experiment_step=True, steps=("wget", "boom")):
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("default")
    wf.set_general("environment", shell)
    wf.pipeline["boom"] = {"key": "ls", "description": "list", "switch": True, "path": "", "outdir": "", "inputdir": "",
                           "input_ext": "", "output_ext": "", "scope": "sample",
                           "structure": "key /definitely/missing" if fail else "key /"}
    wf.pipeline["multiqc"].update(switch=experiment_step, structure="echo summary")
    wf.run_order = list(steps) + (["multiqc"] if experiment_step else [])
    ws.write_workflow(wf, "t")
    fq = tmp_path / "fq"
    fq.mkdir()
    (fq / "s1.fastq").write_text("@r\nACGT\n+\nIIII\n")
    make_tsv(fq, "Homo sapiens", tmp_path / "exp.tsv")
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath="/x.fa"))
    return ws


def build_and_run(ws, tmp_path, scheduler="local"):
    assert main(["-w", str(ws.root), "build", str(tmp_path / "exp.tsv"), "-s", scheduler, "-p", "t"]) == 0
    if scheduler == "local":
        subprocess.run(["bash", str(ws.root / "master.sh")], check=True, capture_output=True)


def states(ws, query=lambda s, i: None):
    return {j.name: (j.state, j.step, j.exit_code) for j in mon.collect(ws, ["exp"], mon.ErrorScanner(), query)}


# --- shells other than bash: tracking, fail-fast, dependencies ------------------------
@pytest.mark.parametrize("shell", SHELLS)
def test_every_supported_shell_reports_failures_and_skips_dependants(tmp_path, shell):
    ws = workspace_with(tmp_path, shell)
    build_and_run(ws, tmp_path)
    got = states(ws)
    assert got["s_s1_0"][0] == mon.FAILED and got["s_s1_0"][1] == "boom" and got["s_s1_0"][2] not in (None, 0)
    assert got["e_exp"][0] == mon.SKIPPED


@pytest.mark.parametrize("shell", SHELLS)
def test_every_supported_shell_reports_success(tmp_path, shell):
    ws = workspace_with(tmp_path, shell, fail=False)
    build_and_run(ws, tmp_path)
    got = states(ws)
    assert got["s_s1_0"][0] == mon.DONE and got["s_s1_0"][1] == "boom"
    assert got["e_exp"][0] == mon.DONE
    assert (ws.output_dir / "exp" / "fastq" / "s1_0_HS.fastq.gz").exists()


def test_untracked_shell_still_gets_its_outcome_from_the_local_scheduler(tmp_path):
    ws = workspace_with(tmp_path, "bash")
    wf = ws.load_workflow("t")
    wf.set_general("environment", "bash")
    ws.write_workflow(wf, "t")
    build_and_run(ws, tmp_path)
    exp = ws.output_dir / "exp"
    # Simulate a shell without tracking: remove what the job wrote and let the local submit record it.
    (exp / "status" / "s_s1_0.status").unlink()
    script = exp / "scripts" / "s_s1_0.sh"
    subprocess.run(["bash", "-c", (ws.root / "master.sh").read_text().split("# ---- exp")[0]
                    + f'\nsubmit {script} > /dev/null'], check=True, capture_output=True)
    assert states(ws)["s_s1_0"][0] == mon.FAILED


@pytest.mark.skipif(not shutil.which("tcsh"), reason="tcsh not installed")
def test_reference_index_commands_use_csh_syntax_under_tcsh(tmp_path):
    ws = workspace_with(tmp_path, "tcsh")
    fasta = tmp_path / "mt.fa"
    fasta.write_text(">MT\nACGT\n")
    for ext in (".fai", ".bwt"):
        (tmp_path / f"mt.fa{ext}").write_text("x")          # already indexed: the guards must skip the tools
    wf = ws.load_workflow("t")
    wf.run_order = ["bwa-align"]
    cmds = index_commands(ReferenceSet(refpath=str(fasta), source="ncbi_mt"), wf)
    assert cmds and all(c.startswith("if ( ! -e") for c in cmds)
    script = tmp_path / "prep.csh"
    script.write_text("\n".join(cmds) + "\necho ok\n")
    result = subprocess.run(["tcsh", str(script)], capture_output=True, text=True)
    assert result.returncode == 0 and result.stdout.strip() == "ok", result.stderr


# --- SGE: dependants must not run after a failed dependency --------------------------
@pytest.mark.parametrize("shell", SHELLS)
def test_dependants_skip_themselves_when_a_dependency_failed(tmp_path, shell):
    ws = workspace_with(tmp_path, shell, fail=False)
    build_and_run(ws, tmp_path, scheduler="sge")
    exp = ws.output_dir / "exp"
    (exp / "status").mkdir(exist_ok=True)
    run = lambda name: subprocess.run([shell, str(exp / "scripts" / f"{name}.sh")], capture_output=True, text=True)  # noqa: E731
    # SGE's -hold_jid releases e_exp once s_s1_0 has ended, whatever its exit code:
    (exp / "status" / "s_s1_0.status").write_text("FAILED boom 1 2 1\n")
    assert run("e_exp").returncode == 0
    assert (exp / "status" / "e_exp.status").read_text().startswith("SKIPPED s_s1_0")
    (exp / "status" / "s_s1_0.status").write_text("DONE boom 1 2 0\n")
    assert run("e_exp").returncode == 0
    assert (exp / "status" / "e_exp.status").read_text().startswith("DONE")


# --- monitor ---------------------------------------------------------------------------
def fake_exp(ws, tracked=True, status="RUNNING count 1 2 -"):
    exp = ws.experiment_dir("E")
    (exp / "status").mkdir(parents=True)
    (exp / "pipemaker_manifest.json").write_text(json.dumps({"scheduler": "sge", "status_tracking": tracked}))
    (exp / "jobs.tsv").write_text("job_id\tjob_name\tscript\tsubmitted\tkind\n11\ts_a_0\t/a.sh\t1\tsample\n")
    if status:
        (exp / "status" / "s_a_0.status").write_text(status + "\n")
    return exp


def test_errors_watch_keeps_watching_while_jobs_run(tmp_path, monkeypatch):
    ws = Workspace(tmp_path / "ws")
    fake_exp(ws)
    sleeps = []

    def fake_sleep(seconds):
        sleeps.append(seconds)
        raise KeyboardInterrupt  # stop after the first refresh

    monkeypatch.setattr(mon.time, "sleep", fake_sleep)
    out = io.StringIO()
    mon.monitor(ws, ["E"], watch=5, errors_only=True, stream=out, color=False, query=lambda s, i: {"11": (mon.RUNNING, "r")})
    assert sleeps == [5] and "All jobs have finished" not in out.getvalue()


def test_finished_untracked_jobs_are_not_rerun(tmp_path):
    ws = Workspace(tmp_path / "ws")
    fake_exp(ws, tracked=False, status=None)
    gone = lambda s, i: {}  # noqa: E731  - qstat no longer lists the job
    (job,) = mon.collect(ws, ["E"], mon.ErrorScanner(), gone)
    assert job.state == mon.ENDED and job.state in mon.FINAL_STATES
    assert plan(ws, ["E"], query=gone) == []
    # With tracking, a job that vanished without ever writing its status did not run: rerun it.
    ws2 = Workspace(tmp_path / "ws2")
    fake_exp(ws2, tracked=True, status=None)
    (job,) = mon.collect(ws2, ["E"], mon.ErrorScanner(), gone)
    assert job.state == mon.UNKNOWN


def test_error_scanner_restarts_on_a_rewritten_log(tmp_path):
    log = tmp_path / "job.err"
    log.write_text("ERROR: first run failed\n" + "progress\n" * 5)
    scanner = mon.ErrorScanner()
    assert scanner.scan(log)[0] == 1
    log.write_text("rerun is fine\n" + "progress line of the rerun\n" * 20)  # same file, now longer
    assert scanner.scan(log) == (0, [])


# --- build / rerun -------------------------------------------------------------------
def test_build_with_only_an_experiment_level_step(tmp_path, capsys):
    ws = workspace_with(tmp_path)
    assert main(["-w", str(ws.root), "build", str(tmp_path / "exp.tsv"), "-p", "t", "--steps", "multiqc"]) == 0
    assert (ws.root / "master.sh").exists()
    assert "1 job script(s)" in capsys.readouterr().out


def test_a_suggested_resume_step_outside_the_build_does_not_block_the_rerun(tmp_path, monkeypatch):
    from pipemaker import rerun as rerun_mod
    from pipemaker.diagnose import Suggestion

    ws = workspace_with(tmp_path)
    build_and_run(ws, tmp_path)
    monkeypatch.setattr(rerun_mod, "diagnose", lambda job, ctx: [Suggestion("x", rerun={"from_step": "STAR-alignment", "clean": True})])
    items = {i.name: i for i in plan(ws, ["exp"], RerunOptions(auto_fix=True))}
    assert items["s_s1_0"].start_step == "boom" and not items["s_s1_0"].clean   # falls back to the failed step
    with pytest.raises(Exception, match="not part of this build"):
        plan(ws, ["exp"], RerunOptions(from_step="STAR-alignment"))           # an explicit --from still errors
