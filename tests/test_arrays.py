import os
import subprocess

from mitopiper import monitor as mon
from mitopiper.cli import main
from mitopiper.config import Workspace
from mitopiper.local_samples import make_tsv
from mitopiper.references import ReferenceRegistry, ReferenceSet
from mitopiper.schedulers import JobResources, header

FAKE_SBATCH = r"""#!/bin/bash
# fake sbatch: runs the job (array tasks one by one) right away and prints a job id
script="${@: -1}"
id=$(( $(cat "$FAKE_STATE/id" 2>/dev/null || echo 100) + 1 )); echo $id > "$FAKE_STATE/id"
echo "$@" >> "$FAKE_STATE/calls"
range=$(grep -m1 -- '--array=' "$script" | sed -E 's/.*--array=1-([0-9]+).*/\1/')
if [ -n "$range" ]; then
    for i in $(seq 1 "$range"); do SLURM_ARRAY_TASK_ID=$i bash "$script" >/dev/null 2>&1; done
else
    bash "$script" >/dev/null 2>&1
fi
echo $id
"""


def test_array_directives():
    res = JobResources(threads=2)
    assert "#SBATCH --array=1-40%5" in header("slurm", "a_E", "/e", res, (40, 5))
    assert "#SBATCH --array=1-40\n" in header("slurm", "a_E", "/e", res, (40, 0))
    sge = header("sge", "a_E", "/e", res, (40, 5))
    assert "#$ -t 1-40" in sge and "#$ -tc 5" in sge
    assert '#BSUB -J "a_E[1-40]%5"' in header("lsf", "a_E", "/e", res, (40, 5))
    assert "--array" not in header("slurm", "s_1", "/e", res)


def test_task_state_falls_back_to_the_array():
    info = {"7_1": (mon.DONE, "COMPLETED"), "7_*": (mon.PENDING, "PENDING"), "9": (mon.RUNNING, "r")}
    assert mon.task_state(info, "7_1") == (mon.DONE, "COMPLETED")
    assert mon.task_state(info, "7_3") == (mon.PENDING, "PENDING")
    assert mon.task_state(info, "9_2") == (mon.RUNNING, "r")
    assert mon.task_state(info, "8") is None
    runner = lambda cmd: "7_1|COMPLETED\n7_[2-4%2]|PENDING\n" if cmd[0] == "sacct" else None  # noqa: E731
    assert mon.query_scheduler("slurm", ["7_1", "7_2"], runner) == {"7_1": (mon.DONE, "COMPLETED"), "7_*": (mon.PENDING, "PENDING")}


def test_samples_run_as_one_limited_array_job(tmp_path, monkeypatch):
    bin_dir, state = tmp_path / "bin", tmp_path / "state"
    bin_dir.mkdir()
    state.mkdir()
    (bin_dir / "sbatch").write_text(FAKE_SBATCH)
    (bin_dir / "sbatch").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_STATE", str(state))
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("default")
    wf.set_general("max_running", 2)
    wf.pipeline["multiqc"].update(switch=True, structure="echo summary")
    wf.run_order = ["wget", "multiqc"]
    ws.write_workflow(wf, "t")
    fq = tmp_path / "fq"
    fq.mkdir()
    for s in ("a", "b", "c"):
        (fq / f"{s}.fastq").write_text("@r\nACGT\n+\nIIII\n")
    make_tsv(fq, "Homo sapiens", tmp_path / "exp.tsv")
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath="/x.fa"))
    assert main(["-w", str(ws.root), "build", "--no-check", str(tmp_path / "exp.tsv"), "-s", "slurm", "-p", "t"]) == 0
    scripts = ws.output_dir / "exp" / "scripts"
    assert "#SBATCH --array=1-3%2" in (scripts / "a_exp.array").read_text()
    assert len((scripts / "a_exp.list").read_text().split()) == 3
    subprocess.run(["bash", str(ws.root / "master.sh")], check=True, capture_output=True)
    calls = (state / "calls").read_text().splitlines()
    assert len(calls) == 3                                    # check job, one array job for 3 samples, experiment job
    assert calls[1].endswith("a_exp.array")
    jobs = {j.name: j for j in mon.collect(ws, ["exp"], mon.ErrorScanner(), lambda s, i: None)}
    assert {n: j.state for n, j in jobs.items()} == {"check_exp": mon.DONE, "s_a_0": mon.DONE, "s_b_1": mon.DONE, "s_c_2": mon.DONE,
                                                           "e_exp": mon.DONE}
    assert jobs["s_b_1"].job_id.endswith("_2") and jobs["s_b_1"].stderr.endswith("s_b_1.err")
    assert "afterok:" in calls[-1]                            # the experiment job waits for the array
