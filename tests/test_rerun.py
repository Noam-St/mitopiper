import io
import subprocess

import pytest

from pipemaker import monitor as mon
from pipemaker.cli import main
from pipemaker.config import Workspace
from pipemaker.diagnose import Context, diagnose, double_memory, double_walltime, merged_rerun_options
from pipemaker.references import ReferenceRegistry, ReferenceSet
from pipemaker.rerun import RerunOptions, plan, rerun, write


# ------------------------------------------------------------------ diagnosis
def job(errors=(), rc=None, sched="", step="bwa-align", state=mon.FAILED):
    return mon.JobStatus("E", "s_x_0", state=state, exit_code=rc, scheduler_state=sched, step=step,
                         errors=[mon.ErrorLine(i, t) for i, t in enumerate(errors, 1)])


CTX = Context(steps=["wget", "trim", "bwa-align", "view"], memory="30G", walltime="1-20:00:00", scheduler="slurm",
              exp_dir="/w/output/E",
              pipeline={"wget": {"outdir": "fastq"}, "trim": {"outdir": "post_trim"},
                        "bwa-align": {"outdir": "alignment", "key": "bwa mem"}, "view": {"outdir": "alignment"}},
              references={"Homo_sapiens": {"refpath": "/refs/hs.fa", "gtfpath": "/refs/hs.gtf"}})


@pytest.mark.parametrize("kwargs,cause,rerun_opts", [
    (dict(sched="OUT_OF_MEMORY", state=mon.KILLED), "out of memory", {"memory": "60G"}),
    (dict(rc=137), "out of memory", {"memory": "60G"}),
    (dict(errors=["slurmstepd: error: *** JOB 5 CANCELLED DUE TO TIME LIMIT ***"], state=mon.KILLED), "time limit",
     {"walltime": "3-16:00:00"}),
    (dict(errors=["s.sh: line 3: bwa: command not found"], rc=127), "'bwa' is not installed",
     {"free_form_hint": "module load bwa"}),
    (dict(errors=["gzip: /w/output/E/fastq/x_1.fastq.gz: unexpected end of file"], step="trim"), "corrupt",
     {"from_step": "wget", "clean": True}),
    (dict(errors=["[E::bwa_idx_load_from_disk] fail to locate the index files"]), "BWA index", {}),
    (dict(errors=["[main_samview] fail to read the header from \"/w/output/E/alignment/x.sam\"."], step="view"),
     "empty or damaged input", {"from_step": "bwa-align", "clean": True}),
    (dict(errors=["[E::hts_open_format] fail to open file '/refs/hs.fa'", "/refs/hs.fa: No such file or directory"]),
     "reference file", {}),
    (dict(errors=["ls: cannot access '/w/output/E/post_trim/x_trimmed.fq.gz': No such file or directory"]),
     "which step 'trim' should have produced", {"from_step": "trim"}),
    (dict(errors=["wget: unable to resolve host address 'ftp.sra.ebi.ac.uk'"], step="wget"), "no network access", {}),
    (dict(errors=["write error: No space left on device"]), "disk", {}),
    (dict(rc=139), "segmentation fault", {"memory": "60G"}),
    (dict(rc=3), "Unrecognised failure", {}),
])
def test_diagnosis_rules(kwargs, cause, rerun_opts):
    suggestions = diagnose(job(**kwargs), CTX)
    assert cause.lower() in suggestions[0].cause.lower(), [s.cause for s in suggestions]
    assert suggestions[0].commands
    assert merged_rerun_options(suggestions[:1]) == rerun_opts


def test_resource_doubling():
    assert double_memory("30G") == "60G"
    assert double_memory("600M") == "1.2G"
    assert double_memory("") == "8G"
    assert double_walltime("1-20:00:00") == "3-16:00:00"
    assert double_walltime("04:00:00") == "08:00:00"
    assert double_walltime("90") == "03:00:00"
    assert double_walltime("junk") == "1-00:00:00"


# ------------------------------------------------------------------- rerun
def make_ws(tmp_path, fail_cmd="mytool"):
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("default")
    base = {"switch": True, "path": "", "outdir": "counts", "inputdir": "fastq", "scope": "sample"}
    wf.pipeline["count"] = dict(base, key=f"{fail_cmd} count", description="count", output_ext=".txt",
                                input_ext=".fastq.gz", structure="key input output_pipe")
    wf.pipeline["multiqc"]["switch"] = True
    wf.pipeline["multiqc"]["structure"] = "echo summary"
    wf.run_order = ["wget", "count", "multiqc"]
    ws.write_workflow(wf, "demo")
    fq = tmp_path / "fq"
    fq.mkdir()
    for s in ("a", "b"):
        (fq / f"{s}.fastq").write_text("@r\nACGT\n+\nIIII\n")
    assert main(["-w", str(ws.root), "make-tsv", str(fq), "--org", "Homo sapiens", "-o", str(tmp_path / "demo.tsv")]) == 0
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath="/x.fa"))
    assert main(["-w", str(ws.root), "build", str(tmp_path / "demo.tsv"), "-s", "local", "-p", "demo"]) == 0
    subprocess.run(["bash", str(ws.root / "master.sh")], check=True, capture_output=True)
    return ws


def states(ws):
    return {j.name: j.state for j in mon.collect(ws, ["demo"], mon.ErrorScanner(), lambda s, i: None)}


def test_rerun_resumes_at_failed_step_with_fix(tmp_path):
    ws = make_ws(tmp_path)
    assert states(ws) == {"s_a_0": mon.FAILED, "s_b_1": mon.FAILED, "e_demo": mon.SKIPPED}
    items = plan(ws)
    assert [(i.name, i.start_step) for i in items] == [("s_a_0", "count"), ("s_b_1", "count"), ("e_demo", "")]
    out = io.StringIO()
    rc = rerun(ws, options=RerunOptions(free_form=["mytool() { zcat \"$2\" | wc -l; }"]), submit=True, stream=out)
    assert rc == 0
    script = (ws.output_dir / "demo" / "scripts" / "s_a_0.sh").read_text()
    assert "pm_step count" in script and "pm_step wget" not in script  # finished steps are not repeated
    assert states(ws) == {"s_a_0": mon.DONE, "s_b_1": mon.DONE, "e_demo": mon.DONE}
    assert (ws.output_dir / "demo" / "counts" / "a_0_HS.txt").read_text().strip() == "4"
    assert list((ws.output_dir / "demo" / "manifests").glob("rerun_*.json"))
    assert rerun(ws, stream=io.StringIO()) == 0  # nothing left to rerun


def test_rerun_options_from_start_clean_and_overrides(tmp_path):
    ws = make_ws(tmp_path)
    opts = RerunOptions(jobs=["s_a_0"], from_start=True, clean=True, memory="64G", walltime="2-00:00:00")
    items = plan(ws, options=opts)
    assert [(i.name, i.start_step) for i in items] == [("s_a_0", "wget"), ("e_demo", "")]
    write(ws, items, RerunOptions(scheduler="slurm"))
    script = (ws.output_dir / "demo" / "scripts" / "s_a_0.sh").read_text()
    assert "#SBATCH --mem=64G" in script and "#SBATCH --time=2-00:00:00" in script
    assert "rm -f" in script and "a_0_HS.fastq.gz" in script
    assert "pm_step wget" in script
    master = (ws.root / "rerun.sh").read_text()
    assert "sbatch --parsable" in master and "status/s_a_0.status" in master
    subprocess.run(["bash", "-n", str(ws.root / "rerun.sh")], check=True)
    with pytest.raises(Exception, match="not part of this build"):
        plan(ws, options=RerunOptions(from_step="STAR-alignment"))


def test_auto_fix_applies_per_job_suggestions_and_cancels_stuck_jobs(tmp_path):
    ws = make_ws(tmp_path)
    exp = ws.output_dir / "demo"
    # pretend this was a SLURM run: a was OOM-killed mid-step, b failed normally, the experiment job is stuck
    rows = (exp / "jobs.tsv").read_text().replace("s_a_0\ts_a_0", "11\ts_a_0").replace("s_b_1\ts_b_1", "12\ts_b_1") \
        .replace("e_demo\te_demo", "13\te_demo")
    (exp / "jobs.tsv").write_text(rows)
    (exp / "status" / "s_a_0.status").write_text("RUNNING count 1 2 -\n")
    (exp / "status" / "e_demo.status").unlink()
    manifest = (exp / "pipemaker_manifest.json").read_text().replace('"scheduler": "local"', '"scheduler": "slurm"')
    (exp / "pipemaker_manifest.json").write_text(manifest)
    query = lambda s, ids: {"11": (mon.KILLED, "OUT_OF_MEMORY"), "12": (mon.FAILED, "FAILED"), "13": (mon.PENDING, "PENDING")}  # noqa: E731
    items = {i.name: i for i in plan(ws, options=RerunOptions(auto_fix=True), query=query)}
    assert items["s_a_0"].overrides == {"memory": "60G"}
    assert items["s_b_1"].overrides == {}
    write(ws, list(items.values()))
    assert "#SBATCH --mem=60G" in (exp / "scripts" / "s_a_0.sh").read_text()
    assert "#SBATCH --mem=30G" in (exp / "scripts" / "s_b_1.sh").read_text()
    assert "scancel 13 || true" in (ws.root / "rerun.sh").read_text()


def test_monitor_shows_diagnosis_and_shortcuts(tmp_path):
    ws = make_ws(tmp_path)
    out = io.StringIO()
    assert mon.monitor(ws, stream=out, color=False) == 1
    text = out.getvalue()
    assert "'mytool' is not installed" in text and "python monitor.py --rerun" in text
    out = io.StringIO()
    mon.monitor(ws, errors_only=True, stream=out, color=False)
    assert "suggested fix:" in out.getvalue() and "$ pipemaker rerun --free-form 'module load mytool'" in out.getvalue()


def test_monitor_rerun_shortcut(tmp_path, capfd):
    ws = make_ws(tmp_path)
    assert main(["-w", str(ws.root), "monitor", "--rerun"]) == 0
    assert "3 job(s) to rerun" in capfd.readouterr().out
    assert states(ws)["s_a_0"] == mon.FAILED  # still failing: the tool is still missing
    assert main(["-w", str(ws.root), "rerun", "--dry-run"]) == 0
