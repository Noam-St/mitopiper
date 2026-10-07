import io
import subprocess

import pytest

from mitopiper import monitor as mon
from mitopiper.cli import main
from mitopiper.config import Workspace
from mitopiper.diagnose import Context, diagnose, double_memory, double_walltime, merged_rerun_options
from mitopiper.references import ReferenceRegistry, ReferenceSet
from mitopiper.rerun import RerunJob, RerunOptions, plan, rerun, write


# ------------------------------------------------------------------ diagnosis
def job(errors=(), rc=None, sched="", step="bwa-align", state=mon.FAILED):
    return mon.JobStatus("E", "s_x_0", state=state, exit_code=rc, scheduler_state=sched, step=step,
                         errors=[mon.ErrorLine(i, t) for i, t in enumerate(errors, 1)])


def _ctx(steps=("wget", "fastqc_pre", "trim", "fastqc", "bwa-align", "view", "sort", "index", "mpileup")):
    from mitopiper.config import Workflow

    wf = Workflow.default()
    return Context(steps=list(steps), pipeline=wf.pipeline, memory="30G", walltime="1-20:00:00", scheduler="slurm",
                   exp_dir="/w/output/E", references={"Homo_sapiens": {"refpath": "/refs/hs.fa", "gtfpath": "/refs/hs.gtf"}})


CTX = _ctx()


@pytest.mark.parametrize("kwargs,cause,rerun_opts", [
    (dict(sched="OUT_OF_MEMORY", state=mon.KILLED), "out of memory", {"memory": "60G"}),
    (dict(rc=137), "out of memory", {"memory": "60G"}),
    (dict(errors=["slurmstepd: error: *** JOB 5 CANCELLED DUE TO TIME LIMIT ***"], state=mon.KILLED), "time limit",
     {"walltime": "3-16:00:00"}),
    (dict(errors=["s.sh: line 3: bwa: command not found"], rc=127), "'bwa' is not installed",
     {"free_form_hint": "module load bwa"}),
    (dict(errors=["gzip: /w/output/E/fastq/x_1.fastq.gz: unexpected end of file"], step="trim"), "corrupt",
     {"from_step": "wget", "clean": True}),
    (dict(errors=["gzip: /w/output/E/post_trim/x_1_val_1.fq.gz: invalid compressed data--crc error"]), "corrupt",
     {"from_step": "trim", "clean": True}),                        # not fastqc, which writes no data
    (dict(errors=["[E::bgzf_read_block] EOF marker is absent. The input /w/output/E/alignment/x.sorted.bam is probably truncated"],
          step="mpileup"), "damaged input from 'sort'", {"from_step": "sort", "clean": True}),   # not index (.bai)
    (dict(errors=["md5sum: /w/output/E/fastq/x_1.fastq.gz: FAILED", "md5sum: WARNING: 1 computed checksum did NOT match"],
          step="wget"), "failed checksum", {"from_step": "wget", "clean": True}),
    (dict(errors=["[E::bwa_idx_load_from_disk] fail to locate the index files"]), "BWA index", {}),
    (dict(errors=["[main_samview] fail to read the header from \"/w/output/E/alignment/x.sam\"."], step="view"),
     "empty or damaged input", {"from_step": "bwa-align", "clean": True}),
    (dict(errors=["[E::hts_open_format] fail to open file '/refs/hs.fa'", "/refs/hs.fa: No such file or directory"]),
     "reference file", {}),
    (dict(errors=["ls: cannot access '/w/output/E/post_trim/x_trimmed.fq.gz': No such file or directory"]),
     "which step 'trim' should have produced", {"from_step": "trim"}),       # not fastqc (same folder)
    (dict(errors=["[main_samview] fail to open \"/w/output/E/alignment/x.sam\": No such file or directory"], step="view"),
     "which step 'bwa-align' should have produced", {"from_step": "bwa-align"}),   # not index (same folder, later)
    (dict(errors=["gzip: /w/output/E/fastq/x_1.fastq.gz: No such file or directory"], step="trim"),
     "which step 'wget' should have produced", {"from_step": "wget"}),        # not fastqc_pre
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
    wf.set_general("tool_check", False)  # these tests are about failures inside the pipeline
    ws.write_workflow(wf, "demo")
    fq = tmp_path / "fq"
    fq.mkdir()
    for s in ("a", "b"):
        (fq / f"{s}.fastq").write_text("@r\nACGT\n+\nIIII\n")
    assert main(["-w", str(ws.root), "make-tsv", str(fq), "--org", "Homo sapiens", "-o", str(tmp_path / "demo.tsv")]) == 0
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath="/x.fa"))
    assert main(["-w", str(ws.root), "build", "--no-check", str(tmp_path / "demo.tsv"), "-s", "local", "-p", "demo"]) == 0
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
    manifest = (exp / "mitopiper_manifest.json").read_text().replace('"scheduler": "local"', '"scheduler": "slurm"')
    (exp / "mitopiper_manifest.json").write_text(manifest)
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
    assert "'mytool' is not installed" in text and "mitopiper monitor --rerun" in text
    out = io.StringIO()
    mon.monitor(ws, errors_only=True, stream=out, color=False)
    assert "suggested fix:" in out.getvalue() and "$ mitopiper rerun --free-form 'module load mytool'" in out.getvalue()


def test_monitor_rerun_shortcut(tmp_path, capfd):
    ws = make_ws(tmp_path)
    assert main(["-w", str(ws.root), "monitor", "--rerun"]) == 0
    assert "3 job(s) to rerun" in capfd.readouterr().out
    assert states(ws)["s_a_0"] == mon.FAILED  # still failing: the tool is still missing
    assert main(["-w", str(ws.root), "rerun", "--dry-run"]) == 0


def test_corrupt_input_from_outside_the_build_gets_no_automatic_resume_point():
    ctx = _ctx(steps=("trim", "fastqc", "bwa-align"))            # built with --steps, no download step
    (sug, *_) = diagnose(job(["gzip: /w/output/E/fastq/x_1.fastq.gz: unexpected end of file"], step="trim"), ctx)
    assert "not produced by any step of this build" in sug.advice[0]
    assert merged_rerun_options([sug]) == {}


def test_fastqc_truncated_file_is_corrupt_input():
    # FastQC is usually the first step to read a truncated FASTQ (seen on Grid Engine)
    ctx = _ctx(steps=("wget", "fastqc_pre", "trim"))
    lines = ["Failed to process file patient_1_SO_2.fastq.gz",
             "uk.ac.babraham.FastQC.Sequence.SequenceFormatException: Ran out of data in the middle of a fastq entry.  "
             "Your file is probably truncated"]
    (sug, *_) = diagnose(job(lines, step="fastqc_pre"), ctx)
    assert sug.cause.startswith("An input file is corrupt or incomplete")


def test_rerun_with_a_fixed_preset_checks_the_preset_tools(registered, data_dir):
    # The build used a tool path that does not exist; the check job of a rerun with the fixed
    # preset must look for the tools of that preset, not those recorded at build time.
    ws = registered
    wf = ws.load_workflow("default")
    wf.pipeline["mpileup"]["path"] = "/nonexistent/samtools"
    ws.write_workflow(wf, "default")
    assert main(["-w", str(ws.root), "build", "--no-check", str(data_dir / "ena_runs.tsv"), "-s", "sge"]) == 0
    check = ws.output_dir / "ena_runs" / "scripts" / "check_ena_runs.sh"
    assert "/nonexistent/samtools" in check.read_text()
    wf.pipeline["mpileup"]["path"] = ""
    ws.write_workflow(wf, "default")
    item = RerunJob("ena_runs", "check_ena_runs", "check", mon.FAILED, "check_tools")
    write(ws, [item], RerunOptions(preset="default"))
    assert "/nonexistent/samtools" not in check.read_text() and "samtools" in check.read_text()
    write(ws, [item], RerunOptions())                    # without --preset: as built
    assert "/nonexistent/samtools" in check.read_text()
