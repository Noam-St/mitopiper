import io
import json
import subprocess

from pipemaker import monitor as mon
from pipemaker.cli import main
from pipemaker.config import Workspace
from pipemaker.references import ReferenceRegistry, ReferenceSet


def local_build(tmp_path, fail=True):
    """Build + run (local scheduler) one sample whose second step fails or succeeds."""
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("default")
    wf.pipeline["boom"] = {"key": "ls", "description": "list a file", "switch": True, "path": "", "outdir": "",
                           "output_ext": "", "inputdir": "", "input_ext": "", "scope": "sample",
                           "structure": "key /definitely/missing/file" if fail else "key /"}
    wf.run_order = ["wget", "boom"]
    wf.set_general("tool_check", False)  # these tests are about failures inside the pipeline
    ws.write_workflow(wf, "t")
    fq = tmp_path / "fq"
    fq.mkdir()
    (fq / "s1.fastq").write_text("@r\nACGT\n+\nIIII\n")
    assert main(["-w", str(ws.root), "make-tsv", str(fq), "--org", "Homo sapiens", "-o", str(tmp_path / "exp.tsv")]) == 0
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath="/x.fa"))
    assert main(["-w", str(ws.root), "build", "--no-check", str(tmp_path / "exp.tsv"), "-s", "local", "-p", "t"]) == 0
    subprocess.run(["bash", str(ws.root / "master.sh")], check=True, capture_output=True)
    return ws


def run_monitor(ws, *args):
    out = io.StringIO()
    rc = mon.monitor(ws, stream=out, color=False, **dict(args))
    return rc, out.getvalue()


def test_failed_job_is_reported_with_step_exit_code_and_stderr(tmp_path):
    ws = local_build(tmp_path, fail=True)
    rc, text = run_monitor(ws)
    assert rc == 1
    assert "s_s1_0" in text and "FAILED rc=2" in text and "boom" in text
    assert "ls: cannot access '/definitely/missing/file'" in text
    rc, detail = run_monitor(ws, ("errors_only", True))
    assert "<-- failed here" in detail and "/definitely/missing/file" in detail
    (job,) = json.loads(run_monitor(ws, ("as_json", True))[1])
    assert job["state"] == "FAILED" and job["step"] == "boom" and job["exit_code"] == 2
    assert [s[0] for s in job["steps"]] == ["wget", "boom"]


def test_successful_job_and_watch_stops_when_done(tmp_path):
    ws = local_build(tmp_path, fail=False)
    rc, text = run_monitor(ws, ("watch", 0.01))
    assert rc == 0 and "DONE" in text and "All jobs have finished." in text
    assert main(["-w", str(ws.root), "monitor", "--no-color"]) == 0


def make_exp(ws, name="E", scheduler="slurm"):
    exp = ws.experiment_dir(name)
    (exp / "status").mkdir(parents=True)
    (exp / "stderr").mkdir()
    (exp / "scripts").mkdir()
    (exp / "pipemaker_manifest.json").write_text(json.dumps({"scheduler": scheduler}))
    rows = ["job_id\tjob_name\tscript\tsubmitted\tkind",
            "11\ts_a_0\t/s/a.sh\t100\tsample", "12\ts_b_1\t/s/b.sh\t100\tsample",
            "13\ts_c_2\t/s/c.sh\t100\tsample", "14\ts_d_3\t/s/d.sh\t100\tsample"]
    (exp / "jobs.tsv").write_text("\n".join(rows) + "\n")
    (exp / "scripts" / "s_e_4.sh").write_text("")
    return exp


def test_states_combine_status_files_and_scheduler(tmp_path):
    ws = Workspace(tmp_path / "ws")
    exp = make_exp(ws)
    (exp / "status" / "s_a_0.status").write_text("RUNNING align 100 200 -\n")   # scheduler: TIMEOUT
    (exp / "status" / "s_b_1.status").write_text("RUNNING trim 100 200 -\n")    # scheduler: RUNNING
    (exp / "status" / "s_c_2.status").write_text("DONE mpileup 100 300 0\n")
    (exp / "stderr" / "s_a_0_11.err").write_text("ok\nslurmstepd: error: *** JOB 11 CANCELLED DUE TO TIME LIMIT ***\n")

    def fake_query(scheduler, ids):
        assert scheduler == "slurm" and sorted(ids) == ["11", "12", "13", "14"]
        return {"11": (mon.KILLED, "TIMEOUT"), "12": (mon.RUNNING, "RUNNING"), "14": (mon.PENDING, "PENDING")}

    jobs = {j.name: j for j in mon.collect(ws, ["E"], mon.ErrorScanner(), fake_query)}
    assert jobs["s_a_0"].state == mon.KILLED and jobs["s_a_0"].scheduler_state == "TIMEOUT"
    assert jobs["s_a_0"].error_count == 1 and "TIME LIMIT" in jobs["s_a_0"].last_error
    assert jobs["s_b_1"].state == mon.RUNNING and jobs["s_b_1"].step == "trim"
    assert jobs["s_c_2"].state == mon.DONE and jobs["s_c_2"].elapsed == 200
    assert jobs["s_d_3"].state == mon.PENDING
    assert jobs["s_e_4"].state == mon.NOT_SUBMITTED
    # Without scheduler access a job that reported RUNNING is trusted.
    jobs = {j.name: j for j in mon.collect(ws, ["E"], mon.ErrorScanner(), lambda s, i: None)}
    assert jobs["s_a_0"].state == mon.RUNNING and jobs["s_d_3"].state == mon.PENDING


def test_find_log_for_each_scheduler(tmp_path):
    exp = tmp_path
    (exp / "stderr").mkdir()
    for name in ("s_a_0.e11", "s_b_1_12.err", "13.err", "s_d_3.err"):
        (exp / "stderr" / name).write_text("")
    assert mon.find_log(exp, "stderr", "s_a_0", "11").name == "s_a_0.e11"      # SGE
    assert mon.find_log(exp, "stderr", "s_b_1", "12").name == "s_b_1_12.err"   # SLURM
    assert mon.find_log(exp, "stderr", "s_c_2", "13").name == "13.err"         # LSF
    assert mon.find_log(exp, "stderr", "s_d_3", "s_d_3").name == "s_d_3.err"   # local
    assert mon.find_log(exp, "stderr", "zzz", "99") is None


def test_error_scanner_is_incremental_and_ignores_benign_lines(tmp_path):
    log = tmp_path / "x.err"
    log.write_text("Processing reads\n0 errors found\nERROR: truncated gzip")
    scanner = mon.ErrorScanner()
    assert scanner.scan(log)[0] == 0              # the error line is not complete yet
    with open(log, "a") as handle:
        handle.write(" input\nsamtools: Segmentation fault\n")
    count, errors = scanner.scan(log)
    assert count == 2 and [e.line for e in errors] == [3, 4]
    assert errors[0].text == "ERROR: truncated gzip input"
    assert mon.ErrorScanner(ignore=()).scan(log)[0] == 3


def test_query_scheduler_parsing():
    outputs = {
        "sacct": "11|TIMEOUT\n12|RUNNING\n13|CANCELLED by 42\n",
        "qstat": "job-ID prior name user state\n----\n  21 0.5 s_a u r 01/01\n  22 0.5 s_b u qw 01/01\n  23 0.5 s_c u Eqw 01/01\n",
        "bjobs": "31|RUN\n32|EXIT\n33|PEND\n",
    }
    runner = lambda cmd: outputs.get(cmd[0])  # noqa: E731
    assert mon.query_scheduler("slurm", ["11", "12", "13"], runner) == {
        "11": (mon.KILLED, "TIMEOUT"), "12": (mon.RUNNING, "RUNNING"), "13": (mon.KILLED, "CANCELLED")}
    assert mon.query_scheduler("sge", ["21", "22", "23"], runner) == {
        "21": (mon.RUNNING, "r"), "22": (mon.PENDING, "qw"), "23": (mon.FAILED, "Eqw")}
    assert mon.query_scheduler("lsf", ["31", "32", "33"], runner) == {
        "31": (mon.RUNNING, "RUN"), "32": (mon.FAILED, "EXIT"), "33": (mon.PENDING, "PEND")}
    assert mon.query_scheduler("slurm", ["1"], lambda cmd: None) is None
    assert mon.query_scheduler("local", ["1"], runner) is None


def test_default_experiments_come_from_last_build(tmp_path):
    ws = local_build(tmp_path)
    assert mon.default_experiments(ws) == ["exp"]
    (ws.root / "last_build.json").unlink()
    assert mon.default_experiments(ws) == ["exp"]


def test_last_error_strips_bash_script_prefix():
    job = mon.JobStatus("E", "j", errors=[mon.ErrorLine(3, "/x/scripts/s_a_0.sh: line 26: fastqc: command not found")])
    assert job.last_error == "line 26: fastqc: command not found"


def test_jobs_not_yet_submitted_by_a_running_master_are_pending(tmp_path):
    ws = Workspace(tmp_path / "ws")
    exp = make_exp(ws)
    (exp / ".submitting").write_text("")
    jobs = {j.name: j for j in mon.collect(ws, ["E"], mon.ErrorScanner(), lambda s, i: None)}
    assert jobs["s_e_4"].state == mon.PENDING and jobs["s_e_4"].scheduler_state == "waiting for master.sh"
