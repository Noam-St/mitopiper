import io
import json
import subprocess

from mitopiper import monitor as mon
from mitopiper.cli import main
from mitopiper.config import Workspace
from mitopiper.references import ReferenceRegistry, ReferenceSet


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
    (exp / "mitopiper_manifest.json").write_text(json.dumps({"scheduler": scheduler}))
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


def test_cutadapt_reports_are_not_errors(tmp_path):
    log = tmp_path / "s.err"
    log.write_text("Maximum trimming error rate: 0.1 (default)\nNo. of allowed errors:\n"
                   "length\tcount\texpect\tmax.err\terror counts\ncutadapt: error: unrecognized arguments\n")
    count, errors = mon.ErrorScanner().scan(log)
    assert count == 1 and "unrecognized" in errors[0].text


# qacct records of a real Son of Grid Engine (8.1.9, daimh/sge) run, shortened.
QACCT = {
    "4": "==============================================================\nqname        all.q\njobname      check_e\n"
         "jobnumber    4\nfailed       0    \nexit_status  0                   \nru_wallclock 6s\n"
         "category     -q all.q -l h_rt=900,mem_free=1G -pe shared 1\n",
    "27": "==============================================================\nqname        all.q\njobname      s_control_0\n"
          "jobnumber    27\nfailed       37  : qmaster enforced h_rt, h_cpu, or h_vmem limit\n"
          "exit_status  137                  (Killed)\nru_wallclock 11s\n"
          "category     -q all.q -l h_rt=10,mem_free=1G -pe shared 1\n",
    "30": "==============================================================\nqname        all.q\njobname      to_cancel\n"
          "jobnumber    30\nfailed       100 : assumedly after job\nexit_status  137                  (Killed)\n"
          "ru_wallclock 4s\ncategory     -l h_rt=600\n",
    "31": "==============================================================\njobnumber    31\n"
          "failed       37  : qmaster enforced h_rt, h_cpu, or h_vmem limit\nexit_status  137\nru_wallclock 50s\n"
          "category     -l h_rt=3600,h_vmem=2G\n",
}


def test_sge_accounting_reads_why_a_job_ended(monkeypatch):
    monkeypatch.setattr(mon, "_SGE_ACCOUNTING", {})
    calls = []

    def qacct(cmd):
        calls.append(cmd)
        return QACCT.get(cmd[2], "error: job id not found\n")

    assert mon.sge_accounting("4", qacct) == (mon.DONE, "COMPLETED")
    assert mon.sge_accounting("27", qacct) == (mon.KILLED, "TIMEOUT")          # ran its whole h_rt of 10 s
    assert mon.sge_accounting("30", qacct) == (mon.KILLED, "CANCELLED")        # qdel
    assert mon.sge_accounting("31", qacct) == (mon.KILLED, "OUT_OF_MEMORY")    # far below h_rt, asked h_vmem
    assert mon.sge_accounting("99", qacct) is None                             # not (yet) in the accounting
    mon.sge_accounting("27", qacct)
    assert [c[2] for c in calls] == ["4", "27", "30", "31", "99"]              # finished jobs are asked once
    assert mon.sge_accounting("12_3", lambda cmd: calls.append(cmd)) is None and calls[-1] == ["qacct", "-j", "12", "-t", "3"]


def test_bwa_info_lines_are_not_a_missing_index():
    from mitopiper.diagnose import diagnose
    job = mon.JobStatus(experiment="E", name="s_a_0", kind="sample", job_id="27", script="")
    job.state, job.step, job.scheduler_state = mon.KILLED, "bwa-align", "TIMEOUT"
    job.errors = [mon.ErrorLine(1, "[M::bwa_idx_load_from_disk] read 0 ALT contigs")] if hasattr(mon, "ErrorLine") else []
    from mitopiper.diagnose import Context
    causes = [s.cause for s in diagnose(job, Context(scheduler="sge"))]
    assert causes == ["The job hit its time limit in step 'bwa-align'"]        # no BWA index, no Eqw
    job.scheduler_state = "Eqw"
    assert any("Eqw" in s.cause for s in diagnose(job, Context(scheduler="sge")))


# qacct -j 17163 on the BGU cluster (N1GE 6.0u8): the job number was used twice before, by other people,
# and a time limit is recorded as failed 100 with no category line. Shortened.
N1GE_17163 = (
    "==============================================================\nqname        barlev.q\nhostname     bhn1080\n"
    "owner        someone\njobname      prf_v146_f\njobnumber    17163\nqsub_time    Thu May 21 10:00:00 2026\n"
    "failed       0    \nexit_status  0\nru_wallclock 1177\nmaxvmem      34.7G\n"
    "==============================================================\nqname        intel_all.q\nhostname     bhn1099\n"
    "owner        other\njobname      run.sh\njobnumber    17163\nqsub_time    Thu Jul 30 09:00:00 2026\n"
    "failed       0    \nexit_status  0\nru_wallclock 159\nmaxvmem      2.1G\n"
    "==============================================================\nqname        bioinfo.q\nhostname     bhn1109\n"
    "owner        labuser\njobname      s_control_0\njobnumber    17163\nqsub_time    Mon Oct 12 15:38:10 2026\n"
    "failed       100 : assumedly after job\nexit_status  137\nru_wallclock 6\nmaxvmem      4.4G\n"
)


def test_sge_accounting_on_n1ge_with_reused_job_numbers(monkeypatch, tmp_path):
    monkeypatch.setattr(mon, "_SGE_ACCOUNTING", {})
    submitted = mon._qacct_time("Mon Oct 12 15:38:09 2026")
    script = tmp_path / "s_control_0.sh"
    script.write_text("#!/bin/bash\n#$ -N s_control_0\n#$ -l mem_free=4G\n#$ -l h_rt=0:00:05\necho\n")
    assert mon.script_walltime(str(script)) == 5
    qacct = lambda cmd: N1GE_17163  # noqa: E731
    # killed by h_rt=5 s after 6 s: a time limit, not a cancellation
    assert mon.sge_accounting("17163", qacct, "s_control_0", submitted, requested=5) == (mon.KILLED, "TIMEOUT")
    # the same signature far below the limit is a qdel
    monkeypatch.setattr(mon, "_SGE_ACCOUNTING", {})
    assert mon.sge_accounting("17163", qacct, "s_control_0", submitted, requested=1800) == (mon.KILLED, "CANCELLED")
    # before today's record is written, the old records of other jobs must not be taken for it
    monkeypatch.setattr(mon, "_SGE_ACCOUNTING", {})
    old_only = N1GE_17163.rsplit("=" * 62, 1)[0]
    assert mon.sge_accounting("17163", lambda cmd: old_only, "s_control_0", submitted, requested=5) is None
    records = mon.own_records(mon.qacct_records(N1GE_17163), "s_control_0", submitted)
    assert [r["qname"] for r in records] == ["bioinfo.q"]


def test_download_retries_are_not_errors(tmp_path):
    """BGU round 5: a DONE job showed ERR 3 for the attempts mitopiper download had recovered from."""
    log = tmp_path / "s.err"
    log.write_text("[mitopiper] download 1/10: ftp.sra.ebi.ac.uk/x.gz\n"
                   "[mitopiper] checksum failed for /o/x.gz (aa, expected bb): deleting it\n"
                   "[mitopiper] wget failed (exit code 4)\n"
                   "[mitopiper] ERROR: download failed 10 times: ftp.sra.ebi.ac.uk/x.gz\n")
    count, lines = mon.ErrorScanner().scan(log)
    assert count == 1 and "download failed 10 times" in lines[0].text


def test_a_huge_log_is_searched_only_at_its_end(tmp_path):
    # BGU: 4 million identical error lines made every monitor refresh, menu and rerun take minutes
    log = tmp_path / "s.err"
    log.write_text("start\n" + "".join(f"[E::fai_retrieve] Failed to retrieve block {i}\n" for i in range(1, 2001)))
    scanner = mon.ErrorScanner()
    scanner.MAX_SCAN, scanner.BLOCK = 4000, 333                   # small, so that the skip and the blocks are used
    count, errors = scanner.scan(log)
    assert scanner.partial(log) and 0 < count < 2000
    assert errors[-1].line == 2001 and errors[-1].text.endswith("block 2000")     # line numbers still right
    assert [e.line for e in errors] == list(range(2001 - len(errors) + 1, 2002))
    with open(log, "a") as handle:
        handle.write("ERROR: the end\n")
    assert scanner.scan(log)[0] == count + 1 and scanner.scan(log)[1][-1].line == 2002     # then incremental
    job = mon.JobStatus("E", "s", error_count=count, errors_partial=True)
    assert mon._error_count(job) == f"{count}+"
    assert mon.ErrorScanner().scan(log)[0] == 2001 and not mon.ErrorScanner().partial(log)   # small: read in full
