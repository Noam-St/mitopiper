import os
import subprocess

import pytest

from mitopiper.schedulers import ExperimentJobs, JobResources, header, master_script, normalise


def test_normalise_aliases():
    assert normalise("qsub") == "sge"
    assert normalise("SBATCH") == "slurm"
    assert normalise("bsub") == "lsf"
    with pytest.raises(ValueError):
        normalise("pbs")


def test_headers():
    res = JobResources(queue="q", threads=4, memory="8G", walltime="10:00:00", account="lab")
    sge = header("sge", "job", "/exp", res)
    assert "#$ -pe shared 4" in sge and "#$ -q q" in sge and "#$ -l mem_free=8G" in sge
    slurm = header("slurm", "job", "/exp", res)
    assert "#SBATCH --account=lab" in slurm and "#SBATCH --time=10:00:00" in slurm
    assert "#SBATCH --output=/exp/stdout/%x_%j.out" in slurm
    assert "#$ -l h_rt=10:00:00" in sge
    lsf = header("lsf", "job", "/exp", res)
    assert "#BSUB -P lab" in lsf and '#BSUB -R "rusage[mem=8G]"' in lsf and "#BSUB -W 10:00" in lsf
    long = JobResources(walltime="1-20:00:00")
    assert "#$ -l h_rt=44:00:00" in header("sge", "j", "/e", long) and "#BSUB -W 44:00" in header("lsf", "j", "/e", long)
    assert "h_rt" not in header("sge", "j", "/e", JobResources(walltime=""))
    assert header("local", "job", "/exp", res) == "#!/bin/bash\n"


def test_no_account_or_queue_when_empty():
    slurm = header("slurm", "job", "/exp", JobResources())
    assert "--account" not in slurm and "--partition" not in slurm


@pytest.mark.parametrize("scheduler,needle", [("sge", "qsub -hold_jid"), ("slurm", "sbatch --parsable"), ("lsf", "bsub")])
def test_master_script_dependencies(scheduler, needle):
    text = master_script(scheduler, [ExperimentJobs("E", "/e/scripts/prep.sh", ["/e/scripts/a.sh"], "/e/scripts/exp.sh")])
    assert needle in text
    assert 'jid=$(submit /e/scripts/prep.sh ${prep[@]+"${prep[@]}"})' in text and 'prep=("$jid")' in text
    assert 'jid=$(submit /e/scripts/exp.sh ${jobs[@]+"${jobs[@]}"})' in text
    # every submission is recorded for the monitor
    assert "jobs_file=/e/jobs.tsv" in text and "marker=/e/.submitting" in text
    assert 'record "$jid" "$script" sample' in text and "record \"$jid\" /e/scripts/exp.sh experiment" in text


def test_master_script_submits_the_tool_check_first():
    text = master_script("slurm", [ExperimentJobs("E", "/e/scripts/prep.sh", ["/e/scripts/a.sh"], "/e/scripts/e.sh",
                                                  "/e", check="/e/scripts/check.sh")])
    order = [text.index(x) for x in ("submit /e/scripts/check.sh", "submit /e/scripts/prep.sh", "for script in", "submit /e/scripts/e.sh")]
    assert order == sorted(order)
    assert "record \"$jid\" /e/scripts/check.sh check" in text
    # without a reference job the samples wait for the check directly
    text = master_script("slurm", [ExperimentJobs("E", None, ["/e/scripts/a.sh"], None, "/e", check="/e/scripts/check.sh")])
    assert 'prep=("$jid")' in text and 'submit "$script" ${prep[@]+"${prep[@]}"}' in text


def test_stage_helpers():
    from mitopiper.schedulers import JobResources, check_resources, job_name, kind_of, stage_dependencies

    assert [kind_of(job_name(k, "E")) for k in ("check", "prepare", "experiment")] == ["check", "prepare", "experiment"]
    assert kind_of("s_SAMN1_0") == "sample"
    assert stage_dependencies("c", "p", ("a", "b")) == {"check": [], "prepare": ["c"], "download": ["p"], "sample": ["p"],
                                                         "experiment": ["a", "b"]}
    assert stage_dependencies("c", None, ()) == {"check": [], "prepare": ["c"], "download": ["c"], "sample": ["c"],
                                                 "experiment": ["c"]}
    assert stage_dependencies() == {"check": [], "prepare": [], "download": [], "sample": [], "experiment": []}
    assert stage_dependencies("c", "p", ("a",), download="d")["sample"] == ["d"]
    small = check_resources(JobResources(queue="q", threads=16, memory="64G", walltime="3-00:00:00"))
    assert (small.queue, small.threads, small.memory, small.walltime) == ("q", 1, "1G", "00:15:00")
    jobs = ExperimentJobs("E", "/p.sh", ["/a.sh"], "/e.sh", check="/c.sh")
    assert jobs.scripts() == ["/c.sh", "/p.sh", "/a.sh", "/e.sh"]


FAKE_OLD_QSUB = r"""#!/bin/sh
# An old Grid Engine: no -terse, the usual answer on stdout.
for a in "$@"; do [ "$a" = -terse ] && { echo 'qsub: ERROR! invalid option argument "-terse"' >&2; exit 1; }; done
echo "$@" >> "$QSUB_LOG"
n=$(($(wc -l < "$QSUB_LOG") + 100))
case "$*" in *a_*) echo "Your job-array $n.1-2:1 (\"x\") has been submitted" ;;
             *) echo "Your job $n (\"x\") has been submitted" ;; esac
"""


def test_sge_master_reads_job_ids_without_terse(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "qsub").write_text(FAKE_OLD_QSUB)
    (bin_dir / "qsub").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    monkeypatch.setenv("QSUB_LOG", str(tmp_path / "qsub.log"))
    from mitopiper.schedulers import _SUBMIT_FUNCTIONS
    script = _SUBMIT_FUNCTIONS["sge"] + "\n" + "\n".join([
        'a=$(submit /x/check_e.sh)', 'b=$(submit /x/s_1.sh "$a")', 'c=$(submit /x/a_e.array "$a" "$b")',
        'echo "$a $b $c"'])
    out = subprocess.run(["bash", "-euo", "pipefail", "-c", script], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.split() == ["101", "102", "103"]
    log = (tmp_path / "qsub.log").read_text().splitlines()
    assert log[1] == "-hold_jid 101 /x/s_1.sh" and log[2] == "-hold_jid 101,102 /x/a_e.array"


def test_sge_submit_fails_loudly_when_no_id(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "qsub").write_text("#!/bin/sh\necho 'Unable to run job: denied'\n")
    (bin_dir / "qsub").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    from mitopiper.schedulers import _SUBMIT_FUNCTIONS
    out = subprocess.run(["bash", "-euo", "pipefail", "-c", _SUBMIT_FUNCTIONS["sge"] + "\njid=$(submit /x/s.sh); echo got"],
                         capture_output=True, text=True)
    assert out.returncode != 0 and "got" not in out.stdout and "Unable to run job: denied" in out.stderr
