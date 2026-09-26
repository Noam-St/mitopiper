import pytest

from pipemaker.schedulers import ExperimentJobs, JobResources, header, master_script, normalise


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
    lsf = header("lsf", "job", "/exp", res)
    assert "#BSUB -P lab" in lsf and '#BSUB -R "rusage[mem=8G]"' in lsf
    assert header("local", "job", "/exp", res) == "#!/bin/bash\n"


def test_no_account_or_queue_when_empty():
    slurm = header("slurm", "job", "/exp", JobResources())
    assert "--account" not in slurm and "--partition" not in slurm


@pytest.mark.parametrize("scheduler,needle", [("sge", "qsub -terse"), ("slurm", "sbatch --parsable"), ("lsf", "bsub")])
def test_master_script_dependencies(scheduler, needle):
    text = master_script(scheduler, [ExperimentJobs("E", "/e/scripts/prep.sh", ["/e/scripts/a.sh"], "/e/scripts/exp.sh")])
    assert needle in text
    assert "jid=$(submit /e/scripts/prep.sh)" in text and 'prep+=("$jid")' in text
    assert 'jid=$(submit /e/scripts/exp.sh ${jobs[@]+"${jobs[@]}"})' in text
    # every submission is recorded for the monitor
    assert "jobs_file=/e/jobs.tsv" in text and "marker=/e/.submitting" in text
    assert 'record "$jid" "$script" sample' in text and "record \"$jid\" /e/scripts/exp.sh experiment" in text
