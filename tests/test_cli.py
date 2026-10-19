import json
import subprocess

from mitopiper.cli import main
from mitopiper.config import Workspace
from mitopiper.references import ReferenceSet


def run(ws, *args):
    return main(["-w", str(ws.root), *args])


def test_build_writes_scripts_master_and_manifest(registered, data_dir):
    ws = registered
    assert run(ws, "build", "--no-check", str(data_dir / "ena_runs.tsv"), "-s", "slurm") == 0
    exp = ws.output_dir / "ena_runs"
    scripts = sorted(p.name for p in (exp / "scripts").iterdir())
    assert scripts == ["check_ena_runs.sh", "s_SAMEA1_0.sh", "s_SAMEA3_1.sh"]
    text = (exp / "scripts" / "s_SAMEA1_0.sh").read_text()
    assert text.startswith("#!/bin/bash\n#SBATCH")
    assert "set -eo pipefail" in text
    assert (exp / "metadata.tsv").exists()
    assert (data_dir / "ena_runs.tsv").exists()  # the user's file is never moved
    manifest = json.loads((exp / "mitopiper_manifest.json").read_text())
    assert manifest["scheduler"] == "slurm" and len(manifest["samples"]) == 2
    master = (ws.root / "master.sh").read_text()
    assert "sbatch --parsable" in master
    subprocess.run(["bash", "-n", str(ws.root / "master.sh")], check=True)
    subprocess.run(["bash", "-n", str(exp / "scripts" / "s_SAMEA1_0.sh")], check=True)
    assert ws.load_workflow("last_used")


def test_dry_run_writes_plan_only(registered, data_dir):
    ws = registered
    assert run(ws, "build", str(data_dir / "ena_runs.tsv"), "--dry-run", "--steps", "wget,trim", "--samples", "0") == 0
    exp = ws.output_dir / "ena_runs"
    plan = json.loads((exp / "plan.json").read_text())
    assert plan["steps"] == ["wget", "trim"]
    assert [s["name"] for s in plan["samples"]] == ["SAMEA1_0_DR"]
    assert not (exp / "scripts").exists()
    assert not (ws.root / "master.sh").exists()


def test_missing_reference_is_reported(workspace, data_dir):
    assert run(workspace, "build", "--no-check", str(data_dir / "ena_runs.tsv")) == 1


def test_reference_override(registered, data_dir):
    ws = registered
    from mitopiper.references import ReferenceRegistry

    ReferenceRegistry(ws).add("Danio_rerio", ReferenceSet(refpath="/refs/dr/GRCz11.fa"), "GRCz11", activate=False)
    assert run(ws, "build", "--no-check", str(data_dir / "ena_runs.tsv"), "--reference", "GRCz11", "--samples", "0") == 0
    text = (ws.output_dir / "ena_runs" / "scripts" / "s_SAMEA1_0.sh").read_text()
    assert "/refs/dr/GRCz11.fa" in text


def test_local_scheduler_runs_jobs(tmp_path, data_dir):
    ws = Workspace(tmp_path / "ws")
    fq = tmp_path / "fastq"
    fq.mkdir()
    (fq / "s1_R1.fastq").write_text("@r\nACGT\n+\nIIII\n")
    (fq / "s1_R2.fastq").write_text("@r\nACGT\n+\nIIII\n")
    assert run(ws, "make-tsv", str(fq), "--org", "Homo sapiens", "-o", str(tmp_path / "local.tsv")) == 0
    assert run(ws, "refs", "add", "Homo sapiens", "--fasta", "/x.fa") == 0
    assert run(ws, "build", "--no-check", str(tmp_path / "local.tsv"), "-s", "local", "--steps", "wget") == 0
    subprocess.run(["bash", str(ws.root / "master.sh")], check=True, capture_output=True)
    out = ws.output_dir / "local" / "fastq"
    assert sorted(p.name for p in out.iterdir()) == ["s1_0_HS_1.fastq.gz", "s1_0_HS_2.fastq.gz"]


def test_presets_and_cleanup(registered, data_dir, capsys):
    ws = registered
    assert run(ws, "presets", "list") == 0
    assert "default" in capsys.readouterr().out
    assert run(ws, "build", "--no-check", str(data_dir / "ena_runs.tsv")) == 0
    exp = ws.output_dir / "ena_runs"
    (exp / "stdout" / "job.out").write_text("log")
    assert run(ws, "cleanup", "--yes") == 0
    assert not (exp / "stdout" / "job.out").exists()
    assert not (exp / "scripts").exists()
    assert not (ws.root / "master.sh").exists()
    assert (exp / "metadata.tsv").exists()


def test_presets_set(tmp_path, capsys):
    ws = Workspace(tmp_path / "ws")
    assert run(ws, "init") == 0
    assert run(ws, "presets", "set", "mt-variants,mt-variants-numt", "qsub_queue=bioinfo.q", "qsub_threads=2",
               "array_jobs=true", "index.path=/opt/bin/samtools", "index.memory=4G") == 0
    assert "mt-variants: qsub_queue = 'bioinfo.q' (was 'bioinfo.q')" in capsys.readouterr().out
    for preset in ("mt-variants", "mt-variants-numt"):
        wf = ws.load_workflow(preset)
        assert (wf.general("qsub_threads"), wf.general("array_jobs")) == (2, True)
        assert wf.pipeline["index"]["path"] == "/opt/bin/samtools" and wf.pipeline["index"]["memory"] == "4G"
    # a typo, or a step one preset lacks, changes nothing in any preset
    assert run(ws, "presets", "set", "mt-variants", "qsub_mem=8G") == 1
    assert run(ws, "presets", "set", "mt-variants,mt-variants-numt", "walltime=1:00:00", "nosuchstep.threads=2") == 1
    assert ws.load_workflow("mt-variants").general("walltime") != "1:00:00"


def test_end_of_input_exits_cleanly(tmp_path, monkeypatch, capsys):
    import io
    import sys

    monkeypatch.setattr(sys, "stdin", io.StringIO(""))  # Ctrl-D / piped input ran out
    assert main(["-w", str(tmp_path / "ws")]) == 130
    assert "Aborted." in capsys.readouterr().out


def test_builds_made_before_the_rename_are_still_read(registered, data_dir):
    from mitopiper import monitor
    from mitopiper.build import load_manifest

    ws = registered
    assert run(ws, "build", "--no-check", str(data_dir / "ena_runs.tsv"), "-s", "slurm") == 0
    exp = ws.output_dir / "ena_runs"
    (exp / "mitopiper_manifest.json").rename(exp / "pipemaker_manifest.json")
    assert load_manifest(exp)["scheduler"] == "slurm"
    assert monitor.experiment_scheduler(exp) == "slurm"


def test_step_resources_reach_the_job_header(registered, data_dir):
    ws = registered
    wf = ws.load_workflow("default")
    wf.pipeline["STAR-alignment"]["switch"] = True
    wf.run_order = ["wget", "STAR-alignment"]
    ws.write_workflow(wf, "star")
    assert run(ws, "build", "--no-check", str(data_dir / "ena_runs.tsv"), "-s", "slurm", "-p", "star") == 0
    text = (ws.output_dir / "ena_runs" / "scripts" / "s_SAMEA1_0.sh").read_text()
    assert "#SBATCH --mem=32G" in text and "#SBATCH --cpus-per-task=8" in text and "--runThreadN 8" in text


def test_presets_reset_restores_the_bundled_version(tmp_path, capsys):
    ws = Workspace(tmp_path / "ws")
    assert run(ws, "init") == 0
    assert run(ws, "presets", "set", "mt-variants-numtparser", "qsub_queue=old.q") == 0
    assert run(ws, "presets", "reset", "mt-variants-numtparser") == 0
    assert ws.load_workflow("mt-variants-numtparser").general("qsub_queue") != "old.q"
    kept = [p for p in ws.list_presets() if p.startswith("mt-variants-numtparser_old_")]
    assert kept and ws.load_workflow(kept[0]).general("qsub_queue") == "old.q"
    assert run(ws, "presets", "reset", "no-such-preset") == 1


def test_version_names_the_commit(monkeypatch, capsys):
    import pytest

    import mitopiper
    from mitopiper import _commit

    monkeypatch.setattr(_commit, "COMMIT", "abc1234")              # as git archive stamps the release
    monkeypatch.setattr(_commit, "DATE", "2026-10-19")
    with pytest.raises(SystemExit):
        main(["--version"])
    assert capsys.readouterr().out.strip() == f"mitopiper {mitopiper.__version__} (commit abc1234, 2026-10-19)"
    monkeypatch.setattr(_commit, "COMMIT", "$Format:%h$")          # not stamped and no git: says so
    monkeypatch.setattr("subprocess.run", lambda *a, **k: (_ for _ in ()).throw(OSError("no git")))
    assert mitopiper.build_info().startswith("commit unknown")
