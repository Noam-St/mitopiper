import io
import os
import shutil
import subprocess

import pytest

from pipemaker import monitor as mon
from pipemaker import toolcheck
from pipemaker.cli import main
from pipemaker.config import Workflow, Workspace
from pipemaker.diagnose import Context, diagnose
from pipemaker.local_samples import make_tsv
from pipemaker.references import ReferenceRegistry, ReferenceSet
from pipemaker.rerun import RerunOptions, rerun

SHELLS = [s for s in ("sh", "tcsh") if shutil.which(s)]


def fake_tool(folder, name):
    folder.mkdir(exist_ok=True)
    path = folder / name
    path.write_text("#!/bin/sh\nexit 0\n")
    path.chmod(0o755)
    return path


def test_tools_are_collected_from_the_generated_commands():
    tools = toolcheck.required_tools({
        "wget": ["wget -c -O /o/x.fastq.gz ftp/x.fastq.gz", 'echo "abc  /o/x" | md5sum -c -', "sleep 30",
                 'cat "/a" "/b" > "/c"', "rm -f /a", 'ln -sf "/in" "/o"'],
        "trim": ["trim_galore -q 20 --paired /a /b"],
        "STAR-alignment": ["/opt/STAR/bin/STAR --readFilesCommand gunzip -c --genomeDir /g"],
        "prepare_references": ['if [ ! -s "/r.fa.fai" ]; then samtools faidx "/r.fa"; fi'],
    })
    assert tools == {
        "/opt/STAR/bin/STAR": ["STAR-alignment"], "cutadapt": ["trim (via trim_galore)"],
        "gunzip": ["STAR-alignment"], "md5sum": ["wget"], "samtools": ["prepare_references"],
        "trim_galore": ["trim"], "wget": ["wget"],
    }


@pytest.mark.parametrize("shell", SHELLS)
def test_check_here_applies_the_free_form_setup(tmp_path, shell):
    fake_tool(tmp_path / "bin", "mytool")
    wf = Workflow.default()
    wf.set_general("environment", shell)
    tools = {"mytool": ["count"], "ls": ["x"]}
    assert toolcheck.check_here(tools, wf).missing == ["mytool"]
    setup = f"setenv PATH {tmp_path}/bin:$PATH" if shell == "tcsh" else f'export PATH="{tmp_path}/bin:$PATH"'
    wf.set_general("free_form", setup)
    assert toolcheck.check_here(tools, wf).ok


def test_check_here_reports_when_it_cannot_run():
    wf = Workflow.default()
    wf.set_general("environment", "fish-that-does-not-exist")
    result = toolcheck.check_here({"x": ["s"]}, wf)
    assert not result.checked and "not available" in result.reason


def test_install_plan_prefers_the_active_conda_environment():
    env = {"CONDA_PREFIX": "/envs/rna", "CONDA_DEFAULT_ENV": "rna"}
    which = lambda name: "/usr/bin/mamba" if name == "mamba" else None  # noqa: E731
    plan = toolcheck.install_plan(["trim_galore", "STAR", "cutadapt", "weirdtool", "/opt/x/tool"], env, which)
    assert plan.commands == [["mamba", "install", "-y", "-p", "/envs/rna", "-c", "conda-forge", "-c", "bioconda",
                              "trim-galore", "star", "cutadapt", "weirdtool"]]
    assert "/opt/x/tool" in plan.manual and "rna" in plan.environment
    plan = toolcheck.install_plan(["multiqc", "bwa"], env={}, which=lambda n: None)
    assert plan.commands[0][1:4] == ["-m", "pip", "install"] and plan.commands[0][-1] == "multiqc"
    assert "bwa" in plan.manual and "module load" in plan.manual["bwa"]


def test_install_retries_package_by_package():
    calls = []

    def runner(cmd, check=False):
        calls.append(cmd[-1])
        return subprocess.CompletedProcess(cmd, 1 if cmd[-1] in ("b", "bad") else 0)

    plan = toolcheck.InstallPlan([["conda", "install", "-y", "-p", "/e", "-c", "conda-forge", "-c", "bioconda", "a", "bad"]], {})
    failed = toolcheck.install(plan, runner=runner)
    assert calls == ["bad", "a", "bad"] and [f[-1] for f in failed] == ["bad"]


# --------------------------------------------------------------- whole pipeline
def workspace(tmp_path, shell="sh"):
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("default")
    wf.set_general("environment", shell)
    wf.pipeline["count"] = {"key": "mytool count", "description": "count", "switch": True, "path": "", "outdir": "counts",
                            "inputdir": "fastq", "input_ext": ".fastq.gz", "output_ext": ".txt", "scope": "sample",
                            "structure": "key input output_pipe"}
    wf.pipeline["multiqc"].update(switch=True, structure="echo summary")
    wf.run_order = ["wget", "count", "multiqc"]
    ws.write_workflow(wf, "t")
    fq = tmp_path / "fq"
    fq.mkdir()
    (fq / "a.fastq").write_text("@r\nACGT\n+\nIIII\n")
    make_tsv(fq, "Homo sapiens", tmp_path / "exp.tsv")
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath="/x.fa"))
    return ws


@pytest.mark.parametrize("shell", SHELLS)
def test_check_job_stops_the_pipeline_and_a_rerun_with_the_fix_completes_it(tmp_path, shell):
    ws = workspace(tmp_path, shell)
    assert main(["-w", str(ws.root), "build", "--no-check", str(tmp_path / "exp.tsv"), "-s", "local", "-p", "t"]) == 0
    assert (ws.output_dir / "exp" / "scripts" / "check_exp.sh").exists()
    subprocess.run(["bash", str(ws.root / "master.sh")], check=True, capture_output=True)
    jobs = {j.name: j for j in mon.collect(ws, ["exp"], mon.ErrorScanner(), lambda s, i: None)}
    assert jobs["check_exp"].state == mon.FAILED and jobs["check_exp"].kind == "check"
    assert jobs["s_a_0"].state == mon.SKIPPED and jobs["e_exp"].state == mon.SKIPPED
    assert not (ws.output_dir / "exp" / "fastq" / "a_0_HS.fastq.gz").exists()        # nothing ran
    (sug, *_) = diagnose(jobs["check_exp"], Context.from_manifest(None))
    assert "missing tool(s): mytool" in sug.cause and sug.rerun["free_form_hint"] == "module load mytool"
    # fix: make mytool available to the jobs, rerun everything that did not run
    fake_tool(tmp_path / "bin", "mytool")
    setup = f"setenv PATH {tmp_path}/bin:$PATH" if shell == "tcsh" else f'export PATH="{tmp_path}/bin:$PATH"'
    assert rerun(ws, ["exp"], RerunOptions(free_form=[setup]), submit=True, stream=io.StringIO()) == 0
    states = {j.name: j.state for j in mon.collect(ws, ["exp"], mon.ErrorScanner(), lambda s, i: None)}
    assert states == {"check_exp": mon.DONE, "s_a_0": mon.DONE, "e_exp": mon.DONE}


def test_build_halts_on_missing_tools_unless_fixed(tmp_path, monkeypatch, capsys):
    ws = workspace(tmp_path)
    args = ["-w", str(ws.root), "build", str(tmp_path / "exp.tsv"), "-s", "slurm", "-p", "t"]
    assert main(args) == 1
    assert "mytool" in capsys.readouterr().out
    # --auto-fix: the (fake) installer puts mytool on the PATH
    bin_dir = tmp_path / "bin"
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    monkeypatch.setattr(toolcheck, "install_plan", lambda missing: toolcheck.InstallPlan(
        [["sh", "-c", f"mkdir -p {bin_dir} && printf '#!/bin/sh\\n' > {bin_dir}/mytool && chmod +x {bin_dir}/mytool"]],
        {}, "a test environment"))
    assert main(args + ["--auto-fix"]) == 0
    assert "All 2 tools" in capsys.readouterr().out     # mytool (+ md5sum is not needed: local files) and ...
    assert main(["-w", str(ws.root), "check"]) == 0


def test_tool_check_can_be_switched_off(tmp_path):
    ws = workspace(tmp_path)
    wf = ws.load_workflow("t")
    wf.set_general("tool_check", False)
    ws.write_workflow(wf, "t")
    assert main(["-w", str(ws.root), "build", "--no-check", str(tmp_path / "exp.tsv"), "-p", "t"]) == 0
    assert not (ws.output_dir / "exp" / "scripts" / "check_exp.sh").exists()


def test_interactive_preflight_adds_a_setup_line_and_rebuilds(tmp_path):
    from pipemaker.interactive import run_interactive

    ws = workspace(tmp_path)
    fake_tool(tmp_path / "bin", "mytool")
    answers = iter(["0", str(tmp_path / "exp.tsv"), "1", "n", "all", "",     # scheduler sge, metadata, preset t
                    "2", f'export PATH="{tmp_path}/bin:$PATH"',              # add a set-up line -> rebuild -> recheck ok
                    "n"])                                                    # don't launch
    assert run_interactive(ws, lambda prompt="": next(answers)) == 0
    assert f'export PATH="{tmp_path}/bin:$PATH"' in (ws.output_dir / "exp" / "scripts" / "s_a_0.sh").read_text()
    assert f"{tmp_path}/bin" in ws.load_workflow("current").general("free_form")
