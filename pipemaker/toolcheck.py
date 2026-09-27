"""Pre-flight check that every tool a build calls can be found.

A build calls executables in three places: the sample steps, the
experiment-level steps and the reference preparation job.  This module

* collects those executables from the generated commands (plus tools they
  depend on, e.g. cutadapt for trim_galore),
* writes ``check_<exp>.sh``, the first job of every build: it runs in the job
  environment on the cluster and stops the pipeline when a tool is missing,
* runs the same check here and now (``pipemaker check``, and before
  launching from the interactive mode), with the job's ``free_form`` set-up
  lines applied, and
* can try to install missing tools into the active environment (conda/mamba
  from bioconda + conda-forge, or pip for Python tools).
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from typing import Callable, Optional

from . import jobscript
from .config import Workflow

#: Shell keywords/builtins and utilities every system has; never reported.
IGNORED = {
    "if", "then", "else", "elif", "fi", "for", "foreach", "do", "done", "end", "endif", "while", "case", "esac",
    "[", "[[", "]", "test", "echo", "printf", "set", "unset", "setenv", "export", "exit", "cd", "true", "false",
    "source", ".", "eval", "exec", "time", "nohup", "env", "module", "conda", "mamba", "trap", "wait", "sleep",
    "mkdir", "rm", "ln", "cat", "cp", "mv", "cut", "date", "ls", "touch", "grep", "sed", "sh", "bash",
    "!", "(", ")", "{", "}",
}
#: Tools that other tools call behind the scenes.
HIDDEN_DEPENDENCIES = {"trim_galore": ["cutadapt"]}
#: Options whose value is itself a command (``STAR --readFilesCommand gunzip -c``).
COMMAND_OPTIONS = ("--readFilesCommand",)

#: executable -> conda package (bioconda / conda-forge)
CONDA_PACKAGES = {
    "trim_galore": "trim-galore", "fastqc": "fastqc", "bwa": "bwa", "samtools": "samtools", "STAR": "star",
    "htseq-count": "htseq", "multiqc": "multiqc", "cutadapt": "cutadapt", "wget": "wget", "md5sum": "coreutils",
    "gzip": "gzip", "gunzip": "gzip", "zcat": "gzip", "pigz": "pigz", "bowtie2": "bowtie2", "hisat2": "hisat2",
    "salmon": "salmon", "kallisto": "kallisto", "bedtools": "bedtools", "featureCounts": "subread",
    "picard": "picard", "gatk": "gatk4", "bcftools": "bcftools", "minimap2": "minimap2",
    "fasterq-dump": "sra-tools", "prefetch": "sra-tools", "curl": "curl",
}
#: executable -> pip package, used when no conda environment is active
PIP_PACKAGES = {"multiqc": "multiqc", "htseq-count": "HTSeq", "cutadapt": "cutadapt"}


@dataclass
class CheckResult:
    tools: dict                      # tool -> sorted list of steps using it
    missing: list = field(default_factory=list)
    checked: bool = True             # False when the check could not run here
    reason: str = ""

    @property
    def ok(self) -> bool:
        return self.checked and not self.missing

    def describe(self) -> str:
        if not self.checked:
            return f"The tool check could not run here: {self.reason}"
        if not self.missing:
            return f"All {len(self.tools)} tools used by the pipeline were found."
        lines = [f"{len(self.missing)} tool(s) used by the pipeline are not available in this environment:"]
        for tool in self.missing:
            steps = ", ".join(self.tools.get(tool, []))
            lines.append(f"  {tool:<16} used by: {steps}")
        return "\n".join(lines)


# ------------------------------------------------------------------ collect
def executables(line: str) -> list:
    """Executables a generated command line runs (first word of each pipeline segment)."""
    wrapped = re.match(r"^\s*(?:\S*/)?(?:ba)?sh -c '(.*)'\s*$", line)
    if wrapped:  # e.g. the reference downloads: sh -c 'set -e; wget ...; md5sum -c -'
        return executables(wrapped.group(1))
    found = []
    for segment in re.split(r"\|\||&&|\||;", line):
        try:
            words = shlex.split(segment, posix=True)
        except ValueError:
            words = segment.split()
        words = [w for w in words if not re.match(r"^\d*[<>]", w)]
        # skip leading keywords (if / then / !) to reach the command itself
        while words and words[0] in ("if", "then", "else", "elif", "!", "do", "(", "[", "[["):
            if words[0] in ("[", "[["):
                words = []
                break
            words = words[1:]
        if words:
            found.append(words[0])
        for i, word in enumerate(words[:-1]):
            if word in COMMAND_OPTIONS:
                found.append(words[i + 1])
    return found


def required_tools(step_lines: dict) -> dict:
    """``{step: [command lines]}`` -> ``{tool: [steps]}`` (keywords and core utilities left out)."""
    tools: dict = {}
    for step, lines in step_lines.items():
        for line in lines:
            for exe in executables(line):
                name = os.path.basename(exe) if not exe.startswith("/") else exe
                if name in IGNORED or os.path.basename(name) in IGNORED or name.startswith("$"):
                    continue
                tools.setdefault(name, set()).add(step)
                for dep in HIDDEN_DEPENDENCIES.get(os.path.basename(name), []):
                    tools.setdefault(dep, set()).add(f"{step} (via {os.path.basename(name)})")
    return {tool: sorted(steps) for tool, steps in sorted(tools.items(), key=lambda kv: kv[0].lower())}


def merge(*tool_maps: dict) -> dict:
    merged: dict = {}
    for tools in tool_maps:
        for tool, steps in tools.items():
            merged.setdefault(tool, set()).update(steps)
    return {tool: sorted(steps) for tool, steps in sorted(merged.items(), key=lambda kv: kv[0].lower())}


# ---------------------------------------------------------------- job script
def check_job_body(tools: dict, workflow: Workflow) -> str:
    """Body of ``check_<exp>.sh``: report every missing tool, fail if any is missing."""
    fam = jobscript.family(jobscript.shell_of(workflow))
    if fam == "csh":
        out = jobscript._csh_step("check_tools") + ['set pm_missing = ""']
        for tool, steps in tools.items():
            out += [
                f"sh -c 'command -v {shlex.quote(tool)}' >& /dev/null",
                "if ( $status != 0 ) then",
                "    sh -c " + shlex.quote(f'echo "MISSING {tool} (used by: {", ".join(steps)})" >&2'),
                f'    set pm_missing = "$pm_missing {tool}"',
                "endif",
            ]
        out += [
            'if ( "$pm_missing" != "" ) then',
            '    sh -c "echo \'[pipemaker] not found:$pm_missing - load or install them (general setting free_form), then rerun\' >&2"',
            '    echo "FAILED check_tools $PM_START `date +%s` 1" > "$PM_STATUS"',
            "    exit 1",
            "endif",
            f'echo "[pipemaker] all {len(tools)} tools found"',
            'echo "DONE check_tools $PM_START `date +%s` 0" > "$PM_STATUS"',
        ]
        return "\n".join(out) + "\n"
    out = ["pm_step check_tools" if fam == "posix" else "", 'pm_missing=""']
    for tool, steps in tools.items():
        out.append(f'if ! command -v {shlex.quote(tool)} >/dev/null 2>&1; then '
                   f'echo "MISSING {tool} (used by: {", ".join(steps)})" >&2; pm_missing="$pm_missing {tool}"; fi')
    out += [
        'if [ -n "$pm_missing" ]; then',
        '    echo "[pipemaker] not found:$pm_missing - load or install them (general setting free_form), then rerun" >&2',
        "    exit 1",
        "fi",
        f'echo "[pipemaker] all {len(tools)} tools found"',
    ]
    return "\n".join(line for line in out if line) + "\n"


# ------------------------------------------------------------- local check
def check_here(tools: dict, workflow: Workflow, runner: Callable = subprocess.run, timeout: int = 120) -> CheckResult:
    """Look the tools up in the current environment, after the jobs' ``free_form`` set-up."""
    if not tools:
        return CheckResult(tools)
    shell = jobscript.shell_of(workflow)
    fam = jobscript.family(shell)
    shell_path = shutil.which(shell)
    if fam is None or shell_path is None:
        return CheckResult(tools, checked=False,
                           reason=f"the job shell '{shell}' is not available on this machine" if shell_path is None
                           else f"unsupported job shell '{shell}'")
    free_form = str(workflow.general("free_form", "") or "")
    names = " ".join(shlex.quote(t) for t in tools)
    if fam == "csh":
        script = f"{free_form}\nforeach pm_t ({names})\n  sh -c \"command -v $pm_t\" >& /dev/null\n" \
                 "  if ( $status != 0 ) echo \"MISSING $pm_t\"\nend\n"
    else:
        script = f"{free_form}\nfor pm_t in {names}; do command -v \"$pm_t\" >/dev/null 2>&1 || echo \"MISSING $pm_t\"; done\n"
    with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as handle:
        handle.write(script)
        path = handle.name
    try:
        result = runner([shell_path, path], capture_output=True, text=True, timeout=timeout, check=False)
    except (OSError, subprocess.SubprocessError) as exc:
        return CheckResult(tools, checked=False, reason=str(exc))
    finally:
        os.unlink(path)
    missing = [line.split(None, 1)[1].strip() for line in result.stdout.splitlines() if line.startswith("MISSING ")]
    return CheckResult(tools, missing=missing)


# ----------------------------------------------------------------- install
@dataclass
class InstallPlan:
    commands: list                     # argv lists to run, in order
    manual: dict                       # tool -> why it cannot be installed automatically
    environment: str = ""

    def describe(self) -> str:
        lines = []
        if self.commands:
            lines.append(f"Will install into {self.environment}:")
            lines += ["  $ " + " ".join(shlex.quote(a) for a in cmd) for cmd in self.commands]
        for tool, why in self.manual.items():
            lines.append(f"  {tool}: {why}")
        return "\n".join(lines)


def install_plan(missing: list, env: Optional[dict] = None, which: Callable = shutil.which) -> InstallPlan:
    env = os.environ if env is None else env
    prefix = env.get("CONDA_PREFIX", "")
    manager = next((m for m in ("mamba", "micromamba", "conda") if which(m)), None) if prefix else None
    manual, conda_pkgs, pip_pkgs = {}, [], []
    for tool in missing:
        if tool.startswith("/"):
            manual[tool] = "an absolute path that does not exist - fix the step's 'path' setting"
        elif manager:
            package = CONDA_PACKAGES.get(tool, tool.lower())
            if package not in conda_pkgs:
                conda_pkgs.append(package)
        elif tool in PIP_PACKAGES:
            if PIP_PACKAGES[tool] not in pip_pkgs:
                pip_pkgs.append(PIP_PACKAGES[tool])
        else:
            manual[tool] = ("no active conda environment to install it into - load it with 'module load ...' "
                            "or install it and add it to your PATH")
    commands = []
    environment = ""
    if conda_pkgs:
        commands.append([manager, "install", "-y", "-p", prefix, "-c", "conda-forge", "-c", "bioconda", *conda_pkgs])
        environment = f"the conda environment {env.get('CONDA_DEFAULT_ENV') or prefix} ({manager})"
    if pip_pkgs:
        in_venv = sys.prefix != getattr(sys, "base_prefix", sys.prefix)
        commands.append([sys.executable, "-m", "pip", "install", *([] if in_venv else ["--user"]), *pip_pkgs])
        environment = environment or f"Python {sys.executable} (pip)"
    return InstallPlan(commands, manual, environment)


def install(plan: InstallPlan, runner: Callable = subprocess.run) -> list:
    """Run the install commands; returns the ones that failed.

    A failing conda command is retried one package at a time so that one
    unknown package does not block the others.
    """
    failed = []
    for cmd in plan.commands:
        if runner(cmd, check=False).returncode == 0:
            continue
        if cmd and os.path.basename(str(cmd[0])) in ("mamba", "micromamba", "conda"):
            head, packages = cmd[:9], cmd[9:]
            for package in packages:
                if runner(head + [package], check=False).returncode != 0:
                    failed.append(head + [package])
        else:
            failed.append(cmd)
    return failed


ENV_NOTE = ("Jobs inherit the environment they are submitted from (SGE -V, SLURM and LSF by default), so submit "
            "from this environment - or add its activation (e.g. 'conda activate NAME') to the general setting "
            "'free_form'.")


def run_preflight(tools: dict, workflow: Workflow, auto_fix: bool = False, stream=None,
                  runner: Callable = subprocess.run, check: Callable = None) -> CheckResult:
    """Check the tools here; with ``auto_fix`` try to install the missing ones and check again."""
    stream = stream or sys.stdout
    check = check or (lambda t, w: check_here(t, w, runner=runner))
    result = check(tools, workflow)
    if result.ok or not result.checked:
        stream.write(result.describe() + "\n")
        return result
    stream.write(result.describe() + "\n")
    if not auto_fix:
        return result
    plan = install_plan(result.missing)
    stream.write(plan.describe() + "\n")
    if plan.commands:
        stream.flush()
        install(plan, runner=runner)
        result = check(tools, workflow)
        stream.write(("After installing: " + result.describe()) + "\n")
        if result.ok:
            stream.write(ENV_NOTE + "\n")
    return result


def experiment_tools(workspace, experiment: str) -> tuple:
    """``(tools, workflow)`` recorded in the manifest of a build."""
    import json

    from .config import Workflow as _Workflow

    manifest = json.loads((workspace.experiment_dir(experiment) / "pipemaker_manifest.json").read_text(encoding="utf-8"))
    return manifest.get("tools", {}), _Workflow(**manifest["workflow"])
