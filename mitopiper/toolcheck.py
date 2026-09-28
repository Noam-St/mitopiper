"""Pre-flight check that every tool a build calls can be found.

A build calls executables in three places: the sample steps, the
experiment-level steps and the reference preparation job.  This module

* collects those executables from the generated commands (plus tools they
  depend on, e.g. cutadapt for trim_galore),
* writes ``check_<exp>.sh``, the first job of every build: it runs in the job
  environment on the cluster and stops the pipeline when a tool is missing,
* runs the same check here and now (``mitopiper check``, and before
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
from functools import partial
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
#: Tools that are not packaged anywhere: how to get them.
MANUAL_INSTALL = {
    "numt_parser.py": "not packaged - 'git clone https://github.com/adeflamingh/NuMt_parser' and add its bin/ "
                      "directory to your PATH (or set the step's 'path' to bin/numt_parser.py)",
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

    @property
    def failed(self) -> bool:
        """The check ran and found missing tools."""
        return self.checked and bool(self.missing)

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
        # skip leading keywords (if / then / !) to reach the command itself; a test ([ ... ]) runs nothing
        while words and words[0] in ("if", "then", "else", "elif", "!", "do", "("):
            words = words[1:]
        if not words or words[0] in ("[", "[["):
            continue
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
                base = os.path.basename(exe)
                if base in IGNORED or exe.startswith("$"):
                    continue
                name = exe if exe.startswith("/") else base
                tools.setdefault(name, set()).add(step)
                for dep in HIDDEN_DEPENDENCIES.get(base, []):
                    tools.setdefault(dep, set()).add(f"{step} (via {base})")
    return {tool: sorted(steps) for tool, steps in sorted(tools.items(), key=lambda kv: kv[0].lower())}


def merge(*tool_maps: dict) -> dict:
    merged: dict = {}
    for tools in tool_maps:
        for tool, steps in tools.items():
            merged.setdefault(tool, set()).update(steps)
    return {tool: sorted(steps) for tool, steps in sorted(merged.items(), key=lambda kv: kv[0].lower())}


# ---------------------------------------------------------------- job script
#: Prints the first line with a version number that ``$1 --version`` (or, for tools such as bwa
#: that only show it in their usage text, ``$1`` alone) writes. No input, time-limited.
_VERSION_OF = ('v() { t=""; command -v timeout >/dev/null 2>&1 && t="timeout 20"; '
               '{ $t "$1" --version </dev/null 2>&1 || $t "$1" </dev/null 2>&1; } '
               "| grep -m1 -E '[0-9]+\\.[0-9]+' | tr '\\t' ' '; }")


def lookup_script(tools: dict, versions_file: Optional[str] = None) -> str:
    """A POSIX ``sh`` script that prints ``MISSING <tool> (used by: ...)`` on stderr for every
    tool it cannot find and exits with status 1 if any is missing.

    With ``versions_file``, the path and reported version of every tool found are written
    to it as ``tool<TAB>path<TAB>version`` lines (for methods sections and provenance).
    """
    lines = ['m=""']
    if versions_file:
        quoted_file = shlex.quote(versions_file)
        lines += [_VERSION_OF, f"printf 'tool\\tpath\\tversion\\n' > {quoted_file}"]
    for tool, steps in tools.items():
        quoted = shlex.quote(tool)
        found = ""
        if versions_file:
            found = (f' && printf "%s\\t%s\\t%s\\n" {quoted} "$(command -v {quoted})" "$(v {quoted})" '
                     f">> {quoted_file}")
        lines.append(f'{{ command -v {quoted} >/dev/null 2>&1{found}; }} || '
                     f'{{ echo "MISSING {tool} (used by: {", ".join(steps)})" >&2; m="$m {tool}"; }}')
    lines.append('if [ -n "$m" ]; then echo "[mitopiper] not found:$m - load or install them (general setting free_form), then rerun" >&2; exit 1; fi')
    lines.append(f'echo "[mitopiper] all {len(tools)} tools found"')
    return "; ".join(lines)


def lookup_command(tools: dict, versions_file: Optional[str] = None) -> str:
    """:func:`lookup_script` as one command line that works in every job shell."""
    return "sh -c " + shlex.quote(lookup_script(tools, versions_file))


VERSIONS_FILE = "tool_versions.tsv"


def check_job_body(tools: dict, workflow: Workflow, exp_dir: Optional[str] = None) -> str:
    """Body of ``check_<exp>.sh``: report every missing tool, fail if any is missing, and record
    the versions of the tools found in ``<exp>/tool_versions.tsv``."""
    from .builder import StepCommands

    versions = f"{exp_dir}/{VERSIONS_FILE}" if exp_dir else None
    step = StepCommands("check_tools", "check that the tools of the pipeline can be found", [lookup_command(tools, versions)])
    return jobscript.body([step], "this build", workflow)


# ------------------------------------------------------------- local check
def check_here(tools: dict, workflow: Workflow, runner: Callable = subprocess.run, timeout: int = 120) -> CheckResult:
    """Look the tools up in the current environment, after the jobs' ``free_form`` set-up."""
    if not tools:
        return CheckResult(tools)
    shell = jobscript.shell_of(workflow)
    shell_path = shutil.which(shell)
    if jobscript.family(shell) is None or shell_path is None:
        return CheckResult(tools, checked=False,
                           reason=f"the job shell '{shell}' is not available on this machine" if shell_path is None
                           else f"unsupported job shell '{shell}'")
    free_form = str(workflow.general("free_form", "") or "")
    with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as handle:
        handle.write(f"{free_form}\n{lookup_command(tools)}\n")
        path = handle.name
    try:
        result = runner([shell_path, path], capture_output=True, text=True, timeout=timeout, check=False)
    except (OSError, subprocess.SubprocessError) as exc:
        return CheckResult(tools, checked=False, reason=str(exc))
    finally:
        os.unlink(path)
    missing = [line.split()[1] for line in result.stderr.splitlines() if line.startswith("MISSING ")]
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
        elif tool in MANUAL_INSTALL:
            manual[tool] = MANUAL_INSTALL[tool]
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
    check = check or partial(check_here, runner=runner)
    result = check(tools, workflow)
    stream.write(result.describe() + "\n")
    if not (result.failed and auto_fix):
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
    from .build import load_manifest  # build imports this module

    manifest = load_manifest(workspace.experiment_dir(experiment))
    return manifest.get("tools", {}), Workflow(**manifest["workflow"])
