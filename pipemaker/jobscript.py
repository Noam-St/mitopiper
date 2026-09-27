"""Shell-aware assembly of job scripts.

Every job script reports its progress to ``<exp>/status/<job>.status`` (one
line ``STATE STEP START UPDATED RC``) and ``<exp>/status/<job>.steps`` for
``pipemaker monitor``, checks that the jobs it depends on finished
successfully, and stops at the first failing command when ``fail_fast`` is on.

Two shell families are supported:

* ``posix`` - bash, zsh, ksh, sh, dash: an ``EXIT`` trap records success or
  failure with the exit code.
* ``csh`` - csh, tcsh: no traps, so the exit status is checked after every
  command.

Any other shell gets no tracking (the monitor then relies on the scheduler).
"""

from __future__ import annotations

from typing import Optional

from .config import Workflow

POSIX_SHELLS = ("bash", "zsh", "ksh", "sh", "dash")
CSH_SHELLS = ("csh", "tcsh")
PIPEFAIL_SHELLS = ("bash", "zsh", "ksh")


def shell_of(workflow: Workflow) -> str:
    return str(workflow.general("environment", "bash") or "bash").strip()


def family(shell: str) -> Optional[str]:
    if shell in POSIX_SHELLS:
        return "posix"
    if shell in CSH_SHELLS:
        return "csh"
    return None


def fail_fast(workflow: Workflow) -> bool:
    return workflow.flag("fail_fast", True)


def append_free_form(workflow: Workflow, lines: list) -> None:
    """Append set-up lines to the ``free_form`` general setting (one per line)."""
    lines = [line for line in lines if line]
    if not lines:
        return
    current = str(workflow.general("free_form", "") or "")
    if current and not current.endswith("\n"):
        current += "\n"
    workflow.set_general("free_form", current + "\n".join(lines) + "\n")


# ------------------------------------------------------------------ posix
_POSIX_TRACKING = r"""PM_STATUS="@STATUS@.status"; PM_STEPS="@STATUS@.steps"; PM_STEP="-"; PM_START=$(date +%s)
mkdir -p "@STATUSDIR@"; : > "$PM_STEPS"
pm_step() {
    PM_STEP="$1"
    echo "RUNNING $1 $PM_START $(date +%s) -" > "$PM_STATUS"
    echo "$1 $(date +%s)" >> "$PM_STEPS"
    echo "[pipemaker] step $1 started $(date '+%F %T')"
}
pm_exit() {
    pm_rc=$?
    if [ "$pm_rc" -eq 0 ]; then pm_state=DONE; else pm_state=FAILED; fi
    echo "$pm_state $PM_STEP $PM_START $(date +%s) $pm_rc" > "$PM_STATUS"
}
trap pm_exit EXIT
echo "RUNNING - $PM_START $PM_START -" > "$PM_STATUS"
"""

_POSIX_GUARD = r"""# Run only if the jobs this one depends on finished successfully
# (SGE's -hold_jid also releases a job when its dependency failed).
for pm_dep in @DEPS@; do
    pm_state=$(cut -d' ' -f1 "@STATUSDIR@/$pm_dep.status" 2>/dev/null || true)
    if [ "$pm_state" != DONE ]; then
        echo "[pipemaker] skipped: $pm_dep did not finish successfully" >&2
        trap - EXIT
        echo "SKIPPED $pm_dep $PM_START $(date +%s) -" > "$PM_STATUS"
        exit 0
    fi
done
"""

# -------------------------------------------------------------------- csh
_CSH_TRACKING = r"""set PM_STATUS = "@STATUS@.status"
set PM_STEPS = "@STATUS@.steps"
set PM_STEP = "-"
set PM_START = `date +%s`
mkdir -p "@STATUSDIR@"
cp /dev/null "$PM_STEPS"
echo "RUNNING - $PM_START $PM_START -" > "$PM_STATUS"
"""

_CSH_GUARD = r"""# Run only if the jobs this one depends on finished successfully
foreach pm_dep (@DEPS@)
    set pm_state = "missing"
    if ( -e "@STATUSDIR@/$pm_dep.status" ) set pm_state = `cut -d' ' -f1 "@STATUSDIR@/$pm_dep.status"`
    if ( "$pm_state" != "DONE" ) then
        echo "[pipemaker] skipped: $pm_dep did not finish successfully"
        echo "SKIPPED $pm_dep $PM_START `date +%s` -" > "$PM_STATUS"
        exit 0
    endif
end
"""


def _csh_step(step: str) -> list:
    return [
        f'set PM_STEP = "{step}"',
        f'echo "RUNNING {step} $PM_START `date +%s` -" > "$PM_STATUS"',
        f'echo "{step} `date +%s`" >> "$PM_STEPS"',
        f'echo "[pipemaker] step {step} started `date \'+%F %T\'`"',
    ]


_CSH_CHECK = [
    "set pm_rc = $status",
    "if ( $pm_rc != 0 ) then",
    '    echo "FAILED $PM_STEP $PM_START `date +%s` $pm_rc" > "$PM_STATUS"',
    "    exit $pm_rc",
    "endif",
]


# ------------------------------------------------------------------ public
def preamble(workflow: Workflow, exp_dir: str, job_name: str, depends_on: Optional[list] = None) -> str:
    """Tracking, dependency guard, ``free_form`` and fail-fast set-up."""
    shell = shell_of(workflow)
    fam = family(shell)
    status_dir = f"{exp_dir}/status"
    base = f"{status_dir}/{job_name}"
    text = ""
    if fam == "posix":
        text += _POSIX_TRACKING.replace("@STATUSDIR@", status_dir).replace("@STATUS@", base)
        if depends_on:
            text += _POSIX_GUARD.replace("@STATUSDIR@", status_dir).replace("@DEPS@", " ".join(depends_on))
    elif fam == "csh":
        text += _CSH_TRACKING.replace("@STATUSDIR@", status_dir).replace("@STATUS@", base)
        if depends_on:
            text += _CSH_GUARD.replace("@STATUSDIR@", status_dir).replace("@DEPS@", " ".join(depends_on))
    free_form = str(workflow.general("free_form", "") or "")
    if free_form:
        text += free_form if free_form.endswith("\n") else free_form + "\n"
    if fail_fast(workflow) and fam == "posix":
        text += "set -eo pipefail\n" if shell in PIPEFAIL_SHELLS else "set -e\n"
    return text


def body(steps: list, label: str, workflow: Workflow, before: Optional[list] = None) -> str:
    """The commands of ``steps`` (StepCommands), with progress reporting.

    ``before`` are extra commands run after the first step marker (e.g. the
    clean-up of a rerun).
    """
    fam = family(shell_of(workflow))
    csh_fail_fast = fam == "csh" and fail_fast(workflow)
    out = ["set pm_rc = 0"] if fam == "csh" else []
    for i, step in enumerate(steps):
        if fam == "posix":
            out.append(f"pm_step {step.step}")
        elif fam == "csh":
            out += _csh_step(step.step)
        out.append(f'echo "Running {step.step}: {step.description} on {label}"')
        lines = list(before or []) + list(step.lines) if i == 0 else list(step.lines)
        for line in lines:
            out.append(line)
            if fam == "csh" and line.strip():
                # Every command is checked - including pipelines that start with echo,
                # such as the checksum test `echo "<md5>  <file>" | md5sum -c -`.
                out += _CSH_CHECK if csh_fail_fast else ["set pm_rc = $status"]
    if fam == "csh":
        # Like the EXIT trap of POSIX shells, report the exit status of the last command.
        out += [
            "if ( $pm_rc == 0 ) then",
            '    echo "DONE $PM_STEP $PM_START `date +%s` 0" > "$PM_STATUS"',
            "else",
            '    echo "FAILED $PM_STEP $PM_START `date +%s` $pm_rc" > "$PM_STATUS"',
            "endif",
            "exit $pm_rc",
        ]
    return "\n".join(out) + ("\n" if out else "")


def tracks_status(workflow: Workflow) -> bool:
    return family(shell_of(workflow)) is not None
