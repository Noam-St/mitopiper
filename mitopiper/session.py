"""Interactive monitor session: launch -> watch -> diagnose -> rerun -> watch ...

Used by the interactive mode right after the job scripts are generated (and
from its start menu for an earlier build).  The table refreshes on its own;
meanwhile you can type a command and press Enter:

    Enter     refresh now            e         diagnose every failed job
    j NAME    details of one job     r         rerun failed jobs (automatic fixes)
    o         rerun with options     i SEC     refresh interval
    h         help                   q         leave (jobs keep running)
"""

from __future__ import annotations

import datetime as _dt
import os
import select
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Callable, Optional

from .config import Workspace
from .diagnose import diagnose
from .monitor import (
    DEFAULT_ERROR_PATTERNS,
    FAILED,
    FINAL_STATES,
    KILLED,
    NOT_SUBMITTED,
    PENDING,
    SKIPPED,
    UNKNOWN,
    ErrorScanner,
    _contexts,
    collect,
    query_scheduler,
    render_detail,
    render_table,
    rerun_hint,
    summary,
    tail,
)

HELP = """Commands (type and press Enter):
  Enter      refresh now                 e        diagnose every failed job
  j NAME     details of a job (name/id)  r        rerun failed jobs with the automatic fixes
  o          rerun with your own options i SEC    change the refresh interval
  h          this help                   q        leave the monitor (jobs keep running)
"""
PROBLEM_STATES = (FAILED, KILLED, SKIPPED, UNKNOWN)


def launch(script: Path, log: Path) -> subprocess.Popen:
    """Run a submit script in the background (it survives leaving mitopiper)."""
    handle = open(log, "w", encoding="utf-8")
    return subprocess.Popen(["bash", str(script)], stdout=handle, stderr=subprocess.STDOUT,
                            cwd=str(script.parent), start_new_session=True)


def stdin_poll(timeout: float) -> Optional[str]:
    """Wait up to ``timeout`` seconds for a line on stdin; ``None`` on timeout."""
    try:
        ready, _, _ = select.select([sys.stdin], [], [], timeout)
    except (ValueError, OSError, TypeError):  # no selectable stdin (e.g. Windows): block instead
        ready = [sys.stdin]
    if not ready:
        return None
    line = sys.stdin.readline()
    return "q" if line == "" else line.strip()


class MonitorSession:
    def __init__(
        self,
        workspace: Workspace,
        experiments: list,
        interval: float = 15.0,
        poll: Callable = stdin_poll,
        ask: Callable = input,
        stream=None,
        query: Callable = query_scheduler,
        color: Optional[bool] = None,
        launcher: Callable = launch,
        stay_when_done: bool = False,
    ):
        self.ws = workspace
        self.stay_when_done = stay_when_done   # 'monitor -i': keep offering j/e/q once every job is done
        self.experiments = experiments
        self.interval = interval
        self.poll = poll
        self.ask = ask
        self.out = stream or sys.stdout
        self.query = query
        self.launcher = launcher
        if color is None:
            color = hasattr(self.out, "isatty") and self.out.isatty() and os.environ.get("NO_COLOR") is None
        self.color = color
        self.scanner = ErrorScanner(DEFAULT_ERROR_PATTERNS)
        self.contexts = _contexts(workspace, experiments)
        self.process: Optional[subprocess.Popen] = None
        self.process_log: Optional[Path] = None
        self.process_reported = False
        self.statuses: list = []

    # ------------------------------------------------------------- helpers
    def write(self, text: str = "") -> None:
        self.out.write(text if text.endswith("\n") else text + "\n")
        self.out.flush()

    def start(self, script: Path, log_name: str) -> None:
        self.process_log = self.ws.root / log_name
        self.process = self.launcher(script, self.process_log)
        self.process_reported = False

    def submitting(self) -> bool:
        return self.process is not None and self.process.poll() is None

    def refresh(self) -> None:
        self.statuses = collect(self.ws, self.experiments, self.scanner, self.query)
        if self.submitting():  # master.sh has not reached these jobs yet
            for job in self.statuses:
                if job.state == NOT_SUBMITTED:
                    job.state, job.scheduler_state = PENDING, "waiting for master.sh"

    def launch_problem(self) -> Optional[str]:
        """Message if the submit script failed (once)."""
        if self.process is None or self.process_reported or self.process.poll() in (None, 0):
            return None
        self.process_reported = True
        lines = tail(self.process_log, 15) if self.process_log else []
        text = "\n".join("   " + l for l in lines)
        hint = ""
        if any("command not found" in l for l in lines):
            hint = ("\nThe scheduler's submit command is not available here: run mitopiper on a login node of the "
                    "cluster, or rebuild with the scheduler your cluster uses (sge/slurm/lsf).")
        return f"Submitting failed (exit code {self.process.returncode}); end of {self.process_log}:\n{text}{hint}\n"

    def draw(self, footer: str = "") -> None:
        width = shutil.get_terminal_size((120, 40)).columns
        stamp = _dt.datetime.now().strftime("%H:%M:%S")
        screen = "\033[H\033[J" if self.color else ""
        state = " - submitting jobs..." if self.submitting() else ""
        body = render_table(self.statuses, self.color, width)
        self.out.write(f"{screen}mitopiper monitor - {stamp} - refreshing every {self.interval:g}s{state}\n\n{body}")
        problem = self.launch_problem()
        if problem:
            self.out.write("\n" + problem)
        self.write(footer)

    def pause(self) -> None:
        self.ask("\nPress Enter to return to the monitor ")
        if not self.color:
            self.out.write("\n")

    def find(self, key: str):
        return [s for s in self.statuses if key in (s.name, s.job_id)]

    # ----------------------------------------------------------- commands
    def show_errors(self) -> None:
        problems = [s for s in self.statuses if s.state in PROBLEM_STATES or s.error_count]
        if not problems:
            self.write("\nNo failed jobs and no errors in stderr.")
        for job in problems:
            sug = diagnose(job, self.contexts.get(job.experiment)) if job.state in PROBLEM_STATES else None
            self.write("\n" + render_detail(job, 10, sug))
        self.pause()

    def show_job(self, key: str) -> None:
        jobs = self.find(key)
        if not jobs:
            self.write(f"No job named '{key}'.")
        for job in jobs:
            sug = diagnose(job, self.contexts.get(job.experiment)) if job.state in PROBLEM_STATES else None
            self.write("\n" + render_detail(job, 20, sug))
        self.pause()

    def rerun(self, custom: bool = False) -> bool:
        """Plan, confirm and submit a rerun; returns True if something was submitted."""
        from .rerun import RerunOptions, describe, plan, write

        options = RerunOptions(auto_fix=not custom)
        if custom:
            options = self.ask_rerun_options()
            if options is None:
                return False
        try:
            items = plan(self.ws, self.experiments, options, query=self.query, scanner=self.scanner)
        except Exception as exc:  # noqa: BLE001 - shown to the user, session continues
            self.write(f"Cannot rerun: {exc}")
            self.pause()
            return False
        if not items:
            self.write("Nothing to rerun: no failed, killed or skipped jobs.")
            self.pause()
            return False
        self.write(f"\n{len(items)} job(s) will be resubmitted:\n{describe(items)}")
        if self.ask("Submit them now? (Y/n) ").strip().lower().startswith("n"):
            return False
        try:
            script = write(self.ws, items, options)
        except Exception as exc:  # noqa: BLE001
            self.write(f"Cannot rerun: {exc}")
            self.pause()
            return False
        self.start(script, "rerun.log")
        self.write(f"Resubmitting with {script} ...")
        time.sleep(0.5)
        return True

    def ask_rerun_options(self):
        from .rerun import RerunOptions

        failed = [s for s in self.statuses if s.state in PROBLEM_STATES]
        self.write("\nRerun options (press Enter to skip a question)")
        opts = RerunOptions()
        names = self.ask(f"Jobs to rerun [all failed: {' '.join(s.name for s in failed) or '-'}]: ").split()
        opts.jobs = names or None
        steps = sorted({st for c in self.contexts.values() for st in c.steps})
        start = self.ask(f"Resume at (failed step | 'start' | one of {', '.join(steps)}): ").strip()
        if start == "start":
            opts.from_start = True
        elif start:
            opts.from_step = start
        opts.clean = self.ask("Delete the outputs of the resume step first? (y/N) ").strip().lower().startswith("y")
        opts.memory = self.ask("Memory per job (e.g. 64G): ").strip() or None
        opts.walltime = self.ask("Time limit (e.g. 3-00:00:00): ").strip() or None
        hints = sorted({h for s in failed for h in self._hints(s)})
        if hints:
            line = self.ask(f"Extra setup line [suggested: {'; '.join(hints)}] (Enter = use it, '-' = none, or type your own): ").strip()
            line = "\n".join(hints) if line == "" else ("" if line == "-" else line)
        else:
            line = self.ask("Extra setup line (e.g. module load samtools): ").strip()
        opts.free_form = [line] if line else []
        opts.auto_fix = self.ask("Also apply the automatic fixes? (Y/n) ").strip().lower() not in ("n", "no")
        return opts

    def _hints(self, job) -> list:
        return [s.rerun["free_form_hint"] for s in diagnose(job, self.contexts.get(job.experiment)) if "free_form_hint" in s.rerun]

    def handle(self, command: str) -> bool:
        """Run a command; returns False when the user wants to leave."""
        cmd, _, arg = command.strip().partition(" ")
        cmd = cmd.lower()
        if cmd in ("q", "quit", "exit"):
            return False
        if cmd in ("h", "help", "?"):
            self.write(HELP)
            self.pause()
        elif cmd == "e":
            self.show_errors()
        elif cmd == "j":
            if arg:
                self.show_job(arg.strip())
            else:
                self.write("Usage: j JOB_NAME")
        elif cmd == "r":
            self.rerun(custom=False)
        elif cmd == "o":
            self.rerun(custom=True)
        elif cmd == "i":
            try:
                self.interval = max(1.0, float(arg))
            except ValueError:
                self.write("Usage: i SECONDS")
        elif cmd and self.find(command.strip()):
            self.show_job(command.strip())
        elif cmd:
            self.write(f"Unknown command '{command}'. Type h for help.")
        return True

    # ---------------------------------------------------------------- main
    def run(self) -> int:
        """Monitor until the user leaves; returns 1 if jobs are left failed."""
        try:
            while True:
                self.refresh()
                finished = all(s.state in FINAL_STATES for s in self.statuses) and not self.submitting()
                if not finished:
                    self.draw("Commands: Enter refresh | j NAME details | e errors | r rerun failed | o rerun with options | q quit")
                    command = self.poll(self.interval)
                    if command is not None and not self.handle(command):
                        return self._leave()
                    continue
                problems = [s for s in self.statuses if s.state in PROBLEM_STATES]
                if self.process is not None and self.process.returncode not in (None, 0) and not problems:
                    self.draw("The jobs could not be submitted: fix the problem above, then submit again "
                              f"with 'bash {self.ws.root / 'master.sh'}'.")
                    return 1
                if not problems:
                    self.draw("All jobs finished successfully.")
                    if not self.stay_when_done:
                        return 0
                    answer = self.ask("j NAME details | e errors | q quit\n> ")
                    if not self.handle(answer or "q"):
                        return self._leave(done=True)
                    continue
                failed = [s for s in problems if s.state in (FAILED, KILLED)]
                hint = rerun_hint(failed, self.contexts, shortcuts=False) if failed else ""
                advice = ""
                if any(self._hints(s) for s in failed):
                    advice = ("Some failures need a change in the job environment: choose 'o' to rerun with the "
                              "suggested setup line (e.g. a module load).\n")
                self.draw(f"All jobs have finished; {len(problems)} need attention.\n{hint}{advice}")
                answer = self.ask("r rerun with automatic fixes | o rerun with options | e errors | j NAME details | q quit\n> ")
                if not self.handle(answer or "q"):
                    return self._leave(done=True)
        except KeyboardInterrupt:
            return self._leave()

    def _leave(self, done: bool = False) -> int:
        counts = ", ".join(f"{n} {s.lower()}" for s, n in summary(self.statuses).items())
        self.write(f"\nLeaving the monitor ({counts}).")
        if not done:
            self.write("Submitted jobs keep running.")
        self.write("Come back any time with: mitopiper monitor -i   (or start mitopiper and choose 'monitor')")
        return 1 if any(s.state in PROBLEM_STATES for s in self.statuses) else 0
