"""mitopiper: generate HPC job scripts for sequencing-data pipelines.

mitopiper turns a sample metadata table (ENA or SRA) and a configurable
workflow into one job script per sample for SGE, SLURM or LSF, plus a master
script that submits them with the right dependencies.
"""

__version__ = "2.0.0.dev0"
__author__ = "Noam Shtolz"


def build_info() -> str:
    """``commit d0245d6, 2026-10-19`` for the commit this copy was made from: stamped into the release archive
    (``_commit.py``), else read from git when running from a checkout, else ``commit unknown``."""
    import subprocess
    from pathlib import Path

    from ._commit import COMMIT, DATE

    if not COMMIT.startswith("$"):
        return f"commit {COMMIT}, {DATE}"
    try:
        out = subprocess.run(["git", "-C", str(Path(__file__).resolve().parent), "log", "-1", "--format=%h %cs"],
                             capture_output=True, text=True, timeout=10).stdout.split()
    except (OSError, subprocess.TimeoutExpired):
        out = []
    if len(out) == 2:
        return f"commit {out[0]}, {out[1]}"
    return "commit unknown: installed from a git clone (see `git log -1` there) or an unstamped copy"
