"""Remove temporary files left behind by pipemaker builds."""

from __future__ import annotations

import logging
import re
import shutil
from typing import Callable, Optional

from .config import Workspace

log = logging.getLogger(__name__)

_LEGACY_SCRIPT = re.compile(r"^run\d+\.csh$")
_MASTER = ("master.sh", "master.csh")


def _confirm_default(question: str) -> bool:
    return input(f"{question} (y/n) ").strip().lower().startswith("y")


def cleanup(
    workspace: Workspace,
    folders: list,
    experiments: Optional[list] = None,
    keep_scripts: bool = False,
    confirm: Optional[Callable[[str], bool]] = _confirm_default,
) -> list:
    """Delete the files inside ``folders`` of every (or the given) experiment.

    Unless ``keep_scripts`` is set, generated job scripts and master scripts
    are removed too.  ``confirm`` is asked once per folder; pass ``None`` to
    skip confirmation.  Returns the list of removed paths.
    """
    removed = []
    if not keep_scripts:
        for name in _MASTER:
            path = workspace.root / name
            if path.is_file():
                path.unlink()
                removed.append(path)
    if not workspace.output_dir.exists():
        return removed
    exp_dirs = sorted(d for d in workspace.output_dir.iterdir() if d.is_dir())
    if experiments:
        exp_dirs = [d for d in exp_dirs if d.name in experiments]
    for folder in folders:
        if confirm is not None and not confirm(f"Clear '{folder}' in {len(exp_dirs)} experiment folder(s)?"):
            continue
        for exp_dir in exp_dirs:
            target = exp_dir / folder
            if not target.is_dir():
                log.info("Folder %s does not exist in %s", folder, exp_dir.name)
                continue
            for path in target.iterdir():
                if path.is_file() or path.is_symlink():
                    path.unlink()
                    removed.append(path)
                else:
                    log.warning("Not removing directory %s", path)
    if not keep_scripts:
        for exp_dir in exp_dirs:
            scripts = exp_dir / "scripts"
            if scripts.is_dir():
                shutil.rmtree(scripts)
                removed.append(scripts)
            for path in exp_dir.iterdir():
                if path.is_file() and _LEGACY_SCRIPT.match(path.name):
                    path.unlink()
                    removed.append(path)
    return removed
