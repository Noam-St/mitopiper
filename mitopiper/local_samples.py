"""Create a metadata table for FASTQ files that are already on disk.

The generated TSV follows the ENA column names, with absolute file paths in
``fastq_ftp``.  When such a table is built, the ``wget`` step links (or
gzips) the local files into the experiment instead of downloading them, so
no files need to be copied or renamed up-front.

Illumina lane files (``..._L001_R1_001.fastq.gz``) of the same sample are
grouped into one sample and concatenated by the download step.

:func:`import_fastq_folder` brings a folder of FASTQ files into a workspace (``fastq/<name>/``, as
links, copies or by moving them) and writes its table to ``metadata/<name>.tsv``, so that ``build`` and
the interactive mode can take a folder instead of a metadata table or an accession.
:func:`browse_for_folder` opens a graphical folder chooser for that, where a display is available.
"""

from __future__ import annotations

import csv
import filecmp
import logging
import os
import re
import shutil
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)

FASTQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$")
MATE_RE = re.compile(r"^(?P<stem>.+?)[_.]R?(?P<mate>[12])(?P<tail>_\d{3})?$")
LANE_RE = re.compile(r"_L\d{3}")


def _strip_ext(name: str) -> str:
    return FASTQ_RE.sub("", name)


def scan_fastq_folder(folder, ignore_paired: bool = False, force_paired: bool = False) -> list:
    """Return a list of ``(run_id, sample_id, [files])`` found in ``folder``."""
    folder = Path(folder).expanduser().resolve()
    files = sorted(p for p in folder.iterdir() if p.is_file() and FASTQ_RE.search(p.name))
    if not files:
        raise FileNotFoundError(f"No FASTQ files found in {folder}")
    runs: "OrderedDict[str, dict]" = OrderedDict()
    for path in files:
        stem = _strip_ext(path.name)
        match = MATE_RE.match(stem)
        if match and not ignore_paired:
            run_id = match.group("stem") + (match.group("tail") or "")
            runs.setdefault(run_id, {})[int(match.group("mate"))] = str(path)
        else:
            runs.setdefault(stem, {})[1] = str(path)
    if force_paired and all(len(v) == 1 for v in runs.values()):
        # Pair consecutive files when names carry no mate information.
        flat = [v[1] for v in runs.values()]
        if len(flat) % 2:
            raise ValueError("--force-paired needs an even number of FASTQ files")
        runs = OrderedDict((_strip_ext(Path(flat[i]).name), {1: flat[i], 2: flat[i + 1]}) for i in range(0, len(flat), 2))
    result = []
    for run_id, mates in runs.items():
        if 2 in mates and 1 not in mates:
            log.warning("Run %s has only a second mate; treating it as single-end", run_id)
        ordered = [mates[k] for k in sorted(mates)]
        sample_id = LANE_RE.sub("", run_id)
        sample_id = re.sub(r"_\d{3}$", "", sample_id)
        result.append((run_id, sample_id, ordered))
    return result


def make_tsv(folder, organism: str, output: Optional[Path] = None, library_strategy: str = "",
             ignore_paired: bool = False, force_paired: bool = False) -> Path:
    """Write ``<folder name>.tsv`` (or ``output``) describing the FASTQ files in ``folder``."""
    folder = Path(folder).expanduser().resolve()
    runs = scan_fastq_folder(folder, ignore_paired, force_paired)
    return write_table(runs, organism, Path(output) if output else Path.cwd() / f"{folder.name}.tsv", library_strategy)


def write_table(runs: list, organism: str, output: Path, library_strategy: str = "") -> Path:
    """Write the ``(run_id, sample_id, [files])`` of :func:`scan_fastq_folder` as an ENA-style table."""
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["run_accession", "sample_accession", "scientific_name", "fastq_ftp", "sample_title", "library_strategy"]
    with open(output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for run_id, sample_id, files in runs:
            writer.writerow(
                {
                    "run_accession": run_id,
                    "sample_accession": sample_id,
                    "scientific_name": organism,
                    "fastq_ftp": ";".join(files),
                    "sample_title": sample_id,
                    "library_strategy": library_strategy,
                }
            )
    log.info("Wrote %d run(s) to %s", len(runs), output)
    return output


# ------------------------------------------------------------- folders into a workspace
#: How import_fastq_folder places the files in the workspace.
FASTQ_MODES = ("link", "copy", "move")


def is_fastq_folder(path) -> bool:
    """Whether ``path`` is a folder holding FASTQ files (.fastq/.fq, gzipped or not)."""
    path = Path(path).expanduser()
    return path.is_dir() and any(p.is_file() and FASTQ_RE.search(p.name) for p in path.iterdir())


def _place(source: Path, target: Path, mode: str) -> None:
    """Link, copy or move ``source`` to ``target``; a target that already is that file (or an identical
    copy of it) is left alone."""
    if target.exists() or target.is_symlink():
        if target.resolve() == source.resolve() or filecmp.cmp(source, target, shallow=True):
            return                                   # placed by an earlier import
        raise FileExistsError(f"{target} exists already and is not {source}: remove it, or import under another "
                              "name")
    if mode == "link":
        target.symlink_to(source.resolve())
    elif mode == "copy":
        shutil.copy2(source, target)
    else:
        shutil.move(str(source), str(target))


def import_fastq_folder(workspace_root, folder, organism: str, mode: str = "link", name: Optional[str] = None,
                        library_strategy: str = "", ignore_paired: bool = False, force_paired: bool = False,
                        dry_run: bool = False) -> Path:
    """Place the FASTQ files of ``folder`` in ``<workspace>/fastq/<name>/`` (``mode``: link, copy or move;
    default ``name``: the folder's) and write their table to ``<workspace>/metadata/<name>.tsv``, which is
    returned. Links keep the files where they are; ``move`` takes them out of ``folder``.

    With ``dry_run`` nothing is placed or written in the workspace: the table, with the paths the files would
    have there, goes to a temporary file."""
    if mode not in FASTQ_MODES:
        raise ValueError(f"unknown mode '{mode}': choose one of {', '.join(FASTQ_MODES)}")
    if not organism.strip():
        raise ValueError("the organism of the FASTQ files is needed (e.g. --org 'Homo sapiens')")
    folder = Path(folder).expanduser().resolve()
    name = name or folder.name
    target = Path(workspace_root) / "fastq" / name
    if mode == "move" and not is_fastq_folder(folder) and is_fastq_folder(target):
        log.info("%s was moved into %s already", folder, target)
        return make_tsv(target, organism, Path(workspace_root) / "metadata" / f"{name}.tsv", library_strategy,
                        ignore_paired, force_paired)
    runs = scan_fastq_folder(folder, ignore_paired, force_paired)
    if dry_run:
        import tempfile

        placed = [(run, sample, [str(target / Path(f).name) for f in files]) for run, sample, files in runs]
        log.info("Would %s %d run(s) of %s into %s (dry run: nothing placed)", mode, len(runs), folder, target)
        return write_table(placed, organism, Path(tempfile.mkdtemp(prefix="mitopiper-")) / f"{name}.tsv",
                           library_strategy)
    target.mkdir(parents=True, exist_ok=True)
    for _, _, files in runs:
        for file in files:
            _place(Path(file), target / Path(file).name, mode)
    log.info("%s %d run(s) of %s into %s", {"link": "Linked", "copy": "Copied", "move": "Moved"}[mode],
             len(runs), folder, target)
    return make_tsv(target, organism, Path(workspace_root) / "metadata" / f"{name}.tsv", library_strategy,
                    ignore_paired, force_paired)


def browse_for_folder(title: str = "Choose the folder with the FASTQ files", start: Optional[str] = None) -> str:
    """A folder chosen in a graphical dialog: Tk, else zenity or kdialog. Returns "" when cancelled.

    Raises RuntimeError when no dialog can be shown (no display, e.g. ssh without -X, or no GUI toolkit)."""
    if sys.platform.startswith("linux") and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")):
        raise RuntimeError("no display to show a file browser on (on a cluster: connect with ssh -X or -Y)")
    start = str(Path(start or Path.cwd()).expanduser())
    try:
        import tkinter
        from tkinter import filedialog

        root = tkinter.Tk()
        root.withdraw()
        try:
            return filedialog.askdirectory(title=title, initialdir=start, mustexist=True) or ""
        finally:
            root.destroy()
    except Exception as exc:                          # no tkinter, or Tk cannot open the display
        tk_error = exc
    for command in (["zenity", "--file-selection", "--directory", f"--title={title}", f"--filename={start}/"],
                    ["kdialog", "--getexistingdirectory", start, "--title", title]):
        if shutil.which(command[0]):
            done = subprocess.run(command, capture_output=True, text=True)
            return done.stdout.strip() if done.returncode == 0 else ""
    raise RuntimeError(f"no file browser available (tkinter: {tk_error}; zenity and kdialog not found)")
