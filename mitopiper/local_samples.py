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
import json
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


def _strip_ext(name: str) -> str:
    return FASTQ_RE.sub("", name)


_SEP = re.compile(r"([_.-])")
_MATE = re.compile(r"^([Rr]?)([12])$")
_LANE = re.compile(r"^L\d{3}$")


def _tokens(stem: str) -> list:
    """``["x", "_", "R1", "_", "001"]``: name parts and the separators between them."""
    return _SEP.split(stem)


def _join_without(tokens: list, drop: set) -> str:
    """The name with the tokens at ``drop`` (and the separator before each) left out."""
    keep = []
    for i, token in enumerate(tokens):
        if i in drop or (i + 1 in drop and i % 2 == 1):
            continue
        keep.append(token)
    return "".join(keep).strip("_.-")


def _sample_name(tokens: list, first_mate: int) -> str:
    """The sample: the name before the mate token, without a lane (``_L001``) at its end."""
    head = [t for t in tokens[:max(first_mate - 1, 0)]]
    if len(head) >= 1 and _LANE.match(head[-1]):
        head = head[:-2]
    return "".join(head).strip("_.-") or "".join(tokens)


def _pair_files(stems: dict) -> tuple:
    """Mates among ``stems`` ({stem: path}): ``(pairs, singles)`` with pairs ``[(run, sample, path1, path2)]``.

    FIX [AUDIT-2026-11-26-H1, review]: one pattern (``<name>_R1[_001]``) found the mates; Trim Galore's
    ``x_1_val_1``/``x_2_val_2``, ``x_R1_trimmed`` and ``y-R1`` became two single-end samples, and single-end files
    numbered ``donor_1``, ``donor_2``, ``donor_3`` had the first two paired. Two files are mates when their names
    differ only in mate tokens (1/2, R1/R2), one all 1, the other all 2; bare numbers are not mates when the same
    name also comes with 3 or more."""
    split = {stem: _tokens(stem) for stem in stems}
    candidates: dict = {}                                   # key -> {mate: [(stem, positions)]}
    for stem, tokens in split.items():
        mates = [i for i in range(0, len(tokens), 2) if _MATE.match(tokens[i])]
        for mask in range(1, 1 << len(mates)):
            chosen = [mates[b] for b in range(len(mates)) if mask >> b & 1]
            values = {_MATE.match(tokens[i]).group(2) for i in chosen}
            if len(values) != 1:
                continue
            key = tuple("\0" + _MATE.match(t).group(1).upper() if i in chosen else t for i, t in enumerate(tokens))
            candidates.setdefault(key, {}).setdefault(values.pop(), []).append((stem, tuple(chosen)))
    numbered = {}                                           # names with a number of 3 or more at a token
    for stem, tokens in split.items():
        for i in range(0, len(tokens), 2):
            if tokens[i].isdigit() and int(tokens[i]) >= 3:
                numbered.setdefault(tuple(tokens[:i] + ["\0"] + tokens[i + 1:]), True)
    pairs, used = [], set()
    # the closest names first: fewest differing tokens (run_1_R1 + run_1_R2 before run_1_R1 + run_2_R2), R1/R2 before
    # bare numbers, the rightmost token first
    def priority(item):
        masked = [i for i, t in enumerate(item[0]) if t.startswith("\0")]
        return len(masked), any(item[0][i] == "\0" for i in masked), -max(masked)

    for key, by_mate in sorted(candidates.items(), key=priority):
        ones = [c for c in by_mate.get("1", []) if c[0] not in used]
        twos = [c for c in by_mate.get("2", []) if c[0] not in used]
        if len(ones) != 1 or len(twos) != 1:
            continue
        (stem1, chosen), (stem2, _) = ones[0], twos[0]
        bare = all(not split[stem1][i][:1].isalpha() for i in chosen)
        if bare and tuple("\0" if t.startswith("\0") else t for t in key) in numbered:
            continue                                        # donor_1, donor_2, donor_3: numbers, not mates
        tokens = split[stem1]
        pairs.append((_join_without(tokens, set(chosen)), _sample_name(tokens, min(chosen)), stems[stem1], stems[stem2]))
        used |= {stem1, stem2}
    return pairs, [stem for stem in stems if stem not in used]


def scan_fastq_folder(folder, ignore_paired: bool = False, force_paired: bool = False) -> list:
    """Return a list of ``(run_id, sample_id, [files])`` found in ``folder``."""
    folder = Path(folder).expanduser().resolve()
    files = sorted(p for p in folder.iterdir() if p.is_file() and FASTQ_RE.search(p.name))
    if not files:
        raise FileNotFoundError(f"No FASTQ files found in {folder}")
    stems = {_strip_ext(p.name): str(p) for p in files}
    runs: "OrderedDict[str, dict]" = OrderedDict()
    sample_of: dict = {}
    if ignore_paired:
        singles = list(stems)
    else:
        pairs, singles = _pair_files(stems)
        for run_id, sample, first, second in pairs:
            runs[run_id] = {1: first, 2: second}
            # FIX [AUDIT-2026-11-26-H1]: the sample is the name before the mate without a lane (_L001); the _NNN
            # removed before belonged to the name itself (patient_101, _102, _103 were pooled into one sample)
            sample_of[run_id] = sample
    for stem in singles:
        tokens = _tokens(stem)
        last = max((i for i in range(0, len(tokens), 2) if re.match(r"^[Rr][12]$", tokens[i])), default=None)
        if last is not None and not ignore_paired:          # x_R1_001 alone: a single-end run of sample x
            run_id, sample, mate = _join_without(tokens, {last}), _sample_name(tokens, last), int(tokens[last][1])
        else:
            run_id, sample, mate = stem, stem, 1
        if run_id in runs:
            # FIX [AUDIT-2026-11-26-H1, review]: C.fq.gz next to C_R1/C_R2 was dropped silently
            log.warning("%s has the name of a paired run (%s): it becomes the sample %s_unpaired", stems[stem], run_id,
                        run_id)
            run_id = sample = f"{run_id}_unpaired"
        if mate == 2:
            log.warning("Run %s has only a second mate; treating it as single-end", run_id)
        runs[run_id] = {mate: stems[stem]}
        sample_of[run_id] = sample
    if force_paired and all(len(v) == 1 for v in runs.values()):
        # Pair consecutive files when names carry no mate information.
        flat = [next(iter(v.values())) for v in runs.values()]
        if len(flat) % 2:
            raise ValueError("--force-paired needs an even number of FASTQ files")
        runs = OrderedDict((_strip_ext(Path(flat[i]).name), {1: flat[i], 2: flat[i + 1]}) for i in range(0, len(flat), 2))
        sample_of = {run_id: run_id for run_id in runs}
    result = []
    for run_id, mates in sorted(runs.items(), key=lambda kv: min(kv[1].values())):
        ordered = [mates[k] for k in sorted(mates)]
        result.append((run_id, sample_of.get(run_id, run_id), ordered))
    pooled = {}
    for run_id, sample_id, _ in result:
        pooled.setdefault(sample_id, []).append(run_id)
    for sample_id, run_ids in pooled.items():
        if len(run_ids) > 1:
            log.info("Sample %s gets %d runs (lanes): %s", sample_id, len(run_ids), ", ".join(run_ids))
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
    if mode == "move" and not source.exists() and target.is_file():
        return                                       # moved by an interrupted earlier import of this folder
    if target.exists() or target.is_symlink():
        if target.resolve() == source.resolve() or filecmp.cmp(source, target, shallow=True):
            return                                   # placed by an earlier import
        raise FileExistsError(f"{target} exists already and is not {source}: remove it, or import under another "
                              "name")
    if mode == "link":
        target.symlink_to(source.resolve())
        return
    if mode == "move":
        try:
            os.rename(source, target)                # same file system: atomic
            return
        except OSError:
            pass
    # FIX [AUDIT-2026-11-26-M3, review]: a copy cut short left a partial file under the final name
    part = target.with_name(target.name + ".part")
    shutil.copy2(source, part)
    part.replace(target)
    if mode == "move":
        source.unlink()


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
    # FIX [AUDIT-2026-11-26-M3, review]: a move cut short left part of the files in the workspace; importing the folder
    # again scanned only what was left there, and the samples (or mates) already moved were silently missing from
    # the table. The runs found are written down first, and an import of the same folder resumes from them.
    journal = target / ".mitopiper_import.json"
    resumed = None
    if not dry_run and journal.is_file():
        try:
            record = json.loads(journal.read_text(encoding="utf-8"))
        except ValueError:
            record = {}
        if record.get("source") == str(folder) and record.get("mode") == mode:
            resumed = [(run, sample, list(files)) for run, sample, files in record.get("runs", [])]
            log.info("Resuming the import of %s into %s, which was interrupted", folder, target)
    if resumed is None and mode == "move" and not is_fastq_folder(folder) and is_fastq_folder(target):
        log.info("%s was moved into %s already", folder, target)
        return make_tsv(target, organism, Path(workspace_root) / "metadata" / f"{name}.tsv", library_strategy,
                        ignore_paired, force_paired)
    runs = resumed if resumed is not None else scan_fastq_folder(folder, ignore_paired, force_paired)
    if dry_run:
        import tempfile

        placed = [(run, sample, [str(target / Path(f).name) for f in files]) for run, sample, files in runs]
        log.info("Would %s %d run(s) of %s into %s (dry run: nothing placed)", mode, len(runs), folder, target)
        return write_table(placed, organism, Path(tempfile.mkdtemp(prefix="mitopiper-")) / f"{name}.tsv",
                           library_strategy)
    target.mkdir(parents=True, exist_ok=True)
    journal.write_text(json.dumps({"source": str(folder), "mode": mode, "runs": runs}), encoding="utf-8")
    for _, _, files in runs:
        for file in files:
            _place(Path(file), target / Path(file).name, mode)
    log.info("%s %d run(s) of %s into %s", {"link": "Linked", "copy": "Copied", "move": "Moved"}[mode],
             len(runs), folder, target)
    # FIX [AUDIT-2026-11-26-M3]: the table lists the runs of this folder only: scanning the target added the samples
    # of an earlier import under the same name (e.g. two folders called fastq/), with this import's organism
    placed = [(run, sample, [str(target / Path(f).name) for f in files]) for run, sample, files in runs]
    table = write_table(placed, organism, Path(workspace_root) / "metadata" / f"{name}.tsv", library_strategy)
    journal.unlink(missing_ok=True)
    return table


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
