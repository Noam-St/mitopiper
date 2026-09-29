"""Mitochondrial-genome specifics: circularity and nuclear mitochondrial segments (NUMTs).

Circularity
    The mtDNA is circular but references are linear, so reads that span the
    origin of the sequence (in the D-loop) align poorly to either end. As in
    the original ``circulize.py`` / ``combine_processed_pileups.py`` scripts,
    reads are also aligned to a *rotated* copy of the reference, whose first
    ``shift`` bases were moved to the end so that the origin lies inside the
    sequence. :func:`combine_pileups` then takes the positions near the origin
    (the first ``window`` bases and, by default, the last ``tail_window``
    bases) from the rotated alignment, shifted back to normal coordinates,
    and every other position from the normal alignment.

NUMTs
    Copies of mitochondrial sequence in the nuclear genome attract reads when
    only the mtDNA is used as reference, which looks like heteroplasmy. Reads
    are therefore aligned to the nuclear genome *plus* the mtDNA
    (:func:`combine_references`) and only reads that map to the mtDNA with a
    high mapping quality are kept. Alternatively, `Numt Parser
    <https://github.com/adeflamingh/NuMt_parser>`_ (de Flamingh et al. 2023,
    J. Hered. 114:120) classifies every read as cytoplasmic mtDNA (``cymt``) or
    NUMT by comparing its alignments to the mtDNA and to a NUMT reference.
    :func:`numt_reference` derives that reference from the nuclear genome
    (where fragments of the mtDNA align) and :func:`numt_read_list` turns
    Numt Parser's table into the list of reads to keep.

Per-position annotation and heteroplasmy of the resulting pileups are in
:mod:`mitopiper.pileup`.

The functions are also available as ``mitopiper mt ...`` commands, which the
bundled mitochondrial presets call from their jobs.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Iterator, Optional, TextIO

MITO_NAMES = {"chrm", "mt", "m", "chrmt", "mito", "mitochondrion", "chrmito"}


def read_fasta(handle: TextIO) -> Iterator[tuple]:
    """Yield ``(header, sequence)`` of every record (header without ``>``)."""
    header, chunks = None, []
    for line in handle:
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def write_record(handle: TextIO, header: str, sequence: str, width: int = 60) -> None:
    handle.write(f">{header}\n")
    for i in range(0, len(sequence), width):
        handle.write(sequence[i:i + width] + "\n")


def contigs(fasta) -> list:
    """Sequence names (first word of each header) of a FASTA file."""
    with open(fasta, encoding="utf-8") as handle:
        return [line[1:].split()[0] for line in handle if line.startswith(">") and line[1:].strip()]


def sequence_length(fasta) -> int:
    """Length of the first sequence of a FASTA file."""
    with open(fasta, encoding="utf-8") as handle:
        for _, sequence in read_fasta(handle):
            return len(sequence)
    raise ValueError(f"{fasta} contains no sequence")


# ------------------------------------------------------------------ circularity
def rotate_fasta(fasta, output, shift: int = 500) -> int:
    """Write ``fasta`` with the first ``shift`` bases of its first sequence moved to the end.

    The header is kept, so alignments to the rotated copy report the same
    contig name. Other records (if any) are copied unchanged. Returns the
    length of the rotated sequence.
    """
    with open(fasta, encoding="utf-8") as handle:
        records = list(read_fasta(handle))
    if not records:
        raise ValueError(f"{fasta} contains no sequence")
    header, sequence = records[0]
    if not 0 < shift < len(sequence):
        raise ValueError(f"shift must be between 1 and {len(sequence) - 1} (the sequence has {len(sequence)} bases)")
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as out:
        write_record(out, header, sequence[shift:] + sequence[:shift])
        for other_header, other in records[1:]:
            write_record(out, other_header, other)
    return len(sequence)


def rotated_to_normal(position: int, length: int, shift: int) -> int:
    """Position on the normal reference of 1-based ``position`` on the rotated one."""
    return position + shift if position <= length - shift else position - (length - shift)


def combine_pileups(normal, rotated, output, length: int, shift: int = 500, window: int = 400,
                    tail_window: Optional[int] = None) -> dict:
    """Merge the pileups of the normal and the rotated alignment into one, in normal coordinates.

    Positions ``1..window`` and ``length-tail_window+1..length`` (``tail_window``
    defaults to ``window``) come from the rotated alignment, where reads that span the
    origin align; all others come from the normal alignment. ``window`` may not exceed
    ``shift``, and ``tail_window`` must stay below ``length - shift``, so that the positions
    taken from the rotated alignment lie inside its sequence rather than at its ends.

    Returns counts of the lines taken from each input.
    """
    tail_window = window if tail_window is None else tail_window
    if window > shift:
        raise ValueError(f"window ({window}) must not exceed the rotation shift ({shift})")
    if tail_window >= length - shift:
        raise ValueError(f"tail window ({tail_window}) must be smaller than length - shift ({length - shift})")

    def from_rotated(position: int) -> bool:
        return position <= window or position > length - tail_window

    contig = None
    rows = []
    counts = {"normal": 0, "rotated": 0, "other": 0}
    with open(normal, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[1].isdigit():
                continue
            contig = contig or fields[0]
            if fields[0] != contig:  # another sequence (not the mtDNA): keep as is
                rows.append(((1, fields[0], int(fields[1])), line))
                counts["other"] += 1
            elif not from_rotated(int(fields[1])):
                rows.append(((0, fields[0], int(fields[1])), line))
                counts["normal"] += 1
    with open(rotated, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[1].isdigit():
                continue
            if contig is not None and fields[0] != contig:
                continue
            contig = contig or fields[0]
            position = rotated_to_normal(int(fields[1]), length, shift)
            if from_rotated(position):
                fields[1] = str(position)
                rows.append(((0, fields[0], position), "\t".join(fields) + "\n"))
                counts["rotated"] += 1
    rows.sort(key=lambda row: row[0])
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as out:
        out.writelines(line for _, line in rows)
    return counts


# ----------------------------------------------------------------------- NUMTs
def is_mitochondrial(header: str, mt_names: set) -> bool:
    name = header.split()[0]
    return (name.lower() in MITO_NAMES or name in mt_names or name.split(".")[0] in mt_names
            or re.search(r"mitochondri", header, re.IGNORECASE) is not None)


def combine_references(nuclear, mt, output, log: TextIO = sys.stderr) -> list:
    """Write the nuclear genome followed by the mtDNA, for NUMT-aware alignment.

    Mitochondrial sequences already present in the nuclear assembly (``chrM``, ``MT``,
    headers mentioning "mitochondrion", or the mtDNA accession) are left out: a
    duplicate copy would make every mtDNA read map twice with mapping quality 0.
    The nuclear FASTA is streamed, so whole genomes need little memory.
    Returns the names of the sequences that were left out.
    """
    with open(mt, encoding="utf-8") as handle:
        mt_records = list(read_fasta(handle))
    if not mt_records:
        raise ValueError(f"{mt} contains no sequence")
    mt_names = {h.split()[0] for h, _ in mt_records} | {h.split()[0].split(".")[0] for h, _ in mt_records}
    dropped = []
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(nuclear, encoding="utf-8") as src, open(output, "w", encoding="utf-8") as out:
        keep = True
        for line in src:
            if line.startswith(">"):
                header = line[1:].strip()
                keep = not is_mitochondrial(header, mt_names)
                if not keep:
                    dropped.append(header.split()[0])
                    print(f"[mitopiper] leaving out {header.split()[0]} (mitochondrial) of {nuclear}", file=log)
            if keep:
                out.write(line if line.endswith("\n") else line + "\n")
        for header, sequence in mt_records:
            write_record(out, header, sequence)
    return dropped


def fragments(mt, output, size: int = 100, step: int = 20) -> int:
    """Write overlapping ``size`` bp fragments of the (circular) mtDNA, every ``step`` bases.

    Aligned to the nuclear genome, they reveal the NUMTs (:func:`numt_regions`).
    Fragments near the end wrap around to the start. Returns the number written.
    """
    with open(mt, encoding="utf-8") as handle:
        records = list(read_fasta(handle))
    if not records:
        raise ValueError(f"{mt} contains no sequence")
    if size < 20 or step < 1:
        raise ValueError("fragment size must be at least 20 and step at least 1")
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(output, "w", encoding="utf-8") as out:
        for header, sequence in records:
            name = header.split()[0]
            circular = sequence * (size // max(len(sequence), 1) + 2)
            for start in range(0, len(sequence), step):
                write_record(out, f"{name}_{start + 1}", circular[start:start + size])
                count += 1
    return count


_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def _reference_span(cigar: str) -> int:
    return sum(int(n) for n, op in _CIGAR.findall(cigar) if op in "MDN=X")


def numt_regions(sam: TextIO, mt_names: set, min_length: int = 50, flank: int = 100, merge: int = 200) -> list:
    """Nuclear regions (``contig:start-end``) where mtDNA fragments align, from SAM text.

    Hits on mitochondrial contigs and shorter than ``min_length`` are ignored;
    hits closer than ``merge`` bases are joined and ``flank`` bases are added on
    both sides (clipped to the contig when its length is in the header).
    """
    lengths, hits = {}, {}
    for line in sam:
        if line.startswith("@SQ"):
            tags = dict(f.split(":", 1) for f in line.rstrip("\n").split("\t")[1:] if ":" in f)
            if "SN" in tags and tags.get("LN", "").isdigit():
                lengths[tags["SN"]] = int(tags["LN"])
            continue
        if line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6 or not fields[1].isdigit() or int(fields[1]) & 4 or fields[2] == "*":
            continue
        contig, start, span = fields[2], int(fields[3]), _reference_span(fields[5])
        if span < min_length or is_mitochondrial(contig, mt_names):
            continue
        hits.setdefault(contig, []).append((start, start + span - 1))
    regions = []
    for contig in sorted(hits):
        merged = []
        for start, end in sorted(hits[contig]):
            if merged and start <= merged[-1][1] + merge:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        for start, end in merged:
            end = min(end + flank, lengths[contig]) if contig in lengths else end + flank
            regions.append(f"{contig}:{max(1, start - flank)}-{end}")
    return regions


def numt_reference(sam, mt, output, **options) -> int:
    """Write the NUMT regions found in ``sam`` (mtDNA fragments aligned to the nuclear genome)
    to ``output``, one ``contig:start-end`` per line for ``samtools faidx -r``. Returns their number."""
    with open(mt, encoding="utf-8") as handle:
        names = {h.split()[0] for h, _ in read_fasta(handle)}
    names |= {n.split(".")[0] for n in names}
    with open(sam, encoding="utf-8") as handle:
        regions = numt_regions(handle, names, **options)
    if not regions:
        raise ValueError(f"no NUMT found in {sam}: the nuclear genome has no copy of the mtDNA fragments; "
                         "use the mt-variants preset (no NUMT filtering) or give a NUMT FASTA with refs add --numt-fasta")
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as out:
        out.writelines(r + "\n" for r in regions)
    return len(regions)


def numt_read_list(table, output, keep: str = "cymt", pairs: str = "all") -> dict:
    """Names of the reads to keep from a Numt Parser output table.

    ``keep`` is ``cymt`` (only reads classified as cytoplasmic mtDNA) or
    ``not-numt`` (also the ``undetermined`` ones). With ``pairs`` ``all`` a read
    is kept only if every row of its name (both mates of a pair) qualifies; with
    ``any`` one qualifying mate is enough. Returns counts per class.
    """
    wanted = {"cymt"} if keep == "cymt" else {"cymt", "undetermined"}
    verdicts, counts = {}, {}
    with open(table, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) < 9:
                continue
            name, candidate = fields[0], fields[8]
            counts[candidate] = counts.get(candidate, 0) + 1
            verdicts.setdefault(name, []).append(candidate in wanted)
    combine = all if pairs == "all" else any
    kept = [name for name, votes in verdicts.items() if combine(votes)]
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as out:
        out.writelines(name + "\n" for name in kept)
    counts["kept"] = len(kept)
    return counts


# ------------------------------------------------------------------------- CLI
def add_parser(sub) -> None:
    """Register ``mitopiper mt ...``."""
    mt = sub.add_parser("mt", help="mitochondrial helpers used by the mt presets (rotation, pileup merge, NUMTs)")
    msub = mt.add_subparsers(dest="action", required=True)
    rot = msub.add_parser("rotate", help="move the first SHIFT bases of the mtDNA to its end")
    rot.add_argument("fasta")
    rot.add_argument("output")
    rot.add_argument("--shift", type=int, default=500)
    comb = msub.add_parser("combine-pileups", help="merge normal and rotated pileups into one circular-aware pileup")
    comb.add_argument("normal")
    comb.add_argument("rotated")
    comb.add_argument("output")
    comb.add_argument("--reference", required=True, help="the normal (not rotated) mtDNA FASTA, for its length")
    comb.add_argument("--shift", type=int, default=500)
    comb.add_argument("--window", type=int, default=400, help="positions at the start taken from the rotated pileup")
    comb.add_argument("--tail-window", type=int, help="positions at the end taken from the rotated pileup (default: --window)")
    refs = msub.add_parser("combine-references", help="nuclear genome + mtDNA (without duplicate mt sequences)")
    refs.add_argument("nuclear")
    refs.add_argument("mt")
    refs.add_argument("output")
    frag = msub.add_parser("fragments", help="overlapping fragments of the mtDNA, to find NUMTs in a nuclear genome")
    frag.add_argument("mt")
    frag.add_argument("output")
    frag.add_argument("--size", type=int, default=100)
    frag.add_argument("--step", type=int, default=20)
    regions = msub.add_parser("numt-regions", help="nuclear regions where mtDNA fragments align (for samtools faidx -r)")
    regions.add_argument("sam", help="the fragments aligned to the nuclear genome (SAM)")
    regions.add_argument("mt", help="the mtDNA FASTA (its contigs are ignored)")
    regions.add_argument("output")
    regions.add_argument("--min-length", type=int, default=50)
    regions.add_argument("--flank", type=int, default=100)
    regions.add_argument("--merge", type=int, default=200)
    reads = msub.add_parser("numt-reads", help="reads to keep from a Numt Parser table (for samtools view -N)")
    reads.add_argument("table")
    reads.add_argument("output")
    reads.add_argument("--keep", choices=("cymt", "not-numt"), default="cymt",
                       help="cymt: reads classified as mtDNA; not-numt: also the undetermined ones")
    reads.add_argument("--pairs", choices=("all", "any"), default="all",
                       help="all: every mate must qualify; any: one qualifying mate keeps the pair")
    from . import pileup

    pileup.add_parsers(msub)
    names = msub.add_parser("contigs", help="print the sequence names of a FASTA")
    names.add_argument("fasta")


def run(args) -> int:
    if args.action in ("annotate", "heteroplasmy"):
        from . import pileup

        return pileup.run(args)
    if args.action == "rotate":
        rotate_fasta(args.fasta, args.output, args.shift)
    elif args.action == "combine-pileups":
        counts = combine_pileups(args.normal, args.rotated, args.output, sequence_length(args.reference),
                                 args.shift, args.window, args.tail_window)
        print(f"[mitopiper] {counts['rotated']} positions from the rotated alignment, {counts['normal']} from the normal one",
              file=sys.stderr)
    elif args.action == "combine-references":
        combine_references(args.nuclear, args.mt, args.output)
    elif args.action == "fragments":
        fragments(args.mt, args.output, args.size, args.step)
    elif args.action == "numt-regions":
        try:
            count = numt_reference(args.sam, args.mt, args.output, min_length=args.min_length,
                                   flank=args.flank, merge=args.merge)
        except ValueError as exc:
            print(f"[mitopiper] {exc}", file=sys.stderr)
            return 1
        print(f"[mitopiper] {count} NUMT regions", file=sys.stderr)
    elif args.action == "numt-reads":
        counts = numt_read_list(args.table, args.output, args.keep, args.pairs)
        print("[mitopiper] Numt Parser: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())), file=sys.stderr)
    elif args.action == "contigs":
        print(" ".join(contigs(args.fasta)))
    return 0
