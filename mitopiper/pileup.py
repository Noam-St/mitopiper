"""Per-position annotation of mtDNA pileups, and heteroplasmy.

A port of the original ``process_pileup.py`` (pandas + Biopython) to the
standard library. For every position of a ``samtools mpileup`` file it
reports the coverage per strand, the reads starting and ending there, z-scores
of the coverage, optional RPM values (with an htseq-count table), and the gene
covering the position (from a GFF3/GTF or GenBank annotation). It adds allele
counts and the minor allele frequency, so heteroplasmic positions can be
called, and :func:`collect_heteroplasmy` gathers them over all samples.

Differences from ``process_pileup.py``, all deliberate:

* The pileup base string is tokenised instead of matched with regular
  expressions. A read end (``$``) belongs to the base *before* it, and the
  mapping-quality character after ``^`` is skipped, so it is never counted as a
  base; indel sequences (``+2AC``) are not counted as coverage.
* Gene intervals cover exactly their GFF/GenBank coordinates (the original
  included one extra position at the end). Overlapping genes (e.g. ATP8/ATP6)
  are all reported, separated by ``;``, instead of the last one winning.
* z-scores of a contig with constant coverage are left empty instead of NaN.

"pos" and "neg" are the forward (``.``/upper case) and reverse (``,``/lower
case) strand of the alignment, swapped with ``flip_strands`` (e.g. PRO-seq), as
in the original.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path
from typing import Iterable, Optional

ALLELES = ("A", "C", "G", "T", "del")
#: locus_tag suffixes of rRNA genes without a gene name (as in process_pileup.py)
RRNA = {"mgr01": "RNR1", "mgr02": "RNR2"}

COLUMNS = [
    "sample", "chrom", "position", "ref", "coverage", "pos_coverage", "neg_coverage",
    "read_start", "pos_read_start", "neg_read_start", "read_end", "pos_read_end", "neg_read_end",
    "end_start_counts", "pos_end_start_counts", "neg_end_start_counts",
    "z", "pos_z", "neg_z", "end_start_ratio", "pos_end_start_ratio", "neg_end_start_ratio",
    "A", "C", "G", "T", "N", "del", "ins", "allele_depth", "major", "minor", "minor_count",
    "minor_pos", "minor_neg", "minor_af", "heteroplasmic", "alt_major",
    "feature", "gene", "length", "strand",
]
#: With a NUMT allele table (``mitopiper mt numt-alleles``): the bases NUMTs carry at the position, and
#: whether the minor allele is one of them (a heteroplasmy that NUMT reads can explain).
NUMT_COLUMNS = ["numt_bases", "numt_minor"]
RPM_COLUMNS = ["RPM", "pos_RPM", "neg_RPM", "read_end_RPM", "pos_read_end_RPM", "neg_read_end_RPM",
               "read_start_RPM", "pos_read_start_RPM", "neg_read_start_RPM"]


# ----------------------------------------------------------------- pileup
_START = re.compile(r"\^.", re.S)       # read start + its mapping-quality character (which can be any character)
_INDEL = re.compile(r"[+-](\d+)")
_STRAND_OF = {**dict.fromkeys(".ACGTN*", "fwd"), **dict.fromkeys(",acgtn#", "rev")}
COUNTED = ("A", "C", "G", "T", "N", "del")


def parse_bases(ref: str, bases: str) -> dict:
    """Counts of one pileup base string, per strand (``fwd``/``rev``).

    Returns ``{"fwd": {...}, "rev": {...}}`` with the alleles A/C/G/T/N/del and
    ``ins``, ``start``, ``end`` and ``cov`` (bases and deletions) per strand.
    Base strings of deep RNA-seq data hold 10^5 characters, so the work is done
    with regular expressions and ``str.count``; Python only loops over indels,
    read starts and read ends.
    """
    counts = {s: dict.fromkeys(COUNTED + ("ins", "start", "end", "cov"), 0) for s in ("fwd", "rev")}
    if not bases:
        return counts
    if "^" in bases:                      # first, so that a mapping quality of '+', '$' or '^' is not misread
        bases = _START.sub("^", bases)
    if "+" in bases or "-" in bases:      # drop indel sequences (+2AC); an insertion counts for the base before it
        pieces, done = [], 0
        for match in _INDEL.finditer(bases):
            if match.start() < done:
                continue
            strand = _STRAND_OF.get(bases[match.start() - 1]) if match.start() else None
            if bases[match.start()] == "+" and strand:
                counts[strand]["ins"] += 1
            pieces.append(bases[done:match.start()])
            done = match.end() + int(match.group(1))
        pieces.append(bases[done:])
        bases = "".join(pieces)
    for marker, key in (("^", "start"), ("$", "end")):
        if marker in bases:               # a start marks the base after it, an end the base before it
            pattern = r"\^(.)" if marker == "^" else r"(.)\$"
            for base in re.findall(pattern, bases, re.S):
                if base in _STRAND_OF:
                    counts[_STRAND_OF[base]][key] += 1
    ref_allele = ref.upper() if ref.upper() in "ACGT" else "N"
    for strand, match_char, deletion, case in (("fwd", ".", "*", str.upper), ("rev", ",", "#", str.lower)):
        c = counts[strand]
        for allele in "ACGTN":
            c[allele] = bases.count(case(allele))
        c[ref_allele] += bases.count(match_char)
        c["del"] = bases.count(deletion)
        c["cov"] = sum(c[a] for a in COUNTED)
    return counts


def read_pileup(path) -> Iterable[tuple]:
    """Yield ``(chrom, position, ref, depth, bases)`` of a pileup file (read byte for byte: quality
    strings may hold bytes that are not UTF-8)."""
    with open(path, encoding="latin-1") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or not fields[1].isdigit():
                continue
            depth = int(fields[3]) if fields[3].isdigit() else 0
            yield fields[0], int(fields[1]), fields[2], depth, fields[4] if depth else ""


# ------------------------------------------------------------- annotation
def _attributes(text: str) -> dict:
    """GFF3 (``key=value;``) or GTF (``key "value";``) attributes."""
    attrs = {}
    for part in text.strip().strip(";").split(";"):
        part = part.strip()
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
        else:
            continue
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def _gene_name(attrs: dict) -> Optional[str]:
    for key in ("gene", "gene_name", "Name"):
        if attrs.get(key):
            return attrs[key]
    tag = attrs.get("locus_tag") or attrs.get("gene_id") or ""
    suffix = tag.split("_")[1] if "_" in tag else ""
    return RRNA.get(suffix) or (tag or None)


def read_gff(path) -> list:
    """``(contig, start, end, feature, gene, strand)`` of the gene features of a GFF3/GTF file (1-based, inclusive)."""
    genes = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) < 9 or fields[2] != "gene":
                continue
            name = _gene_name(_attributes(fields[8]))
            if name and fields[3].isdigit() and fields[4].isdigit():
                genes.append((fields[0], int(fields[3]), int(fields[4]), "gene", name, fields[6] if fields[6] in "+-" else ""))
    return genes


def read_genbank(path) -> list:
    """The same for a GenBank flat file (``join``/``complement`` locations are split into parts)."""
    genes, contig, in_features, current = [], "", False, None

    def finish():
        name = _gene_name(current["attrs"]) if current else None
        if name:
            for start, end in current["parts"]:
                genes.append((contig, start, end, "gene", name, current["strand"]))

    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("VERSION") or (line.startswith("LOCUS") and not contig):
                contig = line.split()[1]
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if line.startswith("ORIGIN"):
                break
            if not in_features:
                continue
            key, rest = line[5:21].strip(), line[21:].strip()
            if key:
                finish()
                current = None
                if key == "gene":
                    parts = [(int(a), int(b or a)) for a, b in re.findall(r"<?(\d+)(?:\.\.>?(\d+))?", rest)]
                    current = {"parts": parts, "strand": "-" if "complement" in rest else "+", "attrs": {}}
            elif current is not None and rest.startswith("/"):
                match = re.match(r'/(\w+)="?([^"]*)"?', rest)
                if match:
                    current["attrs"].setdefault(match.group(1), match.group(2))
        finish()
    return genes


def read_annotation(path) -> list:
    with open(path, encoding="utf-8", errors="replace") as handle:
        genbank = handle.read(5) == "LOCUS"
    return read_genbank(path) if genbank else read_gff(path)


def annotation_contigs(genes: list, rows_by_contig: dict) -> dict:
    """Pileup contig -> annotation sequence whose genes it gets.

    Contigs match by name. An annotation of a single sequence that is not in the
    pileup (e.g. ``NC_012920.1`` vs ``chrM``) is applied to the pileup's main
    contig (the one with most positions), with a warning.
    """
    seqids = {g[0] for g in genes}
    mapping = {c: c for c in rows_by_contig if c in seqids}
    if genes and not mapping and len(seqids) == 1 and rows_by_contig:
        main = max(rows_by_contig, key=lambda c: len(rows_by_contig[c]))
        (seqid,) = seqids
        print(f"[mitopiper] annotation sequence {seqid} is not in the pileup: using it for {main}", file=sys.stderr)
        mapping[main] = seqid
    elif genes and not mapping:
        print(f"[mitopiper] warning: no pileup contig matches the annotation ({', '.join(sorted(seqids))}): "
              "genes are left empty", file=sys.stderr)
    return mapping


def genes_by_position(genes: list, seqid: Optional[str]) -> dict:
    """Position -> list of ``(feature, gene, length, strand)`` of the genes on ``seqid``."""
    table = {}
    for contig, start, end, feature, name, strand in genes:
        if contig == seqid:
            for position in range(start, end + 1):
                table.setdefault(position, []).append((feature, name, end - start + 1, strand))
    return table


# ---------------------------------------------------------------- counts
def htseq_total(path, include_unaligned: bool = False, extra: int = 0) -> int:
    """Total reads for RPM: the gene counts plus ``__no_feature`` of an htseq-count table
    (plus ``__not_aligned`` with ``include_unaligned``), plus ``extra``."""
    total = extra
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[1].strip().isdigit():
                continue
            name, count = fields[0], int(fields[1])
            if not name.startswith("__") or name == "__no_feature" or (include_unaligned and name == "__not_aligned"):
                total += count
    return total


def _stats(values: list) -> tuple:
    if not values:
        return 0.0, 0.0
    mean = sum(values) / len(values)
    return mean, math.sqrt(sum((v - mean) ** 2 for v in values) / len(values))


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


# ---------------------------------------------------------------- annotate
def annotate_pileup(pileup, output, annotation=None, sample: str = "", flip_strands: bool = False,
                    total_reads: Optional[int] = None, min_depth: int = 50, min_af: float = 0.03,
                    min_strand: int = 2, numt_alleles=None) -> dict:
    """Write the per-position table of ``pileup`` to ``output`` (TSV). Returns counts of positions,
    heteroplasmic positions (and of those whose minor allele a NUMT carries, with ``numt_alleles``)
    and positions whose major allele is not the reference."""
    sample = sample or Path(pileup).name.split(".")[0]
    rows_by_contig = {}
    for chrom, position, ref, depth, bases in read_pileup(pileup):
        rows_by_contig.setdefault(chrom, {})[position] = (ref, depth, parse_bases(ref, bases))
    genes = read_annotation(annotation) if annotation else []
    seqid_of = annotation_contigs(genes, rows_by_contig)
    numt = read_numt_alleles(numt_alleles) if numt_alleles else None
    columns = COLUMNS + (NUMT_COLUMNS if numt is not None else []) + (RPM_COLUMNS if total_reads else [])
    pos_s, neg_s = ("rev", "fwd") if flip_strands else ("fwd", "rev")
    thresholds = (min_depth, min_af, min_strand)
    summary = {"positions": 0, "heteroplasmic": 0, "alt_major": 0, **({"numt_minor": 0} if numt is not None else {})}
    empty = ("N", 0, parse_bases("N", ""))
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for chrom, rows in rows_by_contig.items():
            # positions missing from the pileup have no coverage
            table = [(p, *rows.get(p, empty)) for p in range(1, max(rows) + 1)]
            strands_of = [_strand_counts(depth, counts[pos_s], counts[neg_s]) for _, _, depth, counts in table]
            stats = {prefix: _stats([s[prefix]["coverage"] for s in strands_of]) for prefix in STRAND_PREFIXES}
            gene_table = genes_by_position(genes, seqid_of.get(chrom))
            for (position, ref, _, counts), strands in zip(table, strands_of):
                row = {"sample": sample, "chrom": chrom, "position": position, "ref": ref.upper()}
                for prefix, c in strands.items():
                    mean, std = stats[prefix]
                    row.update({f"{prefix}{k}": v for k, v in c.items()})
                    row[f"{prefix}z"] = (c["coverage"] - mean) / std if std else None
                    # as in process_pileup.py: only where coverage is above mean - 1 sd
                    row[f"{prefix}end_start_ratio"] = (c["end_start_counts"] / c["coverage"]
                                                       if c["coverage"] and c["coverage"] > mean - std else None)
                    if total_reads:
                        for k in ("coverage", "read_end", "read_start"):
                            name = "RPM" if k == "coverage" else f"{k}_RPM"
                            row[f"{prefix}{name}"] = c[k] * 1e6 / total_reads
                row.update(_allele_fields(ref.upper(), counts[pos_s], counts[neg_s], thresholds))
                if numt is not None:
                    bases = numt.get((chrom, position), "")
                    row["numt_bases"] = bases
                    row["numt_minor"] = "yes" if row["minor"] in tuple(bases) and row["minor"] != row["major"] else "no"
                    summary["numt_minor"] += row["heteroplasmic"] == "yes" and row["numt_minor"] == "yes"
                summary["heteroplasmic"] += row["heteroplasmic"] == "yes"
                summary["alt_major"] += row["alt_major"] == "yes"
                hits = gene_table.get(position, [])
                for i, col in enumerate(("feature", "gene", "length", "strand")):
                    row[col] = ";".join(str(h[i]) for h in hits)
                writer.writerow([_fmt(row.get(col)) for col in columns])
                summary["positions"] += 1
    return summary


STRAND_PREFIXES = ("", "pos_", "neg_")


def read_numt_alleles(path) -> dict:
    """``{(chrom, position): bases}`` from a ``mitopiper mt numt-alleles`` table."""
    with open(path, encoding="utf-8") as handle:
        return {(r["chrom"], int(r["position"])): r["numt_bases"] for r in csv.DictReader(handle, delimiter="\t")}


def _strand_counts(depth: int, p: dict, n: dict) -> dict:
    """Coverage, read starts/ends and their sum, in total (the pileup depth) and per strand."""
    result = {}
    for prefix, c in (("pos_", p), ("neg_", n)):
        result[prefix] = {"coverage": c["cov"], "read_start": c["start"], "read_end": c["end"],
                          "end_start_counts": c["start"] + c["end"]}
    result[""] = {k: result["pos_"][k] + result["neg_"][k] for k in result["pos_"]}
    result[""]["coverage"] = depth
    return {prefix: result[prefix] for prefix in STRAND_PREFIXES}


def _allele_fields(ref: str, p: dict, n: dict, thresholds: tuple) -> dict:
    """Allele counts, major/minor allele, minor allele frequency and the heteroplasmy call."""
    min_depth, min_af, min_strand = thresholds
    alleles = {a: p[a] + n[a] for a in COUNTED}
    fields = dict(alleles, ins=p["ins"] + n["ins"])
    depth = sum(alleles[a] for a in ALLELES)
    fields["allele_depth"] = depth
    if not depth:
        return dict(fields, major="", minor="", minor_count=0, minor_pos=0, minor_neg=0, minor_af=None,
                    heteroplasmic="no", alt_major="no")
    major, minor = sorted(ALLELES, key=lambda a: -alleles[a])[:2]
    af = alleles[minor] / depth
    het = depth >= min_depth and af >= min_af and p[minor] >= min_strand and n[minor] >= min_strand
    # A variant needs the same depth as a heteroplasmy and a true majority: without them one read with an
    # error at a barely covered position counted as a homoplasmic variant (seen in RNA-seq on BGU data).
    alt = ref in "ACGT" and major != ref and depth >= min_depth and alleles[major] > depth / 2
    return dict(fields, major=major, minor=minor if alleles[minor] else "", minor_count=alleles[minor],
                minor_pos=p[minor], minor_neg=n[minor], minor_af=af,
                heteroplasmic="yes" if het else "no", alt_major="yes" if alt else "no")


# ------------------------------------------------------------ heteroplasmy
SITE_COLUMNS = ["sample", "chrom", "position", "ref", "gene", "type", "major", "minor", "minor_af",
                "allele_depth", "minor_pos", "minor_neg", "numt_minor"]


def _read_tables(tables: list):
    for table in tables:
        with open(table, encoding="utf-8") as handle:
            yield from csv.DictReader(handle, delimiter="\t")


def collect_heteroplasmy(tables: list, outdir) -> dict:
    """Gather the heteroplasmic and alternative-major positions of annotated tables into
    ``outdir/heteroplasmy_sites.tsv``, and the minor allele frequency of every such position in
    every sample into ``outdir/heteroplasmy_matrix.tsv``. Returns the number of sites per sample.

    The tables are read twice (calls, then the frequencies at the called positions), so memory
    does not grow with positions x samples."""
    sites, per_sample, called = [], {}, {}
    for row in _read_tables(tables):
        per_sample.setdefault(row["sample"], 0)
        for kind, col in (("heteroplasmy", "heteroplasmic"), ("variant", "alt_major")):
            if row[col] == "yes":
                sites.append({**{c: row.get(c, "") for c in SITE_COLUMNS}, "type": kind})
                per_sample[row["sample"]] += 1
                called.setdefault((row["chrom"], int(row["position"])), (row["ref"], row["gene"]))
    af = {}
    for row in _read_tables(tables):
        key = (row["chrom"], int(row["position"]))
        if key in called:
            af[key, row["sample"]] = row["minor_af"]
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "heteroplasmy_sites.tsv", "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, SITE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sites)
    samples = list(per_sample)
    with open(outdir / "heteroplasmy_matrix.tsv", "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "position", "ref", "gene"] + samples)
        for key in sorted(called):
            writer.writerow([*key, *called[key]] + [af.get((key, s), "") for s in samples])
    return per_sample


# ---------------------------------------------------------------------- CLI
def add_parsers(msub) -> None:
    """``mitopiper mt annotate`` and ``mitopiper mt heteroplasmy``."""
    ann = msub.add_parser("annotate", help="per-position coverage, read starts/ends, alleles and genes of a pileup")
    ann.add_argument("pileup")
    ann.add_argument("output", help="output table (TSV)")
    ann.add_argument("--annotation", default="", help="GFF3, GTF or GenBank file with the genes")
    ann.add_argument("--sample", default="", help="sample name (default: from the pileup file name)")
    ann.add_argument("--counts", default="", help="htseq-count table, for RPM values")
    ann.add_argument("--include-unaligned", action="store_true", help="count __not_aligned reads in the RPM total")
    ann.add_argument("--extra-reads", type=int, default=0, help="reads added to the RPM total (e.g. mapped elsewhere)")
    ann.add_argument("--total-reads", type=int, default=0, help="RPM total, instead of --counts")
    ann.add_argument("--flip-strands", action="store_true", help="swap the strands (e.g. PRO-seq)")
    ann.add_argument("--min-depth", type=int, default=50, help="heteroplasmy: minimum depth of A/C/G/T/deletions")
    ann.add_argument("--min-af", type=float, default=0.03, help="heteroplasmy: minimum minor allele frequency")
    ann.add_argument("--min-strand", type=int, default=2, help="heteroplasmy: minimum minor allele reads on each strand")
    ann.add_argument("--numt-alleles", default="", help="table of mitopiper mt numt-alleles: flag positions whose minor "
                                                         "allele a NUMT carries (numt_minor)")
    het = msub.add_parser("heteroplasmy", help="collect heteroplasmic sites of annotated tables of many samples")
    het.add_argument("tables", nargs="+")
    het.add_argument("-o", "--outdir", required=True)
    for parser in (ann, het):
        parser.set_defaults(mt_run=run)


def run(args) -> int:
    if args.action == "annotate":
        total = args.total_reads or (htseq_total(args.counts, args.include_unaligned, args.extra_reads) if args.counts else None)
        summary = annotate_pileup(args.pileup, args.output, args.annotation or None, args.sample, args.flip_strands,
                                  total, args.min_depth, args.min_af, args.min_strand, args.numt_alleles or None)
        numt = f" ({summary['numt_minor']} with a minor allele that a NUMT carries)" if args.numt_alleles else ""
        print(f"[mitopiper] {summary['positions']} positions, {summary['heteroplasmic']} heteroplasmic{numt}, "
              f"{summary['alt_major']} with a non-reference major allele", file=sys.stderr)
    else:
        tables = [t for t in args.tables if Path(t).exists()]
        for missing in sorted(set(args.tables) - set(tables)):
            print(f"[mitopiper] warning: {missing} does not exist (failed sample?) - left out", file=sys.stderr)
        if not tables:
            print("[mitopiper] no annotated tables to collect", file=sys.stderr)
            return 1
        per_sample = collect_heteroplasmy(tables, args.outdir)
        print("[mitopiper] sites per sample: " + ", ".join(f"{s} {n}" for s, n in per_sample.items()), file=sys.stderr)
    return 0
