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
RPM_COLUMNS = ["RPM", "pos_RPM", "neg_RPM", "read_end_RPM", "pos_read_end_RPM", "neg_read_end_RPM",
               "read_start_RPM", "pos_read_start_RPM", "neg_read_start_RPM"]


# ----------------------------------------------------------------- pileup
def parse_bases(ref: str, bases: str) -> dict:
    """Counts of one pileup base string, per strand (``fwd``/``rev``).

    Returns ``{"fwd": {...}, "rev": {...}}`` with the alleles A/C/G/T/N/del and
    ``ins``, ``start``, ``end`` and ``cov`` (bases and deletions) per strand.
    """
    ref = ref.upper()
    counts = {s: dict.fromkeys(("A", "C", "G", "T", "N", "del", "ins", "start", "end", "cov"), 0) for s in ("fwd", "rev")}
    i, last, starting = 0, None, False
    while i < len(bases):
        c = bases[i]
        if c == "^":                       # read start; the next character is its mapping quality
            starting, i = True, i + 2
            continue
        if c == "$":                       # the read of the previous base ends here
            if last:
                counts[last]["end"] += 1
            i += 1
            continue
        if c in "+-":
            match = re.match(r"\d+", bases[i + 1:])
            length = int(match.group()) if match else 0
            if c == "+" and last:
                counts[last]["ins"] += 1
            i += 1 + (len(match.group()) if match else 0) + length
            continue
        strand = None
        if c in ".,":
            strand, allele = ("fwd" if c == "." else "rev"), ref if ref in "ACGT" else "N"
        elif c.upper() in "ACGTN":
            strand, allele = ("fwd" if c.isupper() else "rev"), c.upper()
        elif c in "*#":
            strand, allele = ("fwd" if c == "*" else "rev"), "del"
        if strand:
            counts[strand][allele] += 1
            counts[strand]["cov"] += 1
            if starting:
                counts[strand]["start"] += 1
            last = strand
        else:                              # '>'/'<' reference skips, or anything unknown
            last = None
        starting = False
        i += 1
    return counts


def read_pileup(path) -> Iterable[tuple]:
    """Yield ``(chrom, position, ref, depth, bases)`` of a pileup file."""
    with open(path, encoding="utf-8") as handle:
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
        if current and current["name"]:
            for start, end in current["parts"]:
                genes.append((contig, start, end, "gene", current["name"], current["strand"]))

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
                    current = {"parts": parts, "strand": "-" if "complement" in rest else "+", "name": None, "tag": None}
            elif current is not None and rest.startswith("/"):
                match = re.match(r'/(gene|locus_tag)="?([^"]*)"?', rest)
                if match and match.group(1) == "gene":
                    current["name"] = match.group(2)
                elif match and not current["name"]:
                    tag = match.group(2)
                    current["name"] = RRNA.get(tag.split("_")[1] if "_" in tag else "", tag)
        finish()
    return genes


def read_annotation(path) -> list:
    text = Path(path).read_text(encoding="utf-8", errors="replace")[:200]
    return read_genbank(path) if text.startswith("LOCUS") else read_gff(path)


def genes_by_position(genes: list, contig: str, contigs_in_pileup: set) -> dict:
    """Position -> list of ``(feature, gene, length, strand)`` on ``contig``.

    When the annotation names a single sequence that is not in the pileup (e.g.
    ``NC_012920.1`` vs ``chrM``), it is applied to the pileup's main contig.
    """
    seqids = {g[0] for g in genes}
    use = [g for g in genes if g[0] == contig]
    if not use and len(seqids) == 1 and not seqids & contigs_in_pileup:
        use = genes
    table = {}
    for _, start, end, feature, name, strand in use:
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
                    min_strand: int = 2) -> dict:
    """Write the per-position table of ``pileup`` to ``output`` (TSV). Returns counts of positions,
    heteroplasmic positions and positions whose major allele is not the reference."""
    sample = sample or Path(pileup).name.split(".")[0]
    rows_by_contig = {}
    for chrom, position, ref, depth, bases in read_pileup(pileup):
        rows_by_contig.setdefault(chrom, {})[position] = (ref, depth, parse_bases(ref, bases))
    genes = read_annotation(annotation) if annotation else []
    columns = COLUMNS + (RPM_COLUMNS if total_reads else [])
    pos_s, neg_s = ("rev", "fwd") if flip_strands else ("fwd", "rev")
    summary = {"positions": 0, "heteroplasmic": 0, "alt_major": 0}
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for chrom, rows in rows_by_contig.items():
            empty = ("N", 0, parse_bases("N", ""))
            positions = range(1, max(rows) + 1)   # positions missing from the pileup have no coverage
            table = [(p, *rows.get(p, empty)) for p in positions]
            cov = [t[2] for t in table]
            pos_cov = [t[3][pos_s]["cov"] for t in table]
            neg_cov = [t[3][neg_s]["cov"] for t in table]
            stats = {"": _stats(cov), "pos_": _stats(pos_cov), "neg_": _stats(neg_cov)}
            gene_table = genes_by_position(genes, chrom, set(rows_by_contig))
            for (position, ref, depth, counts), c, pc, nc in zip(table, cov, pos_cov, neg_cov):
                p, n = counts[pos_s], counts[neg_s]
                row = {"sample": sample, "chrom": chrom, "position": position, "ref": ref.upper(),
                       "coverage": c, "pos_coverage": pc, "neg_coverage": nc,
                       "pos_read_start": p["start"], "neg_read_start": n["start"],
                       "pos_read_end": p["end"], "neg_read_end": n["end"]}
                row["read_start"] = p["start"] + n["start"]
                row["read_end"] = p["end"] + n["end"]
                for prefix, strand in (("pos_", p), ("neg_", n)):
                    row[f"{prefix}end_start_counts"] = strand["start"] + strand["end"]
                row["end_start_counts"] = row["pos_end_start_counts"] + row["neg_end_start_counts"]
                for prefix, value in (("", c), ("pos_", pc), ("neg_", nc)):
                    mean, std = stats[prefix]
                    row[f"{prefix}z"] = (value - mean) / std if std else None
                    # as in process_pileup.py: only where coverage is above mean - 1 sd
                    row[f"{prefix}end_start_ratio"] = (row[f"{prefix}end_start_counts"] / value
                                                       if value and value > mean - std else None)
                alleles = {a: p[a] + n[a] for a in ("A", "C", "G", "T", "N", "del")}
                row.update(alleles, ins=p["ins"] + n["ins"])
                ranked = sorted(ALLELES, key=lambda a: -alleles[a])
                allele_depth = sum(alleles[a] for a in ALLELES)
                row["allele_depth"] = allele_depth
                if allele_depth:
                    major, minor = ranked[0], ranked[1]
                    minor_count = alleles[minor]
                    row.update(major=major, minor=minor if minor_count else "", minor_count=minor_count,
                               minor_pos=p[minor], minor_neg=n[minor], minor_af=minor_count / allele_depth)
                    het = (allele_depth >= min_depth and row["minor_af"] >= min_af
                           and p[minor] >= min_strand and n[minor] >= min_strand)
                    alt = ref.upper() in "ACGT" and major != ref.upper()
                    row["heteroplasmic"], row["alt_major"] = ("yes" if het else "no"), ("yes" if alt else "no")
                    summary["heteroplasmic"] += het
                    summary["alt_major"] += alt
                else:
                    row.update(major="", minor="", minor_count=0, minor_pos=0, minor_neg=0, minor_af=None,
                               heteroplasmic="no", alt_major="no")
                hits = gene_table.get(position, [])
                row["feature"] = ";".join(h[0] for h in hits)
                row["gene"] = ";".join(h[1] for h in hits)
                row["length"] = ";".join(str(h[2]) for h in hits)
                row["strand"] = ";".join(h[3] for h in hits)
                if total_reads:
                    for name in ("", "pos_", "neg_"):
                        row[f"{name}RPM"] = row[f"{name}coverage"] * 1e6 / total_reads
                        row[f"{name}read_end_RPM"] = row[f"{name}read_end"] * 1e6 / total_reads
                        row[f"{name}read_start_RPM"] = row[f"{name}read_start"] * 1e6 / total_reads
                writer.writerow([_fmt(row.get(col)) for col in columns])
                summary["positions"] += 1
    return summary


# ------------------------------------------------------------ heteroplasmy
SITE_COLUMNS = ["sample", "chrom", "position", "ref", "gene", "type", "major", "minor", "minor_af",
                "allele_depth", "minor_pos", "minor_neg"]


def collect_heteroplasmy(tables: list, outdir) -> dict:
    """Gather the heteroplasmic and alternative-major positions of annotated tables into
    ``outdir/heteroplasmy_sites.tsv``, and the minor allele frequency of every such position in
    every sample into ``outdir/heteroplasmy_matrix.tsv``. Returns the number of sites per sample."""
    sites, af, samples, per_sample = [], {}, [], {}
    for table in tables:
        with open(table, encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                sample = row["sample"]
                if sample not in per_sample:
                    samples.append(sample)
                    per_sample[sample] = 0
                key = (row["chrom"], int(row["position"]))
                af.setdefault(key, {})[sample] = row["minor_af"]
                af[key].setdefault("_ref", row["ref"])
                af[key].setdefault("_gene", row["gene"])
                kinds = [k for k, col in (("heteroplasmy", "heteroplasmic"), ("variant", "alt_major")) if row[col] == "yes"]
                for kind in kinds:
                    sites.append({**{c: row.get(c, "") for c in SITE_COLUMNS}, "type": kind})
                    per_sample[sample] += 1
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "heteroplasmy_sites.tsv", "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, SITE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sites)
    called = sorted({(s["chrom"], int(s["position"])) for s in sites})
    with open(outdir / "heteroplasmy_matrix.tsv", "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "position", "ref", "gene"] + samples)
        for key in called:
            writer.writerow([key[0], key[1], af[key]["_ref"], af[key]["_gene"]] + [af[key].get(s, "") for s in samples])
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
    het = msub.add_parser("heteroplasmy", help="collect heteroplasmic sites of annotated tables of many samples")
    het.add_argument("tables", nargs="+")
    het.add_argument("-o", "--outdir", required=True)


def run(args) -> int:
    if args.action == "annotate":
        total = args.total_reads or (htseq_total(args.counts, args.include_unaligned, args.extra_reads) if args.counts else None)
        summary = annotate_pileup(args.pileup, args.output, args.annotation or None, args.sample, args.flip_strands,
                                  total, args.min_depth, args.min_af, args.min_strand)
        print(f"[mitopiper] {summary['positions']} positions, {summary['heteroplasmic']} heteroplasmic, "
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
