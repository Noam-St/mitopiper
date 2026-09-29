#!/usr/bin/env python3
"""Simulate a small mtDNA dataset with a known answer, for the mitopiper tutorial.

It writes, into OUTDIR:

* ``reference/mt.fa`` and ``reference/mt.gff3``: a synthetic circular mtDNA with a
  few annotated genes (on both strands, two of them overlapping).
* ``reference/nuclear.fa``: a synthetic nuclear genome. It contains a NUMT (a 3%
  diverged copy of part of the mtDNA) and a ``chrM`` copy of the mtDNA, as many
  real assemblies do. mitopiper leaves the latter out of the nuclear + mtDNA
  reference.
* ``fastq/<sample>_1.fastq.gz`` / ``_2``: paired-end reads of each sample. Fragments
  are drawn from the *circular* mtDNA, so some span the origin. A fraction of
  them come from the NUMT, and every base has a small sequencing error rate.
* ``truth.tsv``: every true heteroplasmy and variant of every sample, and every
  position where the NUMT differs from the mtDNA (where NUMT reads fake
  heteroplasmy).

Everything comes from ``--seed``, so the dataset is the same on every machine.
It needs Python 3 only, and takes a few seconds::

    python tutorial/simulate.py tutorial_data
"""

import argparse
import gzip
import io
import random
from pathlib import Path

COMPLEMENT = str.maketrans("ACGT", "TGCA")

#: sample -> {position: allele frequency of a transition there}; frequency 1.0 = homoplasmic variant.
#: Position 60 is next to the origin, 3000 inside the NUMT region, 12000 in a gene on the minus strand.
SAMPLES = {
    "patient": {60: 0.15, 3000: 0.20, 9000: 0.10, 12000: 0.05, 7000: 1.0},
    "control": {7000: 1.0},
}
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}
GENES = [("RNR1", 700, 1600, "+"), ("ND1", 2800, 3800, "+"), ("ATP8", 8400, 8600, "+"), ("ATP6", 8560, 9240, "+"),
         ("ND6", 11500, 12000, "-"), ("CYTB", 13000, 14100, "+")]


def rc(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def random_seq(rng: random.Random, length: int, gc: float = 0.44) -> str:
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights, k=length))


def mutate(rng: random.Random, seq: str, rate: float) -> tuple:
    """``seq`` with a fraction ``rate`` of its bases substituted; also the positions (0-based) that differ."""
    bases, changed = list(seq), []
    for i in rng.sample(range(len(seq)), int(len(seq) * rate)):
        bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
        changed.append(i)
    return "".join(bases), sorted(changed)


def gzip_text(path: Path):
    """A gzip text file without a timestamp in its header, so every run gives identical files."""
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", mtime=0), encoding="ascii")


def write_fasta(path: Path, records: list) -> None:
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            out.writelines(seq[i:i + 70] + "\n" for i in range(0, len(seq), 70))


def reads(rng, template: str, circular: bool, pairs: int, length: int, insert: int, error: float):
    """Paired reads (R1 forward, R2 reverse complement) of fragments of ``template``."""
    for _ in range(pairs):
        size = max(length, int(rng.gauss(insert, insert * 0.1)))
        if circular:
            start = rng.randrange(len(template))
            fragment = (template + template[:size])[start:start + size]
        else:
            start = rng.randrange(len(template) - size)
            fragment = template[start:start + size]
        if rng.random() < 0.5:
            fragment = rc(fragment)
        yield tuple("".join(b if rng.random() > error else rng.choice("ACGT".replace(b, "")) for b in r)
                    for r in (fragment[:length], rc(fragment)[:length]))


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir")
    parser.add_argument("--seed", type=int, default=2024)
    parser.add_argument("--mt-length", type=int, default=16000)
    parser.add_argument("--depth", type=int, default=300, help="mtDNA read depth")
    parser.add_argument("--numt-fraction", type=float, default=0.06, help="fraction of read pairs from the NUMT")
    parser.add_argument("--read-length", type=int, default=100)
    parser.add_argument("--error", type=float, default=0.002, help="sequencing error rate per base")
    args = parser.parse_args(argv)

    rng = random.Random(args.seed)
    out = Path(args.outdir)
    (out / "reference").mkdir(parents=True, exist_ok=True)
    (out / "fastq").mkdir(parents=True, exist_ok=True)

    mt = random_seq(rng, args.mt_length)
    numt_start, numt_end = 2000, 5000                           # the NUMT copies mtDNA 2001..5000
    numt, numt_diffs = mutate(rng, mt[numt_start:numt_end], 0.03)
    chromosomes = [random_seq(rng, 150_000, gc=0.40) for _ in range(3)]
    chromosomes[1] = chromosomes[1][:60_000] + numt + chromosomes[1][60_000:]
    write_fasta(out / "reference" / "mt.fa", [("synthMT", mt)])
    write_fasta(out / "reference" / "nuclear.fa",
                [(f"chr{i + 1}", c) for i, c in enumerate(chromosomes)] + [("chrM", mt)])
    with open(out / "reference" / "mt.gff3", "w") as gff:
        gff.write("##gff-version 3\n")
        gff.writelines(f"synthMT\tsimulated\tgene\t{s}\t{e}\t.\t{strand}\t.\tID=gene-{g};Name={g};gene={g}\n"
                       for g, s, e, strand in GENES)

    truth = [("sample", "position", "ref", "alt", "af", "type")]
    pairs = int(args.depth * args.mt_length / (2 * args.read_length))
    numt_pairs = int(pairs * args.numt_fraction)
    for sample, changes in SAMPLES.items():
        # one haplotype per allele frequency: reads are drawn from each in proportion
        haplotypes = [(1.0, mt)]
        for position, af in sorted(changes.items()):
            truth.append((sample, position, mt[position - 1], TRANSITION[mt[position - 1]], af,
                          "variant" if af == 1.0 else "heteroplasmy"))
        for position, af in changes.items():
            alt = TRANSITION[mt[position - 1]]
            new = []
            for weight, seq in haplotypes:
                mutated = seq[:position - 1] + alt + seq[position:]
                new += [(weight * af, mutated), (weight * (1 - af), seq)] if af < 1 else [(weight, mutated)]
            haplotypes = [(w, s) for w, s in new if w > 0]
        for i in numt_diffs:
            truth.append((sample, numt_start + i + 1, mt[numt_start + i], numt[i], "", "numt"))
        with gzip_text(out / "fastq" / f"{sample}_1.fastq.gz") as r1, gzip_text(out / "fastq" / f"{sample}_2.fastq.gz") as r2:
            n, quality = 0, "I" * args.read_length
            sources = [(seq, True, round(pairs * w)) for w, seq in haplotypes] + [(numt, False, numt_pairs)]
            for template, circular, count in sources:
                for a, b in reads(rng, template, circular, count, args.read_length, 300, args.error):
                    n += 1
                    r1.write(f"@{sample}.{n}/1\n{a}\n+\n{quality}\n")
                    r2.write(f"@{sample}.{n}/2\n{b}\n+\n{quality}\n")
    with open(out / "truth.tsv", "w") as handle:
        handle.writelines("\t".join(map(str, row)) + "\n" for row in truth)
    print(f"Wrote {out}: mtDNA of {args.mt_length} bp, {pairs + numt_pairs} read pairs per sample "
          f"({numt_pairs} from a NUMT with {len(numt_diffs)} differences), samples {', '.join(SAMPLES)}")


if __name__ == "__main__":
    main()
