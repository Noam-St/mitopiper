import io

import pytest

from mitopiper import mtdna
from mitopiper.cli import main


def fasta(path, header, seq):
    path.write_text(f">{header}\n" + "\n".join(seq[i:i + 7] for i in range(0, len(seq), 7)) + "\n")
    return path


def pileup_line(contig, pos, depth, tag):
    return f"{contig}\t{pos}\tN\t{depth}\t{tag}\tI\n"


def test_rotate_like_circulize(tmp_path):
    src = fasta(tmp_path / "mt.fa", "NC_1.1 mitochondrion", "ABCDEFGHIJ")
    assert mtdna.rotate_fasta(src, tmp_path / "rot.fa", shift=3) == 10
    ((header, seq),) = list(mtdna.read_fasta(open(tmp_path / "rot.fa")))
    assert header == "NC_1.1 mitochondrion" and seq == "DEFGHIJABC"     # first 3 bases moved to the end
    for r in range(1, 11):                                            # every rotated base maps back to itself
        assert "ABCDEFGHIJ"[mtdna.rotated_to_normal(r, 10, 3) - 1] == seq[r - 1]
    with pytest.raises(ValueError):
        mtdna.rotate_fasta(src, tmp_path / "x.fa", shift=10)


def test_combine_pileups_takes_the_origin_from_the_rotated_alignment(tmp_path):
    length, shift, window = 20, 6, 4
    normal = tmp_path / "n.pileup"
    rotated = tmp_path / "r.pileup"
    normal.write_text("".join(pileup_line("MT", p, 10, f"n{p}") for p in range(1, length + 1)))
    # the rotated pileup, in rotated coordinates; its tag is the normal position it stands for
    rotated.write_text("".join(pileup_line("MT", r, 30, f"r{mtdna.rotated_to_normal(r, length, shift)}")
                               for r in range(1, length + 1)))
    counts = mtdna.combine_pileups(normal, rotated, tmp_path / "c.pileup", length, shift, window)
    rows = [line.split("\t") for line in (tmp_path / "c.pileup").read_text().splitlines()]
    assert [int(r[1]) for r in rows] == list(range(1, length + 1))
    source = {int(r[1]): r[4] for r in rows}
    for p in range(1, length + 1):
        near_origin = p <= window or p > length - window
        assert source[p] == (f"r{p}" if near_origin else f"n{p}")
    assert counts == {"normal": 12, "rotated": 8, "other": 0}


def test_combine_pileups_start_only_matches_the_original_scripts(tmp_path):
    # combine_processed_pileups.py replaced only the first `to_return` positions
    normal = tmp_path / "n.pileup"
    rotated = tmp_path / "r.pileup"
    normal.write_text("".join(pileup_line("MT", p, 1, f"n{p}") for p in range(1, 101)))
    rotated.write_text("".join(pileup_line("MT", r, 1, f"r{r}") for r in range(1, 101)))
    mtdna.combine_pileups(normal, rotated, tmp_path / "c.pileup", 100, shift=50, window=40, tail_window=0)
    rows = {int(line.split("\t")[1]): line.split("\t")[4] for line in (tmp_path / "c.pileup").read_text().splitlines()}
    assert rows[1] == "r51" and rows[40] == "r90" and rows[41] == "n41" and rows[100] == "n100"
    with pytest.raises(ValueError, match="must not exceed"):
        mtdna.combine_pileups(normal, rotated, tmp_path / "x", 100, shift=50, window=60)


def test_default_rotation_is_half_the_length(tmp_path):
    src = fasta(tmp_path / "mt.fa", "MT", "ABCDEFGHIJ")
    assert main(["mt", "rotate", str(src), str(tmp_path / "rot.fa")]) == 0
    with open(tmp_path / "rot.fa") as handle:
        assert list(mtdna.read_fasta(handle)) == [("MT", "FGHIJABCDE")]
    # positions 1-2 and 9-10 come from the middle of the rotated copy (6-7 and 4-5), not from its ends
    (tmp_path / "n.p").write_text("".join(pileup_line("MT", p, 1, "n") for p in range(1, 11)))
    (tmp_path / "r.p").write_text("".join(pileup_line("MT", p, 1, f"r{p}") for p in range(1, 11)))
    mtdna.combine_pileups(tmp_path / "n.p", tmp_path / "r.p", tmp_path / "c.p", 10, window=2)
    tags = [line.split("\t")[4] for line in (tmp_path / "c.p").read_text().splitlines()]
    assert tags == ["r6", "r7"] + ["n"] * 6 + ["r4", "r5"]
    assert mtdna.resolve_shift("half", 16569) == 8284 and mtdna.resolve_shift("500", 16569) == 500


def test_combine_references_drops_mitochondrial_copies(tmp_path):
    nuclear = tmp_path / "nuc.fa"
    nuclear.write_text(">chr1 Homo sapiens chromosome 1\nAAAA\nAAAA\n>chrM\nCCCC\n>NC_012920.1 Homo sapiens mitochondrion\nGGGG\n"
                       ">chr2\nTTTT\n")
    mt = fasta(tmp_path / "mt.fa", "NC_012920.1 Homo sapiens mitochondrion, complete genome", "ACGTACGTAC")
    log = io.StringIO()
    dropped = mtdna.combine_references(nuclear, mt, tmp_path / "combined.fa", log)
    assert dropped == ["chrM", "NC_012920.1"] and "chrM" in log.getvalue()
    assert mtdna.contigs(tmp_path / "combined.fa") == ["chr1", "chr2", "NC_012920.1"]


def test_mt_commands(tmp_path, capsys):
    mt = fasta(tmp_path / "mt.fa", "MT", "ACGTACGTACGTACGTACGT")
    assert main(["mt", "rotate", str(mt), str(tmp_path / "rot.fa"), "--shift", "5"]) == 0
    assert mtdna.sequence_length(tmp_path / "rot.fa") == 20
    assert main(["mt", "contigs", str(mt)]) == 0
    assert capsys.readouterr().out.strip() == "MT"
    (tmp_path / "n.p").write_text("".join(pileup_line("MT", p, 1, "n") for p in range(1, 21)))
    (tmp_path / "r.p").write_text("".join(pileup_line("MT", p, 1, "r") for p in range(1, 21)))
    assert main(["mt", "combine-pileups", str(tmp_path / "n.p"), str(tmp_path / "r.p"), str(tmp_path / "c.p"),
                 "--reference", str(mt), "--shift", "5", "--window", "3"]) == 0
    assert len((tmp_path / "c.p").read_text().splitlines()) == 20


def test_fragments_wrap_around_the_origin(tmp_path):
    (tmp_path / "mt.fa").write_text(">MT x\nAACCGGTT\n")
    assert mtdna.fragments(tmp_path / "mt.fa", tmp_path / "f.fa", size=20, step=4) == 2
    with open(tmp_path / "f.fa") as handle:
        records = list(mtdna.read_fasta(handle))
    assert records == [("MT_1", ("AACCGGTT" * 3)[:20]), ("MT_5", ("GGTTAACCGGTTAACCGGTTAA")[:20])]


def test_numt_regions_merge_flank_and_ignore_mt_hits():
    sam = io.StringIO("@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chrM\tLN:100\n"
                      "f1\t0\tchr1\t200\t0\t60M\t*\t0\t0\t*\t*\n"
                      "f2\t256\tchr1\t300\t0\t10S50M\t*\t0\t0\t*\t*\n"   # within merge distance
                      "f3\t0\tchr1\t950\t0\t60M\t*\t0\t0\t*\t*\n"        # clipped to the contig end
                      "f4\t0\tchr1\t600\t0\t30M\t*\t0\t0\t*\t*\n"        # too short
                      "f5\t0\tchrM\t1\t0\t60M\t*\t0\t0\t*\t*\n"
                      "f6\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n")
    assert mtdna.numt_regions(sam, set()) == ["chr1:100-449", "chr1:850-1000"]


def test_numt_reference_without_hits_explains_what_to_do(tmp_path):
    (tmp_path / "mt.fa").write_text(">MT\nACGT\n")
    (tmp_path / "hits.sam").write_text("@HD\n")
    with pytest.raises(ValueError, match="mt-variants preset"):
        mtdna.numt_reference(tmp_path / "hits.sam", tmp_path / "mt.fa", tmp_path / "r.txt")


@pytest.mark.parametrize("keep, pairs, expected", [
    ("cymt", "all", ["r1"]), ("cymt", "any", ["r1", "r3", "r4"]),
    ("not-numt", "all", ["r1", "r4"]), ("not-numt", "any", ["r1", "r3", "r4"]),
])
def test_numt_read_list(tmp_path, keep, pairs, expected):
    rows = [("r1", 0, "cymt"), ("r2", 0, "numt"), ("r3", 1, "cymt"), ("r3", 2, "numt"), ("r4", 1, "undetermined"),
            ("r4", 2, "cymt")]
    (tmp_path / "t.tsv").write_text("#read_ID\tread_in_pair\tmt_aln_bp\tmt_mismatch\tmt_identity\tnumt_aln_bp\t"
                                    "numt_mismatch\tnumt_identity\tcandidate\n"
                                    + "".join(f"{r}\t{p}\t1\t0\t1\t1\t0\t1\t{c}\n" for r, p, c in rows))
    counts = mtdna.numt_read_list(tmp_path / "t.tsv", tmp_path / "keep.txt", keep, pairs)
    assert (tmp_path / "keep.txt").read_text().split() == expected
    assert counts["cymt"] == 3 and counts["kept"] == len(expected)


def test_pileups_with_non_utf8_bytes(tmp_path):
    # Seen on real data: "'utf-8' codec can't decode byte 0x81" in circular-merge.
    from mitopiper import pileup
    odd = b"\x81"
    normal = b"".join(b"MT\t%d\tA\t2\t.,\tI" % p + odd + b"\n" for p in range(1, 11))
    (tmp_path / "n.p").write_bytes(normal)
    (tmp_path / "r.p").write_bytes(normal)
    mtdna.combine_pileups(tmp_path / "n.p", tmp_path / "r.p", tmp_path / "c.p", 10, window=2)
    out = (tmp_path / "c.p").read_bytes()
    assert out.count(odd) == 10 and out.splitlines()[4] == normal.splitlines()[4]      # bytes kept as they were
    assert len(list(pileup.read_pileup(tmp_path / "c.p"))) == 10


def test_consensus_takes_the_samples_major_alleles(tmp_path):
    mt = fasta(tmp_path / "mt.fa", "MT human", "AAAAAAAAAA")
    rows = [("MT", 2, "A", 60, "G" * 30 + "g" * 30),        # the sample's allele: substituted
            ("MT", 4, "A", 10, "G" * 10),                     # too little depth
            ("MT", 6, "A", 60, "G" * 30 + "." * 30),          # no majority
            ("MT", 8, "A", 60, "." * 60)]                     # the reference allele
    (tmp_path / "s.pileup").write_text("".join(f"{c}\t{p}\t{r}\t{d}\t{b}\t{'I' * d}\n" for c, p, r, d, b in rows))
    changes = mtdna.consensus(tmp_path / "s.pileup", mt, tmp_path / "c.fa")
    assert changes == [(2, "A", "G")]
    with open(tmp_path / "c.fa") as handle:
        assert list(mtdna.read_fasta(handle)) == [("MT human", "AGAAAAAAAA")]   # same length, same name


MT20 = "ACGTTGCAAGGCTTACGATC"


def _sam(*lines):
    return io.StringIO("@SQ\tSN:numt1\tLN:100\n" + "".join("\t".join(map(str, l)) + "\n" for l in lines))


def test_numt_alleles_reads_the_numt_bases_from_md_tags():
    records = {"mt": MT20}
    found = mtdna.numt_alleles(_sam(
        # secondary alignment (no SEQ, as bwa mem -a writes them): the NUMT has A at mt position 5 (T)
        ("mt_1", 256, "numt1", 1, 0, "20M", "*", 0, 0, "*", "*", "MD:Z:4A15"),
        # reverse strand: the first aligned base is mt position 20 (C); the NUMT's T is an A on the mtDNA
        ("mt_1", 16, "numt1", 30, 0, "20M", "*", 0, 0, "*", "*", "MD:Z:T19"),
        # wraps around the origin: fragment mt_15 of 10 bases ends at position 4 (T); the NUMT has G
        ("mt_15", 0, "numt1", 60, 0, "2S8M", "*", 0, 0, "*", "*", "MD:Z:7G"),
        # ignored: an mtDNA contig, too short, unmapped
        ("mt_1", 0, "chrM", 1, 0, "20M", "*", 0, 0, "*", "*", "MD:Z:0G19"),
        ("mt_1", 0, "numt1", 1, 0, "5M15S", "*", 0, 0, "*", "*", "MD:Z:0G4"),
        ("mt_1", 4, "*", 0, 0, "*", "*", 0, 0, "*", "*"),
    ), records, min_length=8)
    assert found[("mt", 5)] == {"T", "A"} and found[("mt", 20)] == {"C", "A"}
    assert found[("mt", 4)] == {"T", "G"} and found[("mt", 1)] == {"A"}
    assert found[("mt", 17)] == {"G"}                     # the first aligned base of mt_15, after its clipped 2


def test_numt_alleles_command(tmp_path, capsys):
    (tmp_path / "mt.fa").write_text(f">mt test\n{MT20}\n")
    (tmp_path / "a.sam").write_text(_sam(("mt_1", 0, "numt1", 1, 0, "20M", "*", 0, 0, "*", "*", "MD:Z:4A15")).getvalue())
    assert main(["mt", "numt-alleles", str(tmp_path / "a.sam"), str(tmp_path / "mt.fa"), str(tmp_path / "a.tsv"),
                 "--min-length", "10"]) == 0
    lines = (tmp_path / "a.tsv").read_text().splitlines()
    assert lines == ["chrom\tposition\tref\tnumt_bases", "mt\t5\tT\tA"]     # only where a NUMT differs
    assert "1 mtDNA position(s)" in capsys.readouterr().err
