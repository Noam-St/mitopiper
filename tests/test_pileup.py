import csv

import pytest

from mitopiper import pileup
from mitopiper.cli import main


def rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.parametrize("bases, expected", [
    # read start (the mapping quality '.' is not a base), forward match that ends, reverse mismatch
    ("^..$c", {"fwd": {"A": 1, "start": 1, "end": 1, "cov": 1}, "rev": {"C": 1, "cov": 1}}),
    # insertion after a reverse base, deletion placeholders, a deletion's sequence is not coverage
    (",+2ag*#.-3ACG", {"fwd": {"A": 1, "del": 1, "cov": 2}, "rev": {"A": 1, "ins": 1, "del": 1, "cov": 2}}),
    # reference skips (spliced reads) are not coverage; N counts as an allele
    ("><Nn", {"fwd": {"N": 1, "cov": 1}, "rev": {"N": 1, "cov": 1}}),
    ("", {"fwd": {}, "rev": {}}),
])
def test_parse_bases(bases, expected):
    counts = pileup.parse_bases("a", bases)
    for strand in ("fwd", "rev"):
        assert {k: v for k, v in counts[strand].items() if v} == expected[strand]


def test_annotations_from_gff_gtf_and_genbank(tmp_path):
    (tmp_path / "a.gff3").write_text(
        "##gff-version 3\nNC_1\tRefSeq\tregion\t1\t60\t.\t+\t.\tID=r\n"
        "NC_1\tRefSeq\tgene\t5\t20\t.\t+\t.\tID=g1;Name=ND1;gene=ND1\n"
        "NC_1\tRefSeq\tgene\t40\t55\t.\t-\t.\tID=g2;locus_tag=X_mgr02\n")
    (tmp_path / "a.gtf").write_text('NC_1\tx\tgene\t5\t20\t.\t+\t.\tgene_id "g1"; gene_name "ND1";\n')
    (tmp_path / "a.gbk").write_text(
        "LOCUS       NC_1  60 bp    DNA     circular\nVERSION     NC_1.1\nFEATURES             Location/Qualifiers\n"
        "     gene            join(5..10,12..20)\n                     /gene=\"ND1\"\n"
        "     gene            complement(40..55)\n                     /locus_tag=\"X_mgr02\"\nORIGIN\n//\n")
    assert pileup.read_annotation(tmp_path / "a.gff3") == [("NC_1", 5, 20, "gene", "ND1", "+"),
                                                            ("NC_1", 40, 55, "gene", "RNR2", "-")]
    assert pileup.read_annotation(tmp_path / "a.gtf") == [("NC_1", 5, 20, "gene", "ND1", "+")]
    assert pileup.read_annotation(tmp_path / "a.gbk") == [("NC_1.1", 5, 10, "gene", "ND1", "+"),
                                                           ("NC_1.1", 12, 20, "gene", "ND1", "+"),
                                                           ("NC_1.1", 40, 55, "gene", "RNR2", "-")]


def write_pileup(path, lines):
    path.write_text("".join(f"chrM\t{p}\t{ref}\t{depth}\t{bases}\t{'I' * max(depth, 1)}\n" for p, ref, depth, bases in lines))


def test_annotate_reports_coverage_alleles_genes_and_heteroplasmy(tmp_path):
    het = "." * 40 + "," * 40 + "G" * 5 + "g" * 5          # 10% G on both strands
    alt = "T" * 30 + "t" * 30                              # homoplasmic variant
    write_pileup(tmp_path / "s.pileup", [(1, "A", 90, het), (2, "C", 60, alt), (4, "G", 3, "^].$,G")])
    (tmp_path / "a.gff3").write_text("NC_012920.1\tx\tgene\t1\t2\t.\t+\t.\tgene=ATP8\n"
                                     "NC_012920.1\tx\tgene\t2\t3\t.\t+\t.\tgene=ATP6\n")
    summary = pileup.annotate_pileup(tmp_path / "s.pileup", tmp_path / "s.tsv", tmp_path / "a.gff3", sample="S1",
                                     min_depth=50, min_af=0.05, min_strand=2, total_reads=1_000_000)
    assert summary == {"positions": 4, "heteroplasmic": 1, "alt_major": 1}
    r = rows(tmp_path / "s.tsv")
    assert [x["position"] for x in r] == ["1", "2", "3", "4"]       # the missing position 3 is filled in
    assert (r[0]["major"], r[0]["minor"], r[0]["minor_af"], r[0]["heteroplasmic"]) == ("A", "G", "0.111111", "yes")
    assert (r[0]["pos_coverage"], r[0]["neg_coverage"], r[0]["RPM"]) == ("45", "45", "90")
    assert (r[1]["major"], r[1]["alt_major"], r[1]["heteroplasmic"], r[1]["gene"]) == ("T", "yes", "no", "ATP8;ATP6")
    assert r[0]["gene"] == "ATP8"                               # the annotation's NC_012920.1 is this chrM
    assert (r[2]["coverage"], r[2]["major"], r[2]["gene"], r[2]["minor_af"]) == ("0", "", "ATP6", "")
    assert (r[3]["pos_read_start"], r[3]["pos_read_end"], r[3]["pos_coverage"], r[3]["neg_coverage"]) == ("1", "1", "2", "1")
    assert r[0]["sample"] == "S1"


def test_annotate_flip_strands_and_strand_filter(tmp_path):
    write_pileup(tmp_path / "s.pileup", [(1, "A", 100, "." * 50 + "," * 40 + "G" * 10)])   # minor only forward
    pileup.annotate_pileup(tmp_path / "s.pileup", tmp_path / "a.tsv", min_depth=10)
    pileup.annotate_pileup(tmp_path / "s.pileup", tmp_path / "b.tsv", min_depth=10, flip_strands=True)
    a, b = rows(tmp_path / "a.tsv")[0], rows(tmp_path / "b.tsv")[0]
    assert a["heteroplasmic"] == "no" and (a["minor_pos"], a["minor_neg"]) == ("10", "0")   # strand bias
    assert (a["pos_coverage"], b["pos_coverage"]) == ("60", "40")
    assert a["z"] == ""                                          # constant coverage: no z-score


def test_htseq_total(tmp_path):
    (tmp_path / "c.txt").write_text("ND1\t10\nND2\t5\n__no_feature\t3\n__ambiguous\t7\n__not_aligned\t100\n")
    assert pileup.htseq_total(tmp_path / "c.txt") == 18
    assert pileup.htseq_total(tmp_path / "c.txt", include_unaligned=True, extra=2) == 120


def test_heteroplasmy_collects_sites_and_matrix(tmp_path, capsys):
    tables = []
    for name, g in (("S1", 10), ("S2", 0)):
        write_pileup(tmp_path / f"{name}.pileup", [(1, "A", 100, "." * 45 + "," * 45 + "G" * (g // 2) + "g" * (g // 2)
                                                    + "." * (10 - g)), (2, "C", 60, "t" * 60)])
        pileup.annotate_pileup(tmp_path / f"{name}.pileup", tmp_path / f"{name}.tsv", sample=name)
        tables.append(str(tmp_path / f"{name}.tsv"))
    assert main(["mt", "heteroplasmy", *tables, str(tmp_path / "missing.tsv"), "-o", str(tmp_path / "out")]) == 0
    assert "missing.tsv does not exist" in capsys.readouterr().err
    sites = rows(tmp_path / "out" / "heteroplasmy_sites.tsv")
    assert [(s["sample"], s["position"], s["type"]) for s in sites] == [
        ("S1", "1", "heteroplasmy"), ("S1", "2", "variant"), ("S2", "2", "variant")]
    matrix = rows(tmp_path / "out" / "heteroplasmy_matrix.tsv")
    assert [(m["position"], m["S1"], m["S2"]) for m in matrix] == [("1", "0.1", "0"), ("2", "0", "0")]


def test_annotate_cli(tmp_path, capsys):
    write_pileup(tmp_path / "x_0_HS.circular.pileup", [(1, "A", 2, ".,")])
    assert main(["mt", "annotate", str(tmp_path / "x_0_HS.circular.pileup"), str(tmp_path / "o.tsv")]) == 0
    assert rows(tmp_path / "o.tsv")[0]["sample"] == "x_0_HS"
    assert "1 positions" in capsys.readouterr().err


def test_variant_needs_depth_and_a_majority():
    # A single read with an error at a barely covered position was called a homoplasmic variant.
    def counts(**alleles):
        base = {a: 0 for a in pileup.COUNTED}
        base.update(ins=0)
        pos = dict(base, **{a: v // 2 for a, v in alleles.items()})
        neg = dict(base, **{a: v - v // 2 for a, v in alleles.items()})
        return pos, neg

    thresholds = (50, 0.03, 2)
    assert pileup._allele_fields("A", *counts(G=1), thresholds)["alt_major"] == "no"          # depth 1
    assert pileup._allele_fields("A", *counts(G=60), thresholds)["alt_major"] == "yes"        # a real variant
    assert pileup._allele_fields("A", *counts(G=40, A=40), thresholds)["alt_major"] == "no"   # a tie is no majority
    assert pileup._allele_fields("A", *counts(G=45, A=40, C=15), thresholds)["alt_major"] == "no"  # 45%: no majority
