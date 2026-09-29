"""The tutorial simulator and the benchmark comparison (runs without bwa/samtools)."""
import csv
import gzip
import importlib.util
import subprocess
from pathlib import Path

from mitopiper import mtdna, pileup

REPO = Path(__file__).resolve().parent.parent


def load(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


simulate = load(REPO / "tutorial" / "simulate.py")
compare = load(REPO / "benchmark" / "compare.py")


def rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_simulated_dataset_is_consistent(tmp_path):
    simulate.main([str(tmp_path / "d"), "--depth", "4"])
    ref = tmp_path / "d" / "reference"
    mt = dict(mtdna.read_fasta(open(ref / "mt.fa")))["synthMT"]
    nuclear = dict(mtdna.read_fasta(open(ref / "nuclear.fa")))
    assert len(mt) == 16000 and nuclear["chrM"] == mt and len(nuclear["chr2"]) == 153000
    truth = rows(tmp_path / "d" / "truth.tsv")
    assert all(t["ref"] == mt[int(t["position"]) - 1] != t["alt"] for t in truth)
    numt = [t for t in truth if t["sample"] == "patient" and t["type"] == "numt"]
    assert len(numt) == 90 and all(2000 < int(t["position"]) <= 5000 for t in numt)
    assert nuclear["chr2"][60000:63000] != mt[2000:5000]                  # the diverged NUMT copy
    with gzip.open(tmp_path / "d" / "fastq" / "patient_1.fastq.gz", "rt") as r1, \
            gzip.open(tmp_path / "d" / "fastq" / "patient_2.fastq.gz", "rt") as r2:
        a, b = r1.read().splitlines(), r2.read().splitlines()
    assert len(a) == len(b) and len(a) // 4 > 300 and a[0].endswith("/1") and b[0].endswith("/2")
    assert pileup.read_annotation(ref / "mt.gff3")[0] == ("synthMT", 700, 1600, "gene", "RNR1", "+")
    again = tmp_path / "again"
    simulate.main([str(again), "--depth", "4"])
    assert (again / "truth.tsv").read_text() == (tmp_path / "d" / "truth.tsv").read_text()   # deterministic


def experiment(root, label, calls, depth=100, origin_depth=50):
    """An experiment folder with a linear/circular pileup and an annotated table of one sample."""
    exp = root / label
    (exp / "pileup").mkdir(parents=True)
    lines, circular = [], []
    for p in range(1, 1001):
        near = p <= 100 or p > 900
        lines.append(f"chrM\t{p}\tA\t{origin_depth if near else depth}\t{'.' * 2}\tII\n")
        circular.append(f"chrM\t{p}\tA\t{depth}\t{'.' * 2}\tII\n")
    (exp / "pileup" / "s_0_HS.pileup").write_text("".join(lines))
    (exp / "pileup" / "s_0_HS.circular.pileup").write_text("".join(circular))
    pile = exp / "s.pileup"
    pile.write_text("".join(f"chrM\t{p}\tA\t100\t{'.' * (100 - int(100 * calls.get(p, 0)))}"
                            f"{'G' * int(50 * calls.get(p, 0))}{'g' * int(50 * calls.get(p, 0))}\tI\n" for p in range(1, 1001)))
    pileup.annotate_pileup(pile, exp / "annotated" / "s_0_HS.annotated.tsv", sample="s_0_HS")
    return f"{label}={exp}"


def test_compare_reports_origin_concordance_and_accuracy(tmp_path):
    a = experiment(tmp_path, "plain", {10: 0.2, 500: 0.1, 600: 0.3})
    b = experiment(tmp_path, "mapq", {10: 0.2, 500: 0.1})
    (tmp_path / "truth.tsv").write_text("sample\tposition\tref\talt\taf\ttype\n"
                                        "s\t10\tA\tG\t0.2\theteroplasmy\ns\t500\tA\tG\t0.1\theteroplasmy\n"
                                        "s\t600\tA\tG\t\tnumt\n")
    assert compare.main([a, b, "-o", str(tmp_path / "out"), "--truth", str(tmp_path / "truth.tsv"), "--window", "100"]) == 0
    origin = {(r["experiment"], r["pileup"]): float(r["origin_vs_other"]) for r in rows(tmp_path / "out" / "origin_coverage.tsv")}
    assert origin[("plain", "linear")] == 0.5 and origin[("plain", "circular")] == 1.0
    (conc,) = rows(tmp_path / "out" / "concordance.tsv")
    assert (conc["shared"], conc["only_a"], conc["only_b"], conc["jaccard"]) == ("2", "1", "0", "0.6667")
    acc = {r["experiment"]: r for r in rows(tmp_path / "out" / "accuracy.tsv")}
    assert (acc["plain"]["false_positives"], acc["plain"]["false_calls_at_numt_differences"]) == ("1", "1")
    assert (acc["mapq"]["precision"], acc["mapq"]["recall"], acc["mapq"]["af_mean_abs_error"]) == ("1", "1", "0")


def test_run_presets_usage():
    result = subprocess.run(["bash", str(REPO / "benchmark" / "run_presets.sh")], capture_output=True, text=True)
    assert result.returncode == 2 and "run_presets.sh WORKSPACE METADATA SCHEDULER" in result.stdout
