from pipemaker.local_samples import scan_fastq_folder


def touch(folder, *names):
    for n in names:
        (folder / n).write_text("")


def test_pairs_and_lanes(tmp_path):
    touch(tmp_path, "A_L001_R1_001.fastq.gz", "A_L001_R2_001.fastq.gz", "A_L002_R1_001.fastq.gz",
          "A_L002_R2_001.fastq.gz", "B.fq.gz", "notes.txt")
    runs = scan_fastq_folder(tmp_path)
    assert [(r, s, len(f)) for r, s, f in runs] == [
        ("A_L001_001", "A", 2), ("A_L002_001", "A", 2), ("B", "B", 1)
    ]


def test_ignore_and_force_paired(tmp_path):
    touch(tmp_path, "x_1.fastq", "x_2.fastq")
    assert len(scan_fastq_folder(tmp_path, ignore_paired=True)) == 2
    other = tmp_path / "o"
    other.mkdir()
    touch(other, "a.fastq", "b.fastq")
    (run,) = scan_fastq_folder(other, force_paired=True)
    assert len(run[2]) == 2
