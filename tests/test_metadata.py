import pytest

from mitopiper.metadata import (
    MetadataError,
    ena_fastq_urls,
    load_samples,
    organism_initials,
    parse_selection,
    select_read_files,
)


def test_ena_groups_runs_by_sample_and_skips_missing_links(data_dir, caplog):
    samples = load_samples(data_dir / "ena_runs.tsv")
    assert [s.accession for s in samples] == ["SAMEA1", "SAMEA3"]
    fish = samples[0]
    assert fish.run_accessions == ["ERR100", "ERR101"]
    assert fish.paired
    assert fish.name == "SAMEA1_0_DR"
    assert "ERR102" in caplog.text


def test_three_file_runs_keep_only_mates(data_dir):
    fish = load_samples(data_dir / "ena_runs.tsv")[0]
    assert fish.runs[1].files == [
        "ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/ERR101_1.fastq.gz",
        "ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/ERR101_2.fastq.gz",
    ]
    assert fish.runs[1].md5 == ["ccc", "ddd"]


def test_ignore_paired_keeps_first_mate(data_dir):
    fish = load_samples(data_dir / "ena_runs.tsv", ignore_paired=True)[0]
    assert not fish.paired
    assert all(len(r.files) == 1 and r.files[0].endswith("_1.fastq.gz") for r in fish.runs)


def test_sra_runinfo_derives_ena_urls(data_dir):
    samples = load_samples(data_dir / "SraRunInfo.csv")
    assert [s.accession for s in samples] == ["SAMN01", "SAMN02", "SAMN03"]
    assert samples[0].organism == "Mus_musculus"
    assert samples[0].paired and not samples[1].paired
    assert samples[0].runs[0].files[0] == "ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/007/SRR1234567/SRR1234567_1.fastq.gz"
    assert samples[1].runs[0].files == ["ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/078/SRR12345678/SRR12345678.fastq.gz"]
    assert samples[2].runs[0].files == ["ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/SRR123456/SRR123456.fastq.gz"]
    assert samples[0].title == "liver_1"


def test_legacy_sample_level_table(data_dir):
    (sample,) = load_samples(data_dir / "legacy_sample_level.tsv")
    assert sample.run_accessions == ["SRR1", "SRR2"]
    assert sample.paired


def test_unknown_format(tmp_path):
    path = tmp_path / "x.tsv"
    path.write_text("a\tb\n1\t2\n")
    with pytest.raises(MetadataError):
        load_samples(path)


@pytest.mark.parametrize(
    "run,expected",
    [
        ("ERR123456", "ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456/ERR123456.fastq.gz"),
        ("SRR1234567", "ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/007/SRR1234567/SRR1234567.fastq.gz"),
        ("DRR123456789", "ftp.sra.ebi.ac.uk/vol1/fastq/DRR123/789/DRR123456789/DRR123456789.fastq.gz"),
    ],
)
def test_ena_url_layout(run, expected):
    assert ena_fastq_urls(run, paired=False) == [expected]


def test_select_read_files_single():
    assert select_read_files(["a.fastq.gz"]) == (["a.fastq.gz"], [""])


@pytest.mark.parametrize(
    "text,expected",
    [("all", [0, 1, 2, 3, 4]), ("2", [2]), ("1-3", [1, 2, 3]), ("0 4", [0, 4]), ("0,4,4", [0, 4]), ("3-1", [1, 2, 3])],
)
def test_parse_selection(text, expected):
    assert parse_selection(text, 5) == expected


@pytest.mark.parametrize("text", ["7", "a", "1-9"])
def test_parse_selection_errors(text):
    with pytest.raises(ValueError):
        parse_selection(text, 5)


@pytest.mark.parametrize("name,initials", [("Homo_sapiens", "HS"), ("Homo sapiens", "HS"), ("Escherichia", "E")])
def test_initials(name, initials):
    assert organism_initials(name) == initials
