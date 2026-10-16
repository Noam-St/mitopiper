"""mitopiper download: retries (a wrong checksum included), other protocols, the SRA fallback."""
import gzip
import hashlib
import io
import os
import subprocess

import pytest

from mitopiper import download as dl
from mitopiper.toolcheck import required_tools

GOOD = b"@r\nACGT\n+\nIIII\n"


def fake_wget(outcomes, calls):
    """A runner standing in for wget: each call writes the next outcome (bytes, or an exit code)."""
    def run(command, **kwargs):
        calls.append(command)
        outcome = outcomes.pop(0)
        dest = command[command.index("-O") + 1]
        if isinstance(outcome, int):
            return subprocess.CompletedProcess(command, outcome)
        with open(dest, "wb") as out:
            out.write(outcome)
        return subprocess.CompletedProcess(command, 0)
    return run


def test_url_variants_go_round_the_protocols_of_ena():
    assert dl.url_variants("ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1.fastq.gz") == [
        "ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1.fastq.gz", "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1.fastq.gz",
        "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1.fastq.gz"]
    assert dl.url_variants("https://ftp.sra.ebi.ac.uk/x.gz") == ["https://ftp.sra.ebi.ac.uk/x.gz", "ftp://ftp.sra.ebi.ac.uk/x.gz"]
    assert dl.url_variants("https://example.org/x.gz") == ["https://example.org/x.gz"]


def test_wrong_checksums_and_failures_are_retried_with_growing_waits(tmp_path):
    calls, waits = [], []
    dest = tmp_path / "x.fastq.gz"
    runner = fake_wget([b"corrupt", 8, b"corrupt again", GOOD], calls)
    log = io.StringIO()
    how = dl.download("ftp.sra.ebi.ac.uk/a.gz", str(dest), md5=hashlib.md5(GOOD).hexdigest(), downloader_args="-c",
                      tries=10, wait=30, runner=runner, sleep=waits.append, log=log)
    assert how == "archive" and dest.read_bytes() == GOOD
    assert [c[-1] for c in calls] == ["ftp.sra.ebi.ac.uk/a.gz", "https://ftp.sra.ebi.ac.uk/a.gz",
                                      "ftp://ftp.sra.ebi.ac.uk/a.gz", "ftp.sra.ebi.ac.uk/a.gz"]
    assert calls[0][:4] == ["wget", "-c", "-O", str(dest)]
    assert waits == [30, 60, 120] and log.getvalue().count("checksum failed") == 2


def test_waits_are_capped():
    waits = []
    dl.fetch_with_retries("u", "/nonexistent/x", tries=8, wait=100, runner=fake_wget([1] * 8, []),
                          sleep=waits.append, log=io.StringIO())
    assert waits == [100, 200, 400, 600, 600, 600, 600]


def fake_fasterq_dump(tmp_path, monkeypatch, paired=True):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    count = tmp_path / "dumps"
    body = (f'#!/bin/sh\necho x >> "{count}"\n'
            'while [ $# -gt 1 ]; do case "$1" in --outdir) out=$2; shift;; esac; shift; done; run=$1\n')
    body += ('printf "@a/1\\nACGT\\n+\\nIIII\\n" > "$out/${run}_1.fastq"\nprintf "@a/2\\nTTTT\\n+\\nIIII\\n" > "$out/${run}_2.fastq"\n'
             'printf "@u\\nGG\\n+\\nII\\n" > "$out/$run.fastq"\n' if paired else 'printf "@s\\nACGT\\n+\\nIIII\\n" > "$out/$run.fastq"\n')
    (bin_dir / "fasterq-dump").write_text(body)
    (bin_dir / "fasterq-dump").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    return count


def test_after_every_attempt_failed_both_mates_come_from_one_sra_dump(tmp_path, monkeypatch):
    dumps = fake_fasterq_dump(tmp_path, monkeypatch)
    real = subprocess.run

    def runner(command, **kwargs):              # wget always refused; everything else runs
        return subprocess.CompletedProcess(command, 8) if command[0] == "wget" else real(command, **kwargs)

    out = tmp_path / "fastq"
    for mate in (1, 2):
        how = dl.download("ftp.sra.ebi.ac.uk/SRR9_%d.fastq.gz" % mate, str(out / f"S_{mate}.fastq.gz"), md5="0" * 32,
                          tries=3, wait=0, run="SRR9", mate=mate, runner=runner, sleep=lambda s: None, log=io.StringIO())
        assert how == "sra"
    assert gzip.open(out / "S_1.fastq.gz").read() == b"@a/1\nACGT\n+\nIIII\n"
    assert gzip.open(out / "S_2.fastq.gz").read() == b"@a/2\nTTTT\n+\nIIII\n"
    assert dumps.read_text().count("x") == 1 and not (out / ".sra" / "SRR9").exists()


def test_single_end_from_the_sra(tmp_path, monkeypatch):
    fake_fasterq_dump(tmp_path, monkeypatch, paired=False)
    dest = tmp_path / "S.fastq.gz"
    how = dl.download("u", str(dest), tries=1, run="ERR5", mate=0, runner=lambda c, **k: (
        subprocess.CompletedProcess(c, 4) if c[0] == "wget" else subprocess.run(c, **k)), log=io.StringIO())
    assert how == "sra" and gzip.open(dest).read().startswith(b"@s")


@pytest.mark.parametrize("run, sra, message", [
    ("", True, "download failed 2 times"),                       # no run accession: nothing to fall back on
    ("SRR9", False, "SRA fallback is switched off"),
    ("SRR9", True, "fasterq-dump (SRA Toolkit) is not installed"),
])
def test_when_nothing_works_the_error_says_why(tmp_path, monkeypatch, run, sra, message):
    monkeypatch.setenv("PATH", str(tmp_path))                   # no fasterq-dump
    with pytest.raises(dl.DownloadError, match=message.replace("(", r"\(").replace(")", r"\)")):
        dl.download("u", str(tmp_path / "x.gz"), tries=2, wait=0, run=run, sra=sra,
                    runner=fake_wget([8, 8], []), sleep=lambda s: None, log=io.StringIO())


def test_the_tool_check_sees_the_downloader():
    line = 'mitopiper download --downloader wget --downloader-args=-c --md5 abc --tries 10 --wait 30 "u" "d"'
    assert set(required_tools({"wget": [line]})) == {"mitopiper", "wget"}
