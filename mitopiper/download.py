"""Downloads of the sample jobs that survive a flaky archive: ``mitopiper download``.

The sample jobs fetch each FASTQ file with this command instead of calling wget directly:

1. the file is downloaded with the preset's tool (``wget -c``) up to ``tries`` times. A checksum that
   does not match counts as a failed attempt: the file is deleted (``wget -c`` would take a full-size
   corrupt file as complete) and downloaded again. Other failures keep the partial file, which the
   next attempt resumes;
2. ENA serves its files over http(s) and ftp: the attempts go round the URL as given, then https:// and
   ftp://, so a proxy or a block that refuses one protocol (an HTTP 403 seen on the BGU cluster) does
   not fail all of them;
3. the wait between attempts doubles from ``wait`` seconds, up to 10 minutes;
4. when every attempt failed, the run accession is known and the SRA Toolkit is installed (it is
   optional), the run is taken from the SRA itself with ``fasterq-dump`` and compressed. Without it the
   download fails after the attempts, as a plain wget would. Both mates of a run come from one dump, kept
   in ``<folder>/.sra/<run>/`` until the second mate is taken. The archive's MD5 is of its own file,
   so there is no checksum to compare; fasterq-dump checks its download itself.
"""

from __future__ import annotations

import hashlib
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Callable, TextIO

MAX_WAIT = 600
ENA_HOSTS = ("ftp.sra.ebi.ac.uk/", "ftp.ebi.ac.uk/")
RUN = re.compile(r"^[SED]RR\d+$")
SRA_TOOL = "fasterq-dump"


class DownloadError(Exception):
    """Every way of getting the file failed."""


def md5sum(path, chunk: int = 1 << 20) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def url_variants(url: str) -> list:
    """The URL, then the same ENA file over https:// and ftp:// (other URLs only as given)."""
    rest = re.sub(r"^[a-z]+://", "", url)
    if not rest.startswith(ENA_HOSTS):
        return [url]
    variants = [url] + [f"{scheme}://{rest}" for scheme in ("https", "ftp")]
    return list(dict.fromkeys(variants))


def fetch_with_retries(url: str, dest: str, md5: str = "", downloader: str = "wget", downloader_args: str = "",
                       tries: int = 10, wait: int = 30, runner: Callable = subprocess.run,
                       sleep: Callable = time.sleep, log: TextIO = sys.stderr) -> bool:
    """Download ``url`` to ``dest`` (checked against ``md5`` when given). Returns whether it worked."""
    variants = url_variants(url)
    for attempt in range(1, tries + 1):
        source = variants[(attempt - 1) % len(variants)]
        print(f"[mitopiper] download {attempt}/{tries}: {source}", file=log, flush=True)
        rc = runner([downloader, *shlex.split(downloader_args), "-O", dest, source]).returncode
        if rc != 0:
            print(f"[mitopiper] {downloader} failed (exit code {rc})", file=log, flush=True)
        elif not md5:
            return True
        else:
            found = md5sum(dest)
            if found == md5.lower():
                print(f"[mitopiper] {dest}: checksum OK", file=log, flush=True)
                return True
            print(f"[mitopiper] checksum failed for {dest} ({found}, expected {md5}): deleting it", file=log, flush=True)
            Path(dest).unlink(missing_ok=True)
        if attempt < tries:
            pause = min(wait * 2 ** (attempt - 1), MAX_WAIT)
            print(f"[mitopiper] waiting {pause} s before the next attempt", file=log, flush=True)
            sleep(pause)
    return False


def from_sra(run: str, mate: int, dest: str, threads: int = 1, tool: str = SRA_TOOL,
             runner: Callable = subprocess.run, log: TextIO = sys.stderr) -> None:
    """Write mate ``mate`` (0: single-end) of ``run`` from the SRA to ``dest`` (gzipped)."""
    cache = Path(dest).parent / ".sra" / run
    done = cache / ".complete"
    if not done.exists():
        shutil.rmtree(cache, ignore_errors=True)
        cache.mkdir(parents=True)
        print(f"[mitopiper] taking {run} from the SRA ({tool})", file=log, flush=True)
        command = [tool, "--split-3", "--threads", str(threads), "--outdir", str(cache), "--temp", str(cache), run]
        if runner(command).returncode != 0:
            raise DownloadError(f"{tool} {run} failed too")
        done.touch()
    fastq = cache / (f"{run}_{mate}.fastq" if mate else f"{run}.fastq")
    if not fastq.is_file():
        made = ", ".join(sorted(p.name for p in cache.glob("*.fastq"))) or "nothing"
        raise DownloadError(f"{tool} {run} did not write {fastq.name} (it wrote {made})")
    compressor = ["pigz", "-p", str(threads)] if shutil.which("pigz") else ["gzip"]
    with open(dest, "wb") as out:
        if runner([*compressor, "-c", str(fastq)], stdout=out).returncode != 0:
            raise DownloadError(f"could not compress {fastq}")
    fastq.unlink()
    if not list(cache.glob(f"{run}_[12].fastq")):          # every mate taken (the unpaired reads are not used)
        shutil.rmtree(cache, ignore_errors=True)
    print(f"[mitopiper] {dest}: from the SRA ({run}); the archive's MD5 is of its own file, so it was not "
          "compared", file=log, flush=True)


def download(url: str, dest: str, md5: str = "", downloader: str = "wget", downloader_args: str = "",
             tries: int = 10, wait: int = 30, run: str = "", mate: int = 0, threads: int = 1, sra: bool = True,
             runner: Callable = subprocess.run, sleep: Callable = time.sleep, log: TextIO = sys.stderr) -> str:
    """Get one FASTQ file, from the archive or else from the SRA. Returns "archive" or "sra"."""
    Path(dest).parent.mkdir(parents=True, exist_ok=True)
    if fetch_with_retries(url, dest, md5, downloader, downloader_args, tries, wait, runner, sleep, log):
        return "archive"
    if not (sra and RUN.match(run or "")):
        raise DownloadError(f"download failed {tries} times: {url}"
                            + ("" if sra else " (the SRA fallback is switched off: sra_fallback)"))
    if not shutil.which(SRA_TOOL):
        raise DownloadError(f"download failed {tries} times: {url} (the SRA was not tried: {SRA_TOOL} (SRA Toolkit) "
                            "is not installed; it is optional: conda install -c bioconda sra-tools)")
    print(f"[mitopiper] download failed {tries} times: trying the SRA", file=log, flush=True)
    Path(dest).unlink(missing_ok=True)
    from_sra(run, mate, dest, threads, runner=runner, log=log)
    return "sra"


# ---------------------------------------------------------------------- CLI
def add_parser(sub) -> None:
    d = sub.add_parser("download", help="download a FASTQ file with retries, checksum and an SRA fallback "
                                        "(used by the sample jobs)")
    d.add_argument("url")
    d.add_argument("dest")
    d.add_argument("--md5", default="", help="expected MD5 of the file")
    d.add_argument("--downloader", default="wget", help="download tool, called as TOOL ARGS -O DEST URL")
    d.add_argument("--downloader-args", default="", help="its options, e.g. '-c'")
    d.add_argument("--tries", type=int, default=10, help="attempts before the SRA fallback, if the SRA Toolkit is "
                                                          "installed (default 10)")
    d.add_argument("--wait", type=int, default=30, help="seconds before the 2nd attempt, doubling up to 600")
    d.add_argument("--run", default="", help="run accession (SRR/ERR/DRR), for the SRA fallback")
    d.add_argument("--mate", type=int, default=0, choices=(0, 1, 2), help="1 or 2 for paired-end, 0 for single-end")
    d.add_argument("--threads", type=int, default=1)
    d.add_argument("--no-sra", action="store_true", help="no SRA fallback")


def run(args) -> int:
    try:
        download(args.url, args.dest, args.md5, args.downloader, args.downloader_args, max(args.tries, 1), args.wait,
                 args.run, args.mate, args.threads, sra=not args.no_sra)
    except DownloadError as exc:
        print(f"[mitopiper] ERROR: {exc}", file=sys.stderr)
        return 1
    return 0
