"""End-to-end runs of the bundled mtDNA presets with stand-in aligners.

bwa and samtools are replaced by small scripts (the tests are about the wiring, not the
aligners): the stand-in `samtools mpileup` prints one line per reference position with depth 1
for the normal and 2 for the rotated reference, so the merged pileup shows where every position
came from. `mitopiper mt ...` runs for real.
"""
import os
import shutil
import time
import subprocess
import sys
from pathlib import Path

import pytest

from mitopiper import monitor as mon
from mitopiper import mtdna
from mitopiper.cli import main
from mitopiper.config import Workspace
from mitopiper.local_samples import make_tsv
from mitopiper.references import ReferenceRegistry, ReferenceSet

REPO = Path(__file__).resolve().parent.parent
MT = "ACGT" * 7 + "AC"  # 30 bases

FAKE_SAMTOOLS = r"""#!/bin/sh
cmd=$1; shift
case "$cmd" in
  faidx) out=""; f=$1; while [ $# -gt 0 ]; do case "$1" in -o) out=$2; shift;; esac; shift; done
         if [ -n "$out" ]; then printf '>chr1:1-10\nAAAAAAAAAA\n' > "$out"; else touch "$f.fai"; fi ;;
  index) touch "$1.bai" ;;
  fastq) echo "@r"; echo "ACGT"; echo "+"; echo "IIII" ;;
  sort)  out=""; in="-"; while [ $# -gt 0 ]; do case "$1" in -o) out=$2; shift;; *) in=$1;; esac; shift; done
         if [ "$in" = - ]; then cat > "$out"; else cp "$in" "$out"; fi ;;
  view)  out=""; while [ $# -gt 0 ]; do case "$1" in -o) out=$2; shift;; -) cat > /dev/null;; esac; shift; done
         if [ -n "$out" ]; then echo bam > "$out"; else echo bam; fi ;;
  mpileup) ref=""; while [ $# -gt 0 ]; do case "$1" in -f) ref=$2; shift;; esac; shift; done
         depth=1; case "$ref" in *rotated*) depth=2;; esac
         awk -v d=$depth '/^>/{split(substr($0,2),h," "); name=h[1]; next} {seq=seq $0}
                          END{for(i=1;i<=length(seq);i++) printf "%s\t%d\t%s\t%d\t.\tI\n", name, i, substr(seq,i,1), d}' "$ref" ;;
  *) exit 0 ;;
esac
"""
FAKE_BWA = r"""#!/bin/sh
[ "$1" = index ] && touch "$2.bwt"
case " $* " in *" - "*) cat > /dev/null ;; esac   # read the piped reads, as bwa does (no SIGPIPE upstream)
echo "@HD"
case "$*" in *fragments.fa*) printf '@SQ\tSN:chr1\tLN:10\nNC_012920.1_1\t0\tchr1\t3\t0\t60M\t*\t0\t0\t*\t*\n' ;; esac
"""
# Numt Parser's table: r1 is mtDNA, r2 a NUMT, r3 one mate of each.
FAKE_NUMT_PARSER = r"""#!/bin/sh
while [ $# -gt 0 ]; do case "$1" in --outfile) out=$2; shift;; esac; shift; done
printf '#read_ID\tread_in_pair\tmt_aln_bp\tmt_mismatch\tmt_identity\tnumt_aln_bp\tnumt_mismatch\tnumt_identity\tcandidate\n' > "$out"
printf 'r1\t0\t4\t0\t1.0\tNone\tNone\tNone\tcymt\nr2\t0\t4\t1\t0.75\t4\t0\t1.0\tnumt\n' >> "$out"
printf 'r3\t1\t4\t0\t1.0\t4\t1\t0.75\tcymt\nr3\t2\t4\t1\t0.75\t4\t0\t1.0\tnumt\n' >> "$out"
"""


@pytest.fixture
def cluster(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    tools = {"samtools": FAKE_SAMTOOLS, "bwa": FAKE_BWA, "numt_parser.py": FAKE_NUMT_PARSER,
             "mitopiper": f'#!/bin/sh\nPYTHONPATH="{REPO}" exec "{sys.executable}" -m mitopiper "$@"\n'}
    for name, body in tools.items():
        (bin_dir / name).write_text(body)
        (bin_dir / name).chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    fq = tmp_path / "fq"
    fq.mkdir()
    (fq / "a.fastq").write_text("@r\nACGT\n+\nIIII\n")
    make_tsv(fq, "Homo sapiens", tmp_path / "exp.tsv")
    (tmp_path / "mt.fa").write_text(f">NC_012920.1 Homo sapiens mitochondrion\n{MT}\n")
    return ws, tmp_path


def run_preset(ws, tmp_path, preset, steps, shell="bash"):
    wf = ws.load_workflow(preset)
    wf.set_general("environment", shell)
    wf.set_general("mt_rotation", 5)
    wf.set_general("mt_junction_window", 3)
    ws.write_workflow(wf, preset)
    code = main(["-w", str(ws.root), "build", "--no-check", str(tmp_path / "exp.tsv"), "-s", "local", "-p", preset,
                 "--steps", ",".join(steps)])
    if code == 0:
        subprocess.run(["bash", str(ws.root / "master.sh")], check=True, capture_output=True)
    return code


def merged_depths(ws):
    rows = [line.split("\t") for line in (ws.output_dir / "exp" / "pileup" / "a_0_HS.circular.pileup").read_text().splitlines()]
    assert [int(r[1]) for r in rows] == list(range(1, 31))
    return {int(r[1]): int(r[3]) for r in rows}


CIRCULAR = ["align-rotated", "index-rotated", "mpileup-rotated", "circular-merge"]


def test_mt_variants_merges_the_rotated_alignment_around_the_origin(cluster):
    ws, tmp_path = cluster
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath=str(tmp_path / "mt.fa")))
    (tmp_path / "mt.gff3").write_text("NC_012920.1\tRefSeq\tgene\t2\t4\t.\t+\t.\tgene=TRNF\n")
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath=str(tmp_path / "mt.fa"), gtfpath=str(tmp_path / "mt.gff3")))
    steps = ["wget", "bwa-align", "view", "sort", "index", "mpileup"] + CIRCULAR + ["annotate", "heteroplasmy"]
    assert run_preset(ws, tmp_path, "mt-variants", steps) == 0
    states = {j.name: j.state for j in mon.collect(ws, ["exp"], mon.ErrorScanner(), lambda s, i: None)}
    assert all(s == mon.DONE for s in states.values()), states
    rotated = ws.references_dir / "derived" / "Homo_sapiens" / "mt.rotated5.fa"
    assert mtdna.sequence_length(rotated) == 30 and (rotated.parent / "mt.rotated5.fa.bwt").exists()
    depths = merged_depths(ws)
    assert [p for p, d in depths.items() if d == 2] == [1, 2, 3, 28, 29, 30]   # next to the origin: rotated alignment
    with open(ws.output_dir / "exp" / "annotated" / "a_0_HS.annotated.tsv") as handle:
        table = list(__import__("csv").DictReader(handle, delimiter="\t"))
    assert [(r["position"], r["coverage"], r["gene"]) for r in table[:5]] == [
        ("1", "2", ""), ("2", "2", "TRNF"), ("3", "2", "TRNF"), ("4", "1", "TRNF"), ("5", "1", "")]
    assert (ws.output_dir / "exp" / "heteroplasmy" / "heteroplasmy_matrix.tsv").exists()
    script = (ws.output_dir / "exp" / "scripts" / "e_exp.sh").read_text()
    assert f"mitopiper mt heteroplasmy -o {ws.output_dir}/exp/heteroplasmy {ws.output_dir}/exp/annotated/a_0_HS.annotated.tsv" in script


def test_annotate_leaves_out_a_missing_annotation(tmp_path):
    from mitopiper.builder import build_sample_job
    from mitopiper.metadata import Run, Sample

    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("mt-variants")
    sample = Sample(0, "S1", "Homo_sapiens", [Run("r", ["f/x.fastq.gz"])])
    job = build_sample_job(sample, wf, str(tmp_path / "ws/output/E"), ReferenceSet(refpath="/refs/mt.fa"))
    (line,) = {s.step: s.lines for s in job.steps}["annotate"]
    assert "--annotation" not in line and "--sample S1_0_HS --min-depth 50 --min-af 0.03 --min-strand 2" in line


@pytest.mark.parametrize("shell", ["bash"] + (["tcsh"] if __import__("shutil").which("tcsh") else []))
def test_mt_variants_numt_aligns_to_nuclear_plus_mt(cluster, shell):
    ws, tmp_path = cluster
    (tmp_path / "nuc.fa").write_text(">chr1\nAAAAAAAAAA\n>chrM\nCCCCCCCC\n")
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath=str(tmp_path / "mt.fa"), nuc_refpath=str(tmp_path / "nuc.fa")))
    steps = ["wget", "align-combined", "index-combined", "mt-reads", "index-mt", "mpileup"] + CIRCULAR
    assert run_preset(ws, tmp_path, "mt-variants-numt", steps, shell) == 0
    combined = ws.references_dir / "derived" / "Homo_sapiens" / "mt.with_nuclear.fa"
    assert mtdna.contigs(combined) == ["chr1", "NC_012920.1"]               # the nuclear chrM copy was left out
    script = (ws.output_dir / "exp" / "scripts" / "s_a_0.sh").read_text()
    assert "-q 20 -o" in script and f"`mitopiper mt contigs {tmp_path}/mt.fa`" in script
    assert "samtools fastq" in script and "mt.rotated5.fa -" in script    # the rotated realignment uses the mt reads
    depths = merged_depths(ws)
    assert [p for p, d in depths.items() if d == 2] == [1, 2, 3, 28, 29, 30]


def test_numt_preset_needs_a_nuclear_genome(cluster, caplog):
    ws, tmp_path = cluster
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath=str(tmp_path / "mt.fa")))
    assert run_preset(ws, tmp_path, "mt-variants-numt", ["wget", "align-combined"]) == 1
    assert "NUMT-aware" in caplog.text and "refs fetch-genome Homo_sapiens --nuclear-of" in caplog.text


NUMT_PARSER = ["wget", "align-mt", "pileup-mtref", "consensus", "index-consensus", "mt-sam", "align-numt", "numt-parser", "cymt-reads", "mt-reads", "index-mt", "mpileup"]


@pytest.mark.parametrize("shell", ["bash"] + (["tcsh"] if __import__("shutil").which("tcsh") else []))
def test_numt_parser_preset_finds_numts_in_the_nuclear_genome(cluster, shell):
    ws, tmp_path = cluster
    (tmp_path / "nuc.fa").write_text(">chr1\nAAAAAAAAAA\n")
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath=str(tmp_path / "mt.fa"), nuc_refpath=str(tmp_path / "nuc.fa")))
    assert run_preset(ws, tmp_path, "mt-variants-numtparser", NUMT_PARSER + CIRCULAR, shell) == 0
    states = {j.name: j.state for j in mon.collect(ws, ["exp"], mon.ErrorScanner(), lambda s, i: None)}
    assert all(s == mon.DONE for s in states.values()), states
    derived = ws.references_dir / "derived" / "Homo_sapiens"
    assert (derived / "mt.numt_search.regions.txt").read_text() == "chr1:1-10\n"   # hit at 3..62, flanked, clipped
    assert (derived / "mt.numts.fa.bwt").exists() and (derived / "mt.with_nuclear.fa.bwt").exists()
    assert not (tmp_path / "nuc.fa.bwt").exists()        # nothing is written next to the user's genome
    out = ws.output_dir / "exp"
    # intermediates marked temporary are gone once the job finished; what the results need is kept
    assert not list(out.glob("numt/*.sam")) and not list(out.glob("alignment/*.mtref.sorted.bam"))
    assert list(out.glob("alignment/*.mt.bam")) and list(out.glob("numt/*.numt_parser.tsv"))
    numt = ws.output_dir / "exp" / "numt"
    assert (numt / "a_0_HS.cymt_reads.txt").read_text() == "r1\n"                 # both mates of r3 must be cymt
    script = (ws.output_dir / "exp" / "scripts" / "s_a_0.sh").read_text()
    # haplotype-aware: Numt Parser compares the reads with the sample's own mtDNA, which is kept
    assert (f"numt_parser.py --mt-fasta {numt}/a_0_HS.consensus.fa --numt-fasta {derived}/mt.numts.fa "
            f"--mt-sam {numt}/a_0_HS.consensus.sam") in script
    assert (numt / "a_0_HS.consensus.fa").read_text().startswith(">")
    assert f"-N {numt}/a_0_HS.cymt_reads.txt" in script and "-h -F 0x904" in script
    assert [p for p, d in merged_depths(ws).items() if d == 2] == [1, 2, 3, 28, 29, 30]


def test_numt_parser_preset_uses_a_given_numt_reference(cluster, caplog):
    ws, tmp_path = cluster
    registry = ReferenceRegistry(ws)
    registry.add("Homo_sapiens", ReferenceSet(refpath=str(tmp_path / "mt.fa")))
    assert run_preset(ws, tmp_path, "mt-variants-numtparser", ["wget", "align-numt"]) == 1
    assert "needs a NUMT reference" in caplog.text
    (tmp_path / "numts.fa").write_text(">numt1\nACGTACGT\n")
    registry.add("Homo_sapiens", ReferenceSet(refpath=str(tmp_path / "mt.fa"), numt_refpath=str(tmp_path / "numts.fa")))
    assert run_preset(ws, tmp_path, "mt-variants-numtparser", NUMT_PARSER) == 0
    prep = next((ws.output_dir / "exp" / "scripts").glob("prep_*")).read_text()
    assert "numt-regions" not in prep and f'bwa index "{tmp_path}/numts.fa"' in prep


def test_bundled_presets_are_written_once(tmp_path):
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    assert {"mt-variants", "mt-variants-numt", "mt-variants-numtparser", "rnaseq-expression"} <= set(ws.list_presets())
    wf = ws.load_workflow("rnaseq-expression")
    assert wf.active_order() == ["wget", "fastqc_pre", "trim", "fastqc", "STAR-alignment", "featureCounts", "multiqc"]
    wf.set_general("qsub_queue", "mine")
    ws.write_workflow(wf, "rnaseq-expression")
    Workspace(tmp_path / "ws").initialise()                                  # never overwritten
    assert ws.load_workflow("rnaseq-expression").general("qsub_queue") == "mine"


def test_rnaseq_expression_commands(tmp_path):
    from mitopiper.builder import build_sample_job
    from mitopiper.metadata import Run, Sample

    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("rnaseq-expression")
    ref = ReferenceSet(refpath="/refs/GRCh38/genome.fa", gtfpath="/refs/GRCh38/genes.gtf")
    pe = Sample(0, "S1", "Homo_sapiens", [Run("r", ["f/x_1.fastq.gz", "f/x_2.fastq.gz"])])
    se = Sample(1, "S2", "Homo_sapiens", [Run("r", ["f/y.fastq.gz"])])
    steps = {s.step: s.lines for s in build_sample_job(pe, wf, str(tmp_path / "ws/output/E"), ref).steps}
    (star,) = steps["STAR-alignment"]
    assert "--genomeDir /refs/GRCh38" in star and "--runThreadN 8" in star and "--readFilesCommand gunzip -c" in star
    (counts,) = steps["featureCounts"]
    assert counts.startswith("featureCounts -T 4 -a /refs/GRCh38/genes.gtf -t exon -g gene_id -p --countReadPairs -o ")
    assert counts.endswith("S1_0_HSAligned.sortedByCoord.out.bam")
    se_counts = {s.step: s.lines for s in build_sample_job(se, wf, str(tmp_path / "ws/output/E"), ref).steps}["featureCounts"][0]
    assert "--countReadPairs" not in se_counts


def test_references_you_add_get_their_missing_indices(tmp_path):
    from mitopiper.references import prepare_commands

    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("mt-variants")
    fasta = tmp_path / "mt.fa"
    fasta.write_text(">m\nACGT\n")
    (tmp_path / "mt.fa.fai").write_text("m\t4\t3\t4\t5\n")                       # already there
    commands = prepare_commands(ReferenceSet(refpath=str(fasta)), "Homo_sapiens", str(ws.root), wf, csh=False,
                                steps=["bwa-align"])
    assert commands == [f'if [ ! -s "{fasta}.bwt" ]; then bwa index "{fasta}"; fi']
    elsewhere = ReferenceSet(refpath="/cluster/only/mt.fa")                       # not on this machine: left to you
    assert prepare_commands(elsewhere, "Homo_sapiens", str(ws.root), wf, csh=False, steps=["bwa-align"]) == []


def test_present_star_index_does_not_read_the_fasta(tmp_path, monkeypatch):
    from mitopiper import references
    from mitopiper.references import prepare_commands

    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("rnaseq-expression")
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">c\nACGT\n")
    (tmp_path / "SA").write_text("index")
    monkeypatch.setattr(references, "validate_fasta", lambda path: pytest.fail("the FASTA was read"))
    assert prepare_commands(ReferenceSet(refpath=str(fasta)), "Homo_sapiens", str(ws.root), wf, csh=False) == []
    (tmp_path / "SA").unlink()                                                   # now it is needed: STAR's command
    monkeypatch.setattr(references, "validate_fasta", lambda path: 4)
    (command,) = prepare_commands(ReferenceSet(refpath=str(fasta)), "Homo_sapiens", str(ws.root), wf, csh=False)
    assert "genomeGenerate" in command and "--genomeSAindexNbases" in command


SLOW_BWA = r"""#!/bin/sh
# bwa index that takes a moment and notices another indexer running at the same time
[ "$1" = index ] || exit 0
if [ -e "$2.indexing" ]; then echo COLLISION >> "$2.log"; fi
touch "$2.indexing"; echo start >> "$2.log"; sleep 1; echo index > "$2.bwt"; rm -f "$2.indexing"
"""


def test_reference_jobs_of_parallel_builds_take_turns(tmp_path, monkeypatch):
    from mitopiper.references import prepare_job_lines

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, body in {"bwa": SLOW_BWA, "samtools": FAKE_SAMTOOLS,
                       "mitopiper": f'#!/bin/sh\nPYTHONPATH="{REPO}" exec "{sys.executable}" -m mitopiper "$@"\n'}.items():
        (bin_dir / name).write_text(body)
        (bin_dir / name).chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    monkeypatch.setenv("MITOPIPER_LOCK_WAIT", "1")          # look at the lock every second, not every minute
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    (tmp_path / "mt.fa").write_text(f">NC_012920.1\n{MT}\n")
    ref = ReferenceSet(refpath=str(tmp_path / "mt.fa"))
    wf = ws.load_workflow("mt-variants")
    wf.set_general("mt_rotation", 5)
    lines = prepare_job_lines(ref, "Homo_sapiens", str(ws.root), wf, steps=["bwa-align", "align-rotated"])
    assert lines[-1].startswith("sh -c ") and "\n" not in lines[-1]
    shells = ["bash"] + [s for s in ("tcsh", "csh") if shutil.which(s)][:1]     # the job shell may be csh
    scripts = []
    for i, shell in enumerate(shells + ["bash"]):
        (tmp_path / f"prep{i}.sh").write_text("\n".join(lines) + "\n")
        scripts.append(subprocess.Popen([shell, str(tmp_path / f"prep{i}.sh")], stdout=subprocess.DEVNULL))
    assert [p.wait() for p in scripts] == [0] * len(scripts)
    assert not (ws.references_dir / "derived" / "Homo_sapiens" / ".prepare.lock.d").exists()      # released
    rotated = ws.references_dir / "derived" / "Homo_sapiens" / "mt.rotated5.fa"
    for log in (tmp_path / "mt.fa.log", rotated.with_name(rotated.name + ".log")):
        assert log.read_text().split() == ["start"], log      # indexed once, never two at the same time


def test_a_stale_reference_lock_is_taken_over(tmp_path, monkeypatch):
    # A reference job killed by the scheduler (SIGKILL: no trap) leaves its lock behind.
    from mitopiper.references import prepare_job_lines

    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    (tmp_path / "mt.fa").write_text(f">NC_012920.1\n{MT}\n")
    wf = ws.load_workflow("mt-variants")
    wf.set_general("mt_rotation", 5)
    lines = prepare_job_lines(ReferenceSet(refpath=str(tmp_path / "mt.fa")), "Homo_sapiens", str(ws.root), wf,
                              steps=["align-rotated"])
    lock = ws.references_dir / "derived" / "Homo_sapiens" / ".prepare.lock.d"
    lock.mkdir(parents=True)
    old = time.time() - 20 * 60
    os.utime(lock, (old, old))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, body in {"bwa": "#!/bin/sh\n[ \"$1\" = index ] && echo x > \"$2.bwt\"\nexit 0\n", "samtools": FAKE_SAMTOOLS,
                       "mitopiper": f'#!/bin/sh\nPYTHONPATH="{REPO}" exec "{sys.executable}" -m mitopiper "$@"\n'}.items():
        (bin_dir / name).write_text(body)
        (bin_dir / name).chmod(0o755)
    env = dict(os.environ, PATH=f"{bin_dir}:{os.environ['PATH']}", MITOPIPER_LOCK_WAIT="1")
    (tmp_path / "prep.sh").write_text("\n".join(lines) + "\n")
    out = subprocess.run(["bash", str(tmp_path / "prep.sh")], env=env, capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, out.stderr
    assert "taking over the stale lock" in out.stderr and not lock.exists()


def test_temporary_outputs_are_removed_last_and_only_on_success():
    from mitopiper.builder import TEMPORARY_STEP, build_sample_job
    from mitopiper.metadata import Run, Sample

    from mitopiper.defaults import bundled_presets
    from mitopiper.config import Workflow
    preset = Workflow(**bundled_presets()["mt-variants"])
    sample = Sample(index=0, accession="S", organism="Homo_sapiens", runs=[Run("R", ["/f_1.fq.gz", "/f_2.fq.gz"])])
    job = build_sample_job(sample, preset, "/w/output/E", ReferenceSet(refpath="/r/mt.fa"))
    last = job.steps[-1]
    assert last.step == TEMPORARY_STEP and [s.step for s in job.steps].index("annotate") < len(job.steps) - 1
    assert last.lines == ['rm -f "/w/output/E/alignment/S_0_HS.sam"', 'rm -f "/w/output/E/alignment/S_0_HS.bam"']
    # a rerun resumed late still removes what the first attempt left
    late = build_sample_job(sample, preset, "/w/output/E", ReferenceSet(refpath="/r/mt.fa"), steps=["annotate"])
    assert late.steps[-1].lines == last.lines
    preset.set_general("fail_fast", False)               # the job goes on after errors: delete nothing
    assert build_sample_job(sample, preset, "/w/output/E", ReferenceSet(refpath="/r/mt.fa")).steps[-1].step == "annotate"


def test_existing_numt_reference_is_reused_without_the_genome_index(tmp_path):
    # After a run, the 3 GB nuclear + mtDNA copy and its index can be deleted; the NUMT reference is kept.
    from mitopiper.references import derived_commands, derived_paths

    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    ref = ReferenceSet(refpath=str(tmp_path / "mt.fa"), nuc_refpath=str(tmp_path / "nuc.fa"))
    wf = ws.load_workflow("mt-variants-numtparser")
    before = " ".join(derived_commands(ref, "Homo_sapiens", str(ws.root), wf))
    assert "combine-references" in before and "numt-regions" in before
    numts = derived_paths(ref, "Homo_sapiens", str(ws.root), "half")["numt"]
    Path(numts).parent.mkdir(parents=True)
    Path(numts).write_text(">chr1:1-10\nACGT\n")
    after = " ".join(derived_commands(ref, "Homo_sapiens", str(ws.root), wf))
    assert "combine-references" not in after and "numt-regions" not in after and "numts.fa" in after
