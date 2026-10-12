from mitopiper.builder import (
    build_experiment_steps,
    build_sample_job,
    guess_threads,
    render_settings,
    render_structure,
    selected_steps,
)
from mitopiper.config import Workflow
from mitopiper.metadata import Run, Sample
from mitopiper.references import ReferenceSet

REF = ReferenceSet(refpath="/refs/chrM.fa", gtfpath="/refs/chrM.gff3")
EXP = "/data/output/EXP"


def se_sample(**kw):
    return Sample(index=3, accession="SRR1", organism="Homo_sapiens", runs=[Run("SRR1", ["ftp/SRR1.fastq.gz"], ["abc"])], **kw)


def pe_sample():
    return Sample(index=0, accession="SAMN1", organism="Danio_rerio",
                  runs=[Run("SRR1", ["ftp/SRR1_1.fastq.gz", "ftp/SRR1_2.fastq.gz"]),
                        Run("SRR2", ["ftp/SRR2_1.fastq.gz", "ftp/SRR2_2.fastq.gz"])])


def lines(job):
    return {s.step: s.lines for s in job.steps}


def test_single_end_default_pipeline():
    job = build_sample_job(se_sample(), Workflow.default(), EXP, REF)
    got = lines(job)
    assert job.job_name == "s_SRR1_3"
    assert got["wget"][0] == "wget -c -O /data/output/EXP/fastq/SRR1_3_HS.fastq.gz ftp/SRR1.fastq.gz"
    assert got["wget"][1].startswith("sh -c ") and 'echo "abc  /data/output/EXP/fastq/SRR1_3_HS.fastq.gz" | md5sum -c -' \
        in got["wget"][1]
    assert got["wget"][2:] == ["sleep 30"]
    assert got["trim"] == ["trim_galore -q 20 -o /data/output/EXP/post_trim /data/output/EXP/fastq/SRR1_3_HS.fastq.gz"]
    assert got["bwa-align"] == [
        "bwa mem -k 19 /refs/chrM.fa /data/output/EXP/post_trim/SRR1_3_HS_trimmed.fq.gz > /data/output/EXP/alignment/SRR1_3_HS.sam"
    ]
    assert got["view"] == ["samtools view -h -b /data/output/EXP/alignment/SRR1_3_HS.sam > /data/output/EXP/alignment/SRR1_3_HS.bam"]
    assert got["sort"] == [
        "samtools sort -o /data/output/EXP/alignment/SRR1_3_HS.sorted.bam /data/output/EXP/alignment/SRR1_3_HS.bam"
    ]
    assert got["index"] == ["samtools index /data/output/EXP/alignment/SRR1_3_HS.sorted.bam"]
    assert got["mpileup"] == [
        "samtools mpileup -Q 30 -f /refs/chrM.fa /data/output/EXP/alignment/SRR1_3_HS.sorted.bam > /data/output/EXP/pileup/SRR1_3_HS.pileup"
    ]


def test_paired_end_multi_run_concatenates():
    job = build_sample_job(pe_sample(), Workflow.default(), EXP, REF)
    wget = lines(job)["wget"]
    assert 'cat "/data/output/EXP/fastq/SAMN1_0_DR_1.part0.fastq.gz" "/data/output/EXP/fastq/SAMN1_0_DR_1.part1.fastq.gz" > "/data/output/EXP/fastq/SAMN1_0_DR_1.fastq.gz"' in wget
    assert any(l.endswith("ftp/SRR2_2.fastq.gz") for l in wget)
    assert not any("md5sum" in l for l in wget)  # no md5 known
    trim = lines(job)["trim"][0]
    assert "--paired /data/output/EXP/fastq/SAMN1_0_DR_1.fastq.gz /data/output/EXP/fastq/SAMN1_0_DR_2.fastq.gz" in trim
    assert len(lines(job)["fastqc"]) == 2


def test_local_files_are_linked_or_gzipped():
    sample = Sample(index=0, accession="S1", organism="Homo_sapiens",
                    runs=[Run("r1", ["/in/a_R1.fastq.gz", "/in/a_R2.fastq"])])
    wget = lines(build_sample_job(sample, Workflow.default(), EXP, REF, steps=["wget"]))["wget"]
    assert wget == [
        'ln -sf "/in/a_R1.fastq.gz" "/data/output/EXP/fastq/S1_0_HS_1.fastq.gz"',
        'gzip -c "/in/a_R2.fastq" > "/data/output/EXP/fastq/S1_0_HS_2.fastq.gz"',
    ]


def test_if_se_if_pe_flags():
    params = {
        "se": {"key": "--single", "value": "_if_SE", "switch": True, "order": 0},
        "pe": {"key": "--paired", "value": "_if_PE", "switch": True, "order": 1},
    }
    assert render_settings(params, {"sample": "x"}, {}, paired=True, output_ext="") == "--paired"
    assert render_settings(params, {"sample": "x"}, {}, paired=False, output_ext="") == "--single"


def test_settings_order_switch_and_placeholders():
    params = {
        "b": {"key": "-b", "value": "_gtf", "switch": True, "order": 1},
        "a": {"key": "-a", "value": "{sample}.txt", "switch": True, "order": 0},
        "off": {"key": "-z", "value": "", "switch": False, "order": 2},
        "o": {"key": "-o", "value": "dir_trim", "switch": 1, "order": 3},
        "out": {"key": "--output_file", "value": "_fastq_dir", "switch": True, "order": 4},
    }
    rendered = render_settings(
        {k: v for k, v in params.items() if k != "out"} | {"output_file": params["out"]},
        {"sample": "S", "gtf": "/g.gff", "fastq_dir": "/fq"}, {"dir_trim": "/exp/post_trim"}, False, ".txt",
    )
    assert rendered == "-a S.txt -b /g.gff -o /exp/post_trim --output_file /fq/S.txt"


def test_structure_is_not_evaluated_as_python():
    out = render_structure("key settings __import__('os') -t 1,2 {input} output_pipe",
                           {"command_key": "tool", "command_settings": "-x", "command_input": "/in",
                            "output_pipe": "> /out"}, {})
    assert out == "tool -x __import__('os') -t 1,2 /in > /out"


def test_steps_selection_and_skip():
    wf = Workflow.default()
    assert selected_steps(wf, ["trim", "fastqc"]) == ["trim", "fastqc"]
    assert "wget" not in selected_steps(wf, skip=["wget"])


def test_experiment_scope_steps():
    wf = Workflow.default()
    wf.pipeline["multiqc"]["switch"] = True
    wf.run_order.append("multiqc")
    samples = [se_sample(), pe_sample()]
    job = build_sample_job(samples[0], wf, EXP, REF)
    assert "multiqc" not in lines(job)
    (step,) = build_experiment_steps(samples, wf, EXP)
    assert step.lines == ["multiqc -o /data/output/EXP/multiqc /data/output/EXP"]


def test_guess_threads():
    wf = Workflow.default()
    assert guess_threads(wf) == 4  # nothing active declares threads
    wf.pipeline["STAR-alignment"]["switch"] = True
    wf.run_order.append("STAR-alignment")
    wf.settings["STAR-alignment"]["runThreadN"]["value"] = 12
    assert guess_threads(wf) == 12


def test_output_ext_forms():
    from mitopiper.builder import output_exts

    assert output_exts({"output_ext": ".bam"}) == {"0": ".bam"}
    assert output_exts({"output_ext": ""}) == {"0": ""}
    assert output_exts({"output_ext": ".bam;.bam.bai"}) == {"0": ".bam", "1": ".bam.bai"}
    assert output_exts({"output_ext": "bam=.bam; log=.log"}) == {"bam": ".bam", "log": ".log"}
    assert output_exts({"output_ext": [".a", ".b"]}) == {"0": ".a", "1": ".b"}


def test_multiple_outputs_are_referencable():
    wf = Workflow.default()
    for step in ("STAR-alignment", "htseq-count"):
        wf.pipeline[step]["switch"] = True
    # htseq counts STAR's sorted BAM and writes the count table plus a SAM of assignments
    wf.pipeline["htseq-count"]["structure"] = "key settings -o output_htseq-count:assign output_STAR-alignment:bam gtf > output_htseq-count"
    wf.pipeline["htseq-count"]["output_ext"] = {"counts": ".counts.txt", "assign": ".assign.sam"}
    wf.run_order = ["STAR-alignment", "htseq-count", "multiqc"]
    wf.pipeline["multiqc"]["switch"] = True
    wf.pipeline["multiqc"]["structure"] = "echo outputs_STAR-alignment:log"
    sample = se_sample()
    got = lines(build_sample_job(sample, wf, EXP, REF))
    assert got["htseq-count"][0].endswith(
        "-o /data/output/EXP/count/SRR1_3_HS.assign.sam "
        "/data/output/EXP/alignment/SRR1_3_HS" "Aligned.sortedByCoord.out.bam /refs/chrM.gff3 "
        "> /data/output/EXP/count/SRR1_3_HS.counts.txt"
    )
    assert "--outFileNamePrefix /data/output/EXP/alignment/SRR1_3_HS" in got["STAR-alignment"][0]
    (exp,) = build_experiment_steps([sample, pe_sample()], wf, EXP)
    assert exp.lines == [
        "echo /data/output/EXP/alignment/SRR1_3_HSLog.final.out /data/output/EXP/alignment/SAMN1_0_DRLog.final.out"
    ]


def test_output_all_and_embedded():
    from mitopiper.builder import step_io, substitute

    wf = Workflow.default()
    wf.pipeline["index"]["output_ext"] = ".bai;.csi"
    outputs, _ = step_io(wf, EXP, "S")
    assert outputs["output_index"] == "/data/output/EXP/alignment/S.bai"
    assert outputs["output_index:all"] == "/data/output/EXP/alignment/S.bai /data/output/EXP/alignment/S.csi"
    assert substitute("x={output_index:1}", outputs) == "x=/data/output/EXP/alignment/S.csi"


def test_per_step_resources():
    from mitopiper.build import job_resources

    wf = Workflow.default()
    wf.set_general("qsub_threads", 2)
    wf.set_general("qsub_memory", "8G")
    wf.set_general("walltime", "02:00:00")
    base = job_resources(wf, ["trim", "fastqc"])
    assert (base.threads, base.memory, base.walltime) == (2, "8G", "02:00:00")
    wf.pipeline["trim"]["walltime"] = "01:30:00"
    wf.pipeline["bwa-align"].update(threads="6", memory="12G", walltime="01:00:00")
    res = job_resources(wf, ["trim", "bwa-align"])
    assert (res.threads, res.memory, res.walltime) == (6, "12G", "02:30:00")
    wf.pipeline["view"]["memory"] = "1500M"                      # smaller than the general setting: ignored
    assert job_resources(wf, ["view"]).memory == "8G"


def test_star_uses_the_step_threads():
    wf = Workflow.default()
    wf.pipeline["STAR-alignment"]["switch"] = True
    wf.run_order = ["STAR-alignment"]
    (cmd,) = lines(build_sample_job(se_sample(), wf, EXP, REF))["STAR-alignment"]
    assert "--runThreadN 8" in cmd


def test_a_corrupt_download_is_downloaded_once_more(tmp_path):
    # Seen on the BGU cluster: full size, wrong MD5; 'wget -c' then kept the corrupt file.
    import hashlib
    import shutil
    import subprocess

    good = b"good data\n"
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "wget").write_text("#!/bin/sh\n# first download corrupt, then right\n"
                                  'while [ "$1" != -O ]; do shift; done; out=$2\n'
                                  f'if [ -e "{tmp_path}/once" ]; then printf "good data\\n" > "$out"; '
                                  f'else touch "{tmp_path}/once"; printf "bad! data\\n" > "$out"; fi\n')
    (bin_dir / "wget").chmod(0o755)
    sample = Sample(index=0, accession="SAMN1", organism="Danio_rerio",
                    runs=[Run("SRR1", ["ftp/SRR1.fastq.gz"], md5=[hashlib.md5(good).hexdigest()])])
    wf = Workflow.default()
    wf.set_general("download_delay", 0)
    for shell in ["bash"] + [s for s in ("tcsh",) if shutil.which(s)]:
        (tmp_path / "once").unlink(missing_ok=True)
        job = build_sample_job(sample, wf, str(tmp_path / "exp"), REF, steps=["wget"])
        (tmp_path / "exp" / "fastq").mkdir(parents=True, exist_ok=True)
        script = tmp_path / "dl.sh"
        script.write_text("\n".join(lines(job)["wget"]) + "\n")
        env = dict(__import__("os").environ, PATH=f"{bin_dir}:" + __import__("os").environ["PATH"])
        out = subprocess.run([shell, str(script)], env=env, capture_output=True, text=True)
        assert out.returncode == 0 and "downloading" in out.stderr, (shell, out.stderr)
        assert (tmp_path / "exp" / "fastq" / "SAMN1_0_DR.fastq.gz").read_bytes() == good
