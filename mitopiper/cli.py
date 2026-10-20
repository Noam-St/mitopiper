"""Command line interface.

Run ``mitopiper`` without arguments for the interactive menus, or use the
sub-commands for scripted, reproducible builds::

    mitopiper build PRJNA12345.tsv --preset rnaseq --scheduler slurm --samples all
    mitopiper build PRJEB12345 --scheduler slurm     # metadata fetched from ENA
    mitopiper build PRJNA12345.tsv --dry-run --steps wget,trim,fastqc
    mitopiper refs fetch-mt Danio_rerio
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from pathlib import Path
from typing import Optional

from . import __version__, build_info
from .build import BuildError, BuildOptions, build_experiment, resolve_metadata_source, write_master
from .cleanup import cleanup
from .config import ConfigError, Workspace
from .metadata import MetadataError, load_samples, parse_selection
from .references import ReferenceRegistry, RefError, ReferenceSet, fetch_mt_reference
from .schedulers import SCHEDULERS, normalise
from .ui import table

log = logging.getLogger("mitopiper")


def _csv_list(value: Optional[str]) -> Optional[list]:
    if not value:
        return None
    return [v.strip() for v in value.replace(" ", ",").split(",") if v.strip()]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mitopiper", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__} ({build_info()})")
    parser.add_argument("-w", "--workdir", default=".", help="workspace directory (default: current directory)")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="more logging (-vv for debug)")
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("interactive", help="menu driven interface (default)")
    sub.add_parser("init", help="create the workspace skeleton and default preset")

    b = sub.add_parser("build", help="generate job scripts non-interactively")
    b.add_argument("metadata", nargs="*", help="metadata file(s), previous experiment name(s), ENA/SRA accession(s) such as "
                                               "PRJEB12345, or folder(s) of FASTQ files (with --org)")
    b.add_argument("-p", "--preset", default="default", help="workflow preset to use (default: %(default)s)")
    b.add_argument("-s", "--scheduler", default="sge", help=f"one of {', '.join(SCHEDULERS)} (default: %(default)s)")
    b.add_argument("--samples", default="all", help="rows to process: all, 3, 1-4 (inclusive), '0 2 5' (default: all)")
    b.add_argument("--name", help="experiment name (default: metadata file name); only with a single metadata file")
    b.add_argument("--steps", help="comma separated steps to run instead of the preset's run order")
    b.add_argument("--skip", help="comma separated steps to leave out")
    b.add_argument("--reference", help="name of the reference set to use for every organism that has it")
    fq = b.add_argument_group("FASTQ folders", "a folder of .fastq/.fq(.gz) files is brought into the workspace "
                              "(fastq/<folder name>/) and described in metadata/<folder name>.tsv, as make-tsv does")
    fq.add_argument("--org", help="scientific name of the organism of the FASTQ files, e.g. 'Homo sapiens'")
    fq.add_argument("--lib-strat", default="", help="library strategy of the FASTQ files, e.g. RNA-Seq")
    fq.add_argument("--fastq-mode", choices=("link", "copy", "move"), default="link",
                    help="link the files (default; they stay where they are), copy them, or move them into the workspace")
    fq.add_argument("--browse", action="store_true", help="choose the FASTQ folder in a file browser (needs a display)")
    b.add_argument("--fetch-mt", action="store_true", help="download missing references (mtDNA) from NCBI instead of failing")
    b.add_argument("--dry-run", action="store_true", help="resolve everything and write plan.json, but no scripts")
    b.add_argument("--refresh", action="store_true", help="download the metadata of an accession again even if it was used before")
    b.add_argument("--auto-fix", action="store_true",
                   help="if tools the jobs call are missing here, try to install them into the current environment")
    b.add_argument("--no-check", action="store_true",
                   help="do not check here that the tools exist (e.g. when generating scripts on another machine)")

    ck = sub.add_parser("check", help="check that every tool the jobs call can be found (optionally install them)")
    ck.add_argument("experiments", nargs="*", help="experiments (default: those of the last build)")
    ck.add_argument("--auto-fix", action="store_true", help="try to install missing tools into the current environment")

    f = sub.add_parser("fetch-metadata", help="download the ENA run report of a project/study/sample accession")
    f.add_argument("accession", help="e.g. PRJEB12345, PRJNA67890, SRP012345, SAMN01234567")
    f.add_argument("-o", "--output", help="output TSV (default: metadata/<ACCESSION>.tsv in the workspace)")

    p = sub.add_parser("presets", help="list, show, change or save workflow presets")
    psub = p.add_subparsers(dest="action", required=True)
    psub.add_parser("list")
    ps = psub.add_parser("show")
    ps.add_argument("name", nargs="?", default="default")
    ps.add_argument("--json", action="store_true", help="dump the full preset as JSON")
    pset = psub.add_parser("set", help="change general settings or step fields of presets",
                           description="Examples: presets set mt-variants qsub_queue=bioinfo.q walltime=02:00:00; "
                                       "presets set mt-variants,mt-variants-numt numt-parser.memory=32G mpileup.path=/opt/bin/samtools")
    pset.add_argument("names", help="preset, or several separated by commas")
    pset.add_argument("assignments", nargs="+", metavar="KEY=VALUE",
                      help="a general setting (qsub_memory=8G), or STEP.FIELD=VALUE for a step (trim.threads=4)")
    preset_reset = psub.add_parser("reset", help="restore bundled presets to the version of this mitopiper (after an "
                                                 "upgrade); the old copy is kept as NAME_old_DATE")
    preset_reset.add_argument("names", nargs="+", help="bundled presets, e.g. mt-variants-numtparser")
    psave = psub.add_parser("save", help="copy a preset (default: last_used) to a new name")
    psave.add_argument("name")
    psave.add_argument("--from", dest="source", default="last_used")

    r = sub.add_parser("refs", help="manage reference genomes")
    rsub = r.add_subparsers(dest="action", required=True)
    rsub.add_parser("list")
    ra = rsub.add_parser("add", help="register an existing reference")
    ra.add_argument("organism")
    ra.add_argument("--fasta", required=True)
    ra.add_argument("--gtf", default="")
    ra.add_argument("--nuc-fasta", default="")
    ra.add_argument("--numt-fasta", default="", help="NUMT sequences for Numt Parser (default: found in --nuc-fasta)")
    ra.add_argument("--chrominfo", default="")
    ra.add_argument("--name", default="default")
    ra.add_argument("--no-activate", action="store_true")
    ru = rsub.add_parser("use", help="switch the active reference of an organism")
    ru.add_argument("organism")
    ru.add_argument("name")
    rr = rsub.add_parser("remove")
    rr.add_argument("organism")
    rr.add_argument("name")
    rf = rsub.add_parser("fetch-mt", help="download a mitochondrial genome from NCBI")
    rf.add_argument("organism")
    rf.add_argument("--accession", help="RefSeq accession (default: search NCBI)")
    rf.add_argument("--name", default="mt")
    rf.add_argument("--no-activate", action="store_true")
    rg = rsub.add_parser("fetch-genome", help="plan the download of the nuclear genome from NCBI Genome (large!)")
    rg.add_argument("organism")
    rg.add_argument("--accession", help="assembly accession, e.g. GCF_000002035.6 (default: the reference genome)")
    rg.add_argument("--name", default="genome", help="name of the new reference set (default: %(default)s)")
    rg.add_argument("--nuclear-of", metavar="SET",
                    help="instead of a new set, use the genome as the nuclear FASTA (_fasta_nuc) of reference set SET")
    rg.add_argument("--no-activate", action="store_true")
    rg.add_argument("-y", "--yes", action="store_true", help="do not ask for confirmation after the size warning")

    c = sub.add_parser("cleanup", help="remove temporary files (logs, scripts) from experiments")
    c.add_argument("folders", nargs="*", default=["stdout", "stderr"], help="folders to empty (default: stdout stderr)")
    c.add_argument("-e", "--experiment", action="append", help="limit to this experiment (repeatable)")
    c.add_argument("-y", "--yes", action="store_true", help="do not ask for confirmation")
    c.add_argument("--keep-scripts", action="store_true", help="keep generated job scripts")

    mo = sub.add_parser("monitor", help="follow submitted jobs and show errors from their stderr")
    mo.add_argument("experiments", nargs="*", help="experiments to follow (default: those of the last build)")
    mo.add_argument("-w", "--watch", nargs="?", type=float, const=30.0, metavar="SECONDS",
                    help="refresh every SECONDS (default 30) until all jobs have finished")
    mo.add_argument("-j", "--job", help="details of one job (name or scheduler id): steps, error lines, stderr tail")
    mo.add_argument("-e", "--errors", action="store_true", help="details of every failed job and every job with errors in stderr")
    mo.add_argument("--tail", type=int, default=20, help="stderr lines to show in details (default: %(default)s)")
    mo.add_argument("--pattern", action="append", help="extra regular expression that counts as an error (repeatable)")
    mo.add_argument("--ignore", action="append", help="regular expression for stderr lines that are not errors (repeatable)")
    mo.add_argument("--json", action="store_true", help="machine readable output")
    mo.add_argument("--no-color", action="store_true")
    mo.add_argument("-i", "--interactive", action="store_true",
                    help="interactive monitor: live table plus commands to inspect, diagnose and rerun jobs")
    mo.add_argument("--rerun", action="store_true",
                    help="shortcut: resubmit failed jobs from their failed step with the automatic fixes (= rerun --auto-fix --submit)")

    rr = sub.add_parser("rerun", help="resubmit failed jobs, resuming at the step that failed")
    rr.add_argument("experiments", nargs="*", help="experiments (default: those of the last build)")
    rr.add_argument("-j", "--job", action="append", help="only this job (name or id; repeatable), whatever its state")
    rr.add_argument("--from", dest="from_step", metavar="STEP", help="resume at STEP instead of the step that failed")
    rr.add_argument("--from-start", action="store_true", help="run the whole job again")
    rr.add_argument("--clean", action="store_true", help="delete the outputs of the step the job resumes at first")
    rr.add_argument("--auto-fix", action="store_true",
                    help="apply the suggested automatic fixes per job (more memory/time, clean corrupt input, earlier resume point)")
    rr.add_argument("--memory", help="memory per job for the rerun, e.g. 64G")
    rr.add_argument("--walltime", help="time limit for the rerun, e.g. 3-00:00:00")
    rr.add_argument("--threads", type=int, help="threads per job for the rerun")
    rr.add_argument("--queue", help="queue/partition for the rerun")
    rr.add_argument("--free-form", action="append", default=[], metavar="LINE",
                    help="extra line for the rerun jobs, e.g. 'module load samtools' (repeatable)")
    rr.add_argument("-p", "--preset", help="use this (fixed) preset instead of the workflow recorded at build time")
    rr.add_argument("-s", "--scheduler", help="submit with another scheduler (e.g. local)")
    rr.add_argument("--include-pending", action="store_true", help="also resubmit jobs that are still pending (they are cancelled first)")
    rr.add_argument("--dry-run", action="store_true", help="only show what would be rerun")
    rr.add_argument("--submit", action="store_true", help="run rerun.sh right away")

    from . import download, mtdna

    mtdna.add_parser(sub)
    download.add_parser(sub)

    m = sub.add_parser("make-tsv", help="create a metadata table for local FASTQ files")
    m.add_argument("folder", help="folder containing the FASTQ files")
    m.add_argument("--org", required=True, help="scientific name of the organism")
    m.add_argument("-o", "--output", help="output TSV (default: <folder name>.tsv)")
    m.add_argument("-l", "--lib-strat", default="", help="library strategy, e.g. RNA-Seq")
    m.add_argument("--ignore-paired", action="store_true")
    m.add_argument("--force-paired", action="store_true")
    return parser


def _metadata_items(ws: Workspace, args) -> list:
    """The metadata tables of the build: files, experiments and accessions as given, FASTQ folders
    (and the one chosen with --browse) imported into the workspace first."""
    from .local_samples import browse_for_folder, import_fastq_folder

    items = list(args.metadata)
    if args.browse:
        try:
            folder = browse_for_folder()
        except RuntimeError as exc:
            raise BuildError(f"{exc}; give the folder's path instead") from None
        if not folder:
            raise BuildError("no folder was chosen")
        print(f"FASTQ folder: {folder}")
        items.append(folder)
    if not items:
        raise BuildError("give a metadata file, an experiment, an accession or a FASTQ folder (or --browse)")
    tables = []
    for item in items:
        if Path(item).expanduser().is_dir():
            if not args.org:
                raise BuildError(f"{item} is a folder of FASTQ files: give their organism with --org 'Genus species'")
            tables.append(import_fastq_folder(ws.root, item, args.org, args.fastq_mode, library_strategy=args.lib_strat,
                                              dry_run=args.dry_run))
            done = {"link": "linked", "copy": "copied", "move": "moved"}[args.fastq_mode]
            target = ws.root / "fastq" / Path(item).expanduser().resolve().name
            print(f"{item}: files would be {done} into {target} (dry run)" if args.dry_run else
                  f"{item}: files {done} into {target}, sample table {tables[-1]}")
        else:
            tables.append(resolve_metadata_source(ws, item, refresh=args.refresh))
    return tables


def cmd_build(ws: Workspace, args) -> int:
    ws.initialise()
    items = _metadata_items(ws, args)
    if args.name and len(items) > 1:
        raise BuildError("--name can only be used with a single metadata file")
    workflow = ws.activate(args.preset)
    scheduler = normalise(args.scheduler)
    options = BuildOptions(
        scheduler=scheduler,
        steps=_csv_list(args.steps),
        skip=_csv_list(args.skip),
        reference=args.reference,
        dry_run=args.dry_run,
    )
    registry = ReferenceRegistry(ws)

    def resolver(organism: str) -> Optional[ReferenceSet]:
        if not args.fetch_mt:
            return None
        ref = fetch_mt_reference(ws, organism)
        registry.add(organism, ref, "mt")
        return ref

    results = []
    for metadata in items:
        ignore_paired = str(workflow.general("ignore_paired_end", False)).lower() in ("true", "1", "yes")
        samples = load_samples(metadata, ignore_paired=ignore_paired)
        selection = parse_selection(args.samples, len(samples))
        result = build_experiment(ws, workflow, metadata, selection, options, resolver, name=args.name, samples=samples)
        if result.skipped_organisms:
            log.error(
                "No reference for: %s. Register one with 'mitopiper refs add' or rerun with --fetch-mt.",
                ", ".join(result.skipped_organisms),
            )
        results.append(result)
        if args.dry_run:
            print(f"[dry-run] {result.name}: {len(result.plan['samples'])} sample(s), steps: {' -> '.join(result.plan['steps'])}")
            print(f"[dry-run] plan written to {result.directory / 'plan.json'}")
        else:
            count = len(result.jobs.scripts())
            print(f"{result.name}: {count} job script(s) in {result.directory / 'scripts'}")
    if args.dry_run:
        return 0
    if not any(r.jobs.samples or r.jobs.experiment for r in results):
        log.error("No job scripts were generated")
        return 1
    master = write_master(ws, scheduler, results, shell=str(workflow.general("environment", "bash")))
    ws.copy_preset("current", "last_used")
    from .build import tool_check_enabled

    if not args.no_check and tool_check_enabled(workflow):
        from . import toolcheck

        tools = toolcheck.merge(*[r.tools for r in results])
        result = toolcheck.run_preflight(tools, workflow, auto_fix=args.auto_fix)
        if result.failed:
            log.error("Halted: %d tool(s) the jobs call were not found here (%s). Make them available "
                      "(module load / conda activate in 'free_form'), rerun with --auto-fix to install them, "
                      "or --no-check if the jobs run in a different environment. The scripts are written; "
                      "the check job would stop the pipeline on the cluster anyway.",
                      len(result.missing), ", ".join(result.missing))
            return 1
    print(f"Submit everything with: bash {master}")
    return 0


def cmd_presets(ws: Workspace, args) -> int:
    ws.initialise()
    if args.action == "list":
        for name in ws.list_presets():
            print(name)
    elif args.action == "show":
        wf = ws.load_workflow(args.name)
        if args.json:
            print(json.dumps(wf.as_dict(), indent=2))
            return 0
        rows = []
        for step in wf.run_order + [s for s in wf.pipeline if s not in wf.run_order]:
            conf = wf.pipeline.get(step, {})
            rows.append([step, "on" if step in wf.active_order() else "off", conf.get("scope", "sample"), conf.get("key", ""), conf.get("description", "")])
        print(table(rows, ["Step", "State", "Scope", "Command", "Description"]))
        print(table([[k, v[0] if isinstance(v, list) else v] for k, v in wf.general_settings.items()], ["General setting", "Value"]))
    elif args.action == "set":
        return _set_preset_values(ws, args.names.split(","), args.assignments)
    elif args.action == "reset":
        from .defaults import bundled_presets
        from .config import Workflow

        bundled = bundled_presets()
        unknown = [n for n in args.names if n not in bundled]
        if unknown:
            log.error("Not a bundled preset: %s (bundled: %s)", ", ".join(unknown), ", ".join(sorted(bundled)))
            return 1
        stamp = __import__("datetime").date.today().strftime("%Y%m%d")
        for name in args.names:
            old = f"{name}_old_{stamp}"
            if name in ws.list_presets():
                ws.save_preset(old, name)
            ws.write_workflow(Workflow(**bundled[name]), name)
            print(f"Reset preset '{name}' to this version of mitopiper" + (f" (previous copy: '{old}')" if name in ws.list_presets() else ""))
    elif args.action == "save":
        ws.save_preset(args.name, args.source)
        print(f"Saved preset '{args.name}' from '{args.source}'")
    return 0


def _setting_value(text: str):
    """``true``/``false`` -> bool and whole numbers -> int, like the values the presets store; else the text."""
    if text.lower() in ("true", "false"):
        return text.lower() == "true"
    return int(text) if re.fullmatch(r"-?\d+", text) else text


def _set_preset_values(ws: Workspace, presets: list, assignments: list) -> int:
    """``presets set``: every assignment is checked against every preset before anything is written."""
    workflows = {name: ws.load_workflow(name) for name in presets}
    changes = []
    for item in assignments:
        key, sep, value = item.partition("=")
        if not sep or not key:
            log.error("'%s' is not KEY=VALUE", item)
            return 1
        for name, wf in workflows.items():
            step, dot, field = key.partition(".")
            if dot:
                if step not in wf.pipeline:
                    log.error("preset '%s' has no step '%s' (see 'mitopiper presets show %s')", name, step, name)
                    return 1
            elif key not in wf.general_settings:
                log.error("'%s' is not a general setting of preset '%s' (see 'mitopiper presets show %s')", key, name, name)
                return 1
            changes.append((name, wf, step if dot else None, field if dot else key, _setting_value(value)))
    for name, wf, step, field, value in changes:
        old = wf.pipeline[step].get(field, "") if step else wf.general(field)
        if step:
            wf.pipeline[step][field] = str(value) if field in ("memory", "threads") else value
        else:
            wf.set_general(field, value)
        off = " - this step is off in this preset" if step and step not in wf.active_order() else ""
        print(f"{name}: {step + '.' if step else ''}{field} = {value!r} (was {old!r}){off}")
    for name, wf in workflows.items():
        ws.write_workflow(wf, name)
    return 0


def cmd_refs(ws: Workspace, args) -> int:
    ws.initialise()
    registry = ReferenceRegistry(ws)
    if args.action == "list":
        rows = []
        for org in registry.organisms():
            for name in registry.names(org):
                ref = registry.get(org, name)
                rows.append([org, name, "*" if registry.active_name(org) == name else "", ref.refpath, ref.gtfpath,
                             ref.nuc_refpath or "", ref.numt_refpath or "", ref.source, ref.accession])
        print(table(rows, ["Organism", "Name", "Active", "FASTA", "GTF/GFF", "Nuclear FASTA", "NUMT FASTA", "Source",
                           "Accession"]) if rows else "No references registered.")
    elif args.action == "add":
        ref = ReferenceSet(refpath=args.fasta, gtfpath=args.gtf, nuc_refpath=args.nuc_fasta,
                           numt_refpath=args.numt_fasta, chrominfopath=args.chrominfo)
        if not Path(args.fasta).exists():
            log.warning("%s does not exist on this machine (fine if it exists on the cluster)", args.fasta)
        registry.add(args.organism.replace(" ", "_"), ref, args.name, activate=not args.no_activate)
        print(f"Registered reference '{args.name}' for {args.organism}")
    elif args.action == "use":
        registry.use(args.organism.replace(" ", "_"), args.name)
        print(f"{args.organism} now uses reference '{args.name}'")
    elif args.action == "remove":
        registry.remove(args.organism.replace(" ", "_"), args.name)
    elif args.action == "fetch-genome":
        from .references import attach_nuclear, fetch_genome_reference

        organism = args.organism.replace(" ", "_")
        try:
            ref, genome = fetch_genome_reference(ws, organism, args.accession)
        except RefError as exc:
            raise RefError(f"{exc}. Register an existing genome instead: "
                           f"mitopiper refs add {organism} --fasta /path/genome.fa --gtf /path/genes.gtf") from exc
        print(genome.warning())
        if not args.yes:
            if not sys.stdin.isatty():
                raise RefError("Confirm the download with --yes")
            if not _confirm("Continue?"):
                print("Nothing was changed.")
                return 1
        if args.nuclear_of:
            if args.nuclear_of not in registry.names(organism):
                raise RefError(f"{organism} has no reference set '{args.nuclear_of}'")
            registry.add(organism, attach_nuclear(registry.get(organism, args.nuclear_of), ref), args.nuclear_of, activate=False)
            print(f"'{args.nuclear_of}' now has the nuclear genome {genome.accession} as _fasta_nuc (its FASTA only: "
                  f"the set keeps its own annotation)")
        else:
            registry.add(organism, ref, args.name, activate=not args.no_activate)
            print(f"Registered '{args.name}' = {genome.accession} ({genome.assembly_name})")
        print(f"The files will be downloaded to {ref.refpath.rsplit('/', 1)[0]} by the reference job of the next build.")
    elif args.action == "fetch-mt":
        organism = args.organism.replace(" ", "_")
        ref = fetch_mt_reference(ws, organism, args.accession)
        registry.add(organism, ref, args.name, activate=not args.no_activate)
        print(f"Downloaded {ref.accession} to {ref.refpath}")
    return 0


def cmd_cleanup(ws: Workspace, args) -> int:
    removed = cleanup(ws, args.folders, args.experiment, args.keep_scripts, confirm=None if args.yes else _confirm)
    print(f"Removed {len(removed)} item(s).")
    return 0


def _confirm(question: str) -> bool:
    return input(f"{question} (y/n) ").strip().lower().startswith("y")


def cmd_fetch_metadata(ws: Workspace, args) -> int:
    from .ena import fetch_run_report, is_accession

    if not is_accession(args.accession):
        raise BuildError(f"'{args.accession}' does not look like an ENA/SRA/DDBJ accession")
    dest = Path(args.output) if args.output else ws.root / "metadata" / f"{args.accession.upper()}.tsv"
    fetch_run_report(args.accession, dest)
    print(f"Metadata written to {dest}")
    return 0


def cmd_check(ws: Workspace, args) -> int:
    from . import toolcheck
    from .monitor import default_experiments

    experiments = args.experiments or default_experiments(ws)
    if not experiments:
        print("No builds found. Run 'mitopiper build' first.")
        return 1
    # Experiments that share a job environment (shell + set-up lines) are checked in one go.
    groups: dict = {}
    for exp in experiments:
        tools, workflow = toolcheck.experiment_tools(ws, exp)
        key = (toolcheck.jobscript.shell_of(workflow), str(workflow.general("free_form", "") or ""))
        group = groups.setdefault(key, {"experiments": [], "tools": [], "workflow": workflow})
        group["experiments"].append(exp)
        group["tools"].append(tools)
    status = 0
    for group in groups.values():
        print(f"== {', '.join(group['experiments'])}")
        result = toolcheck.run_preflight(toolcheck.merge(*group["tools"]), group["workflow"], auto_fix=args.auto_fix)
        status |= int(result.failed)
    return status


def cmd_rerun(ws: Workspace, args) -> int:
    from .rerun import RerunOptions, rerun

    options = RerunOptions(
        jobs=args.job, from_step=args.from_step, from_start=args.from_start, clean=args.clean,
        include_pending=args.include_pending, auto_fix=args.auto_fix, preset=args.preset,
        scheduler=normalise(args.scheduler) if args.scheduler else None, memory=args.memory,
        walltime=args.walltime, threads=args.threads, queue=args.queue, free_form=args.free_form,
    )
    return rerun(ws, args.experiments, options, submit=args.submit, dry_run=args.dry_run)


def cmd_monitor(ws: Workspace, args) -> int:
    from .monitor import monitor

    if args.interactive:
        from .monitor import default_experiments
        from .session import MonitorSession

        experiments = args.experiments or default_experiments(ws)
        if not experiments:
            print("No submitted experiments found. Build and run master.sh first.")
            return 1
        return MonitorSession(ws, experiments, interval=args.watch or 15.0, color=False if args.no_color else None,
                              stay_when_done=sys.stdin.isatty()).run()
    if args.rerun:
        from .rerun import RerunOptions, rerun

        return rerun(ws, args.experiments, RerunOptions(auto_fix=True), submit=True)

    return monitor(
        ws, args.experiments, watch=args.watch, job=args.job, errors_only=args.errors, as_json=args.json,
        color=False if args.no_color else None, tail_lines=args.tail,
        extra_patterns=args.pattern, extra_ignore=args.ignore,
    )


def cmd_make_tsv(ws: Workspace, args) -> int:
    from .local_samples import make_tsv

    out = make_tsv(args.folder, args.org, Path(args.output) if args.output else None, args.lib_strat, args.ignore_paired, args.force_paired)
    print(f"Metadata written to {out}\nNext: mitopiper build {out} --steps <steps> (the wget step links local files)")
    return 0


def main(argv: Optional[list] = None) -> int:
    args = build_parser().parse_args(argv)
    level = logging.WARNING - 10 * min(args.verbose + 1, 2)
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")
    ws = Workspace(Path(args.workdir))
    handlers = {
        "build": cmd_build,
        "presets": cmd_presets,
        "refs": cmd_refs,
        "cleanup": cmd_cleanup,
        "make-tsv": cmd_make_tsv,
        "fetch-metadata": cmd_fetch_metadata,
        "monitor": cmd_monitor,
        "rerun": cmd_rerun,
        "check": cmd_check,
        "mt": lambda ws, args: __import__("mitopiper.mtdna", fromlist=["run"]).run(args),
        "download": lambda ws, args: __import__("mitopiper.download", fromlist=["run"]).run(args),
    }
    try:
        if args.command in (None, "interactive"):
            from .interactive import run_interactive

            return run_interactive(ws)
        if args.command == "init":
            ws.initialise()
            print(f"Workspace initialised in {ws.root}")
            return 0
        return handlers[args.command](ws, args)
    except (BuildError, ConfigError, MetadataError, RefError, KeyError, ValueError, FileNotFoundError) as exc:
        message = exc.args[0] if isinstance(exc, KeyError) and exc.args else exc
        log.error("%s", message)
        return 2
    except (KeyboardInterrupt, EOFError):  # Ctrl-C, or Ctrl-D / end of piped input
        print("\nAborted.")
        return 130
    except BrokenPipeError:  # output piped into e.g. `head`, which stopped reading
        import os

        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 0


if __name__ == "__main__":
    sys.exit(main())
