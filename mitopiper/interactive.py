"""Menu-driven interface (the original ``maker.py`` experience)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Callable, Optional

from . import jobscript, toolcheck
from .build import (
    BuildError,
    BuildOptions,
    build_experiment,
    find_experiment_metadata,
    list_experiments,
    resolve_metadata_source,
    tool_check_enabled,
    write_master,
)
from .builder import PLACEHOLDERS, output_exts, selected_steps
from .config import INTERNAL_PRESETS, ConfigError, Workflow, Workspace
from .ena import is_accession
from .metadata import MetadataError, load_samples, parse_selection
from .references import ReferenceRegistry, RefError, ReferenceSet, fetch_mt_reference
from .schedulers import SCHEDULERS
from .ui import ask, banner, choose, choose_many, table, yes_no

log = logging.getLogger(__name__)

PLACEHOLDER_HELP = (
    "Values may contain placeholders that are replaced per sample:\n  "
    + ", ".join("_" + p for p in PLACEHOLDERS)
    + "\nas well as output_<step> / dir_<step> (output file / directory of another step;\n"
    "output_<step>:<name> or output_<step>:all for steps with several outputs),\n"
    "{name} anywhere inside a value, and _if_SE / _if_PE for flags that only apply to\n"
    "single- or paired-end samples.\n"
)


def key_fixer(key: str) -> str:
    """Add ``-``/``--`` to a bare option name according to unix conventions."""
    key = key.strip()
    if not key or key.startswith("-"):
        return key
    return ("--" + key.replace(" ", "-")) if len(key) > 1 else "-" + key


class Menus:
    """Configuration menus operating on the ``current`` preset."""

    def __init__(self, workspace: Workspace, input_fn: Callable = input):
        self.ws = workspace
        self.input = input_fn
        self.wf: Workflow = workspace.load_workflow("current")

    def save(self) -> None:
        self.ws.write_workflow(self.wf, "current")

    # ----------------------------------------------------------- main menu
    def run(self) -> None:
        banner("Welcome to the configuration menu!")
        entries = [
            ("change general settings", self.change_general_settings),
            ("change pipeline step settings", self.change_pipeline_settings),
            ("activate/deactivate or reorder steps", self.modify_pipeline),
            ("add a new setting to a step", self.add_new_pipeline_setting),
            ("add a new pipeline step", self.add_new_pipeline_step),
            ("modify a pipeline step", self.modify_pipeline_step),
            ("manage references", self.manage_references),
            ("save current settings as a preset", self.save_current_settings),
            ("copy steps/settings between presets", self.copy_between_presets),
        ]
        while True:
            idx = choose("\nWhat would you like to do?", [e[0] for e in entries], input_fn=self.input)
            if idx is None:
                return
            entries[idx][1]()

    # ----------------------------------------------------- general settings
    def change_general_settings(self) -> None:
        gs = self.wf.general_settings
        while True:
            names = list(gs)
            print(table([[i, n, gs[n][0], gs[n][1]] for i, n in enumerate(names)], ["#", "Name", "Value", "Description"]))
            idx = choose("Which setting would you like to change?", names, input_fn=self.input)
            if idx is None:
                break
            name = names[idx]
            current = gs[name][0]
            if isinstance(current, bool):
                gs[name][0] = not current
                print(f"{name} set to {gs[name][0]}\n")
            elif name == "free_form":
                value = self.input(f"Enter a line to append to the free-form block (empty to wipe it)\nCurrent:\n{current}\n> ")
                gs[name][0] = "" if value == "" else f"{current}{value}\n"
            else:
                gs[name][0] = ask(f"New value for {name} ({gs[name][1]})", str(current), self.input)
            self.save()

    # --------------------------------------------------- pipeline settings
    def change_pipeline_settings(self) -> None:
        active = self.wf.active_order()
        while True:
            idx = choose("Choose an active step to configure:", active, input_fn=self.input)
            if idx is None:
                return
            step = active[idx]
            params = self.wf.settings.setdefault(step, {})
            if not params:
                print(f"No settings were added to {step} yet; use 'add a new setting to a step'.\n")
                continue
            names = sorted(params, key=lambda n: params[n].get("order", 0))
            rows = [[i, n, params[n].get("key"), params[n].get("value"), params[n].get("switch"), params[n].get("description")] for i, n in enumerate(names)]
            print(table(rows, ["#", "Name", "Key", "Value", "On", "Description"]))
            answer = self.input("Enter the number of a setting to modify, 'r' to reorder, 'exit' to go back\n> ").strip()
            if answer.lower() == "r":
                self.change_order(step)
                continue
            if not answer.isdigit() or int(answer) >= len(names):
                continue
            name = names[int(answer)]
            param = params[name]
            action = choose(f"{name}:", ["toggle on/off", "change value"], input_fn=self.input)
            if action == 0:
                param["switch"] = not bool(param.get("switch", True))
            elif action == 1:
                print(PLACEHOLDER_HELP)
                order = self.wf.run_order
                if step in order:
                    print("Steps before this one: " + ", ".join(order[: order.index(step)]))
                param["value"] = ask("New value", str(param.get("value", "")), self.input)
            self.save()

    def change_order(self, step: str) -> None:
        params = self.wf.settings[step]
        names = sorted(params, key=lambda n: params[n].get("order", 0))
        for i, n in enumerate(names):
            print(f"  {i}\t{n}\t{params[n].get('key')}")
        chosen = choose_many("Enter the new order (e.g. '2 0 1'); unlisted parameters keep their relative order", len(names), self.input)
        if chosen is None:
            return
        new = [names[i] for i in chosen] + [n for i, n in enumerate(names) if i not in chosen]
        for i, n in enumerate(new):
            params[n]["order"] = i
        self.save()

    # -------------------------------------------------------- run order
    def modify_pipeline(self) -> None:
        while True:
            print("Current run order:\n" + "\n".join(f"  {s}" for s in self.wf.active_order()) + "\n")
            idx = choose("", ["activate/deactivate steps", "reorder the run order"], input_fn=self.input)
            if idx is None:
                return
            steps = list(self.wf.pipeline)
            if idx == 0:
                rows = [[i, s, "on" if self.wf.pipeline[s].get("switch") else "off", self.wf.pipeline[s].get("scope", "sample")] for i, s in enumerate(steps)]
                print(table(rows, ["#", "Step", "State", "Scope"]))
                chosen = choose_many("Enter the step(s) to toggle", len(steps), self.input)
                if chosen is None:
                    continue
                for i in chosen:
                    step = steps[i]
                    self.wf.pipeline[step]["switch"] = not bool(self.wf.pipeline[step].get("switch"))
                    if self.wf.pipeline[step]["switch"] and step not in self.wf.run_order:
                        self.wf.run_order.append(step)
            else:
                order = self.wf.active_order()
                for i, s in enumerate(order):
                    print(f"  {i}\t{s}")
                chosen = choose_many("Enter the new order (e.g. '0 2 1 3'); steps left out are deactivated", len(order), self.input)
                if chosen is None:
                    continue
                new = [order[i] for i in chosen]
                for step in order:
                    if step not in new:
                        self.wf.pipeline[step]["switch"] = False
                self.wf.run_order = new
            self.wf.run_order = [s for s in self.wf.run_order if s in self.wf.pipeline and self.wf.pipeline[s].get("switch")]
            self.save()

    # ---------------------------------------------------- new settings
    def add_new_pipeline_setting(self) -> None:
        steps = list(self.wf.pipeline)
        idx = choose("Which step should get a new setting?", steps, input_fn=self.input)
        if idx is None:
            return
        step = steps[idx]
        params = self.wf.settings.setdefault(step, {})
        name = ask("Setting name", input_fn=self.input)
        if not name:
            return
        description = ask("Description", input_fn=self.input)
        key = key_fixer(ask(f"Option key as documented for {step} (e.g. -q or --quality)", input_fn=self.input))
        print(PLACEHOLDER_HELP)
        value = ask("Value (leave empty for an on/off flag)", input_fn=self.input)
        params[name] = {"description": description, "key": key, "value": value, "switch": True, "order": len(params)}
        self.save()
        print(f"Setting added: {name}\t{key} {value}\n")

    def add_new_pipeline_step(self) -> None:
        explanations = {
            "key": "The command that runs the step (e.g. samtools view)",
            "description": "What the step does",
            "path": "Absolute path of the executable (leave empty to use $PATH)",
            "outdir": "Output directory (relative to the experiment directory)",
            "output_ext": "Extension of the produced file (e.g. .bam). Several outputs: name=ext pairs separated by ';' "
            "(e.g. bam=.bam;log=.log), referenced as output_<step>:<name>",
            "inputdir": "Directory of the input file (e.g. alignment)",
            "input_ext": "Extension of the input file (e.g. .sam)",
            "structure": "Command template; empty for 'key settings input'. Tokens: key, settings, input, output, "
            "output_pipe, output_<step>, dir_<step> and the placeholders listed in the manual",
            "scope": "'sample' (one job per sample) or 'experiment' (runs once after all samples)",
        }
        name = ask("Name of the new step (must be unique)", input_fn=self.input).replace(" ", "_")
        if not name or name in self.wf.pipeline:
            print("Invalid or duplicate step name.\n")
            return
        step = {}
        for field, text in explanations.items():
            default = "sample" if field == "scope" else ""
            step[field] = ask(f"{field}: {text}", default, self.input)
        if step["scope"] not in ("sample", "experiment"):
            step["scope"] = "sample"
        step["switch"] = True
        self.wf.pipeline[name] = step
        self.wf.settings.setdefault(name, {})
        self.wf.run_order.append(name)
        self.save()
        print("Pipeline step added.\n")

    def modify_pipeline_step(self) -> None:
        steps = list(self.wf.pipeline)
        idx = choose("Which step would you like to modify?", steps, input_fn=self.input)
        if idx is None:
            return
        step = self.wf.pipeline[steps[idx]]
        fields = [f for f in step if f != "switch"]
        print(table([[i, f, step[f]] for i, f in enumerate(fields)], ["#", "Field", "Value"]))
        fidx = choose("Which field?", fields, input_fn=self.input)
        if fidx is None:
            return
        field = fields[fidx]
        current = step[field]
        if field == "output_ext" and isinstance(current, (dict, list)):
            current = ";".join(f"{k}={v}" for k, v in output_exts(step).items())
        step[field] = ask(f"New value for {field}", str(current), self.input)
        self.save()

    # ------------------------------------------------------ references
    def manage_references(self) -> None:
        registry = ReferenceRegistry(self.ws)
        while True:
            rows = []
            for org in registry.organisms():
                for name in registry.names(org):
                    ref = registry.get(org, name)
                    active = "*" if registry.active_name(org) == name else ""
                    rows.append([org, name + active, ref.refpath, ref.gtfpath, ref.source])
            print(table(rows, ["Organism", "Name (*active)", "FASTA", "GTF/GFF", "Source"]) if rows else "No references registered yet.\n")
            action = choose(
                "References:",
                ["add/replace a reference", "switch active reference", "download mtDNA reference from NCBI", "remove a reference",
                 "download the nuclear genome from NCBI Genome (not recommended: very large)"],
                input_fn=self.input,
            )
            if action is None:
                return
            try:
                if action == 0:
                    org = ask("Organism (e.g. Homo_sapiens)", input_fn=self.input).replace(" ", "_")
                    name = ask("Name of this reference set", "default", self.input)
                    registry.add(org, prompt_reference_paths(self.input), name, activate=True)
                elif action in (1, 3):
                    orgs = registry.organisms()
                    oi = choose("Organism:", orgs, input_fn=self.input)
                    if oi is None:
                        continue
                    names = registry.names(orgs[oi])
                    ni = choose("Reference:", names, input_fn=self.input)
                    if ni is None:
                        continue
                    (registry.use if action == 1 else registry.remove)(orgs[oi], names[ni])
                elif action == 4:
                    org = ask("Organism (e.g. Homo_sapiens)", input_fn=self.input).replace(" ", "_")
                    ref = fetch_genome_interactive(self.ws, org, self.input)
                    if ref is not None:
                        name = ask("Name of this reference set", "genome" if ref.source == "ncbi_genome" else "default", self.input)
                        registry.add(org, ref, name, activate=True)
                elif action == 2:
                    org = ask("Organism (e.g. Homo_sapiens)", input_fn=self.input).replace(" ", "_")
                    acc = ask("RefSeq accession (leave empty to search NCBI)", input_fn=self.input)
                    name = ask("Name of this reference set", "mt", self.input)
                    registry.add(org, fetch_mt_reference(self.ws, org, acc or None), name, activate=True)
            except RefError as exc:
                print(f"Error: {exc}\n")

    # --------------------------------------------------------- presets
    def save_current_settings(self) -> None:
        name = ask("Name for the new preset", input_fn=self.input)
        try:
            self.ws.save_preset(name)
            print(f"Preset '{name}' saved.\n")
        except ConfigError as exc:
            print(f"Error: {exc}\n")

    def copy_between_presets(self) -> None:
        presets = [p for p in self.ws.list_presets() if p not in INTERNAL_PRESETS]
        src = choose("Copy steps and settings FROM:", presets, input_fn=self.input)
        if src is None:
            return
        dst = choose("... TO:", presets, input_fn=self.input)
        if dst is None:
            return
        source = self.ws.load_workflow(presets[src])
        target = self.ws.load_workflow(presets[dst])
        for step, conf in source.pipeline.items():
            target.pipeline.setdefault(step, conf)
        for step, params in source.settings.items():
            target.settings.setdefault(step, params)
        self.ws.write_workflow(target, presets[dst])
        print(f"Copied missing steps/settings from {presets[src]} to {presets[dst]} (nothing was overwritten).\n")


def prompt_reference_paths(input_fn: Callable = input) -> ReferenceSet:
    fasta = ask("Full path to the reference FASTA (pre-indexed for your aligner)", input_fn=input_fn)
    nuc = ask("Full path to a nuclear reference FASTA (optional)", input_fn=input_fn)
    gtf = ask("Full path to the GTF/GFF annotation (optional)", input_fn=input_fn)
    chrom = ask("Full path to a ChromInfo file (optional)", input_fn=input_fn)
    return ReferenceSet(refpath=fasta.strip('"'), nuc_refpath=nuc.strip('"'), gtfpath=gtf.strip('"'), chrominfopath=chrom.strip('"'))


def fetch_genome_interactive(workspace: Workspace, organism: str, input_fn: Callable = input, client=None):
    """Find the nuclear genome in NCBI Genome, warn about its size and confirm.

    If NCBI has none (or it cannot be reached), the user can enter a path instead.
    Returns the reference set, or ``None`` to go back.
    """
    from .references import fetch_genome_reference

    print(f"Looking up the {organism.replace('_', ' ')} genome in the NCBI Genome database ...")
    try:
        ref, genome = fetch_genome_reference(workspace, organism, client=client)
    except RefError as exc:
        print(f"Error: {exc}")
        if yes_no("Enter the path to a genome FASTA manually instead?", True, input_fn):
            return prompt_reference_paths(input_fn)
        return None
    print("\n" + genome.warning() + "\n")
    if not yes_no(f"Download {genome.accession} ({genome.assembly_name}) anyway?", False, input_fn):
        return None
    return ref


def interactive_resolver(workspace: Workspace, input_fn: Callable = input):
    """Ask the user what to do when an organism has no reference."""
    registry = ReferenceRegistry(workspace)

    def resolve(organism: str) -> Optional[ReferenceSet]:
        while True:
            idx = choose(
                f"No reference is registered for {organism}.",
                [
                    "download the mitochondrial genome from NCBI automatically",
                    "enter paths to an existing reference (recommended)",
                    "skip this organism",
                    "download the whole nuclear genome from NCBI Genome (not recommended: very large)",
                ],
                allow_exit=False,
                input_fn=input_fn,
            )
            try:
                if idx == 0:
                    ref = fetch_mt_reference(workspace, organism)
                    registry.add(organism, ref, "mt")
                    return ref
                if idx == 1:
                    ref = prompt_reference_paths(input_fn)
                    registry.add(organism, ref, "default")
                    return ref
                if idx == 3:
                    ref = fetch_genome_interactive(workspace, organism, input_fn)
                    if ref is not None:
                        registry.add(organism, ref, "genome" if ref.source == "ncbi_genome" else "default")
                        return ref
                    continue
                return None
            except RefError as exc:
                print(f"Error: {exc}\n")

    return resolve


def _choose_metadata(workspace: Workspace, input_fn: Callable) -> list:
    experiments = list_experiments(workspace)
    print("\nPrevious experiments:")
    for i, name in enumerate(experiments):
        print(f"  {i}\t{name}")
    print(
        "\nEnter the path of a new ENA/SRA metadata file, a project/study accession to fetch from ENA\n"
        "(e.g. PRJEB12345, PRJNA67890, SRP012345), or the number(s) of previous experiments\n"
        "(e.g. '2', '0-3', '0 2 5' or 'all')."
    )
    while True:
        answer = input_fn("> ").strip().strip('"')
        if not answer:
            continue
        candidate = Path(answer).expanduser()
        for path in (candidate, candidate.with_suffix(".tsv"), candidate.with_suffix(".txt"), candidate.with_suffix(".csv")):
            if path.is_file():
                return [path]
        if is_accession(answer):
            try:
                return [resolve_metadata_source(workspace, answer)]
            except (BuildError, MetadataError) as exc:
                print(f"Error: {exc}\n")
                continue
        try:
            chosen = parse_selection(answer, len(experiments)) if experiments else []
        except ValueError:
            chosen = []
        if chosen:
            return [find_experiment_metadata(workspace.experiment_dir(experiments[i])) for i in chosen]
        print(f"'{answer}' is not an existing file, an accession or a valid experiment selection.\n")


def preflight(workspace: Workspace, results: list, input_fn: Callable, rebuild: Callable,
              check: Callable = None, run_install: Callable = None) -> tuple:
    """Make sure the tools the jobs call exist before launching.

    Returns ``(results, proceed)``; ``results`` changes when a set-up line was
    added and the scripts were rebuilt.
    """
    check = check or toolcheck.check_here
    run_install = run_install or toolcheck.install
    while True:
        workflow = workspace.load_workflow("current")
        tools = toolcheck.merge(*[r.tools for r in results])
        result = check(tools, workflow)
        if result.ok:
            print(result.describe())
            return results, True
        if not result.checked:
            print(result.describe() + "\nThe check job will verify the tools on the cluster before anything else runs.")
            return results, True
        print("\n" + result.describe())
        idx = choose(
            "How do you want to continue?",
            [
                "I will install / load them myself (then check again)",
                "let mitopiper try to install them into the current environment",
                "add a set-up line to every job (e.g. 'module load ...') and regenerate the scripts",
                "launch anyway (the check job stops the pipeline on the cluster if they are really missing)",
                "stop here (the scripts are kept)",
            ],
            allow_exit=False, input_fn=input_fn,
        )
        if idx == 0:
            input_fn("Install or load the tools in this environment, then press Enter to check again ")
        elif idx == 1:
            plan = toolcheck.install_plan(result.missing)
            print(plan.describe())
            if not plan.commands:
                print("Nothing can be installed automatically here; choose another option.")
                continue
            if yes_no("Run these commands now?", True, input_fn):
                failed = run_install(plan)
                if failed:
                    print(f"{len(failed)} install command(s) failed; see the output above.")
                print(toolcheck.ENV_NOTE)
        elif idx == 2:
            suggestion = f"module load {' '.join(result.missing)}"
            line = ask("Set-up line added to every job", suggestion, input_fn)
            if line:
                results = rebuild(line)
                print("The scripts were regenerated with the new set-up line (saved in the current preset).")
        elif idx == 3:
            return results, True
        else:
            return results, False


def _last_build(workspace: Workspace) -> list:
    """Experiments of the last build that were actually launched."""
    import json

    try:
        last = json.loads((workspace.root / "last_build.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    return [e for e in last.get("experiments", []) if (workspace.experiment_dir(e) / "jobs.tsv").exists()]


def _default_session(workspace: Workspace, experiments: list, input_fn: Callable):
    from .session import MonitorSession

    return MonitorSession(workspace, experiments, ask=input_fn)


def run_interactive(workspace: Workspace, input_fn: Callable = input, session_factory: Callable = _default_session) -> int:
    """Build -> launch -> monitor -> rerun, all from one interactive session.

    ``session_factory(workspace, experiments, input_fn)`` creates the monitor
    session (replaceable for testing).
    """
    workspace.initialise()
    banner("Welcome to mitopiper - job script generator for sequencing pipelines")

    previous = _last_build(workspace)
    if previous:
        what = choose("What would you like to do?",
                      ["build and launch new jobs", f"monitor / rerun the jobs of the last build ({', '.join(previous)})"],
                      allow_exit=False, input_fn=input_fn)
        if what == 1:
            return session_factory(workspace, previous, input_fn).run()

    idx = choose("Choose your job scheduler:", list(SCHEDULERS), allow_exit=False, input_fn=input_fn)
    scheduler = SCHEDULERS[idx]
    metadata_files = _choose_metadata(workspace, input_fn)

    presets = [p for p in workspace.list_presets() if p != "current"]
    pidx = choose("Choose a workflow preset:", presets, allow_exit=False, input_fn=input_fn)
    workspace.activate(presets[pidx])

    if yes_no("Would you like to change the configuration first?", False, input_fn):
        Menus(workspace, input_fn).run()
    workflow = workspace.load_workflow("current")

    resolver = interactive_resolver(workspace, input_fn)
    results = []
    builds = []  # (metadata, selection, options, samples) to rebuild after a set-up change
    for metadata in metadata_files:
        ignore_paired = bool(workflow.general("ignore_paired_end", False))
        try:
            samples = load_samples(metadata, ignore_paired=ignore_paired)
        except MetadataError as exc:
            print(f"Error: {exc}\nSkipping {metadata}.\n")
            continue
        rows = [[s.index, s.organism, " ".join(s.run_accessions), s.accession, s.title, s.library_strategy, "PE" if s.paired else "SE"] for s in samples]
        print(table(rows, ["Index", "Organism", "Runs", "Sample", "Title", "Strategy", "Layout"]))
        selection = choose_many("Which samples should be processed? (e.g. 'all', '2', '1-4', '0 3 5')", len(samples), input_fn, allow_exit=False)

        steps = selected_steps(workflow)
        print("Steps that will run: " + " -> ".join(steps))
        answer = input_fn("Press enter to keep them all, or list the steps to run (space separated)\n> ").strip()
        chosen_steps = answer.split() if answer else None
        options = BuildOptions(scheduler=scheduler, steps=chosen_steps)
        try:
            result = build_experiment(workspace, workflow, metadata, selection, options, resolver, samples=samples)
        except (KeyError, MetadataError) as exc:
            print(f"Error: {exc}\n")
            continue
        results.append(result)
        builds.append((metadata, selection, options, samples))
        count = len(result.jobs.scripts())
        print(f"{count} job script(s) written to {result.directory / 'scripts'}\n")

    if not results:
        print("Nothing was generated.")
        return 1

    def rebuild(setup_line: str) -> list:
        """Add a set-up line to every job (saved in the current preset) and regenerate the scripts."""
        wf = workspace.load_workflow("current")
        jobscript.append_free_form(wf, [setup_line])
        workspace.write_workflow(wf)
        return [build_experiment(workspace, wf, m, sel, opt, resolver, samples=smp) for m, sel, opt, smp in builds]

    proceed = True
    if tool_check_enabled(workflow):
        results, proceed = preflight(workspace, results, input_fn, rebuild)
        workflow = workspace.load_workflow("current")
    master = write_master(workspace, scheduler, results, shell=str(workflow.general("environment", "bash")))
    workspace.copy_preset("current", "last_used")
    if not proceed:
        print("Stopped before launching; the scripts are kept. Fix the environment and run 'bash master.sh', "
              "or start mitopiper again.")
        return 1
    total = sum(len(r.jobs.samples) for r in results)
    print(f"\nAll scripts are ready: {total} sample job(s), submitted by {master}.")
    if scheduler == "local":
        print("NOTE: the 'local' scheduler runs the jobs one after another on this machine.")
    else:
        print("IMPORTANT: every sample is submitted at once - make sure your cluster allows that many jobs!")
    if not yes_no("Launch the jobs now and follow them in the monitor?", True, input_fn):
        print(f"Submit later with 'bash {master}', then follow the jobs with 'mitopiper monitor -i'.")
        return 0
    session = session_factory(workspace, [r.name for r in results if r.jobs.samples or r.jobs.experiment], input_fn)
    session.start(master, "master.log")
    return session.run()
