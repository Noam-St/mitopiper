"""Workspace and preset (workflow) configuration handling.

A *workspace* is the directory pipemaker operates in (the current directory by
default).  Its layout is::

    <workspace>/
        config/
            default/        built-in defaults (regenerated keys are merged in)
            current/        the working copy used for the current invocation
            last_used/      snapshot of the configuration of the last build
            <preset>/       user-saved presets
            global/         settings shared by every preset (refpaths.json)
            credentials.json  optional NCBI e-mail / API key
        output/<experiment>/  generated scripts and pipeline outputs
        references/           automatically downloaded references

Every preset directory holds four JSON files: ``pipeline``, ``settings``,
``run_order`` and ``general_settings``.  The layout is identical to the one
used by the original ``maker.py`` so existing workspaces keep working.
"""

from __future__ import annotations

import copy
import json
import logging
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .defaults import CONFIG_FILES, default_config

log = logging.getLogger(__name__)

RESERVED_PRESETS = ("current", "global")
INTERNAL_PRESETS = ("default", "current", "global", "last_used")


class ConfigError(Exception):
    """Raised for missing or malformed configuration."""


@dataclass
class Workflow:
    """The four configuration files that define a workflow."""

    pipeline: dict
    settings: dict
    run_order: list
    general_settings: dict

    def general(self, key: str, default: Any = None) -> Any:
        """Return the value of a general setting (stored as ``[value, description]``)."""
        entry = self.general_settings.get(key)
        if entry is None:
            return default
        if isinstance(entry, list):
            return entry[0] if entry else default
        return entry

    def set_general(self, key: str, value: Any) -> None:
        entry = self.general_settings.get(key)
        if isinstance(entry, list) and len(entry) > 1:
            entry[0] = value
        else:
            self.general_settings[key] = [value, ""]

    def active_order(self) -> list:
        """Run order restricted to steps that exist and are switched on."""
        return [s for s in self.run_order if s in self.pipeline and _truthy(self.pipeline[s].get("switch", True))]

    def as_dict(self) -> dict:
        return {
            "pipeline": self.pipeline,
            "settings": self.settings,
            "run_order": self.run_order,
            "general_settings": self.general_settings,
        }

    def copy(self) -> "Workflow":
        return Workflow(**copy.deepcopy(self.as_dict()))

    @classmethod
    def default(cls) -> "Workflow":
        return cls(**default_config())


def _truthy(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "y", "on")
    return bool(value)


@dataclass
class Workspace:
    """Access to the configuration and output directories of a workspace."""

    root: Path
    _initialised: bool = field(default=False, repr=False)

    def __post_init__(self) -> None:
        self.root = Path(self.root).expanduser().resolve()

    # ------------------------------------------------------------------ paths
    @property
    def config_dir(self) -> Path:
        return self.root / "config"

    @property
    def output_dir(self) -> Path:
        return self.root / "output"

    @property
    def references_dir(self) -> Path:
        return self.root / "references"

    def preset_dir(self, preset: str) -> Path:
        return self.config_dir / preset

    def experiment_dir(self, name: str) -> Path:
        return self.output_dir / name

    # ------------------------------------------------------------- raw files
    def read_json(self, name: str, preset: str = "current") -> Any:
        path = self.preset_dir(preset) / f"{name}.json"
        try:
            with open(path, "r", encoding="utf-8") as handle:
                return json.load(handle)
        except FileNotFoundError as exc:
            raise ConfigError(f"Missing config file {path}") from exc
        except json.JSONDecodeError as exc:
            raise ConfigError(f"Malformed JSON in {path}: {exc}") from exc

    def write_json(self, name: str, data: Any, preset: str = "current") -> Path:
        directory = self.preset_dir(preset)
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"{name}.json"
        tmp = path.with_suffix(".json.tmp")
        with open(tmp, "w", encoding="utf-8") as handle:
            json.dump(data, handle, sort_keys=True, indent=4)
        os.replace(tmp, path)
        return path

    # ---------------------------------------------------------------- set up
    def initialise(self) -> None:
        """Create the directory skeleton and default preset, then migrate presets."""
        if self._initialised:
            return
        for directory in (self.config_dir, self.output_dir, self.preset_dir("global")):
            directory.mkdir(parents=True, exist_ok=True)
        if not (self.preset_dir("default") / "pipeline.json").exists():
            self.write_workflow(Workflow.default(), "default")
            log.info("Created default preset in %s", self.preset_dir("default"))
        self.migrate_presets()
        self._initialised = True

    def migrate_presets(self) -> None:
        """Merge keys that were added to the built-in defaults into every preset.

        Existing values are never overwritten.  Settings parameters without an
        ``order`` field are given one based on their current position.
        """
        defaults = default_config()
        for preset in self.list_presets(include_internal=True):
            if preset == "global":
                continue
            for name in CONFIG_FILES:
                path = self.preset_dir(preset) / f"{name}.json"
                if not path.exists():
                    continue
                try:
                    data = self.read_json(name, preset)
                except ConfigError as exc:
                    log.warning("%s", exc)
                    continue
                changed = False
                if name in ("pipeline", "settings", "general_settings") and isinstance(data, dict):
                    for key, value in defaults[name].items():
                        if key not in data:
                            data[key] = copy.deepcopy(value)
                            changed = True
                if name == "pipeline":
                    for step in data.values():
                        if isinstance(step, dict) and "scope" not in step:
                            step["scope"] = "sample"
                            changed = True
                if name == "settings":
                    for params in data.values():
                        if not isinstance(params, dict):
                            continue
                        for i, param in enumerate(params.values()):
                            if isinstance(param, dict) and "order" not in param:
                                param["order"] = i
                                changed = True
                if changed:
                    self.write_json(name, data, preset)
                    log.debug("Migrated %s/%s.json", preset, name)

    # --------------------------------------------------------------- presets
    def list_presets(self, include_internal: bool = False) -> list:
        if not self.config_dir.exists():
            return []
        presets = sorted(p.name for p in self.config_dir.iterdir() if p.is_dir())
        if include_internal:
            return presets
        return [p for p in presets if p not in RESERVED_PRESETS]

    def load_workflow(self, preset: str = "current") -> Workflow:
        if not self.preset_dir(preset).is_dir():
            raise ConfigError(f"Preset '{preset}' does not exist in {self.config_dir}")
        data = {name: self.read_json(name, preset) for name in CONFIG_FILES}
        return Workflow(**data)

    def write_workflow(self, workflow: Workflow, preset: str = "current") -> None:
        for name, data in workflow.as_dict().items():
            self.write_json(name, data, preset)

    def copy_preset(self, source: str, target: str) -> None:
        """Copy every JSON file of ``source`` into ``target`` (overwriting)."""
        src = self.preset_dir(source)
        if not src.is_dir():
            raise ConfigError(f"Preset '{source}' does not exist")
        dst = self.preset_dir(target)
        dst.mkdir(parents=True, exist_ok=True)
        for path in src.glob("*.json"):
            shutil.copyfile(path, dst / path.name)

    def activate(self, preset: str) -> Workflow:
        """Make ``preset`` the working copy (``current``) and return it."""
        if preset in RESERVED_PRESETS:
            raise ConfigError(f"'{preset}' is reserved and cannot be used as a preset")
        self.copy_preset(preset, "current")
        return self.load_workflow("current")

    def save_preset(self, name: str, source: str = "current") -> None:
        if not name or name in RESERVED_PRESETS or name == "default" or os.sep in name:
            raise ConfigError(f"Invalid preset name '{name}'")
        self.copy_preset(source, name)

    # ----------------------------------------------------------- credentials
    def ncbi_credentials(self) -> tuple:
        """Return ``(email, api_key)`` from the environment or credentials.json."""
        email = os.environ.get("NCBI_EMAIL", "")
        api_key = os.environ.get("NCBI_API_KEY", "")
        path = self.config_dir / "credentials.json"
        if path.exists():
            try:
                with open(path, "r", encoding="utf-8") as handle:
                    creds = json.load(handle)
                email = email or creds.get("Email", "") or creds.get("email", "")
                api_key = api_key or creds.get("API", "") or creds.get("api_key", "")
            except (OSError, json.JSONDecodeError) as exc:
                log.warning("Could not read %s: %s", path, exc)
        return email, api_key
