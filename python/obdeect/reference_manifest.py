"""Create and verify reproducible optical reference manifests.

The configuration is intentionally explicit: an absent executable, convention,
or input file is an error, never an implicit local default.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any

from obdeect.model_import import ImportError as ModelImportError
from obdeect.model_import import resolve_model, sha256


class ManifestError(ValueError):
    """A reference manifest is incomplete or no longer matches its inputs."""


_REQUIRED = (
    "sim_telarray_release",
    "hessio_decoder",
    "site",
    "production_version",
    "telescope_variants",
    "configuration_overrides",
    "seeds",
    "atmosphere_extinction",
    "photon_blocks",
    "coordinate_frames",
    "software_revisions",
    "commands",
)


def _object(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ManifestError(f"{label} must be an object")
    return value


def _git_revision(root: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode or len(result.stdout.strip()) != 40:
        raise ManifestError(f"cannot determine Git revision for {root}")
    return result.stdout.strip()


def _file(path: str, label: str) -> Path:
    if not isinstance(path, str) or not Path(path).is_absolute():
        raise ManifestError(f"{label} must be an absolute file path")
    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise ManifestError(f"{label} is missing: {path}")
    return resolved


def _validate(config: dict[str, Any]) -> None:
    missing = set(_REQUIRED) - config.keys()
    if missing:
        raise ManifestError(f"missing required fields: {', '.join(sorted(missing))}")
    for key in ("sim_telarray_release", "hessio_decoder", "site", "production_version"):
        if not isinstance(config[key], str) or not config[key].strip():
            raise ManifestError(f"{key} must be a nonempty string")
    for key in (
        "configuration_overrides",
        "seeds",
        "atmosphere_extinction",
        "coordinate_frames",
        "software_revisions",
    ):
        if not _object(config[key], key):
            raise ManifestError(f"{key} must be nonempty")
    variants = config["telescope_variants"]
    if (
        not isinstance(variants, list)
        or not variants
        or any(not isinstance(item, str) or not item for item in variants)
        or len(set(variants)) != len(variants)
    ):
        raise ManifestError("telescope_variants must be distinct model names")
    blocks = config["photon_blocks"]
    if not isinstance(blocks, list) or not blocks:
        raise ManifestError("photon_blocks must list input files")
    for block in blocks:
        if not isinstance(block, str):
            raise ManifestError("photon_blocks entries must be paths")
    commands = config["commands"]
    if not isinstance(commands, list) or not commands:
        raise ManifestError("commands must list reference invocations")
    for command in commands:
        command = _object(command, "command")
        if set(command) != {"name", "argv", "cwd", "environment", "inputs"}:
            raise ManifestError("command requires name, argv, cwd, environment, inputs")
        if not isinstance(command["name"], str) or not command["name"]:
            raise ManifestError("command name must be nonempty")
        if (
            not isinstance(command["argv"], list)
            or not command["argv"]
            or any(not isinstance(arg, str) or not arg for arg in command["argv"])
        ):
            raise ManifestError("command argv must be a nonempty string list")
        if not isinstance(command["cwd"], str) or not command["cwd"]:
            raise ManifestError("command cwd must be a path")
        if not Path(command["cwd"]).is_absolute() or not Path(command["cwd"]).is_dir():
            raise ManifestError(f"command cwd is missing: {command['cwd']}")
        _file(command["argv"][0], "command executable")
        if not isinstance(command["environment"], dict) or any(
            not isinstance(key, str) or not isinstance(value, str)
            for key, value in command["environment"].items()
        ):
            raise ManifestError("command environment must map strings to strings")
        if not isinstance(command["inputs"], list) or any(
            not isinstance(path, str) for path in command["inputs"]
        ):
            raise ManifestError("command inputs must list file paths")


def freeze(config: dict[str, Any], model_root: Path) -> dict[str, Any]:
    """Resolve selected model versions and hash every declared reference input."""
    _validate(config)
    version = config["production_version"]
    try:
        models = {
            name: resolve_model(model_root, name, version) for name in config["telescope_variants"]
        }
    except ModelImportError as error:
        raise ManifestError(str(error)) from error
    paths = set(config["photon_blocks"])
    for command in config["commands"]:
        paths.add(command["argv"][0])
        paths.update(command["inputs"])
    hashes = {path: sha256(_file(path, "reference input")) for path in sorted(paths)}
    return {
        "format": "obdeect.reference-manifest.v1",
        "configuration": config,
        "model_checkout_revision": _git_revision(model_root),
        "models": models,
        "file_sha256": hashes,
    }


def verify(manifest: dict[str, Any], model_root: Path) -> None:
    """Fail if a configured file, selected production, or checkout changed."""
    if manifest.get("format") != "obdeect.reference-manifest.v1":
        raise ManifestError("unsupported reference manifest format")
    config = _object(manifest.get("configuration"), "configuration")
    expected = freeze(config, model_root)
    if manifest != expected:
        raise ManifestError("reference manifest differs from current revisions or file hashes")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze", "verify"))
    parser.add_argument("path", type=Path, help="configuration or frozen manifest JSON")
    parser.add_argument("--model-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, help="required for freeze")
    args = parser.parse_args()
    try:
        data = _object(json.loads(args.path.read_text(encoding="utf-8")), "JSON root")
        if args.action == "freeze":
            if args.output is None:
                raise ManifestError("freeze requires --output")
            result = freeze(data, args.model_root)
            args.output.write_text(
                json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
        else:
            verify(data, args.model_root)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise SystemExit(f"reference manifest failed: {error}") from error


if __name__ == "__main__":
    main()
