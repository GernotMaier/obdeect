"""Create and verify reproducible reference-run records for optical comparisons.

The configuration is intentionally explicit: an absent executable, convention,
or input file is an error, never an implicit local default.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import subprocess
import time
from pathlib import Path
from typing import Any

from obdeect.model_import import ImportError as ModelImportError
from obdeect.model_import import resolve_model, sha256


class ReferenceRunError(ValueError):
    """A reference-run record is incomplete or no longer matches its inputs."""


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
        raise ReferenceRunError(f"{label} must be an object")
    return value


def _git_revision(root: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode or len(result.stdout.strip()) != 40:
        raise ReferenceRunError(f"cannot determine Git revision for {root}")
    return result.stdout.strip()


def _file(path: str, label: str) -> Path:
    if not isinstance(path, str) or not Path(path).is_absolute():
        raise ReferenceRunError(f"{label} must be an absolute file path")
    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise ReferenceRunError(f"{label} is missing: {path}")
    return resolved


def _validate(config: dict[str, Any]) -> None:
    missing = set(_REQUIRED) - config.keys()
    if missing:
        raise ReferenceRunError(f"missing required fields: {', '.join(sorted(missing))}")
    for key in ("sim_telarray_release", "hessio_decoder", "site", "production_version"):
        if not isinstance(config[key], str) or not config[key].strip():
            raise ReferenceRunError(f"{key} must be a nonempty string")
    for key in (
        "configuration_overrides",
        "seeds",
        "atmosphere_extinction",
        "coordinate_frames",
        "software_revisions",
    ):
        if not _object(config[key], key):
            raise ReferenceRunError(f"{key} must be nonempty")
    variants = config["telescope_variants"]
    if (
        not isinstance(variants, list)
        or not variants
        or any(not isinstance(item, str) or not item for item in variants)
        or len(set(variants)) != len(variants)
    ):
        raise ReferenceRunError("telescope_variants must be distinct model names")
    blocks = config["photon_blocks"]
    if not isinstance(blocks, list) or not blocks:
        raise ReferenceRunError("photon_blocks must list input files")
    for block in blocks:
        if not isinstance(block, str):
            raise ReferenceRunError("photon_blocks entries must be paths")
    commands = config["commands"]
    if not isinstance(commands, list) or not commands:
        raise ReferenceRunError("commands must list reference invocations")
    for command in commands:
        command = _object(command, "command")
        if set(command) != {"name", "argv", "cwd", "environment", "inputs"}:
            raise ReferenceRunError("command requires name, argv, cwd, environment, inputs")
        if not isinstance(command["name"], str) or not command["name"]:
            raise ReferenceRunError("command name must be nonempty")
        if (
            not isinstance(command["argv"], list)
            or not command["argv"]
            or any(not isinstance(arg, str) or not arg for arg in command["argv"])
        ):
            raise ReferenceRunError("command argv must be a nonempty string list")
        if not isinstance(command["cwd"], str) or not command["cwd"]:
            raise ReferenceRunError("command cwd must be a path")
        if not Path(command["cwd"]).is_absolute() or not Path(command["cwd"]).is_dir():
            raise ReferenceRunError(f"command cwd is missing: {command['cwd']}")
        _file(command["argv"][0], "command executable")
        if not isinstance(command["environment"], dict) or any(
            not isinstance(key, str)
            or not key
            or "=" in key
            or "\0" in key
            or not isinstance(value, str)
            or "\0" in value
            for key, value in command["environment"].items()
        ):
            raise ReferenceRunError("command environment must map strings to strings")
        if not isinstance(command["inputs"], list) or any(
            not isinstance(path, str) for path in command["inputs"]
        ):
            raise ReferenceRunError("command inputs must list file paths")


def freeze(config: dict[str, Any], model_root: Path) -> dict[str, Any]:
    """Resolve selected model versions and hash every declared reference input."""
    _validate(config)
    version = config["production_version"]
    try:
        models = {
            name: resolve_model(model_root, name, version) for name in config["telescope_variants"]
        }
    except ModelImportError as error:
        raise ReferenceRunError(str(error)) from error
    paths = set(config["photon_blocks"])
    for command in config["commands"]:
        paths.add(command["argv"][0])
        paths.update(command["inputs"])
    hashes = {path: sha256(_file(path, "reference input")) for path in sorted(paths)}
    return {
        "format": "obdeect.reference-run.v1",
        "configuration": copy.deepcopy(config),
        "model_checkout_revision": _git_revision(model_root),
        "models": models,
        "file_sha256": hashes,
    }


def verify(reference_run: dict[str, Any], model_root: Path) -> None:
    """Fail if a configured file, selected production, or checkout changed."""
    if reference_run.get("format") != "obdeect.reference-run.v1":
        raise ReferenceRunError("unsupported reference-run record format")
    config = _object(reference_run.get("configuration"), "configuration")
    expected = freeze(config, model_root)
    if reference_run != expected:
        raise ReferenceRunError(
            "reference-run record differs from current revisions or file hashes"
        )


def execute(
    reference_run: dict[str, Any], model_root: Path, output: Path, timeout_s: float = 3600.0
) -> dict[str, Any]:
    """Verify frozen inputs, then execute recorded argv and retain logs on failure.

    Commands run without a shell. The declared environment overrides the current
    process environment; the declared overrides are recorded for reproduction.
    Output directories must be new so previous evidence is never overwritten.
    """
    if not math.isfinite(timeout_s) or timeout_s <= 0:
        raise ReferenceRunError("timeout must be positive")
    verify(reference_run, model_root)
    output.mkdir(parents=True, exist_ok=False)
    results: dict[str, Any] = {
        "format": "obdeect.reference-execution.v1",
        "passed": False,
        "commands": [],
    }
    summary = output / "execution.json"
    for index, command in enumerate(reference_run["configuration"]["commands"]):
        stdout_path = output / f"{index:03d}.stdout.txt"
        stderr_path = output / f"{index:03d}.stderr.txt"
        # Preserve only declared environment in the public report; inherited
        # credentials are neither logged nor needed as evidence for optical physics.
        environment = {**os.environ, **command["environment"]}
        start = time.monotonic()
        record = {
            "name": command["name"],
            "argv": command["argv"],
            "cwd": command["cwd"],
            "environment": command["environment"],
            "returncode": None,
            "stdout": stdout_path.name,
            "stderr": stderr_path.name,
        }
        try:
            with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
                verify(reference_run, model_root)
                completed = subprocess.run(
                    command["argv"],
                    cwd=command["cwd"],
                    env=environment,
                    stdout=stdout,
                    stderr=stderr,
                    timeout=timeout_s,
                    check=False,
                )
            record["returncode"] = completed.returncode
            verify(reference_run, model_root)
        except (OSError, subprocess.TimeoutExpired, ReferenceRunError) as error:
            record["error"] = str(error)
        record["elapsed_s"] = time.monotonic() - start
        record["stdout_sha256"] = sha256(stdout_path)
        record["stderr_sha256"] = sha256(stderr_path)
        results["commands"].append(record)
        results["passed"] = record["returncode"] == 0 and "error" not in record
        summary.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
        if not results["passed"]:
            raise ReferenceRunError(f"recorded command {command['name']} failed; see {summary}")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--action", choices=("freeze", "verify", "execute"), required=True)
    parser.add_argument(
        "--input", type=Path, required=True, help="configuration or frozen reference-run JSON"
    )
    parser.add_argument("--model-root", type=Path, required=True)
    parser.add_argument(
        "--output", type=Path, help="file for freeze; new log directory for execute"
    )
    args = parser.parse_args()
    try:
        data = _object(json.loads(args.input.read_text(encoding="utf-8")), "JSON root")
        if args.action == "freeze":
            if args.output is None:
                raise ReferenceRunError("freeze requires --output")
            result = freeze(data, args.model_root)
            args.output.write_text(
                json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            print(f"Wrote frozen reference-run record to {args.output}")
        elif args.action == "execute":
            if args.output is None:
                raise ReferenceRunError("execute requires --output")
            execute(data, args.model_root, args.output)
            print(f"Executed frozen reference commands; logs in {args.output}")
        else:
            verify(data, args.model_root)
            print(f"Verified reference-run record {args.input}")
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise SystemExit(f"reference-run check failed: {error}") from error


if __name__ == "__main__":
    main()
