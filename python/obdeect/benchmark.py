"""Benchmark explicit optical commands only after matched photon equivalence.

Adapters must write normalized comparison CSVs. No simulator conventions or
kernel/I/O timings are inferred from executable names or total wall time.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import platform
import signal
import statistics
import subprocess
import sys
import time
from pathlib import Path

from obdeect.model_import import sha256
from obdeect.production_validation import (
    ProductionValidationError,
    compare,
    read_comparison_table,
    validate_comparison_fixture,
    validate_optical_model,
)


class BenchmarkError(ValueError):
    """Benchmark provenance, optical equivalence or execution failed."""


def _file(value, description):
    if not isinstance(value, str) or not Path(value).is_absolute() or not Path(value).is_file():
        raise BenchmarkError(f"{description} must name an existing absolute file")
    return Path(value)


def freeze(configuration: dict) -> dict:
    """Bind commands, optical model and source to immutable input hashes."""
    required = {
        "optical_model",
        "source",
        "fixture",
        "tolerances",
        "provenance",
        "repeats",
        "timeout_s",
        "variants",
        "engines",
        "inputs",
    }
    if not isinstance(configuration, dict) or set(configuration) != required:
        raise BenchmarkError(f"configuration requires exactly {sorted(required)}")
    config = copy.deepcopy(configuration)
    if (
        not isinstance(config["repeats"], int)
        or isinstance(config["repeats"], bool)
        or config["repeats"] < 10
    ):
        raise BenchmarkError("at least ten repeats are required")
    if (
        isinstance(config["timeout_s"], bool)
        or not isinstance(config["timeout_s"], (float, int))
        or not math.isfinite(config["timeout_s"])
        or config["timeout_s"] <= 0
    ):
        raise BenchmarkError("timeout_s must be finite and positive")
    if not isinstance(config["provenance"], dict) or not config["provenance"]:
        raise BenchmarkError("explicit provenance is required")
    model_path = _file(config["optical_model"], "optical_model")
    source_path = _file(config["source"], "source")
    model = json.loads(model_path.read_text())
    if not isinstance(model, dict):
        raise BenchmarkError("optical_model must contain a JSON object")
    try:
        validate_optical_model(model, "generic")
        validate_comparison_fixture(config["fixture"])
    except ProductionValidationError as error:
        raise BenchmarkError(str(error)) from error
    fixture = config["fixture"]
    if (
        not isinstance(fixture, dict)
        or fixture.get("optical_model_sha256") != model["optical_model_sha256"]
        or fixture.get("source_sha256") != sha256(source_path)
    ):
        raise BenchmarkError("fixture must bind the exact optical model and source hashes")
    tolerances = config["tolerances"]
    if (
        not isinstance(tolerances, dict)
        or set(tolerances)
        != {"focal_position_m", "path_length_m", "arrival_time_ns", "incidence_deg"}
        or any(
            isinstance(value, bool)
            or not isinstance(value, (float, int))
            or not math.isfinite(value)
            or value < 0
            for value in tolerances.values()
        )
    ):
        raise BenchmarkError("tolerances must declare finite nonnegative optical bounds")
    variants = config["variants"]
    if not isinstance(variants, list) or not variants:
        raise BenchmarkError("variants must declare thread/block configurations")
    seen = set()
    for variant in variants:
        if (
            not isinstance(variant, dict)
            or set(variant) != {"threads", "block_size"}
            or any(
                not isinstance(v, int) or isinstance(v, bool) or v <= 0 for v in variant.values()
            )
        ):
            raise BenchmarkError("variants require positive integer threads and block_size")
        key = (variant["threads"], variant["block_size"])
        if key in seen:
            raise BenchmarkError("duplicate thread/block variant")
        seen.add(key)
    if len({v[0] for v in seen}) < 2 or len({v[1] for v in seen}) < 2:
        raise BenchmarkError("reproducibility requires multiple thread counts and block sizes")
    if not isinstance(config["inputs"], list):
        raise BenchmarkError("inputs must list all additional command input files")
    paths = {str(model_path), str(source_path)}
    paths.update(str(_file(value, "input")) for value in config["inputs"])
    engines = config["engines"]
    if not isinstance(engines, dict) or set(engines) != {"baseline", "candidate"}:
        raise BenchmarkError("engines must declare baseline and candidate")
    for name, engine in engines.items():
        if not isinstance(engine, dict) or set(engine) != {
            "argv",
            "cwd",
            "environment",
            "timing_metrics",
        }:
            raise BenchmarkError("engine requires argv, cwd, environment and timing_metrics")
        argv = engine["argv"]
        if (
            not isinstance(argv, list)
            or not argv
            or any(not isinstance(arg, str) or not arg for arg in argv)
        ):
            raise BenchmarkError("engine argv must be nonempty strings")
        if name == "candidate" and not all(
            any(token in arg for arg in argv) for token in ("{threads}", "{block_size}")
        ):
            raise BenchmarkError("candidate argv must apply {threads} and {block_size}")
        paths.add(str(_file(argv[0], "engine executable")))
        if not any("{output}" in arg for arg in argv):
            raise BenchmarkError("engine argv must use {output} for a new comparison CSV")
        if not isinstance(engine["timing_metrics"], bool) or (
            engine["timing_metrics"] and not any("{metrics}" in arg for arg in argv)
        ):
            raise BenchmarkError("timing_metrics requires a {metrics} output argument")
        if (
            not isinstance(engine["cwd"], str)
            or not Path(engine["cwd"]).is_absolute()
            or not Path(engine["cwd"]).is_dir()
        ):
            raise BenchmarkError("engine cwd must be an existing absolute directory")
        env = engine["environment"]
        if not isinstance(env, dict) or any(
            not isinstance(k, str)
            or not k
            or "=" in k
            or "\0" in k
            or not isinstance(v, str)
            or "\0" in v
            for k, v in env.items()
        ):
            raise BenchmarkError("engine environment must map valid strings")
    return {
        "format": "obdeect.benchmark-input.v1",
        "configuration": config,
        "file_sha256": {path: sha256(Path(path)) for path in sorted(paths)},
    }


def verify(record: dict) -> None:
    """Reject altered executables, model, sources or adapter inputs."""
    if (
        not isinstance(record, dict)
        or record.get("format") != "obdeect.benchmark-input.v1"
        or freeze(record.get("configuration")) != record
    ):
        raise BenchmarkError("benchmark inputs differ from the frozen record")


def _run(argv, cwd, environment, directory, timeout_s):
    """Measure one direct child using OS wait4 resource accounting.

    RSS is the maximum resident set of the measured process, not simultaneous
    summed memory across a process tree. macOS reports bytes; Linux reports KiB.
    """
    if not hasattr(os, "wait4") or sys.platform not in {"darwin", "linux"}:
        raise BenchmarkError("per-process RSS measurement requires macOS or Linux wait4")
    with (
        (directory / "stdout.txt").open("wb") as stdout,
        (directory / "stderr.txt").open("wb") as stderr,
    ):
        start = time.perf_counter()
        process = subprocess.Popen(
            argv,
            cwd=cwd,
            env={**os.environ, **environment},
            stdout=stdout,
            stderr=stderr,
            start_new_session=True,
        )
        timed_out = False
        while True:
            pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid:
                process.returncode = os.waitstatus_to_exitcode(status)
                break
            if not timed_out and time.perf_counter() - start > timeout_s:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass  # Child exited between wait4 and the timeout check.
                timed_out = True
            time.sleep(0.001)
        elapsed = time.perf_counter() - start
    measurement = {
        "elapsed_s": elapsed,
        "user_cpu_s": usage.ru_utime,
        "system_cpu_s": usage.ru_stime,
        "peak_process_rss_bytes": int(usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)),
        "returncode": process.returncode,
        "timed_out": timed_out,
    }
    if timed_out or process.returncode:
        raise BenchmarkError(f"command failed: {measurement}; logs in {directory}")
    return measurement


def _spread(values):
    return {
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
        "stdev": statistics.stdev(values),
    }


def execute(record: dict, output: Path) -> dict:
    """Validate every timed result, then publish speedup and resource summaries."""
    verify(record)
    output.mkdir(parents=True, exist_ok=False)
    config = record["configuration"]
    report = {
        "format": "obdeect.benchmark-report.v1",
        "passed": False,
        "runs": [],
        "platform": platform.platform(),
        "python": sys.version,
        "input_record": record,
        "memory_semantics": "peak RSS of measured process; not summed process-tree RSS",
        "timing_semantics": (
            "process wall time includes startup and output I/O; "
            "kernel/I/O only when explicitly instrumented"
        ),
    }
    summary_path = output / "benchmark.json"
    references = {}
    zero = dict.fromkeys(
        ("focal_position_m", "path_length_m", "arrival_time_ns", "incidence_deg"), 0.0
    )
    try:
        for variant_index, variant in enumerate(config["variants"]):
            for repeat in range(config["repeats"]):
                pair = {}
                # Alternate execution order to reduce a consistent first-run bias.
                names = ("baseline", "candidate") if repeat % 2 == 0 else ("candidate", "baseline")
                for name in names:
                    verify(record)
                    engine = config["engines"][name]
                    directory = output / f"v{variant_index:03d}-r{repeat:03d}-{name}"
                    directory.mkdir()
                    table = directory / "arrivals.csv"
                    metrics = directory / "timing.json"
                    tokens = {
                        "{output}": str(table.resolve()),
                        "{metrics}": str(metrics.resolve()),
                        "{threads}": str(variant["threads"]),
                        "{block_size}": str(variant["block_size"]),
                    }
                    argv = []
                    for arg in engine["argv"]:
                        for token, value in tokens.items():
                            arg = arg.replace(token, value)
                        argv.append(arg)
                    measurement = _run(
                        argv, engine["cwd"], engine["environment"], directory, config["timeout_s"]
                    )
                    verify(record)
                    rows = read_comparison_table(table)
                    pair[name] = rows
                    prior = references.get(name)
                    if prior is None:
                        references[name] = rows
                    elif rows != prior:
                        # Unlike cross-engine equivalence, reproducibility also
                        # covers loss vertices and every normalized diagnostic.
                        compare(rows, prior, zero, config["fixture"])
                        raise BenchmarkError(
                            f"{name} normalized outputs are not exactly reproducible"
                        )
                    measurement.update(
                        engine=name,
                        variant=variant,
                        repeat=repeat,
                        argv=argv,
                        output_sha256=sha256(table),
                    )
                    if engine["timing_metrics"]:
                        instrumented = json.loads(metrics.read_text())
                        if (
                            not isinstance(instrumented, dict)
                            or set(instrumented) != {"kernel_s", "io_s"}
                            or any(
                                isinstance(v, bool)
                                or not isinstance(v, (float, int))
                                or not math.isfinite(v)
                                or v < 0
                                for v in instrumented.values()
                            )
                            or sum(instrumented.values()) > measurement["elapsed_s"]
                        ):
                            raise BenchmarkError(
                                "timing metrics require nonnegative kernel_s/io_s within wall time"
                            )
                        measurement.update(instrumented)
                    report["runs"].append(measurement)
                residuals = compare(
                    pair["candidate"], pair["baseline"], config["tolerances"], config["fixture"]
                )
                report["runs"][-1]["matched_residuals"] = residuals
        aggregates = {}
        for variant in config["variants"]:
            key = f"threads={variant['threads']},block_size={variant['block_size']}"
            aggregates[key] = {}
            for name in ("baseline", "candidate"):
                runs = [
                    run
                    for run in report["runs"]
                    if run["engine"] == name and run["variant"] == variant
                ]
                fields = ["elapsed_s", "peak_process_rss_bytes", "user_cpu_s", "system_cpu_s"]
                if config["engines"][name]["timing_metrics"]:
                    fields += ["kernel_s", "io_s"]
                aggregates[key][name] = {
                    field: _spread([run[field] for run in runs]) for field in fields
                }
            aggregates[key]["wall_speedup"] = (
                aggregates[key]["baseline"]["elapsed_s"]["median"]
                / aggregates[key]["candidate"]["elapsed_s"]["median"]
            )
            if all(config["engines"][name]["timing_metrics"] for name in ("baseline", "candidate")):
                candidate = aggregates[key]["candidate"]["kernel_s"]["median"]
                if candidate > 0:
                    aggregates[key]["kernel_speedup"] = (
                        aggregates[key]["baseline"]["kernel_s"]["median"] / candidate
                    )
        report.update(
            passed=True,
            summaries=aggregates,
            reproducibility=(
                "exact normalized per-photon comparison across repetitions "
                "and thread/block variants"
            ),
        )
    except (OSError, ValueError) as error:
        report["error"] = str(error)
        raise BenchmarkError(f"benchmark failed: {error}; see {summary_path}") from error
    finally:
        summary_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--action", choices=("freeze", "verify", "execute"), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        data = json.loads(args.input.read_text())
        if args.action == "freeze":
            if args.output is None:
                raise BenchmarkError("freeze requires --output")
            args.output.write_text(json.dumps(freeze(data), indent=2, sort_keys=True) + "\n")
        elif args.action == "verify":
            verify(data)
        else:
            if args.output is None:
                raise BenchmarkError("execute requires a new --output directory")
            execute(data, args.output)
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error


if __name__ == "__main__":
    main()
