#!/usr/bin/env python3
"""Run reproducible analytic/source examples and render their recorded paths."""

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path


def validate_rows(
    name: str, rows: list[dict[str, str]], photons: int, require_detected: bool
) -> dict[str, int]:
    if len(rows) != photons or any(
        None in row or any(value is None for value in row.values()) for row in rows
    ):
        raise ValueError(f"{name}: malformed CSV or missing photons")
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
        count = int(row["point_count"])
        if not 1 <= count <= 4:
            raise ValueError(f"{name}: invalid path length")
    if require_detected and not counts.get("detected"):
        raise ValueError(f"{name}: no focal-plane hits")
    return counts


def run_examples(build: Path, output: Path, photons: int) -> None:
    root = Path(__file__).resolve().parents[1]
    output.mkdir(parents=True, exist_ok=True)
    suffix = ".exe" if os.name == "nt" else ""
    cases = [
        ("toy_star", "obdeect_toy", ["--source", "star"], "toy-mst", True),
        (
            "toy_flasher",
            "obdeect_toy",
            ["--source", "illuminator", "--distance-m", "50"],
            "toy-mst",
            False,
        ),
        (
            "toy_laser",
            "obdeect_toy",
            ["--source", "laser", "--divergence-deg", "0.1"],
            "toy-mst",
            True,
        ),
        ("lst_on_axis", "obdeect_ctao", ["--telescope", "LST"], "LST", True),
        (
            "lst_off_axis",
            "obdeect_ctao",
            ["--telescope", "LST", "--field-x-deg", "0.5"],
            "LST",
            True,
        ),
        ("mst_sphere", "obdeect_ctao", ["--telescope", "MST"], "MST", True),
    ]
    summary = {}
    for name, executable, flags, telescope, require_detected in cases:
        csv_path = output / f"{name}.csv"
        command = [
            str((build / (executable + suffix)).resolve()),
            *flags,
            "--photons",
            str(photons),
            "--output",
            str(csv_path),
        ]
        subprocess.run(command, check=True)
        with csv_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        counts = validate_rows(name, rows, photons, require_detected)
        subprocess.run(
            [
                sys.executable,
                str(root / "python/obdeect/plotting.py"),
                str(csv_path),
                "--telescope",
                telescope,
                "--output",
                str(output / f"{name}.png"),
            ],
            check=True,
        )
        summary[name] = {"status_counts": counts, "photons": photons}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, default=Path("build"))
    parser.add_argument("--output", type=Path, default=Path("out/examples"))
    parser.add_argument("--photons", type=int, default=1000)
    args = parser.parse_args()
    if args.photons <= 0:
        parser.error("--photons must be positive")
    run_examples(args.build, args.output, args.photons)
