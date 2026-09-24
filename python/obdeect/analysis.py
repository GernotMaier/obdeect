"""Weighted focal-plane analysis for recorded obdeect trace CSV files.

This module deliberately analyses the detector intersections recorded by the
tracer.  It does not re-trace rays or infer losses from a focal-plane image.
"""

import argparse
import csv
import json
import math
import subprocess
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class FocalPlaneHit:
    """One detected focal-surface intersection, in telescope-frame metres."""

    x_m: float
    y_m: float
    weight: float


@dataclass(frozen=True)
class PsfResult:
    """Weighted PSF and trace-loss summary for a single field direction.

    ``d80_m`` is the diameter of the smallest centroid-centred circle that
    contains at least 80 percent of detected optical weight.  It is an
    empirical weighted quantile, so no histogram binning is involved.
    """

    input_weight: float
    detected_weight: float
    optical_throughput: float
    detected_count: int
    centroid_x_m: float
    centroid_y_m: float
    r80_m: float
    d80_m: float
    terminal_count: dict[str, int]
    terminal_input_weight: dict[str, float]


def _finite_nonnegative(value: str | None, field: str, path: Path, row_number: int) -> float:
    try:
        result = float(value if value not in (None, "") else "1")
    except ValueError as error:
        raise ValueError(f"{path}:{row_number}: invalid {field}") from error
    if not math.isfinite(result) or result < 0.0:
        raise ValueError(f"{path}:{row_number}: {field} must be finite and non-negative")
    return result


def _final_hit(row: dict[str, str], path: Path, row_number: int) -> tuple[float, float]:
    try:
        point_count = int(row["point_count"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}:{row_number}: invalid point_count") from error
    if point_count < 1:
        raise ValueError(f"{path}:{row_number}: detected photon has no path vertices")
    point = point_count - 1
    try:
        x, y = float(row[f"x{point}_m"]), float(row[f"y{point}_m"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}:{row_number}: missing final focal-plane vertex") from error
    if not math.isfinite(x) or not math.isfinite(y):
        raise ValueError(f"{path}:{row_number}: focal-plane coordinates must be finite")
    return x, y


def analyse_trace_csv(path: Path) -> PsfResult:
    """Compute a weighted focal-plane PSF and loss closure from one trace CSV."""
    input_weight = 0.0
    detected_weight = 0.0
    terminal_count: dict[str, int] = defaultdict(int)
    terminal_input_weight: dict[str, float] = defaultdict(float)
    hits: list[FocalPlaneHit] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"status", "point_count"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: CSV must contain status and point_count columns")
        for row_number, row in enumerate(reader, start=2):
            status = row["status"]
            if not status:
                raise ValueError(f"{path}:{row_number}: terminal status is required")
            source_weight = _finite_nonnegative(
                row.get("source_weight"), "source_weight", path, row_number
            )
            input_weight += source_weight
            terminal_count[status] += 1
            terminal_input_weight[status] += source_weight
            if status != "detected":
                continue
            throughput = _finite_nonnegative(row.get("throughput"), "throughput", path, row_number)
            x, y = _final_hit(row, path, row_number)
            weight = source_weight * throughput
            detected_weight += weight
            if weight > 0.0:
                hits.append(FocalPlaneHit(x, y, weight))
    if input_weight <= 0.0:
        raise ValueError(f"{path}: total input weight must be positive")
    if not hits:
        raise ValueError(f"{path}: no detected focal-plane weight")
    centroid_x = sum(hit.x_m * hit.weight for hit in hits) / detected_weight
    centroid_y = sum(hit.y_m * hit.weight for hit in hits) / detected_weight
    radii = sorted(
        (math.hypot(hit.x_m - centroid_x, hit.y_m - centroid_y), hit.weight) for hit in hits
    )
    target_weight = 0.8 * detected_weight
    cumulative_weight = 0.0
    r80 = 0.0
    for radius, weight in radii:
        cumulative_weight += weight
        if cumulative_weight >= target_weight:
            r80 = radius
            break
    return PsfResult(
        input_weight=input_weight,
        detected_weight=detected_weight,
        optical_throughput=detected_weight / input_weight,
        detected_count=terminal_count["detected"],
        centroid_x_m=centroid_x,
        centroid_y_m=centroid_y,
        r80_m=r80,
        d80_m=2.0 * r80,
        terminal_count=dict(sorted(terminal_count.items())),
        terminal_input_weight=dict(sorted(terminal_input_weight.items())),
    )


def write_scan_csv(results: Iterable[tuple[float, PsfResult]], path: Path) -> None:
    """Write field-angle PSF metrics in a portable, plot-ready CSV table."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "field_x_deg",
                "input_weight",
                "detected_weight",
                "optical_throughput",
                "detected_count",
                "centroid_x_m",
                "centroid_y_m",
                "r80_m",
                "d80_m",
            ),
        )
        writer.writeheader()
        for field_x_deg, result in results:
            writer.writerow({
                "field_x_deg": field_x_deg,
                **{
                    key: value
                    for key, value in asdict(result).items()
                    if key not in {"terminal_count", "terminal_input_weight"}
                },
            })


def _plot_scan(results: list[tuple[float, PsfResult]], output: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError("Install matplotlib to request a scan plot") from error
    angle = [entry[0] for entry in results]
    d80 = [entry[1].d80_m for entry in results]
    throughput = [entry[1].optical_throughput for entry in results]
    figure, (psf_axis, throughput_axis) = plt.subplots(2, 1, sharex=True, figsize=(6, 6))
    psf_axis.plot(angle, d80, marker="o")
    psf_axis.set(ylabel="D80 [m]", title="Weighted focal-plane PSF scan")
    throughput_axis.plot(angle, throughput, marker="o")
    throughput_axis.set(xlabel="field x angle [deg]", ylabel="optical throughput")
    figure.tight_layout()
    figure.savefig(output, dpi=160)
    plt.close(figure)


def _derive(args: argparse.Namespace) -> int:
    if len(args.field_x_deg) not in (0, len(args.paths)):
        raise ValueError("--field-x-deg must be supplied once per input CSV")
    results = [analyse_trace_csv(path) for path in args.paths]
    payload = [asdict(result) for result in results]
    if args.field_x_deg:
        payload = [
            {"field_x_deg": angle, **result}
            for angle, result in zip(args.field_x_deg, payload, strict=True)
        ]
    args.output.write_text(
        json.dumps(payload[0] if len(payload) == 1 else payload, indent=2) + "\n",
        encoding="utf-8",
    )
    if args.scan_csv:
        if not args.field_x_deg:
            raise ValueError("--scan-csv requires --field-x-deg")
        write_scan_csv(list(zip(args.field_x_deg, results, strict=True)), args.scan_csv)
    if args.plot:
        if not args.field_x_deg:
            raise ValueError("--plot requires --field-x-deg")
        _plot_scan(list(zip(args.field_x_deg, results, strict=True)), args.plot)
    return 0


def _scan(args: argparse.Namespace) -> int:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results: list[tuple[float, PsfResult]] = []
    for index, field_x_deg in enumerate(args.field_x_deg):
        safe_angle = f"{field_x_deg:+.6f}".replace("+", "p").replace("-", "m")
        trace_path = args.output_dir / f"field_x_{index:04d}_{safe_angle}_deg.csv"
        command = [
            str(args.executable),
            "--photons",
            str(args.photons),
            "--field-x-deg",
            str(field_x_deg),
            "--field-y-deg",
            str(args.field_y_deg),
            "--output",
            str(trace_path),
            *args.extra_argument,
        ]
        subprocess.run(command, check=True)
        results.append((field_x_deg, analyse_trace_csv(trace_path)))
    write_scan_csv(results, args.output_dir / "psf_scan.csv")
    (args.output_dir / "psf_scan.json").write_text(
        json.dumps(
            [{"field_x_deg": angle, **asdict(result)} for angle, result in results], indent=2
        )
        + "\n",
        encoding="utf-8",
    )
    if args.plot:
        _plot_scan(results, args.plot)
    return 0


def main() -> int:
    """Entry point for deriving or simulating deterministic off-axis PSF scans."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    derive = commands.add_parser("derive", help="analyse existing trace CSV output")
    derive.add_argument("paths", type=Path, nargs="+")
    derive.add_argument("--field-x-deg", type=float, nargs="*", default=[])
    derive.add_argument("--output", type=Path, required=True)
    derive.add_argument("--scan-csv", type=Path)
    derive.add_argument("--plot", type=Path)
    derive.set_defaults(handler=_derive)
    scan = commands.add_parser("scan", help="trace and analyse a field-angle scan")
    scan.add_argument("--executable", type=Path, required=True)
    scan.add_argument("--field-x-deg", type=float, nargs="+", required=True)
    scan.add_argument("--field-y-deg", type=float, default=0.0)
    scan.add_argument("--photons", type=int, required=True)
    scan.add_argument("--output-dir", type=Path, required=True)
    scan.add_argument("--extra-argument", action="append", default=[])
    scan.add_argument("--plot", type=Path)
    scan.set_defaults(handler=_scan)
    args = parser.parse_args()
    if getattr(args, "photons", 1) <= 0:
        parser.error("--photons must be positive")
    for angle in (*args.field_x_deg, getattr(args, "field_y_deg", 0.0)):
        if not math.isfinite(angle):
            parser.error("field angles must be finite")
    try:
        return args.handler(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
