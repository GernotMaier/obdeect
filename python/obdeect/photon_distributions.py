"""Summarise and render weighted photon distributions in focal space and time.

The input is the strict ``obdeect-arrival-v1`` optical-boundary contract.  A
time value is the vacuum geometric arrival time (emission time plus geometric
path divided by ``c``); it is intentionally labelled as such, rather than
being presented as a detector/electronics timestamp.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from xml.sax.saxutils import escape

from obdeect.result_contract import ArrivalContractError, OpticalArrival, read_arrivals

_SPEED_OF_LIGHT_M_PER_NS = 0.299792458


@dataclass(frozen=True)
class WeightedHistogram:
    """One fixed-width histogram with SI edges and optical weights."""

    lower: float
    upper: float
    bins: tuple[float, ...]


@dataclass(frozen=True)
class PhotonDistributions:
    """Detected-photon focal map and vacuum-geometric time distribution."""

    detected_count: int
    detected_weight: float
    focal_x_m: WeightedHistogram
    focal_y_m: WeightedHistogram
    arrival_time_ns: WeightedHistogram
    focal_plane_bins: tuple[tuple[float, ...], ...]


def _bounds(values: list[float]) -> tuple[float, float]:
    lower, upper = min(values), max(values)
    if lower == upper:
        padding = max(abs(lower) * 0.05, 1.0e-9)
        return lower - padding, upper + padding
    padding = 0.05 * (upper - lower)
    return lower - padding, upper + padding


def _histogram(values: list[float], weights: list[float], bins: int) -> WeightedHistogram:
    lower, upper = _bounds(values)
    result = [0.0] * bins
    scale = bins / (upper - lower)
    for value, weight in zip(values, weights, strict=True):
        index = min(bins - 1, max(0, int((value - lower) * scale)))
        result[index] += weight
    return WeightedHistogram(lower, upper, tuple(result))


def photon_distributions(arrivals: list[OpticalArrival], bins: int = 64) -> PhotonDistributions:
    """Calculate weighted focal-plane and arrival-time distributions.

    Only detected photons have a focal position and take part in the spatial
    map.  Their time includes their source emission time and traced geometric
    path, preserving pulses and source-time offsets in the arrival contract.
    """
    if bins < 2 or bins > 512:
        raise ValueError("bins must be between 2 and 512")
    detected = [
        arrival for arrival in arrivals if arrival.detected and arrival.optical_weight > 0.0
    ]
    if not detected:
        raise ValueError("arrival table contains no detected optical weight")
    x = [float(arrival.focal_x_m) for arrival in detected]
    y = [float(arrival.focal_y_m) for arrival in detected]
    times = [
        arrival.emission_time_ns + arrival.path_length_m / _SPEED_OF_LIGHT_M_PER_NS
        for arrival in detected
    ]
    weights = [arrival.optical_weight for arrival in detected]
    x_histogram, y_histogram, time_histogram = (
        _histogram(values, weights, bins) for values in (x, y, times)
    )
    plane = [[0.0] * bins for _ in range(bins)]
    x_scale = bins / (x_histogram.upper - x_histogram.lower)
    y_scale = bins / (y_histogram.upper - y_histogram.lower)
    for x_value, y_value, weight in zip(x, y, weights, strict=True):
        x_index = min(bins - 1, max(0, int((x_value - x_histogram.lower) * x_scale)))
        y_index = min(bins - 1, max(0, int((y_value - y_histogram.lower) * y_scale)))
        plane[y_index][x_index] += weight
    return PhotonDistributions(
        detected_count=len(detected),
        detected_weight=sum(weights),
        focal_x_m=x_histogram,
        focal_y_m=y_histogram,
        arrival_time_ns=time_histogram,
        focal_plane_bins=tuple(tuple(row) for row in plane),
    )


def _histogram_svg(
    histogram: WeightedHistogram, *, x: int, y: int, width: int, height: int, title: str
) -> str:
    maximum = max(histogram.bins) or 1.0
    bar_width = width / len(histogram.bins)
    bars = []
    for index, weight in enumerate(histogram.bins):
        bar_height = height * weight / maximum
        bars.append(
            f'<rect x="{x + index * bar_width:.3f}" y="{y + height - bar_height:.3f}" '
            f'width="{bar_width:.3f}" height="{bar_height:.3f}" fill="#276fbf"/>'
        )
    return "\n".join([
        f'<text x="{x}" y="{y - 10}" font-size="16">{escape(title)}</text>',
        f'<rect x="{x}" y="{y}" width="{width}" height="{height}" fill="white" stroke="#334155"/>',
        *bars,
        f'<text x="{x}" y="{y + height + 20}" font-size="12">{histogram.lower:.6g}</text>',
        f'<text x="{x + width - 45}" y="{y + height + 20}" font-size="12">'
        f"{histogram.upper:.6g}</text>",
    ])


def write_distribution_svg(distributions: PhotonDistributions, output: Path) -> None:
    """Render a dependency-free, weighted focal map and time histogram SVG."""
    bins = len(distributions.focal_x_m.bins)
    maximum = max(max(row) for row in distributions.focal_plane_bins) or 1.0
    map_x, map_y, size = 65, 65, 420
    cell = size / bins
    cells = []
    for y_index, row in enumerate(distributions.focal_plane_bins):
        for x_index, weight in enumerate(row):
            if weight <= 0.0:
                continue
            intensity = math.sqrt(weight / maximum)
            red = round(25 + 220 * intensity)
            green = round(35 + 170 * intensity)
            blue = round(70 + 80 * (1 - intensity))
            cells.append(
                f'<rect x="{map_x + x_index * cell:.3f}" '
                f'y="{map_y + (bins - 1 - y_index) * cell:.3f}" '
                f'width="{cell:.3f}" height="{cell:.3f}" '
                f'fill="#{red:02x}{green:02x}{blue:02x}"/>'
            )
    svg = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1120 610">',
        '<rect width="1120" height="610" fill="white"/>',
        "<title>Weighted detected-photon distributions</title>",
        '<text x="65" y="32" font-size="24">Detected photon distributions</text>',
        f'<text x="65" y="53" font-size="14">{distributions.detected_count} photons; '
        f"weighted optical total {distributions.detected_weight:.6g}</text>",
        f'<text x="65" y="505" font-size="14">focal x [m]: '
        f"{distributions.focal_x_m.lower:.6g} to {distributions.focal_x_m.upper:.6g}</text>",
        f'<text x="65" y="527" font-size="14">focal y [m]: '
        f"{distributions.focal_y_m.lower:.6g} to {distributions.focal_y_m.upper:.6g}</text>",
        f'<rect x="{map_x}" y="{map_y}" width="{size}" height="{size}" fill="#101827"/>',
        *cells,
        '<text x="65" y="565" font-size="16">weighted focal-plane distribution [m]</text>',
        _histogram_svg(
            distributions.arrival_time_ns,
            x=570,
            y=105,
            width=470,
            height=160,
            title="vacuum geometric arrival time [ns]",
        ),
        _histogram_svg(
            distributions.focal_x_m,
            x=570,
            y=355,
            width=470,
            height=120,
            title="focal x distribution [m]",
        ),
        "</svg>",
    ]
    try:
        output.write_text("\n".join(svg) + "\n", encoding="utf-8")
    except OSError as error:
        raise ArrivalContractError(f"cannot write {output}: {error}") from error


def main() -> None:
    """Create a self-contained photon-space/time SVG from arrival records."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="obdeect-arrival-v1 CSV")
    parser.add_argument("--output", type=Path, required=True, help="distribution SVG")
    parser.add_argument("--bins", type=int, default=64, help="bins per distribution axis (2-512)")
    args = parser.parse_args()
    try:
        distributions = photon_distributions(read_arrivals(args.input), args.bins)
        write_distribution_svg(distributions, args.output)
    except (ArrivalContractError, ValueError) as error:
        raise SystemExit(f"photon distribution failed: {error}") from error
    print(
        f"Wrote weighted photon distributions for {distributions.detected_count} detections "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
