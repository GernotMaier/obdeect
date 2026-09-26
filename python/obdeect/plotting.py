#!/usr/bin/env python3
"""Optional visualization for CSV paths from the reference and CTAO reference tracers."""

import argparse
import csv
import json
import math
from pathlib import Path
from xml.sax.saxutils import escape


def _parse_path_row(path: Path, row: dict[str, str]) -> list[tuple[float, float, float]]:
    """Validate and return the ragged path vertices from one CSV record."""
    points = []
    try:
        count = int(row["point_count"])
        if not 1 <= count <= 4:
            raise ValueError("point_count must be between 1 and 4")
        for index in range(count):
            point = tuple(float(row[f"{axis}{index}_m"]) for axis in "xyz")
            if not all(math.isfinite(value) for value in point):
                raise ValueError("non-finite path vertex")
            points.append(point)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: invalid path vertices: {error}") from error
    return points


def read_paths(path: Path):
    """Yield terminal status and validated vertices from a trace CSV."""
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            yield row["status"], _parse_path_row(path, row)


def read_trace_rows(path: Path):
    """Yield validated path records with source and spectral metadata."""
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            yield row, row["status"], _parse_path_row(path, row)


def focal_plane_hits(path: Path):
    """Yield detected focal-plane ``(x, y, weight)`` samples from a trace CSV.

    The final ragged path vertex is the detector-surface intersection.  Modern
    CSV output supplies source weights and throughput; older files remain
    usable with unit weight.
    """
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("status") != "detected":
                continue
            try:
                point_index = int(row["point_count"]) - 1
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"{path}: detected row has an invalid point_count") from error
            if point_index < 0:
                raise ValueError(f"{path}: detected row has no path vertices")
            try:
                source_weight = float(row.get("source_weight", "1") or "1")
                throughput = float(row.get("throughput", "1") or "1")
                x = float(row[f"x{point_index}_m"])
                y = float(row[f"y{point_index}_m"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"{path}: detected row has invalid focal-plane data") from error
            if not all(math.isfinite(value) for value in (source_weight, throughput, x, y)):
                raise ValueError(f"{path}: detected row has non-finite focal-plane data")
            if source_weight < 0.0 or throughput < 0.0:
                raise ValueError(f"{path}: detected row has negative optical weight")
            yield (
                x,
                y,
                source_weight * throughput,
            )


TELESCOPE_NAMES = ("reference-mst", "LST", "MST", "SST", "SCT")


def draw_reference_telescope(axis, telescope: str):
    """Draw a side projection; reference models show optical surfaces only."""
    from matplotlib.patches import Rectangle

    if telescope == "reference-mst":
        mirror_x = [(-6.0 + 12.0 * index / 200.0) for index in range(201)]
        mirror_z = [9.75 - math.sqrt(9.75**2 - x**2) for x in mirror_x]
        axis.plot(mirror_x, mirror_z, color="0.15", linewidth=2.0, label="spherical primary")
        axis.add_patch(
            Rectangle(
                (-0.55, 4.625), 1.1, 0.5, color="black", alpha=0.25, label="camera body projection"
            )
        )
        for base_x, top_x in ((-4.2, -0.7), (0.0, 0.0), (4.2, 0.7)):
            axis.plot([base_x, top_x], [0.30, 4.625], color="0.3", linewidth=1)
        axis.plot([-2.0, 2.0], [4.875, 4.875], color="black", linewidth=0.7, label="focal screen")
        return
    # Only validated analytic reference outlines are drawn here.
    # SC geometry must come from a compiled-scene export.
    if telescope in ("SST", "SCT"):
        axis.text(
            0.02,
            0.98,
            "Unvalidated SC paths; geometry outline unavailable",
            transform=axis.transAxes,
            va="top",
        )
        return
    radius, camera_z, focal_radius = {
        "LST": (11.5, 28.0, 1.5),
        "MST": (6.0, 16.0, 1.0),
    }[telescope]
    primary_x = [(-radius + 2.0 * radius * index / 200.0) for index in range(201)]
    # An outline deliberately avoids reimplementing C++ prescription math.
    primary_z = (
        [x * x / (4 * camera_z) for x in primary_x]
        if telescope == "LST"
        else [2 * camera_z - math.sqrt((2 * camera_z) ** 2 - x * x) for x in primary_x]
    )
    axis.plot(primary_x, primary_z, color="0.15", linewidth=2.0, label="M1 reference aperture")
    axis.plot(
        [-focal_radius, focal_radius],
        [camera_z, camera_z],
        color="black",
        linewidth=0.9,
        label="focal plane",
    )


def draw_structure(plt, telescope: str, output: Path):
    """Render both side projections of the available telescope geometry."""
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for axis, horizontal in zip(axes, "xy", strict=True):
        draw_reference_telescope(axis, telescope)
        axis.set(xlabel=f"telescope {horizontal} [m]", ylabel="telescope z [m]")
        axis.set_aspect("equal", adjustable="box")
        axis.set_ylim(-0.5, {"reference-mst": 6.0, "LST": 30.0, "MST": 18.0}.get(telescope, 8.0))
    axes[0].legend(loc="best")
    description = (
        "reference mirror, camera and supports"
        if telescope == "reference-mst"
        else "reference optical outline"
    )
    figure.suptitle(f"{telescope} {description}")
    figure.tight_layout()
    figure.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def draw_compiled_structure(plt, scene_path: Path, output: Path):
    """Plot the explicit compiled primary facets and focal boundary.

    The plot consumes only geometry present in the compiled scene JSON. Any
    unavailable secondary, support, or obscuration geometry is listed in the
    figure annotation instead of being replaced by an illustrative outline.
    """
    scene = json.loads(scene_path.read_text(encoding="utf-8"))
    facets = scene.get("primary", {}).get("facets", [])
    if not facets:
        raise SystemExit("compiled scene contains no primary facets")
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for axis, coordinate in zip(axes, (0, 1), strict=True):
        for facet in facets:
            centre = facet.get("nominal_centre_m")
            if not isinstance(centre, list) or len(centre) != 3:
                continue
            axis.scatter(centre[coordinate], centre[2], s=8, c="tab:blue", alpha=0.65)
        camera = scene.get("camera", {})
        elements = camera.get("pixels", [])
        if elements:
            extent = max(
                math.hypot(float(item["centre_xy_m"][0]), float(item["centre_xy_m"][1]))
                for item in elements
            )
            focal = float(scene.get("focal_length_m", facets[0].get("focal_length_m", 0.0)))
            axis.plot(
                [-extent, extent],
                [focal, focal],
                color="tab:orange",
                linewidth=2,
                label="focal boundary",
            )
        axis.set(xlabel=f"compiled telescope {'xy'[coordinate]} [m]", ylabel="compiled z [m]")
        axis.set_aspect("equal", adjustable="box")
    axes[0].scatter([], [], s=8, c="tab:blue", label="mirror facet centres")
    axes[0].legend(loc="best")
    blockers = scene.get("report", {}).get("trace_blockers", [])
    subtitle = "compiled model geometry"
    if blockers:
        subtitle += "; unresolved: " + ", ".join(str(item) for item in blockers[:2])
    figure.suptitle(subtitle)
    figure.tight_layout()
    figure.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def draw_rays(plt, path: Path, telescope: str, output: Path, max_paths: int):
    """Render recorded paths in the two telescope side projections."""
    figure, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=True)
    colours = {
        "detected": "tab:blue",
        "blocked_camera": "tab:red",
        "blocked_mast": "tab:orange",
        "missed_primary": "0.5",
        "missed_screen": "tab:purple",
    }
    source_colours = {"star": "tab:blue", "illuminator": "tab:green", "laser": "tab:red"}
    source_seen = set()
    incidence = []
    for count, (row, status, points) in enumerate(read_trace_rows(path)):
        if count >= max_paths:
            break
        source = row.get("source_kind", "unknown")
        source_seen.add(source)
        colour = source_colours.get(source, colours.get(status, "black"))
        if row.get("incidence_focal_deg"):
            incidence.append(float(row["incidence_focal_deg"]))
        for axis, coordinate in zip(axes, (0, 1), strict=True):
            axis.plot(
                [point[coordinate] for point in points],
                [point[2] for point in points],
                color=colour,
                alpha=0.25,
                linewidth=0.7,
            )
    for axis, horizontal in zip(axes, "xy", strict=True):
        draw_reference_telescope(axis, telescope)
        axis.set(xlabel=f"telescope {horizontal} [m]", ylabel="telescope z [m]")
        axis.set_aspect("equal", adjustable="box")
        axis.set_ylim(-0.5, {"reference-mst": 8.0, "LST": 32.0, "MST": 20.0}.get(telescope, 8.0))
    axes[0].legend(loc="best")
    labels = ", ".join(sorted(source_seen)) or "unknown source"
    angle = f"; focal incidence mean {sum(incidence) / len(incidence):.3g}°" if incidence else ""
    figure.suptitle(f"{telescope} recorded ray paths — {labels}{angle}")
    figure.tight_layout()
    figure.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def draw_focal_plane(plt, path: Path, output: Path, bins: int, telescope: str):
    """Render surviving-photon intensity and its Cartesian projections."""
    import numpy as np

    samples = list(focal_plane_hits(path))
    if not samples:
        raise SystemExit("No surviving photons reached the focal plane in this CSV")
    values = np.asarray(samples, dtype=float)
    x, y, weights = values.T
    extent = max(float(np.max(np.abs(x))), float(np.max(np.abs(y))), 1.0e-6)
    extent *= 1.05

    figure = plt.figure(figsize=(9, 8))
    grid = figure.add_gridspec(
        2, 2, width_ratios=(4, 1.25), height_ratios=(1.25, 4), hspace=0.06, wspace=0.06
    )
    top = figure.add_subplot(grid[0, 0])
    image = figure.add_subplot(grid[1, 0])
    right = figure.add_subplot(grid[1, 1], sharey=image)
    histogram = image.hist2d(
        x,
        y,
        bins=bins,
        range=[[-extent, extent], [-extent, extent]],
        weights=weights,
        cmap="viridis",
    )
    figure.colorbar(histogram[3], ax=image, label="weighted surviving photons / bin")
    bin_edges = np.linspace(-extent, extent, bins + 1)
    top.hist(x, bins=bin_edges, weights=weights, color="tab:blue")
    right.hist(y, bins=bin_edges, weights=weights, orientation="horizontal", color="tab:blue")
    image.set(
        xlabel="focal-plane x [m]",
        ylabel="focal-plane y [m]",
        title=f"{telescope} focal plane",
    )
    image.set_aspect("equal", adjustable="box")
    top.set(ylabel="weighted surviving photons")
    right.set(xlabel="weighted surviving photons")
    top.tick_params(labelbottom=False)
    right.tick_params(labelleft=False)
    figure.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def draw_focal_plane_svg(path: Path, output: Path, bins: int, telescope: str):
    """Write a weighted PSF map and measured summary without plotting dependencies."""
    from obdeect.analysis import analyse_trace_csv

    samples = list(focal_plane_hits(path))
    result = analyse_trace_csv(path)
    extent = max(max(abs(x), abs(y)) for x, y, _ in samples) * 1.05
    extent = max(extent, 1e-6)
    histogram = [[0.0] * bins for _ in range(bins)]
    for x, y, weight in samples:
        ix = min(bins - 1, max(0, int((x + extent) * bins / (2 * extent))))
        iy = min(bins - 1, max(0, int((y + extent) * bins / (2 * extent))))
        histogram[iy][ix] += weight
    maximum = max(max(row) for row in histogram)
    left = 95
    top = 90
    size = 520
    cell = size / bins
    elements = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 720 760">',
        '<rect width="720" height="760" fill="white"/>',
        f"<title>{escape(telescope)} weighted focal-plane PSF</title>",
        f'<text x="95" y="38" font-size="24">{escape(telescope)} focal-plane PSF</text>',
        '<rect x="95" y="90" width="520" height="520" fill="#101827"/>',
    ]
    for iy, row in enumerate(histogram):
        for ix, weight in enumerate(row):
            if weight <= 0:
                continue
            fraction = math.sqrt(weight / maximum)
            red = round(25 + 230 * fraction)
            green = round(35 + 190 * fraction)
            blue = round(70 + 85 * (1 - fraction))
            elements.append(
                f'<rect x="{left + ix * cell:.4f}" y="{top + (bins - 1 - iy) * cell:.4f}" '
                f'width="{cell:.4f}" height="{cell:.4f}" fill="#{red:02x}{green:02x}{blue:02x}"/>'
            )
    centroid_x = left + (result.centroid_x_m + extent) * size / (2 * extent)
    centroid_y = top + (extent - result.centroid_y_m) * size / (2 * extent)
    elements.extend([
        f'<circle cx="{centroid_x:.4f}" cy="{centroid_y:.4f}" r="5" '
        'fill="none" stroke="white" stroke-width="2"/>',
        '<text x="95" y="645" font-size="16">x and y [m]; '
        "colour: weighted detected photons per bin</text>",
        f'<text x="95" y="676" font-size="16">Centroid: '
        f"({result.centroid_x_m:.6g}, {result.centroid_y_m:.6g}) m</text>",
        f'<text x="95" y="702" font-size="16">D80: {result.d80_m:.6g} m; '
        f"throughput: {result.optical_throughput:.4g}</text>",
        f'<text x="95" y="728" font-size="14">{bins} x {bins} bins; '
        f"axis extent: ±{extent:.6g} m; input: {escape(path.name)}</text>",
        "</svg>",
    ])
    output.write_text("\n".join(elements) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Plot an obdeect telescope reference outline and traced paths."
    )
    parser.add_argument(
        "paths", type=Path, nargs="?", help="trace CSV (required for rays and focal plane)"
    )
    parser.add_argument("--output", type=Path, default=Path("artificial_mst_paths.png"))
    parser.add_argument("--max-paths", type=int, default=300)
    parser.add_argument(
        "--view",
        choices=("structure", "compiled-structure", "rays", "focal-plane"),
        default="rays",
        help="structure outline, recorded ray paths, or weighted focal-plane image",
    )
    parser.add_argument(
        "--focal-plane",
        action="store_true",
        help=(
            "plot surviving focal-plane photons and weighted x/y projections "
            "instead of photon paths"
        ),
    )
    parser.add_argument("--bins", type=int, default=64, help="focal-plane histogram bins per axis")
    parser.add_argument(
        "--scene-json",
        type=Path,
        help="compiled scene JSON for --view compiled-structure",
    )
    parser.add_argument(
        "--telescope",
        choices=TELESCOPE_NAMES,
        default="reference-mst",
        help="outline only; does not modify trace coordinates",
    )
    args = parser.parse_args()
    if args.max_paths < 1 or args.bins < 1:
        parser.error("--max-paths and --bins must be positive")
    if args.paths is None and args.view not in {"structure", "compiled-structure"}:
        parser.error("paths CSV is required for rays and focal plane")
    if args.focal_plane and args.view != "rays":
        parser.error("--focal-plane cannot be combined with --view")
    if args.view == "compiled-structure" and args.scene_json is None:
        parser.error("--scene-json is required for --view compiled-structure")

    if args.focal_plane and args.output.suffix.lower() == ".svg":
        draw_focal_plane_svg(args.paths, args.output, args.bins, args.telescope)
        return

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise SystemExit(
            "Install optional visualization dependency: python -m pip install matplotlib"
        ) from error

    if args.view == "structure":
        draw_structure(plt, args.telescope, args.output)
        return
    if args.view == "compiled-structure":
        draw_compiled_structure(plt, args.scene_json, args.output)
        return
    if args.focal_plane or args.view == "focal-plane":
        draw_focal_plane(plt, args.paths, args.output, args.bins, args.telescope)
        return
    draw_rays(plt, args.paths, args.telescope, args.output, args.max_paths)


if __name__ == "__main__":
    main()
