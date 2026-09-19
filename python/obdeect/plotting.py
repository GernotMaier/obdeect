#!/usr/bin/env python3
"""Optional visualization for CSV paths from the toy and CTAO reference tracers."""

import argparse
import csv
import math
from pathlib import Path


def read_paths(path: Path):
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            points = []
            for index in range(int(row["point_count"])):
                points.append((
                    float(row[f"x{index}_m"]),
                    float(row[f"y{index}_m"]),
                    float(row[f"z{index}_m"]),
                ))
            yield row["status"], points


TELESCOPE_NAMES = ("toy-mst", "LST", "MST", "SST", "SCT")


def draw_reference_telescope(axis, telescope: str):
    """Draw only a labelled diagnostic outline; rays always originate in CSV."""
    from matplotlib.patches import Rectangle

    if telescope == "toy-mst":
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
        axis.axhline(4.875, color="black", linewidth=0.7, label="focal screen")
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
    specification = {
        "LST": (11.5, 28.0, None),
        "MST": (6.0, 16.0, None),
        "SST": (2.1205, 2.15, 3.1084),
        "SCT": (4.8319, 5.5863, 1.5),
    }[telescope]
    radius, camera_z, secondary_z = specification
    primary_x = [(-radius + 2.0 * radius * index / 200.0) for index in range(201)]
    # An outline deliberately avoids reimplementing C++ prescription math.
    primary_z = (
        [x * x / (4 * camera_z) for x in primary_x]
        if telescope == "LST"
        else [2 * camera_z - math.sqrt((2 * camera_z) ** 2 - x * x) for x in primary_x]
    )
    axis.plot(primary_x, primary_z, color="0.15", linewidth=2.0, label="M1 reference aperture")
    if secondary_z is not None:
        secondary_radius = 0.90 if telescope == "SST" else 2.71
        axis.plot(
            [-secondary_radius, secondary_radius],
            [secondary_z, secondary_z],
            color="0.3",
            linewidth=2.0,
            label="M2 reference aperture",
        )
    axis.axhline(camera_z, color="black", linewidth=0.9, label="focal plane")


def main():
    parser = argparse.ArgumentParser(
        description="Plot an obdeect telescope reference outline and traced paths."
    )
    parser.add_argument("paths", type=Path)
    parser.add_argument("--output", type=Path, default=Path("toy_mst_paths.png"))
    parser.add_argument("--max-paths", type=int, default=300)
    parser.add_argument(
        "--telescope",
        choices=TELESCOPE_NAMES,
        default="toy-mst",
        help="outline only; does not modify trace coordinates",
    )
    args = parser.parse_args()

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise SystemExit(
            "Install optional visualization dependency: python -m pip install matplotlib"
        ) from error

    fig, axis = plt.subplots(figsize=(8, 8))
    colours = {
        "detected": "tab:blue",
        "blocked_camera": "tab:red",
        "blocked_mast": "tab:orange",
        "missed_primary": "0.5",
        "missed_screen": "tab:purple",
    }
    for count, (status, points) in enumerate(read_paths(args.paths)):
        if count >= args.max_paths:
            break
        axis.plot(
            [point[0] for point in points],
            [point[2] for point in points],
            color=colours.get(status, "black"),
            alpha=0.20,
            linewidth=0.6,
        )

    draw_reference_telescope(axis, args.telescope)
    axis.set(
        xlabel="telescope x [m]",
        ylabel="telescope z [m]",
        title=f"obdeect: 400-nm {args.telescope} reference paths",
    )
    axis.set_aspect("equal", adjustable="box")
    axis.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(args.output, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
