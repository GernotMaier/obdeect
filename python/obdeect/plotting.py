#!/usr/bin/env python3
"""Optional visualization for the dependency-free MST baseline CSV output."""

import argparse
import csv
import math
from pathlib import Path


def read_paths(path: Path):
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            points = []
            for index in range(int(row["point_count"])):
                points.append((float(row[f"x{index}_m"]), float(row[f"y{index}_m"]), float(row[f"z{index}_m"])))
            yield row["status"], points


def main():
    parser = argparse.ArgumentParser(description="Plot simple MST structure and sampled obdeect paths.")
    parser.add_argument("paths", type=Path)
    parser.add_argument("--output", type=Path, default=Path("toy_mst_paths.png"))
    parser.add_argument("--max-paths", type=int, default=300)
    args = parser.parse_args()

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise SystemExit("Install optional visualization dependency: python -m pip install matplotlib") from error

    fig, axis = plt.subplots(figsize=(8, 8))
    colours = {"detected": "tab:blue", "blocked_camera": "tab:red", "blocked_mast": "tab:orange",
               "missed_primary": "0.5", "missed_screen": "tab:purple"}
    for count, (status, points) in enumerate(read_paths(args.paths)):
        if count >= args.max_paths:
            break
        axis.plot([point[0] for point in points], [point[2] for point in points],
                  color=colours.get(status, "black"), alpha=0.20, linewidth=0.6)

    # The drawing is a diagnostic sketch only; path coordinates always come
    # from the C++ trace output and are never re-derived here.
    mirror_x = [(-6.0 + 12.0 * index / 200.0) for index in range(201)]
    mirror_z = [9.75 - math.sqrt(9.75**2 - x**2) for x in mirror_x]
    axis.plot(mirror_x, mirror_z, color="0.15", linewidth=2.0, label="spherical primary")
    for base_x in (-4.2, 0.0, 4.2):
        # x-z projections of the four finite support legs; the two y-plane
        # legs overlap at x=0 in this view.
        axis.plot([base_x, math.copysign(0.7, base_x) if base_x else 0.0], [0.30, 4.625],
                  color="0.25", linewidth=1.2)
    axis.add_patch(plt.Circle((0, 4.875), 0.55, color="black", alpha=0.25, label="camera body"))
    axis.axhline(4.875, color="black", linewidth=0.7, label="focal screen")
    axis.set(xlabel="telescope x [m]", ylabel="telescope z [m]", title="obdeect: 400-nm MST baseline paths")
    axis.set_aspect("equal", adjustable="box")
    axis.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(args.output, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
