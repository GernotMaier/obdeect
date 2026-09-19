#!/usr/bin/env python3
"""Cross-check the obdeect toy spherical primary against iactrace.

This deliberately compares only the shared analytic problem: parallel rays,
one continuous spherical primary and a plane at the paraxial focal distance.
Camera and mast shadows are disabled because this iactrace reference scene has
no matching obstruction definition. It is not a validation of the CTAO MST
configuration, whose 86 segmented facets, camera and structure the current
obdeect prototype does not yet implement.
"""

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

import jax

# The C++ kernel is binary64.  Make the reference binary64 too: iactrace
# defaults to float32 for throughput, which is a different numerical contract
# rather than an optical disagreement.
jax.config.update("jax_enable_x64", True)

# JAX must be configured before importing jax.numpy; these imports therefore
# intentionally follow the configuration call.
import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
from iactrace import Telescope  # noqa: E402
from iactrace.telescope.mirrors import spherical  # noqa: E402


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def points(rows: list[dict[str, str]], index: int) -> np.ndarray:
    return np.asarray([
        [float(row[f"x{index}_m"]), float(row[f"y{index}_m"]), float(row[f"z{index}_m"])]
        for row in rows
    ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--photons", type=int, default=256)
    parser.add_argument("--tolerance-m", type=float, default=1.0e-10)
    args = parser.parse_args()
    if args.photons <= 0 or args.tolerance_m <= 0.0:
        parser.error("--photons and --tolerance-m must be positive")

    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "paths.csv"
        subprocess.run(
            [
                args.executable,
                "--photons",
                str(args.photons),
                "--no-structure",
                "--output",
                str(output),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        rows = load_csv(output)

    if len(rows) != args.photons:
        raise AssertionError(f"obdeect emitted {len(rows)} rows for {args.photons} photons")

    origins = jnp.asarray(points(rows, 0))
    directions = jnp.broadcast_to(jnp.array([0.0, 0.0, -1.0]), origins.shape)
    mirror = spherical(
        position=(0.0, 0.0, 0.0),
        focal_length=4.875,
        radius=6.0,
        key=jax.random.key(0),
    )
    telescope = Telescope(
        mirror_groups=[mirror],
        camera_position=jnp.array([0.0, 0.0, 4.875]),
        camera_rotation=jnp.array([0.0, 0.0, 0.0]),
        name="matched_spherical_reference",
    )
    trace = telescope.trace(origins, directions, jnp.ones(args.photons), record_trajectory=True)
    if trace.trajectory is None:
        raise AssertionError("iactrace did not return the requested trajectory")
    reference = np.asarray(trace.trajectory.points)
    if reference.shape != (3, args.photons, 3):
        raise AssertionError(f"unexpected iactrace trajectory shape: {reference.shape}")

    mirror_delta = np.linalg.norm(points(rows, 1) - reference[1], axis=1)
    screen_delta = np.linalg.norm(points(rows, 2) - reference[2], axis=1)
    maximum = max(float(mirror_delta.max()), float(screen_delta.max()))
    print(f"iactrace version: {__import__('iactrace').__version__}")
    print(f"photons: {args.photons}")
    print(f"max mirror-point residual [m]: {mirror_delta.max():.9g}")
    print(f"max screen-point residual [m]: {screen_delta.max():.9g}")
    print(f"acceptance tolerance [m]: {args.tolerance_m:.9g}")
    if maximum > args.tolerance_m:
        raise SystemExit("FAILED: analytic spherical paths disagree with iactrace")
    print("PASS: matched analytic spherical trace agrees with iactrace")


if __name__ == "__main__":
    main()
