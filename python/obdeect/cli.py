"""The ``obdeect`` command: compile optical models, trace sources, plot."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

from obdeect.formats import ModelError
from obdeect.model import STATUS_NAMES, OpticalModel
from obdeect.psf import summary


def _compile(args: argparse.Namespace) -> None:
    model = OpticalModel.from_simulation_models(
        args.simulation_models,
        args.model,
        args.version,
        configuration_seed=args.seed,
        zenith_deg=args.zenith_deg,
    )
    model.save(args.output)
    print(
        f"{args.model} {args.version}: optical model {model.sha256[:12]} written to {args.output}"
    )
    for note in model.document["notes"]:
        print(f"note: {note}")


def _trace(args: argparse.Namespace) -> None:
    model = OpticalModel.load(args.model)
    fx, fy = (math.radians(a) for a in args.field_deg)
    towards = np.array([math.tan(fx), math.tan(fy), 1.0])
    photons = model.star_photons(
        args.photons,
        seed=args.seed,
        towards_source=tuple(towards / np.linalg.norm(towards)),
        distance_m=args.distance_m,
        wavelength_nm=args.wavelength_nm,
    )
    arrivals = model.trace(
        photons, seed=args.seed, apply_responses=not args.geometry_only, threads=args.threads
    )
    if args.output:
        _write(args.output, photons, arrivals)
    print(
        json.dumps(
            {"model_sha256": model.sha256, **summary(arrivals, model.sampling_area_m2)}, indent=1
        )
    )


def _write(path: Path, photons: dict, arrivals: dict) -> None:
    columns = {"photon_id": photons["photon_id"], **arrivals}
    if path.suffix == ".npz":
        np.savez_compressed(path, status_names=np.array(STATUS_NAMES), **columns)
        return
    names, data = [], []
    for name, values in columns.items():
        if values.ndim == 2:
            names += [f"{name}_{axis}" for axis in "xyz"]
            data += list(values.T)
        else:
            names.append(name)
            data.append(values)
    np.savetxt(
        path, np.column_stack(data), delimiter=",", header=",".join(names), comments="", fmt="%.17g"
    )


def _plot(args: argparse.Namespace) -> None:
    from obdeect.plot import plot

    arrivals = dict(np.load(args.arrivals)) if args.arrivals else None
    plot(OpticalModel.load(args.model), arrivals, args.output)
    print(f"plot written to {args.output}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="obdeect", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    c = commands.add_parser("compile", help="compile an optical model from simulation-models")
    c.add_argument(
        "--simulation-models", type=Path, required=True, help="simulation-models checkout"
    )
    c.add_argument("--model", required=True, help="telescope model, e.g. LSTN-design")
    c.add_argument("--version", required=True, help="production version, e.g. 7.0.0")
    c.add_argument("--seed", type=int, required=True, help="seed of the random mirror errors")
    c.add_argument(
        "--zenith-deg", type=float, help="telescope zenith angle for zenith-dependent errors"
    )
    c.add_argument("--output", type=Path, required=True)
    c.set_defaults(run=_compile)

    t = commands.add_parser("trace", help="trace a star or point source")
    t.add_argument("--model", type=Path, required=True)
    t.add_argument(
        "--seed", type=int, required=True, help="seed of source sampling and mirror scatter"
    )
    t.add_argument("--photons", type=int, default=100000)
    t.add_argument(
        "--field-deg",
        type=float,
        nargs=2,
        default=(0.0, 0.0),
        metavar=("X", "Y"),
        help="source direction: field angles from the optical axis",
    )
    t.add_argument(
        "--distance-m", type=float, default=math.inf, help="source distance (default: star)"
    )
    t.add_argument("--wavelength-nm", type=float, default=400.0)
    t.add_argument("--geometry-only", action="store_true", help="do not apply optical responses")
    t.add_argument("--threads", type=int, default=0, help="0 uses every core")
    t.add_argument("--output", type=Path, help="per-photon results (.npz or .csv)")
    t.set_defaults(run=_trace)

    p = commands.add_parser("plot", help="plot the dish layout and the focal-plane image")
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--arrivals", type=Path, help="results of 'obdeect trace' (.npz)")
    p.add_argument("--output", type=Path, required=True)
    p.set_defaults(run=_plot)

    args = parser.parse_args(argv)
    try:
        args.run(args)
    except (ModelError, OSError) as error:
        print(f"obdeect: error: {error}", file=sys.stderr)
        return 1
    return 0
