"""Launch the C++ executables shipped in the obdeect wheel."""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def executable_path(name: str = "obdeect-simtools-raytrace") -> Path:
    """Return the native executable from a wheel or editable package installation.

    Parameters
    ----------
    name : str
        Executable name, without a Windows ``.exe`` suffix.

    Raises
    ------
    FileNotFoundError
        If the wheel does not contain the requested executable.
    """

    from obdeect import __path__ as package_paths

    # Editable installs can put Python sources and native files in separate
    # package locations. Search the package paths in their declared order.
    native_dirs = [Path(location) / "_native" for location in package_paths]
    candidates = []
    for native_dir in native_dirs:
        candidates.append(native_dir / name)
        if os.name == "nt" and not name.endswith(".exe"):
            candidates.append(native_dir / f"{name}.exe")
    for candidate in candidates:
        if candidate.is_file():
            if os.name != "nt" and not os.access(candidate, os.X_OK):
                raise PermissionError(f"Packaged obdeect executable is not executable: {candidate}")
            return candidate
    raise FileNotFoundError(
        f"Packaged obdeect executable {name!r} was not found in "
        f"{', '.join(str(directory) for directory in native_dirs)}. "
        "Reinstall obdeect-dev with a wheel for this platform."
    )


def _run(name: str) -> int:
    """Forward command-line arguments to a packaged C++ executable."""

    command = [str(executable_path(name)), *sys.argv[1:]]
    if os.name == "nt":
        return subprocess.call(command)
    os.execv(command[0], command)
    return 1


def demo_mst_main() -> int:
    """Run the packaged developer MST demonstration."""

    return _run("obdeect-demo-mst")


def analytic_optics_main() -> int:
    """Run the packaged analytic-optics developer diagnostic."""

    return _run("obdeect-analytic-optics")


def simtools_raytrace_main() -> int:
    """Run the packaged simtools-compatible ray tracer."""
    arguments = sys.argv[1:]
    if "--photon-input-frame" not in arguments:
        return _run("obdeect-simtools-raytrace")

    def take(name, default=None):
        if name not in arguments:
            return default
        index = arguments.index(name)
        if index + 1 >= len(arguments):
            raise ValueError(f"missing value for {name}")
        value = arguments[index + 1]
        del arguments[index : index + 2]
        return value

    frame = take("--photon-input-frame")
    if frame == "telescope":
        return subprocess.call([str(executable_path()), *arguments])
    if frame != "ground":
        raise ValueError("photon input frame must be ground or telescope")
    zenith, azimuth = take("--pointing-zenith-deg"), take("--pointing-azimuth-deg")
    if zenith is None or azimuth is None:
        raise ValueError("ground photon input requires explicit pointing zenith and azimuth")
    seed = int(take("--pointing-seed", "0"))
    from obdeect.observing_geometry import prepare_observing_model
    from obdeect.telescope_frame import (
        prepare_ground_photons,
        resolve_telescope_frame,
        single_telescope_identity,
    )

    if "--optical-model" not in arguments or "--photon-input" not in arguments:
        raise ValueError("ground replay requires an optical model and photon input")
    model_index, input_index = (
        arguments.index("--optical-model") + 1,
        arguments.index("--photon-input") + 1,
    )
    original_model = Path(arguments[model_index])
    original_input = Path(arguments[input_index])
    for name in ("--output", "--photon-output", "--interactions-output"):
        if name in arguments:
            output = Path(arguments[arguments.index(name) + 1]).resolve()
            if output in (original_input.resolve(), original_model.resolve()):
                raise ValueError("ground replay output must not overwrite its input")
    context = json.loads(original_model.read_text())["telescope_frame"]
    if context["axis_origin_m"] is None:
        with original_input.open(newline="") as stream:
            first = next(csv.DictReader(stream), {})
        if not all(f"telescope_{axis}_m" in first for axis in "xyz"):
            raise ValueError("ground replay requires elevation-axis coordinates in model or input")
        context["axis_origin_m"] = [float(first[f"telescope_{axis}_m"]) for axis in "xyz"]
    telescope_id = single_telescope_identity(original_input)
    pointing = resolve_telescope_frame(
        context,
        float(azimuth),
        float(zenith),
        seed=seed,
        telescope_id=telescope_id,
    )
    with tempfile.TemporaryDirectory(prefix="obdeect-ground-replay-") as temporary:
        directory = Path(temporary)
        prepared_model = prepare_observing_model(
            original_model, directory / "optical-model.json", pointing.known_zenith_deg
        )
        prepared_input = directory / "photons.csv"
        prepare_ground_photons(
            original_input,
            prepared_input,
            context,
            float(azimuth),
            float(zenith),
            seed=seed,
            telescope_id=telescope_id,
        )
        arguments[model_index], arguments[input_index] = str(prepared_model), str(prepared_input)
        return subprocess.call([str(executable_path()), *arguments])
