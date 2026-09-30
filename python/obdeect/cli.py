"""Launch the C++ executables shipped in the obdeect wheel."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def executable_path(name: str = "obdeect-simtools-raytrace") -> Path:
    """Return the packaged native executable named ``name``.

    Parameters
    ----------
    name : str
        Executable name, without a Windows ``.exe`` suffix.

    Raises
    ------
    FileNotFoundError
        If the wheel does not contain the requested executable.
    """

    # Editable scikit-build installations import Python directly from the
    # checkout, while CMake installs executables under site-packages.  Check
    # both locations so `pip install -e .` has the same launcher behaviour as
    # a regular wheel.
    native_dirs = [Path(__file__).resolve().parent / "_native"]
    native_dirs.extend(Path(entry) / "obdeect" / "_native" for entry in sys.path if entry)
    candidates = [directory / name for directory in native_dirs]
    if os.name == "nt" and not name.endswith(".exe"):
        candidates.extend(directory / f"{name}.exe" for directory in native_dirs)
    for candidate in candidates:
        if candidate.is_file():
            if os.name != "nt" and not os.access(candidate, os.X_OK):
                raise PermissionError(f"Packaged obdeect executable is not executable: {candidate}")
            return candidate
    raise FileNotFoundError(
        f"Packaged obdeect executable {name!r} was not found in any obdeect/_native directory. "
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

    return _run("obdeect-simtools-raytrace")
