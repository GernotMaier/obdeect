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

    native_dir = Path(__file__).resolve().parent / "_native"
    candidates = [native_dir / name]
    if os.name == "nt" and not name.endswith(".exe"):
        candidates.append(native_dir / f"{name}.exe")
    for candidate in candidates:
        if candidate.is_file():
            if os.name != "nt" and not os.access(candidate, os.X_OK):
                raise PermissionError(f"Packaged obdeect executable is not executable: {candidate}")
            return candidate
    raise FileNotFoundError(
        f"Packaged obdeect executable {name!r} was not found in {native_dir}. "
        "Reinstall obdeect-dev with a wheel for this platform."
    )


def _run(name: str) -> int:
    """Forward command-line arguments to a packaged C++ executable."""

    command = [str(executable_path(name)), *sys.argv[1:]]
    if os.name == "nt":
        return subprocess.call(command)
    os.execv(command[0], command)
    return 1


def reference_main() -> int:
    """Run the packaged  reference ray tracer."""

    return _run("obdeect_reference")


def ctao_main() -> int:
    """Run the packaged CTAO reference ray tracer."""

    return _run("obdeect_ctao")


def simtools_raytrace_main() -> int:
    """Run the packaged simtools-compatible ray tracer."""

    return _run("obdeect-simtools-raytrace")
