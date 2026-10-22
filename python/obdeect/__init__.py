"""obdeect: optical ray tracing for imaging atmospheric Cherenkov telescopes."""

import importlib
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("obdeect-dev")
except PackageNotFoundError:
    __version__ = "0.0.0"

# Imported on first use, so that reading the version needs no NumPy.
_EXPORTS = {
    "DETECTED": "obdeect.model",
    "STATUS_NAMES": "obdeect.model",
    "OpticalModel": "obdeect.model",
    "ModelError": "obdeect.formats",
}
__all__ = list(_EXPORTS)


def __getattr__(name: str):
    if name in _EXPORTS:
        return getattr(importlib.import_module(_EXPORTS[name]), name)
    raise AttributeError(f"module 'obdeect' has no attribute {name!r}")
