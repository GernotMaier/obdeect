"""Python helpers and native executable discovery for obdeect."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("obdeect-dev")
except PackageNotFoundError:
    __version__ = "0.1.0"

from .cli import executable_path


def read_paths(*args, **kwargs):
    """Read traced paths, importing plotting dependencies only when needed."""

    from .plotting import read_paths as _read_paths

    return _read_paths(*args, **kwargs)


__all__ = ["executable_path", "read_paths"]
