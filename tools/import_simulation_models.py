#!/usr/bin/env python3
"""Backward-compatible executable wrapper for the installed importer."""

from __future__ import annotations

import sys
from pathlib import Path

# Keep ``python tools/import_simulation_models.py`` usable from a source tree
# while the supported installed entry point remains ``obdeect-import-simulation-models``.
sys.path.insert(0, str(Path(__file__).parents[1] / "python"))

# Re-export the public importer API for callers that historically loaded this
# source-tree script as a module. The installed console entry point uses the
# package module directly.
from obdeect.model_import import ImportError, main, resolve_model

__all__ = ["ImportError", "main", "resolve_model"]

if __name__ == "__main__":
    main()
