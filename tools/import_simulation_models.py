#!/usr/bin/env python3
"""Backward-compatible executable wrapper for the installed importer."""

from __future__ import annotations

import sys
from pathlib import Path

# Keep ``python tools/import_simulation_models.py`` usable from a source tree
# while the supported installed entry point remains ``obdeect-import-simulation-models``.
sys.path.insert(0, str(Path(__file__).parents[1] / "python"))

from obdeect.model_import import main

if __name__ == "__main__":
    main()
