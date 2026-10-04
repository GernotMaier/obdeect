"""Check installed entry points and trace a small, generic optical model."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sysconfig
from importlib.metadata import distribution
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from obdeect.cli import executable_path
from obdeect.result_contract import read_arrivals


def tiny_model() -> dict:
    """Return a nominal spherical reflector and a finite detector in SI units."""
    surface = {
        "id": 0,
        "shape": "circle",
        "centre_m": [0.0, 0.0, 0.0],
        "normal": [0.0, 0.0, 1.0],
        "tangent": [1.0, 0.0, 0.0],
        "diameter_m": 1.0,
        "focal_length_m": 10.0,
    }
    detector = {**surface, "id": 1, "centre_m": [0.0, 0.0, 10.0]}
    detector.pop("focal_length_m")
    model = {
        "format": "obdeect.compiled-optical-model.v1",
        "provenance": {"model": "installation-check", "model_version": "1.0.0"},
        "report": {
            "native_trace_ready": True,
            "production_trace_ready": False,
            "trace_blockers": ["generic installation check; no production qualification"],
        },
        "trace_model": {
            "kind": "segmented",
            "primary_facets": [surface],
            "detector_surfaces": [detector],
            "cylinder_obscurers": [],
            "primary_reflectivity": [
                {"wavelength_nm": 300.0, "response": 0.8},
                {"wavelength_nm": 500.0, "response": 0.8},
            ],
        },
    }
    canonical = json.dumps(model, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    model["optical_model_sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    return model


def main() -> None:
    """Run without relying on the source checkout or telescope model data."""
    from obdeect import _core

    for entry in distribution("obdeect-dev").entry_points:
        if entry.group == "console_scripts":
            command = Path(sysconfig.get_path("scripts")) / entry.name
            subprocess.run([str(command), "--help"], check=True, stdout=subprocess.DEVNULL)
    with TemporaryDirectory() as directory:
        root = Path(directory)
        model_path = root / "model.json"
        model_path.write_text(json.dumps(tiny_model(), ensure_ascii=False, sort_keys=True))
        output = root / "arrivals.csv"
        subprocess.run(
            [
                str(executable_path()),
                "--optical-model",
                str(model_path),
                "--photons",
                "8",
                "--distance-m",
                "100",
                "--output",
                str(output),
            ],
            check=True,
        )
        arrivals = read_arrivals(output)
        assert len(arrivals) == 8 and any(row.detected for row in arrivals)
        assert all(abs(row.throughput - 0.8) < 1e-12 for row in arrivals if row.detected)
        with output.open(newline="") as handle:
            assert len(list(csv.DictReader(handle))) == 8
        bulk = _core.OpticalModel(str(model_path)).trace(
            np.array([[0.0, 0.0, 50.0]]),
            np.array([[0.0, 0.0, -1.0]]),
            np.array([400.0]),
            np.array([0.0]),
            np.array([1.0]),
            np.array([0], dtype=np.uint64),
        )
        assert _core.STATUS_NAMES[bulk["status"][0]] == "detected"
        np.testing.assert_array_equal(bulk["position_m"], [[0.0, 0.0, 10.0]])
        np.testing.assert_allclose(bulk["optical_weight"], [0.8], rtol=0, atol=1e-12)
    print("Installed commands and generic compiled-model trace passed.")


if __name__ == "__main__":
    main()
