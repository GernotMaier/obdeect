"""Installed bulk tracing owns results and borrows validated contiguous input."""

import csv
import gc
import hashlib
import importlib
import json
import runpy
import subprocess
import sys
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pytest
from obdeect.cli import executable_path
from obdeect.result_contract import read_arrivals

try:
    from obdeect import _core as core
except ImportError:
    build_module = Path(__file__).resolve().parents[2] / "build/python-api"
    sys.path.insert(0, str(build_module))
    try:
        core = importlib.import_module("_core")
    except ImportError as error:
        raise unittest.SkipTest(
            "bulk tracing requires the optional native Python module"
        ) from error


@pytest.fixture
def model_path(tmp_path):
    helper = Path(__file__).resolve().parents[2] / "tools/installed_smoke.py"
    model = runpy.run_path(str(helper))["tiny_model"]()
    path = tmp_path / "optical-model.json"
    path.write_text(json.dumps(model, ensure_ascii=False, sort_keys=True))
    return path


def photons(count=8):
    positions = np.zeros((count, 3), dtype=np.float64)
    positions[:, 0] = np.linspace(-0.2, 0.2, count)
    positions[:, 2] = 50
    directions = np.zeros_like(positions)
    directions[:, 2] = -1
    return (
        positions,
        directions,
        np.full(count, 400.0),
        np.arange(count, dtype=np.float64),
        np.linspace(0.5, 2.0, count),
        np.arange(count, dtype=np.uint64),
    )


def test_bulk_loss_closure_and_output_lifetime(model_path):
    tracer = core.OpticalModel(str(model_path))
    arrays = photons()
    arrays[0][-1, 0] = 2  # misses the finite primary
    before = tuple(array.copy() for array in arrays)
    for array in arrays:
        array.flags.writeable = False
    result = tracer.trace(*arrays)
    np.testing.assert_array_equal(result["status"][:-1], 0)
    assert core.STATUS_NAMES[result["status"][-1]] == "missed_primary"
    np.testing.assert_allclose(result["optical_weight"][:-1], arrays[4][:-1] * 0.8)
    np.testing.assert_allclose(
        result["optical_weight"] + result["response_loss_weight"] + result["terminal_loss_weight"],
        arrays[4],
    )
    for original, saved in zip(arrays, before, strict=True):
        np.testing.assert_array_equal(original, saved)
    retained = result["position_m"]
    expected = retained.copy()
    del result, tracer, arrays
    gc.collect()
    np.testing.assert_array_equal(retained, expected)


def test_bulk_order_blocks_and_concurrent_calls(model_path):
    tracer = core.OpticalModel(str(model_path))
    arrays = photons(19)
    expected = tracer.trace(*arrays)
    reversed_result = tracer.trace(*(np.ascontiguousarray(array[::-1]) for array in arrays))
    for name, values in expected.items():
        np.testing.assert_array_equal(values, reversed_result[name][::-1])
    chunks = [
        tracer.trace(*(array[start : start + 4] for array in arrays)) for start in range(0, 19, 4)
    ]
    for name, values in expected.items():
        np.testing.assert_array_equal(values, np.concatenate([chunk[name] for chunk in chunks]))
    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(lambda _: tracer.trace(*arrays), range(2)))
    for result in results:
        for name, values in expected.items():
            np.testing.assert_array_equal(values, result[name])


def test_bulk_rejects_invalid_input_and_tampered_model(model_path):
    tracer = core.OpticalModel(str(model_path))
    arrays = list(photons())
    with pytest.raises(TypeError):
        tracer.trace(arrays[0].astype(np.float32), *arrays[1:])
    with pytest.raises(TypeError):
        tracer.trace(arrays[0][:, ::-1], *arrays[1:])
    with pytest.raises(ValueError, match="same length"):
        tracer.trace(arrays[0][:-1], *arrays[1:])
    arrays[5][1] = arrays[5][0]
    with pytest.raises(ValueError, match="unique"):
        tracer.trace(*arrays)
    arrays = list(photons())
    arrays[1][0, 2] = -2
    with pytest.raises(ValueError, match="unit directions"):
        tracer.trace(*arrays)
    data = json.loads(model_path.read_text())
    data["trace_model"]["detector_surfaces"][0]["diameter_m"] *= 2
    model_path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="content hash"):
        core.OpticalModel(str(model_path))


def test_empty_bulk_batch(model_path):
    result = core.OpticalModel(str(model_path)).trace(*photons(0))
    assert result["position_m"].shape == (0, 3)
    assert all(len(value) == 0 for value in result.values())


def test_bulk_material_window_and_strict_fields(tmp_path):
    """Synthetic analytic window: phase/group timing and native loader parity."""

    def surface(identifier, z, role, **extra):
        return {
            "id": identifier,
            "shape": "circle",
            "diameter_m": 4,
            "role": role,
            "frame": {
                "origin_m": [0, 0, z],
                "x_axis": [1, 0, 0],
                "y_axis": [0, 1, 0],
                "z_axis": [0, 0, 1],
            },
            **extra,
        }

    def curve(value):
        return [{"wavelength_nm": wavelength, "value": value} for wavelength in (300, 500)]

    model = {
        "format": "obdeect.compiled-optical-model.v1",
        "provenance": {"model": "synthetic-window", "model_version": "1"},
        "trace_model": {
            "kind": "nonsequential",
            "max_interactions": 8,
            "entrance_medium_id": 4294967295,
            "materials": [
                {
                    "id": 1,
                    "phase_index": curve(1.5),
                    "group_index": curve(1.8),
                    "absorption_per_m": curve(0.2),
                }
            ],
            "surfaces": [
                surface(1, 2, "refractive_interface", front_medium_id=4294967295, back_medium_id=1),
                surface(2, 1, "refractive_interface", front_medium_id=1, back_medium_id=4294967295),
                surface(3, 0, "detector"),
            ],
        },
    }
    path = tmp_path / "window.json"

    def write():
        data = {key: value for key, value in model.items() if key != "optical_model_sha256"}
        model["optical_model_sha256"] = hashlib.sha256(
            json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        path.write_text(json.dumps(model, sort_keys=True))

    write()
    arrays = list(photons(1))
    arrays[0][0] = [0, 0, 3]
    arrays[3][0] = 5
    result = core.OpticalModel(str(path)).trace(*arrays)
    assert core.STATUS_NAMES[result["status"][0]] == "detected"
    assert result["path_length_m"][0] == 3
    assert result["optical_path_m"][0] == 3.5
    assert result["arrival_time_ns"][0] == pytest.approx(5 + 3.8 / 0.299792458)
    assert result["throughput"][0] == pytest.approx(0.96**2 * np.exp(-0.2))
    model["trace_model"]["surfaces"][0]["unimplemented_coating"] = 0.8
    write()
    with pytest.raises(ValueError, match="invalid schema"):
        core.OpticalModel(str(path))


def test_bulk_measured_camera_response_and_keyed_scatter(model_path):
    data = json.loads(model_path.read_text())
    data.pop("optical_model_sha256")
    trace = data["trace_model"]
    trace["camera_response"] = {
        "semantics": "measured_complete_response",
        "camera_transmission": 0.5,
        "camera_filter": [{"wavelength_nm": w, "response": 0.8} for w in (300, 500)],
        "lightguide_efficiency": [{"incidence_angle_deg": a, "response": 0.9} for a in (0, 90)],
    }
    trace["primary_scatter"] = {
        "sigma1_rad": 0.001,
        "sigma2_rad": 0.002,
        "fraction2": 0.2,
        "method": "outgoing_angles",
        "seed": 123,
    }
    data["optical_model_sha256"] = hashlib.sha256(
        json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    model_path.write_text(json.dumps(data, sort_keys=True))
    tracer = core.OpticalModel(str(model_path))
    arrays = photons()
    result = tracer.trace(*arrays)
    np.testing.assert_allclose(result["throughput"], 0.8 * 0.5 * 0.8 * 0.9)
    reversed_result = tracer.trace(*(np.ascontiguousarray(array[::-1]) for array in arrays))
    for key, values in result.items():
        np.testing.assert_array_equal(values, reversed_result[key][::-1])
    pieces = [
        tracer.trace(*(array[start : start + 3] for array in arrays))
        for start in range(0, len(arrays[0]), 3)
    ]
    for key, values in result.items():
        np.testing.assert_array_equal(values, np.concatenate([piece[key] for piece in pieces]))


def test_bulk_matches_installed_native_photon_replay(model_path, tmp_path):
    arrays = photons()
    source = tmp_path / "photons.csv"
    fields = (
        "run_id,event_id,array_id,telescope_id,photon_id,x_m,y_m,z_m,dx,dy,dz,"
        "wavelength_nm,time_ns,weight"
    )
    with source.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields.split(","))
        for i in range(len(arrays[0])):
            writer.writerow((
                0,
                0,
                0,
                0,
                i,
                *arrays[0][i],
                *arrays[1][i],
                400,
                arrays[3][i],
                arrays[4][i],
            ))
    output = tmp_path / "arrivals.csv"
    try:
        executable = executable_path()
    except FileNotFoundError:
        executable = (
            Path(__file__).resolve().parents[2] / "build/python-api/obdeect-simtools-raytrace"
        )
    assert executable.is_file()
    subprocess.run(
        [
            str(executable),
            "--optical-model",
            str(model_path),
            "--photon-input",
            str(source),
            "--input-block-size",
            "3",
            "--output",
            str(output),
        ],
        check=True,
    )
    native = read_arrivals(output)
    bulk = core.OpticalModel(str(model_path)).trace(*arrays)
    for i, arrival in enumerate(native):
        assert core.STATUS_NAMES[bulk["status"][i]] == arrival.status
        assert bulk["photon_id"][i] == arrival.photon_id
        assert bulk["terminal_surface_id"][i] == arrival.terminal_surface_id
        assert bulk["arrival_time_ns"][i] == arrival.arrival_time_ns
        assert bulk["optical_weight"][i] == arrival.optical_weight
        np.testing.assert_array_equal(bulk["position_m"][i], arrival.interaction_points_m[-1])
