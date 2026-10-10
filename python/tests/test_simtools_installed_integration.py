"""Exercise both simtools workflows against the packaged native executable."""

import json
import logging
import math
import runpy
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

try:
    import astropy.units as u
    import simtools  # noqa: F401
except ImportError as error:
    raise unittest.SkipTest("simtools integration requires simtools and astropy") from error
from obdeect import executable_path
from obdeect.imaging_list import load_imaging_metadata
from obdeect.result_contract import read_arrivals
from simtools import settings
from simtools.ray_tracing.incident_angles import IncidentAnglesCalculator
from simtools.ray_tracing.ray_tracing import RayTracing


@pytest.fixture
def installed_model(tmp_path, monkeypatch):
    """Create a generic reflector without a database or production model assets."""
    smoke = runpy.run_path(str(Path(__file__).parents[2] / "tools" / "installed_smoke.py"))
    model = tmp_path / "model.json"
    model.write_text(json.dumps(smoke["tiny_model"]()))
    monkeypatch.setattr(settings.config, "_ray_tracing_backend", "obdeect")
    monkeypatch.setattr(settings.config, "_args", {"obdeect_optical_model_file": model})
    monkeypatch.setattr(settings.config, "_obdeect_exe", None)
    assert settings.config.obdeect_exe == executable_path("obdeect-simtools-raytrace")
    return model


def test_ray_tracing_executes_packaged_model_for_unique_offsets(installed_model, tmp_path):
    """Run the actual backend factory, file naming and native transport."""
    ray = object.__new__(RayTracing)
    ray._logger = logging.getLogger(__name__)
    ray.telescope_model = SimpleNamespace(site="North", name="LSTN-01", label="smoke")
    ray.label = "smoke"
    ray.output_directory = tmp_path
    ray.single_mirror_mode = False
    ray.zenith_angle = 0.0
    ray.mirrors = [{"source_distance": 10.0}]
    files = []
    for offset in (0.0, 0.1):
        simulator = ray._create_simulator(offset, 0, 0, ray.mirrors[0], True, False)
        simulator.run()
        photons = np.loadtxt(simulator.output_file)
        assert 0 < len(photons) < simulator.photons_per_run
        metadata = load_imaging_metadata(installed_model, 10000, offset, 0, 0)
        image = ray._create_psf_image(simulator.output_file, 1000, 0.8)
        expected_area_m2 = (
            math.pi * metadata.sampling_radius_m**2 * len(photons) / simulator.photons_per_run
        )
        assert image.get_effective_area() == pytest.approx(expected_area_m2)
        result = ray._analyze_image(image, offset, 0, offset, 0.8, 1)
        assert result[5].to_value(u.m**2) == pytest.approx(expected_area_m2)
        assert np.all(photons[:, 17] == 1)
        files.append(simulator.output_file)
        # Reuse preserves an existing result; force regenerates a damaged file.
        original = simulator.output_file.read_text()
        simulator.output_file.write_text("existing result")
        simulator.run()
        assert simulator.output_file.read_text() == "existing result"
        simulator.force_simulate = True
        simulator.run()
        assert simulator.output_file.read_text() == original
    assert files[0] != files[1]


def test_incident_angles_executes_packaged_model(installed_model, tmp_path, monkeypatch):
    """Read native incidence measurements and write the real ECSV distribution."""
    calculator = object.__new__(IncidentAnglesCalculator)
    calculator.logger = logging.getLogger(__name__)
    calculator.label = "smoke"
    calculator.results_dir = tmp_path
    calculator.calculate_primary_secondary_angles = False
    calculator.config_data = {
        "telescope": "LSTN-01",
        "source_distance": 10 * u.km,
        "number_of_photons": 100,
        "off_axis_angle": 0 * u.deg,
        "obdeect_optical_model_file": installed_model,
    }
    monkeypatch.setattr(
        "simtools.ray_tracing.incident_angles.MetadataCollector.dump", lambda **_: None
    )
    result = calculator.run()
    assert len(result) > 0
    assert result["angle_incidence_focal"].unit == u.deg
    assert (tmp_path / "incident_angles_smoke_LSTN-01_off0.ecsv").is_file()
    assert all(
        arrival.incidence_focal_deg is not None
        for arrival in read_arrivals(tmp_path / "arrivals_smoke_LSTN-01_off0.csv")
        if arrival.detected
    )
