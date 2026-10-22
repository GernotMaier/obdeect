"""Optical models end to end: compile, store, trace, analyse, command line."""

import json
import math

import numpy as np
import pytest
from obdeect import ModelError, OpticalModel
from obdeect.cli import main
from obdeect.model import DETECTED, STATUS_NAMES
from obdeect.psf import containment, summary
from obdeect.repository import read_telescope
from telescopes import dual, record, segmented, write_simulation_models


def model_of(records_files, **options):
    records, files = records_files
    return OpticalModel.from_records(
        records, files.__getitem__, configuration_seed=options.pop("seed", 1), **options
    )


def test_store_and_reload(tmp_path):
    model = model_of(segmented())
    path = tmp_path / "model.json"
    model.save(path)
    assert OpticalModel.load(path).document == model.document
    tampered = json.loads(path.read_text())
    tampered["optics"]["focal_surface"]["plane_z_m"] += 0.01
    with pytest.raises(ModelError, match="sha256"):
        OpticalModel(tampered)


def test_model_document_is_isolated_from_mutation(tmp_path):
    document = model_of(segmented()).document
    model = OpticalModel(document)
    original_area = model.sampling_area_m2
    original_hash = model.sha256
    original_photons = model.star_photons(8, seed=3)

    document["sampling"]["disk_radius_m"] = 1.0
    document["optics"]["focal_surface"]["plane_z_m"] = -99.0
    exposed = model.document
    exposed["sampling"]["disk_radius_m"] = 2.0
    exposed["optics"]["focal_surface"]["plane_z_m"] = -50.0

    assert model.sampling_area_m2 == original_area
    assert model.sha256 == original_hash
    assert np.array_equal(
        model.star_photons(8, seed=3)["position_m"], original_photons["position_m"]
    )
    assert model.document["optics"]["focal_surface"]["plane_z_m"] != -50.0
    saved = tmp_path / "model.json"
    model.save(saved)
    assert json.loads(saved.read_text()) == OpticalModel.load(saved).document
    assert json.loads(saved.read_text())["sampling"]["disk_radius_m"] != 1.0


def test_native_errors_name_the_json_path():
    document = model_of(segmented()).document
    document["optics"]["primary"]["facets"][2]["normal"] = [0.0, 0.0, 2.0]
    from obdeect.compiler import model_hash

    document["sha256"] = model_hash(document)
    with pytest.raises(ModelError, match=r"facets\[2\]\.normal: expected a unit vector"):
        OpticalModel(document)


def test_on_axis_star():
    model = model_of(segmented())
    photons = model.star_photons(200000, seed=4)
    arrivals = model.trace(photons, seed=4, apply_responses=False)
    result = summary(arrivals, model.sampling_area_m2)
    # Seven 1 m circles, projected along the axis, minus the 0.6 m square
    # camera shadow (inside the central mirror).
    normals = [f["normal"][2] for f in model.document["optics"]["primary"]["facets"]]
    expected = math.pi * 0.25 * sum(normals) - 0.6**2
    assert result["effective_area_m2"] == pytest.approx(expected, rel=0.005)
    assert abs(result["centroid_x_m"]) < 1e-4 and abs(result["centroid_y_m"]) < 1e-4
    detected = arrivals["status"] == DETECTED
    assert np.allclose(arrivals["position_m"][detected, 2], 8.5)
    assert np.all(arrivals["primary_surface"][detected] < 7)
    assert np.all(np.isnan(arrivals["secondary_point_m"]))
    blocked = arrivals["status"] == STATUS_NAMES.index("blocked")
    assert np.all(
        arrivals["terminal_surface"][blocked] == model.document["optics"]["obscurers"][0]["id"]
    )
    # Optical responses: the telescope transmission is 0.9 on axis.
    weighted = model.trace(photons, seed=4)
    assert weighted["weight"][detected] == pytest.approx(0.9)


def test_results_do_not_depend_on_threads_or_order():
    model = model_of(
        segmented(
            mirror_reflection_random_angle=record(
                "mirror_reflection_random_angle", [0.01, 0.2, 0.05], ["deg", "null", "deg"]
            )
        )
    )
    photons = model.star_photons(50000, seed=9, towards_source=(0.01, 0.0, math.sqrt(1 - 1e-4)))
    one = model.trace(photons, seed=2, threads=1)
    many = model.trace(photons, seed=2, threads=8)
    for name in one:
        assert np.array_equal(one[name], many[name], equal_nan=True), name
    reversed_photons = {name: values[::-1] for name, values in photons.items()}
    backwards = model.trace(reversed_photons, seed=2, threads=3)
    assert np.array_equal(backwards["position_m"][::-1], one["position_m"])
    assert not np.array_equal(model.trace(photons, seed=3)["direction"], one["direction"])


def test_star_photons_depend_only_on_their_ids():
    model = model_of(segmented())
    many = model.star_photons(20000, seed=5, distance_m=50.0)
    part = model.star_photons(7000, seed=5, distance_m=50.0, first_id=9000)
    for name, values in part.items():
        assert np.array_equal(values, many[name][9000:16000]), name


def test_emitter_photons_cross_the_disk_in_the_expected_number():
    model = model_of(segmented())
    sampling = model.document["sampling"]
    height, radius = 500.0, sampling["disk_radius_m"]
    position = (0.0, 0.0, sampling["disk_z_m"] + height)
    emitted = 2e8
    photons = model.emitter_photons(
        emitted, seed=4, position_m=position, intensity=([0.0, 180.0], [1.0, 1.0])
    )
    expected = emitted * (1.0 - height / math.hypot(height, radius)) / 2.0
    assert abs(len(photons["weight"]) - expected) < 4.0 * math.sqrt(expected)
    assert np.all((photons["weight"] > 0.99) & (photons["weight"] <= 1.0))
    assert np.all(photons["position_m"] == position)
    again = model.emitter_photons(
        emitted, seed=4, position_m=position, intensity=([0.0, 180.0], [1.0, 1.0])
    )
    assert np.array_equal(again["direction"], photons["direction"])
    assert np.all(photons["time_ns"] == 0.0)


def test_flasher_light_goes_straight_to_the_camera():
    model = model_of(segmented())
    camera = model.focal_surface_z_m
    photons = model.emitter_photons(
        1e7,
        seed=2,
        position_m=(0.0, 0.0, 0.0),
        intensity=([0.0, 10.0, 10.001], [1.0, 1.0, 0.0]),
        pulse=([-1.0, 1.0], [1.0, 1.0]),
        disk=((0.0, 0.0, camera), (0.0, 0.0, 1.0), 0.5),
    )
    arrivals = model.trace(photons, seed=2, direct=True)
    assert np.all(arrivals["status"] == DETECTED)
    assert np.allclose(arrivals["position_m"][:, 2], camera)
    assert np.array_equal(arrivals["weight"], photons["weight"])
    assert np.all(arrivals["primary_surface"] == np.iinfo(np.uint32).max)  # No mirror.
    assert np.all(np.abs(photons["time_ns"]) <= 1.0)


def test_off_axis_image_moves_by_the_plate_scale():
    model = model_of(segmented(camera_body_diameter=record("camera_body_diameter", 0.0, "cm")))
    theta = math.radians(1.0)
    photons = model.star_photons(
        50000, seed=1, towards_source=(math.sin(theta), 0.0, math.cos(theta))
    )
    result = summary(model.trace(photons, seed=1, apply_responses=False), model.sampling_area_m2)
    assert result["centroid_x_m"] == pytest.approx(-8.0 * math.tan(theta), rel=0.01)


def test_dual_mirror_focus():
    model = model_of(dual(), zenith_deg=20.0)
    photons = model.star_photons(100000, seed=6, distance_m=1e7)
    arrivals = model.trace(photons, seed=6, apply_responses=False)
    result = summary(arrivals, model.sampling_area_m2)
    # The SST design images a star on-axis into ~6 mm (sim_telarray: 7.05 mm with mirror errors).
    assert 0.005 < result["d80_m"] < 0.007
    detected = arrivals["status"] == DETECTED
    assert np.all(np.isfinite(arrivals["secondary_point_m"][detected]))
    assert np.count_nonzero(arrivals["status"] == STATUS_NAMES.index("numerical_failure")) == 0
    assert 6.0 < result["effective_area_m2"] < 10.0


def test_containment():
    x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    result = containment(x, np.zeros(5), fraction=0.8)
    assert result["centroid_x_m"] == 2.0 and result["radius_m"] == 2.0


def test_repository_and_command_line(tmp_path, capsys):
    root = tmp_path / "simulation-models"
    records, files = segmented()
    write_simulation_models(root, "TEST-design", records, files)
    telescope = read_telescope(root, "TEST-design", "9.9.9")
    assert set(telescope.records) == set(records) and telescope.propagation_group_index == 1.0
    model, arrivals = tmp_path / "model.json", tmp_path / "arrivals.npz"
    compile_args = f"compile --simulation-models {root} --version 9.9.9 --seed 3 --output {model}"
    assert main([*compile_args.split(), "--model", "TEST-design"]) == 0
    trace_args = f"trace --model {model} --seed 1 --photons 2000 --field-deg 0.5 0"
    assert main([*trace_args.split(), "--output", str(arrivals)]) == 0
    report = json.loads(capsys.readouterr().out.split("\n", 1)[1])
    assert report["status"]["detected"] > 0
    assert report["model_sha256"] == OpticalModel.load(model).sha256
    assert np.load(arrivals)["position_m"].shape == (2000, 3)
    csv = tmp_path / "arrivals.csv"
    assert main(f"trace --model {model} --seed 1 --photons 10 --output {csv}".split()) == 0
    assert csv.read_text().startswith("photon_id,status,")
    plot = tmp_path / "plot.png"
    assert main(f"plot --model {model} --arrivals {arrivals} --output {plot}".split()) == 0
    assert plot.stat().st_size > 0
    assert main([*compile_args.split(), "--model", "MISSING"]) == 1


def test_patch_production_overlays_its_history(tmp_path):
    root = tmp_path / "simulation-models"
    records, files = segmented()
    write_simulation_models(root, "TEST-design", records, files)
    patch = root / "productions" / "9.9.10"
    patch.mkdir()
    (patch / "info.yml").write_text(
        'model_version: "9.9.10"\nmodel_update: "patch_update"\n'
        'model_version_history:\n  - "9.9.9"\nmodel_status: "production"\n'
    )
    table = {"focal_length": "2.0.0"}
    manifest = {"model_version": "9.9.10", "parameters": {"TEST-design": table}}
    manifest["deprecated_parameters"] = ["camera_depth"]
    (patch / "TEST-design.json").write_text(json.dumps(manifest))
    focal = records["focal_length"] | {"parameter_version": "2.0.0", "value": 810.0}
    directory = root / "model_parameters" / "TEST-design" / "focal_length"
    (directory / "focal_length-2.0.0.json").write_text(json.dumps(focal))

    patched = read_telescope(root, "TEST-design", "9.9.10")
    assert patched.records["focal_length"]["value"] == 810.0
    assert "camera_depth" not in patched.records
    assert set(patched.records) == set(records) - {"camera_depth"}
    assert read_telescope(root, "TEST-design", "9.9.9").records["focal_length"]["value"] == 800.0
    with pytest.raises(ModelError, match="no production table MISSING in 9.9.10"):
        read_telescope(root, "MISSING", "9.9.10")
