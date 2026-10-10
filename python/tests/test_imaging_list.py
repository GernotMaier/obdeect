import csv
import json
import math
from pathlib import Path

import pytest
from obdeect.imaging_list import load_imaging_metadata, write_imaging_list
from obdeect.result_contract import ArrivalContractError


def model_file(tmp_path: Path, kind: str = "segmented") -> Path:
    trace = {
        "kind": kind,
        "primary_facets": [
            {"centre_m": [3, 4, 0], "diameter_m": 2, "shape": "circle", "normal": [0, 0, 1]}
        ],
    }
    if kind == "axisymmetric":
        trace = {
            "kind": kind,
            "primary": {"outer_radius_m": 6, "vertex_z_m": 0},
            "secondary": {"vertex_z_m": 2},
            "detector": {"vertex_z_m": 1},
        }
    path = tmp_path / "model.json"
    path.write_text(
        json.dumps({
            "trace_model": trace,
            "random_seeds": {"detector_configuration_seed": 12},
            "focal_length_m": 10,
            "camera": {"rotation_deg": 30},
            "optical_model_sha256": "a" * 64,
        })
    )
    return path


def arrivals_file(tmp_path: Path) -> Path:
    path = tmp_path / "arrivals.csv"
    fields = (
        "contract_version,photon_id,source_kind,wavelength_nm,emission_time_ns,source_weight,"
        "throughput,status,point_count,path_length_m,x0_m,y0_m,z0_m,x1_m,y1_m,z1_m,"
        "x2_m,y2_m,z2_m,final_dx,final_dy,final_dz,incidence_primary_deg,"
        "incidence_focal_deg,arrival_time_ns"
    ).split(",")
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        for index, response in enumerate([0, 0.5, 0]):
            writer.writerow(
                dict(
                    zip(
                        fields,
                        [
                            "obdeect-arrival-v1",
                            index,
                            "star",
                            400,
                            0,
                            1,
                            response,
                            "detected" if index < 2 else "missed_primary",
                            3 if index < 2 else 1,
                            10,
                            0,
                            0,
                            6,
                            0,
                            0,
                            0,
                            1 + index,
                            2,
                            10,
                            0,
                            0,
                            1,
                            20,
                            0,
                            7,
                        ],
                        strict=True,
                    )
                )
            )
    return path


@pytest.mark.parametrize("kind", ["segmented", "axisymmetric"])
def test_model_metadata_matches_simtel_launch_disk(tmp_path: Path, kind: str) -> None:
    metadata = load_imaging_metadata(model_file(tmp_path, kind), 10000, 0, 2.5)
    assert metadata.sampling_radius_m == pytest.approx(7.2)
    assert metadata.prime_focus == (kind == "segmented")
    assert metadata.camera_rotation_deg == 30


@pytest.mark.parametrize(
    ("shape", "factor"),
    [
        ("circle", 1),
        ("square", math.sqrt(2)),
        ("hexagon", 2 / math.sqrt(3)),
        ("hexagon_flat_x", 2 / math.sqrt(3)),
        ("hexagon_flat_y", 2 / math.sqrt(3)),
    ],
)
def test_sampling_radius_contains_non_circular_facet_corners(tmp_path, shape, factor) -> None:
    path = model_file(tmp_path)
    model = json.loads(path.read_text())
    model["trace_model"]["primary_facets"][0]["shape"] = shape
    model["trace_model"]["primary_facets"][0]["normal"] = [0, 0, 1]
    path.write_text(json.dumps(model))
    metadata = load_imaging_metadata(path, 10000, 0, 0)
    assert metadata.sampling_radius_m == pytest.approx(1.2 * (5 + factor))


def test_alignment_zenith_is_independent_of_source_offset(tmp_path: Path) -> None:
    path = model_file(tmp_path)
    model = json.loads(path.read_text())
    model["primary"] = {"alignment": {"zenith_angle_deg": 20}}
    path.write_text(json.dumps(model))
    assert load_imaging_metadata(path, 10000, 0.5, 0, 20).sampling_radius_m == pytest.approx(7.2)
    with pytest.raises(ValueError, match="different telescope zenith"):
        load_imaging_metadata(path, 10000, 0.5, 0, 20.5)


def test_imaging_list_preserves_geometric_rows_and_header(tmp_path: Path) -> None:
    metadata = load_imaging_metadata(model_file(tmp_path), 10000, 0, 2.5)
    output = tmp_path / "image.lis"
    assert write_imaging_list(arrivals_file(tmp_path), output, metadata, 3) == 2
    lines = output.read_text().splitlines()
    header = next(line for line in lines if "falling on an area" in line).split()
    assert int(header[4]) == 3
    assert float(header[14]) == pytest.approx(math.pi * 7.2**2)
    rows = [list(map(float, line.split())) for line in lines if not line.startswith("#")]
    assert len(rows[0]) == 33
    assert rows[0][17] == 0  # A zero-throughput geometric hit is still an imaging row.
    theta = math.radians(30)
    assert rows[0][2] * math.cos(theta) - rows[0][3] * math.sin(theta) == pytest.approx(100)
    assert rows[0][2] * math.sin(theta) + rows[0][3] * math.cos(theta) == pytest.approx(200)
    assert "# Focal_length = 1000 cm" in lines
    assert "# Camera rotation angle = 30 deg" in lines


def test_imaging_list_rejects_incorrect_normalization(tmp_path: Path) -> None:
    metadata = load_imaging_metadata(model_file(tmp_path), 10000, 0, 0)
    with pytest.raises(ArrivalContractError, match="count"):
        write_imaging_list(arrivals_file(tmp_path), tmp_path / "image.lis", metadata, 4)
    with pytest.raises(ValueError, match="overwrite"):
        path = arrivals_file(tmp_path)
        write_imaging_list(path, path, metadata, 3)


@pytest.mark.parametrize(
    "field,value,message",
    [
        ("source_weight", "2", "unit-weight"),
        ("source_kind", "laser", "unit-weight"),
        ("sampling_area_m2", "1", "sampling area"),
    ],
)
def test_imaging_list_rejects_incompatible_source(
    tmp_path: Path, field: str, value: str, message: str
) -> None:
    metadata = load_imaging_metadata(model_file(tmp_path), 10000, 0, 0)
    path = arrivals_file(tmp_path)
    with path.open() as source:
        reader = csv.DictReader(source)
        fields = list(reader.fieldnames or ())
        rows = list(reader)
    if field not in fields:
        fields.append(field)
    for row in rows:
        row[field] = value
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ArrivalContractError, match=message):
        write_imaging_list(path, tmp_path / "image.lis", metadata, 3)


def test_launch_plane_is_above_obscurers(tmp_path: Path) -> None:
    path = model_file(tmp_path)
    model = json.loads(path.read_text())
    model["trace_model"]["cylinder_obscurers"] = [
        {
            "first_endpoint_m": [0, 0, 20],
            "second_endpoint_m": [0, 0, 21],
            "diameter_m": 2,
        }
    ]
    path.write_text(json.dumps(model))
    metadata = load_imaging_metadata(path, 10000, 0, 0)
    assert metadata.entrance_z_m > 22
    with pytest.raises(ValueError, match="source distance"):
        load_imaging_metadata(path, 10, 0, 0)


def test_missing_mirror_effects_rejected_and_seeds_recorded(tmp_path: Path) -> None:
    path = model_file(tmp_path)
    model = json.loads(path.read_text())
    model["report"] = {
        "trace_blockers": ["run-specific panel alignment and distance are not compiled"]
    }
    path.write_text(json.dumps(model))
    with pytest.raises(ValueError, match="recompile"):
        load_imaging_metadata(path, 10000, 0, 0)
    model["report"]["trace_blockers"] = []
    path.write_text(json.dumps(model))
    metadata = load_imaging_metadata(path, 10000, 0, 0)
    output = tmp_path / "image.lis"
    write_imaging_list(arrivals_file(tmp_path), output, metadata, 3, ray_tracing_seed=13)
    assert "# detector_configuration_seed = 12" in output.read_text()
    assert "# ray_tracing_seed = 13" in output.read_text()


def test_legacy_scatter_blocker_and_alignment_seed_are_preserved(tmp_path: Path) -> None:
    path = model_file(tmp_path)
    model = json.loads(path.read_text())
    model["report"] = {"trace_blockers": ["mirror scatter requires an explicit --scatter-seed"]}
    model["random_seeds"] = {}
    model["primary"] = {"alignment": {"seed": 41}}
    path.write_text(json.dumps(model))
    with pytest.raises(ValueError, match="recompile"):
        load_imaging_metadata(path, 10000, 0, 0)
    model["report"]["trace_blockers"] = []
    path.write_text(json.dumps(model))
    metadata = load_imaging_metadata(path, 10000, 0, 0)
    assert metadata.detector_configuration_seed == 41
