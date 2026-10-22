"""Readers of simulation-models records and sim_telarray files."""

import math

import pytest
from obdeect.formats import (
    ModelError,
    components,
    length_m,
    lengths_m,
    read_mirror_list,
    read_response,
    read_segmentation,
    read_spatial_map,
    read_structure_table,
    secondary_segment_frame,
    segmentation_from_record,
)
from obdeect.simtel_tables import evaluate, parse_rpol_table
from telescopes import MIRROR_LIST, STRUCTURE


def test_lengths_convert_to_metres():
    assert length_m({"value": 250.0, "unit": "cm"}, "x") == 2.5
    with pytest.raises(ModelError, match="unsupported length unit"):
        length_m({"value": 1.0, "unit": "inch"}, "x")
    with pytest.raises(ModelError, match="must be positive"):
        length_m({"value": 0.0, "unit": "m"}, "x")


def test_list_parameters_take_one_or_per_value_units():
    # simtools writes a single unit when all components share it.
    assert lengths_m({"value": [1.0, 20.0], "unit": "cm"}, "x", 2) == [0.01, 0.2]
    assert lengths_m({"value": [1.0, 2.0], "unit": ["m", "mm"]}, "x", 2) == [1.0, 0.002]
    assert components({"value": [0.1, 0.0], "unit": None}, "x") == ([0.1, 0.0], [None, None])
    with pytest.raises(ModelError, match="units do not match"):
        components({"value": [1.0, 2.0], "unit": ["cm"]}, "x")
    with pytest.raises(ModelError, match="requires 3 components"):
        components({"value": [1.0, 2.0], "unit": "cm"}, "x", 3)


def test_mirror_list():
    facets = read_mirror_list(MIRROR_LIST + "10 20 30 -1500 1 -40 # trailing comment\n")
    assert len(facets) == 8
    assert facets[1] == {
        "id": 1,
        "centre_m": [1.25, 0.0, 0.0],
        "shape": "circle",
        "diameter_m": 1.0,
        "catalogue_focal_length_m": 0.0,
    }
    assert facets[-1]["shape"] == "hexagon_flat_x"
    assert facets[-1]["centre_m"] == [0.1, 0.2, -0.4]
    assert facets[-1]["catalogue_focal_length_m"] == -15.0
    with pytest.raises(ModelError, match="invalid shape"):
        read_mirror_list("0 0 120 0 7\n")


def test_segmentation():
    segments = read_segmentation("hex 1 10 20 80\nyhex 1 0 0 80 5\nring 3 100 200 120 0 2\n")
    assert [s["shape"] for s in segments] == ["hexagon", "hexagon"] + ["annular_sector"] * 3
    assert segments[1]["rotation_deg"] == 95.0
    assert [s["start_deg"] for s in segments[2:]] == [0.0, 120.0, 240.0]
    assert segments[2]["gap_m"] == 0.02
    polygon = read_segmentation("polygon 1 90 0 0 100 0 0 100\n")[0]
    assert polygon["vertices_xy_m"][1] == pytest.approx([0.0, -1.0])
    with pytest.raises(ModelError, match="convex"):
        read_segmentation("polygon 1 0 0 0 100 0 50 10 100 100 0 100\n")
    structured = segmentation_from_record({
        "value": [
            {
                "kind": "hex",
                "count": 1,
                "x": {"value": 10.0, "unit": "cm"},
                "y": {"value": 20.0, "unit": "cm"},
                "diameter": {"value": 80.0, "unit": "cm"},
            }
        ]
    })
    assert structured == read_segmentation("hex 1 10 20 80")
    flipped = secondary_segment_frame(segments)
    assert flipped[0]["centre_xy_m"] == [-0.1, 0.2]
    assert flipped[0]["rotation_deg"] == 180.0
    assert flipped[2]["start_deg"] == 60.0 and flipped[2]["gap_at_start"]


def test_response_tables():
    table = read_response("wavelength reflectivity\n300 0.5\n400 0.7\n400 0.9\n500 0.6\n", "r")
    assert table["x"] == [300.0, 400.0, 500.0]
    assert table["values"] == [0.5, 0.8, 0.6]  # Repeated wavelengths are averaged.
    rpol = read_response("#@RPOL@ 1\n300 0.5\n400 0.7\n500 0.6\n", "r")
    assert evaluate(rpol["x"], rpol["values"], rpol["interpolation"], 350.0) == pytest.approx(0.6)
    grid = read_response(
        "wavelength incidence_angle reflectivity\n300 0 0.8\n300 60 0.4\n400 0 0.9\n400 60 0.5\n",
        "r",
    )
    assert grid["y"] == [0.0, 60.0] and grid["values"] == [0.8, 0.4, 0.9, 0.5]
    with pytest.raises(ModelError, match="values in"):
        read_response("wavelength reflectivity\n300 0.5\n400 1.5\n", "r")


def test_rpol_spline_and_zero_boundary():
    table = parse_rpol_table("#@RPOL@ 1 scheme=3,clip=1\n0 0\n1 1\n2 4\n3 9\n4 16\n")
    assert evaluate(table["x"], table["response"], table["interpolation"], 2.5) == pytest.approx(
        6.25, rel=0.02
    )
    assert evaluate(table["x"], table["response"], table["interpolation"], 5.0) == 0.0


def test_spatial_map():
    table = read_spatial_map("#@RPOL@ 2\n 0 100\n0 1 0.9\n50 0.8 0.7\n", "map")
    assert table["x_m"] == [0.0, 0.5] and table["y_m"] == [0.0, 1.0]
    assert table["values"] == [[1.0, 0.9], [0.8, 0.7]]


def test_structure_tables():
    table = read_structure_table(STRUCTURE, "s")
    assert table["frame"] == "primary_mirror_frame"
    assert table["columns"][:5] == ["id", "group", "shape", "diameter", "x1"]
    mast, plate = table["rows"]
    assert mast["id"] == "mast-0" and mast["shape"] == "cylinder" and mast["diameter"] == 0.1
    assert math.isnan(mast["x3"]) and math.isnan(plate["diameter"])
    assert math.isclose(plate["z4"], 3.2)
    with pytest.raises(ModelError, match="metres"):
        read_structure_table(STRUCTURE.replace("diameter, unit: m", "diameter, unit: cm"), "s")
    with pytest.raises(ModelError, match="coordinate_frame"):
        read_structure_table(
            STRUCTURE.replace("#   coordinate_frame: primary_mirror_frame\n", ""), "s"
        )
    with pytest.raises(ModelError, match="infinite"):
        read_structure_table(STRUCTURE.replace("0.1 0.0 1.6", "inf 0.0 1.6"), "s")
