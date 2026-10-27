"""Compilation of parameter records into optical models."""

import math

import numpy as np
import pytest
from obdeect.compiler import compile_model, convex_hull, model_hash, response_grid
from obdeect.formats import ModelError
from telescopes import STRUCTURE, dual, record, segmented


def build(records_files, **options):
    records, files = records_files
    options.setdefault("configuration_seed", 1)
    return compile_model(records, files.__getitem__, **options)


def test_davies_cotton_chief_rays_meet_at_the_focus():
    optics = build(segmented())["optics"]
    focus = np.array([0.0, 0.0, 8.0 + 0.5])  # f - mirror_offset
    assert optics["focal_surface"]["plane_z_m"] == pytest.approx(focus[2])
    for facet in optics["primary"]["facets"]:
        centre, normal = np.array(facet["centre_m"]), np.array(facet["normal"])
        reflected = np.array([0, 0, -1.0]) - 2 * np.dot([0, 0, -1.0], normal) * normal
        to_focus = (focus - centre) / np.linalg.norm(focus - centre)
        assert np.allclose(reflected, to_focus, atol=1e-12)
        # A Davies-Cotton mirror's focal length is its distance to the focus.
        assert facet["radius_of_curvature_m"] == pytest.approx(2 * np.linalg.norm(focus - centre))
        assert abs(np.dot(facet["tangent"], normal)) < 1e-12


def test_random_errors_follow_the_seed():
    errors = {
        "mirror_align_random_horizontal": record(
            "mirror_align_random_horizontal", [0.01, 20.0, 0.0, 0.0], ["deg", "deg", "null", "null"]
        ),
        "random_focal_length": record("random_focal_length", [5.0, 0.0], ["cm", "cm"]),
    }
    first = build(segmented(**errors), configuration_seed=7)
    assert first == build(segmented(**errors), configuration_seed=7)
    other = build(segmented(**errors), configuration_seed=8)
    assert first["sha256"] != other["sha256"]
    normals = [f["normal"] for f in first["optics"]["primary"]["facets"]]
    assert normals != [f["normal"] for f in other["optics"]["primary"]["facets"]]
    assert first["provenance"]["configuration_seed"] == 7
    assert first["sha256"] == model_hash(first)


def test_focus_offset_zenith_terms_need_both_components():
    def plane(values, zenith=60.0):
        focus = record("focus_offset", values, ["cm", "deg", "cm", "cm"])
        return build(segmented(focus_offset=focus), zenith_deg=zenith)["optics"]["focal_surface"][
            "plane_z_m"
        ]

    assert plane([2.0, 20.0, 0.0, 0.0]) == pytest.approx(8.52)
    assert plane([2.0, 20.0, 3.0, 0.0]) == pytest.approx(8.52)  # Both zenith terms are required.
    expected = 8.52 + 0.01 * (
        3.0 * (math.cos(math.radians(60)) - math.cos(math.radians(20)))
        + 4.0 * (math.sin(math.radians(60)) - math.sin(math.radians(20)))
    )
    assert plane([2.0, 20.0, 3.0, 4.0]) == pytest.approx(expected)
    with pytest.raises(ModelError, match="zenith"):
        plane([2.0, 20.0, 3.0, 4.0], zenith=None)


def test_camera_body_and_transmission():
    model = build(segmented())
    camera = model["optics"]["obscurers"]
    assert len(camera) == 1 and camera[0]["shape"] == "plate" and camera[0]["legs"] == ["incoming"]
    assert camera[0]["centre_m"][2] == pytest.approx(8.5) and camera[0]["aperture"] == "square"
    transmission = model["optics"]["telescope_transmission"]
    assert transmission == {
        "on_axis": 0.9,
        "amplitude": 0.1,
        "angular_scale_rad": math.radians(5.0),
        "power": 2.0,
        "outer_power": 1.0,
    }
    assert not model["provenance"]["structure"]


def test_structure_replaces_the_telescope_transmission():
    model = build(segmented(telescope_structure=STRUCTURE))
    optics = model["optics"]
    assert optics["telescope_transmission"] is None
    assert model["provenance"]["structure"]
    assert any("telescope_transmission is not applied" in note for note in model["notes"])
    mast, plate = optics["obscurers"][1:]
    assert mast["name"] == "mast-0" and mast["group"] == "camera_support"
    assert mast["legs"] == ["incoming", "after_primary", "after_secondary"]
    assert mast["start_m"][2] == pytest.approx(0.1 + 0.5)  # Dish frame -> telescope frame.
    assert plate["shape"] == "polyhedron" and len(plate["planes"]) == 4
    ids = [o["id"] for o in optics["obscurers"]] + [f["id"] for f in optics["primary"]["facets"]]
    assert len(ids) == len(set(ids)) and optics["focal_surface"]["id"] not in ids


@pytest.mark.parametrize(
    ("old", "new", "message"),
    [
        ("cylinder 0.1", "cylinder nan", "two end points and a diameter"),
        ("0.0 0.5 7.9 nan nan nan", "nan nan nan 0.0 0.5 7.9", "unused points must come last"),
        ("0.0 0.5 7.9 nan nan nan", "0.0 0.5 7.9 1.0 nan nan", "nan and numeric"),
        ("convex_polyhedron nan", "convex_polyhedron 0.2", "has no diameter"),
        ("convex_polyhedron", "box", "unsupported shape box"),
    ],
)
def test_structure_rows_are_checked(old, new, message):
    with pytest.raises(ModelError, match=message):
        build(segmented(telescope_structure=STRUCTURE.replace(old, new)))


def test_convex_hull():
    cube = [[x, y, z] for x in (0, 1) for y in (0, 1) for z in (0, 1)]
    planes = convex_hull(cube)
    assert len(planes) == 6
    for nx, ny, nz, d in planes:
        assert max(nx * x + ny * y + nz * z for x, y, z in cube) == pytest.approx(d)
    with pytest.raises(ModelError, match="coplanar"):
        convex_hull([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])


def test_response_grid():
    table = {
        "x": [300.0, 400.0],
        "y": [],
        "values": [0.5, 0.9],
        "interpolation": {"boundary": "clamp", "scheme": 1},
    }
    grid = response_grid(table, factor=0.5)
    assert grid["wavelength_min_nm"] == 300.0 and len(grid["values"]) == 101
    assert grid["values"][50] == [pytest.approx(0.35)]


def test_dual_mirror_telescope():
    model = build(dual(), zenith_deg=20.0)
    optics = model["optics"]
    primary, secondary = optics["primary"]["asphere"], optics["secondary"]["asphere"]
    assert primary["vertex_z_m"] == 0.0 and secondary["vertex_z_m"] == pytest.approx(3.1084)
    # The secondary frame is rotated by pi: its curvature flips sign.
    assert secondary["coefficients_m"][1] == pytest.approx(-0.00229358 * 100)
    assert optics["focal_surface"]["asphere"]["vertex_z_m"] == pytest.approx(2.5888)
    assert not optics["secondary_blocks_incoming"]
    names = [o["name"] for o in optics["obscurers"]]
    assert names == ["secondary-shadow", "camera-back", "camera-side-wall", "secondary-baffle"]
    side = optics["obscurers"][2]
    assert side["start_m"][2] == pytest.approx(2.0888) and side["end_m"][2] == pytest.approx(2.5888)
    assert optics["telescope_transmission"]["outer_power"] == 1.0
    assert model["sampling"]["disk_radius_m"] == pytest.approx(1.2 * 2.1205)
    single_unit = dict(dual()[0]["primary_mirror_parameters"], unit="cm")
    same = build(dual(primary_mirror_parameters=single_unit), zenith_deg=20.0)["optics"]
    assert same["primary"] == optics["primary"]


def test_secondary_baffle_wall_thickness_or_outer_radius():
    def baffle(values):
        optics = build(
            dual(secondary_mirror_baffle=record("secondary_mirror_baffle", values, "cm"))
        )
        return optics["optics"]["obscurers"][-1]

    wall = baffle([294.8, 273.6, 97.0, 10.0, 0.0])
    assert wall["thickness_m"] == pytest.approx(0.10)
    assert wall["start_radius_m"] == wall["end_radius_m"] == pytest.approx(0.97)
    # A fourth value not below R1 is the outer radius.
    assert baffle([294.8, 273.6, 97.0, 107.0, 0.0])["thickness_m"] == pytest.approx(0.10)
    cone = baffle([294.8, 273.6, 97.0, 5.0, 90.0])
    assert cone["end_radius_m"] == pytest.approx(0.90) and cone["thickness_m"] == pytest.approx(
        0.05
    )


def test_secondary_reflects_like_the_primary_without_its_own_table():
    table = "wavelength reflectivity\n300 0.8\n500 0.9\n"
    reflectivity = record("mirror_reflectivity", "ref.dat", file=True)
    records, _ = dual(mirror_reflectivity=reflectivity)
    files = {"mirror_reflectivity": table}
    for secondary in (None, record("secondary_mirror_reflectivity", "same", file=True)):
        if secondary:
            records["secondary_mirror_reflectivity"] = secondary
        optics = build((records, files), zenith_deg=20.0)["optics"]
        assert optics["secondary"]["reflectivity"] == optics["primary"]["reflectivity"]
        assert optics["primary"]["reflectivity"]["values"][100] == [pytest.approx(0.85)]


def test_degradation_maps_accept_filename_options():
    spatial_map = "#@RPOL@ 2\n 0 100\n0 1 0.9\n50 0.8 0.7\n"
    primary_record = record("primary_mirror_degraded_map", "map.dat#rpol:clip=1", file=True)
    single_records, files = segmented(primary_mirror_degraded_map=primary_record)
    files["primary_mirror_degraded_map"] = spatial_map
    model = build((single_records, files))
    assert model["optics"]["primary"]["degradation"]["outside"] == "zero"

    primary_record = record("primary_mirror_degraded_map", "map.dat#rpol:clip=1", file=True)
    secondary_record = record("secondary_mirror_degraded_map", "map.dat#rpol:clip=1", file=True)
    dual_records, files = dual(
        primary_mirror_degraded_map=primary_record,
        secondary_mirror_degraded_map=secondary_record,
    )
    files["primary_mirror_degraded_map"] = spatial_map
    files["secondary_mirror_degraded_map"] = spatial_map
    optics = build((dual_records, files), zenith_deg=20.0)["optics"]
    assert optics["primary"]["degradation"]["outside"] == "zero"
    assert optics["secondary"]["degradation"]["outside"] == "zero"


def test_rejected_models():
    with pytest.raises(ModelError, match="mirror_class 1"):
        build(segmented(mirror_class=record("mirror_class", 1)))
    with pytest.raises(ModelError, match="missing parameter focal_length"):
        build(segmented(focal_length=None))
    with pytest.raises(ModelError, match="seed"):
        build(segmented(), configuration_seed=-1)
