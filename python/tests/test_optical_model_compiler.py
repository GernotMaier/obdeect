"""Tests for the generic simulation-models optical model compiler."""

import json
import math
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from obdeect.model_import import resolve_model
from obdeect.observing_geometry import resolve_observing_geometry
from obdeect.optical_model_compiler import (
    OpticalModelCompileError,
    _dual_reflector_surfaces,
    _secondary_segment_frame,
    apply_panel_alignment,
    build_plot_geometry,
    build_trace_model,
    compile_optical_model,
    derive_nominal_single_reflector,
    parse_model_segmentation,
    parse_obscuration_cylinders,
    parse_simtel_mirror_list,
    parse_simtel_segmentation,
    parse_wavelength_response,
    require_trace_ready,
    trace_surface_rows,
    write_native_optical_model,
)


class TestOpticalModelCompiler(unittest.TestCase):
    def test_observing_geometry_matches_fresh_compile_without_accumulating_errors(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_ir(root)
            production = root / "productions/1.0.0/GENERIC.json"
            manifest = json.loads(production.read_text())
            additions = {
                "focal_length": (1600, "cm", False),
                "dish_shape_length": (1600, "cm", False),
                "mirror_offset": (0, "cm", False),
                "mirror_class": (0, None, False),
                "parabolic_dish": (True, None, False),
                "mirror_align_random_horizontal": ([0.01, 20, 0.02, 0.03], ["deg"] * 4, False),
                "mirror_align_random_vertical": ([0.02, 20, 0.01, 0.04], ["deg"] * 4, False),
                "focus_offset": ([1, 20, 2, 3], ["cm", "deg", "", "null"], False),
                "camera_config_file": ("camera.dat", None, True),
                "camera_pixels": (1, None, False),
                "camera_degraded_map": ("camera-map.dat", None, True),
            }
            (root / "model_parameters/Files/camera.dat").write_text(
                "PixType 1 0 0 1 0 1 0 0.9 0.7\nPixel 0 1 0 0\n"
            )
            (root / "model_parameters/Files/camera-map.dat").write_text(
                "#@RPOL@ 3\n-1 -1 0.2\n-1 1 0.4\n1 -1 0.6\n1 1 0.8\n"
            )
            for name, (value, unit, is_file) in additions.items():
                path = root / f"model_parameters/GENERIC/{name}/{name}-1.0.0.json"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(
                    json.dumps({
                        "instrument": "GENERIC",
                        "parameter": name,
                        "parameter_version": "1.0.0",
                        "type": "string" if is_file else "float64",
                        "file": is_file,
                        "value": value,
                        "unit": unit,
                    })
                )
                manifest["parameters"]["GENERIC"][name] = "1.0.0"
            production.write_text(json.dumps(manifest))
            ir = resolve_model(root, "GENERIC", "1.0.0")
            base = compile_optical_model(ir, root, alignment_seed=19, alignment_zenith_deg=20)
            self.assertTrue(base["trace_model"]["camera_degradation_in_detector_frame"])
            self.assertEqual(base["trace_model"]["camera_degradation"]["x_basis"], [1, 0, 0])
            self.assertEqual(base["trace_model"]["detector_surfaces"][0]["response_y_sign"], -1)
            original = json.dumps(base, sort_keys=True)
            for angle in (17, 23, 0, -3):
                with self.subTest(angle=angle):
                    resolved = resolve_observing_geometry(base, angle)
                    fresh = compile_optical_model(
                        ir, root, alignment_seed=19, alignment_zenith_deg=angle
                    )
                    self.assertEqual(resolved["trace_model"], fresh["trace_model"])
                    self.assertEqual(json.dumps(base, sort_keys=True), original)
                    round_trip = resolve_observing_geometry(resolved, 20)
                    self.assertEqual(round_trip["trace_model"], base["trace_model"])

    def test_derived_reference_parameters_are_accounted_without_changing_transport(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = compile_optical_model(self.make_ir(root), root)
            production = root / "productions/1.0.0/GENERIC.json"
            manifest = json.loads(production.read_text())
            manifest["parameters"]["GENERIC"]["effective_focal_length"] = "1.0.0"
            production.write_text(json.dumps(manifest))
            path = root / (
                "model_parameters/GENERIC/effective_focal_length/effective_focal_length-1.0.0.json"
            )
            path.parent.mkdir()
            path.write_text(
                json.dumps({
                    "instrument": "GENERIC",
                    "parameter": "effective_focal_length",
                    "parameter_version": "1.0.0",
                    "type": "float64",
                    "unit": "cm",
                    "file": False,
                    "value": 1601,
                })
            )
            model = compile_optical_model(resolve_model(root, "GENERIC", "1.0.0"), root)
            report = model["report"]
            self.assertNotIn("effective_focal_length", report["deferred"])
            self.assertNotIn("effective_focal_length", report["consumed"])
            self.assertEqual(
                report["field_coverage"]["effective_focal_length"]["disposition"],
                "reference_diagnostic",
            )
            self.assertIn(
                "sha256",
                report["reference_diagnostics"]["effective_focal_length"]["parameter_record"],
            )
            self.assertEqual(model["primary"], baseline["primary"])

    def test_parses_model_obscuration_cylinders_in_metres(self):
        cylinders = parse_obscuration_cylinders(
            "# %ECSV 1.0\nid group x1 y1 z1 x2 y2 z2 diameter\nmast-1 mast 0 0 1 0 0 4 0.2\n"
        )
        self.assertEqual(cylinders[0]["id"], "mast-1")
        self.assertEqual(cylinders[0]["first_endpoint_m"], [0.0, 0.0, 1.0])
        self.assertEqual(cylinders[0]["diameter_m"], 0.2)

    def test_trace_model_detector_bound_includes_outer_pixel_entrance(self):
        optical_model = {
            "report": {"facet_geometry_evidence": {"normal_status": "nominal_unperturbed"}},
            "primary": {
                "facets": [
                    {
                        "id": 0,
                        "shape": "circle",
                        "diameter_m": 1.0,
                        "focal_length_m": 16.0,
                        "nominal_centre_m": [0.0, 0.0, 0.0],
                        "nominal_normal": [0.0, 0.0, 1.0],
                        "nominal_tangent": [1.0, 0.0, 0.0],
                    }
                ]
            },
            "camera": {
                "pixel_types": [{"id": 2, "funnel_diameter_m": 0.12}],
                "pixels": [{"type_id": 2, "centre_xy_m": [1.0, 0.0]}],
            },
            "focal_length_m": 16.0,
        }
        rows, _, _ = trace_surface_rows(optical_model)
        detector = rows[-1]
        self.assertAlmostEqual(detector["diameter_m"], 2.12)

    def test_dual_reflector_surfaces_preserve_si_aspheres_and_holes(self):
        coefficient = {"value": [0.0, 0.01], "unit": ["cm", "cm"]}
        parameters = {
            "primary_mirror_parameters": coefficient,
            "primary_mirror_diameter": {"value": 400.0, "unit": "cm"},
            "primary_mirror_hole_diameter": {"value": 20.0, "unit": "cm"},
            "secondary_mirror_parameters": {"value": [300.0, 0.02], "unit": ["cm", "cm"]},
            "secondary_mirror_diameter": {"value": 180.0, "unit": "cm"},
            "secondary_mirror_hole_diameter": {"value": 0.0, "unit": "cm"},
            "focal_surface_parameters": {"value": [200.0], "unit": ["cm"]},
        }
        surfaces = _dual_reflector_surfaces(parameters)
        self.assertIsNotNone(surfaces)
        self.assertEqual(surfaces["primary"]["inner_radius_m"], 0.1)
        self.assertEqual(surfaces["primary"]["outer_radius_m"], 2.0)
        self.assertAlmostEqual(surfaces["primary"]["coefficient_m"][1], 1.0)
        self.assertAlmostEqual(surfaces["secondary"]["coefficient_m"][0], 3.0)
        self.assertAlmostEqual(surfaces["secondary"]["coefficient_m"][1], -2.0)
        # The local secondary normal is reversed by sim_telarray's pi rotation.
        parameters["secondary_mirror_parameters"]["value"][0] = -300.0
        surfaces = _dual_reflector_surfaces(parameters)
        self.assertEqual(surfaces["secondary"]["coefficient_m"][:2], [3.0, -2.0])

    def test_dual_reflector_preserves_simtel_reference_radius_convention(self):
        parameters = {
            "primary_mirror_parameters": {"value": [0.0, 0.1], "unit": ["cm", "cm"]},
            "primary_mirror_ref_radius": {"value": 500.0, "unit": "cm"},
            "primary_mirror_diameter": {"value": 800.0, "unit": "cm"},
            "primary_mirror_hole_diameter": {"value": 0.0, "unit": "cm"},
            "secondary_mirror_parameters": {"value": [1.0], "unit": ["cm"]},
            "secondary_mirror_ref_radius": {"value": 500.0, "unit": "cm"},
            "secondary_mirror_diameter": {"value": 400.0, "unit": "cm"},
            "secondary_mirror_hole_diameter": {"value": 0.0, "unit": "cm"},
            "focal_surface_parameters": {"value": [2.0], "unit": ["cm"]},
            "focal_surface_ref_radius": {"value": 500.0, "unit": "cm"},
        }
        surfaces = _dual_reflector_surfaces(parameters)
        self.assertIsNotNone(surfaces)
        self.assertEqual(surfaces["primary"]["radial_scale_m"], 5.0)
        self.assertEqual(surfaces["primary"]["coefficient_m"][:2], [0.0, 0.5])
        self.assertEqual(surfaces["focal_surface"]["coefficient_m"][0], 10.0)

    def test_trace_model_contains_provenance_bound_surface(self):
        optical_model = {
            "provenance": {"model": "GENERIC", "model_version": "1.0.0", "input_records": {}},
            "optical_model_sha256": "a" * 64,
            "report": {"facet_geometry_evidence": {"normal_status": "nominal_unperturbed"}},
            "primary": {
                "facets": [
                    {
                        "id": 0,
                        "shape": "hexagon_flat_y",
                        "diameter_m": 1.2,
                        "focal_length_m": 16.0,
                        "nominal_centre_m": [0.0, 0.0, 0.0],
                        "nominal_normal": [0.0, 0.0, 1.0],
                        "nominal_tangent": [1.0, 0.0, 0.0],
                    }
                ]
            },
            "camera": {"pixels": [{"centre_xy_m": [0.0, 0.1]}]},
            "focal_length_m": 16.0,
        }
        trace_model = build_trace_model(optical_model)
        self.assertEqual(trace_model["kind"], "segmented")
        self.assertEqual(trace_model["primary_facets"][0]["shape"], "hexagon_flat_y")
        self.assertEqual(trace_model["detector_surfaces"][0]["shape"], "circle")

    def test_alignment_exact_frame_seed_and_order_independence(self):
        import copy

        facets = [{"id": 1, "centre_m": [2.0, 3.0, 0.0]}, {"id": 2, "centre_m": [-3.0, 1.0, 0.0]}]
        parameters = {
            "focal_length": {"value": 16.0, "unit": "m"},
            "dish_shape_length": {"value": 16.0, "unit": "m"},
            "mirror_offset": {"value": 0.0, "unit": "m"},
            "parabolic_dish": {"value": True},
        }
        derive_nominal_single_reflector(facets, parameters)
        facet = facets[0]
        phi = math.atan2(3, 2)
        inclination = math.acos(facet["nominal_normal"][2])
        expected = [
            math.cos(phi) ** 2 * math.cos(inclination) + math.sin(phi) ** 2,
            math.cos(phi) * math.sin(phi) * (math.cos(inclination) - 1),
            math.cos(phi) * math.sin(inclination),
        ]
        for actual, value in zip(facet["nominal_tangent"], expected):
            self.assertAlmostEqual(actual, value, places=14)
        errors = {
            name: {"value": [0.01, 28.0, 0.0, 0.0], "unit": ["deg", "deg", "null", "null"]}
            for name in ("mirror_align_random_horizontal", "mirror_align_random_vertical")
        }
        reversed_facets = copy.deepcopy(facets[::-1])
        nominal = copy.deepcopy(facets)
        apply_panel_alignment(facets, errors, 4, None)
        apply_panel_alignment(reversed_facets, errors, 4, None)
        self.assertEqual(facets, reversed_facets[::-1])
        self.assertNotEqual(facets[0]["nominal_normal"], nominal[0]["nominal_normal"])
        for facet in facets:
            self.assertAlmostEqual(
                sum(a * b for a, b in zip(facet["nominal_normal"], facet["nominal_tangent"])),
                0,
                places=14,
            )
        errors["mirror_align_random_horizontal"]["value"][2] = 0.01
        with self.assertRaisesRegex(OpticalModelCompileError, "alignment zenith"):
            apply_panel_alignment(nominal, errors, 4, None)
        apply_panel_alignment(nominal, errors, 4, 20)

    def test_segmented_housing_uses_incoming_planes_and_telescope_frame(self):
        optical_model = {
            "report": {"facet_geometry_evidence": {"normal_status": "nominal_unperturbed"}},
            "primary": {
                "facets": [
                    {
                        "id": 0,
                        "shape": "circle",
                        "diameter_m": 1,
                        "focal_length_m": 16,
                        "nominal_centre_m": [0, 0, 0],
                        "nominal_normal": [0, 0, 1],
                        "nominal_tangent": [1, 0, 0],
                    }
                ]
            },
            "camera": {
                "pixel_types": [{"id": 0, "funnel_diameter_m": 0.1}],
                "pixels": [{"type_id": 0, "centre_xy_m": [0, 0]}],
            },
            "focal_length_m": 16,
        }
        optical_model["camera"]["rotation_deg"] = 17
        optical_model["camera"]["housing"] = {
            "shape": "square",
            "diameter_m": 2,
            "front_z_m": 15,
            "depth_m": 1,
        }
        with TemporaryDirectory() as directory:
            with self.assertRaisesRegex(OpticalModelCompileError, "use JSON"):
                write_native_optical_model(optical_model, Path(directory) / "surfaces.csv")
        trace = build_trace_model(optical_model)
        planes = trace["incoming_obscurer_planes"]
        self.assertEqual([p["centre_m"][2] for p in planes], [15, 16])
        self.assertEqual(planes[0]["tangent"], [1, 0, 0])
        self.assertEqual(planes[0]["shape"], "square")
        components = build_plot_geometry(trace, {"trace_blockers": []})["components"]
        self.assertEqual(
            [p["id"] for p in planes], [c["id"] for c in components if c["role"] == "obscurer"]
        )

    def test_plot_index_contains_only_native_finite_components(self):
        trace = {
            "kind": "segmented",
            "primary_facets": [{"id": 4}],
            "detector_surfaces": [{"id": 7}],
            "cylinder_obscurers": [{"id": 9}],
        }
        index = build_plot_geometry(
            trace, {"trace_blockers": ["quadrilateral obscurers remain deferred"]}
        )
        self.assertEqual([item["id"] for item in index["components"]], [4, 7, 9])
        self.assertEqual(index["components"][2]["source"], "trace_model.cylinder_obscurers")
        self.assertEqual(index["unavailable_roles"], ["quadrilateral_obscurers"])
        dual = build_plot_geometry({"kind": "axisymmetric"}, {"trace_blockers": []})
        self.assertEqual(
            [item["role"] for item in dual["components"]], ["primary", "secondary", "detector"]
        )
        dual = build_plot_geometry(
            {
                "kind": "axisymmetric",
                "primary_to_secondary_planes": [{"id": 10}],
                "incoming_obscurer_planes": [{"id": 11}],
                "primary_to_secondary_cylinders": [{"id": 12}],
            },
            {"trace_blockers": []},
        )
        self.assertEqual(
            [(item["id"], item["role"]) for item in dual["components"][-3:]],
            [(10, "obscurer"), (11, "obscurer"), (12, "opaque_cylinder")],
        )

    def test_trace_readiness_requires_resolved_optical_model(self):
        require_trace_ready({"report": {"trace_blockers": [], "trace_ready": True}})
        with self.assertRaisesRegex(OpticalModelCompileError, "detector surfaces"):
            require_trace_ready({
                "report": {"trace_blockers": ["physical detector surfaces are unresolved"]}
            })
        with self.assertRaisesRegex(OpticalModelCompileError, "trace-model binding is unavailable"):
            require_trace_ready({"report": {"trace_blockers": []}})

    def test_nominal_panel_normal_points_chief_ray_at_focal_plane(self):
        # A panel at r=2 m on a 16 m DC dish retains its signed placement and
        # its nominal normal bisects the incident ray and focal-point ray.
        parameters = {
            "focal_length": {"value": 1600, "unit": "cm"},
            "dish_shape_length": {"value": 1600, "unit": "cm"},
            "mirror_offset": {"value": -100, "unit": "cm"},
            "parabolic_dish": {"value": False},
        }
        facets = [{"centre_m": [2.0, 0.0, 0.0]}]
        derive_nominal_single_reflector(facets, parameters)
        sag = 16 - math.sqrt(16**2 - 2**2)
        centre_z = sag + 1
        self.assertAlmostEqual(facets[0]["nominal_centre_m"][2], centre_z)
        focal_distance = math.hypot(2, 17 - centre_z)
        outgoing = [-2 / focal_distance, 0.0, (17 - centre_z) / focal_distance]
        normal_length = math.hypot(outgoing[0], outgoing[2] + 1)
        self.assertAlmostEqual(facets[0]["nominal_normal"][0], outgoing[0] / normal_length)
        self.assertAlmostEqual(facets[0]["nominal_normal"][2], (outgoing[2] + 1) / normal_length)

    def make_ir(self, root: Path, *, mirror_contents: str | None = None, focal_cm=1600.0) -> dict:
        asset = root / "model_parameters/Files/mirrors.dat"
        asset.parent.mkdir(parents=True)
        asset.write_text(
            mirror_contents
            or "# x y diameter focal shape z source metadata\n0 100 120 0 1 20 # id=M01\n"
        )
        (asset.parent / "filter.dat").write_text("wavelength transmission\n300 0.8\n400 0.9\n")
        production = root / "productions/1.0.0/GENERIC.json"
        production.parent.mkdir(parents=True)
        versions = {
            name: "1.0.0"
            for name in (
                "mirror_list",
                "mirror_focal_length",
                "camera_body_diameter",
                "camera_filter",
            )
        }
        production.write_text(
            json.dumps({
                "model_version": "1.0.0",
                "production_table_name": "GENERIC",
                "parameters": {"GENERIC": versions},
            })
        )
        values = {
            "mirror_list": ("mirrors.dat", None, True),
            "mirror_focal_length": (focal_cm, "cm", False),
            "camera_body_diameter": (200.0, "cm", False),
            "camera_filter": ("filter.dat", None, True),
        }
        for name, (value, unit, is_file) in values.items():
            parameter = root / f"model_parameters/GENERIC/{name}/{name}-1.0.0.json"
            parameter.parent.mkdir(parents=True)
            parameter.write_text(
                json.dumps({
                    "instrument": "GENERIC",
                    "parameter": name,
                    "parameter_version": "1.0.0",
                    "type": "string" if is_file else "float64",
                    "unit": unit,
                    "value": value,
                    "file": is_file,
                })
            )
        return resolve_model(root, "GENERIC", "1.0.0")

    def test_measured_scatter_is_enabled_by_default_and_preserves_units(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_ir(root)
            production = root / "productions/1.0.0/GENERIC.json"
            data = json.loads(production.read_text())
            name = "mirror_reflection_random_angle"
            data["parameters"]["GENERIC"][name] = "1.0.0"
            production.write_text(json.dumps(data))
            record = root / f"model_parameters/GENERIC/{name}/{name}-1.0.0.json"
            record.parent.mkdir(parents=True)
            parameter = {
                "instrument": "GENERIC",
                "parameter": name,
                "parameter_version": "1.0.0",
                "type": "float64",
                "unit": ["deg", "null", "deg"],
                "value": [0.0255, 0.2, 0.05],
                "file": False,
            }
            record.write_text(json.dumps(parameter))
            ir = resolve_model(root, "GENERIC", "1.0.0")
            nominal = compile_optical_model(ir, root)
            self.assertIn(name, nominal["report"]["consumed"])
            self.assertEqual(nominal["primary"]["scatter"]["seed"], 0)
            compiled = compile_optical_model(ir, root, scatter_seed=21)
            scatter = compiled["primary"]["scatter"]
            self.assertEqual(scatter["seed"], 21)
            self.assertEqual(scatter["method"], "outgoing_angles")
            self.assertAlmostEqual(scatter["sigma1_rad"], math.radians(0.0255))
            self.assertAlmostEqual(scatter["sigma2_rad"], math.radians(0.05))
            self.assertEqual(scatter["fraction2"], 0.2)
            self.assertNotIn(name, compiled["report"]["deferred"])
            parameter["unit"] = ["deg", "", "deg"]
            record.write_text(json.dumps(parameter))
            empty_unit = compile_optical_model(
                resolve_model(root, "GENERIC", "1.0.0"), root, scatter_seed=21
            )
            self.assertEqual(empty_unit["primary"]["scatter"], compiled["primary"]["scatter"])
            ir = resolve_model(root, "GENERIC", "1.0.0")
            for seed in (True, -1, 2**64):
                with self.subTest(seed=seed), self.assertRaises(OpticalModelCompileError):
                    compile_optical_model(ir, root, scatter_seed=seed)
            for values in ([90, 0, 0], [0, 1, 0], [-1, 0, 0], [0, 0, True]):
                parameter["value"] = values
                record.write_text(json.dumps(parameter))
                ir = resolve_model(root, "GENERIC", "1.0.0")
                with self.subTest(values=values), self.assertRaises(OpticalModelCompileError):
                    compile_optical_model(ir, root, scatter_seed=21)

    def test_axisymmetric_trace_preserves_scatter_housing_and_shadow(self):
        surface = {
            "inner_radius_m": 0,
            "outer_radius_m": 2,
            "radial_scale_m": 1,
            "coefficient_m": [0.0] * 13,
        }
        scatter = {
            "sigma1_rad": 0.001,
            "fraction2": 0,
            "sigma2_rad": 0,
            "method": "surface_slopes",
            "seed": 21,
        }
        optical_model = {
            "primary": {
                "aspheric_surface": surface,
                "scatter": scatter,
                "reflectivity": [
                    {"wavelength_nm": 300, "response": 0.8},
                    {"wavelength_nm": 500, "response": 0.8},
                ],
            },
            "secondary": {
                "kind": "aspheric_mirror",
                "reflectivity": [
                    {"wavelength_nm": 300, "response": 0.9},
                    {"wavelength_nm": 500, "response": 0.9},
                ],
                **surface,
                "scatter": scatter,
                "incoming_shadow": {"diameter_m": 3, "z_m": 3},
            },
            "focal_surface": surface,
            "camera": {
                "housing": {"shape": "square", "diameter_m": 0.5, "front_z_m": 1.5, "depth_m": 0.3},
            },
        }
        trace = build_trace_model(optical_model)
        self.assertEqual(trace["primary_scatter"], scatter)
        self.assertEqual(trace["secondary_scatter"], scatter)
        self.assertFalse(trace["block_incoming_secondary"])
        plane = trace["primary_to_secondary_planes"][0]
        self.assertEqual(plane["shape"], "square")
        self.assertEqual(plane["centre_m"], [0, 0, 1.2])
        cylinder = trace["primary_to_secondary_cylinders"][0]
        self.assertEqual(cylinder["second_endpoint_m"], [0, 0, 1.5])
        self.assertEqual(trace["incoming_obscurer_planes"][0]["diameter_m"], 3)
        ids = [
            trace[key]
            for key in ("primary_surface_id", "secondary_surface_id", "detector_surface_id")
        ] + [plane["id"], cylinder["id"], trace["incoming_obscurer_planes"][0]["id"]]
        self.assertEqual(len(ids), len(set(ids)))

    def test_imports_dual_mirror_scatter_and_obscuration_from_verified_records(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_ir(root)
            (root / "model_parameters/Files/segments.dat").write_text("hex 1 0 0 400 0\n")
            (root / "model_parameters/Files/reflectivity.dat").write_text(
                "wavelength reflectivity\n300 0.8\n400 0.9\n"
            )
            (root / "model_parameters/Files/camera.dat").write_text(
                'PixType 1 0 0 1 0 1 0 "filter.dat"\nPixel 0 1 0 0\n'
            )
            production = root / "productions/1.0.0/GENERIC.json"
            data = json.loads(production.read_text())
            parameters = {
                "primary_mirror_segmentation": ("segments.dat", None, True),
                "primary_mirror_parameters": ([0.0], ["cm"], False),
                "primary_mirror_diameter": (400, "cm", False),
                "primary_mirror_hole_diameter": (0, "cm", False),
                "secondary_mirror_parameters": ([250.0], ["cm"], False),
                "secondary_mirror_diameter": (180, "cm", False),
                "secondary_mirror_hole_diameter": (0, "cm", False),
                "focal_surface_parameters": ([200.0], ["cm"], False),
                "mirror_reflectivity": ("reflectivity.dat", None, True),
                "secondary_mirror_reflectivity": ("reflectivity.dat", None, True),
                "mirror_reflection_random_angle": ([0.0255, 0, 0], ["deg", "null", "deg"], False),
                "camera_config_file": ("camera.dat", None, True),
                "camera_pixels": (1, None, False),
                "camera_body_shape": (2, None, False),
                "camera_depth": (20, "cm", False),
                "secondary_mirror_shadow_diameter": (214, "cm", False),
                "secondary_mirror_shadow_offset": (0, "cm", False),
            }
            for name, (value, unit, is_file) in parameters.items():
                record = root / f"model_parameters/GENERIC/{name}/{name}-1.0.0.json"
                record.parent.mkdir(parents=True)
                record.write_text(
                    json.dumps({
                        "instrument": "GENERIC",
                        "parameter": name,
                        "parameter_version": "1.0.0",
                        "type": "string" if is_file else "float64",
                        "unit": unit,
                        "value": value,
                        "file": is_file,
                    })
                )
                data["parameters"]["GENERIC"][name] = "1.0.0"
            production.write_text(json.dumps(data))
            optical_model = compile_optical_model(
                resolve_model(root, "GENERIC", "1.0.0"), root, scatter_seed=21
            )
            trace = optical_model["trace_model"]
            self.assertEqual(trace["primary_scatter"]["method"], "surface_slopes")
            self.assertEqual(trace["secondary_scatter"], trace["primary_scatter"])
            self.assertAlmostEqual(trace["incoming_obscurer_planes"][0]["diameter_m"], 2.14)
            self.assertEqual(trace["incoming_obscurer_planes"][0]["centre_m"][2], 2.5)
            self.assertEqual(trace["primary_to_secondary_planes"][0]["shape"], "square")
            self.assertEqual(trace["primary_to_secondary_planes"][0]["centre_m"][2], 1.8)
            self.assertEqual(len(trace["primary_to_secondary_planes"]), 1)
            self.assertEqual(
                trace["primary_to_secondary_cylinders"][0]["second_endpoint_m"][2], 2.0
            )
            self.assertNotIn("camera_depth", optical_model["report"]["deferred"])
            shadow_path = root / (
                "model_parameters/GENERIC/secondary_mirror_shadow_diameter/"
                "secondary_mirror_shadow_diameter-1.0.0.json"
            )
            shadow_parameter = json.loads(shadow_path.read_text())
            shadow_parameter["value"] = -1
            shadow_path.write_text(json.dumps(shadow_parameter))
            automatic = compile_optical_model(resolve_model(root, "GENERIC", "1.0.0"), root)
            self.assertAlmostEqual(automatic["secondary"]["incoming_shadow"]["diameter_m"], 1.8)
            surface_path = root / (
                "model_parameters/GENERIC/secondary_mirror_parameters/"
                "secondary_mirror_parameters-1.0.0.json"
            )
            surface_parameter = json.loads(surface_path.read_text())
            surface_parameter.update(value=[250, -0.01], unit=["cm", "cm"])
            surface_path.write_text(json.dumps(surface_parameter))
            convex = compile_optical_model(resolve_model(root, "GENERIC", "1.0.0"), root)
            self.assertEqual(convex["secondary"]["incoming_shadow"]["z_m"], 2.5)
            with self.assertRaisesRegex(OpticalModelCompileError, "use compiled JSON"):
                write_native_optical_model(optical_model, root / "unsupported.csv")

    def test_compiles_tracked_mirror_list_without_dropping_deferred_fields(self):
        # T-IR-004: source units, shape, position and zero-focal fallback survive compilation.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            optical_model = compile_optical_model(self.make_ir(root), root)
        facet = optical_model["primary"]["facets"][0]
        self.assertEqual(optical_model["format"], "obdeect.compiled-optical-model.v1")
        self.assertEqual(facet["shape"], "hexagon_flat_x")
        self.assertEqual(facet["centre_m"], [0.0, 1.0, 0.2])
        self.assertEqual(facet["diameter_m"], 1.2)
        self.assertEqual(facet["focal_length_m"], 16.0)
        self.assertEqual(
            optical_model["report"]["deferred"], ["camera_body_diameter", "camera_filter"]
        )
        self.assertEqual(
            optical_model["report"]["facet_geometry_evidence"]["normal_status"], "unavailable"
        )
        self.assertIn(
            "No normals", optical_model["report"]["facet_geometry_evidence"]["interpretation"]
        )
        self.assertEqual(len(optical_model["optical_model_sha256"]), 64)

    def test_trace_model_preserves_complete_segmented_geometry(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            # Keep this focused on the trace-model contract; the production
            # compiler's camera parser is tested separately.
            optical_model = compile_optical_model(self.make_ir(root), root)
            optical_model["camera"] = {"pixels": [{"centre_xy_m": [0.1, 0.0]}]}
            optical_model["focal_length_m"] = 16.0
            optical_model["report"]["facet_geometry_evidence"]["normal_status"] = (
                "nominal_unperturbed"
            )
            for facet in optical_model["primary"]["facets"]:
                facet["nominal_centre_m"] = [*facet["centre_m"][:2], 0.0]
                facet["nominal_normal"] = [0.0, 0.0, 1.0]
                facet["nominal_tangent"] = [1.0, 0.0, 0.0]
            trace_model = build_trace_model(optical_model)
        self.assertEqual(trace_model["kind"], "segmented")
        self.assertEqual(trace_model["primary_facets"][0]["id"], 0)
        self.assertEqual(trace_model["detector_surfaces"][0]["shape"], "circle")

    def test_axisymmetric_export_preserves_outer_primary_reflectivity(self):
        surface = {
            "coefficient_m": [0.0] * 13,
            "inner_radius_m": 0.0,
            "outer_radius_m": 2.0,
            "radial_scale_m": 1.0,
        }
        optical_model = {
            "provenance": {"model": "GENERIC", "model_version": "1.0.0"},
            "optical_model_sha256": "a" * 64,
            "primary": {
                "aspheric_surface": surface,
                "reflectivity": [
                    {"wavelength_nm": 300.0, "response": 0.8},
                    {"wavelength_nm": 500.0, "response": 0.9},
                ],
            },
            "secondary": {
                "kind": "aspheric_mirror",
                **surface,
                "reflectivity": [
                    {"wavelength_nm": 300.0, "response": 0.7},
                    {"wavelength_nm": 500.0, "response": 0.8},
                ],
            },
            "focal_surface": surface,
        }
        with TemporaryDirectory() as directory:
            output = Path(directory) / "axisymmetric.csv"
            write_native_optical_model(optical_model, output)
            lines = output.read_text().splitlines()
        self.assertTrue(any(line.startswith("primary_reflectivity,300,0.8") for line in lines))
        self.assertTrue(any(line.startswith("primary_reflectivity,500,0.9") for line in lines))

        trace_model = build_trace_model(optical_model)
        self.assertEqual(
            [
                {k: v for k, v in row.items() if k not in ("interpolation", "wavelength_sampling")}
                for row in trace_model["primary_reflectivity"]
            ],
            optical_model["primary"]["reflectivity"],
        )
        self.assertEqual(
            [
                {k: v for k, v in row.items() if k not in ("interpolation", "wavelength_sampling")}
                for row in trace_model["secondary_reflectivity"]
            ],
            optical_model["secondary"]["reflectivity"],
        )

    def test_structured_segments_preserve_units_counts_and_secondary_gap_edge(self):
        parameter = {
            "value": [
                {
                    "kind": "ring",
                    "count": 2,
                    "r_min": {"value": 100, "unit": "cm"},
                    "r_max": {"value": 2, "unit": "m"},
                    "dphi": {"value": 40, "unit": "deg"},
                    "phi0": {"value": -10, "unit": "deg"},
                    "gap": {"value": 10, "unit": "mm"},
                }
            ]
        }
        segments = parse_model_segmentation(parameter)
        self.assertEqual(segments, parse_simtel_segmentation("ring 2 100 200 40 -10 1"))
        secondary = _secondary_segment_frame(segments)
        self.assertEqual([item["start_deg"] for item in secondary], [150, 110])
        self.assertTrue(all(item["gap_at_start"] for item in secondary))
        hexagon = parse_simtel_segmentation("hex 1 100 20 40 15")
        reflected = _secondary_segment_frame(hexagon)[0]
        self.assertEqual(reflected["centre_xy_m"], [-1.0, 0.2])
        self.assertEqual(reflected["rotation_deg"], 165)
        for key, value in (("count", True), ("extra", 1), ("kind", "unsupported")):
            invalid = {"value": [{**parameter["value"][0], key: value}]}
            with self.subTest(key=key), self.assertRaises(OpticalModelCompileError):
                parse_model_segmentation(invalid)

    def test_incidence_response_preserves_full_grid_and_rejects_missing_knots(self):
        header = "wavelength incidence_angle reflectivity\n"
        body = "500 20 0.6\n300 0 0.8\n500 0 0.9\n300 20 0.7\n"
        response = parse_wavelength_response(header + body, "mirror")
        self.assertEqual(
            [(row["wavelength_nm"], row["incidence_angle_deg"]) for row in response],
            [(300, 0), (300, 20), (500, 0), (500, 20)],
        )
        surface = {
            "coefficient_m": [0.0] * 13,
            "inner_radius_m": 0.0,
            "outer_radius_m": 2.0,
            "radial_scale_m": 1.0,
        }
        model = {
            "primary": {"aspheric_surface": surface, "reflectivity": response},
            "secondary": {"kind": "aspheric_mirror", **surface, "reflectivity": response},
            "focal_surface": surface,
        }
        trace = build_trace_model(model)
        self.assertEqual(
            [
                {k: v for k, v in row.items() if k not in ("interpolation", "wavelength_sampling")}
                for row in trace["primary_reflectivity"]
            ],
            response,
        )
        self.assertEqual(
            [
                {k: v for k, v in row.items() if k not in ("interpolation", "wavelength_sampling")}
                for row in trace["secondary_reflectivity"]
            ],
            response,
        )
        self.assertEqual(trace["secondary_reflectivity"][0]["interpolation"]["boundary"], "clamp")
        self.assertEqual(
            trace["primary_reflectivity"][0]["wavelength_sampling"],
            {
                "width_nm": 1.0,
                "offset_nm": 0,
                "first_bin": 200,
                "last_bin": 999,
                "projection_angle_deg": 0.0,
            },
        )
        measured = [{**row, "reflectivity_rms": 0.02} for row in response]
        model["primary"]["reflectivity"] = measured
        model["secondary"]["reflectivity"] = measured
        trace = build_trace_model(model)
        self.assertEqual(
            [
                {k: v for k, v in row.items() if k not in ("interpolation", "wavelength_sampling")}
                for row in trace["primary_reflectivity"]
            ],
            response,
        )
        self.assertEqual(trace["primary_reflectivity"][0]["interpolation"]["boundary"], "clamp")
        self.assertEqual(
            [
                {k: v for k, v in row.items() if k not in ("interpolation", "wavelength_sampling")}
                for row in trace["secondary_reflectivity"]
            ],
            response,
        )
        self.assertEqual(trace["secondary_reflectivity"][0]["interpolation"]["boundary"], "clamp")
        self.assertEqual(model["primary"]["reflectivity"][0]["reflectivity_rms"], 0.02)
        model["camera"] = {"filter_response": response}
        camera_response = build_trace_model(model)["camera_response"]
        self.assertEqual(
            camera_response["camera_filter"][0]["wavelength_sampling"],
            {"width_nm": 1.0, "offset_nm": 0.5, "first_bin": 0, "last_bin": 999},
        )
        self.assertEqual(
            [
                {k: v for k, v in row.items() if k not in ("interpolation", "wavelength_sampling")}
                for row in camera_response["camera_filter"]
            ],
            response,
        )
        self.assertEqual(camera_response["camera_filter"][0]["interpolation"]["boundary"], "clamp")
        self.assertEqual(camera_response["camera_transmission"], 1.0)
        self.assertNotIn("lightguide_efficiency", camera_response)
        model["camera"] = {"transmission": 0.9}
        self.assertEqual(build_trace_model(model)["camera_response"]["camera_transmission"], 0.9)
        for invalid in (
            body.replace("500 20 0.6\n", ""),
            body + "500 20 0.6\n",
            body.replace("20", "91"),
            body.replace("500", "300"),
        ):
            with self.subTest(invalid=invalid), self.assertRaises(OpticalModelCompileError):
                parse_wavelength_response(header + invalid, "mirror")
        with self.assertRaises(OpticalModelCompileError):
            parse_wavelength_response("wavelength reflectivity unsupported\n300 0.8 1\n", "mirror")

    def test_response_preserves_measurement_metadata(self):
        response = parse_wavelength_response(
            "wavelength reflectivity reflectivity_rms reflectivity_min reflectivity_max\n"
            "300 0.8 0.02 0.7 0.9\n500 0.6 0.03 0.5 0.7\n",
            "mirror",
        )
        self.assertEqual(response[0]["reflectivity_rms"], 0.02)
        self.assertEqual(response[1]["reflectivity_min"], 0.5)
        self.assertEqual(response[1]["reflectivity_max"], 0.7)
        with self.assertRaises(OpticalModelCompileError):
            parse_wavelength_response(
                "wavelength reflectivity reflectivity_rms\n300 0.8 invalid\n500 0.6 0.03\n",
                "mirror",
            )

    def test_axisymmetric_export_rejects_unrepresented_obscurers(self):
        surface = {
            "coefficient_m": [0.0] * 13,
            "inner_radius_m": 0.0,
            "outer_radius_m": 2.0,
            "radial_scale_m": 1.0,
        }
        optical_model = {
            "provenance": {"model": "GENERIC", "model_version": "1.0.0"},
            "optical_model_sha256": "a" * 64,
            "primary": {
                "aspheric_surface": surface,
                "cylinder_obscurers": [{"id": "mast"}],
            },
            "secondary": {"kind": "aspheric_mirror", **surface},
            "focal_surface": surface,
        }
        with TemporaryDirectory() as directory:
            with self.assertRaisesRegex(
                OpticalModelCompileError, "cannot represent cylinder obscurers"
            ):
                write_native_optical_model(optical_model, Path(directory) / "axisymmetric.csv")

    def test_parses_real_simtel_comment_suffix_without_turning_it_into_alignment(self):
        # T-IR-006: LST mirror-list rows carry an optional z followed by a
        # sim_telarray comment/panel ID; none of that is a facet orientation.
        facets = parse_simtel_mirror_list(
            "  1022.49 -462.00 151.00 2912.50 3 0.0 #% id=198\n"
            " -620.80 0.00 120.00 0.00 1 # no z supplied\n",
            fallback_focal_length_m=16.0,
        )
        self.assertEqual(facets[0]["centre_m"], [10.2249, -4.62, 0.0])
        self.assertEqual(facets[0]["shape"], "hexagon_flat_y")
        self.assertAlmostEqual(facets[1]["centre_m"][0], -6.208)
        self.assertEqual(facets[1]["centre_m"][1:], [0.0, 0.0])
        self.assertEqual(facets[1]["focal_length_m"], 16.0)
        self.assertNotIn("unit_normal", facets[0])
        self.assertNotIn("rotation_deg", facets[0])

    def test_parses_ecsv_mirror_list_with_float_shape_code(self):
        facets = parse_simtel_mirror_list(
            "# %ECSV 1.0\n"
            "mirror_x mirror_y mirror_diameter focal_length shape_type mirror_z mirror_panel_id\n"
            "0.0 100.0 120.0 1600.0 3.0 0.0 7\n",
            fallback_focal_length_m=None,
        )
        self.assertEqual(facets[0]["id"], 0)
        self.assertEqual(facets[0]["shape"], "hexagon_flat_y")
        self.assertEqual(facets[0]["focal_length_m"], 16.0)

    def test_rejects_invalid_optional_mirror_height(self):
        with self.assertRaisesRegex(OpticalModelCompileError, "invalid numeric field"):
            parse_simtel_mirror_list("0 0 120 1600 1 missing\n", fallback_focal_length_m=16.0)

    def test_zero_catalogue_fallback_allows_real_lst_rows_with_panel_focal_lengths(self):
        # T-IR-007: the LST catalogue sets mirror_focal_length to zero while
        # its mirror-list provides each panel's focal length.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            ir = self.make_ir(
                root,
                mirror_contents="1022.49 -462.00 151.00 2912.50 3 0.0 #% id=198\n",
                focal_cm=0.0,
            )
            optical_model = compile_optical_model(ir, root)
        self.assertEqual(optical_model["primary"]["facets"][0]["focal_length_m"], 29.125)

    def test_accepts_repository_root_containing_simulation_models_data_package(self):
        # T-IR-008: importer and compiler accept the same released-checkout layout.
        with TemporaryDirectory() as directory:
            checkout = Path(directory) / "simulation-models-repository"
            data_root = checkout / "simulation-models"
            optical_model = compile_optical_model(self.make_ir(data_root), checkout)
        self.assertEqual(len(optical_model["primary"]["facets"]), 1)

    def test_hash_mismatch_and_malformed_records_fail_closed(self):
        # T-IR-005: asset provenance and input syntax are validated before use.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            ir = self.make_ir(root)
            (root / "model_parameters/Files/mirrors.dat").write_text("0 0 120 0 1\n")
            with self.assertRaisesRegex(OpticalModelCompileError, "IR assets differs"):
                compile_optical_model(ir, root)
        with self.assertRaisesRegex(OpticalModelCompileError, "unsupported shape"):
            parse_simtel_mirror_list("0 0 120 1600 9\n", fallback_focal_length_m=16.0)

    def test_deferred_asset_and_parameter_records_are_verified(self):
        # T-IR-009: deferred input cannot change without invalidating provenance.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            ir = self.make_ir(root)
            (root / "model_parameters/Files/filter.dat").write_text("300 0.1\n")
            with self.assertRaisesRegex(OpticalModelCompileError, "IR assets differs"):
                compile_optical_model(ir, root)
            (root / "model_parameters/Files/filter.dat").write_text(
                "wavelength transmission\n300 0.8\n400 0.9\n"
            )
            ir["parameters"]["camera_body_diameter"]["value"] = 999.0
            with self.assertRaisesRegex(OpticalModelCompileError, "IR parameters differs"):
                compile_optical_model(ir, root)

    def test_parses_explicit_hex_and_ring_footprints(self):
        # T-IR-010: ring groups expand to stable IDs; geometry remains 2D.
        segments = parse_simtel_segmentation(
            "hex 1 -85.6 0 84.6 0\nRING 2 100 200 180 -90 1.4 # cm and degrees\n"
        )
        self.assertEqual([segment["id"] for segment in segments], [0, 1, 2])
        self.assertEqual(segments[0]["centre_xy_m"], [-0.856, 0.0])
        self.assertEqual(segments[1]["inner_radius_m"], 1.0)
        self.assertEqual(segments[2]["start_deg"], 90.0)
        self.assertAlmostEqual(segments[1]["gap_m"], 0.014)
        with self.assertRaisesRegex(OpticalModelCompileError, "polygon segment"):
            parse_simtel_segmentation("polygon 1 0 0 1 0\n")

    def test_defaults_omitted_segmentation_rotation_start_and_gap_to_zero(self):
        segments = parse_simtel_segmentation("hex 1 -85.6 0 84.6\nring 2 100 200 180\n")
        self.assertEqual(segments[0]["rotation_deg"], 0.0)
        self.assertEqual(segments[1]["start_deg"], 0.0)
        self.assertEqual(segments[2]["start_deg"], 180.0)
        self.assertEqual(segments[1]["gap_m"], 0.0)

    def test_partial_ring_spacing_and_yhex_orientation_follow_reference_parser(self):
        segments = parse_simtel_segmentation("ring 2 100 200 60 10 1\nyhex 1 0 0 80 5\n")
        self.assertEqual(segments[0]["start_deg"], 10.0)
        self.assertEqual(segments[1]["start_deg"], 70.0)
        self.assertEqual(segments[2]["rotation_deg"], 95.0)

    def test_compiles_dual_mirror_segmentation_without_invented_normals(self):
        # T-IR-011: nullable mirror_list uses explicit segmentation assets.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_ir(root)
            parameter = root / "model_parameters/GENERIC/mirror_list/mirror_list-1.0.0.json"
            data = json.loads(parameter.read_text())
            data["value"] = None
            parameter.write_text(json.dumps(data))
            production = root / "productions/1.0.0/GENERIC.json"
            data = json.loads(production.read_text())
            for name, contents in (
                ("primary_mirror_segmentation", "hex 1 0 0 80 0\n"),
                ("secondary_mirror_segmentation", "RING 2 10 20 180 0 0\n"),
            ):
                (root / f"model_parameters/Files/{name}.dat").write_text(contents)
                record = root / f"model_parameters/GENERIC/{name}/{name}-1.0.0.json"
                record.parent.mkdir(parents=True)
                record.write_text(
                    json.dumps({
                        "instrument": "GENERIC",
                        "parameter": name,
                        "parameter_version": "1.0.0",
                        "type": "string",
                        "unit": None,
                        "value": f"{name}.dat",
                        "file": True,
                    })
                )
                data["parameters"]["GENERIC"][name] = "1.0.0"
            production.write_text(json.dumps(data))
            ir = resolve_model(root, "GENERIC", "1.0.0")
            optical_model = compile_optical_model(ir, root)
        self.assertEqual(optical_model["primary"]["kind"], "segmented_footprints")
        self.assertEqual(len(optical_model["primary"]["segments"]), 1)
        self.assertEqual(len(optical_model["secondary"]["segments"]), 2)
        self.assertEqual(
            optical_model["report"]["facet_geometry_evidence"]["normal_status"], "unavailable"
        )

    def test_camera_layout_count_and_nested_response_provenance(self):
        # T-IR-012: camera channels are counted from the file, not only the record.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_ir(root)
            files = root / "model_parameters/Files"
            (files / "response.dat").write_text("300 0.8\n")
            (files / "camera.dat").write_text(
                'PixType 1 0 2 0.6 2 0.7 0.1 "response.dat"\nPixel 0 1 0 0\nPixel 1 1 1 0\n'
            )
            production = root / "productions/1.0.0/GENERIC.json"
            data = json.loads(production.read_text())
            for name, value, is_file in (
                ("camera_config_file", "camera.dat", True),
                ("camera_pixels", 2, False),
            ):
                record = root / f"model_parameters/GENERIC/{name}/{name}-1.0.0.json"
                record.parent.mkdir(parents=True)
                record.write_text(
                    json.dumps({
                        "instrument": "GENERIC",
                        "parameter": name,
                        "parameter_version": "1.0.0",
                        "type": "string" if is_file else "int64",
                        "unit": None,
                        "value": value,
                        "file": is_file,
                    })
                )
                data["parameters"]["GENERIC"][name] = "1.0.0"
            production.write_text(json.dumps(data))
            ir = resolve_model(root, "GENERIC", "1.0.0")
            optical_model = compile_optical_model(ir, root)
            self.assertEqual(optical_model["report"]["camera_layout_evidence"]["pixel_count"], 2)
            self.assertIn("response.dat", optical_model["provenance"]["nested_assets"])
            original_asset = optical_model["provenance"]["nested_assets"]["response.dat"]
            # Physical dimensions override the layout dimensions, but the
            # response file remains bound to the same pixel type and provenance.
            physical_types = root / (
                "model_parameters/GENERIC/camera_pixel_types/camera_pixel_types-1.0.0.json"
            )
            physical_types.parent.mkdir(parents=True)
            physical_types.write_text(
                json.dumps({
                    "instrument": "GENERIC",
                    "parameter": "camera_pixel_types",
                    "parameter_version": "1.0.0",
                    "type": "dict",
                    "unit": None,
                    "file": False,
                    "value": [
                        {
                            "type_id": 1,
                            "pmt_type": 0,
                            "cathode_shape": 0,
                            "cathode_diameter_cm": 2,
                            "funnel_shape": 2,
                            "funnel_diameter_cm": 4,
                            "funnel_depth_cm": 10,
                        }
                    ],
                })
            )
            data["parameters"]["GENERIC"]["camera_pixel_types"] = "1.0.0"
            production.write_text(json.dumps(data))
            ir = resolve_model(root, "GENERIC", "1.0.0")
            overridden = compile_optical_model(ir, root)
            pixel_type = overridden["camera"]["pixel_types"][0]
            self.assertEqual(pixel_type["funnel_diameter_m"], 0.04)
            self.assertEqual(pixel_type["response_files"], ["response.dat"])
            self.assertEqual(
                overridden["provenance"]["nested_assets"]["response.dat"], original_asset
            )
            (files / "response.dat").write_text("300 0.5\n")
            changed = compile_optical_model(ir, root)
            self.assertNotEqual(
                changed["provenance"]["nested_assets"]["response.dat"]["sha256"],
                original_asset["sha256"],
            )
            self.assertNotEqual(changed["optical_model_sha256"], overridden["optical_model_sha256"])
            (files / "response.dat").unlink()
            with self.assertRaisesRegex(OpticalModelCompileError, "missing or ambiguous"):
                compile_optical_model(ir, root)
            (files / "response.dat").write_text("300 0.5\n")
            count_record = root / "model_parameters/GENERIC/camera_pixels/camera_pixels-1.0.0.json"
            data = json.loads(count_record.read_text())
            data["value"] = 3
            count_record.write_text(json.dumps(data))
            with self.assertRaisesRegex(
                OpticalModelCompileError, "focal-plane element count differs"
            ):
                compile_optical_model(resolve_model(root, "GENERIC", "1.0.0"), root)

    def test_rejects_nonfinite_fallback_focal_length(self):
        with self.assertRaisesRegex(OpticalModelCompileError, "finite and positive"):
            parse_simtel_mirror_list("0 0 120 0 1\n", fallback_focal_length_m=float("nan"))


if __name__ == "__main__":
    unittest.main()
