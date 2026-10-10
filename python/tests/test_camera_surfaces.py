"""Analytic finite camera-plane placement and frame conventions."""

import math
import unittest

from obdeect.camera_config import CameraConfigError, parse_camera_pixel_types
from obdeect.camera_surfaces import compile_camera_surfaces
from obdeect.optical_model_compiler import build_trace_model


class TestCameraSurfaces(unittest.TestCase):
    def test_projected_assignment_grid_uses_source_bounds_and_focal_plane(self):
        camera = self.camera([1])
        camera["rotation_deg"] = 90
        camera["pixels"][0]["z_offset_m"] = 0.2
        planes = compile_camera_surfaces(
            camera, dict(coefficient_m=[10], radial_scale_m=1), 1, reflected=True
        )
        grid = planes["assignment_grid"]
        self.assertEqual((grid["nx"], grid["ny"]), (4, 4))
        self.assertEqual((grid["x_low_m"], grid["x_high_m"]), (0, 1.02))
        self.assertEqual((grid["y_low_m"], grid["y_high_m"]), (-0.02, 0.02))
        self.assertEqual(grid["reference_plane_z_m"], 10)
        self.assertAlmostEqual(grid["x_basis"][1], 1)
        self.assertAlmostEqual(grid["y_basis"][0], -1)
        self.assertEqual(planes["entrance_surfaces"][0]["assignment_radius_m"], 0.02)

    def test_global_camera_rotation_follows_inverse_reference_frame(self):
        for reflected, expected_y in ((False, 1), (True, 1)):
            camera = self.camera([1])
            camera["rotation_deg"] = 90
            plane = compile_camera_surfaces(
                camera, dict(coefficient_m=[10], radial_scale_m=1), 1, reflected=reflected
            )["entrance_surfaces"][0]
            self.assertAlmostEqual(plane["centre_m"][0], 0)
            self.assertAlmostEqual(plane["centre_m"][1], expected_y)
            self.assertAlmostEqual(plane["tangent"][0], 0)
            self.assertAlmostEqual(plane["tangent"][1], expected_y)

    def camera(self, positions):
        types = parse_camera_pixel_types({
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
            ]
        })
        pixels = [
            dict(
                id=index,
                type_id=1,
                centre_xy_m=[x, 0],
                z_offset_m=0,
                rotation_deg=0,
                normal_slopes=[0, 0],
                module=2,
                enabled=True,
            )
            for index, x in enumerate(positions)
        ]
        return dict(pixel_types=types, pixels=pixels)

    def test_prime_focus_frame_preserves_transverse_position_and_physical_depth(self):
        camera = self.camera([1])
        camera["pixels"][0]["z_offset_m"] = 0.2
        planes = compile_camera_surfaces(
            camera, dict(coefficient_m=[10], radial_scale_m=1), 1, reflected=True
        )
        entrance, cathode = planes["entrance_surfaces"][0], planes["cathode_surfaces"][0]
        self.assertEqual(entrance["centre_m"], [1, 0, 9.8])
        self.assertEqual(entrance["normal"], [0, 0, -1])
        self.assertEqual(entrance["response_y_sign"], -1)
        self.assertAlmostEqual(cathode["centre_m"][2], 9.9)
        self.assertEqual(entrance["shape"], "square")
        self.assertEqual(cathode["shape"], "circle")

    def test_individual_and_common_module_planes_follow_analytic_paraboloid(self):
        camera = self.camera([0, 2, 4])
        focal = dict(coefficient_m=[10, 0.5], radial_scale_m=1)
        individual = compile_camera_surfaces(camera, focal, 0, reflected=False)["entrance_surfaces"]
        self.assertEqual([p["centre_m"][2] for p in individual], [10, 12, 18])
        self.assertAlmostEqual(individual[1]["normal"][0], -2 / math.sqrt(5))
        common = compile_camera_surfaces(camera, focal, 2, reflected=False)["entrance_surfaces"]
        self.assertTrue(all(p["normal"] == common[0]["normal"] for p in common))
        for plane in common:
            x, y, z = plane["centre_m"]
            self.assertAlmostEqual(-2 * x + z, 40 / 3 - 4)
            self.assertAlmostEqual(
                sum(a * b for a, b in zip(plane["normal"], plane["tangent"], strict=True)), 0
            )
        parallel = compile_camera_surfaces(camera, focal, 3, reflected=False)["entrance_surfaces"]
        self.assertTrue(all(p["normal"] == [0, 0, 1] for p in parallel))
        self.assertTrue(all(math.isclose(p["centre_m"][2], 40 / 3) for p in parallel))

    def test_incomplete_or_conflicting_alignment_is_rejected(self):
        camera = self.camera([0])
        camera["pixels"][0]["normal_slopes"] = [1, 0]
        focal = dict(coefficient_m=[10], radial_scale_m=1)
        with self.assertRaisesRegex(CameraConfigError, "conflicts"):
            compile_camera_surfaces(camera, focal, 1, reflected=False)
        del camera["pixels"][0]["z_offset_m"]
        with self.assertRaisesRegex(CameraConfigError, "requires"):
            compile_camera_surfaces(camera, focal, 0, reflected=False)

    def test_disabled_pixels_keep_module_placement_but_emit_no_surfaces(self):
        focal = dict(coefficient_m=[10, 0.5], radial_scale_m=1)
        for mode in range(4):
            for reflected in (False, True):
                with self.subTest(mode=mode, reflected=reflected):
                    camera = self.camera([0, 2, 4])
                    complete = compile_camera_surfaces(camera, focal, mode, reflected=reflected)
                    camera["pixels"][1]["enabled"] = False
                    active = compile_camera_surfaces(camera, focal, mode, reflected=reflected)
                    for key in ("entrance_surfaces", "cathode_surfaces"):
                        self.assertEqual(active[key], [complete[key][0], complete[key][2]])
                    self.assertFalse(camera["pixels"][1]["enabled"])

    def test_no_enabled_entrances_is_rejected(self):
        camera = self.camera([0])
        camera["pixels"][0]["enabled"] = False
        for mode in range(4):
            with (
                self.subTest(mode=mode),
                self.assertRaisesRegex(CameraConfigError, "no enabled pixel entrances"),
            ):
                compile_camera_surfaces(
                    camera, dict(coefficient_m=[10], radial_scale_m=1), mode, reflected=False
                )

    def test_trace_models_use_only_enabled_physical_entrances(self):
        camera = self.camera([0, 2])
        camera["pixels"][0]["enabled"] = False
        camera.update(
            compile_camera_surfaces(
                camera, dict(coefficient_m=[10], radial_scale_m=1), 1, reflected=False
            )
        )
        surface = dict(
            coefficient_m=[0.0] * 13,
            inner_radius_m=0.0,
            outer_radius_m=4.0,
            radial_scale_m=1.0,
        )
        optical_model = {
            "camera": camera,
            "report": {"facet_geometry_evidence": {"normal_status": "nominal_unperturbed"}},
            "primary": {
                "facets": [
                    dict(
                        id=0,
                        shape="circle",
                        diameter_m=1.0,
                        focal_length_m=10.0,
                        nominal_centre_m=[0.0, 0.0, 0.0],
                        nominal_normal=[0.0, 0.0, 1.0],
                    )
                ]
            },
        }
        segmented = build_trace_model(optical_model)
        response = [
            {"wavelength_nm": 300.0, "response": 0.8},
            {"wavelength_nm": 500.0, "response": 0.8},
        ]
        optical_model["primary"] = {"aspheric_surface": surface, "reflectivity": response}
        optical_model["secondary"] = {
            "kind": "aspheric_mirror",
            "reflectivity": response,
            **surface,
        }
        optical_model["focal_surface"] = surface
        axisymmetric = build_trace_model(optical_model)
        for trace in (segmented, axisymmetric):
            with self.subTest(kind=trace["kind"]):
                detectors = trace["detector_surfaces"]
                self.assertEqual(len(detectors), 1)
                self.assertEqual(detectors[0]["centre_m"], [2, 0, 10])
                self.assertTrue(detectors[0].get("enabled", True))
