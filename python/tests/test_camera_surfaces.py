"""Analytic finite camera-plane placement and frame conventions."""

import math
import unittest

from obdeect.camera_config import CameraConfigError, parse_camera_pixel_types
from obdeect.camera_surfaces import compile_camera_surfaces


class TestCameraSurfaces(unittest.TestCase):
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

    def test_prime_focus_pi_rotation_preserves_physical_depth(self):
        camera = self.camera([1])
        camera["pixels"][0]["z_offset_m"] = 0.2
        planes = compile_camera_surfaces(
            camera, dict(coefficient_m=[10], radial_scale_m=1), 1, reflected=True
        )
        entrance, cathode = planes["entrance_surfaces"][0], planes["cathode_surfaces"][0]
        self.assertEqual(entrance["centre_m"], [-1, 0, 9.8])
        self.assertEqual(entrance["normal"], [0, 0, -1])
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
