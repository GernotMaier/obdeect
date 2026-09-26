"""Camera layout parser preserves pixel identity and rejects malformed geometry."""

import unittest

from obdeect.camera_config import CameraConfigError, parse_camera_layout, parse_camera_layout_ecsv


class TestCameraConfig(unittest.TestCase):
    def test_ecsv_layout_preserves_focal_plane_coordinates(self):
        layout = parse_camera_layout_ecsv(
            "# %ECSV 1.0\npixel_id type_id x_cm y_cm enabled\n3 1 -2.0 1.5 1\n4 2 2.0 -1.5 1\n"
        )
        self.assertEqual(layout["kind"], "focal_plane_layout")
        self.assertEqual(layout["pixels"][0]["centre_xy_m"], [-0.02, 0.015])

    def test_pixel_layout_and_response_reference(self):
        layout = parse_camera_layout(
            'PixType 1 0 2 0.6 2 0.7 0.1 "response.dat"\n'
            "Rotate 10\nPixel 0 1 -2 3 7 1 2 0x01 1\n"
            "Pixel 1 1 4 -5 7 1 3 0x01 0\nMajorityTrigger 0 1\n"
        )
        self.assertEqual([pixel["id"] for pixel in layout["pixels"]], [0, 1])
        self.assertEqual(layout["pixels"][0]["centre_xy_m"], [-0.02, 0.03])
        self.assertEqual(layout["pixels"][1]["source_columns"][-1], "0")
        self.assertEqual(layout["pixel_types"][0]["response_files"], ["response.dat"])
        self.assertEqual(layout["rotation_deg"], 10.0)
        self.assertEqual(layout["deferred_directives"], {"MajorityTrigger": 1})

    def test_duplicate_and_unknown_directive_fail(self):
        prefix = "PixType 1 0 2 0.6 2 0.7 0.1 1.0\n"
        with self.assertRaisesRegex(CameraConfigError, "duplicate Pixel"):
            parse_camera_layout(prefix + "Pixel 0 1 0 0\nPixel 0 1 1 1\n")
        with self.assertRaisesRegex(CameraConfigError, "unsupported directive"):
            parse_camera_layout(prefix + "Mystery 1\nPixel 0 1 0 0\n")
        with self.assertRaisesRegex(CameraConfigError, "undefined PixType"):
            parse_camera_layout(prefix + "Pixel 0 2 0 0\n")


if __name__ == "__main__":
    unittest.main()
