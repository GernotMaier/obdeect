"""Checks for reproducible sim_telarray radial PSF products."""

import csv
import gzip
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))

from obdeect.simtel_reference_psf import (
    analyse_imaging_list,
    integration_radius,
    update_summary,
    write_profile,
)


class TestSimtelReferencePsf(unittest.TestCase):
    def test_rotation_counts_pixel_misses_and_exact_containment(self):
        data = (
            b"# Camera rotation angle = 90 deg\n"
            b"# Telescope 1 with 4 photons from 1 star(s) falling on an area of 10 m^2\n"
            b"1 -1 1 0\n"
            b"1 2 0 1\n"
        )
        result = analyse_imaging_list(data)
        self.assertEqual(result["focal_plane_crossings"], 2)
        self.assertEqual(result["launched_photons"], 4)
        self.assertAlmostEqual(result["effective_area_m2"], 5)
        self.assertAlmostEqual(result["centroid_x_m"], -0.005)
        self.assertAlmostEqual(result["centroid_y_m"], 0.005)
        self.assertAlmostEqual(integration_radius(result["radii_cm"], 0.8), 2**-0.5)

    def test_profile_reaches_every_crossing(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "profile.csv"
            write_profile(path, [1.0, 2.0, 3.0, 4.0])
            with path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1001)
        self.assertEqual(rows[0]["enclosed_fraction"], "0.0")
        self.assertEqual(rows[-1]["integration_radius_cm"], "4.0")
        self.assertEqual(rows[-1]["enclosed_fraction"], "1.0")

    def test_check_rejects_modified_cumulative_product(self):
        raw = (
            b"# Camera rotation angle = 0 deg\n"
            b"# Telescope 1 with 4 photons from 1 star(s) falling on an area of 10 m^2\n"
            b"1 -1 0 0\n1 1 1 0\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / "photons.gz"
            archive.write_bytes(gzip.compress(raw, mtime=0))
            result = analyse_imaging_list(raw)
            summary = root / "summary.json"
            summary.write_text(
                json.dumps({
                    "rows": [
                        {
                            "imaging_list_archive": archive.name,
                            "cumulative_profile": "profile.csv",
                            "focal_plane_crossings": result["focal_plane_crossings"],
                            "centroid_x_m": result["centroid_x_m"],
                            "centroid_y_m": result["centroid_y_m"],
                            "effective_area_m2": result["effective_area_m2"],
                            "sha256": {
                                "imaging_list": hashlib.sha256(raw).hexdigest(),
                                "imaging_list_archive": hashlib.sha256(
                                    archive.read_bytes()
                                ).hexdigest(),
                            },
                        }
                    ]
                })
            )
            update_summary(summary, root)
            update_summary(summary, root, check=True)
            (root / "profile.csv").write_text("corrupt\n")
            with self.assertRaisesRegex(ValueError, "radial PSF product differs"):
                update_summary(summary, root, check=True)


if __name__ == "__main__":
    unittest.main()
