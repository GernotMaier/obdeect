"""Tests for weighted focal-space and geometric-time distributions."""

import tempfile
import unittest
from pathlib import Path

from obdeect.photon_distributions import photon_distributions, write_distribution_svg
from obdeect.result_contract import OpticalArrival


class TestPhotonDistributions(unittest.TestCase):
    def test_weighted_space_and_time_are_conserved(self):
        arrivals = [
            OpticalArrival(1, "star", 400, 2, 2, 0.5, "detected", 3, 1, -1, 0, ((1, -1, 0),)),
            OpticalArrival(2, "laser", 400, 4, 1, 1, "detected", 6, -1, 1, 0, ((-1, 1, 0),)),
            OpticalArrival(
                3,
                "illuminator",
                400,
                7,
                10,
                0.2,
                "missed_primary",
                9,
                None,
                None,
                None,
                ((0, 0, 0),),
            ),
        ]
        result = photon_distributions(arrivals, bins=4)

        self.assertEqual(result.detected_count, 2)
        self.assertEqual(result.detected_weight, 2.0)
        self.assertEqual(sum(result.focal_x_m.bins), 2.0)
        self.assertEqual(sum(result.focal_y_m.bins), 2.0)
        self.assertEqual(sum(result.arrival_time_ns.bins), 2.0)
        self.assertEqual(sum(sum(row) for row in result.focal_plane_bins), 2.0)
        self.assertGreater(result.arrival_time_ns.upper, 4 + 6 / 0.299792458)

    def test_svg_is_self_contained_and_labels_the_time_convention(self):
        arrival = OpticalArrival(1, "star", 400, 0, 1, 1, "detected", 1, 0, 0, 0, ((0, 0, 0),))
        result = photon_distributions([arrival], bins=2)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "distributions.svg"
            write_distribution_svg(result, output)
            contents = output.read_text(encoding="utf-8")
        self.assertIn("weighted focal-plane distribution", contents)
        self.assertIn("vacuum geometric arrival time", contents)
        self.assertIn("<svg", contents)
        self.assertNotIn("/>width=", contents)


if __name__ == "__main__":
    unittest.main()
