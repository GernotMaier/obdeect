"""Analytic ambient propagation tests independent of telescope geometry."""

import math
import unittest

from obdeect.propagation import ambient_group_index


class TestPropagation(unittest.TestCase):
    def test_exponential_refractivity_matches_reference_lookup(self):
        contents = "altitude refractive_index\n" + "\n".join(
            f"{height} {0.0003 * math.exp(-height / 8):.17g}" for height in (0, 1, 2, 3)
        )
        for altitude in (0, 123, 1000, 1789):
            self.assertAlmostEqual(
                ambient_group_index(contents, altitude),
                1 + 0.0003 * math.exp(-altitude / 8000),
                places=15,
            )
        self.assertEqual(ambient_group_index(contents, 3000), 1)
        with self.assertRaisesRegex(ValueError, "upper layer"):
            ambient_group_index(contents, 2500)

    def test_invalid_profile_rejects(self):
        for contents in (
            "altitude index\n0 1\n1 2",
            "altitude refractive_index\n0 0\n1 0.1",
            "altitude refractive_index\n0 0.1\n0 0.2",
        ):
            with self.assertRaises(ValueError):
                ambient_group_index(contents, 100)
