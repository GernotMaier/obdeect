"""Analytic input-boundary tests, independent of optical image analysis."""

import csv
import math
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from obdeect.optical_model_compiler import compile_telescope_frame, compile_telescope_transmission
from obdeect.telescope_frame import (
    prepare_ground_photons,
    resolve_telescope_frame,
    single_telescope_identity,
)


class TestTelescopeFrame(unittest.TestCase):
    def context(self):
        return dict(
            axis_origin_m=[10, 20, 30],
            axes_offsets_m=[2, 3],
            known_error_deg=0,
            unknown_error_deg=0,
        )

    def test_mount_offsets_at_zenith_and_horizon(self):
        frame = resolve_telescope_frame(self.context(), 0, 0)
        self.assertEqual(frame.to_local((10, 20, 30), (0, 0, -1)), ((1, 0, 0), (0, 0, -1)))
        frame = resolve_telescope_frame(self.context(), 0, 90)
        position, direction = frame.to_local((10, 20, 30), (-1, 0, 0))
        for actual, expected in zip(position, (3, 0, -2), strict=True):
            self.assertAlmostEqual(actual, expected)
        for actual, expected in zip(direction, (0, 0, -1), strict=True):
            self.assertAlmostEqual(actual, expected)

    def test_ground_round_trip_including_both_pointing_errors(self):
        context = self.context()
        context.update(known_error_deg=0.1, unknown_error_deg=0.2)
        first = resolve_telescope_frame(context, 37, 23, seed=17, telescope_id=4)
        self.assertEqual(first, resolve_telescope_frame(context, 37, 23, seed=17, telescope_id=4))
        self.assertNotEqual(
            first, resolve_telescope_frame(context, 37, 23, seed=17, telescope_id=5)
        )
        point, direction = (1, 2, 3), (0.6, 0, -0.8)
        restored = first.to_ground(*first.to_local(point, direction))
        for actual, expected in zip(
            (*restored[0], *restored[1]), (*point, *direction), strict=True
        ):
            self.assertAlmostEqual(actual, expected)
        disabled = resolve_telescope_frame(context, 37, 23, pointing_errors=False)
        self.assertEqual(disabled.known_zenith_deg, 23)
        self.assertAlmostEqual(disabled.zenith_rad, math.radians(23))

    def test_axis_height_added_once_to_ground_element_origin(self):
        context = compile_telescope_frame({
            "array_element_position_ground": dict(value=[1, 2, 4], unit="m"),
            "telescope_axis_height": dict(value=300, unit="cm"),
            "axes_offsets": dict(value=[10, 20], unit="cm"),
        })
        self.assertEqual(context["axis_origin_m"], [1, 2, 7])
        self.assertEqual(context["axes_offsets_m"], [0.1, 0.2])

    def test_csv_preserves_photon_identity_time_and_weights(self):
        with TemporaryDirectory() as directory:
            source, output = Path(directory) / "ground.csv", Path(directory) / "local.csv"
            source.write_text(
                "photon_id,x_m,y_m,z_m,dx,dy,dz,time_ns,weight\n19,10,20,31,0,0,-1,42,0.7\n"
            )
            prepare_ground_photons(source, output, self.context(), 0, 0)
            with output.open() as stream:
                row = next(csv.DictReader(stream))
            self.assertEqual((row["photon_id"], row["time_ns"], row["weight"]), ("19", "42", "0.7"))
            self.assertEqual([float(row[f"{axis}_m"]) for axis in "xyz"], [1, 0, 1])

    def test_csv_rejects_same_source_and_output_without_truncating(self):
        with TemporaryDirectory() as directory:
            source = Path(directory) / "ground.csv"
            original = "photon_id,x_m,y_m,z_m,dx,dy,dz\n19,10,20,31,0,0,-1\n"
            source.write_text(original)
            with self.assertRaisesRegex(ValueError, "must be different files"):
                prepare_ground_photons(source, source, self.context(), 0, 0)
            self.assertEqual(source.read_text(), original)

    def test_csv_pointing_seed_uses_telescope_identity(self):
        with TemporaryDirectory() as directory:
            source, output = Path(directory) / "ground.csv", Path(directory) / "local.csv"
            source.write_text(
                "telescope_id,x_m,y_m,z_m,dx,dy,dz\n"
                "4,10,20,31,0,0,-1\n4,10,20,31,0,0,-1\n5,10,20,31,0,0,-1\n"
            )
            context = self.context()
            context.update(known_error_deg=1.0, unknown_error_deg=1.0)
            prepare_ground_photons(source, output, context, 0, 0, seed=17)
            with output.open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows[0], rows[1])
            self.assertNotEqual(rows[0]["dx"], rows[2]["dx"])

    def test_single_telescope_identity_rejects_mixed_inputs(self):
        with TemporaryDirectory() as directory:
            source = Path(directory) / "ground.csv"
            source.write_text("telescope_id,x_m\n4,0\n4,1\n")
            self.assertEqual(single_telescope_identity(source), "4")
            source.write_text("telescope_id,x_m\n4,0\n5,1\n")
            with self.assertRaisesRegex(ValueError, "one consistent telescope identity"):
                single_telescope_identity(source)

    def test_transmission_reference_defaults_and_class_two_outer_power(self):
        parameter = dict(value=[0.96, 1, 0.2, 0, 0, 3])
        single = compile_telescope_transmission(parameter, 2, 10)
        dual = compile_telescope_transmission(parameter, 2, 10, dual=True)
        self.assertEqual(single["angular_scale_rad"], 0.1)
        self.assertEqual(single["power"], 2)
        self.assertEqual(single["outer_power"], 3)
        self.assertEqual(dual["outer_power"], 1)
