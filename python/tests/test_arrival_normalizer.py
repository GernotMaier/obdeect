"""Tests for conversion from the arrival contract to comparison rows."""

import csv
import os
import tempfile
import unittest
from pathlib import Path

from obdeect.arrival_normalizer import normalize_arrivals
from obdeect.result_contract import ArrivalContractError


class TestArrivalNormalizer(unittest.TestCase):
    def test_output_aliases_cannot_overwrite_frozen_arrivals(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "arrivals.csv"
            source.write_text("frozen input must remain untouched\n")
            symlink, hardlink = source.with_name("symlink.csv"), source.with_name("hardlink.csv")
            symlink.symlink_to(source)
            os.link(source, hardlink)
            for output in (source, symlink, hardlink):
                with self.subTest(output=output.name):
                    with self.assertRaisesRegex(ArrivalContractError, "alias"):
                        normalize_arrivals(source, output)
                    self.assertEqual(source.read_text(), "frozen input must remain untouched\n")

    def test_normalizes_detected_and_loss_rows(self):
        """Detected positions and terminal loss positions remain distinguishable."""
        header = (
            "contract_version,photon_id,source_kind,wavelength_nm,emission_time_ns,"
            "source_weight,throughput,status,point_count,path_length_m,"
            "incidence_primary_deg,incidence_focal_deg,"
            "x0_m,y0_m,z0_m,x1_m,y1_m,z1_m\n"
        )
        detected = "obdeect-arrival-v1,4,star,400,0,1,1,detected,2,12,3.5,1.5,1,2,3,4,5,6\n"
        lost = "obdeect-arrival-v1,5,star,400,0,1,0,missed_primary,1,7,,,8,9,10,,,\n"
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "arrivals.csv"
            normalized = Path(directory) / "normalized.csv"
            source.write_text(header + detected + lost, encoding="utf-8")

            self.assertEqual(normalize_arrivals(source, normalized), 2)
            with normalized.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(rows[0]["status"], "detected")
        self.assertEqual(rows[0]["focal_x_m"], "4.0")
        self.assertEqual(rows[0]["focal_y_m"], "5.0")
        self.assertEqual(rows[0]["incidence_primary_deg"], "3.5")
        self.assertAlmostEqual(float(rows[0]["arrival_time_ns"]), 12 / 0.299792458)
        self.assertEqual(rows[0]["incidence_secondary_deg"], "0.0")
        self.assertEqual(rows[1]["status"], "missed_primary")
        self.assertEqual(rows[1]["focal_x_m"], "8.0")
        self.assertEqual(rows[1]["focal_y_m"], "9.0")
        self.assertEqual(rows[1]["incidence_primary_deg"], "0.0")

    def test_provenance_requires_recorded_time_and_preserves_it(self):
        header = (
            "contract_version,photon_id,source_kind,wavelength_nm,emission_time_ns,"
            "source_weight,throughput,status,point_count,path_length_m,"
            "x0_m,y0_m,z0_m,x1_m,y1_m,z1_m"
        )
        row = "obdeect-arrival-v1,4,star,400,0,1,1,detected,2,12,1,2,3,4,5,6"
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "arrivals.csv"
            normalized = Path(directory) / "normalized.csv"
            source.write_text(header + "\n" + row + "\n")
            with self.assertRaisesRegex(ArrivalContractError, "recorded arrival time"):
                normalize_arrivals(source, normalized, source_sha256="b" * 64)
            self.assertFalse(normalized.exists())
            source.write_text(header + ",arrival_time_ns\n" + row + ",99\n")
            normalize_arrivals(source, normalized, source_sha256="b" * 64)
            with normalized.open(newline="") as handle:
                result = next(csv.DictReader(handle))
            self.assertEqual(float(result["arrival_time_ns"]), 99)
            self.assertNotIn("interaction_surface_ids", result)
            with self.assertRaisesRegex(ArrivalContractError, "invalid source_sha256"):
                normalize_arrivals(source, normalized, source_sha256="invalid")


if __name__ == "__main__":
    unittest.main()
