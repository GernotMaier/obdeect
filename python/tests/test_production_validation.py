"""Tests for the fail-closed production comparison gate."""

import csv
import tempfile
import unittest
from pathlib import Path

from obdeect.production_validation import (
    ProductionValidationError,
    compare,
    read_comparison_table,
    validate_optical_model,
)


class TestProductionValidation(unittest.TestCase):
    def _optical_model(self, *, sst=False):
        optical_model = {
            "report": {"native_trace_ready": True, "trace_blockers": []},
            "primary": {"facets": [{"id": 1}]},
            "detector": {"surface_id": 3},
            "materials": {"mirror": {"response": "tabulated"}},
        }
        if sst:
            optical_model["secondary"] = {"facets": [{"id": 2}]}
        return optical_model

    def _rows(self, path: Path, *, status="detected", x=0.0):
        fields = (
            "photon_id",
            "status",
            "focal_x_m",
            "focal_y_m",
            "path_length_m",
            "arrival_time_ns",
            "incidence_primary_deg",
            "incidence_secondary_deg",
            "incidence_focal_deg",
        )
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerow({
                "photon_id": 7,
                "status": status,
                "focal_x_m": x,
                "focal_y_m": 0,
                "path_length_m": 12,
                "arrival_time_ns": 42,
                "incidence_primary_deg": 10,
                "incidence_secondary_deg": 0,
                "incidence_focal_deg": 4,
            })

    def test_optical_model_requires_every_production_component(self):
        with self.assertRaisesRegex(ProductionValidationError, "production ready"):
            validate_optical_model({"report": {"trace_blockers": ["missing mirrors"]}}, "LST")
        with self.assertRaisesRegex(ProductionValidationError, "secondary"):
            validate_optical_model(self._optical_model(), "SST")
        validate_optical_model(self._optical_model(sst=True), "SST")

    def test_comparison_accepts_identical_resolved_photons(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "first.csv", Path(directory) / "second.csv"
            self._rows(first)
            self._rows(second)
            result = compare(
                read_comparison_table(first),
                read_comparison_table(second),
                {
                    "focal_position_m": 1e-4,
                    "path_length_m": 1e-4,
                    "arrival_time_ns": 1e-4,
                    "incidence_deg": 1e-4,
                },
            )
        self.assertEqual(result["photon_count"], 1)
        self.assertIsNone(result["first_divergence"])

    def test_comparison_rejects_first_optical_divergence(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "first.csv", Path(directory) / "second.csv"
            self._rows(first, x=0.0)
            self._rows(second, x=0.1)
            with self.assertRaisesRegex(ProductionValidationError, '"photon_id": 7'):
                compare(
                    read_comparison_table(first),
                    read_comparison_table(second),
                    {
                        "focal_position_m": 1e-3,
                        "path_length_m": 1e-3,
                        "arrival_time_ns": 1e-3,
                        "incidence_deg": 1e-3,
                    },
                )

    def test_comparison_only_checks_loss_categories(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "first.csv", Path(directory) / "second.csv"
            self._rows(first, status="missed_primary", x=0.0)
            self._rows(second, status="missed_primary", x=100.0)
            result = compare(
                read_comparison_table(first),
                read_comparison_table(second),
                {
                    "focal_position_m": 1e-3,
                    "path_length_m": 1e-3,
                    "arrival_time_ns": 1e-3,
                    "incidence_deg": 1e-3,
                },
            )
        self.assertIsNone(result["first_divergence"])
        self.assertEqual(result["max_focal_position_residual_m"], 0.0)

    def test_comparison_rejects_arrival_time_divergence(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "first.csv", Path(directory) / "second.csv"
            self._rows(first)
            self._rows(second)
            with second.open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle))
            row["arrival_time_ns"] = "42.1"
            with second.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=row.keys())
                writer.writeheader()
                writer.writerow(row)
            with self.assertRaisesRegex(ProductionValidationError, '"arrival_time_ns"'):
                compare(
                    read_comparison_table(first),
                    read_comparison_table(second),
                    {
                        "focal_position_m": 1e-3,
                        "path_length_m": 1e-3,
                        "arrival_time_ns": 1e-3,
                        "incidence_deg": 1e-3,
                    },
                )

    def test_comparison_rejects_unknown_status(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "arrivals.csv"
            self._rows(path, status="typo_status")
            with self.assertRaisesRegex(ProductionValidationError, "invalid status"):
                read_comparison_table(path)

    def test_comparison_rejects_malformed_tolerances(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "first.csv", Path(directory) / "second.csv"
            self._rows(first)
            self._rows(second)
            for value in ("bad", None, [], {}, True):
                tolerances = {
                    "focal_position_m": 1e-4,
                    "path_length_m": 1e-4,
                    "arrival_time_ns": 1e-4,
                    "incidence_deg": value,
                }
                with (
                    self.subTest(value=value),
                    self.assertRaisesRegex(ProductionValidationError, "tolerances must define"),
                ):
                    compare(read_comparison_table(first), read_comparison_table(second), tolerances)


if __name__ == "__main__":
    unittest.main()
