"""Tests for the fail-closed production comparison gate."""

import csv
import hashlib
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import pytest
from obdeect.production_validation import (
    ComparisonRow,
    ProductionValidationError,
    compare,
    main,
    read_comparison_table,
    validate_optical_model,
)


class TestProductionValidation(unittest.TestCase):
    def _optical_model(self, *, sst=False):
        optical_model = {
            "format": "obdeect.compiled-optical-model.v1",
            "report": {
                "native_trace_ready": True,
                "production_trace_ready": True,
                "trace_blockers": [],
            },
            "trace_model": {
                "kind": "segmented",
                "primary_facets": [{"id": 1}],
                "detector_surfaces": [{"id": 3}],
            },
        }
        if sst:
            optical_model["trace_model"] = {
                "kind": "axisymmetric",
                "primary": {},
                "secondary": {},
                "detector": {},
            }
        canonical = json.dumps(
            optical_model, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        ).encode()
        optical_model["optical_model_sha256"] = hashlib.sha256(canonical).hexdigest()
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

    def test_optical_model_accepts_embedded_trace_schema(self):
        validate_optical_model(self._optical_model(), "LST")
        validate_optical_model(self._optical_model(sst=True), "SST")

    def test_optical_model_rejects_tampered_hash_and_nominal_readiness(self):
        model = self._optical_model()
        model["trace_model"]["primary_facets"][0]["id"] = 99
        with self.assertRaisesRegex(ProductionValidationError, "content hash"):
            validate_optical_model(model, "LST")
        model = self._optical_model()
        del model["report"]["production_trace_ready"]
        model["report"]["trace_ready"] = True
        with self.assertRaisesRegex(ProductionValidationError, "not production ready"):
            validate_optical_model(model, "LST")

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

    def test_comparison_rejects_undeclared_all_loss(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "first.csv", Path(directory) / "second.csv"
            self._rows(first, status="missed_primary", x=0.0)
            self._rows(second, status="missed_primary", x=100.0)
            with self.assertRaisesRegex(ProductionValidationError, "detection coverage"):
                compare(read_comparison_table(first), read_comparison_table(second), TOLERANCES)

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


TOLERANCES = dict(
    focal_position_m=1e-3, path_length_m=1e-3, arrival_time_ns=1e-3, incidence_deg=1e-3
)
FIXTURE = dict(
    minimum_detected=1,
    required_surfaces=["3"],
    optical_model_sha256="a" * 64,
    source_sha256="b" * 64,
)
PHOTON = ComparisonRow(
    7, "detected", 0, 0, 12, 42, 10, 0, 4, 400, 1, 0.8, "3", "a" * 64, "b" * 64, (0, 0, 1)
)


@pytest.mark.parametrize(
    "field,value",
    [
        ("source_weight", 2),
        ("throughput", 0.7),
        ("wavelength_nm", 401),
        ("terminal_surface_id", "4"),
        ("source_sha256", "c" * 64),
        ("optical_model_sha256", "c" * 64),
        ("final_direction", (1, 0, 0)),
    ],
)
def test_declared_fixture_rejects_wrong_observables(field, value):
    with pytest.raises(ProductionValidationError):
        compare({7: PHOTON}, {7: replace(PHOTON, **{field: value})}, TOLERANCES, FIXTURE)


def test_declared_negative_fixture_ignores_arbitrary_loss_endpoint_but_checks_weight():
    fixture = FIXTURE | dict(minimum_detected=0, allow_all_loss=True, required_surfaces=[])
    lost = replace(PHOTON, status="missed_primary", throughput=0)
    assert (
        compare(
            {7: lost}, {7: replace(lost, path_length_m=500, focal_x_m=100)}, TOLERANCES, fixture
        )["detected_count"]
        == 0
    )
    with pytest.raises(ProductionValidationError):
        compare({7: lost}, {7: replace(lost, source_weight=2)}, TOLERANCES, fixture)


def test_comparator_persists_failed_report(tmp_path, monkeypatch):
    output = tmp_path / "nested" / "report.json"
    missing = tmp_path / "missing.json"
    monkeypatch.setattr(
        "sys.argv",
        [
            "production-validation",
            "--optical-model",
            str(missing),
            "--telescope-family",
            "LST",
            "--obdeect-arrivals",
            str(missing),
            "--simtel-arrivals",
            str(missing),
            "--tolerances",
            str(missing),
            "--fixture",
            str(missing),
            "--output",
            str(output),
        ],
    )
    with pytest.raises(SystemExit):
        main()
    report = json.loads(output.read_text())
    assert report["passed"] is False
    assert "missing.json" in report["error"]


def test_matching_invalid_status_cannot_pass():
    invalid = replace(PHOTON, status="unknown")
    with pytest.raises(ProductionValidationError, match="invalid terminal status"):
        compare({7: invalid}, {7: invalid}, TOLERANCES)


@pytest.mark.parametrize(
    "fixture",
    [
        FIXTURE | dict(minimum_detected=2),
        FIXTURE | dict(required_surfaces=["4"]),
        FIXTURE | dict(minimum_detected=0),
    ],
)
def test_declared_coverage_is_required(fixture):
    with pytest.raises(ProductionValidationError):
        compare({7: PHOTON}, {7: PHOTON}, TOLERANCES, fixture)


def test_first_different_interaction_surface_is_reported():
    first = replace(PHOTON, interaction_surface_ids=(3, 4, 5))
    second = replace(PHOTON, interaction_surface_ids=(3, 9, 5))
    with pytest.raises(ProductionValidationError, match='"interaction_index": 1'):
        compare({7: first}, {7: second}, TOLERANCES, FIXTURE)


@pytest.mark.parametrize(
    "response,terminal",
    [(0.1, 0), (0.1, 0.1), (float("nan"), 0), (None, 0.2), (-0.1, 0.3)],
)
def test_direct_adapter_cannot_bypass_optical_loss_ledger(response, terminal):
    invalid = replace(PHOTON, response_loss_fraction=response, terminal_loss_fraction=terminal)
    with pytest.raises(ProductionValidationError, match="loss"):
        compare({7: invalid}, {7: invalid}, TOLERANCES, FIXTURE)


def test_direct_adapter_accepts_closed_optical_loss_ledger():
    photon = replace(PHOTON, response_loss_fraction=0.2, terminal_loss_fraction=0)
    summary = compare({7: photon}, {7: photon}, TOLERANCES, FIXTURE)
    assert summary["response_loss_weight"] == 0.2
    assert summary["terminal_loss_weight"] == 0


def test_comparison_table_rejects_detected_terminal_loss(tmp_path):
    path = tmp_path / "invalid.csv"
    fields = [*PHOTON.__dataclass_fields__]
    fields.remove("final_direction")
    fields.remove("interaction_surface_ids")
    values = {field: getattr(PHOTON, field) for field in fields}
    values.update(response_loss_fraction=0.1, terminal_loss_fraction=0.1)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(values)
    with pytest.raises(ProductionValidationError, match="terminal loss"):
        read_comparison_table(path)


@pytest.mark.parametrize("identity", [8, -7, True, (1, 7), (1, 2, 3, 4, -5, 7)])
def test_direct_adapter_requires_correct_identity_binding(identity):
    with pytest.raises(ProductionValidationError, match="identity"):
        compare({identity: PHOTON}, {identity: PHOTON}, TOLERANCES, FIXTURE)


def test_direct_adapter_rejects_mixed_identity_context():
    rows = {7: PHOTON, (0, 0, 0, 0, 0, 7): PHOTON}
    with pytest.raises(ProductionValidationError, match="identity"):
        compare(rows, rows, TOLERANCES, FIXTURE)
