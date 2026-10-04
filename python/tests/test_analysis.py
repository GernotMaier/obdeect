"""Analytic contracts for focal-plane PSF analysis."""

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).parents[1]))

from obdeect.analysis import PsfResult, _scan, analyse_trace_csv, main, write_scan_csv


class TestPsfAnalysis(unittest.TestCase):
    def _trace(self, directory: Path) -> Path:
        path = directory / "trace.csv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "status",
                    "point_count",
                    "source_weight",
                    "throughput",
                    "x0_m",
                    "y0_m",
                    "x1_m",
                    "y1_m",
                ),
            )
            writer.writeheader()
            writer.writerows((
                {
                    "status": "detected",
                    "point_count": 2,
                    "source_weight": 0.2,
                    "throughput": 1,
                    "x0_m": 0,
                    "y0_m": 0,
                    "x1_m": 0,
                    "y1_m": 0,
                },
                {
                    "status": "detected",
                    "point_count": 2,
                    "source_weight": 1.6,
                    "throughput": 0.5,
                    "x0_m": 0,
                    "y0_m": 0,
                    "x1_m": 1,
                    "y1_m": 0,
                },
                {
                    "status": "blocked_camera",
                    "point_count": 1,
                    "source_weight": 1,
                    "throughput": 0,
                    "x0_m": 0,
                    "y0_m": 0,
                    "x1_m": "",
                    "y1_m": "",
                },
            ))
        return path

    def test_uses_photon_optical_response_and_unbinned_d80(self):
        # T-VIS-015: D80 uses each detected photon's optical response, not a
        # photon-count quantile or a histogram-dependent approximation.
        with tempfile.TemporaryDirectory() as directory:
            result = analyse_trace_csv(self._trace(Path(directory)))
        self.assertEqual(result.detected_count, 2)
        self.assertAlmostEqual(result.input_weight, 2.8)
        self.assertAlmostEqual(result.detected_weight, 1.0)
        self.assertAlmostEqual(result.optical_throughput, 1.0 / 2.8)
        self.assertAlmostEqual(result.centroid_x_m, 0.8)
        self.assertAlmostEqual(result.r80_m, 0.2)
        self.assertAlmostEqual(result.d80_m, 0.4)
        self.assertEqual(result.terminal_count, {"blocked_camera": 1, "detected": 2})
        self.assertEqual(result.terminal_input_weight, {"blocked_camera": 1.0, "detected": 1.8})

    def test_arrival_contract_preserves_psf_and_validates_loss_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            legacy = self._trace(root)
            expected = analyse_trace_csv(legacy)
            with legacy.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            for photon_id, row in enumerate(rows):
                row.update(
                    contract_version="obdeect-arrival-v1",
                    photon_id=photon_id,
                    source_kind="replay",
                    wavelength_nm=400,
                    emission_time_ns=0,
                    path_length_m=1,
                    z0_m=0,
                    z1_m=1 if row["status"] == "detected" else "",
                )
            modern = root / "arrivals.csv"

            def write_arrivals():
                with modern.open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=rows[0])
                    writer.writeheader()
                    writer.writerows(rows)

            write_arrivals()
            self.assertEqual(analyse_trace_csv(modern), expected)
            # Contract validation covers terminal losses too, after the detections.
            rows[-1]["z0_m"] = "nan"
            write_arrivals()
            with self.assertRaisesRegex(ValueError, "non-finite z0_m"):
                analyse_trace_csv(modern)
            rows[-1]["z0_m"] = 0
            rows[-1]["photon_id"] = 0
            write_arrivals()
            with self.assertRaisesRegex(ValueError, "duplicate photon_id"):
                analyse_trace_csv(modern)

    def test_derivation_writes_field_angle_metrics_without_retracing(self):
        # T-VIS-016: an off-axis metric carries its configured field direction
        # and the analysis consumes the recorded detector intersection only.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            trace = self._trace(root)
            output = root / "metrics.json"
            scan_csv = root / "scan.csv"
            with mock.patch(
                "sys.argv",
                [
                    "obdeect-psf",
                    "derive",
                    "--input",
                    str(trace),
                    "--field-x-deg",
                    "0.5",
                    "--output",
                    str(output),
                    "--scan-csv",
                    str(scan_csv),
                ],
            ):
                self.assertEqual(main(), 0)
            payload = json.loads(output.read_text())
            self.assertEqual(payload["field_x_deg"], 0.5)
            self.assertAlmostEqual(payload["d80_m"], 0.4)
            self.assertIn("field_x_deg", scan_csv.read_text().splitlines()[0])

    def test_scan_csv_has_one_row_per_angle(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = analyse_trace_csv(self._trace(root))
            output = root / "scan.csv"
            write_scan_csv([(0.0, result), (1.0, result)], output)
            self.assertEqual(len(output.read_text().splitlines()), 3)

    def test_scan_trace_names_are_unique_for_close_angles(self):
        result = PsfResult(1.0, 1.0, 1.0, 1, 0.0, 0.0, 0.0, 0.0, {}, {})
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = type(
                "ScanArgs",
                (),
                {
                    "output_dir": root,
                    "field_x_deg": [0.0000001, 0.0000002],
                    "field_y_deg": 0.0,
                    "photons": 4,
                    "executable": root / "obdeect-demo-mst",
                    "extra_argument": [],
                    "plot": None,
                },
            )()
            with (
                mock.patch("obdeect.analysis.subprocess.run") as run,
                mock.patch("obdeect.analysis.analyse_trace_csv", return_value=result),
            ):
                _scan(args)
            output_paths = [
                call.args[0][call.args[0].index("--output") + 1] for call in run.call_args_list
            ]
        self.assertEqual(len(output_paths), 2)
        self.assertNotEqual(output_paths[0], output_paths[1])

    def test_missing_final_vertex_is_reported_as_value_error(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "missing_vertex.csv"
            path.write_text(
                "status,point_count,source_weight,throughput,x0_m,y0_m,x1_m,y1_m\n"
                "detected,2,1,1,0,0\n"
            )
            with self.assertRaisesRegex(ValueError, "missing final focal-plane vertex"):
                analyse_trace_csv(path)

    def test_no_detected_photons_is_reported_as_a_trace_outcome(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "misses.csv"
            path.write_text(
                "status,point_count,source_weight,throughput,x0_m,y0_m\nmissed_primary,1,1,0,0,0\n"
            )
            with self.assertRaisesRegex(ValueError, "no detected focal-plane photons"):
                analyse_trace_csv(path)

    def test_missing_or_invalid_weights_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self._trace(Path(directory))
            contents = path.read_text()
            for corrupted, message in (
                (contents.replace("detected,2,0.2,1,", "detected,2,,1,"), "missing source_weight"),
                (
                    contents.replace("blocked_camera,1,1,0,", "unknown,1,1,0,"),
                    "invalid terminal status",
                ),
                (
                    contents.replace("detected,2,0.2,1,", "detected,2,0.2,1.1,"),
                    "throughput must not exceed 1",
                ),
                (
                    contents.replace("blocked_camera,1,1,0,", "blocked_camera,1,1,0.5,"),
                    "lost photon has nonzero throughput",
                ),
            ):
                path.write_text(corrupted)
                with self.assertRaisesRegex(ValueError, message):
                    analyse_trace_csv(path)


if __name__ == "__main__":
    unittest.main()
