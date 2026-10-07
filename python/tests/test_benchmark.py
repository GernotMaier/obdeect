"""Execution and fail-closed equivalence checks for generic benchmark adapters."""

import copy
import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

from obdeect.benchmark import BenchmarkError, execute, freeze, verify
from obdeect.model_import import sha256


@unittest.skipUnless(hasattr(os, "wait4") and sys.platform in {"linux", "darwin"}, "wait4 RSS")
class BenchmarkTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        model = {
            "format": "obdeect.compiled-optical-model.v1",
            "report": {
                "production_trace_ready": True,
                "native_trace_ready": True,
                "trace_blockers": [],
            },
            "trace_model": {"kind": "segmented", "primary_facets": [{}], "detector_surfaces": [{}]},
        }
        canonical = json.dumps(model, sort_keys=True, separators=(",", ":")).encode()
        model["optical_model_sha256"] = hashlib.sha256(canonical).hexdigest()
        model_path = self.directory / "model.json"
        model_path.write_text(json.dumps(model))
        source = self.directory / "source.txt"
        source.write_text("explicit resolved photons fixture\n")
        fixture = {
            "minimum_detected": 1,
            "required_surfaces": ["5"],
            "optical_model_sha256": model["optical_model_sha256"],
            "source_sha256": sha256(source),
        }
        table = self.directory / "template.csv"
        table.write_text(
            "photon_id,status,focal_x_m,focal_y_m,path_length_m,arrival_time_ns,"
            "incidence_primary_deg,incidence_secondary_deg,incidence_focal_deg,wavelength_nm,"
            "source_weight,throughput,terminal_surface_id,optical_model_sha256,source_sha256,"
            "final_dx,final_dy,final_dz\n"
            f"1,detected,0,0,6,20,0,0,0,400,2,0.8,5,{fixture['optical_model_sha256']},"
            f"{fixture['source_sha256']},0,0,1\n"
        )
        script = self.directory / "adapter.py"
        script.write_text(
            "import json, sys\nfrom pathlib import Path\n"
            "text=Path(sys.argv[1]).read_text()\n"
            "if sys.argv[-1]=='drift': text=text.replace('detected,0,0','detected,1,0')\n"
            "if sys.argv[-1]=='thread-drift' and sys.argv[4]=='2': "
            "text=text.replace('detected,0,0','detected,0.0001,0')\n"
            "Path(sys.argv[2]).write_text(text)\n"
            "Path(sys.argv[3]).write_text(json.dumps({'kernel_s':0.0001,'io_s':0.0001}))\n"
        )
        engine = {
            "argv": [
                sys.executable,
                str(script),
                str(table),
                "{output}",
                "{metrics}",
                "{threads}",
                "{block_size}",
                "same",
            ],
            "cwd": str(self.directory),
            "environment": {},
            "timing_metrics": True,
        }
        self.config = {
            "optical_model": str(model_path),
            "source": str(source),
            "fixture": fixture,
            "tolerances": {
                "focal_position_m": 0.001,
                "path_length_m": 0,
                "arrival_time_ns": 0,
                "incidence_deg": 0,
            },
            "provenance": {"software_revision": "test fixture", "physics": "matched fixture only"},
            "repeats": 10,
            "timeout_s": 5,
            "variants": [{"threads": 1, "block_size": 1}, {"threads": 2, "block_size": 7}],
            "engines": {"baseline": copy.deepcopy(engine), "candidate": copy.deepcopy(engine)},
            "inputs": [str(script), str(table)],
        }

    def test_ten_matched_repeats_with_resource_and_kernel_separation(self):
        report = execute(freeze(self.config), self.directory / "result")
        self.assertTrue(report["passed"])
        self.assertEqual(len(report["runs"]), 40)
        self.assertTrue(all(run["peak_process_rss_bytes"] > 0 for run in report["runs"]))
        for summary in report["summaries"].values():
            self.assertGreater(summary["wall_speedup"], 0)
            self.assertEqual(summary["kernel_speedup"], 1)
            self.assertIn("stdev", summary["baseline"]["io_s"])
        residuals = [
            run["matched_residuals"] for run in report["runs"] if "matched_residuals" in run
        ]
        self.assertEqual(len(residuals), 20)
        self.assertTrue(all(row["max_focal_position_residual_m"] == 0 for row in residuals))
        self.assertEqual(json.loads((self.directory / "result/benchmark.json").read_text()), report)

    def test_no_speedup_after_optical_failure(self):
        self.config["engines"]["candidate"]["argv"][-1] = "drift"
        output = self.directory / "failure"
        with self.assertRaisesRegex(BenchmarkError, "comparison failed"):
            execute(freeze(self.config), output)
        report = json.loads((output / "benchmark.json").read_text())
        self.assertFalse(report["passed"])
        self.assertNotIn("summaries", report)
        self.assertIn("first_divergence", report["error"])

    def test_thread_drift_fails_even_inside_cross_engine_tolerance(self):
        self.config["engines"]["candidate"]["argv"][-1] = "thread-drift"
        with self.assertRaisesRegex(BenchmarkError, "comparison failed"):
            execute(freeze(self.config), self.directory / "drift")

    def test_uninstrumented_runs_do_not_claim_kernel_speedup(self):
        for engine in self.config["engines"].values():
            engine["timing_metrics"] = False
        report = execute(freeze(self.config), self.directory / "wall-only")
        for summary in report["summaries"].values():
            self.assertNotIn("kernel_speedup", summary)
            self.assertNotIn("kernel_s", summary["candidate"])

    def test_bad_instrumentation_retains_failure_without_speedup(self):
        script = Path(self.config["inputs"][0])
        script.write_text(script.read_text().replace("'kernel_s':0.0001", "'kernel_s':-1"))
        output = self.directory / "bad-metrics"
        with self.assertRaisesRegex(BenchmarkError, "timing metrics"):
            execute(freeze(self.config), output)
        report = json.loads((output / "benchmark.json").read_text())
        self.assertFalse(report["passed"])
        self.assertNotIn("summaries", report)

    def test_timeout_retains_logs_and_withholds_speedup(self):
        script = Path(self.config["inputs"][0])
        script.write_text("import time\ntime.sleep(30)\n")
        self.config["timeout_s"] = 0.05
        output = self.directory / "timeout"
        with self.assertRaisesRegex(BenchmarkError, "timed_out.*True"):
            execute(freeze(self.config), output)
        report = json.loads((output / "benchmark.json").read_text())
        self.assertFalse(report["passed"])
        self.assertNotIn("summaries", report)
        self.assertTrue((output / "v000-r000-baseline/stderr.txt").exists())

    def test_provenance_and_repetition_constraints(self):
        record = freeze(self.config)
        verify(record)
        Path(self.config["source"]).write_text("changed\n")
        with self.assertRaisesRegex(BenchmarkError, "hashes"):
            verify(record)
        self.config["repeats"] = 9
        with self.assertRaisesRegex(BenchmarkError, "ten repeats"):
            freeze(self.config)

    def test_rejects_fake_reproducibility_variants(self):
        self.config["variants"] = [{"threads": 1, "block_size": 1}]
        with self.assertRaisesRegex(BenchmarkError, "multiple thread"):
            freeze(self.config)

    def test_rejects_invalid_acceptance_before_execution(self):
        for value in (None, True, -1, 0):
            with self.subTest(minimum=value):
                self.config["fixture"]["minimum_detected"] = value
                with self.assertRaisesRegex(BenchmarkError, "coverage"):
                    freeze(self.config)

    def test_instrumentation_requires_an_actual_boolean(self):
        for value in (0, 1, "true", None):
            with self.subTest(value=value):
                self.config["engines"]["candidate"]["timing_metrics"] = value
                with self.assertRaisesRegex(BenchmarkError, "timing_metrics"):
                    freeze(self.config)


if __name__ == "__main__":
    unittest.main()
