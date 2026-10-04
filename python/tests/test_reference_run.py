"""Reference-run records reject missing conventions and changed inputs."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from obdeect.reference_run import ReferenceRunError, execute, freeze, verify


class TestReferenceRun(unittest.TestCase):
    def test_freeze_and_verify(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model_root = root / "checkout"
            production = model_root / "productions" / "1.0"
            production.mkdir(parents=True)
            (production / "TEST.json").write_text(
                json.dumps({
                    "model_version": "1.0",
                    "production_table_name": "TEST",
                    "parameters": {"TEST": {}},
                })
            )
            subprocess.run(["git", "init", "-q", str(model_root)], check=True)
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(model_root),
                    "-c",
                    "user.name=Test",
                    "-c",
                    "user.email=test@example.invalid",
                    "commit",
                    "-qm",
                    "initial",
                    "--allow-empty",
                ],
                check=True,
            )
            block = root / "photons.csv"
            block.write_text("photon_id\n1\n")
            config = {
                "sim_telarray_release": "test-release",
                "hessio_decoder": "test-decoder",
                "site": "test-site",
                "production_version": "1.0",
                "telescope_variants": ["TEST"],
                "configuration_overrides": {"focus": 0},
                "seeds": {"source": 1},
                "atmosphere_extinction": {"mode": "off"},
                "photon_blocks": [str(block)],
                "coordinate_frames": {"telescope": "+z"},
                "software_revisions": {"reference": "test-revision"},
                "commands": [
                    {
                        "name": "reference",
                        "argv": [
                            sys.executable,
                            "-c",
                            "import os; print(os.environ['OPTICAL_TEST'])",
                        ],
                        "cwd": str(root),
                        "environment": {"OPTICAL_TEST": "shared-photons"},
                        "inputs": [str(block)],
                    }
                ],
            }
            reference_run = freeze(config, model_root)
            verify(reference_run, model_root)
            self.assertEqual(len(reference_run["file_sha256"]), 2)
            result = execute(reference_run, model_root, root / "success")
            self.assertTrue(result["passed"])
            self.assertEqual((root / "success/000.stdout.txt").read_text(), "shared-photons\n")
            self.assertEqual(result["commands"][0]["returncode"], 0)
            with self.assertRaises(FileExistsError):
                execute(reference_run, model_root, root / "success")
            config["commands"][0]["argv"] = [sys.executable, "-c", "import sys; sys.exit(3)"]
            failing = freeze(config, model_root)
            with self.assertRaisesRegex(ReferenceRunError, "failed; see"):
                execute(failing, model_root, root / "failure")
            failure = json.loads((root / "failure/execution.json").read_text())
            self.assertFalse(failure["passed"])
            self.assertEqual(failure["commands"][0]["returncode"], 3)
            config["commands"][0]["argv"] = [
                sys.executable,
                "-c",
                "from pathlib import Path; Path(" + repr(str(block)) + ").write_text('changed')",
            ]
            mutating = freeze(config, model_root)
            with self.assertRaisesRegex(ReferenceRunError, "failed; see"):
                execute(mutating, model_root, root / "mutating")
            mutation = json.loads((root / "mutating/execution.json").read_text())
            self.assertFalse(mutation["passed"])
            self.assertEqual(mutation["commands"][0]["returncode"], 0)
            self.assertIn("differs", mutation["commands"][0]["error"])
            # Freezing detaches the record from subsequent caller configuration edits.
            self.assertNotEqual(reference_run["configuration"]["commands"], config["commands"])
            for timeout in (0, -1, float("nan"), float("inf")):
                with (
                    self.subTest(timeout=timeout),
                    self.assertRaisesRegex(ReferenceRunError, "timeout"),
                ):
                    execute(reference_run, model_root, root / "bad-timeout", timeout)
            block.write_text("photon_id\n2\n")
            with self.assertRaisesRegex(ReferenceRunError, "differs"):
                verify(reference_run, model_root)
            block.write_text("photon_id\n1\n")
            (production / "TEST.json").write_text(
                json.dumps({
                    "model_version": "1.0",
                    "production_table_name": "TEST",
                    "parameters": {"TEST": {}},
                    "comment": "changed without a new commit",
                })
            )
            with self.assertRaisesRegex(ReferenceRunError, "differs"):
                verify(reference_run, model_root)
            del config["coordinate_frames"]
            with self.assertRaisesRegex(ReferenceRunError, "coordinate_frames"):
                freeze(config, model_root)


if __name__ == "__main__":
    unittest.main()
