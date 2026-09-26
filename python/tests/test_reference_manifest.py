"""Reference manifests reject missing conventions and changed inputs."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from obdeect.reference_manifest import ManifestError, freeze, verify


class TestReferenceManifest(unittest.TestCase):
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
                        "argv": [sys.executable],
                        "cwd": str(root),
                        "environment": {},
                        "inputs": [str(block)],
                    }
                ],
            }
            manifest = freeze(config, model_root)
            verify(manifest, model_root)
            self.assertEqual(len(manifest["file_sha256"]), 2)
            block.write_text("photon_id\n2\n")
            with self.assertRaisesRegex(ManifestError, "differs"):
                verify(manifest, model_root)
            block.write_text("photon_id\n1\n")
            (production / "TEST.json").write_text(
                json.dumps({
                    "model_version": "1.0",
                    "production_table_name": "TEST",
                    "parameters": {"TEST": {}},
                    "comment": "changed without a new commit",
                })
            )
            with self.assertRaisesRegex(ManifestError, "differs"):
                verify(manifest, model_root)
            del config["coordinate_frames"]
            with self.assertRaisesRegex(ManifestError, "coordinate_frames"):
                freeze(config, model_root)


if __name__ == "__main__":
    unittest.main()
