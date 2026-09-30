"""Tests for simulation-models production resolution."""

import json
import tempfile
import unittest
from pathlib import Path

from obdeect import model_import as IMPORTER


class TestSimulationModelsImport(unittest.TestCase):
    def make_tree(self, root: Path, *, asset_exists: bool = True):
        production = root / "productions" / "1.2.3"
        parameter = root / "model_parameters" / "TEST" / "focal_length"
        production.mkdir(parents=True)
        parameter.mkdir(parents=True)
        (production / "TEST.json").write_text(
            json.dumps({
                "model_version": "1.2.3",
                "production_table_name": "TEST",
                "parameters": {"TEST": {"focal_length": "1.0.0", "mirror_list": "1.0.0"}},
            })
        )
        (parameter / "focal_length-1.0.0.json").write_text(
            json.dumps({
                "instrument": "TEST",
                "parameter": "focal_length",
                "parameter_version": "1.0.0",
                "type": "float64",
                "unit": "cm",
                "value": 123.0,
                "file": False,
            })
        )
        mirror = root / "model_parameters" / "TEST" / "mirror_list"
        mirror.mkdir()
        (mirror / "mirror_list-1.0.0.json").write_text(
            json.dumps({
                "instrument": "TEST",
                "parameter": "mirror_list",
                "parameter_version": "1.0.0",
                "type": "string",
                "unit": None,
                "value": "mirrors.dat",
                "file": True,
            })
        )
        if asset_exists:
            assets = root / "model_parameters" / "Files"
            assets.mkdir()
            (assets / "mirrors.dat").write_text("one facet\n")

    def test_resolves_every_parameter_and_asset_with_hashes(self):
        # T-IMPORT-001: the IR retains every selected ray-tracing record and asset.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "simulation-models"
            self.make_tree(root)
            optical_model = IMPORTER.resolve_model(root, "TEST", "1.2.3")
        self.assertEqual(optical_model["format"], "obdeect.simulation-models-optical-model-ir.v1")
        self.assertEqual(optical_model["parameters"]["focal_length"]["value"], 123.0)
        self.assertEqual(
            set(optical_model["input_records"]),
            {"production_manifest", "parameter:focal_length", "parameter:mirror_list"},
        )
        self.assertIn("mirror_list", optical_model["assets"])
        self.assertEqual(len(optical_model["assets"]["mirror_list"]["sha256"]), 64)

    def test_excludes_camera_electronics_from_ray_tracing_ir(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "simulation-models"
            self.make_tree(root)
            manifest = root / "productions/1.2.3/TEST.json"
            data = json.loads(manifest.read_text())
            data["parameters"]["TEST"]["fadc_noise"] = "1.0.0"
            manifest.write_text(json.dumps(data))
            electronics = root / "model_parameters/TEST/fadc_noise"
            electronics.mkdir()
            (electronics / "fadc_noise-1.0.0.json").write_text(
                json.dumps({
                    "instrument": "TEST",
                    "parameter": "fadc_noise",
                    "parameter_version": "1.0.0",
                    "type": "float64",
                    "unit": "ct",
                    "value": 1.0,
                    "file": False,
                })
            )
            optical_model = IMPORTER.resolve_model(root, "TEST", "1.2.3")
        self.assertNotIn("fadc_noise", optical_model["parameters"])
        self.assertNotIn("parameter:fadc_noise", optical_model["input_records"])

    def test_keeps_dual_reflector_optical_parameters(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "simulation-models"
            self.make_tree(root)
            manifest = root / "productions/1.2.3/TEST.json"
            data = json.loads(manifest.read_text())
            data["parameters"]["TEST"]["secondary_mirror_parameters"] = "1.0.0"
            manifest.write_text(json.dumps(data))
            parameter = root / "model_parameters/TEST/secondary_mirror_parameters"
            parameter.mkdir()
            (parameter / "secondary_mirror_parameters-1.0.0.json").write_text(
                json.dumps({
                    "instrument": "TEST",
                    "parameter": "secondary_mirror_parameters",
                    "parameter_version": "1.0.0",
                    "type": "float64",
                    "unit": ["cm"],
                    "value": [1.0],
                    "file": False,
                })
            )
            optical_model = IMPORTER.resolve_model(root, "TEST", "1.2.3")
        self.assertIn("secondary_mirror_parameters", optical_model["parameters"])

    def test_missing_declared_asset_fails_closed(self):
        # T-IMPORT-002: an asset reference can never become an untracked path.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "simulation-models"
            self.make_tree(root, asset_exists=False)
            with self.assertRaises(IMPORTER.ImportError):
                IMPORTER.resolve_model(root, "TEST", "1.2.3")

    def test_resolves_documented_parameter_local_asset(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "simulation-models"
            self.make_tree(root, asset_exists=False)
            local_asset = root / "model_parameters/TEST/mirror_list/mirrors.dat"
            local_asset.write_text("one facet\n")
            optical_model = IMPORTER.resolve_model(root, "TEST", "1.2.3")
        self.assertEqual(
            optical_model["assets"]["mirror_list"]["path"],
            "model_parameters/TEST/mirror_list/mirrors.dat",
        )

    def test_path_components_fail_closed(self):
        # T-IMPORT-003: caller controlled selections never escape the source root.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "simulation-models"
            self.make_tree(root)
            with self.assertRaises(IMPORTER.ImportError):
                IMPORTER.resolve_model(root, "../TEST", "1.2.3")

    def test_manifest_parameter_traversal_and_missing_type(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "simulation-models"
            self.make_tree(root)
            parameter = root / "model_parameters/TEST/focal_length/focal_length-1.0.0.json"
            data = json.loads(parameter.read_text())
            del data["type"]
            parameter.write_text(json.dumps(data))
            with self.assertRaisesRegex(IMPORTER.ImportError, "incomplete"):
                IMPORTER.resolve_model(root, "TEST", "1.2.3")
            manifest = root / "productions/1.2.3/TEST.json"
            data = json.loads(manifest.read_text())
            data["parameters"]["TEST"] = {"../outside": "1.0.0"}
            manifest.write_text(json.dumps(data))
            with self.assertRaises(IMPORTER.ImportError):
                IMPORTER.resolve_model(root, "TEST", "1.2.3")

    def test_disabled_asset_and_metadata_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "simulation-models"
            self.make_tree(root)
            parameter = root / "model_parameters/TEST/mirror_list/mirror_list-1.0.0.json"
            data = json.loads(parameter.read_text())
            data.update(value=None, site="North", schema_version="0.4.0")
            parameter.write_text(json.dumps(data))
            result = IMPORTER.resolve_model(root, "TEST", "1.2.3")
            self.assertEqual(result["parameters"]["mirror_list"]["site"], "North")
            self.assertNotIn("mirror_list", result["assets"])

    def test_accepts_repository_root_containing_data_package(self):
        with tempfile.TemporaryDirectory() as directory:
            checkout = Path(directory) / "simulation-models-repository"
            root = checkout / "simulation-models"
            self.make_tree(root)
            result = IMPORTER.resolve_model(checkout, "TEST", "1.2.3")
            self.assertEqual(result["model"], "TEST")


if __name__ == "__main__":
    unittest.main()
