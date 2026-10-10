"""Tests for simulation-models production resolution."""

import json
import tempfile
import unittest
from pathlib import Path

from obdeect import model_import as IMPORTER


class TestSimulationModelsImport(unittest.TestCase):
    def test_declared_site_requires_environment_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_tree(root)
            parameter = root / "model_parameters/TEST/focal_length/focal_length-1.0.0.json"
            record = json.loads(parameter.read_text())
            record["site"] = "North"
            parameter.write_text(json.dumps(record))
            with self.assertRaisesRegex(
                IMPORTER.ImportError, "site environment manifest is missing"
            ):
                IMPORTER.resolve_model(root, "TEST", "1.2.3")

    def test_explicit_site_environment_is_resolved_with_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_tree(root)
            focal = root / "model_parameters/TEST/focal_length/focal_length-1.0.0.json"
            record = json.loads(focal.read_text())
            record["site"] = "Example"
            focal.write_text(json.dumps(record))
            manifest = root / "productions/1.2.3/OBS-Example.json"
            names = ("atmospheric_profile", "corsika_observation_level")
            manifest.write_text(
                json.dumps({
                    "model_version": "1.2.3",
                    "parameters": {"OBS-Example": dict.fromkeys(names, "1.0.0")},
                })
            )
            for name, value, unit, file in (
                ("atmospheric_profile", "profile.ecsv", None, True),
                ("corsika_observation_level", 1000, "m", False),
            ):
                path = root / "model_parameters/OBS-Example" / name / f"{name}-1.0.0.json"
                path.parent.mkdir(parents=True)
                path.write_text(
                    json.dumps(
                        dict(
                            instrument="OBS-Example",
                            parameter=name,
                            parameter_version="1.0.0",
                            value=value,
                            unit=unit,
                            file=file,
                        )
                    )
                )
            profile = root / "model_parameters/OBS-Example/atmospheric_profile/profile.ecsv"
            profile.write_text("altitude refractive_index\n0 0.0002\n1 0.0002\n2 0.0002\n")
            model = IMPORTER.resolve_model(root, "TEST", "1.2.3")
            self.assertEqual(model["environment"]["propagation_group_index"], 1.0002)
            self.assertEqual(model["environment"]["observation_level_m"], 1000)
            self.assertEqual(model["parameters"]["focal_length"]["site"], "Example")
            self.assertEqual(
                model["input_records"]["environment_atmospheric_profile"]["sha256"],
                IMPORTER.sha256(profile),
            )
            profile_parameter = (
                root
                / "model_parameters/OBS-Example/atmospheric_profile"
                / "atmospheric_profile-1.0.0.json"
            )
            profile_record = json.loads(profile_parameter.read_text())
            for invalid_value in (None, "missing.ecsv", "../outside.ecsv", "bad\0.ecsv"):
                profile_record["value"] = invalid_value
                profile_parameter.write_text(json.dumps(profile_record))
                with self.assertRaisesRegex(IMPORTER.ImportError, "site atmospheric profile"):
                    IMPORTER.resolve_model(root, "TEST", "1.2.3")
            profile_record["value"] = "profile.ecsv"
            profile_parameter.write_text(json.dumps(profile_record))
            site_parameter = (
                root
                / "model_parameters/OBS-Example/corsika_observation_level"
                / "corsika_observation_level-1.0.0.json"
            )
            site_record = json.loads(site_parameter.read_text())
            site_record["parameter_version"] = "9.9.9"
            site_parameter.write_text(json.dumps(site_record))
            with self.assertRaisesRegex(
                IMPORTER.ImportError, "site environment parameter identity"
            ):
                IMPORTER.resolve_model(root, "TEST", "1.2.3")
            site_record["parameter_version"] = "1.0.0"
            site_parameter.write_text(json.dumps(site_record))
            profile.write_text("altitude refractive_index\n0\n1 0.0002\n2 0.0002\n")
            with self.assertRaisesRegex(IMPORTER.ImportError, "invalid site environment"):
                IMPORTER.resolve_model(root, "TEST", "1.2.3")
            profile.write_text("altitude refractive_index\n0 nan\n1 0.0002\n2 0.0002\n")
            with self.assertRaisesRegex(IMPORTER.ImportError, "invalid site environment"):
                IMPORTER.resolve_model(root, "TEST", "1.2.3")

    def test_instance_inherits_design_and_overrides_with_complete_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_tree(root)
            instance = root / "productions/1.2.3/INSTANCE.json"
            instance.write_text(
                json.dumps({
                    "model_version": "1.2.3",
                    "production_table_name": "INSTANCE",
                    "design_model": {"INSTANCE": "TEST"},
                    "parameters": {"INSTANCE": {"focal_length": "1.0.0"}},
                })
            )
            path = root / "model_parameters/INSTANCE/focal_length/focal_length-1.0.0.json"
            path.parent.mkdir(parents=True)
            path.write_text(
                json.dumps({
                    "instrument": "INSTANCE",
                    "parameter": "focal_length",
                    "parameter_version": "1.0.0",
                    "type": "float64",
                    "unit": "cm",
                    "value": 456,
                    "file": False,
                })
            )
            model = IMPORTER.resolve_model(root, "INSTANCE", "1.2.3")
            self.assertEqual(model["parameters"]["focal_length"]["value"], 456)
            self.assertEqual(model["parameters"]["mirror_list"]["instrument"], "TEST")
            self.assertIn("design_manifest:TEST", model["input_records"])
            self.assertIn("source_parameter:TEST:focal_length", model["input_records"])
            self.assertIn("source_parameter:INSTANCE:focal_length", model["input_records"])
            design = root / "productions/1.2.3/TEST.json"
            data = json.loads(design.read_text())
            data["design_model"] = {"TEST": "INSTANCE"}
            design.write_text(json.dumps(data))
            with self.assertRaisesRegex(IMPORTER.ImportError, "cyclic"):
                IMPORTER.resolve_model(root, "INSTANCE", "1.2.3")

    def test_table_filename_options_preserve_asset_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_tree(root)
            production = root / "productions/1.2.3/TEST.json"
            manifest = json.loads(production.read_text())
            manifest["parameters"]["TEST"]["mirror_reflectivity"] = "1.0.0"
            production.write_text(json.dumps(manifest))
            path = root / "model_parameters/TEST/mirror_list/mirror_list-1.0.0.json"
            parameter = json.loads(path.read_text())
            parameter["parameter"] = "mirror_reflectivity"
            parameter["value"] = "mirrors.dat#rpol:clip,scheme=0"
            path = root / "model_parameters/TEST/mirror_reflectivity/mirror_reflectivity-1.0.0.json"
            path.parent.mkdir()
            path.write_text(json.dumps(parameter))
            model = IMPORTER.resolve_model(root, "TEST", "1.2.3")
            self.assertEqual(
                model["parameters"]["mirror_reflectivity"]["value"], parameter["value"]
            )
            self.assertTrue(model["assets"]["mirror_reflectivity"]["path"].endswith("mirrors.dat"))

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

    def test_unknown_required_optical_parameter_is_not_silently_discarded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_tree(root)
            manifest = root / "productions/1.2.3/TEST.json"
            data = json.loads(manifest.read_text())
            data["parameters"]["TEST"]["new_window_curvature"] = "1.0.0"
            manifest.write_text(json.dumps(data))
            folder = root / "model_parameters/TEST/new_window_curvature"
            folder.mkdir()
            parameter = folder / "new_window_curvature-1.0.0.json"
            parameter.write_text(
                json.dumps({
                    "instrument": "TEST",
                    "parameter": "new_window_curvature",
                    "parameter_version": "1.0.0",
                    "type": "float64",
                    "unit": "m",
                    "value": 1.0,
                    "file": False,
                    "required_for_trace": True,
                })
            )
            with self.assertRaisesRegex(IMPORTER.ImportError, "required optical parameter"):
                IMPORTER.resolve_model(root, "TEST", "1.2.3")
            data = json.loads(parameter.read_text())
            data["required_for_trace"] = False
            parameter.write_text(json.dumps(data))
            coverage = IMPORTER.resolve_model(root, "TEST", "1.2.3")["source_parameter_coverage"]
            self.assertEqual(coverage["new_window_curvature"]["disposition"], "unsupported")
            self.assertEqual(
                coverage["new_window_curvature"]["record"]["sha256"], IMPORTER.sha256(parameter)
            )

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
            data.update(value=None, schema_version="0.4.0")
            parameter.write_text(json.dumps(data))
            result = IMPORTER.resolve_model(root, "TEST", "1.2.3")
            self.assertEqual(result["parameters"]["mirror_list"]["schema_version"], "0.4.0")
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
