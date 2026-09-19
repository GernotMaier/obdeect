"""Tests for the generic simulation-models scene compiler."""

import hashlib
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from obdeect.scene_compiler import SceneCompileError, compile_scene, parse_simtel_mirror_list


class TestSceneCompiler(unittest.TestCase):
    def make_ir(self, root: Path) -> dict:
        asset = root / "model_parameters/Files/mirrors.dat"
        asset.parent.mkdir(parents=True)
        asset.write_text(
            "# x y diameter focal shape z source metadata\n0 100 120 0 1 20 # id=M01\n"
        )
        return {
            "format": "obdeect.simulation-models-ir.v1",
            "model": "GENERIC",
            "model_version": "1.0.0",
            "input_records": {},
            "assets": {
                "mirror_list": {
                    "path": "model_parameters/Files/mirrors.dat",
                    "sha256": hashlib.sha256(asset.read_bytes()).hexdigest(),
                }
            },
            "parameters": {
                "mirror_list": {"file": True, "value": "mirrors.dat"},
                "mirror_focal_length": {"file": False, "value": 1600.0, "unit": "cm"},
                "camera_body_diameter": {"file": False, "value": 200.0, "unit": "cm"},
            },
        }

    def test_compiles_tracked_mirror_list_without_dropping_deferred_fields(self):
        # T-IR-004: source units, shape, position and zero-focal fallback survive compilation.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            scene = compile_scene(self.make_ir(root), root)
        facet = scene["primary"]["facets"][0]
        self.assertEqual(scene["format"], "obdeect.compiled-scene.v1")
        self.assertEqual(facet["shape"], "hexagon_flat_y")
        self.assertEqual(facet["centre_m"], [0.0, 1.0, 0.2])
        self.assertEqual(facet["diameter_m"], 1.2)
        self.assertEqual(facet["focal_length_m"], 16.0)
        self.assertEqual(scene["report"]["deferred"], ["camera_body_diameter"])
        self.assertEqual(scene["report"]["facet_geometry_evidence"]["normal_status"], "unavailable")
        self.assertIn("No normals", scene["report"]["facet_geometry_evidence"]["interpretation"])
        self.assertEqual(len(scene["scene_sha256"]), 64)

    def test_parses_real_simtel_comment_suffix_without_turning_it_into_alignment(self):
        # T-IR-006: LST mirror-list rows carry an optional z followed by a
        # sim_telarray comment/panel ID; none of that is a facet orientation.
        facets = parse_simtel_mirror_list(
            "  1022.49 -462.00 151.00 2912.50 3 0.0 #% id=198\n"
            " -620.80 0.00 120.00 0.00 1 # no z supplied\n",
            fallback_focal_length_m=16.0,
        )
        self.assertEqual(facets[0]["centre_m"], [10.2249, -4.62, 0.0])
        self.assertEqual(facets[0]["shape"], "hexagon_flat_x")
        self.assertAlmostEqual(facets[1]["centre_m"][0], -6.208)
        self.assertEqual(facets[1]["centre_m"][1:], [0.0, 0.0])
        self.assertEqual(facets[1]["focal_length_m"], 16.0)
        self.assertNotIn("unit_normal", facets[0])
        self.assertNotIn("rotation_deg", facets[0])

    def test_zero_catalogue_fallback_allows_real_lst_rows_with_panel_focal_lengths(self):
        # T-IR-007: the LST catalogue sets mirror_focal_length to zero while
        # its mirror-list provides each panel's focal length.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            ir = self.make_ir(root)
            asset = root / "model_parameters/Files/mirrors.dat"
            asset.write_text("1022.49 -462.00 151.00 2912.50 3 0.0 #% id=198\n")
            ir["assets"]["mirror_list"]["sha256"] = hashlib.sha256(asset.read_bytes()).hexdigest()
            ir["parameters"]["mirror_focal_length"]["value"] = 0.0
            scene = compile_scene(ir, root)
        self.assertEqual(scene["primary"]["facets"][0]["focal_length_m"], 29.125)

    def test_accepts_repository_root_containing_simulation_models_data_package(self):
        # T-IR-008: importer and compiler accept the same released-checkout layout.
        with TemporaryDirectory() as directory:
            checkout = Path(directory) / "simulation-models-repository"
            data_root = checkout / "simulation-models"
            scene = compile_scene(self.make_ir(data_root), checkout)
        self.assertEqual(len(scene["primary"]["facets"]), 1)

    def test_hash_mismatch_and_malformed_records_fail_closed(self):
        # T-IR-005: asset provenance and input syntax are validated before use.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            ir = self.make_ir(root)
            (root / "model_parameters/Files/mirrors.dat").write_text("0 0 120 0 1\n")
            with self.assertRaisesRegex(SceneCompileError, "hash"):
                compile_scene(ir, root)
        with self.assertRaisesRegex(SceneCompileError, "unsupported shape"):
            parse_simtel_mirror_list("0 0 120 1600 9\n", fallback_focal_length_m=16.0)


if __name__ == "__main__":
    unittest.main()
