"""Tests for the generic simulation-models scene compiler."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from obdeect.model_import import resolve_model
from obdeect.scene_compiler import (
    SceneCompileError,
    compile_scene,
    parse_simtel_mirror_list,
    parse_simtel_segmentation,
)


class TestSceneCompiler(unittest.TestCase):
    def make_ir(self, root: Path, *, mirror_contents: str | None = None, focal_cm=1600.0) -> dict:
        asset = root / "model_parameters/Files/mirrors.dat"
        asset.parent.mkdir(parents=True)
        asset.write_text(
            mirror_contents
            or "# x y diameter focal shape z source metadata\n0 100 120 0 1 20 # id=M01\n"
        )
        (asset.parent / "filter.dat").write_text("300 0.8\n400 0.9\n")
        production = root / "productions/1.0.0/GENERIC.json"
        production.parent.mkdir(parents=True)
        versions = {
            name: "1.0.0"
            for name in (
                "mirror_list",
                "mirror_focal_length",
                "camera_body_diameter",
                "camera_filter",
            )
        }
        production.write_text(
            json.dumps({
                "model_version": "1.0.0",
                "production_table_name": "GENERIC",
                "parameters": {"GENERIC": versions},
            })
        )
        values = {
            "mirror_list": ("mirrors.dat", None, True),
            "mirror_focal_length": (focal_cm, "cm", False),
            "camera_body_diameter": (200.0, "cm", False),
            "camera_filter": ("filter.dat", None, True),
        }
        for name, (value, unit, is_file) in values.items():
            parameter = root / f"model_parameters/GENERIC/{name}/{name}-1.0.0.json"
            parameter.parent.mkdir(parents=True)
            parameter.write_text(
                json.dumps({
                    "instrument": "GENERIC",
                    "parameter": name,
                    "parameter_version": "1.0.0",
                    "type": "string" if is_file else "float64",
                    "unit": unit,
                    "value": value,
                    "file": is_file,
                })
            )
        return resolve_model(root, "GENERIC", "1.0.0")

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
        self.assertEqual(scene["report"]["deferred"], ["camera_body_diameter", "camera_filter"])
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

    def test_rejects_invalid_optional_mirror_height(self):
        with self.assertRaisesRegex(SceneCompileError, "invalid numeric field"):
            parse_simtel_mirror_list("0 0 120 1600 1 missing\n", fallback_focal_length_m=16.0)

    def test_zero_catalogue_fallback_allows_real_lst_rows_with_panel_focal_lengths(self):
        # T-IR-007: the LST catalogue sets mirror_focal_length to zero while
        # its mirror-list provides each panel's focal length.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            ir = self.make_ir(
                root,
                mirror_contents="1022.49 -462.00 151.00 2912.50 3 0.0 #% id=198\n",
                focal_cm=0.0,
            )
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
            with self.assertRaisesRegex(SceneCompileError, "IR assets differs"):
                compile_scene(ir, root)
        with self.assertRaisesRegex(SceneCompileError, "unsupported shape"):
            parse_simtel_mirror_list("0 0 120 1600 9\n", fallback_focal_length_m=16.0)

    def test_deferred_asset_and_parameter_records_are_verified(self):
        # T-IR-009: deferred input cannot change without invalidating provenance.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            ir = self.make_ir(root)
            (root / "model_parameters/Files/filter.dat").write_text("300 0.1\n")
            with self.assertRaisesRegex(SceneCompileError, "IR assets differs"):
                compile_scene(ir, root)
            (root / "model_parameters/Files/filter.dat").write_text("300 0.8\n400 0.9\n")
            ir["parameters"]["camera_body_diameter"]["value"] = 999.0
            with self.assertRaisesRegex(SceneCompileError, "IR parameters differs"):
                compile_scene(ir, root)

    def test_parses_explicit_hex_and_ring_footprints(self):
        # T-IR-010: ring groups expand to stable IDs; geometry remains 2D.
        segments = parse_simtel_segmentation(
            "hex 1 -85.6 0 84.6 0\nRING 2 100 200 180 -90 1.4 # cm and degrees\n"
        )
        self.assertEqual([segment["id"] for segment in segments], [0, 1, 2])
        self.assertEqual(segments[0]["centre_xy_m"], [-0.856, 0.0])
        self.assertEqual(segments[1]["inner_radius_m"], 1.0)
        self.assertEqual(segments[2]["start_deg"], 90.0)
        self.assertAlmostEqual(segments[1]["gap_m"], 0.014)
        with self.assertRaisesRegex(SceneCompileError, "unsupported type"):
            parse_simtel_segmentation("polygon 1 0 0 1 0\n")

    def test_defaults_omitted_segmentation_rotation_start_and_gap_to_zero(self):
        segments = parse_simtel_segmentation("hex 1 -85.6 0 84.6\nring 2 100 200 180\n")
        self.assertEqual(segments[0]["rotation_deg"], 0.0)
        self.assertEqual(segments[1]["start_deg"], 0.0)
        self.assertEqual(segments[2]["start_deg"], 180.0)
        self.assertEqual(segments[1]["gap_m"], 0.0)

    def test_compiles_dual_mirror_segmentation_without_invented_normals(self):
        # T-IR-011: nullable mirror_list uses explicit segmentation assets.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_ir(root)
            parameter = root / "model_parameters/GENERIC/mirror_list/mirror_list-1.0.0.json"
            data = json.loads(parameter.read_text())
            data["value"] = None
            parameter.write_text(json.dumps(data))
            production = root / "productions/1.0.0/GENERIC.json"
            data = json.loads(production.read_text())
            for name, contents in (
                ("primary_mirror_segmentation", "hex 1 0 0 80 0\n"),
                ("secondary_mirror_segmentation", "RING 2 10 20 180 0 0\n"),
            ):
                (root / f"model_parameters/Files/{name}.dat").write_text(contents)
                record = root / f"model_parameters/GENERIC/{name}/{name}-1.0.0.json"
                record.parent.mkdir(parents=True)
                record.write_text(
                    json.dumps({
                        "instrument": "GENERIC",
                        "parameter": name,
                        "parameter_version": "1.0.0",
                        "type": "string",
                        "unit": None,
                        "value": f"{name}.dat",
                        "file": True,
                    })
                )
                data["parameters"]["GENERIC"][name] = "1.0.0"
            production.write_text(json.dumps(data))
            ir = resolve_model(root, "GENERIC", "1.0.0")
            scene = compile_scene(ir, root)
        self.assertEqual(scene["primary"]["kind"], "segmented_footprints")
        self.assertEqual(len(scene["primary"]["segments"]), 1)
        self.assertEqual(len(scene["secondary"]["segments"]), 2)
        self.assertEqual(scene["report"]["facet_geometry_evidence"]["normal_status"], "unavailable")

    def test_camera_layout_count_and_nested_response_provenance(self):
        # T-IR-012: camera pixels are counted from the file, not only the record.
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.make_ir(root)
            files = root / "model_parameters/Files"
            (files / "response.dat").write_text("300 0.8\n")
            (files / "camera.dat").write_text(
                'PixType 1 0 2 0.6 2 0.7 0.1 "response.dat"\nPixel 0 1 0 0\nPixel 1 1 1 0\n'
            )
            production = root / "productions/1.0.0/GENERIC.json"
            data = json.loads(production.read_text())
            for name, value, is_file in (
                ("camera_config_file", "camera.dat", True),
                ("camera_pixels", 2, False),
            ):
                record = root / f"model_parameters/GENERIC/{name}/{name}-1.0.0.json"
                record.parent.mkdir(parents=True)
                record.write_text(
                    json.dumps({
                        "instrument": "GENERIC",
                        "parameter": name,
                        "parameter_version": "1.0.0",
                        "type": "string" if is_file else "int64",
                        "unit": None,
                        "value": value,
                        "file": is_file,
                    })
                )
                data["parameters"]["GENERIC"][name] = "1.0.0"
            production.write_text(json.dumps(data))
            ir = resolve_model(root, "GENERIC", "1.0.0")
            scene = compile_scene(ir, root)
            self.assertEqual(scene["report"]["camera_layout_evidence"]["pixel_count"], 2)
            self.assertIn("response.dat", scene["provenance"]["nested_assets"])
            count_record = root / "model_parameters/GENERIC/camera_pixels/camera_pixels-1.0.0.json"
            data = json.loads(count_record.read_text())
            data["value"] = 3
            count_record.write_text(json.dumps(data))
            with self.assertRaisesRegex(SceneCompileError, "pixel count differs"):
                compile_scene(resolve_model(root, "GENERIC", "1.0.0"), root)

    def test_rejects_nonfinite_fallback_focal_length(self):
        with self.assertRaisesRegex(SceneCompileError, "finite and positive"):
            parse_simtel_mirror_list("0 0 120 0 1\n", fallback_focal_length_m=float("nan"))


if __name__ == "__main__":
    unittest.main()
