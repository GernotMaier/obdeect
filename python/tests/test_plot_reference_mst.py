"""Dependency-free checks for the optional path-visualization data reader."""

import importlib.util
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / "obdeect" / "plotting.py"
SPEC = importlib.util.spec_from_file_location("plot_artificial_mst", SCRIPT)
PLOT = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(PLOT)


class TestTracePathReader(unittest.TestCase):
    def test_reads_actual_ragged_trace_vertices(self):
        # T-VIS-009: no optical calculation is performed by the plot reader.
        csv_text = (
            "photon_id,wavelength_nm,emission_time_ns,source_weight,throughput,status,point_count,path_length_m,x0_m,y0_m,z0_m,x1_m,y1_m,z1_m,x2_m,y2_m,z2_m\n"
            "4,400,1.5,2,1,detected,3,20,1,2,20,1,2,0,0.1,0.2,4.875\n"
            "5,400,0,2,0,blocked_mast,2,10,3,4,20,3,4,8,0,0,0\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "paths.csv"
            path.write_text(csv_text)
            paths = list(PLOT.read_paths(path))
        self.assertEqual(
            paths[0], ("detected", [(1.0, 2.0, 20.0), (1.0, 2.0, 0.0), (0.1, 0.2, 4.875)])
        )
        self.assertEqual(paths[1], ("blocked_mast", [(3.0, 4.0, 20.0), (3.0, 4.0, 8.0)]))

    def test_all_reference_telescope_outlines_are_declared(self):
        # T-VIS-010: every C++ catalogue name is selectable by the plot CLI.
        self.assertEqual(PLOT.TELESCOPE_NAMES, ("reference-mst", "LST", "MST", "SST", "SCT"))

    def test_extracts_weighted_detected_focal_plane_hits(self):
        csv_text = (
            "photon_id,status,point_count,source_weight,throughput,x0_m,y0_m,z0_m,x1_m,y1_m,z1_m\n"
            "1,detected,2,4,0.5,0,0,20,0.1,-0.2,4.875\n"
            "2,missed_screen,2,8,0,0,0,20,2,2,4.875\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "paths.csv"
            path.write_text(csv_text)
            hits = list(PLOT.focal_plane_hits(path))
        self.assertEqual(hits, [(0.1, -0.2, 2.0)])

    def test_old_csv_defaults_focal_plane_weight_to_one(self):
        csv_text = "status,point_count,x0_m,y0_m,z0_m\ndetected,1,0.3,0.4,4.875\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "paths.csv"
            path.write_text(csv_text)
            hits = list(PLOT.focal_plane_hits(path))
        self.assertEqual(hits, [(0.3, 0.4, 1.0)])

    def test_compiled_model_sets_full_focal_plane_extent(self):
        model = {
            "trace_model": {
                "kind": "segmented",
                "detector_surfaces": [
                    {
                        "centre_m": [0.2, -0.3, 28.0],
                        "normal": [0.0, 0.0, 1.0],
                        "tangent": [1.0, 0.0, 0.0],
                        "shape": "circle",
                        "diameter_m": 2.0,
                    },
                    {
                        "centre_m": [0.0, 0.0, 28.0],
                        "normal": [0.0, 0.0, 1.0],
                        "tangent": [1.0, 0.0, 0.0],
                        "shape": "circle",
                        "diameter_m": 1.0,
                    },
                ],
            }
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.json"
            path.write_text(json.dumps(model))
            extent = PLOT.compiled_focal_plane_extent(path)
        self.assertAlmostEqual(extent, math.hypot(0.2, -0.3) + 1.0)

    def test_compiled_mirror_polygons_follow_serialised_facet_frame(self):
        model = {
            "trace_model": {
                "kind": "segmented",
                "primary_facets": [
                    {
                        "centre_m": [1.0, 2.0, 3.0],
                        "normal": [0.0, 0.0, 1.0],
                        "tangent": [1.0, 0.0, 0.0],
                        "diameter_m": 2.0,
                        "shape": "square",
                    }
                ],
            }
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.json"
            path.write_text(json.dumps(model))
            facets = PLOT.compiled_facet_polygons(path)
        vertices, height = facets[0]
        self.assertEqual(height, 3.0)
        self.assertEqual(len(vertices), 4)
        self.assertEqual(vertices[0], (0.0, 1.0, 3.0))
        self.assertEqual(vertices[2], (2.0, 3.0, 3.0))

    def test_compiled_mirror_rejects_non_segmented_model(self):
        model = {"trace_model": {"kind": "axisymmetric"}}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.json"
            path.write_text(json.dumps(model))
            with self.assertRaisesRegex(ValueError, "segmented"):
                PLOT.compiled_facet_polygons(path)

    def test_loads_finite_segmented_scene_without_invented_hardware(self):
        model = {
            "provenance": {"model": "mini-dish", "model_version": "1"},
            "optical_model_sha256": "a" * 64,
            "trace_model": {
                "kind": "segmented",
                "primary_facets": [
                    {
                        "id": 4,
                        "shape": "square",
                        "centre_m": [0, 0, 0],
                        "normal": [0, 0, 1],
                        "tangent": [1, 0, 0],
                        "diameter_m": 2,
                    }
                ],
                "detector_surfaces": [
                    {
                        "id": 5,
                        "shape": "circle",
                        "centre_m": [0, 0, 4],
                        "normal": [0, 0, 1],
                        "tangent": [1, 0, 0],
                        "diameter_m": 1,
                    }
                ],
                "cylinder_obscurers": [
                    {
                        "id": 6,
                        "first_endpoint_m": [-1, 0, 1],
                        "second_endpoint_m": [1, 0, 1],
                        "diameter_m": 0.2,
                    }
                ],
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.json"
            path.write_text(json.dumps(model))
            scene = PLOT.load_plot_scene(path)
        self.assertEqual(scene.kind, "segmented")
        self.assertEqual(scene.model_label, "mini-dish/1")
        self.assertEqual(
            [(item.identifier, item.role) for item in scene.polygons],
            [(4, "primary"), (5, "detector")],
        )
        self.assertEqual(scene.obscurers[0].identifier, 6)
        self.assertEqual(scene.unavailable_roles, ())

    def test_axisymmetric_profile_follows_compiled_polynomial(self):
        model = {
            "trace_model": {
                "kind": "axisymmetric",
                "primary": {
                    "vertex_z_m": 2,
                    "inner_radius_m": 0,
                    "outer_radius_m": 3,
                    "radial_scale_m": 2,
                    "coefficient_m": [0, 4] + [0] * 11,
                },
                "secondary": {
                    "vertex_z_m": 4,
                    "inner_radius_m": 1,
                    "outer_radius_m": 2,
                    "radial_scale_m": 1,
                    "coefficient_m": [0] * 13,
                },
                "detector": {
                    "vertex_z_m": 1,
                    "inner_radius_m": 0,
                    "outer_radius_m": 1,
                    "radial_scale_m": 1,
                    "coefficient_m": [0] * 13,
                },
            }
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.json"
            path.write_text(json.dumps(model))
            scene = PLOT.load_plot_scene(path)
        primary = scene.axisymmetric_surfaces[0]
        profile = PLOT._axisymmetric_profile(primary, samples=3)
        self.assertEqual(profile[0], (0.0, 2.0))
        self.assertEqual(profile[-1], (3.0, 11.0))

    def test_scene_loader_rejects_duplicate_component_identifiers(self):
        surface = {
            "id": 4,
            "shape": "circle",
            "centre_m": [0, 0, 0],
            "normal": [0, 0, 1],
            "tangent": [1, 0, 0],
            "diameter_m": 1,
        }
        model = {
            "trace_model": {
                "kind": "segmented",
                "primary_facets": [surface],
                "detector_surfaces": [{**surface, "centre_m": [0, 0, 1]}],
            }
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.json"
            path.write_text(json.dumps(model))
            with self.assertRaisesRegex(ValueError, "unique"):
                PLOT.load_plot_scene(path)

    def test_scene_bounds_include_the_full_cylinder_diameter(self):
        scene = PLOT.PlotScene(
            "segmented",
            "fixture",
            None,
            (),
            (),
            (PLOT.PlotObscurer(1, (1, 2, 3), (2, 3, 4), 0.5),),
            (),
        )
        self.assertEqual(PLOT._scene_bounds(scene), ((0.75, 2.25), (1.75, 3.25), (2.75, 4.25)))

    def test_rejects_invalid_detected_focal_plane_data(self):
        csv_text = "status,point_count,x0_m,y0_m\ndetected,1,nan,0\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "paths.csv"
            path.write_text(csv_text)
            with self.assertRaisesRegex(ValueError, "non-finite"):
                list(PLOT.focal_plane_hits(path))

    def test_svg_psf_uses_weighted_hits_and_analysis(self):
        csv_text = (
            "status,point_count,source_weight,throughput,x0_m,y0_m\n"
            "detected,1,2,1,0,0\n"
            "detected,1,1,1,1,0\n"
            "missed_screen,1,1,0,5,5\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trace.csv"
            output = Path(directory) / "psf.svg"
            path.write_text(csv_text)
            env = {**os.environ, "PYTHONPATH": str(SCRIPT.parents[1])}
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "obdeect.plotting",
                    "--input",
                    str(path),
                    "--focal-plane",
                    "--telescope",
                    "LST",
                    "--bins",
                    "8",
                    "--output",
                    str(output),
                ],
                check=True,
                env=env,
            )
            svg = output.read_text()
        self.assertIn("LST focal-plane PSF", svg)
        self.assertIn("Centroid: (0.333333, 0) m", svg)
        self.assertIn("D80: 1.33333 m", svg)
        self.assertIn("throughput: 0.75", svg)

    def test_preserves_zero_recorded_throughput(self):
        csv_text = "status,point_count,source_weight,throughput,x0_m,y0_m\ndetected,1,3,0,0,0\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "paths.csv"
            path.write_text(csv_text)
            self.assertEqual(list(PLOT.focal_plane_hits(path)), [(0.0, 0.0, 0.0)])

    def test_rejects_nonfinite_path_vertex(self):
        csv_text = "status,point_count,x0_m,y0_m,z0_m\ndetected,1,nan,0,0\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "paths.csv"
            path.write_text(csv_text)
            with self.assertRaisesRegex(ValueError, "non-finite"):
                list(PLOT.read_paths(path))


if __name__ == "__main__":
    unittest.main()
