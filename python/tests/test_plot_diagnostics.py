"""Model geometry, selection and loss observables without optical inference."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from obdeect import plotting as plot
from obdeect.camera_surfaces import compile_camera_surfaces


def model_fixture():
    surface = {
        "id": 4,
        "shape": "square",
        "centre_m": [0, 0, 0],
        "normal": [0, 0, 1],
        "tangent": [1, 0, 0],
        "diameter_m": 2,
    }
    return {
        "provenance": {"model": "finite-fixture", "model_version": "1"},
        "optical_model_sha256": "a" * 64,
        "trace_model": {
            "kind": "segmented",
            "primary_facets": [surface],
            "detector_surfaces": [{**surface, "id": 5, "centre_m": [0, 0, 4], "diameter_m": 0.5}],
            "cylinder_obscurers": [
                {
                    "id": 6,
                    "first_endpoint_m": [-1, 0, 2],
                    "second_endpoint_m": [1, 0, 2],
                    "diameter_m": 0.2,
                }
            ],
        },
        "plot_geometry": {
            "schema_version": 1,
            "frame": {"origin": "primary_vertex", "axes": "+x,+y,+z", "unit": "m"},
            "components": [
                {"id": 4, "role": "primary", "source": "trace_model.primary_facets"},
                {"id": 5, "role": "detector", "source": "trace_model.detector_surfaces"},
                {"id": 6, "role": "opaque_cylinder", "source": "trace_model.cylinder_obscurers"},
            ],
            "unavailable_roles": ["window"],
        },
    }


class DiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.model = self.directory / "model.json"
        self.model.write_text(json.dumps(model_fixture()))
        self.csv = self.directory / "trace.csv"
        self.csv.write_text(
            "photon_id,source_kind,status,source_weight,throughput,point_count,terminal_surface_id,"
            "x0_m,y0_m,z0_m,x1_m,y1_m,z1_m\n"
            "1,star,detected,3,1,2,5,0,0,10,0.1,0,4\n"
            "2,laser,blocked_obscurer,1,0,2,6,0,0,10,0,0,2\n"
            "3,star,missed_primary,2,0,1,4294967295,0.5,0.5,10,0,0,0\n"
        )

    def test_index_references_only_available_components(self):
        optical_model = plot.load_plot_optical_model(self.model)
        self.assertFalse(optical_model.legacy_geometry)
        self.assertEqual(optical_model.frame_origin, "primary_vertex")
        self.assertIn("+z optical axis", plot._optical_model_footer(optical_model))
        self.assertIn("window", plot._optical_model_footer(optical_model))
        broken = model_fixture()
        broken["plot_geometry"]["components"][0]["id"] = 999
        self.model.write_text(json.dumps(broken))
        with self.assertRaisesRegex(ValueError, "absent component"):
            plot.load_plot_optical_model(self.model)

    def test_index_schema_and_units_rejected(self):
        for field, value in (("schema_version", 2), ("frame", {"unit": "cm"})):
            with self.subTest(field=field):
                broken = model_fixture()
                broken["plot_geometry"][field] = value
                self.model.write_text(json.dumps(broken))
                with self.assertRaises(ValueError):
                    plot.load_plot_optical_model(self.model)

    def test_dual_camera_draws_actual_compiled_entrances(self):
        from obdeect.optical_model_compiler import build_plot_geometry

        surface = {
            "vertex_z_m": 0,
            "inner_radius_m": 0,
            "outer_radius_m": 100,
            "radial_scale_m": 1,
            "coefficient_m": [0] * 13,
        }
        plane = {
            "id": 10,
            "source_pixel_id": 42,
            "shape": "square",
            "centre_m": [1, 2, 3],
            "normal": [0, 0, 1],
            "tangent": [0, 1, 0],
            "diameter_m": 0.4,
        }
        fixture = model_fixture()
        trace = {
            "kind": "axisymmetric",
            "primary": surface,
            "secondary": surface,
            "detector": surface,
            "detector_surfaces": [plane],
        }
        fixture["trace_model"] = trace
        fixture["camera"] = {"entrance_surfaces": [plane]}
        fixture["plot_geometry"] = build_plot_geometry(trace, {"trace_blockers": []})
        self.model.write_text(json.dumps(fixture))
        compiled = plot.load_plot_optical_model(self.model)
        self.assertEqual(
            [surface.role for surface in compiled.axisymmetric_surfaces], ["primary", "secondary"]
        )
        self.assertEqual(compiled.polygons[0].identifier, 10)
        self.assertEqual(compiled.camera_pixels[0].identifier, 42)
        self.assertEqual(
            compiled.camera_pixels[0].vertices_xy_m,
            tuple((p[0], p[1]) for p in compiled.polygons[0].vertices_m),
        )
        self.assertLess(plot.compiled_focal_plane_extent(self.model), 3)

    def test_surface_ids_reject_fractional_boolean_and_reserved_values(self):
        for identifier in (4.5, True, -1, 4294967295):
            with self.subTest(identifier=identifier):
                fixture = model_fixture()
                fixture.pop("plot_geometry")
                fixture["trace_model"]["primary_facets"][0]["id"] = identifier
                self.model.write_text(json.dumps(fixture))
                with self.assertRaisesRegex(ValueError, "component id"):
                    plot.load_plot_optical_model(self.model)

    def test_nonorthogonal_frame_rejected(self):
        broken = model_fixture()
        broken["trace_model"]["primary_facets"][0]["tangent"] = [1, 0, 0.1]
        self.model.write_text(json.dumps(broken))
        with self.assertRaisesRegex(ValueError, "perpendicular"):
            plot.load_plot_optical_model(self.model)

    def test_loss_is_input_weight_fraction_empty_bins_masked(self):
        fractions = plot.weighted_loss_bins(list(plot.entrance_samples(self.csv)), 2, 1)
        self.assertEqual(fractions[1][1], 0.5)
        self.assertIsNone(fractions[0][0])
        # Zero detector response still means the photon reached the detector.
        fractions = plot.weighted_loss_bins([(0, 0, 2, "detected", None)], 1, 1)
        self.assertEqual(fractions, [[0.0]])

    def test_component_and_status_filters_use_recorded_fields(self):
        chosen = plot.select_trace_rows(self.csv, statuses=["blocked_obscurer"], component_id=6)
        self.assertEqual([item[0]["photon_id"] for item in chosen], ["2"])
        self.csv.write_text("status,point_count,x0_m,y0_m,z0_m\nblocked_mast,1,0,0,1\n")
        with self.assertRaisesRegex(ValueError, "requires recorded"):
            plot.select_trace_rows(self.csv, component_id=6)

    def test_stratification_is_deterministic_and_retains_statuses(self):
        first = plot.select_trace_rows(self.csv, selection="stratified", seed=42, max_paths=3)
        second = plot.select_trace_rows(self.csv, selection="stratified", seed=42, max_paths=3)
        self.assertEqual(first, second)
        self.assertEqual(
            {item[1] for item in first}, {"detected", "blocked_obscurer", "missed_primary"}
        )

    def test_clipping_retains_recorded_direction(self):
        segments = plot.clip_recorded_path([(0, 0, 100), (0, 0, 0), (1, 0, 4)], 5)
        self.assertEqual(
            segments, [((0.0, 0.0, 5.0), (0.0, 0.0, 0.0)), ((0.0, 0.0, 0.0), (1.0, 0.0, 4.0))]
        )
        self.assertEqual(plot.clip_recorded_path([(100, 0, 1), (100, 0, 2)], 5), [])

    def test_cylinder_footprint_is_finite_and_has_physical_diameter(self):
        cylinder = plot.PlotObscurer(6, (-1, 0, 2), (1, 0, 2), 0.2)
        footprint = plot._cylinder_footprint(cylinder, 0, 1)
        self.assertEqual((min(x for x, _ in footprint), max(x for x, _ in footprint)), (-1.0, 1.0))
        self.assertAlmostEqual(max(y for _, y in footprint) - min(y for _, y in footprint), 0.2)

    def test_hexagon_orientation_and_circle_diameter(self):
        horizontal = plot._aperture_local_vertices("hexagon_flat_x", 2)
        vertical = plot._aperture_local_vertices("hexagon_flat_y", 2)
        self.assertAlmostEqual(max(x for x, _ in horizontal), 1)
        self.assertAlmostEqual(max(y for _, y in vertical), 1)
        circle = plot._aperture_local_vertices("circle", 2)
        self.assertEqual(len(circle), 32)
        self.assertAlmostEqual(max(x for x, _ in circle) - min(x for x, _ in circle), 2)

    def test_export_svg_is_byte_stable(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest("Matplotlib optional")
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        first, second = self.directory / "first.svg", self.directory / "second.svg"
        plot.draw_pupil(plt, self.model, first)
        plot.draw_pupil(plt, self.model, second)
        self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_axisymmetric_section_does_not_bridge_aperture_hole(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest("Matplotlib optional")
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        surface = plot.AxisymmetricSurface("secondary", 4, 1, 2, 1, (0,) * 13)
        optical_model = plot.PlotOpticalModel(
            "axisymmetric", "fixture", None, (), (surface,), (), ()
        )
        figure, axis = plt.subplots()
        plot._draw_optical_model_section(axis, optical_model)
        curves = axis.lines[:2]
        self.assertEqual(len(curves), 2)
        self.assertTrue(all(x <= -1 for x in curves[0].get_xdata()))
        self.assertTrue(all(x >= 1 for x in curves[1].get_xdata()))
        plt.close(figure)

    def test_native_version_runs_shared_contract_before_any_plot_reader(self):
        # Versioned files must not silently fall through the permissive legacy adapter.
        self.csv.write_text(
            "contract_version,status,point_count,x0_m,y0_m,z0_m\n"
            "obdeect-arrival-v1,detected,1,0,0,0\n"
        )
        for reader in (plot.read_paths, plot.read_trace_rows, plot.focal_plane_hits):
            with self.subTest(reader=reader.__name__):
                with self.assertRaisesRegex(ValueError, "missing columns"):
                    list(reader(self.csv))
        from obdeect.analysis import analyse_trace_csv

        with self.assertRaisesRegex(ValueError, "missing columns"):
            analyse_trace_csv(self.csv)

    def test_bounded_material_diagnostics_keep_all_recorded_vertices(self):
        row = {"point_count": "65"}
        row.update({f"{axis}{index}_m": str(index) for index in range(65) for axis in "xyz"})
        points = plot._parse_path_row(self.csv, row)
        self.assertEqual(len(points), 65)
        self.assertEqual(points[-1], (64.0, 64.0, 64.0))
        row["point_count"] = "66"
        with self.assertRaisesRegex(ValueError, "between 1 and 65"):
            plot._parse_path_row(self.csv, row)

    def test_native_readers_use_validated_arrivals_and_keep_extra_metadata(self):
        self.csv.write_text(
            "contract_version,photon_id,source_kind,wavelength_nm,emission_time_ns,source_weight,"
            "throughput,status,point_count,path_length_m,x0_m,y0_m,z0_m,x1_m,y1_m,z1_m,"
            "custom_diagnostic\n"
            "obdeect-arrival-v1,1,laser,400,0,4,0.5,detected,2,6,0,0,10,0.1,0,4,retained\n"
        )
        with patch.object(plot, "_parse_path_row", side_effect=AssertionError("CSV reparsed")):
            self.assertEqual(
                list(plot.read_paths(self.csv)),
                [("detected", [(0.0, 0.0, 10.0), (0.1, 0.0, 4.0)])],
            )
            self.assertEqual(list(plot.focal_plane_hits(self.csv)), [(0.1, 0.0, 2.0)])
        self.assertEqual(
            list(plot.read_trace_rows(self.csv))[0][0]["custom_diagnostic"], "retained"
        )
        # Full iteration must still reject a later invalid native record.
        self.csv.write_text(self.csv.read_text() + self.csv.read_text().splitlines()[1] + "\n")
        for reader in (plot.read_paths, plot.read_trace_rows, plot.focal_plane_hits):
            with self.subTest(reader=reader.__name__):
                with self.assertRaisesRegex(ValueError, "duplicate photon_id"):
                    list(reader(self.csv))

    def test_camera_footprints_remain_two_dimensional_and_require_types(self):
        fixture = model_fixture()
        fixture["camera"] = {
            "rotation_deg": 90,
            "pixel_types": [{"id": 1, "funnel_shape_code": 2, "funnel_diameter_m": 0.1}],
            "pixels": [{"id": 9, "type_id": 1, "centre_xy_m": [0.2, 0]}],
        }
        self.model.write_text(json.dumps(fixture))
        optical_model = plot.load_plot_optical_model(self.model)
        pixel = optical_model.camera_pixels[0]
        self.assertEqual(len(pixel.vertices_xy_m), 4)
        self.assertAlmostEqual(sum(y for _, y in pixel.vertices_xy_m) / 4, 0.2)
        self.assertEqual(len(optical_model.polygons), 2)
        fixture["camera"].pop("pixel_types")
        self.model.write_text(json.dumps(fixture))
        self.assertEqual(plot.load_plot_optical_model(self.model).camera_pixels, ())

    def test_legacy_hexagonal_footprints_agree_with_compiled_entrances(self):
        for shape in (1, 3):
            with self.subTest(shape=shape):
                camera = {
                    "pixel_types": [
                        {
                            "id": 1,
                            "funnel_shape_code": shape,
                            "cathode_shape_code": shape,
                            "funnel_diameter_m": 0.1,
                            "cathode_diameter_m": 0.1,
                            "funnel_depth_m": 0,
                        }
                    ],
                    "pixels": [
                        {
                            "id": 9,
                            "type_id": 1,
                            "centre_xy_m": [0, 0],
                            "z_offset_m": 0,
                            "rotation_deg": 0,
                            "normal_slopes": [0, 0],
                            "module": 0,
                            "enabled": True,
                        }
                    ],
                }
                legacy = plot._camera_pixel_polygons({"camera": camera})
                camera.update(
                    compile_camera_surfaces(
                        camera, {"coefficient_m": [0], "radial_scale_m": 1}, 1, reflected=False
                    )
                )
                compiled = plot._camera_pixel_polygons({"camera": camera})
                self.assertEqual(legacy, compiled)

    def test_finite_aspheric_masks_keep_gap_and_tangent_orientation(self):
        surface = plot.AxisymmetricSurface("primary", 2, 0, 3, 2, (0, 4) + (0,) * 11)
        hexagon = plot._aspheric_mask_polygon(
            {
                "id": 1,
                "shape": "hexagon",
                "centre_xy_m": [1, 0],
                "diameter_m": 0.2,
                "rotation_deg": 0,
            },
            surface,
        )
        self.assertEqual(len(hexagon.vertices_m), 6)
        # z = 2 + r^2; the centre tangent has slope dz/dr=2 at r=1.
        for x, _, z in hexagon.vertices_m:
            self.assertAlmostEqual(z - 3, 2 * (x - 1))
        ring = plot._aspheric_mask_polygon(
            {
                "id": 2,
                "shape": "annular_sector",
                "inner_radius_m": 1,
                "outer_radius_m": 2,
                "start_deg": 0,
                "span_deg": 90,
                "gap_m": 0.1,
            },
            surface,
        )
        import math

        for x, y, z in ring.vertices_m:
            radius = math.hypot(x, y)
            self.assertGreaterEqual(radius, 1 - 1e-12)
            self.assertLessEqual(radius, 2 + 1e-12)
            self.assertLessEqual(math.atan2(y, x), math.pi / 2 - 0.1 / radius + 1e-12)
            self.assertAlmostEqual(z, 2 + radius**2)

    def test_loss_filter_keeps_incident_denominator(self):
        samples = list(plot.entrance_samples(self.csv))
        fractions = plot.weighted_loss_bins(samples, 2, 1, component_id=6)
        self.assertAlmostEqual(fractions[1][1], 1 / 6)

    def test_numeric_path_colours_require_recorded_fields(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest("Matplotlib optional")
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        records = plot.select_trace_rows(self.csv)
        with self.assertRaisesRegex(ValueError, "recorded arrival_time_ns"):
            plot._trace_colour_scale(plt, records, "arrival-time")
        for row, _, _ in records:
            row["arrival_time_ns"] = "42"
        _, normaliser, label = plot._trace_colour_scale(plt, records, "arrival-time")
        self.assertEqual(normaliser.vmin, 42)
        self.assertIn("[ns]", label)

    def test_cli_rejects_invalid_input_dependent_options(self):
        options = (
            ["--selection-seed", "42"],
            ["--selection-seed", "0"],
            ["--path-selection", "first"],
            ["--show-paths"],
            ["--context-radius-m", "5"],
            ["--status", "invented-status"],
            ["--component-id", "-1"],
            ["--view", "section", "--panel", "assembly"],
        )
        for option in options:
            with self.subTest(option=option):
                command = [
                    sys.executable,
                    "-m",
                    "obdeect.plotting",
                    "--view",
                    "telescope",
                    "--optical-model-json",
                    str(self.model),
                    "--output",
                    str(self.directory / "bad.png"),
                ]
                result = subprocess.run(
                    command + option,
                    capture_output=True,
                    text=True,
                    env={**os.environ, "PYTHONPATH": str(Path(plot.__file__).parents[1])},
                )
                self.assertEqual(result.returncode, 2)
                self.assertFalse((self.directory / "bad.png").exists())

    def test_mask_bounds_do_not_expand_to_uncovered_parent_annulus(self):
        surface = plot.AxisymmetricSurface("primary", 0, 0, 10, 1, (0,) * 13)
        footprint = plot._aspheric_mask_polygon(
            {
                "id": 1,
                "shape": "hexagon",
                "centre_xy_m": [1, 0],
                "diameter_m": 0.2,
                "rotation_deg": 0,
            },
            surface,
        )
        optical_model = plot.PlotOpticalModel(
            "axisymmetric", "fixture", None, (footprint,), (surface,), (), ()
        )
        bounds = plot._optical_model_bounds(optical_model)
        self.assertLess(bounds[0][1], 1.2)
        self.assertGreater(bounds[0][0], 0.8)

    def test_panel_selection_overlays_and_pdf_exports(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest("Matplotlib optional")
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for panel in ("assembly", "section", "pupil", "focal-geometry"):
            with self.subTest(panel=panel):
                output = self.directory / f"{panel}.svg"
                plot.draw_telescope(
                    plt, self.model, output, panel=panel, input_path=self.csv, context_radius_m=5
                )
                self.assertIn("finite-fixture", output.read_text())
                self.assertIn("[m]", output.read_text())
        first, second = self.directory / "plate-first.svg", self.directory / "plate-second.svg"
        plot.draw_telescope(plt, self.model, first, input_path=self.csv, context_radius_m=5)
        plot.draw_telescope(plt, self.model, second, input_path=self.csv, context_radius_m=5)
        self.assertEqual(first.read_bytes(), second.read_bytes())
        pdf = self.directory / "plate.pdf"
        plot.draw_telescope(plt, self.model, pdf)
        self.assertTrue(pdf.read_bytes().startswith(b"%PDF"))
        # Focal aspheric circles must autoscale to include both negative axes.
        surface = plot.AxisymmetricSurface("detector", 4, 0, 0.3, 1, (0,) * 13)
        optical_model = plot.PlotOpticalModel(
            "axisymmetric", "fixture", None, (), (surface,), (), ()
        )
        figure, axis = plt.subplots()
        plot._draw_focal_geometry(axis, optical_model)
        self.assertLess(axis.get_xlim()[0], -0.3)
        self.assertGreater(axis.get_ylim()[1], 0.3)
        plt.close(figure)

    def test_compiled_focal_svg_uses_model_provenance(self):
        output = self.directory / "focal.svg"
        plot.draw_focal_plane_svg(self.csv, output, 8, "incorrect-user-label", self.model)
        svg = output.read_text()
        self.assertIn("finite-fixture/1 focal-plane PSF", svg)
        self.assertNotIn("incorrect-user-label", svg)
        self.assertIn("Compiled model SHA256: " + "a" * 64, svg)

    def test_focal_svg_excludes_hits_outside_compiled_extent(self):
        self.csv.write_text(
            "status,point_count,source_weight,throughput,x0_m,y0_m\n"
            "detected,1,1,1,0,0\n"
            "detected,1,100,1,5,5\n"
        )
        output = self.directory / "focal.svg"
        plot.draw_focal_plane_svg(self.csv, output, 8, "fixture", self.model)
        # An out-of-frame hit must not be clamped into a false edge-bin hit.
        self.assertEqual(output.read_text().count('height="65.0000" fill='), 1)

    def test_rays_share_recorded_overlay_and_model_footer(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        figures = []
        with patch.object(
            plot, "_save_figure", side_effect=lambda figure, *_a, **_kw: figures.append(figure)
        ):
            plot.draw_rays(
                plt,
                self.csv,
                "incorrect-user-label",
                self.directory / "rays.svg",
                3,
                self.model,
                context_radius_m=5,
            )
        figure = figures[0]
        self.assertTrue(any("SHA256: aaaaaaaaaaaa" in text.get_text() for text in figure.texts))
        first, second = figure.axes
        for axis in (first, second):
            labels = axis.get_legend_handles_labels()[1]
            self.assertIn("detected", labels)
            self.assertIn("incoming path clipped", labels)
        terminal_x = [line.get_xdata()[0] for line in first.lines if line.get_marker() == "x"]
        terminal_y = [line.get_xdata()[0] for line in second.lines if line.get_marker() == "x"]
        self.assertEqual(terminal_x, [0.1, 0.0])
        self.assertEqual(terminal_y, [0.0, 0.0])

    def test_cylinder_along_plane_normal_does_not_gain_axial_radius(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # The cylinder occupies y=[0.05, 1], so its flat end never touches y=0.
        cylinder = plot.PlotObscurer(6, (0, 0.05, 2), (0, 1, 2), 0.2)
        model = plot.PlotOpticalModel("segmented", "fixture", None, (), (), (cylinder,), ())
        figure, axis = plt.subplots()
        plot._draw_optical_model_section(axis, model)
        self.assertEqual(len(axis.patches), 0)
        self.assertEqual(len(axis.artists), 1)  # Physical metre scale bar.
        plt.close(figure)

    def test_all_geometry_views_render_headless_with_real_footprints(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest("Matplotlib optional")
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for name, renderer in (
            ("plate", plot.draw_telescope),
            ("pupil", plot.draw_pupil),
            ("section", plot.draw_section),
            ("assembly", plot.draw_assembly_3d),
            ("legacy-mirror", plot.draw_compiled_mirror),
        ):
            with self.subTest(view=name):
                output = self.directory / f"{name}.png"
                renderer(plt, self.model, output)
                self.assertTrue(output.read_bytes().startswith(b"\x89PNG"))
        output = self.directory / "loss.svg"
        plot.draw_loss_map(plt, self.csv, self.model, output, bins=4)
        self.assertIn("lost input weight / incident-bin weight", output.read_text())
        plot.draw_telescope(
            plt,
            self.model,
            self.directory / "overlay.png",
            input_path=self.csv,
            panel="focal-plane-hits",
            context_radius_m=5,
        )
        self.assertFalse(plt.get_fignums())

    def test_cli_writes_files_without_pyplot_with_interactive_backend(self):
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            self.skipTest("Matplotlib optional")
        script = """
import importlib.abc
import sys
from pathlib import Path

class NoGuiImports(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "matplotlib.pyplot" or fullname.startswith((
            "matplotlib.backends.backend_macosx", "matplotlib.backends._macosx",
            "matplotlib.backends.backend_qt", "matplotlib.backends.backend_tk",
            "matplotlib.backends._backend_tk", "matplotlib.backends.backend_gtk",
            "matplotlib.backends.backend_wx",
        )):
            raise AssertionError(f"GUI import attempted: {fullname}")

sys.meta_path.insert(0, NoGuiImports())
import matplotlib
matplotlib.rcParams.update({"backend": "MacOSX", "interactive": True})
from obdeect.plotting import main

model, trace, directory = sys.argv[1:]
for extension in ("png", "pdf", "svg"):
    sys.argv = ["obdeect-plot-reference", "--view", "telescope",
                "--optical-model-json", model, "--input", trace,
                "--colour-by", "wavelength", "--output",
                str(Path(directory) / f"headless.{extension}")]
    main()
assert "matplotlib.pyplot" not in sys.modules
assert matplotlib.rcParams["backend"] == "MacOSX"
"""
        rows = self.csv.read_text().splitlines()
        self.csv.write_text(
            rows[0] + ",wavelength_nm\n" + "\n".join(row + ",400" for row in rows[1:]) + "\n"
        )
        subprocess.run(
            [sys.executable, "-c", script, str(self.model), str(self.csv), str(self.directory)],
            check=True,
            env={**os.environ, "PYTHONPATH": str(Path(plot.__file__).parents[1])},
        )
        self.assertTrue((self.directory / "headless.png").read_bytes().startswith(b"\x89PNG"))
        self.assertTrue((self.directory / "headless.pdf").read_bytes().startswith(b"%PDF"))
        self.assertIn("<svg", (self.directory / "headless.svg").read_text())
        self.assertNotIn("A. orthographic assembly", (self.directory / "headless.svg").read_text())

    def test_structure_cli_uses_model_geometry_and_selected_section_plane(self):
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            self.skipTest("Matplotlib optional")
        output = self.directory / "structure.svg"
        command = [
            sys.executable,
            "-m",
            "obdeect.plotting",
            "--view",
            "structure",
            "--section-plane",
            "yz",
            "--output",
            str(output),
        ]
        env = {**os.environ, "PYTHONPATH": str(Path(plot.__file__).parents[1])}
        subprocess.run(command + ["--optical-model-json", str(self.model)], check=True, env=env)
        svg = output.read_text()
        self.assertIn("finite-fixture/1 compiled optical section", svg)
        self.assertIn("telescope y [m]", svg)
        self.assertIn("a" * 12, svg)
        self.assertIn("opaque cylinder #6", svg)
        self.assertNotIn("spherical primary", svg)
        self.assertNotIn("camera body projection", svg)

        # A supplied invalid model must not fall back to a reference outline.
        self.model.write_text("{}")
        result = subprocess.run(
            command + ["--optical-model-json", str(self.model)], capture_output=True, env=env
        )
        self.assertNotEqual(result.returncode, 0)

        subprocess.run(command, check=True, env=env)
        self.assertIn("reference-mst reference mirror, camera and supports", output.read_text())


if __name__ == "__main__":
    unittest.main()
