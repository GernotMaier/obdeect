"""Analytic prescriptions from sim_telarray manual sections 11 and 12.3.3."""

import copy
import math
import unittest

from obdeect.optical_model_compiler import (
    OpticalModelCompileError,
    _apply_mirror_degradation,
    compile_secondary_baffle,
    derive_nominal_single_reflector,
    parse_obscuration_quadrilaterals,
    parse_rectangular_response,
    parse_simtel_mirror_list,
    parse_spatial_response,
    resolve_focus_offset,
    resolve_panel_prescription,
)
from obdeect.pixel_response_compiler import measured_pixel_response
from obdeect.simtel_tables import TableImportError, parse_rpol_table


class TestManualPanelPrescriptions(unittest.TestCase):
    def test_degradation_scales_polynomials_and_logarithms_once(self):
        for logarithmic in (False, True):
            source = {
                "reflectivity": [
                    {
                        "wavelength_nm": 300,
                        "response": 0.4,
                        "interpolation": {
                            "value_log": logarithmic,
                            "coefficients": [
                                [math.log(0.4), 1, 2, 3] if logarithmic else [0.4, 1, 2, 3]
                            ],
                        },
                    },
                    {"wavelength_nm": 500, "response": 0.8},
                ]
            }
            _apply_mirror_degradation(source, dict(value=0.5, unit=""), "mirror_degradation")
            self.assertEqual(source["reflectivity"][0]["response"], 0.2)
            actual = source["reflectivity"][0]["interpolation"]["coefficients"][0]
            expected = [math.log(0.2), 1, 2, 3] if logarithmic else [0.2, 0.5, 1, 1.5]
            for value, target in zip(actual, expected):
                self.assertAlmostEqual(value, target)
            _apply_mirror_degradation(source, dict(value=0, unit="null"), "mirror_degradation")
            self.assertFalse(source["reflectivity"][0]["interpolation"]["value_log"])
            self.assertEqual(source["reflectivity"][0]["interpolation"]["coefficients"], [[0] * 4])

    def test_constant_pixel_efficiency_and_absolute_spectral_normalization(self):
        table = measured_pixel_response(4, "0 0.5\n", "#THETA: 0 30\n400 0.2\n")
        self.assertEqual(table["angular_efficiency"], [0.5] * 1000)
        for factor in table["spectral_correction"]:
            self.assertAlmostEqual(factor, 0.4)
        self.assertEqual(table["wavelength_bin_origin_nm"], 0.5)
        with self.assertRaises(TableImportError):
            measured_pixel_response(4, "0 nan\n")

    def test_measured_angular_bins_sample_tangent_midpoints(self):
        table = measured_pixel_response(7, "0 0\n90 1\n")
        for index in (0, 19, 999):
            expected = math.degrees(math.atan(0.01 * (index + 0.5))) / 90
            self.assertAlmostEqual(table["angular_efficiency"][index], expected)

    def test_table_options_and_quadratic_coefficients(self):
        table = parse_rpol_table("#@RPOL@ 1 scheme=2 clip\n0 0\n1 1\n2 4\n3 9\n")
        self.assertEqual(table["interpolation"]["boundary"], "zero")
        for index, coefficients in enumerate(table["interpolation"]["coefficients"]):
            for fraction in (0, 0.25, 0.8, 1):
                value = sum(c * fraction**power for power, c in enumerate(coefficients))
                self.assertAlmostEqual(value, (index + fraction) ** 2)
        log_table = parse_rpol_table(
            "#@RPOL@ 1 xlog ylog clip\n1 0.1\n10 1\n", filename_options="noclip"
        )
        self.assertTrue(log_table["interpolation"]["x_log"])
        self.assertTrue(log_table["interpolation"]["value_log"])
        self.assertEqual(log_table["interpolation"]["boundary"], "clamp")

    def parameters(self):
        return {
            "focal_length": {"value": 16, "unit": "m"},
            "dish_shape_length": {"value": 16, "unit": "m"},
            "mirror_offset": {"value": 0, "unit": "m"},
            "parabolic_dish": {"value": True},
        }

    def test_reference_table_scaling_and_option_precedence(self):
        table = parse_rpol_table(
            "#@RPOL@ 1 yscale=100 zscale=2\n1 0.2 0.3\n2 0.4 0.5\n",
            api_options="zcol=3,zscale=3",
            filename_options="zscale=1",
        )
        self.assertEqual(table["response"], [0.3, 0.5])
        logarithmic = parse_rpol_table("#@RPOL@ 1 xlog zlog xscale=2 zscale=2\n1 0.2\n2 0.4\n")
        self.assertEqual(logarithmic["x"], [1, 4])
        for actual, expected in zip(logarithmic["response"], [0.04, 0.16]):
            self.assertAlmostEqual(actual, expected)
        self.assertEqual(
            parse_rpol_table("#@RPOL@ 1 zscale=0\n1 0.2\n2 0.4\n")["response"],
            [0.2, 0.4],
        )

    def test_automatic_and_signed_focal_lengths_have_distinct_semantics(self):
        facets = parse_simtel_mirror_list(
            "200 0 100 0 % auto\n300 0 100 -1600\n400 0 100 1600\n",
            fallback_focal_length_m=None,
            allow_automatic=True,
        )
        parameters = self.parameters()
        parameters["random_focal_length"] = {"value": [0.1, 0.2], "unit": ["m", "m"]}
        before = copy.deepcopy(facets)
        resolve_panel_prescription(facets, parameters, 4)
        automatic_distance = math.hypot(16 - 2**2 / 64, 2)
        self.assertLessEqual(abs(facets[0]["focal_length_m"] - automatic_distance), 0.2)
        self.assertLessEqual(abs(facets[1]["focal_length_m"] - 16), 0.2)
        self.assertNotEqual(facets[1]["focal_length_m"], 16)
        self.assertEqual(facets[2]["focal_length_m"], 16)
        resolve_panel_prescription(before, parameters, 4)
        self.assertEqual(facets, before)

    def test_positive_height_prevents_random_distance_but_negative_height_does_not(self):
        facets = parse_simtel_mirror_list(
            "200 0 100 1600 0 100\n200 0 100 1600 0 -100\n",
            fallback_focal_length_m=None,
        )
        parameters = self.parameters()
        parameters["mirror_align_random_distance"] = {"value": 0.1, "unit": "m"}
        resolve_panel_prescription(facets, parameters, 3)
        self.assertEqual(facets[0]["resolved_focus_distance_m"], math.hypot(15, 2))
        self.assertNotEqual(facets[1]["resolved_focus_distance_m"], math.hypot(15, 2))
        derive_nominal_single_reflector(facets, parameters)
        self.assertEqual(facets[0]["nominal_centre_m"][2], 1)
        self.assertNotEqual(facets[1]["nominal_centre_m"][2], 1)

    def test_grading_flip_and_explicit_scale(self):
        facets = parse_simtel_mirror_list(
            "100 0 100 0 1\n300 0 100 0 3\n500 0 100 1000 2\n",
            fallback_focal_length_m=None,
            allow_automatic=True,
        )
        parameters = self.parameters()
        parameters["flip_mirrors"] = {"value": True}
        parameters["grading_of_focal_length"] = {"value": 0.4, "unit": "m"}
        parameters["mirror_f_scale"] = {"value": [1, 1, 1.1]}
        resolve_panel_prescription(facets, parameters, 0)
        self.assertEqual(facets[0]["centre_m"][:2], [0, 1])
        self.assertEqual(facets[0]["shape"], "hexagon_flat_y")
        self.assertAlmostEqual(facets[0]["focal_length_m"], math.hypot(16 - 1 / 64, 1) - 0.2)
        self.assertAlmostEqual(facets[2]["focal_length_m"], 11)

    def test_focus_deformation_is_linear_and_requires_zenith_context(self):
        parameter = {"value": [2, 30, 3, 4], "unit": ["cm", "deg", "cm", "cm"]}
        self.assertAlmostEqual(resolve_focus_offset(parameter, 30), 0.02)
        expected = 0.02 + 0.03 * (1 - math.cos(math.radians(30))) - 0.04 * 0.5
        self.assertAlmostEqual(resolve_focus_offset(parameter, 0), expected)
        with self.assertRaisesRegex(OpticalModelCompileError, "zenith"):
            resolve_focus_offset(parameter, None)

    def test_rectangular_maps_and_baffle_units(self):
        table = "#@RPOL@ 3\n-100 -100 0.2\n-100 100 0.4\n100 -100 0.6\n100 100 0.8\n"
        spatial = parse_spatial_response(table, "degradation")
        self.assertEqual(spatial["x_m"], [-1, 1])
        self.assertEqual(spatial["response"], [0.2, 0.4, 0.6, 0.8])
        matrix = "#@RPOL@[#Y=] 2\n#Y= -100 100\n-100 0.2 0.4\n100 0.6 0.8\n"
        self.assertEqual(
            parse_rectangular_response(table, "map"), parse_rectangular_response(matrix, "map")
        )
        with self.assertRaisesRegex(OpticalModelCompileError, "complete grid"):
            parse_spatial_response(table.rsplit("100 100", 1)[0], "map")
        baffle = compile_secondary_baffle({"value": [100, 200, 10, 1, 20], "unit": ["cm"] * 5}, 0.5)
        self.assertEqual(baffle["first_endpoint_m"], [0, 0, 0.5])
        self.assertEqual(baffle["second_radius_m"], 0.2)
        self.assertEqual(baffle["thickness_m"], 0.01)
        self.assertIsNone(compile_secondary_baffle({"value": [0, 0, 0], "unit": "cm"}, 0))

    def test_quadrilateral_vertex_import_keeps_metres(self):
        text = "id group x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4\n"
        text += "panel structure -1 -1 2 1 -1 2 1 1 2 -1 1 2\n"
        plate = parse_obscuration_quadrilaterals(text)[0]
        self.assertEqual(plate["vertices_m"][0], [-1, -1, 2])
        with self.assertRaisesRegex(OpticalModelCompileError, "planar"):
            parse_obscuration_quadrilaterals(text.replace("-1 1 2", "-1 1 3"))
        with self.assertRaisesRegex(OpticalModelCompileError, "convex"):
            parse_obscuration_quadrilaterals(text.replace("1 1 2", "0 -0.5 2", 1))
