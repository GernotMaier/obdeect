// Optical-model loading: valid documents load; every error names its path.
#include "check.hpp"

#include "obdeect/model_json.hpp"

#include <string>

using namespace obdeect;

namespace {

const std::string kValid = R"({
  "format": "obdeect.optical-model.v2",
  "sha256": "ignored by the loader",
  "provenance": {"model": "test"},
  "notes": [],
  "optics": {
    "propagation_group_index": 1.0003,
    "telescope_transmission": {"on_axis": 0.9, "amplitude": 0.0, "angular_scale_rad": 0.1,
                               "power": 2.0, "outer_power": 1.0},
    "primary": {
      "facets": [{"id": 0, "centre_m": [0, 0, 0], "normal": [0, 0, 1], "tangent": [1, 0, 0],
                  "shape": "hexagon_flat_y", "diameter_m": 1.2, "radius_of_curvature_m": 32.0}],
      "reflectivity": {"wavelength_min_nm": 300, "wavelength_step_nm": 1, "incidence_deg": [0],
                       "values": [[0.8], [0.9]]},
      "scatter": {"sigma1_rad": 0.0001, "fraction2": 0.2, "sigma2_rad": 0.0004,
                  "method": "outgoing_angles"},
      "degradation": {"x_m": [-1, 1], "y_m": [-1, 1], "values": [[1, 1], [1, 0.9]], "outside": "clamp"}
    },
    "secondary": null,
    "secondary_blocks_incoming": false,
    "focal_surface": {"id": 1, "plane_z_m": 16.0},
    "obscurers": [
      {"id": 2, "name": "mast-0", "group": "camera_support", "legs": ["incoming", "after_primary"],
       "shape": "cylinder", "start_m": [3, 0, 0], "end_m": [0, 0, 16], "diameter_m": 0.2},
      {"id": 3, "name": "box", "group": "dish_structure", "legs": ["incoming"], "shape": "polyhedron",
       "vertices_m": [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
       "planes": [[0, 0, -1, 0], [0, -1, 0, 0], [-1, 0, 0, 0],
                  [0.5773502691896258, 0.5773502691896258, 0.5773502691896258, 0.5773502691896258]]},
      {"id": 4, "name": "camera", "group": "camera", "legs": ["incoming"], "shape": "plate",
       "centre_m": [0, 0, 16], "normal": [0, 0, 1], "tangent": [1, 0, 0], "aperture": "square",
       "diameter_m": 2.0}
    ]
  }
})";

// Load `text` with `from` replaced by `to`; return the error message, if any.
std::string error_after(const std::string &from, const std::string &to) {
  std::string text = kValid;
  const auto at = text.find(from);
  if (at == std::string::npos)
    return "test error: '" + from + "' not found";
  text.replace(at, from.size(), to);
  try {
    (void)load_model(text);
  } catch (const ModelError &error) {
    return error.what();
  }
  return "";
}

bool contains(const std::string &text, const std::string &part) {
  return text.find(part) != std::string::npos;
}

} // namespace

void all_tests() {
  const auto model = load_model(kValid);
  CHECK(model.primary.facets.size() == 1);
  CHECK_NEAR(model.primary.facets[0].half_width, 0.6, 0.0);
  CHECK(model.primary.reflectivity && model.primary.reflectivity->wavelength_count == 2);
  CHECK(model.obscurers.items().size() == 3 && model.obscurers.planes().size() == 4);
  CHECK(model.obscurers.items()[0].legs == (Leg::incoming | Leg::after_primary));
  CHECK(!model.secondary && model.telescope_transmission);

  CHECK(contains(
      error_after(R"("normal": [0, 0, 1], "tangent")", R"("normal": [0, 0, 2], "tangent")"),
      "model.optics.primary.facets[0].normal: expected a unit vector"));
  CHECK(contains(error_after(R"("format": "obdeect.optical-model.v2")", R"("format": "v1")"),
                 "model.format"));
  CHECK(contains(error_after(R"("notes": [],)", R"("notes": [], "extra": 1,)"),
                 "unknown key 'extra'"));
  CHECK(contains(error_after(R"("id": 3, "name": "box")", R"("id": 2, "name": "box")"),
                 "used twice"));
  CHECK(contains(error_after(R"("values": [[0.8], [0.9]])", R"("values": [[0.8], [1.9]])"),
                 "reflectivity.values[1][0]: expected a value in [0, 1]"));
  CHECK(contains(error_after(R"("legs": ["incoming"], "shape": "plate")",
                             R"("legs": ["sideways"], "shape": "plate")"),
                 "obscurers[2].legs[0]: unknown leg"));
  CHECK(contains(error_after("0.5773502691896258]]", "0.1]]"), "vertex lies outside a face plane"));
  CHECK(contains(error_after(R"("plane_z_m": 16.0)", R"("asphere_x": 1)"),
                 "focal_surface: unknown key"));
  CHECK(contains(error_after("\"optics\": {", "\"optics\": {,"), "JSON syntax error"));
}

int main() { return run_checks(all_tests); }
