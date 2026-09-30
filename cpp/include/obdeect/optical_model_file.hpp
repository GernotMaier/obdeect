#pragma once

#include "obdeect/axisymmetric_optics.hpp"
#include "obdeect/json.hpp"
#include "obdeect/photon_buffer.hpp"
#include "obdeect/segmented_optical_model.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <numbers>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace obdeect {

struct AxisymmetricOpticalModel {
  ModelProvenance provenance;
  AxisymmetricMirror primary;
  AxisymmetricMirror secondary;
  AxisymmetricMirror detector;
  std::optional<SpectralResponse> primary_reflectivity;
  std::optional<SpectralResponse> secondary_reflectivity;
};

namespace detail {

[[nodiscard]] inline const json::Value* field(const json::Value& value, std::string_view name) {
  return value.find(name);
}

[[nodiscard]] inline const std::string* string_field(const json::Value& value, std::string_view name) {
  const auto* result = field(value, name);
  return result && result->kind == json::Value::Kind::string ? &result->string : nullptr;
}

[[nodiscard]] inline std::optional<double> number_field(const json::Value& value, std::string_view name) {
  const auto* result = field(value, name);
  if (!result || result->kind != json::Value::Kind::number || !std::isfinite(result->number)) return std::nullopt;
  return result->number;
}

[[nodiscard]] inline std::optional<std::uint32_t> uint_field(const json::Value& value, std::string_view name) {
  const auto number = number_field(value, name);
  if (!number || *number < 0.0 || *number > std::numeric_limits<std::uint32_t>::max() ||
      std::floor(*number) != *number) return std::nullopt;
  return static_cast<std::uint32_t>(*number);
}

[[nodiscard]] inline std::optional<Vec3> vec3_field(const json::Value& value, std::string_view name) {
  const auto* result = field(value, name);
  if (!result || result->kind != json::Value::Kind::array || result->array.size() != 3) return std::nullopt;
  Vec3 vector{};
  double* components[] = {&vector.x, &vector.y, &vector.z};
  for (std::size_t index = 0; index < 3; ++index) {
    if (result->array[index].kind != json::Value::Kind::number || !std::isfinite(result->array[index].number)) return std::nullopt;
    *components[index] = result->array[index].number;
  }
  return vector;
}

[[nodiscard]] inline std::optional<FacetShape> facet_shape(const json::Value& value) {
  const auto* shape = string_field(value, "shape");
  if (!shape) return std::nullopt;
  if (*shape == "circle") return FacetShape::circle;
  if (*shape == "hexagon_flat_y" || *shape == "hexagon") return FacetShape::hexagon_flat_y;
  if (*shape == "square") return FacetShape::square;
  if (*shape == "hexagon_flat_x") return FacetShape::hexagon_flat_x;
  return std::nullopt;
}

[[nodiscard]] inline std::optional<ModelProvenance> provenance(const json::Value& root) {
  const auto* source = field(root, "provenance");
  if (!source || source->kind != json::Value::Kind::object) return std::nullopt;
  const auto* model = string_field(*source, "model");
  const auto* version = string_field(*source, "model_version");
  const auto* hash = string_field(root, "optical_model_sha256");
  if (!model || !version || !hash) return std::nullopt;
  ModelProvenance result{*model, *version, *hash};
  return has_valid_provenance(result) ? std::optional{std::move(result)} : std::nullopt;
}

[[nodiscard]] inline std::optional<SpectralResponse> response_field(const json::Value& value, std::string_view name) {
  const auto* response = field(value, name);
  if (!response || response->kind != json::Value::Kind::array) return std::nullopt;
  SpectralResponse result;
  for (const auto& entry : response->array) {
    if (entry.kind != json::Value::Kind::object) return std::nullopt;
    const auto wavelength = number_field(entry, "wavelength_nm");
    const auto response_value = number_field(entry, "response");
    if (!wavelength || !response_value) return std::nullopt;
    result.wavelength_nm.push_back(*wavelength);
    result.response.push_back(*response_value);
  }
  return result.is_valid() ? std::optional{std::move(result)} : std::nullopt;
}

[[nodiscard]] inline std::optional<AxisymmetricMirror> axisymmetric_surface(const json::Value& value) {
  const auto vertex = number_field(value, "vertex_z_m");
  const auto inner = number_field(value, "inner_radius_m");
  const auto outer = number_field(value, "outer_radius_m");
  const auto scale = number_field(value, "radial_scale_m");
  const auto* coefficients = field(value, "coefficient_m");
  if (!vertex || !inner || !outer || !scale || !coefficients || coefficients->kind != json::Value::Kind::array ||
      coefficients->array.size() != 13) return std::nullopt;
  AxisymmetricMirror result{};
  result.vertex_z_m = *vertex;
  result.inner_radius_m = *inner;
  result.outer_radius_m = *outer;
  result.surface.radial_scale_m = *scale;
  for (std::size_t index = 0; index < result.surface.coefficient_m.size(); ++index) {
    const auto& coefficient = coefficients->array[index];
    if (coefficient.kind != json::Value::Kind::number || !std::isfinite(coefficient.number)) return std::nullopt;
    result.surface.coefficient_m[index] = coefficient.number;
  }
  return is_valid(result) ? std::optional{result} : std::nullopt;
}

[[nodiscard]] inline std::optional<json::Value> read_json(const std::string& path) {
  std::ifstream input(path);
  if (!input) return std::nullopt;
  std::stringstream buffer;
  buffer << input.rdbuf();
  const auto root = json::Parser{buffer.str()}.parse();
  if (!root || root->kind != json::Value::Kind::object || !string_field(*root, "format") ||
      *string_field(*root, "format") != "obdeect.compiled-optical-model.v1") return std::nullopt;
  return root;
}

}  // namespace detail

// Read the single versioned JSON artifact emitted by obdeect-compile-optical-model.
// Its trace_model section is directly consumable; the remainder is retained for
// provenance, validation and inspection.
inline std::optional<CompiledSegmentedOpticalModel> read_segmented_optical_model(const std::string& path) {
  const auto root = detail::read_json(path);
  if (!root) return std::nullopt;
  const auto model_provenance = detail::provenance(*root);
  const auto* trace = detail::field(*root, "trace_model");
  if (!model_provenance || !trace || trace->kind != json::Value::Kind::object ||
      !detail::string_field(*trace, "kind") || *detail::string_field(*trace, "kind") != "segmented") return std::nullopt;
  const auto* facet_values = detail::field(*trace, "primary_facets");
  const auto* detector_values = detail::field(*trace, "detector_surfaces");
  const auto* obscurer_values = detail::field(*trace, "cylinder_obscurers");
  if (!facet_values || !detector_values || !obscurer_values || facet_values->kind != json::Value::Kind::array ||
      detector_values->kind != json::Value::Kind::array || obscurer_values->kind != json::Value::Kind::array) return std::nullopt;
  std::vector<ImportedFacet> facets;
  for (const auto& value : facet_values->array) {
    if (value.kind != json::Value::Kind::object) return std::nullopt;
    const auto id = detail::uint_field(value, "id");
    const auto centre = detail::vec3_field(value, "centre_m");
    const auto normal = detail::vec3_field(value, "normal");
    const auto tangent = detail::vec3_field(value, "tangent");
    const auto diameter = detail::number_field(value, "diameter_m");
    const auto focal_length = detail::number_field(value, "focal_length_m");
    const auto shape = detail::facet_shape(value);
    if (!id || !centre || !normal || !tangent || !diameter || !focal_length || !shape) return std::nullopt;
    facets.push_back({*id, *centre, *normal, *diameter, *focal_length, *shape, *tangent, 2.0 * *focal_length});
  }
  std::vector<ImportedDetectorSurface> detectors;
  for (const auto& value : detector_values->array) {
    if (value.kind != json::Value::Kind::object) return std::nullopt;
    const auto id = detail::uint_field(value, "id");
    const auto centre = detail::vec3_field(value, "centre_m");
    const auto normal = detail::vec3_field(value, "normal");
    const auto tangent = detail::vec3_field(value, "tangent");
    const auto diameter = detail::number_field(value, "diameter_m");
    const auto shape = detail::facet_shape(value);
    if (!id || !centre || !normal || !tangent || !diameter || !shape) return std::nullopt;
    detectors.push_back({*id, *centre, *normal, *diameter, *shape, *tangent});
  }
  std::vector<ImportedCylinderObscurer> obscurers;
  for (const auto& value : obscurer_values->array) {
    if (value.kind != json::Value::Kind::object) return std::nullopt;
    const auto id = detail::uint_field(value, "id");
    const auto first = detail::vec3_field(value, "first_endpoint_m");
    const auto second = detail::vec3_field(value, "second_endpoint_m");
    const auto diameter = detail::number_field(value, "diameter_m");
    if (!id || !first || !second || !diameter) return std::nullopt;
    obscurers.push_back({*id, *first, *second, *diameter});
  }
  const auto reflectivity = detail::response_field(*trace, "primary_reflectivity");
  if (!reflectivity || facets.empty() || detectors.empty()) return std::nullopt;
  return compile_segmented_optical_model({*model_provenance, std::move(facets), std::move(detectors),
                                          std::move(obscurers), *reflectivity});
}

[[nodiscard]] inline std::optional<AxisymmetricOpticalModel> read_axisymmetric_optical_model(const std::string& path) {
  const auto root = detail::read_json(path);
  if (!root) return std::nullopt;
  const auto model_provenance = detail::provenance(*root);
  const auto* trace = detail::field(*root, "trace_model");
  if (!model_provenance || !trace || trace->kind != json::Value::Kind::object ||
      !detail::string_field(*trace, "kind") || *detail::string_field(*trace, "kind") != "axisymmetric") return std::nullopt;
  const auto* primary = detail::field(*trace, "primary");
  const auto* secondary = detail::field(*trace, "secondary");
  const auto* detector = detail::field(*trace, "detector");
  if (!primary || !secondary || !detector) return std::nullopt;
  const auto primary_surface = detail::axisymmetric_surface(*primary);
  const auto secondary_surface = detail::axisymmetric_surface(*secondary);
  const auto detector_surface = detail::axisymmetric_surface(*detector);
  const auto primary_reflectivity = detail::response_field(*trace, "primary_reflectivity");
  const auto secondary_reflectivity = detail::response_field(*trace, "secondary_reflectivity");
  if (!primary_surface || !secondary_surface || !detector_surface || !primary_reflectivity || !secondary_reflectivity)
    return std::nullopt;
  return AxisymmetricOpticalModel{*model_provenance, *primary_surface, *secondary_surface, *detector_surface,
                                  *primary_reflectivity, *secondary_reflectivity};
}

[[nodiscard]] inline PathRecord trace_axisymmetric_optical_model(const Ray& input, std::uint64_t photon_id,
                                                          const AxisymmetricOpticalModel& optical_model) {
  PathRecord record{};
  record.photon_id = photon_id;
  record.points_m[0] = input.position_m;
  record.point_count = 1;
  const auto direction = normalised_checked(input.direction);
  if (!direction) { record.status = PhotonStatus::invalid_input; return record; }
  Ray ray{input.position_m, *direction};
  const auto primary = intersect_axisymmetric_mirror(ray, optical_model.primary);
  if (!primary) { record.status = PhotonStatus::missed_primary; record.final_direction = ray.direction; return record; }
  record.points_m[1] = primary->point_m;
  record.point_count = 2;
  record.path_length_m = primary->distance_m;
  record.incidence_primary_deg = std::acos(std::clamp(std::abs(dot(ray.direction, primary->unit_normal)), 0.0, 1.0)) * 180.0 / std::numbers::pi;
  const auto after_primary = reflect(ray, *primary);
  if (!after_primary) { record.status = PhotonStatus::invalid_input; return record; }
  ray = *after_primary;
  const auto secondary = intersect_axisymmetric_mirror(ray, optical_model.secondary);
  if (!secondary) { record.status = PhotonStatus::missed_screen; record.final_direction = ray.direction; return record; }
  record.points_m[2] = secondary->point_m;
  record.point_count = 3;
  record.path_length_m += secondary->distance_m;
  record.incidence_secondary_deg = std::acos(std::clamp(std::abs(dot(ray.direction, secondary->unit_normal)), 0.0, 1.0)) * 180.0 / std::numbers::pi;
  const auto after_secondary = reflect(ray, *secondary);
  if (!after_secondary) { record.status = PhotonStatus::invalid_input; return record; }
  ray = *after_secondary;
  const auto detector = intersect_axisymmetric_mirror(ray, optical_model.detector);
  if (!detector) { record.status = PhotonStatus::missed_screen; record.final_direction = ray.direction; return record; }
  record.points_m[3] = detector->point_m;
  record.point_count = 4;
  record.path_length_m += detector->distance_m;
  record.incidence_focal_deg = std::acos(std::clamp(std::abs(dot(ray.direction, detector->unit_normal)), 0.0, 1.0)) * 180.0 / std::numbers::pi;
  record.final_direction = ray.direction;
  record.status = PhotonStatus::detected;
  return record;
}

}  // namespace obdeect
