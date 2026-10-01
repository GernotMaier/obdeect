#pragma once

#include "obdeect/axisymmetric_optics.hpp"
#include "obdeect/photon_buffer.hpp"
#include "obdeect/segmented_scene.hpp"

#include <fstream>
#include <cmath>
#include <limits>
#include <optional>
#include <numbers>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace obdeect {

struct AxisymmetricScene {
  ModelProvenance provenance;
  AxisymmetricMirror primary;
  AxisymmetricMirror secondary;
  AxisymmetricMirror detector;
  std::optional<SpectralResponse> primary_reflectivity;
  std::optional<SpectralResponse> secondary_reflectivity;
};

// Dependency-free interchange format emitted by the Python model adapter.
// Keeping parsing here deliberately small makes the wheel usable without a
// JSON library while retaining a strict, versioned boundary between model
// assets and the C++ tracing kernel.
inline std::vector<std::string> split_scene_fields(const std::string& line) {
  std::vector<std::string> fields;
  std::stringstream stream(line);
  std::string field;
  while (std::getline(stream, field, ',')) fields.push_back(field);
  return fields;
}

inline bool scene_number(std::string_view text, double& value) {
  try {
    std::size_t consumed = 0;
    value = std::stod(std::string{text}, &consumed);
    return consumed == text.size() && std::isfinite(value);
  } catch (...) {
    return false;
  }
}

inline bool scene_uint(std::string_view text, std::uint32_t& value) {
  try {
    std::size_t consumed = 0;
    const auto parsed = std::stoul(std::string{text}, &consumed);
    if (consumed != text.size() || parsed > std::numeric_limits<std::uint32_t>::max()) return false;
    value = static_cast<std::uint32_t>(parsed);
    return true;
  } catch (...) {
    return false;
  }
}

inline std::optional<FacetShape> scene_shape(std::string_view text) {
  if (text == "circle") return FacetShape::circle;
  if (text == "hexagon_flat_y" || text == "hexagon") return FacetShape::hexagon_flat_y;
  if (text == "square") return FacetShape::square;
  if (text == "hexagon_flat_x") return FacetShape::hexagon_flat_x;
  return std::nullopt;
}

inline std::optional<CompiledSegmentedScene> read_native_scene(const std::string& path) {
  std::ifstream input(path);
  if (!input) return std::nullopt;
  std::string line;
  const auto read_line = [&]() {
    if (!std::getline(input, line)) return false;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return true;
  };
  if (!read_line() || line != "obdeect-scene-v1") return std::nullopt;
  if (!read_line()) return std::nullopt;
  const auto provenance_fields = split_scene_fields(line);
  if (provenance_fields.size() != 4 || provenance_fields[0] != "provenance") return std::nullopt;
  ModelProvenance provenance{provenance_fields[1], provenance_fields[2], provenance_fields[3]};
  if (!has_valid_provenance(provenance)) return std::nullopt;
  if (!read_line() || (line != "surface_id,role,shape,cx_m,cy_m,cz_m,nx,ny,nz,tx,ty,tz,diameter_m" &&
                       line != "surface_id,role,shape,cx_m,cy_m,cz_m,nx,ny,nz,tx,ty,tz,diameter_m,focal_length_m"))
    return std::nullopt;
  const bool has_focal_length = line.ends_with(",focal_length_m");

  std::vector<ImportedFacet> facets;
  std::vector<ImportedDetectorSurface> detectors;
  std::vector<ImportedCylinderObscurer> obscurers;
  SpectralResponse primary_reflectivity;
  while (read_line()) {
    if (line.empty()) continue;
    const auto fields = split_scene_fields(line);
    if (fields.size() == 3 && fields[0] == "primary_reflectivity") {
      double wavelength = 0.0, response = 0.0;
      if (!scene_number(fields[1], wavelength) || !scene_number(fields[2], response)) return std::nullopt;
      primary_reflectivity.wavelength_nm.push_back(wavelength);
      primary_reflectivity.response.push_back(response);
      continue;
    }
    if (fields.size() == 9 && fields[0] == "obscurer_cylinder") {
      std::uint32_t id = 0;
      if (!scene_uint(fields[1], id)) return std::nullopt;
      double values[7]{};
      for (std::size_t index = 0; index < 7; ++index)
        if (!scene_number(fields[index + 2], values[index])) return std::nullopt;
      obscurers.push_back({id, {values[0], values[1], values[2]}, {values[3], values[4], values[5]},
                           values[6]});
      continue;
    }
    if (fields.size() != (has_focal_length ? 14U : 13U)) return std::nullopt;
    std::uint32_t id = 0;
    if (!scene_uint(fields[0], id)) return std::nullopt;
    const auto shape = scene_shape(fields[2]);
    if (!shape) return std::nullopt;
    double values[10]{};
    for (std::size_t index = 0; index < 10; ++index)
      if (!scene_number(fields[index + 3], values[index])) return std::nullopt;
    const Vec3 centre{values[0], values[1], values[2]};
    const Vec3 normal{values[3], values[4], values[5]};
    const Vec3 tangent{values[6], values[7], values[8]};
    if (fields[1] == "mirror") {
      // The focal length is only used for import validation. The native
      // segmented kernel traces the explicit plane placement and therefore
      // does not infer a prescription from it.
      const double focal_length_m = has_focal_length ? [&] {
        double value = 0.0;
        return scene_number(fields[13], value) ? value : std::numeric_limits<double>::quiet_NaN();
      }() : values[9];
      facets.push_back({id, centre, normal, values[9], focal_length_m, *shape, tangent,
                        has_focal_length ? 2.0 * focal_length_m : 0.0});
    } else if (fields[1] == "detector") {
      detectors.push_back({id, centre, normal, values[9], *shape, tangent});
    } else {
      return std::nullopt;
    }
  }
  if (facets.empty() || detectors.empty()) return std::nullopt;
  const std::optional<SpectralResponse> reflectivity = primary_reflectivity.wavelength_nm.empty()
                                                          ? std::nullopt
                                                          : std::optional{std::move(primary_reflectivity)};
  return compile_segmented_scene(ImportedSegmentedScene{std::move(provenance), std::move(facets),
                                                         std::move(detectors), std::move(obscurers), reflectivity});
}

[[nodiscard]] inline std::optional<AxisymmetricScene> read_native_axisymmetric_scene(const std::string& path) {
  std::ifstream input(path);
  if (!input) return std::nullopt;
  std::string line;
  const auto read_line = [&]() {
    if (!std::getline(input, line)) return false;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return true;
  };
  if (!read_line() || line != "obdeect-axisymmetric-scene-v1") return std::nullopt;
  if (!read_line()) return std::nullopt;
  const auto provenance_fields = split_scene_fields(line);
  if (provenance_fields.size() != 4 || provenance_fields[0] != "provenance") return std::nullopt;
  ModelProvenance provenance{provenance_fields[1], provenance_fields[2], provenance_fields[3]};
  if (!has_valid_provenance(provenance) || !read_line() ||
      line != "role,vertex_z_m,inner_radius_m,outer_radius_m,radial_scale_m,c0_m,c1_m,c2_m,c3_m,c4_m,c5_m,c6_m,c7_m,c8_m,c9_m,c10_m,c11_m,c12_m")
    return std::nullopt;
  std::optional<AxisymmetricMirror> primary, secondary, detector;
  SpectralResponse primary_reflectivity, secondary_reflectivity;
  while (read_line()) {
    const auto fields = split_scene_fields(line);
    if (fields.size() == 3 && (fields[0] == "primary_reflectivity" || fields[0] == "secondary_reflectivity")) {
      double wavelength = 0.0, response = 0.0;
      if (!scene_number(fields[1], wavelength) || !scene_number(fields[2], response)) return std::nullopt;
      auto& table = fields[0] == "primary_reflectivity" ? primary_reflectivity : secondary_reflectivity;
      table.wavelength_nm.push_back(wavelength);
      table.response.push_back(response);
      continue;
    }
    if (fields.size() != 18) return std::nullopt;
    double values[17]{};
    for (std::size_t index = 0; index < 17; ++index)
      if (!scene_number(fields[index + 1], values[index])) return std::nullopt;
    AxisymmetricMirror surface{};
    surface.vertex_z_m = values[0];
    surface.inner_radius_m = values[1];
    surface.outer_radius_m = values[2];
    surface.surface.radial_scale_m = values[3];
    for (std::size_t index = 0; index < surface.surface.coefficient_m.size(); ++index)
      surface.surface.coefficient_m[index] = values[index + 4];
    if (!is_valid(surface)) return std::nullopt;
    if (fields[0] == "primary" && !primary) primary = surface;
    else if (fields[0] == "secondary" && !secondary) secondary = surface;
    else if (fields[0] == "detector" && !detector) detector = surface;
    else return std::nullopt;
  }
  if (!primary || !secondary || !detector) return std::nullopt;
  const auto optional_response = [](SpectralResponse response) -> std::optional<SpectralResponse> {
    return response.wavelength_nm.empty() ? std::nullopt : std::optional{std::move(response)};
  };
  return AxisymmetricScene{std::move(provenance), *primary, *secondary, *detector,
                           optional_response(std::move(primary_reflectivity)),
                           optional_response(std::move(secondary_reflectivity))};
}

[[nodiscard]] inline PathRecord trace_axisymmetric_scene(const Ray& input, std::uint64_t photon_id,
                                                          const AxisymmetricScene& scene) {
  PathRecord record{};
  record.photon_id = photon_id;
  record.points_m[0] = input.position_m;
  record.point_count = 1;
  const auto direction = normalised_checked(input.direction);
  if (!direction) {
    record.status = PhotonStatus::invalid_input;
    return record;
  }
  Ray ray{input.position_m, *direction};
  const auto primary = intersect_axisymmetric_mirror(ray, scene.primary);
  if (!primary) {
    record.status = PhotonStatus::missed_primary;
    record.final_direction = ray.direction;
    return record;
  }
  record.points_m[1] = primary->point_m;
  record.point_count = 2;
  record.path_length_m = primary->distance_m;
  record.incidence_primary_deg = std::acos(std::clamp(std::abs(dot(ray.direction, primary->unit_normal)), 0.0, 1.0)) * 180.0 / std::numbers::pi;
  const auto after_primary = reflect(ray, *primary);
  if (!after_primary) { record.status = PhotonStatus::invalid_input; return record; }
  ray = *after_primary;
  const auto secondary = intersect_axisymmetric_mirror(ray, scene.secondary);
  if (!secondary) { record.status = PhotonStatus::missed_screen; record.final_direction = ray.direction; return record; }
  record.points_m[2] = secondary->point_m;
  record.point_count = 3;
  record.path_length_m += secondary->distance_m;
  record.incidence_secondary_deg = std::acos(std::clamp(std::abs(dot(ray.direction, secondary->unit_normal)), 0.0, 1.0)) * 180.0 / std::numbers::pi;
  const auto after_secondary = reflect(ray, *secondary);
  if (!after_secondary) { record.status = PhotonStatus::invalid_input; return record; }
  ray = *after_secondary;
  const auto detector = intersect_axisymmetric_mirror(ray, scene.detector);
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
