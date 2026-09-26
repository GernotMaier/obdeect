#pragma once

#include "obdeect/segmented_scene.hpp"

#include <fstream>
#include <cmath>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace obdeect {

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
  if (!read_line() ||
      line != "surface_id,role,shape,cx_m,cy_m,cz_m,nx,ny,nz,tx,ty,tz,diameter_m")
    return std::nullopt;

  std::vector<ImportedFacet> facets;
  std::vector<ImportedDetectorSurface> detectors;
  while (read_line()) {
    if (line.empty()) continue;
    const auto fields = split_scene_fields(line);
    if (fields.size() != 13) return std::nullopt;
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
      facets.push_back({id, centre, normal, values[9], values[9], *shape, tangent});
    } else if (fields[1] == "detector") {
      detectors.push_back({id, centre, normal, values[9], *shape, tangent});
    } else {
      return std::nullopt;
    }
  }
  if (facets.empty() || detectors.empty()) return std::nullopt;
  return compile_segmented_scene(ImportedSegmentedScene{std::move(provenance), std::move(facets),
                                                         std::move(detectors)});
}

}  // namespace obdeect
