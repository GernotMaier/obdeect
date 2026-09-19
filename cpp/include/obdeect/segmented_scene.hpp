#pragma once

#include "obdeect/math.hpp"
#include "obdeect/model_import.hpp"

#include <cstdint>
#include <cmath>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace obdeect {

// Observatory-neutral facet data after an external catalogue adapter has
// resolved all units and coordinate transforms. Shapes are retained exactly;
// no polygon is silently approximated as a circle by scene compilation.
enum class FacetShape : std::uint8_t { circle, hexagon_flat_y, square, hexagon_flat_x };

struct ImportedFacet {
  std::uint32_t id{};
  Vec3 centre_m{};
  Vec3 unit_normal{};
  double diameter_m{};
  double focal_length_m{};
  FacetShape shape{FacetShape::circle};
  // The in-plane orientation is observable for square and hexagonal panels.
  // It is supplied by the model adapter; tracing must not manufacture a
  // telescope-specific panel rotation from a dish prescription.
  Vec3 unit_tangent_u{};
};

// A generic finite planar optical-arrival surface. Its placement and aperture
// are supplied by the adapter; scene compilation never infers a camera shape.
// Detector IDs share the scene surface-ID namespace with facet IDs.
struct ImportedDetectorSurface {
  std::uint32_t id{};
  Vec3 centre_m{};
  Vec3 unit_normal{};
  double diameter_m{};
  FacetShape shape{FacetShape::circle};
  Vec3 unit_tangent_u{};
};

struct ImportedSegmentedScene {
  ModelProvenance provenance;
  std::vector<ImportedFacet> primary_facets;
  std::vector<ImportedDetectorSurface> detector_surfaces;

  ImportedSegmentedScene(ModelProvenance imported_provenance, std::vector<ImportedFacet> imported_facets,
                         std::vector<ImportedDetectorSurface> imported_detectors = {})
      : provenance(std::move(imported_provenance)), primary_facets(std::move(imported_facets)),
        detector_surfaces(std::move(imported_detectors)) {}
};

struct CompiledSegmentedScene {
  ModelProvenance provenance;
  std::vector<ImportedFacet> primary_facets;
  std::vector<ImportedDetectorSurface> detector_surfaces;

  CompiledSegmentedScene(ModelProvenance compiled_provenance, std::vector<ImportedFacet> compiled_facets,
                         std::vector<ImportedDetectorSurface> compiled_detectors = {})
      : provenance(std::move(compiled_provenance)), primary_facets(std::move(compiled_facets)),
        detector_surfaces(std::move(compiled_detectors)) {}
};

[[nodiscard]] inline bool has_valid_provenance(const ModelProvenance& provenance) {
  if (provenance.model_name.empty() || provenance.model_version.empty() || provenance.content_hash.size() != 64) {
    return false;
  }
  for (const char character : provenance.content_hash) {
    if (!((character >= '0' && character <= '9') || (character >= 'a' && character <= 'f'))) return false;
  }
  return true;
}

[[nodiscard]] inline bool is_valid(const ImportedFacet& facet) {
  const auto normal = normalised_checked(facet.unit_normal);
  const auto tangent = normalised_checked(facet.unit_tangent_u);
  const bool orientation_is_valid = facet.shape == FacetShape::circle ||
                                    (normal.has_value() && tangent.has_value() &&
                                     std::abs(dot(*normal, *tangent)) <= kEpsilon);
  return normal.has_value() && orientation_is_valid && std::isfinite(facet.centre_m.x) &&
         std::isfinite(facet.centre_m.y) &&
         std::isfinite(facet.centre_m.z) && std::isfinite(facet.diameter_m) &&
         std::isfinite(facet.focal_length_m) && facet.diameter_m > kEpsilon &&
         facet.focal_length_m > kEpsilon;
}

[[nodiscard]] inline bool is_valid(const ImportedDetectorSurface& surface) {
  const auto normal = normalised_checked(surface.unit_normal);
  const auto tangent = normalised_checked(surface.unit_tangent_u);
  const bool orientation_is_valid = surface.shape == FacetShape::circle ||
                                    (normal.has_value() && tangent.has_value() &&
                                     std::abs(dot(*normal, *tangent)) <= kEpsilon);
  return normal.has_value() && orientation_is_valid && std::isfinite(surface.centre_m.x) &&
         std::isfinite(surface.centre_m.y) && std::isfinite(surface.centre_m.z) &&
         std::isfinite(surface.diameter_m) && surface.diameter_m > kEpsilon;
}

[[nodiscard]] inline bool is_valid(const CompiledSegmentedScene& scene) {
  if (!has_valid_provenance(scene.provenance) || scene.primary_facets.empty()) {
    return false;
  }
  std::unordered_set<std::uint32_t> ids;
  ids.reserve(scene.primary_facets.size());
  for (const auto& facet : scene.primary_facets) {
    if (!is_valid(facet) || !ids.insert(facet.id).second) return false;
  }
  for (const auto& detector : scene.detector_surfaces) {
    if (!is_valid(detector) || !ids.insert(detector.id).second) return false;
  }
  return true;
}

// The native compiler is deliberately strict: an adapter must provide a
// normal for every facet rather than allowing the tracing kernel to infer a
// telescope-specific dish prescription.
[[nodiscard]] inline std::optional<CompiledSegmentedScene> compile_segmented_scene(
    const ImportedSegmentedScene& input) {
  if (!has_valid_provenance(input.provenance) || input.primary_facets.empty()) {
    return std::nullopt;
  }
  std::unordered_set<std::uint32_t> ids;
  ids.reserve(input.primary_facets.size());
  for (const auto& facet : input.primary_facets) {
    if (!is_valid(facet) || !ids.insert(facet.id).second) return std::nullopt;
  }
  for (const auto& detector : input.detector_surfaces) {
    if (!is_valid(detector) || !ids.insert(detector.id).second) return std::nullopt;
  }
  return CompiledSegmentedScene{input.provenance, input.primary_facets, input.detector_surfaces};
}

// A planar, finite panel hit.  The mirror-list diameter convention is kept
// explicit: circles use it as their diameter; squares as their side length;
// regular hexagons as their flat-to-flat distance.
struct SegmentedFacetHit {
  std::uint32_t facet_id{};
  double distance_m{};
  Vec3 point_m{};
  Vec3 unit_normal{};
};

[[nodiscard]] inline bool contains_facet_point(const ImportedFacet& facet, const Vec3& point_m) {
  const auto normal = normalised_checked(facet.unit_normal);
  if (!normal) return false;
  const Vec3 displacement = point_m - facet.centre_m;
  if (std::abs(dot(displacement, *normal)) > kEpsilon) return false;
  if (facet.shape == FacetShape::circle) return norm(displacement) <= facet.diameter_m * 0.5 + kEpsilon;

  const auto tangent_u = normalised_checked(facet.unit_tangent_u);
  if (!tangent_u || std::abs(dot(*normal, *tangent_u)) > kEpsilon) return false;
  const auto tangent_v = normalised_checked(cross(*normal, *tangent_u));
  if (!tangent_v) return false;
  const double u = dot(displacement, *tangent_u);
  const double v = dot(displacement, *tangent_v);
  const double apothem = facet.diameter_m * 0.5;
  constexpr double sqrt_three = 1.7320508075688772935;
  switch (facet.shape) {
    case FacetShape::square:
      return std::abs(u) <= apothem + kEpsilon && std::abs(v) <= apothem + kEpsilon;
    case FacetShape::hexagon_flat_y:
      return std::abs(v) <= apothem + kEpsilon &&
             std::abs(sqrt_three * u + v) <= 2.0 * apothem + kEpsilon &&
             std::abs(sqrt_three * u - v) <= 2.0 * apothem + kEpsilon;
    case FacetShape::hexagon_flat_x:
      return std::abs(u) <= apothem + kEpsilon &&
             std::abs(u + sqrt_three * v) <= 2.0 * apothem + kEpsilon &&
             std::abs(u - sqrt_three * v) <= 2.0 * apothem + kEpsilon;
    case FacetShape::circle: break;
  }
  return false;
}

[[nodiscard]] inline bool contains_detector_point(const ImportedDetectorSurface& surface, const Vec3& point_m) {
  const auto normal = normalised_checked(surface.unit_normal);
  if (!normal) return false;
  const Vec3 displacement = point_m - surface.centre_m;
  if (std::abs(dot(displacement, *normal)) > kEpsilon) return false;
  if (surface.shape == FacetShape::circle) return norm(displacement) <= surface.diameter_m * 0.5 + kEpsilon;

  const auto tangent_u = normalised_checked(surface.unit_tangent_u);
  if (!tangent_u || std::abs(dot(*normal, *tangent_u)) > kEpsilon) return false;
  const auto tangent_v = normalised_checked(cross(*normal, *tangent_u));
  if (!tangent_v) return false;
  const double u = dot(displacement, *tangent_u);
  const double v = dot(displacement, *tangent_v);
  const double apothem = surface.diameter_m * 0.5;
  constexpr double sqrt_three = 1.7320508075688772935;
  switch (surface.shape) {
    case FacetShape::square:
      return std::abs(u) <= apothem + kEpsilon && std::abs(v) <= apothem + kEpsilon;
    case FacetShape::hexagon_flat_y:
      return std::abs(v) <= apothem + kEpsilon &&
             std::abs(sqrt_three * u + v) <= 2.0 * apothem + kEpsilon &&
             std::abs(sqrt_three * u - v) <= 2.0 * apothem + kEpsilon;
    case FacetShape::hexagon_flat_x:
      return std::abs(u) <= apothem + kEpsilon &&
             std::abs(u + sqrt_three * v) <= 2.0 * apothem + kEpsilon &&
             std::abs(u - sqrt_three * v) <= 2.0 * apothem + kEpsilon;
    case FacetShape::circle: break;
  }
  return false;
}

[[nodiscard]] inline std::optional<SegmentedFacetHit> intersect_segmented_facet(
    const Ray& ray, const ImportedFacet& facet, double minimum_t_m = kEpsilon) {
  if (!is_valid(facet) || !std::isfinite(minimum_t_m) || minimum_t_m < 0.0) return std::nullopt;
  const auto direction = normalised_checked(ray.direction);
  const auto normal = normalised_checked(facet.unit_normal);
  if (!direction || !normal) return std::nullopt;
  const double denominator = dot(*direction, *normal);
  if (std::abs(denominator) <= kEpsilon) return std::nullopt;
  const double distance_m = dot(facet.centre_m - ray.position_m, *normal) / denominator;
  if (!std::isfinite(distance_m) || distance_m <= minimum_t_m) return std::nullopt;
  const Vec3 point_m = ray.position_m + *direction * distance_m;
  if (!contains_facet_point(facet, point_m)) return std::nullopt;
  return SegmentedFacetHit{facet.id, distance_m, point_m, *normal};
}

[[nodiscard]] inline std::optional<SegmentedFacetHit> intersect_segmented_primary(
    const Ray& ray, const CompiledSegmentedScene& scene) {
  if (!is_valid(scene)) return std::nullopt;
  std::optional<SegmentedFacetHit> nearest;
  for (const auto& facet : scene.primary_facets) {
    const auto candidate = intersect_segmented_facet(ray, facet);
    if (candidate && (!nearest || candidate->distance_m < nearest->distance_m)) nearest = candidate;
  }
  return nearest;
}

struct DetectorSurfaceHit {
  std::uint32_t surface_id{};
  double distance_m{};
  Vec3 point_m{};
  Vec3 unit_normal{};
};

[[nodiscard]] inline std::optional<DetectorSurfaceHit> intersect_detector_surface(
    const Ray& ray, const ImportedDetectorSurface& surface, double minimum_t_m = kEpsilon) {
  if (!is_valid(surface) || !std::isfinite(minimum_t_m) || minimum_t_m < 0.0) return std::nullopt;
  const auto direction = normalised_checked(ray.direction);
  const auto normal = normalised_checked(surface.unit_normal);
  if (!direction || !normal) return std::nullopt;
  const double denominator = dot(*direction, *normal);
  if (std::abs(denominator) <= kEpsilon) return std::nullopt;
  const double distance_m = dot(surface.centre_m - ray.position_m, *normal) / denominator;
  if (!std::isfinite(distance_m) || distance_m <= minimum_t_m) return std::nullopt;
  const Vec3 point_m = ray.position_m + *direction * distance_m;
  if (!contains_detector_point(surface, point_m)) return std::nullopt;
  return DetectorSurfaceHit{surface.id, distance_m, point_m, *normal};
}

[[nodiscard]] inline std::optional<DetectorSurfaceHit> intersect_detector_surfaces(
    const Ray& ray, const CompiledSegmentedScene& scene) {
  if (!is_valid(scene)) return std::nullopt;
  std::optional<DetectorSurfaceHit> nearest;
  for (const auto& surface : scene.detector_surfaces) {
    const auto candidate = intersect_detector_surface(ray, surface);
    if (candidate && (!nearest || candidate->distance_m < nearest->distance_m)) nearest = candidate;
  }
  return nearest;
}

}  // namespace obdeect
