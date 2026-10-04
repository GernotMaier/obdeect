#pragma once

#include "obdeect/intersections.hpp"
#include "obdeect/math.hpp"
#include "obdeect/mirror_scatter.hpp"
#include "obdeect/model_import.hpp"
#include "obdeect/optical_response.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace obdeect {

// Observatory-neutral facet data after an external catalogue adapter has
// resolved all units and coordinate transforms. Shapes are retained exactly;
// no polygon is silently approximated as a circle by optical model compilation.
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
  // Positive radius of the panel's spherical optical surface. A zero value
  // represents a planar panel and is retained for generic plane optical_models.
  double curvature_radius_m{};
};

// A generic finite planar optical-arrival surface. Its placement and aperture
// are supplied by the adapter; optical model compilation never infers a camera shape.
// Detector IDs share the optical_model surface-ID namespace with facet IDs.
struct ImportedDetectorSurface {
  std::uint32_t id{};
  Vec3 centre_m{};
  Vec3 unit_normal{};
  double diameter_m{};
  FacetShape shape{FacetShape::circle};
  Vec3 unit_tangent_u{};
};

// An opaque finite cylinder supplied by a model adapter.  It is deliberately
// separate from optical surfaces because it has no reflective or transmissive
// branch: the first intersection terminates the photon as a component loss.
struct ImportedCylinderObscurer {
  std::uint32_t id{};
  Vec3 first_endpoint_m{};
  Vec3 second_endpoint_m{};
  double diameter_m{};
};

struct ImportedSegmentedOpticalModel {
  ModelProvenance provenance;
  std::vector<ImportedFacet> primary_facets;
  std::vector<ImportedDetectorSurface> detector_surfaces;
  std::vector<ImportedCylinderObscurer> cylinder_obscurers;
  std::optional<SpectralResponse> primary_reflectivity;
  std::optional<MirrorScatter> primary_scatter{};
  std::optional<CameraResponse> camera_response{};

  ImportedSegmentedOpticalModel(
      ModelProvenance imported_provenance, std::vector<ImportedFacet> imported_facets,
      std::vector<ImportedDetectorSurface> imported_detectors = {},
      std::vector<ImportedCylinderObscurer> imported_obscurers = {},
      std::optional<SpectralResponse> imported_reflectivity = std::nullopt)
      : provenance(std::move(imported_provenance)), primary_facets(std::move(imported_facets)),
        detector_surfaces(std::move(imported_detectors)),
        cylinder_obscurers(std::move(imported_obscurers)),
        primary_reflectivity(std::move(imported_reflectivity)) {}
};

class CompiledDetectorPlanes;

struct CompiledSegmentedOpticalModel {
  ModelProvenance provenance;
  std::vector<ImportedFacet> primary_facets;
  std::vector<ImportedDetectorSurface> detector_surfaces;
  std::vector<ImportedCylinderObscurer> cylinder_obscurers;
  std::optional<SpectralResponse> primary_reflectivity;
  std::optional<MirrorScatter> primary_scatter{};
  std::optional<CameraResponse> camera_response{};

  std::shared_ptr<const CompiledDetectorPlanes> detector_planes{};

  CompiledSegmentedOpticalModel(
      ModelProvenance compiled_provenance, std::vector<ImportedFacet> compiled_facets,
      std::vector<ImportedDetectorSurface> compiled_detectors = {},
      std::vector<ImportedCylinderObscurer> compiled_obscurers = {},
      std::optional<SpectralResponse> compiled_reflectivity = std::nullopt)
      : provenance(std::move(compiled_provenance)), primary_facets(std::move(compiled_facets)),
        detector_surfaces(std::move(compiled_detectors)),
        cylinder_obscurers(std::move(compiled_obscurers)),
        primary_reflectivity(std::move(compiled_reflectivity)) {}
};

[[nodiscard]] inline bool has_valid_provenance(const ModelProvenance &provenance) {
  if (provenance.model_name.empty() || provenance.model_version.empty() ||
      provenance.content_hash.size() != 64) {
    return false;
  }
  for (const char character : provenance.content_hash) {
    if (!((character >= '0' && character <= '9') || (character >= 'a' && character <= 'f')))
      return false;
  }
  return true;
}

[[nodiscard]] inline bool is_valid(const ImportedFacet &facet) {
  const auto normal = normalised_checked(facet.unit_normal);
  const auto tangent = normalised_checked(facet.unit_tangent_u);
  const bool orientation_is_valid =
      facet.shape == FacetShape::circle ||
      (normal.has_value() && tangent.has_value() && std::abs(dot(*normal, *tangent)) <= kEpsilon);
  return normal.has_value() && orientation_is_valid && std::isfinite(facet.centre_m.x) &&
         std::isfinite(facet.centre_m.y) && std::isfinite(facet.centre_m.z) &&
         std::isfinite(facet.diameter_m) && std::isfinite(facet.focal_length_m) &&
         facet.diameter_m > kEpsilon && facet.focal_length_m > kEpsilon &&
         std::isfinite(facet.curvature_radius_m) && facet.curvature_radius_m >= 0.0;
}

[[nodiscard]] inline bool is_valid(const ImportedDetectorSurface &surface) {
  const auto normal = normalised_checked(surface.unit_normal);
  const auto tangent = normalised_checked(surface.unit_tangent_u);
  const bool orientation_is_valid =
      surface.shape == FacetShape::circle ||
      (normal.has_value() && tangent.has_value() && std::abs(dot(*normal, *tangent)) <= kEpsilon);
  return normal.has_value() && orientation_is_valid && std::isfinite(surface.centre_m.x) &&
         std::isfinite(surface.centre_m.y) && std::isfinite(surface.centre_m.z) &&
         std::isfinite(surface.diameter_m) && surface.diameter_m > kEpsilon;
}

[[nodiscard]] inline bool is_valid(const ImportedCylinderObscurer &obscurer) {
  return std::isfinite(obscurer.first_endpoint_m.x) && std::isfinite(obscurer.first_endpoint_m.y) &&
         std::isfinite(obscurer.first_endpoint_m.z) &&
         std::isfinite(obscurer.second_endpoint_m.x) &&
         std::isfinite(obscurer.second_endpoint_m.y) &&
         std::isfinite(obscurer.second_endpoint_m.z) && std::isfinite(obscurer.diameter_m) &&
         obscurer.diameter_m > kEpsilon &&
         norm(obscurer.second_endpoint_m - obscurer.first_endpoint_m) > kEpsilon;
}

[[nodiscard]] inline bool is_valid(const CompiledSegmentedOpticalModel &optical_model) {
  if (!has_valid_provenance(optical_model.provenance) || optical_model.primary_facets.empty()) {
    return false;
  }
  std::unordered_set<std::uint32_t> ids;
  ids.reserve(optical_model.primary_facets.size());
  for (const auto &facet : optical_model.primary_facets) {
    if (!is_valid(facet) || !ids.insert(facet.id).second)
      return false;
  }
  for (const auto &detector : optical_model.detector_surfaces) {
    if (!is_valid(detector) || !ids.insert(detector.id).second)
      return false;
  }
  for (const auto &obscurer : optical_model.cylinder_obscurers) {
    if (!is_valid(obscurer) || !ids.insert(obscurer.id).second)
      return false;
  }
  if ((optical_model.primary_reflectivity && !optical_model.primary_reflectivity->is_valid()) ||
      (optical_model.primary_scatter && !optical_model.primary_scatter->is_valid()) ||
      (optical_model.camera_response && !optical_model.camera_response->is_valid()))
    return false;
  return true;
}

// The native compiler is deliberately strict: an adapter must provide a
// normal for every facet rather than allowing the tracing kernel to infer a
// telescope-specific dish prescription.
[[nodiscard]] inline std::optional<CompiledSegmentedOpticalModel>
compile_segmented_optical_model(const ImportedSegmentedOpticalModel &input) {
  if (!has_valid_provenance(input.provenance) || input.primary_facets.empty()) {
    return std::nullopt;
  }
  std::unordered_set<std::uint32_t> ids;
  ids.reserve(input.primary_facets.size());
  for (const auto &facet : input.primary_facets) {
    if (!is_valid(facet) || !ids.insert(facet.id).second)
      return std::nullopt;
  }
  for (const auto &detector : input.detector_surfaces) {
    if (!is_valid(detector) || !ids.insert(detector.id).second)
      return std::nullopt;
  }
  for (const auto &obscurer : input.cylinder_obscurers) {
    if (!is_valid(obscurer) || !ids.insert(obscurer.id).second)
      return std::nullopt;
  }
  if ((input.primary_reflectivity && !input.primary_reflectivity->is_valid()) ||
      (input.primary_scatter && !input.primary_scatter->is_valid()) ||
      (input.camera_response && !input.camera_response->is_valid()))
    return std::nullopt;
  auto result =
      CompiledSegmentedOpticalModel{input.provenance, input.primary_facets, input.detector_surfaces,
                                    input.cylinder_obscurers, input.primary_reflectivity};
  result.primary_scatter = input.primary_scatter;
  result.camera_response = input.camera_response;
  return result;
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

[[nodiscard]] inline bool contains_planar_aperture_point(const Vec3 &centre_m,
                                                         const Vec3 &unit_normal,
                                                         const Vec3 &unit_tangent_u,
                                                         double diameter_m, FacetShape shape,
                                                         const Vec3 &point_m) {
  const auto normal = normalised_checked(unit_normal);
  if (!normal)
    return false;
  const Vec3 displacement = point_m - centre_m;
  if (std::abs(dot(displacement, *normal)) > kEpsilon)
    return false;
  if (shape == FacetShape::circle) {
    const double radius_m = diameter_m * 0.5 + kEpsilon;
    return dot(displacement, displacement) <= radius_m * radius_m;
  }

  const auto tangent_u = normalised_checked(unit_tangent_u);
  if (!tangent_u || std::abs(dot(*normal, *tangent_u)) > kEpsilon)
    return false;
  const auto tangent_v = normalised_checked(cross(*normal, *tangent_u));
  if (!tangent_v)
    return false;
  const double u = dot(displacement, *tangent_u);
  const double v = dot(displacement, *tangent_v);
  const double apothem = diameter_m * 0.5;
  constexpr double sqrt_three = 1.7320508075688772935;
  switch (shape) {
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
  case FacetShape::circle:
    break;
  }
  return false;
}

[[nodiscard]] inline bool contains_facet_point(const ImportedFacet &facet, const Vec3 &point_m) {
  return contains_planar_aperture_point(facet.centre_m, facet.unit_normal, facet.unit_tangent_u,
                                        facet.diameter_m, facet.shape, point_m);
}

[[nodiscard]] inline bool contains_detector_point(const ImportedDetectorSurface &surface,
                                                  const Vec3 &point_m) {
  return contains_planar_aperture_point(surface.centre_m, surface.unit_normal,
                                        surface.unit_tangent_u, surface.diameter_m, surface.shape,
                                        point_m);
}

[[nodiscard]] inline std::optional<SegmentedFacetHit>
intersect_segmented_facet_unchecked(const Ray &ray, const ImportedFacet &facet,
                                    double minimum_t_m = kEpsilon) {
  if (!std::isfinite(minimum_t_m) || minimum_t_m < 0.0)
    return std::nullopt;
  const auto direction = normalised_checked(ray.direction);
  const auto normal = normalised_checked(facet.unit_normal);
  if (!direction || !normal)
    return std::nullopt;
  if (facet.curvature_radius_m == 0.0) {
    const double denominator = dot(*direction, *normal);
    if (std::abs(denominator) <= kEpsilon)
      return std::nullopt;
    const double distance_m = dot(facet.centre_m - ray.position_m, *normal) / denominator;
    if (!std::isfinite(distance_m) || distance_m <= minimum_t_m)
      return std::nullopt;
    const Vec3 point_m = ray.position_m + *direction * distance_m;
    if (!contains_facet_point(facet, point_m))
      return std::nullopt;
    return SegmentedFacetHit{facet.id, distance_m, point_m, *normal};
  }
  // The optical face is concave towards the incoming beam.  With a nominal
  // normal pointing from the primary towards the focal surface, the sphere
  // centre lies on that focal-surface side of the facet.  Putting it behind
  // the facet produces a convex mirror and sends off-centre rays outward.
  const Vec3 sphere_centre = facet.centre_m + *normal * facet.curvature_radius_m;
  const Vec3 offset = ray.position_m - sphere_centre;
  const double projection = dot(offset, *direction);
  const double discriminant =
      projection * projection -
      (dot(offset, offset) - facet.curvature_radius_m * facet.curvature_radius_m);
  if (!std::isfinite(discriminant) || discriminant < 0.0)
    return std::nullopt;
  const double root = std::sqrt(std::max(0.0, discriminant));
  for (const double distance_m : {-projection - root, -projection + root}) {
    if (!std::isfinite(distance_m) || distance_m <= minimum_t_m)
      continue;
    const Vec3 point_m = ray.position_m + *direction * distance_m;
    const auto surface_normal = normalised_checked(point_m - sphere_centre);
    // The geometric sphere normal points away from its centre and therefore
    // opposite the nominal facet normal on the illuminated concave face.
    if (!surface_normal || dot(*surface_normal, *normal) >= 0.0)
      continue;
    const Vec3 tangent_displacement =
        point_m - facet.centre_m - *normal * dot(point_m - facet.centre_m, *normal);
    const Vec3 aperture_point = facet.centre_m + tangent_displacement;
    if (contains_facet_point(facet, aperture_point))
      return SegmentedFacetHit{facet.id, distance_m, point_m, *surface_normal};
  }
  return std::nullopt;
}

[[nodiscard]] inline std::optional<SegmentedFacetHit>
intersect_segmented_facet(const Ray &ray, const ImportedFacet &facet,
                          double minimum_t_m = kEpsilon) {
  if (!is_valid(facet))
    return std::nullopt;
  return intersect_segmented_facet_unchecked(ray, facet, minimum_t_m);
}

[[nodiscard]] inline std::optional<SegmentedFacetHit>
intersect_segmented_primary_unchecked(const Ray &ray,
                                      const CompiledSegmentedOpticalModel &optical_model) {
  std::optional<SegmentedFacetHit> nearest;
  for (const auto &facet : optical_model.primary_facets) {
    const auto candidate = intersect_segmented_facet_unchecked(ray, facet);
    if (candidate &&
        (!nearest || candidate->distance_m < nearest->distance_m ||
         (candidate->distance_m == nearest->distance_m && candidate->facet_id < nearest->facet_id)))
      nearest = candidate;
  }
  return nearest;
}

[[nodiscard]] inline std::optional<SegmentedFacetHit>
intersect_segmented_primary(const Ray &ray, const CompiledSegmentedOpticalModel &optical_model) {
  if (!is_valid(optical_model))
    return std::nullopt;
  return intersect_segmented_primary_unchecked(ray, optical_model);
}

struct DetectorSurfaceHit {
  std::uint32_t surface_id{};
  double distance_m{};
  Vec3 point_m{};
  Vec3 unit_normal{};
};

struct CylinderObscurerHit {
  std::uint32_t surface_id{};
  double distance_m{};
};

[[nodiscard]] inline std::optional<CylinderObscurerHit>
intersect_cylinder_obscurers_unchecked(const Ray &ray,
                                       const CompiledSegmentedOpticalModel &optical_model) {
  std::optional<CylinderObscurerHit> nearest;
  for (const auto &obscurer : optical_model.cylinder_obscurers) {
    const auto distance_m = intersect_closed_finite_cylinder(
        ray, obscurer.first_endpoint_m, obscurer.second_endpoint_m, obscurer.diameter_m * 0.5);
    if (distance_m && (!nearest || *distance_m < nearest->distance_m ||
                       (*distance_m == nearest->distance_m && obscurer.id < nearest->surface_id)))
      nearest = CylinderObscurerHit{obscurer.id, *distance_m};
  }
  return nearest;
}

[[nodiscard]] inline std::optional<DetectorSurfaceHit>
intersect_detector_surface_unchecked(const Ray &ray, const ImportedDetectorSurface &surface,
                                     double minimum_t_m = kEpsilon) {
  if (!std::isfinite(minimum_t_m) || minimum_t_m < 0.0)
    return std::nullopt;
  const auto direction = normalised_checked(ray.direction);
  const auto normal = normalised_checked(surface.unit_normal);
  if (!direction || !normal)
    return std::nullopt;
  const double denominator = dot(*direction, *normal);
  if (std::abs(denominator) <= kEpsilon)
    return std::nullopt;
  const double distance_m = dot(surface.centre_m - ray.position_m, *normal) / denominator;
  if (!std::isfinite(distance_m) || distance_m <= minimum_t_m)
    return std::nullopt;
  const Vec3 point_m = ray.position_m + *direction * distance_m;
  if (!contains_detector_point(surface, point_m))
    return std::nullopt;
  return DetectorSurfaceHit{surface.id, distance_m, point_m, *normal};
}

[[nodiscard]] inline std::optional<DetectorSurfaceHit>
intersect_detector_surface(const Ray &ray, const ImportedDetectorSurface &surface,
                           double minimum_t_m = kEpsilon) {
  if (!is_valid(surface))
    return std::nullopt;
  return intersect_detector_surface_unchecked(ray, surface, minimum_t_m);
}

[[nodiscard]] inline std::optional<DetectorSurfaceHit>
intersect_detector_surfaces_unchecked(const Ray &ray,
                                      const CompiledSegmentedOpticalModel &optical_model) {
  std::optional<DetectorSurfaceHit> nearest;
  for (const auto &surface : optical_model.detector_surfaces) {
    const auto candidate = intersect_detector_surface_unchecked(ray, surface);
    if (candidate && (!nearest || candidate->distance_m < nearest->distance_m ||
                      (candidate->distance_m == nearest->distance_m &&
                       candidate->surface_id < nearest->surface_id)))
      nearest = candidate;
  }
  return nearest;
}

[[nodiscard]] inline std::optional<DetectorSurfaceHit>
intersect_detector_surfaces(const Ray &ray, const CompiledSegmentedOpticalModel &optical_model) {
  if (!is_valid(optical_model))
    return std::nullopt;
  return intersect_detector_surfaces_unchecked(ray, optical_model);
}

} // namespace obdeect
