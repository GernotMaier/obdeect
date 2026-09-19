#pragma once

#include "obdeect/interactions.hpp"
#include "obdeect/tables.hpp"

#include <limits>
#include <cstdint>
#include <vector>

namespace obdeect {

// A compact physical facet representation used after CTAO model import.
// The facet is deliberately a finite circular plane rather than a global dish
// approximation: each imported centre and normal remains observable.
struct CircularFacet {
  std::uint32_t id{};
  Vec3 centre_m{};
  Vec3 unit_normal{};
  double radius_m{};
  Table1DView reflectivity{};
};

struct FacetHit {
  std::uint32_t facet_id{};
  double distance_m{};
  Vec3 point_m{};
  Vec3 unit_normal{};
  double reflectivity{};
};

[[nodiscard]] inline bool is_valid(const CircularFacet& facet) {
  for (const double value : facet.reflectivity.value) {
    if (value > 1.0) return false;
  }
  const auto normal = normalised_checked(facet.unit_normal);
  return normal.has_value() && std::isfinite(facet.centre_m.x) && std::isfinite(facet.centre_m.y) &&
         std::isfinite(facet.centre_m.z) && std::isfinite(facet.radius_m) && facet.radius_m > kEpsilon &&
         facet.reflectivity.is_valid();
}

[[nodiscard]] inline std::optional<FacetHit> intersect_facet(const Ray& ray, const CircularFacet& facet,
                                                              double wavelength_nm,
                                                              double minimum_t_m = kEpsilon) {
  if (!is_valid(facet) || !std::isfinite(wavelength_nm) || !std::isfinite(minimum_t_m) ||
      minimum_t_m < 0.0) return std::nullopt;
  const auto direction = normalised_checked(ray.direction);
  const auto normal = normalised_checked(facet.unit_normal);
  const auto reflectivity = facet.reflectivity.interpolate(wavelength_nm);
  if (!direction || !normal || !reflectivity) return std::nullopt;
  const double denominator = dot(*direction, *normal);
  if (std::abs(denominator) <= kEpsilon) return std::nullopt;
  const double distance_m = dot(facet.centre_m - ray.position_m, *normal) / denominator;
  if (!std::isfinite(distance_m) || distance_m <= minimum_t_m) return std::nullopt;
  const Vec3 point_m = ray.position_m + *direction * distance_m;
  if (norm(point_m - facet.centre_m) > facet.radius_m) return std::nullopt;
  return FacetHit{facet.id, distance_m, point_m, *normal, *reflectivity};
}

struct FacetedMirror {
  std::vector<CircularFacet> facets;
};

[[nodiscard]] inline std::optional<FacetHit> intersect_faceted_mirror(const Ray& ray,
                                                                        const FacetedMirror& mirror,
                                                                        double wavelength_nm) {
  std::optional<FacetHit> nearest;
  for (const auto& facet : mirror.facets) {
    const auto candidate = intersect_facet(ray, facet, wavelength_nm);
    if (candidate && (!nearest || candidate->distance_m < nearest->distance_m)) nearest = candidate;
  }
  return nearest;
}

[[nodiscard]] inline std::optional<Ray> reflect(const Ray& incident, const FacetHit& hit) {
  const auto direction = reflect_specular(incident.direction, hit.unit_normal);
  return direction ? std::optional<Ray>{{hit.point_m, *direction}} : std::nullopt;
}

}  // namespace obdeect
