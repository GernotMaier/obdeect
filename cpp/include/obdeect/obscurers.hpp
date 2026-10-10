#pragma once

#include "obdeect/math.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>

namespace obdeect {

struct OpaqueSurface {
  enum class Shape { quadrilateral, hollow_frustum, solid_frustum };
  std::uint32_t id{};
  Shape shape{Shape::quadrilateral};
  std::array<Vec3, 4> vertices{};
  Vec3 first{}, second{};
  double first_radius_m{}, second_radius_m{}, thickness_m{};
};
struct OpaqueHit {
  std::uint32_t surface_id{};
  double distance_m{};
  Vec3 point_m{}, unit_normal{};
};

[[nodiscard]] inline bool is_valid(const OpaqueSurface &surface) {
  const auto finite = [](const Vec3 &p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
  };
  if (surface.id == std::numeric_limits<std::uint32_t>::max())
    return false;
  if (surface.shape != OpaqueSurface::Shape::quadrilateral &&
      surface.shape != OpaqueSurface::Shape::hollow_frustum &&
      surface.shape != OpaqueSurface::Shape::solid_frustum)
    return false;
  if (surface.shape != OpaqueSurface::Shape::quadrilateral)
    return finite(surface.first) && finite(surface.second) &&
           norm(surface.second - surface.first) > kEpsilon &&
           std::isfinite(surface.first_radius_m) && std::isfinite(surface.second_radius_m) &&
           std::isfinite(surface.thickness_m) && surface.first_radius_m > 0 &&
           surface.second_radius_m > 0 && surface.thickness_m >= 0 &&
           (surface.shape != OpaqueSurface::Shape::solid_frustum || surface.thickness_m == 0);
  for (const auto &p : surface.vertices)
    if (!finite(p))
      return false;
  const auto normal = normalised_checked(
      cross(surface.vertices[1] - surface.vertices[0], surface.vertices[2] - surface.vertices[0]));
  if (!normal)
    return false;
  double sign = 0;
  for (std::size_t i = 0; i < 4; ++i) {
    if (std::abs(dot(surface.vertices[i] - surface.vertices[0], *normal)) > kEpsilon)
      return false;
    const auto a = surface.vertices[(i + 1) % 4] - surface.vertices[i];
    const auto b = surface.vertices[(i + 2) % 4] - surface.vertices[(i + 1) % 4];
    const double turn = dot(cross(a, b), *normal);
    if (std::abs(turn) <= kEpsilon || (sign && turn * sign < 0))
      return false;
    sign = turn;
  }
  return true;
}

// Finite convex plate or an open hollow cylinder/cone, including annular ends.
// All storage and geometry validation happen before the photon kernel.
[[nodiscard]] inline std::optional<OpaqueHit>
intersect_opaque_surface(const Ray &ray, const OpaqueSurface &surface) {
  std::optional<OpaqueHit> nearest;
  const auto consider = [&](double t, Vec3 normal) {
    if (std::isfinite(t) && t > kEpsilon && (!nearest || t < nearest->distance_m))
      nearest = OpaqueHit{surface.id, t, ray.position_m + ray.direction * t, normal};
  };
  if (surface.shape == OpaqueSurface::Shape::quadrilateral) {
    const Vec3 normal = *normalised_checked(cross(surface.vertices[1] - surface.vertices[0],
                                                  surface.vertices[2] - surface.vertices[0]));
    const double denominator = dot(normal, ray.direction);
    if (std::abs(denominator) <= kEpsilon)
      return {};
    const double t = dot(surface.vertices[0] - ray.position_m, normal) / denominator;
    const Vec3 point = ray.position_m + ray.direction * t;
    for (std::size_t i = 0; i < 4; ++i)
      if (dot(cross(surface.vertices[(i + 1) % 4] - surface.vertices[i],
                    point - surface.vertices[i]),
              normal) < -kEpsilon)
        return {};
    consider(t, normal);
    return nearest;
  }
  const Vec3 delta = surface.second - surface.first;
  const double length = norm(delta);
  const Vec3 axis = delta * (1.0 / length), offset = ray.position_m - surface.first;
  const double position_axis = dot(offset, axis), direction_axis = dot(ray.direction, axis);
  const Vec3 radial_position = offset - axis * position_axis,
             radial_direction = ray.direction - axis * direction_axis;
  const double slope = (surface.second_radius_m - surface.first_radius_m) / length;
  for (const double thickness : {0.0, surface.thickness_m}) {
    const double radius = surface.first_radius_m + thickness + slope * position_axis;
    const double a =
        dot(radial_direction, radial_direction) - slope * slope * direction_axis * direction_axis;
    const double b = dot(radial_direction, radial_position) - radius * slope * direction_axis;
    const double c = dot(radial_position, radial_position) - radius * radius;
    const auto root = [&](double t) {
      const double along = position_axis + t * direction_axis;
      if (along < 0 || along > length)
        return;
      const Vec3 radial = radial_position + radial_direction * t;
      const auto normal = normalised_checked(
          radial - axis * (slope * (surface.first_radius_m + thickness + slope * along)));
      if (normal)
        consider(t, thickness == 0 && surface.shape != OpaqueSurface::Shape::solid_frustum
                        ? *normal * -1
                        : *normal);
    };
    if (std::abs(a) > kEpsilon) {
      const double discriminant = b * b - a * c;
      if (discriminant >= 0) {
        const double square = std::sqrt(discriminant);
        root((-b - square) / a);
        root((-b + square) / a);
      }
    } else if (std::abs(b) > kEpsilon)
      root(-c / (2 * b));
    if (surface.thickness_m == 0)
      break;
  }
  if ((surface.thickness_m > 0 || surface.shape == OpaqueSurface::Shape::solid_frustum) &&
      std::abs(direction_axis) > kEpsilon)
    for (double along : {0.0, length}) {
      const double t = (along - position_axis) / direction_axis;
      const Vec3 radial = radial_position + radial_direction * t;
      const double radius = surface.first_radius_m + slope * along, r2 = dot(radial, radial);
      if ((surface.shape == OpaqueSurface::Shape::solid_frustum || r2 >= radius * radius) &&
          r2 <= (radius + surface.thickness_m) * (radius + surface.thickness_m))
        consider(t, axis * (along == 0 ? -1 : 1));
    }
  return nearest;
}

[[nodiscard]] inline std::optional<OpaqueHit>
intersect_opaque_surfaces(const Ray &ray, std::span<const OpaqueSurface> surfaces) {
  std::optional<OpaqueHit> nearest;
  for (const auto &surface : surfaces) {
    const auto hit = intersect_opaque_surface(ray, surface);
    if (hit && (!nearest || hit->distance_m < nearest->distance_m ||
                (hit->distance_m == nearest->distance_m && hit->surface_id < nearest->surface_id)))
      nearest = hit;
  }
  return nearest;
}
} // namespace obdeect
