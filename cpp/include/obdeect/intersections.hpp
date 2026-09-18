#pragma once

#include "obdeect/math.hpp"

#include <algorithm>
#include <optional>

namespace obdeect {

inline std::optional<double> intersect_sphere(const Ray& ray, const Vec3& center, double radius) {
  if (!std::isfinite(radius) || radius <= kEpsilon) return std::nullopt;
  const Vec3 oc = ray.position_m - center;
  const double b = dot(oc, ray.direction);
  const double c = dot(oc, oc) - radius * radius;
  const double discriminant = b * b - c;
  if (discriminant < 0.0) return std::nullopt;
  const double root = std::sqrt(std::max(0.0, discriminant));
  for (const double distance_m : {-b - root, -b + root}) {
    if (distance_m > kEpsilon) return distance_m;
  }
  return std::nullopt;
}

inline std::optional<double> intersect_lower_spherical_cap(const Ray& ray, const Vec3& center,
                                                           double radius) {
  const auto hit = intersect_sphere(ray, center, radius);
  if (!hit) return std::nullopt;
  const Vec3 first_point = ray.position_m + ray.direction * *hit;
  if (first_point.z <= center.z + kEpsilon) return hit;
  const Vec3 oc = ray.position_m - center;
  const double root = std::sqrt(std::max(0.0, dot(oc, ray.direction) * dot(oc, ray.direction) -
                                                   (dot(oc, oc) - radius * radius)));
  const double second = -dot(oc, ray.direction) + root;
  return second > kEpsilon && (ray.position_m + ray.direction * second).z <= center.z + kEpsilon
             ? std::optional<double>{second}
             : std::nullopt;
}

inline std::optional<double> intersect_plane_z(const Ray& ray, double z) {
  if (!std::isfinite(z) || !std::isfinite(ray.position_m.z) || !std::isfinite(ray.direction.z) ||
      std::abs(ray.direction.z) < kEpsilon) {
    return std::nullopt;
  }
  const double distance_m = (z - ray.position_m.z) / ray.direction.z;
  return std::isfinite(distance_m) && distance_m > kEpsilon ? std::optional<double>{distance_m}
                                                              : std::nullopt;
}

inline std::optional<double> intersect_disk_z(const Ray& ray, double z, double radius) {
  if (!std::isfinite(radius) || radius <= kEpsilon) return std::nullopt;
  const auto distance_m = intersect_plane_z(ray, z);
  if (!distance_m) return std::nullopt;
  const Vec3 point = ray.position_m + ray.direction * *distance_m;
  return point.x * point.x + point.y * point.y <= radius * radius ? distance_m : std::nullopt;
}

inline std::optional<double> intersect_finite_cylinder(const Ray& ray, const Vec3& a, const Vec3& b,
                                                        double radius) {
  if (!std::isfinite(radius) || radius <= kEpsilon) return std::nullopt;
  const Vec3 axis = b - a;
  const Vec3 offset = ray.position_m - a;
  const double axis_squared = dot(axis, axis);
  if (!std::isfinite(axis_squared) || axis_squared <= kEpsilon) return std::nullopt;
  const double direction_axis = dot(ray.direction, axis);
  const double offset_axis = dot(offset, axis);
  const double quadratic_a = dot(ray.direction, ray.direction) - direction_axis * direction_axis / axis_squared;
  const double quadratic_b = dot(ray.direction, offset) - direction_axis * offset_axis / axis_squared;
  const double quadratic_c = dot(offset, offset) - offset_axis * offset_axis / axis_squared - radius * radius;
  const double discriminant = quadratic_b * quadratic_b - quadratic_a * quadratic_c;
  if (quadratic_a <= kEpsilon || discriminant < 0.0) return std::nullopt;
  const double root = std::sqrt(std::max(0.0, discriminant));
  for (const double distance_m : {(-quadratic_b - root) / quadratic_a,
                                  (-quadratic_b + root) / quadratic_a}) {
    const double along_axis = offset_axis + distance_m * direction_axis;
    if (distance_m > kEpsilon && along_axis >= 0.0 && along_axis <= axis_squared) return distance_m;
  }
  return std::nullopt;
}

}  // namespace obdeect
