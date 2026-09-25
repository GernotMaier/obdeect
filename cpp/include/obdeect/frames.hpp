#pragma once

#include "obdeect/math.hpp"

#include <cmath>

namespace obdeect {

// A right-handed rigid frame. Axes and origin are expressed in the parent
// frame, with metres for positions. No scale or implicit unit conversion.
struct RigidFrame {
  Vec3 origin_m{};
  Vec3 x_axis{1.0, 0.0, 0.0};
  Vec3 y_axis{0.0, 1.0, 0.0};
  Vec3 z_axis{0.0, 0.0, 1.0};

  [[nodiscard]] bool is_valid() const {
    constexpr double tolerance = 1e-12;
    const auto finite = [](Vec3 vector) {
      return std::isfinite(vector.x) && std::isfinite(vector.y) && std::isfinite(vector.z);
    };
    return finite(origin_m) && finite(x_axis) && finite(y_axis) && finite(z_axis) &&
           std::abs(dot(x_axis, x_axis) - 1.0) <= tolerance &&
           std::abs(dot(y_axis, y_axis) - 1.0) <= tolerance &&
           std::abs(dot(z_axis, z_axis) - 1.0) <= tolerance &&
           std::abs(dot(x_axis, y_axis)) <= tolerance &&
           std::abs(dot(x_axis, z_axis)) <= tolerance &&
           std::abs(dot(y_axis, z_axis)) <= tolerance &&
           norm(cross(x_axis, y_axis) - z_axis) <= tolerance;
  }

  [[nodiscard]] Vec3 point_to_parent(Vec3 point_m) const {
    return origin_m + x_axis * point_m.x + y_axis * point_m.y + z_axis * point_m.z;
  }

  [[nodiscard]] Vec3 point_from_parent(Vec3 point_m) const {
    const Vec3 delta = point_m - origin_m;
    return {dot(delta, x_axis), dot(delta, y_axis), dot(delta, z_axis)};
  }

  [[nodiscard]] Vec3 direction_to_parent(Vec3 direction) const {
    return x_axis * direction.x + y_axis * direction.y + z_axis * direction.z;
  }

  [[nodiscard]] Vec3 direction_from_parent(Vec3 direction) const {
    return {dot(direction, x_axis), dot(direction, y_axis), dot(direction, z_axis)};
  }

  [[nodiscard]] Ray ray_to_parent(Ray ray) const {
    return {point_to_parent(ray.position_m), direction_to_parent(ray.direction)};
  }

  [[nodiscard]] Ray ray_from_parent(Ray ray) const {
    return {point_from_parent(ray.position_m), direction_from_parent(ray.direction)};
  }
};

}  // namespace obdeect
