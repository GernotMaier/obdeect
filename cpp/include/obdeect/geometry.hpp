#pragma once

#include "obdeect/axisymmetric_optics.hpp"
#include "obdeect/intersections.hpp"

#include <variant>

namespace obdeect {

struct DiskSurface {
  double z_m{};
  double radius_m{};
};

struct CylinderObstruction {
  Vec3 start_m{};
  Vec3 end_m{};
  double radius_m{};
};

using OpticalSurface = std::variant<AxisymmetricMirror, DiskSurface>;

[[nodiscard]] inline std::optional<double> intersect(const Ray& ray, const DiskSurface& surface) {
  return intersect_disk_z(ray, surface.z_m, surface.radius_m);
}

[[nodiscard]] inline std::optional<double> intersect(const Ray& ray, const CylinderObstruction& obstruction) {
  return intersect_finite_cylinder(ray, obstruction.start_m, obstruction.end_m, obstruction.radius_m);
}

}  // namespace obdeect
