#pragma once

#include "obdeect/math.hpp"

#include <optional>

namespace obdeect {

[[nodiscard]] inline std::optional<Vec3> reflect_specular(const Vec3& incident, const Vec3& normal) {
  const auto unit_incident = normalised_checked(incident);
  const auto unit_normal = normalised_checked(normal);
  if (!unit_incident || !unit_normal) return std::nullopt;
  return normalised_checked(*unit_incident - *unit_normal * (2.0 * dot(*unit_incident, *unit_normal)));
}

}  // namespace obdeect
