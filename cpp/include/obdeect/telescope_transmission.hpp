#pragma once

#include "obdeect/math.hpp"
#include <cmath>
#include <optional>

namespace obdeect {
struct TelescopeTransmission {
  double on_axis{1}, amplitude{}, angular_scale_rad{1}, power{2}, outer_power{1};
  [[nodiscard]] bool is_valid() const {
    return std::isfinite(on_axis) && on_axis >= 0 && on_axis <= 1 && std::isfinite(amplitude) &&
           amplitude >= 0 && std::isfinite(angular_scale_rad) && angular_scale_rad > 0 &&
           std::isfinite(power) && power > 0 && std::isfinite(outer_power) && outer_power > 0;
  }
  [[nodiscard]] std::optional<double> at_unchecked(const Vec3 &direction) const {
    const auto unit = normalised_checked(direction);
    if (!unit)
      return {};
    const double sine = std::hypot(unit->x, unit->y);
    return on_axis /
           std::pow(1 + amplitude * std::pow(sine / angular_scale_rad, power), outer_power);
  }
};
} // namespace obdeect
