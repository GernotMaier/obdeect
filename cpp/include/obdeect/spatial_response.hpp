#pragma once

#include "obdeect/interpolation.hpp"
#include "obdeect/math.hpp"
#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

namespace obdeect {
// Rectangular physical-coordinate degradation map; values are fractions.
struct SpatialResponse {
  std::vector<double> x_m, y_m, response;
  bool clip{};
  TableInterpolation interpolation{TableBoundary::clamp};
  Vec3 x_basis{1, 0, 0}, y_basis{0, 1, 0};
  [[nodiscard]] bool is_valid() const {
    const auto valid_axis = [](const auto &axis) {
      if (axis.empty())
        return false;
      for (std::size_t i = 0; i < axis.size(); ++i)
        if (!std::isfinite(axis[i]) || (i && axis[i] <= axis[i - 1]))
          return false;
      return true;
    };
    if (!interpolation.is_valid(0, true))
      return false;
    if ((interpolation.x_log && (x_m.empty() || x_m.front() <= 0)) ||
        (interpolation.y_log && (y_m.empty() || y_m.front() <= 0)) ||
        (interpolation.value_log &&
         std::any_of(response.begin(), response.end(), [](double v) { return v <= 0; })))
      return false;
    const auto u = normalised_checked(x_basis), v = normalised_checked(y_basis);
    if (!u || !v || std::abs(dot(*u, *v)) > kEpsilon || std::abs(norm(x_basis) - 1) > kEpsilon ||
        std::abs(norm(y_basis) - 1) > kEpsilon)
      return false;
    return valid_axis(x_m) && valid_axis(y_m) && response.size() == x_m.size() * y_m.size() &&
           std::all_of(response.begin(), response.end(),
                       [](double v) { return std::isfinite(v) && v >= 0 && v <= 1; });
  }
  [[nodiscard]] std::optional<double> at_point_unchecked(const Vec3 &point) const {
    return at_unchecked(dot(point, x_basis), dot(point, y_basis));
  }
  [[nodiscard]] std::optional<double> at_unchecked(double x, double y) const {
    if (!std::isfinite(x) || !std::isfinite(y))
      return {};
    if (x < x_m.front() || x > x_m.back() || y < y_m.front() || y > y_m.back()) {
      if (clip || interpolation.boundary == TableBoundary::zero)
        return 0;
      if (interpolation.boundary == TableBoundary::reject)
        return {};
    }
    x = std::clamp(x, x_m.front(), x_m.back());
    y = std::clamp(y, y_m.front(), y_m.back());
    const auto bracket = [](const auto &axis, double value) {
      return axis.size() == 1
                 ? std::size_t{0}
                 : std::min(std::size_t(std::upper_bound(axis.begin(), axis.end(), value) -
                                        axis.begin() - 1),
                            axis.size() - 2);
    };
    const auto i = bracket(x_m, x), j = bracket(y_m, y), ip = std::min(i + 1, x_m.size() - 1),
               jp = std::min(j + 1, y_m.size() - 1);
    const double a = ip == i ? 0 : table_fraction(x, x_m[i], x_m[ip], interpolation.x_log);
    const double b = jp == j ? 0 : table_fraction(y, y_m[j], y_m[jp], interpolation.y_log);
    const auto value = [&](std::size_t xx, std::size_t yy) {
      return interpolation.value_log ? std::log(response[xx * y_m.size() + yy])
                                     : response[xx * y_m.size() + yy];
    };
    const double result = (1 - a) * ((1 - b) * value(i, j) + b * value(i, jp)) +
                          a * ((1 - b) * value(ip, j) + b * value(ip, jp));
    return interpolation.value_log ? std::exp(result) : result;
  }
};
} // namespace obdeect
