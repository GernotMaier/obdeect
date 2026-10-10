#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <span>
#include <vector>

namespace obdeect {
enum class TableBoundary { reject, clamp, zero };
struct TableInterpolation {
  TableBoundary boundary{TableBoundary::reject};
  unsigned scheme{1};
  bool x_log{}, y_log{}, value_log{};
  std::vector<std::array<double, 4>> coefficients{};
  [[nodiscard]] bool is_valid(std::size_t intervals, bool two_dimensional = false) const {
    if (scheme > 4 || (two_dimensional && (scheme != 1 || !coefficients.empty())))
      return false;
    if (scheme < 2 && !coefficients.empty())
      return false;
    if (scheme >= 2 && coefficients.size() != intervals)
      return false;
    for (const auto &row : coefficients)
      for (double v : row)
        if (!std::isfinite(v))
          return false;
    return true;
  }
};
[[nodiscard]] inline double table_fraction(double x, double lower, double upper, bool logarithmic) {
  return logarithmic ? std::log(x / lower) / std::log(upper / lower)
                     : (x - lower) / (upper - lower);
}
[[nodiscard]] inline std::optional<double>
interpolate_curve_unchecked(std::span<const double> axis, std::span<const double> values, double x,
                            const TableInterpolation &options) {
  if (!std::isfinite(x) || (options.x_log && x <= 0))
    return std::nullopt;
  if (x < axis.front() || x > axis.back()) {
    if (options.boundary == TableBoundary::reject)
      return std::nullopt;
    if (options.boundary == TableBoundary::zero)
      return 0;
    x = std::clamp(x, axis.front(), axis.back());
  }
  const auto upper = std::upper_bound(axis.begin(), axis.end(), x);
  const std::size_t high = std::clamp(std::size_t(upper - axis.begin()), std::size_t{1},
                                      axis.size() - 1),
                    low = high - 1;
  const double t = table_fraction(x, axis[low], axis[high], options.x_log);
  if (options.scheme == 0)
    return values[t < 0.5 ? low : high];
  if (options.scheme >= 2) {
    const auto &c = options.coefficients[low];
    const double value = ((c[3] * t + c[2]) * t + c[1]) * t + c[0];
    return options.value_log ? std::exp(value) : value;
  }
  if (options.value_log)
    return std::exp((1 - t) * std::log(values[low]) + t * std::log(values[high]));
  return (1 - t) * values[low] + t * values[high];
}
} // namespace obdeect
