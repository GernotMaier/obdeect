#pragma once

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <optional>
#include <span>

namespace obdeect {

// Immutable monotonic tabulation for reflectivity, transmission or efficiency.
// It has no extrapolation: callers receive no value outside the physical table.
struct Table1DView {
  std::span<const double> axis;
  std::span<const double> value;

  [[nodiscard]] bool is_valid() const {
    if (axis.size() < 2 || axis.size() != value.size()) return false;
    for (std::size_t index = 0; index < axis.size(); ++index) {
      if (!std::isfinite(axis[index]) || !std::isfinite(value[index]) || value[index] < 0.0) return false;
      if (index > 0 && axis[index] <= axis[index - 1]) return false;
    }
    return true;
  }

  [[nodiscard]] std::optional<double> interpolate(double coordinate) const {
    if (!is_valid() || !std::isfinite(coordinate) || coordinate < axis.front() || coordinate > axis.back()) {
      return std::nullopt;
    }
    const auto upper = std::upper_bound(axis.begin(), axis.end(), coordinate);
    if (upper == axis.begin()) return value.front();
    if (upper == axis.end()) return value.back();
    const std::size_t high = static_cast<std::size_t>(upper - axis.begin());
    const std::size_t low = high - 1;
    const double fraction = (coordinate - axis[low]) / (axis[high] - axis[low]);
    return value[low] + fraction * (value[high] - value[low]);
  }
};

}  // namespace obdeect
