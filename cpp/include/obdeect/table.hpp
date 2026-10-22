#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

namespace obdeect {

namespace detail {
// Index of the lower node and the fraction towards the next one, clamped to
// the axis. A single-node axis is constant.
struct Bracket {
  std::size_t index;
  double fraction;
};

inline Bracket bracket(const std::vector<double> &axis, double value) {
  if (axis.size() < 2 || value <= axis.front())
    return {0, 0.0};
  if (value >= axis.back())
    return {axis.size() - 2, 1.0};
  const auto upper = std::upper_bound(axis.begin(), axis.end(), value);
  const auto index = static_cast<std::size_t>(upper - axis.begin()) - 1;
  return {index, (value - axis[index]) / (axis[index + 1] - axis[index])};
}
} // namespace detail

// Optical response versus wavelength (uniform grid) and incidence angle.
// The compiler resamples the source table with its declared interpolation
// onto the wavelength grid; outside the grid the edge value applies.
struct Response {
  double wavelength_min_nm{};
  double wavelength_step_nm{1.0};
  std::size_t wavelength_count{};
  std::vector<double> incidence_deg; // One entry for a wavelength-only table.
  std::vector<double> values;        // [wavelength][incidence], row-major.

  [[nodiscard]] double at(double wavelength_nm, double incidence_deg_value) const {
    const double position = std::clamp((wavelength_nm - wavelength_min_nm) / wavelength_step_nm,
                                       0.0, static_cast<double>(wavelength_count - 1));
    const auto w = std::min(static_cast<std::size_t>(position), wavelength_count - 1);
    const double fw = position - static_cast<double>(w);
    const std::size_t w1 = std::min(w + 1, wavelength_count - 1);
    const std::size_t na = incidence_deg.size();
    const auto a = detail::bracket(incidence_deg, incidence_deg_value);
    const std::size_t a1 = std::min(a.index + 1, na - 1);
    const auto row = [&](std::size_t i) {
      return values[i * na + a.index] * (1 - a.fraction) + values[i * na + a1] * a.fraction;
    };
    return row(w) * (1 - fw) + row(w1) * fw;
  }
};

// Relative response map on a rectangular grid of surface coordinates.
struct Map2D {
  std::vector<double> x_m, y_m;
  std::vector<double> values; // [x][y], row-major.
  bool zero_outside{};

  [[nodiscard]] double at(double x, double y) const {
    if (zero_outside && (x < x_m.front() || x > x_m.back() || y < y_m.front() || y > y_m.back()))
      return 0.0;
    const auto bx = detail::bracket(x_m, x), by = detail::bracket(y_m, y);
    const std::size_t ny = y_m.size();
    const std::size_t x1 = std::min(bx.index + 1, x_m.size() - 1);
    const std::size_t y1 = std::min(by.index + 1, ny - 1);
    const auto v = [&](std::size_t i, std::size_t j) { return values[i * ny + j]; };
    return (1 - bx.fraction) *
               ((1 - by.fraction) * v(bx.index, by.index) + by.fraction * v(bx.index, y1)) +
           bx.fraction * ((1 - by.fraction) * v(x1, by.index) + by.fraction * v(x1, y1));
  }
};

} // namespace obdeect
