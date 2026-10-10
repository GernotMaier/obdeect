#pragma once

#include "obdeect/segmented_optical_model.hpp"

#include <algorithm>
#include <numbers>
#include <optional>
#include <vector>

namespace obdeect {
// Lookup is deterministic, allocation-free, and tied to the actual detector.
inline std::optional<double> PixelResponses::at_unchecked(std::uint32_t detector_id,
                                                          double wavelength_nm,
                                                          double incidence_deg,
                                                          const Ray &at_entrance) const {
  const auto binding = std::lower_bound(bindings.begin(), bindings.end(), detector_id,
                                        [](const auto &b, auto id) { return b.detector_id < id; });
  if (binding == bindings.end() || binding->detector_id != detector_id)
    return std::nullopt;
  const auto table = std::lower_bound(tables.begin(), tables.end(), binding->table_id,
                                      [](const auto &t, auto id) { return t.id < id; });
  if (!std::isfinite(wavelength_nm) || wavelength_nm <= 0 || !std::isfinite(incidence_deg) ||
      incidence_deg < 0 || incidence_deg >= 90)
    return 0;
  if (table->method == PixelResponseTable::Method::single_reflection) {
    const bool on_cathode = std::abs(dot(at_entrance.position_m - binding->cathode->centre_m,
                                         binding->cathode->unit_normal)) <= kEpsilon &&
                            contains_detector_point(*binding->cathode, at_entrance.position_m);
    return table->transparency *
           (on_cathode || intersect_detector_surface_unchecked(at_entrance, *binding->cathode)
                ? 1
                : table->wall_reflectivity);
  }
  const double index = std::tan(incidence_deg * std::numbers::pi / 180) / table->tangent_bin_width;
  if (index < 0 || index >= double(table->angular_efficiency.size()))
    return 0;
  double value = table->angular_efficiency[std::size_t(index)];
  if (!table->spectral_correction.empty()) {
    const double index = std::floor(
        (wavelength_nm - table->wavelength_bin_origin_nm) / table->wavelength_bin_width_nm + 0.5);
    if (index < 0 || index >= double(table->spectral_correction.size()))
      return 0;
    value *= table->spectral_correction[std::size_t(index)];
  }
  return value;
}
} // namespace obdeect
