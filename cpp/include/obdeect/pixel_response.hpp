#pragma once

#include "obdeect/segmented_optical_model.hpp"

#include <algorithm>
#include <numbers>
#include <optional>
#include <vector>

namespace obdeect {
struct PixelResponseTable {
  enum class Method { measured, single_reflection };
  std::uint32_t id{};
  Method method{Method::measured};
  double tangent_bin_width{}, wavelength_bin_width_nm{}, wavelength_bin_origin_nm{};
  std::vector<double> angular_efficiency{}, spectral_correction{};
  double transparency{1}, wall_reflectivity{1};
};
struct PixelResponseBinding {
  std::uint32_t detector_id{}, table_id{};
  std::optional<ImportedDetectorSurface> cathode{};
};
struct PixelResponses {
  std::vector<PixelResponseTable> tables;
  std::vector<PixelResponseBinding> bindings;

  [[nodiscard]] bool is_valid() const {
    for (std::size_t i = 0; i < tables.size(); ++i) {
      const auto &table = tables[i];
      if (table.method != PixelResponseTable::Method::measured &&
          table.method != PixelResponseTable::Method::single_reflection)
        return false;
      if (i && table.id <= tables[i - 1].id)
        return false;
      const auto fractional = [](double v) { return std::isfinite(v) && v >= 0 && v <= 1; };
      if (!fractional(table.transparency) || !fractional(table.wall_reflectivity))
        return false;
      if (table.method == PixelResponseTable::Method::measured) {
        if (table.angular_efficiency.empty() || !std::isfinite(table.tangent_bin_width) ||
            table.tangent_bin_width <= 0 ||
            !std::all_of(table.angular_efficiency.begin(), table.angular_efficiency.end(),
                         fractional))
          return false;
        if (!table.spectral_correction.empty()) {
          if (!std::isfinite(table.wavelength_bin_width_nm) || table.wavelength_bin_width_nm <= 0 ||
              !std::isfinite(table.wavelength_bin_origin_nm) ||
              !std::all_of(table.spectral_correction.begin(), table.spectral_correction.end(),
                           [](double v) { return std::isfinite(v) && v >= 0; }))
            return false;
        }
      } else if (!table.angular_efficiency.empty() || !table.spectral_correction.empty())
        return false;
    }
    for (std::size_t i = 0; i < bindings.size(); ++i) {
      const auto &binding = bindings[i];
      if (i && binding.detector_id <= bindings[i - 1].detector_id)
        return false;
      const auto table = std::lower_bound(tables.begin(), tables.end(), binding.table_id,
                                          [](const auto &t, auto id) { return t.id < id; });
      if (table == tables.end() || table->id != binding.table_id ||
          (table->method == PixelResponseTable::Method::single_reflection && !binding.cathode) ||
          (binding.cathode && !obdeect::is_valid(*binding.cathode)))
        return false;
    }
    return !tables.empty() && !bindings.empty();
  }

  // Lookup is deterministic, allocation-free, and tied to the actual detector.
  [[nodiscard]] std::optional<double> at_unchecked(std::uint32_t detector_id, double wavelength_nm,
                                                   double incidence_deg,
                                                   const Ray &at_entrance) const {
    const auto binding =
        std::lower_bound(bindings.begin(), bindings.end(), detector_id,
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
    const double index =
        std::tan(incidence_deg * std::numbers::pi / 180) / table->tangent_bin_width;
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
};
} // namespace obdeect
