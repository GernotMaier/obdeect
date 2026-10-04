#pragma once

#include "obdeect/interactions.hpp"
#include "obdeect/tables.hpp"

#include <cmath>
#include <optional>

namespace obdeect {

// Terminal outcome of the deterministic, transmitted branch of a window.
// Fresnel-reflected power is represented in ``transmission`` rather than
// spawning a second ray.  A caller that needs non-sequential reflected paths
// can use dielectric_interface directly and create that branch explicitly.
enum class SlabTransportStatus : unsigned char {
  transmitted,
  total_internal_reflection,
  absorbed,
  invalid_input,
  wavelength_out_of_range,
};

// Declare which losses are already included in a measured transmission table.
// Included losses replace their analytic counterpart; they never multiply it.
enum class SlabTransmissionSemantics : unsigned char {
  additional_coating,
  interfaces_only,
  bulk_only,
  complete_window,
};

struct SlabMaterial {
  // Phase refractive index and Beer--Lambert absorption coefficient, both as
  // a function of vacuum wavelength in nm.  The response is intentionally
  // data-only: an optical model importer may supply any observatory's material table.
  Table1DView refractive_index;
  Table1DView absorption_per_m;

  // Optional measured throughput, with a mandatory declaration of included
  // losses. ``interfaces_only`` covers both faces; ``bulk_only`` covers the
  // slab thickness for which that response was measured. These 1D responses
  // cannot supply unrecorded incidence or thickness dependence.
  std::optional<Table1DView> transmission{};
  std::optional<SlabTransmissionSemantics> transmission_semantics{};

  [[nodiscard]] bool is_valid() const {
    if (!refractive_index.is_valid() || !absorption_per_m.is_valid())
      return false;
    for (const double value : refractive_index.value)
      if (value <= 0.0)
        return false;
    if (transmission.has_value() != transmission_semantics.has_value())
      return false;
    if (transmission) {
      if (!transmission->is_valid())
        return false;
      for (const double value : transmission->value)
        if (value > 1.0)
          return false;
      switch (*transmission_semantics) {
      case SlabTransmissionSemantics::additional_coating:
      case SlabTransmissionSemantics::interfaces_only:
      case SlabTransmissionSemantics::bulk_only:
      case SlabTransmissionSemantics::complete_window:
        break;
      default:
        return false;
      }
    }
    return true;
  }
};

struct SlabTransportResult {
  Vec3 direction{};
  Vec3 position_m{};
  double geometric_path_m{};
  double optical_path_m{};
  double transmission{};
  // Bare dielectric reflectances are diagnostics. They are not applied loss
  // fractions when a measured table replaces the interface response.
  double entry_reflectance{};
  double exit_reflectance{};
  SlabTransportStatus status{SlabTransportStatus::invalid_input};
};

// Trace the transmitted branch through a parallel-sided slab.  ``normal``
// points into the incident medium and ``entry_position_m`` is on the first
// face.  The same incident medium is assumed beyond the second face; use two
// dielectric_interface calls directly for an asymmetric medium stack.
[[nodiscard]] inline SlabTransportResult
transmit_through_slab(Vec3 incident, Vec3 normal, Vec3 entry_position_m, double n_incident,
                      double thickness_m, double wavelength_nm, const SlabMaterial &material) {
  SlabTransportResult result{};
  const auto d = normalised_checked(incident);
  const auto n = normalised_checked(normal);
  if (!d || !n || !material.is_valid() || !std::isfinite(n_incident) || n_incident <= 0.0 ||
      !std::isfinite(thickness_m) || thickness_m < 0.0 || !std::isfinite(wavelength_nm) ||
      wavelength_nm <= 0.0 || !std::isfinite(entry_position_m.x) ||
      !std::isfinite(entry_position_m.y) || !std::isfinite(entry_position_m.z)) {
    return result;
  }

  const auto n_inside = material.refractive_index.interpolate(wavelength_nm);
  const auto absorption = material.absorption_per_m.interpolate(wavelength_nm);
  const auto response = material.transmission ? material.transmission->interpolate(wavelength_nm)
                                              : std::optional<double>{1.0};
  if (!n_inside || !absorption || !response) {
    result.status = SlabTransportStatus::wavelength_out_of_range;
    return result;
  }

  const auto entry = dielectric_interface(*d, *n, n_incident, *n_inside);
  if (!entry)
    return result;
  result.direction = entry->reflected;
  result.position_m = entry_position_m;
  result.entry_reflectance = 0.5 * (entry->reflectance_s + entry->reflectance_p);
  if (!entry->transmitted) {
    result.status = SlabTransportStatus::total_internal_reflection;
    return result;
  }

  const double cos_inside = -dot(*entry->transmitted, *n);
  if (!std::isfinite(cos_inside) || cos_inside <= 0.0)
    return result;
  result.geometric_path_m = thickness_m / cos_inside;
  result.position_m = entry_position_m + *entry->transmitted * result.geometric_path_m;

  // At the far face the incident medium is still the slab, so its outward
  // normal points back toward the entry face, i.e. the same normal vector.
  const auto exit = dielectric_interface(*entry->transmitted, *n, *n_inside, n_incident);
  if (!exit)
    return SlabTransportResult{};
  result.exit_reflectance = 0.5 * (exit->reflectance_s + exit->reflectance_p);
  result.direction = exit->reflected;
  if (!exit->transmitted) {
    result.status = SlabTransportStatus::total_internal_reflection;
    return result;
  }

  result.direction = *exit->transmitted;
  result.optical_path_m = *n_inside * result.geometric_path_m;
  const double entry_t = 1.0 - result.entry_reflectance;
  const double exit_t = 1.0 - result.exit_reflectance;
  const bool includes_interfaces =
      material.transmission_semantics == SlabTransmissionSemantics::interfaces_only ||
      material.transmission_semantics == SlabTransmissionSemantics::complete_window;
  const bool includes_bulk =
      material.transmission_semantics == SlabTransmissionSemantics::bulk_only ||
      material.transmission_semantics == SlabTransmissionSemantics::complete_window;
  const double interfaces = includes_interfaces ? 1.0 : entry_t * exit_t;
  const double bulk = includes_bulk ? 1.0 : std::exp(-*absorption * result.geometric_path_m);
  result.transmission = interfaces * bulk * *response;
  if (!std::isfinite(result.transmission) || result.transmission < 0.0 ||
      result.transmission > 1.0) {
    return SlabTransportResult{};
  }
  result.status =
      result.transmission == 0.0 ? SlabTransportStatus::absorbed : SlabTransportStatus::transmitted;
  return result;
}

} // namespace obdeect
