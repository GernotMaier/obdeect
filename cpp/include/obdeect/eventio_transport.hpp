#pragma once

#include "obdeect/atmosphere.hpp"
#include "obdeect/eventio_photon_input.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace obdeect {
namespace detail {
[[nodiscard]] inline std::uint64_t mix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}
}  // namespace detail

// Spectral children retain the bunch's arrival ray and time and divide its
// photon count exactly. A deterministic stratified draw approximates the
// CORSIKA 1/lambda^2 spectrum without assigning a fake monochromatic colour.
// The caller selects children per bunch to control spectral quadrature error.
[[nodiscard]] inline std::vector<OpticalPhoton> resolve_eventio_spectrum(
    const OpticalPhoton& bunch, const EventioRunInfo& run, std::size_t children,
    std::uint64_t seed = 0) {
  if (!run.has_event_header || run.ceffic || !valid_photon(bunch) || children == 0 ||
      children > 65536 || run.wavelength_lower_nm <= 0 ||
      run.wavelength_upper_nm <= run.wavelength_lower_nm)
    throw std::invalid_argument("invalid EventIO spectral resolution request");
  if (bunch.wavelength_nm > 0) return {bunch};
  std::vector<OpticalPhoton> result;
  result.reserve(children);
  const double inverse_lower = 1.0 / run.wavelength_lower_nm;
  const double inverse_span = inverse_lower - 1.0 / run.wavelength_upper_nm;
  for (std::size_t index = 0; index < children; ++index) {
    const auto bits = detail::mix64(seed ^ bunch.photon_id ^
                                    (static_cast<std::uint64_t>(index) * 0x9e3779b97f4a7c15ULL));
    const double jitter = static_cast<double>(bits >> 11) * 0x1.0p-53;
    const double u = (static_cast<double>(index) + jitter) / static_cast<double>(children);
    auto child = bunch;
    child.wavelength_nm = 1.0 / (inverse_lower - u * inverse_span);
    child.weight = bunch.weight / static_cast<double>(children);
    child.photon_id = detail::mix64(bunch.photon_id ^
                                    (static_cast<std::uint64_t>(index + 1) * 0x9e3779b97f4a7c15ULL));
    result.push_back(child);
  }
  return result;
}

// Apply direct-beam extinction once between source and telescope entrance.
// The return value is weighted light removed from the direct beam. The caller
// records this separately from losses inside the optical scene.
[[nodiscard]] inline double attenuate_eventio_direct_beam(
    std::span<OpticalPhoton> photons, const PhotonBatchContext& context,
    const EventioRunInfo& run, const AtmosphereTransmissionTable& table) {
  if (!run.has_event_header || run.ceffic)
    throw std::invalid_argument("EventIO atmosphere stage cannot process CEFFIC input");
  if (run.curved || run.refraction)
    throw std::invalid_argument("curved/refracted EventIO rays need a path-integral atmosphere model");
  if (!std::isfinite(run.observation_altitude_m) ||
      std::abs(run.observation_altitude_m - table.observation_altitude_m()) > 25.0)
    throw std::invalid_argument("EventIO observation altitude differs from atmosphere table");
  std::vector<double> survival;
  survival.reserve(photons.size());
  for (const auto& photon : photons) survival.push_back(table.direct_survival(photon, context));
  double removed = 0;
  for (std::size_t index = 0; index < photons.size(); ++index) {
    auto& photon = photons[index];
    const double old_weight = photon.weight;
    photon.weight *= survival[index];
    removed += old_weight - photon.weight;
  }
  return removed;
}

}  // namespace obdeect
