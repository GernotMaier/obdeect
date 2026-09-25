#pragma once

#include "obdeect/photon_buffer.hpp"

namespace obdeect {

[[nodiscard]] inline bool validate_photon_block(const PhotonBlockView& photons) {
  if (!photons.is_consistent()) return false;
  for (std::size_t index = 0; index < photons.position_m.size(); ++index) {
    const double direction_length = norm(photons.direction[index]);
    if (!std::isfinite(direction_length) || std::abs(direction_length - 1.0) > 1e-12 ||
        !std::isfinite(photons.position_m[index].x) ||
        !std::isfinite(photons.position_m[index].y) || !std::isfinite(photons.position_m[index].z) ||
        !std::isfinite(photons.wavelength_nm[index]) || photons.wavelength_nm[index] <= 0.0 ||
        !std::isfinite(photons.time_ns[index]) || !std::isfinite(photons.weight[index]) ||
        photons.weight[index] < 0.0) {
      return false;
    }
  }
  return true;
}

}  // namespace obdeect
