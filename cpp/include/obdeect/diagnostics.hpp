#pragma once

#include "obdeect/photon_buffer.hpp"

#include <array>

namespace obdeect {

struct TraceSummary {
  std::array<std::size_t, kPhotonStatusCount> status_count{};
  double input_weight{};
  double detected_weight{};

  void add(PhotonStatus status, double weight) {
    ++status_count[static_cast<std::size_t>(status)];
    input_weight += weight;
    if (status == PhotonStatus::detected) detected_weight += weight;
  }
};

}  // namespace obdeect
