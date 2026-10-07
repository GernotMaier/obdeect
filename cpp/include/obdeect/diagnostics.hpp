#pragma once

#include "obdeect/photon_buffer.hpp"

#include <array>

namespace obdeect {

struct TraceSummary {
  std::array<std::size_t, kPhotonStatusCount> status_count{};
  double input_weight{};
  double detected_weight{};
  double lost_weight{};
  double response_loss_weight{};
  std::array<double, kPhotonStatusCount> terminal_loss_weight{};

  // response_loss is optical absorption before the terminal interaction.
  // remaining_weight is the weight reaching that terminal, whether accepted
  // or lost there. Counts and incident weight remain independent ledgers.
  void add(PhotonStatus status, double weight, double remaining_weight,
           double response_loss = 0.0) {
    ++status_count[static_cast<std::size_t>(status)];
    input_weight += weight;
    response_loss_weight += response_loss;
    const double terminal_loss = status == PhotonStatus::detected ? 0.0 : remaining_weight;
    terminal_loss_weight[static_cast<std::size_t>(status)] += terminal_loss;
    lost_weight += response_loss + terminal_loss;
    if (status == PhotonStatus::detected)
      detected_weight += remaining_weight;
  }

  void add(PhotonStatus status, double weight) { add(status, weight, weight); }
};

} // namespace obdeect
