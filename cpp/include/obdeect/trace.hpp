#pragma once

#include "obdeect/abi.hpp"
#include "obdeect/diagnostics.hpp"
#include "obdeect/scene.hpp"

namespace obdeect {

struct TraceResult {
  PhotonResultBlock photons;
  TraceSummary summary;
};

// Scalar reference block runner. It is intentionally allocation-free inside
// the photon loop: output is sized before tracing and each record is written by
// index. SIMD/threaded kernels must reproduce this contract exactly.
[[nodiscard]] inline TraceResult trace(const CompiledToyScene& scene, const PhotonBlockView& input) {
  TraceResult result{PhotonResultBlock{input.position_m.size()}, {}};
  if (!scene.is_valid() || !validate_photon_block(input)) {
    for (std::size_t index = 0; index < input.position_m.size(); ++index) {
      result.photons.status[index] = PhotonStatus::invalid_input;
    }
    return result;
  }
  constexpr double speed_of_light_m_per_ns = 0.299792458;
  for (std::size_t index = 0; index < input.position_m.size(); ++index) {
    const PathRecord record = trace_toy_mst({input.position_m[index], input.direction[index]}, input.photon_id[index],
                                            scene.configuration);
    const std::size_t final_index = record.point_count == 0 ? 0 : record.point_count - 1;
    result.photons.position_m[index] = record.points_m[final_index];
    result.photons.direction[index] = input.direction[index];
    result.photons.optical_path_m[index] = record.path_length_m;
    result.photons.time_ns[index] = input.time_ns[index] + record.path_length_m / speed_of_light_m_per_ns;
    result.photons.weight[index] = record.status == PhotonStatus::detected ? input.weight[index] : 0.0;
    result.photons.status[index] = record.status;
    result.summary.add(record.status, input.weight[index]);
  }
  return result;
}

}  // namespace obdeect
