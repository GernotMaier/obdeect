#pragma once

#include "obdeect/abi.hpp"
#include "obdeect/diagnostics.hpp"
#include "obdeect/material_path.hpp"
#include "obdeect/optical_model.hpp"
#include "obdeect/reference_optical_model.hpp"
#include "obdeect/segmented_optical_model.hpp"
#include "obdeect/segmented_path.hpp"

namespace obdeect {

struct TraceResult {
  PhotonResultBlock photons;
  TraceSummary summary;
};

// Scalar reference block runner. It is intentionally allocation-free inside
// the photon loop: output is sized before tracing and each record is written by
// index. SIMD/threaded kernels must reproduce this contract exactly.
[[nodiscard]] inline TraceResult trace(const CompiledReferenceOpticalModel &optical_model,
                                       const PhotonBlockView &input) {
  TraceResult result{PhotonResultBlock{input.position_m.size()}, {}};
  if (!optical_model.is_valid() || !validate_photon_block(input)) {
    for (std::size_t index = 0; index < input.position_m.size(); ++index) {
      result.photons.status[index] = PhotonStatus::invalid_input;
      ++result.summary.status_count[static_cast<std::size_t>(PhotonStatus::invalid_input)];
    }
    return result;
  }
  for (std::size_t index = 0; index < input.position_m.size(); ++index) {
    const PathRecord record =
        trace_artificial_mst({input.position_m[index], input.direction[index]},
                             input.photon_id[index], optical_model.configuration);
    const std::size_t final_index = record.point_count == 0 ? 0 : record.point_count - 1;
    result.photons.position_m[index] = record.points_m[final_index];
    result.photons.direction[index] = record.final_direction;
    result.photons.optical_path_m[index] = record.path_length_m;
    result.photons.time_ns[index] =
        input.time_ns[index] + record.path_length_m / kSpeedOfLightMPerNs;
    result.photons.weight[index] =
        record.status == PhotonStatus::detected ? input.weight[index] : 0.0;
    result.photons.status[index] = record.status;
    result.summary.add(record.status, input.weight[index]);
  }
  return result;
}

// A segmented optical_model may supply finite detector surfaces. A detected photon
// records the nearest physical post-reflection intersection; otherwise it
// retains the explicit no-detector terminal status.
[[nodiscard]] inline TraceResult trace(const CompiledSegmentedOpticalModel &optical_model,
                                       const PhotonBlockView &input) {
  TraceResult result{PhotonResultBlock{input.position_m.size()}, {}};
  if (!is_valid(optical_model) || !validate_photon_block(input)) {
    for (std::size_t index = 0; index < input.position_m.size(); ++index) {
      result.photons.status[index] = PhotonStatus::invalid_input;
      ++result.summary.status_count[static_cast<std::size_t>(PhotonStatus::invalid_input)];
    }
    return result;
  }
  for (std::size_t index = 0; index < input.position_m.size(); ++index) {
    const auto path =
        trace_segmented_path({input.position_m[index], input.direction[index]},
                             input.photon_id[index], input.wavelength_nm[index], optical_model);
    result.photons.position_m[index] = path.points_m[path.point_count - 1];
    result.photons.direction[index] = path.final_direction;
    result.photons.optical_path_m[index] = path.path_length_m;
    result.photons.time_ns[index] = input.time_ns[index] + path.path_length_m / kSpeedOfLightMPerNs;
    const double remaining_weight = input.weight[index] * path.surviving_throughput;
    result.photons.weight[index] = path.status == PhotonStatus::detected ? remaining_weight : 0;
    result.photons.surface_id[index] = path.terminal_surface_id;
    result.photons.status[index] = path.status;
    result.summary.add(path.status, input.weight[index], remaining_weight,
                       input.weight[index] - remaining_weight);
  }
  return result;
}

// Nonsequential finite-surface transport with explicit dielectric media.
[[nodiscard]] inline TraceResult trace(const CompiledOpticalModel &optical_model,
                                       const PhotonBlockView &input) {
  TraceResult result{PhotonResultBlock{input.position_m.size()}, {}};
  if (!validate_photon_block(input)) {
    for (std::size_t index = 0; index < input.position_m.size(); ++index) {
      result.photons.status[index] = PhotonStatus::invalid_input;
      ++result.summary.status_count[static_cast<std::size_t>(PhotonStatus::invalid_input)];
    }
    return result;
  }
  for (std::size_t index = 0; index < input.position_m.size(); ++index) {
    const auto path = trace_material_path<false>({input.position_m[index], input.direction[index]},
                                                 input.photon_id[index], input.wavelength_nm[index],
                                                 optical_model);
    result.photons.position_m[index] = path.points_m[path.point_count - 1];
    result.photons.direction[index] = path.final_direction;
    result.photons.optical_path_m[index] = path.optical_path_m;
    result.photons.geometric_path_m[index] = path.path_length_m;
    result.photons.time_ns[index] = input.time_ns[index] + path.group_delay_ns;
    result.photons.surface_id[index] = path.terminal_surface_id;
    result.photons.status[index] = path.status;
    const double remaining = input.weight[index] * path.surviving_throughput;
    result.photons.weight[index] = path.status == PhotonStatus::detected ? remaining : 0;
    result.summary.add(path.status, input.weight[index], remaining,
                       input.weight[index] - remaining);
  }
  return result;
}

} // namespace obdeect
