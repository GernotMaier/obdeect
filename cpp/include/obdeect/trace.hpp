#pragma once

#include "obdeect/abi.hpp"
#include "obdeect/diagnostics.hpp"
#include "obdeect/scene.hpp"
#include "obdeect/segmented_scene.hpp"

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
      ++result.summary.status_count[static_cast<std::size_t>(PhotonStatus::invalid_input)];
    }
    return result;
  }
  for (std::size_t index = 0; index < input.position_m.size(); ++index) {
    const PathRecord record = trace_toy_mst({input.position_m[index], input.direction[index]}, input.photon_id[index],
                                            scene.configuration);
    const std::size_t final_index = record.point_count == 0 ? 0 : record.point_count - 1;
    result.photons.position_m[index] = record.points_m[final_index];
    result.photons.direction[index] = record.final_direction;
    result.photons.optical_path_m[index] = record.path_length_m;
    result.photons.time_ns[index] = input.time_ns[index] + record.path_length_m / kSpeedOfLightMPerNs;
    result.photons.weight[index] = record.status == PhotonStatus::detected ? input.weight[index] : 0.0;
    result.photons.status[index] = record.status;
    result.summary.add(record.status, input.weight[index]);
  }
  return result;
}

// A segmented scene may supply finite detector surfaces. A detected photon
// records the nearest physical post-reflection intersection; otherwise it
// retains the explicit no-detector terminal status.
[[nodiscard]] inline TraceResult trace(const CompiledSegmentedScene& scene, const PhotonBlockView& input) {
  TraceResult result{PhotonResultBlock{input.position_m.size()}, {}};
  if (!is_valid(scene) || !validate_photon_block(input)) {
    for (std::size_t index = 0; index < input.position_m.size(); ++index) {
      result.photons.status[index] = PhotonStatus::invalid_input;
      ++result.summary.status_count[static_cast<std::size_t>(PhotonStatus::invalid_input)];
    }
    return result;
  }
  for (std::size_t index = 0; index < input.position_m.size(); ++index) {
    const auto direction = normalised_checked(input.direction[index]);
    result.photons.position_m[index] = input.position_m[index];
    result.photons.direction[index] = direction.value_or(Vec3{});
    result.photons.time_ns[index] = input.time_ns[index];
    if (!direction) {
      result.photons.status[index] = PhotonStatus::invalid_input;
      result.summary.add(PhotonStatus::invalid_input, input.weight[index]);
      continue;
    }
    const Ray ray{input.position_m[index], *direction};
    const auto hit = intersect_segmented_primary_unchecked(ray, scene);
    if (!hit) {
      result.photons.status[index] = PhotonStatus::missed_primary;
      result.summary.add(PhotonStatus::missed_primary, input.weight[index]);
      continue;
    }
    const auto reflected = reflect_specular(ray.direction, hit->unit_normal);
    if (!reflected) {
      result.photons.status[index] = PhotonStatus::invalid_input;
      result.summary.add(PhotonStatus::invalid_input, input.weight[index]);
      continue;
    }
    result.photons.position_m[index] = hit->point_m;
    result.photons.direction[index] = *reflected;
    result.photons.optical_path_m[index] = hit->distance_m;
    result.photons.time_ns[index] += hit->distance_m / kSpeedOfLightMPerNs;
    const Ray reflected_ray{hit->point_m, *reflected};
    const auto detector_hit = intersect_detector_surfaces_unchecked(reflected_ray, scene);
    if (!detector_hit) {
      result.photons.status[index] = PhotonStatus::no_detector;
      result.summary.add(PhotonStatus::no_detector, input.weight[index]);
      continue;
    }
    result.photons.position_m[index] = detector_hit->point_m;
    result.photons.optical_path_m[index] += detector_hit->distance_m;
    result.photons.time_ns[index] += detector_hit->distance_m / kSpeedOfLightMPerNs;
    result.photons.weight[index] = input.weight[index];
    result.photons.surface_id[index] = detector_hit->surface_id;
    result.photons.status[index] = PhotonStatus::detected;
    result.summary.add(PhotonStatus::detected, input.weight[index]);
  }
  return result;
}

}  // namespace obdeect
