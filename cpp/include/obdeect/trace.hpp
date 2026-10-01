#pragma once

#include "obdeect/abi.hpp"
#include "obdeect/diagnostics.hpp"
#include "obdeect/reference_optical_model.hpp"
#include "obdeect/optical_model.hpp"
#include "obdeect/segmented_optical_model.hpp"

namespace obdeect {

struct TraceResult {
  PhotonResultBlock photons;
  TraceSummary summary;
};

// Scalar reference block runner. It is intentionally allocation-free inside
// the photon loop: output is sized before tracing and each record is written by
// index. SIMD/threaded kernels must reproduce this contract exactly.
[[nodiscard]] inline TraceResult trace(const CompiledReferenceOpticalModel& optical_model, const PhotonBlockView& input) {
  TraceResult result{PhotonResultBlock{input.position_m.size()}, {}};
  if (!optical_model.is_valid() || !validate_photon_block(input)) {
    for (std::size_t index = 0; index < input.position_m.size(); ++index) {
      result.photons.status[index] = PhotonStatus::invalid_input;
      ++result.summary.status_count[static_cast<std::size_t>(PhotonStatus::invalid_input)];
    }
    return result;
  }
  for (std::size_t index = 0; index < input.position_m.size(); ++index) {
    const PathRecord record = trace_artificial_mst({input.position_m[index], input.direction[index]}, input.photon_id[index],
                                            optical_model.configuration);
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

// A segmented optical_model may supply finite detector surfaces. A detected photon
// records the nearest physical post-reflection intersection; otherwise it
// retains the explicit no-detector terminal status.
[[nodiscard]] inline TraceResult trace(const CompiledSegmentedOpticalModel& optical_model, const PhotonBlockView& input) {
  TraceResult result{PhotonResultBlock{input.position_m.size()}, {}};
  if (!is_valid(optical_model) || !validate_photon_block(input)) {
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
    const auto incoming_obscurer = intersect_cylinder_obscurers_unchecked(ray, optical_model);
    const auto hit = intersect_segmented_primary_unchecked(ray, optical_model);
    if (incoming_obscurer && (!hit || incoming_obscurer->distance_m < hit->distance_m)) {
      result.photons.position_m[index] =
          ray.position_m + ray.direction * incoming_obscurer->distance_m;
      result.photons.optical_path_m[index] = incoming_obscurer->distance_m;
      result.photons.time_ns[index] += incoming_obscurer->distance_m / kSpeedOfLightMPerNs;
      result.photons.status[index] = PhotonStatus::blocked_obscurer;
      result.photons.surface_id[index] = incoming_obscurer->surface_id;
      result.summary.add(PhotonStatus::blocked_obscurer, input.weight[index]);
      continue;
    }
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
    const auto outgoing_obscurer = intersect_cylinder_obscurers_unchecked(reflected_ray, optical_model);
    const auto detector_hit = intersect_detector_surfaces_unchecked(reflected_ray, optical_model);
    if (outgoing_obscurer && (!detector_hit || outgoing_obscurer->distance_m < detector_hit->distance_m)) {
      result.photons.position_m[index] =
          reflected_ray.position_m + reflected_ray.direction * outgoing_obscurer->distance_m;
      result.photons.optical_path_m[index] += outgoing_obscurer->distance_m;
      result.photons.time_ns[index] += outgoing_obscurer->distance_m / kSpeedOfLightMPerNs;
      result.photons.surface_id[index] = outgoing_obscurer->surface_id;
      result.photons.status[index] = PhotonStatus::blocked_obscurer;
      result.summary.add(PhotonStatus::blocked_obscurer, input.weight[index]);
      continue;
    }
    if (!detector_hit) {
      result.photons.status[index] = PhotonStatus::no_detector;
      result.summary.add(PhotonStatus::no_detector, input.weight[index]);
      continue;
    }
    result.photons.position_m[index] = detector_hit->point_m;
    result.photons.optical_path_m[index] += detector_hit->distance_m;
    result.photons.time_ns[index] += detector_hit->distance_m / kSpeedOfLightMPerNs;
    const auto reflectivity = optical_model.primary_reflectivity
                                  ? optical_model.primary_reflectivity->at(input.wavelength_nm[index])
                                  : std::optional<double>{1.0};
    if (!reflectivity) {
      result.photons.status[index] = PhotonStatus::invalid_input;
      result.summary.add(PhotonStatus::invalid_input, input.weight[index]);
      continue;
    }
    result.photons.weight[index] = input.weight[index] * *reflectivity;
    result.photons.surface_id[index] = detector_hit->surface_id;
    result.photons.status[index] = PhotonStatus::detected;
    result.summary.add(PhotonStatus::detected, input.weight[index]);
  }
  return result;
}

// Geometry-only non-sequential reference. Every segment considers every
// surface; material_id is retained for later response binding.
[[nodiscard]] inline TraceResult trace(const CompiledOpticalModel& optical_model, const PhotonBlockView& input) {
  TraceResult result{PhotonResultBlock{input.position_m.size()}, {}};
  if (!validate_photon_block(input)) {
    for (std::size_t index = 0; index < input.position_m.size(); ++index) {
      result.photons.status[index] = PhotonStatus::invalid_input;
      ++result.summary.status_count[static_cast<std::size_t>(PhotonStatus::invalid_input)];
    }
    return result;
  }
  for (std::size_t index = 0; index < input.position_m.size(); ++index) {
    Ray ray{input.position_m[index], input.direction[index]};
    result.photons.position_m[index] = ray.position_m;
    result.photons.direction[index] = ray.direction;
    result.photons.time_ns[index] = input.time_ns[index];
    PhotonStatus status = PhotonStatus::interaction_limit;
    for (std::uint32_t interaction = 0; interaction < optical_model.max_interactions(); ++interaction) {
      const auto hit = intersect_nearest_surface(ray, optical_model);
      if (!hit) {
        status = PhotonStatus::escaped_optical_model;
        break;
      }
      result.photons.optical_path_m[index] += hit->distance_m;
      result.photons.time_ns[index] += hit->distance_m / kSpeedOfLightMPerNs;
      result.photons.position_m[index] = hit->point_m;
      if (hit->role == SurfaceRole::detector) {
        status = PhotonStatus::detected;
        result.photons.surface_id[index] = hit->surface_id;
        break;
      }
      if (hit->role == SurfaceRole::obscurer) {
        status = PhotonStatus::blocked_obscurer;
        result.photons.surface_id[index] = hit->surface_id;
        break;
      }
      const auto reflected = reflect_specular(ray.direction, hit->normal);
      if (!reflected) {
        status = PhotonStatus::invalid_input;
        break;
      }
      result.photons.direction[index] = *reflected;
      ray = {hit->point_m, *reflected};
      if (interaction + 1 == optical_model.max_interactions())
        result.photons.surface_id[index] = hit->surface_id;
    }
    result.photons.status[index] = status;
    result.photons.weight[index] = status == PhotonStatus::detected ? input.weight[index] : 0.0;
    result.summary.add(status, input.weight[index]);
  }
  return result;
}

}  // namespace obdeect
