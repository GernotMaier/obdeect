#pragma once

#include "obdeect/interactions.hpp"
#include "obdeect/photon_buffer.hpp"
#include "obdeect/segmented_optical_model.hpp"

#include <numbers>

namespace obdeect {

[[nodiscard]] inline Vec3 cylinder_normal(const CompiledSegmentedOpticalModel &model,
                                          std::uint32_t id, const Vec3 &point) {
  for (const auto &cylinder : model.cylinder_obscurers) {
    if (cylinder.id != id)
      continue;
    const Vec3 delta = cylinder.second_endpoint_m - cylinder.first_endpoint_m;
    const auto axis = normalised_checked(delta);
    if (!axis)
      return {};
    const double projection = dot(point - cylinder.first_endpoint_m, *axis);
    if (projection <= kEpsilon)
      return *axis * -1;
    if (projection >= norm(delta) - kEpsilon)
      return *axis;
    return normalised_checked(point - cylinder.first_endpoint_m - *axis * projection)
        .value_or(Vec3{});
  }
  return {};
}

// Scalar transport shared by the native command and the bulk library runner.
// The immutable model is validated before this allocation-free photon kernel.
[[nodiscard]] inline PathRecord
trace_segmented_path(const Ray &input, std::uint64_t photon_id, double wavelength_nm,
                     const CompiledSegmentedOpticalModel &optical_model) {
  PathRecord path{};
  path.photon_id = photon_id;
  path.wavelength_nm = wavelength_nm;
  path.points_m[0] = input.position_m;
  path.point_count = 1;

  const auto direction = obdeect::normalised_checked(input.direction);
  if (!direction || !std::isfinite(input.position_m.x) || !std::isfinite(input.position_m.y) ||
      !std::isfinite(input.position_m.z) || !std::isfinite(wavelength_nm) || wavelength_nm <= 0) {
    path.status = obdeect::PhotonStatus::invalid_input;
  } else {
    const obdeect::Ray ray{input.position_m, *direction};
    const auto incoming_obscurer =
        obdeect::intersect_cylinder_obscurers_unchecked(ray, optical_model);
    const auto primary_hit = obdeect::intersect_segmented_primary_unchecked(ray, optical_model);
    if (incoming_obscurer &&
        (!primary_hit || incoming_obscurer->distance_m < primary_hit->distance_m)) {
      path.points_m[1] = ray.position_m + ray.direction * incoming_obscurer->distance_m;
      path.point_count = 2;
      path.path_length_m = incoming_obscurer->distance_m;
      path.terminal_surface_id = incoming_obscurer->surface_id;
      path.interaction_surface_ids[0] = incoming_obscurer->surface_id;
      path.interaction_kinds[0] = OpticalInteractionKind::obscurer;
      path.interaction_normals[0] =
          cylinder_normal(optical_model, incoming_obscurer->surface_id, path.points_m[1]);
      path.interaction_incoming_directions[0] = *direction;
      path.interaction_outgoing_directions[0] = *direction;
      path.status = obdeect::PhotonStatus::blocked_obscurer;
      path.final_direction = *direction;
    } else if (!primary_hit) {
      path.status = obdeect::PhotonStatus::missed_primary;
      path.final_direction = *direction;
    } else {
      path.points_m[1] = primary_hit->point_m;
      path.terminal_surface_id = primary_hit->facet_id;
      path.interaction_surface_ids[0] = primary_hit->facet_id;
      path.interaction_normals[0] = primary_hit->unit_normal;
      path.interaction_incoming_directions[0] = *direction;
      path.point_count = 2;
      path.path_length_m = primary_hit->distance_m;
      path.incidence_primary_deg =
          std::acos(
              std::clamp(std::abs(obdeect::dot(*direction, primary_hit->unit_normal)), 0.0, 1.0)) *
          180.0 / std::numbers::pi;
      const auto reflected = obdeect::reflect_specular(*direction, primary_hit->unit_normal);
      if (!reflected) {
        path.status = obdeect::PhotonStatus::invalid_input;
      } else {
        path.interaction_outgoing_directions[0] = *reflected;
        path.final_direction = *reflected;
        const auto response = optical_model.primary_reflectivity
                                  ? optical_model.primary_reflectivity->at_unchecked(
                                        wavelength_nm, path.incidence_primary_deg)
                                  : std::optional<double>{1.0};
        if (!response) {
          path.status = PhotonStatus::invalid_input;
          return path;
        }
        path.surviving_throughput = *response;
        path.interaction_throughput[0] = *response;
        const auto detector_hit = obdeect::intersect_detector_surfaces_unchecked(
            obdeect::Ray{primary_hit->point_m, *reflected}, optical_model);
        const obdeect::Ray reflected_ray{primary_hit->point_m, *reflected};
        const auto outgoing_obscurer =
            obdeect::intersect_cylinder_obscurers_unchecked(reflected_ray, optical_model);
        if (outgoing_obscurer &&
            (!detector_hit || outgoing_obscurer->distance_m < detector_hit->distance_m)) {
          path.points_m[2] =
              reflected_ray.position_m + reflected_ray.direction * outgoing_obscurer->distance_m;
          path.point_count = 3;
          path.terminal_surface_id = outgoing_obscurer->surface_id;
          path.interaction_surface_ids[1] = outgoing_obscurer->surface_id;
          path.interaction_kinds[1] = OpticalInteractionKind::obscurer;
          path.interaction_normals[1] =
              cylinder_normal(optical_model, outgoing_obscurer->surface_id, path.points_m[2]);
          path.interaction_incoming_directions[1] = *reflected;
          path.interaction_outgoing_directions[1] = *reflected;
          path.interaction_throughput[1] = path.surviving_throughput;
          path.status = obdeect::PhotonStatus::blocked_obscurer;
          path.path_length_m += outgoing_obscurer->distance_m;
          path.final_direction = *reflected;
        } else if (!detector_hit) {
          path.status = obdeect::PhotonStatus::no_detector;
          path.final_direction = *reflected;
        } else {
          path.points_m[2] = detector_hit->point_m;
          path.terminal_surface_id = detector_hit->surface_id;
          path.interaction_surface_ids[1] = detector_hit->surface_id;
          path.interaction_kinds[1] = OpticalInteractionKind::detector;
          path.interaction_normals[1] = detector_hit->unit_normal;
          path.interaction_incoming_directions[1] = *reflected;
          path.interaction_outgoing_directions[1] = *reflected;
          path.interaction_throughput[1] = path.surviving_throughput;
          path.point_count = 3;
          path.path_length_m += detector_hit->distance_m;
          path.final_direction = *reflected;
          path.incidence_focal_deg =
              std::acos(std::clamp(std::abs(obdeect::dot(*reflected, detector_hit->unit_normal)),
                                   0.0, 1.0)) *
              180.0 / std::numbers::pi;
          path.status = obdeect::PhotonStatus::detected;
        }
      }
    }
  }

  return path;
}

} // namespace obdeect
