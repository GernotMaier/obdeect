#pragma once

#include "obdeect/detector_assignment.hpp"
#include "obdeect/detector_planes.hpp"
#include "obdeect/interactions.hpp"
#include "obdeect/photon_buffer.hpp"
#include "obdeect/pixel_response.hpp"
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
    if (optical_model.telescope_transmission)
      path.surviving_throughput = *optical_model.telescope_transmission->at_unchecked(*direction);
    const obdeect::Ray ray{input.position_m, *direction};
    auto incoming_obscurer = obdeect::intersect_cylinder_obscurers_unchecked(ray, optical_model);
    const auto opaque_incoming = intersect_opaque_surfaces(ray, optical_model.opaque_obscurers);
    if (opaque_incoming &&
        (!incoming_obscurer || opaque_incoming->distance_m < incoming_obscurer->distance_m))
      incoming_obscurer =
          CylinderObscurerHit{opaque_incoming->surface_id, opaque_incoming->distance_m};
    const auto primary_hit = obdeect::intersect_segmented_primary_unchecked(ray, optical_model);
    const auto camera_shadow = optical_model.incoming_obscurer_planes
                                   ? optical_model.incoming_obscurer_planes->intersect(ray)
                                   : std::optional<DetectorSurfaceHit>{};
    if (camera_shadow && (!primary_hit || camera_shadow->distance_m < primary_hit->distance_m) &&
        (!incoming_obscurer || camera_shadow->distance_m < incoming_obscurer->distance_m)) {
      path.points_m[1] = camera_shadow->point_m;
      path.point_count = 2;
      path.path_length_m = camera_shadow->distance_m;
      path.terminal_surface_id = camera_shadow->surface_id;
      path.interaction_surface_ids[0] = camera_shadow->surface_id;
      path.interaction_kinds[0] = OpticalInteractionKind::obscurer;
      path.interaction_normals[0] = camera_shadow->unit_normal;
      path.interaction_incoming_directions[0] = *direction;
      path.interaction_outgoing_directions[0] = *direction;
      path.status = PhotonStatus::blocked_obscurer;
      path.final_direction = *direction;
    } else if (incoming_obscurer &&
               (!primary_hit || incoming_obscurer->distance_m < primary_hit->distance_m)) {
      path.points_m[1] = ray.position_m + ray.direction * incoming_obscurer->distance_m;
      path.point_count = 2;
      path.path_length_m = incoming_obscurer->distance_m;
      path.terminal_surface_id = incoming_obscurer->surface_id;
      path.interaction_surface_ids[0] = incoming_obscurer->surface_id;
      path.interaction_kinds[0] = OpticalInteractionKind::obscurer;
      path.interaction_normals[0] =
          opaque_incoming && opaque_incoming->surface_id == incoming_obscurer->surface_id
              ? opaque_incoming->unit_normal
              : cylinder_normal(optical_model, incoming_obscurer->surface_id, path.points_m[1]);
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
      const auto reflected = optical_model.primary_scatter
                                 ? reflect_with_scatter(*direction, primary_hit->unit_normal,
                                                        {1, 0, 0}, *optical_model.primary_scatter,
                                                        photon_id, primary_hit->facet_id)
                                 : obdeect::reflect_specular(*direction, primary_hit->unit_normal);
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
        path.surviving_throughput *= *response;
        if (optical_model.primary_degradation) {
          Vec3 point = primary_hit->point_m;
          if (optical_model.primary_degradation_in_facet_frame)
            for (const auto &facet : optical_model.primary_facets)
              if (facet.id == primary_hit->facet_id) {
                const Vec3 offset = point - facet.centre_m;
                point = {dot(offset, facet.response_basis_u), dot(offset, facet.response_basis_v),
                         0};
                break;
              }
          const auto degradation = optical_model.primary_degradation->at_point_unchecked(point);
          if (!degradation) {
            path.status = PhotonStatus::invalid_input;
            return path;
          }
          path.surviving_throughput *= *degradation;
        }
        path.interaction_throughput[0] = path.surviving_throughput;
        const obdeect::Ray detector_ray{primary_hit->point_m, *reflected};
        std::optional<DetectorSurfaceHit> detector_hit;
        if (optical_model.detector_assignment) {
          const auto z = optical_model.detector_assignment->description().reference_plane_z_m;
          if (!z || std::abs(detector_ray.direction.z) <= kEpsilon)
            detector_hit.reset();
          else {
            const double distance = (*z - detector_ray.position_m.z) / detector_ray.direction.z;
            detector_hit =
                distance > kEpsilon
                    ? optical_model.detector_assignment->intersect(
                          detector_ray, detector_ray.position_m + detector_ray.direction * distance)
                    : std::nullopt;
          }
        } else
          detector_hit =
              optical_model.detector_planes
                  ? optical_model.detector_planes->intersect(detector_ray)
                  : obdeect::intersect_detector_surfaces_unchecked(detector_ray, optical_model);
        const obdeect::Ray reflected_ray{primary_hit->point_m, *reflected};
        auto outgoing_obscurer =
            obdeect::intersect_cylinder_obscurers_unchecked(reflected_ray, optical_model);
        const auto opaque_outgoing =
            intersect_opaque_surfaces(reflected_ray, optical_model.opaque_obscurers);
        if (opaque_outgoing &&
            (!outgoing_obscurer || opaque_outgoing->distance_m < outgoing_obscurer->distance_m))
          outgoing_obscurer =
              CylinderObscurerHit{opaque_outgoing->surface_id, opaque_outgoing->distance_m};
        if (outgoing_obscurer &&
            (!detector_hit || outgoing_obscurer->distance_m < detector_hit->distance_m)) {
          path.points_m[2] =
              reflected_ray.position_m + reflected_ray.direction * outgoing_obscurer->distance_m;
          path.point_count = 3;
          path.terminal_surface_id = outgoing_obscurer->surface_id;
          path.interaction_surface_ids[1] = outgoing_obscurer->surface_id;
          path.interaction_kinds[1] = OpticalInteractionKind::obscurer;
          path.interaction_normals[1] =
              opaque_outgoing && opaque_outgoing->surface_id == outgoing_obscurer->surface_id
                  ? opaque_outgoing->unit_normal
                  : cylinder_normal(optical_model, outgoing_obscurer->surface_id, path.points_m[2]);
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
          if (optical_model.camera_degradation) {
            const auto point = optical_model.camera_degradation_in_detector_frame
                                   ? detector_hit->local_position_m
                                   : std::optional{detector_hit->point_m};
            const auto degradation =
                !point ? std::nullopt
                       : optical_model.camera_degradation->at_point_unchecked(*point);
            if (!degradation) {
              path.status = PhotonStatus::invalid_input;
              return path;
            }
            path.surviving_throughput *= *degradation;
          }
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
          if (optical_model.pixel_responses) {
            const auto response = optical_model.pixel_responses->at_unchecked(
                detector_hit->surface_id, wavelength_nm, path.incidence_focal_deg,
                {detector_hit->point_m, *reflected});
            if (!response) {
              path.status = PhotonStatus::invalid_input;
              return path;
            }
            path.surviving_throughput *= *response;
            path.interaction_throughput[1] = path.surviving_throughput;
          }
          path.status = obdeect::PhotonStatus::detected;
          if (optical_model.camera_response) {
            const auto response = optical_model.camera_response->at_unchecked(
                wavelength_nm, path.incidence_focal_deg);
            if (!response)
              path.status = PhotonStatus::invalid_input;
            else
              path.surviving_throughput *= *response;
            path.interaction_throughput[1] = path.surviving_throughput;
          }
        }
      }
    }
  }

  return path;
}

} // namespace obdeect
