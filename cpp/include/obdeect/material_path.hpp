#pragma once

#include "obdeect/optical_model.hpp"
#include "obdeect/photon_buffer.hpp"

#include <numbers>

namespace obdeect {

// Diagnostics have fixed storage bounded by the model's validated interaction
// cap. The bulk specialization stores only entrance and terminal points.
template <bool RecordInteractions = true>
[[nodiscard]] inline auto trace_material_path(const Ray &input, std::uint64_t photon_id,
                                              double wavelength_nm,
                                              const CompiledOpticalModel &model) {
  BasicPathRecord<RecordInteractions ? kMaximumOpticalModelInteractions + 1 : 2> path{};
  path.material_transport = true;
  path.photon_id = photon_id;
  path.wavelength_nm = wavelength_nm;
  path.points_m[0] = input.position_m;
  path.point_count = 1;
  const auto direction = normalised_checked(input.direction);
  if (!direction || !std::isfinite(input.position_m.x) || !std::isfinite(input.position_m.y) ||
      !std::isfinite(input.position_m.z) || !std::isfinite(wavelength_nm) || wavelength_nm <= 0) {
    path.status = PhotonStatus::invalid_input;
    return path;
  }
  Ray ray{input.position_m, *direction};
  std::uint32_t medium = model.entrance_medium_id();
  path.final_direction = ray.direction;
  path.status = PhotonStatus::interaction_limit;
  for (std::uint32_t interaction = 0; interaction < model.max_interactions(); ++interaction) {
    const auto current = model.material_at(medium, wavelength_nm);
    if (!current) {
      path.status = PhotonStatus::invalid_input;
      break;
    }
    const auto hit = intersect_nearest_surface(ray, model);
    if (!hit) {
      path.status = medium == kVacuumMediumId ? PhotonStatus::escaped_optical_model
                                              : PhotonStatus::escaped_material;
      break;
    }
    path.path_length_m += hit->distance_m;
    path.optical_path_m += current->phase_index * hit->distance_m;
    path.group_delay_ns += current->group_index * hit->distance_m / kSpeedOfLightMPerNs;
    path.surviving_throughput *= std::exp(-current->absorption_per_m * hit->distance_m);
    const std::size_t vertex = RecordInteractions ? interaction + 1 : 1;
    path.points_m[vertex] = hit->point_m;
    path.point_count = static_cast<std::uint8_t>(vertex + 1);
    path.terminal_surface_id = hit->surface_id;
    const auto &surface = model.surfaces()[hit->surface_index];
    const double incidence =
        std::acos(std::clamp(std::abs(dot(ray.direction, hit->normal)), 0.0, 1.0)) * 180.0 /
        std::numbers::pi;
    OpticalInteractionKind kind = OpticalInteractionKind::mirror;
    if (hit->role == SurfaceRole::detector) {
      path.status = PhotonStatus::detected;
      path.incidence_focal_deg = incidence;
      kind = OpticalInteractionKind::detector;
    } else if (hit->role == SurfaceRole::obscurer) {
      path.status = PhotonStatus::blocked_obscurer;
      kind = OpticalInteractionKind::obscurer;
    } else if (hit->role == SurfaceRole::mirror) {
      const auto reflected = reflect_specular(ray.direction, hit->normal);
      const auto response = surface.reflectivity
                                ? surface.reflectivity->at_unchecked(wavelength_nm, incidence)
                                : std::optional<double>{1};
      if (!reflected || !response)
        path.status = PhotonStatus::invalid_input;
      else {
        path.final_direction = *reflected;
        path.surviving_throughput *= *response;
      }
    } else {
      kind = OpticalInteractionKind::refractive_interface;
      const bool entering_back = dot(ray.direction, hit->normal) < 0;
      const auto incident_medium = entering_back ? surface.front_medium_id : surface.back_medium_id;
      const auto transmitted_medium =
          entering_back ? surface.back_medium_id : surface.front_medium_id;
      const auto next = model.material_at(transmitted_medium, wavelength_nm);
      if (incident_medium != medium || !next) {
        path.status = PhotonStatus::invalid_input;
      } else {
        const Vec3 normal = entering_back ? hit->normal : hit->normal * -1;
        const auto result =
            dielectric_interface(ray.direction, normal, current->phase_index, next->phase_index);
        if (!result) {
          path.status = PhotonStatus::invalid_input;
        } else if (!result->transmitted) {
          // Total internal reflection carries all power and remains in the
          // incident medium. It participates in the same bounded transport.
          path.final_direction = result->reflected;
        } else {
          const auto response = surface.transmission
                                    ? surface.transmission->at_unchecked(wavelength_nm, incidence)
                                    : std::optional<double>{1};
          if (!response) {
            path.status = PhotonStatus::invalid_input;
          } else {
            const double bare = 1 - (result->reflectance_s + result->reflectance_p) / 2;
            const double interfaces =
                surface.transmission_semantics == InterfaceTransmissionSemantics::complete_interface
                    ? 1
                    : bare;
            path.surviving_throughput *= interfaces * *response;
            path.final_direction = *result->transmitted;
            medium = transmitted_medium;
          }
        }
      }
    }
    if constexpr (RecordInteractions) {
      path.interaction_surface_ids[interaction] = hit->surface_id;
      path.interaction_kinds[interaction] = kind;
      path.interaction_normals[interaction] = hit->normal;
      path.interaction_incoming_directions[interaction] = ray.direction;
      path.interaction_outgoing_directions[interaction] = path.final_direction;
      path.interaction_throughput[interaction] = path.surviving_throughput;
      path.interaction_optical_path_m[interaction] = path.optical_path_m;
      path.interaction_group_delay_ns[interaction] = path.group_delay_ns;
    }
    if (path.surviving_throughput == 0 && path.status == PhotonStatus::interaction_limit)
      path.status = PhotonStatus::absorbed_material;
    if (path.status != PhotonStatus::interaction_limit)
      break;
    ray = {hit->point_m, path.final_direction};
  }
  return path;
}

} // namespace obdeect
