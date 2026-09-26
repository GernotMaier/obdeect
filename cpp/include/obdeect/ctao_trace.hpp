#pragma once

#include "obdeect/ctao_models.hpp"
#include "obdeect/photon_buffer.hpp"

#include <numbers>

namespace obdeect {

[[nodiscard]] inline std::optional<AxisymmetricHit> intersect_primary(const Ray& ray,
                                                                        const CtaoReferenceModel& model) {
  if (model.primary_sphere_radius_m) {
    const Vec3 centre{0.0, 0.0, model.primary_vertex_z_m + *model.primary_sphere_radius_m};
    const auto distance = intersect_lower_spherical_cap(ray, centre, *model.primary_sphere_radius_m);
    if (!distance) return std::nullopt;
    const Vec3 point = ray.position_m + *normalised_checked(ray.direction) * *distance;
    const double radius = std::hypot(point.x, point.y);
    if (radius < model.primary_inner_radius_m || radius > model.primary_outer_radius_m) return std::nullopt;
    const auto normal = normalised_checked(point - centre);
    return normal ? std::optional<AxisymmetricHit>{{*distance, point, *normal}} : std::nullopt;
  }
  return intersect_axisymmetric_mirror(
      ray, {model.primary_vertex_z_m, model.primary_inner_radius_m, model.primary_outer_radius_m,
            model.primary_surface});
}

// Reference optical-chain tracer. It is intentionally geometric only: no
// wavelength response, facet boundaries, obscurations or detector conversion are
// silently invented. Such features are applied only after model import.
[[nodiscard]] inline PathRecord trace_ctao_reference(const Ray& input, std::uint64_t photon_id,
                                                      const CtaoReferenceModel& model) {
  PathRecord record{};
  record.photon_id = photon_id;
  record.points_m[0] = input.position_m;
  record.point_count = 1;
  const auto unit_direction = normalised_checked(input.direction);
  if (!unit_direction || !is_valid(model) || !std::isfinite(input.position_m.x) ||
      !std::isfinite(input.position_m.y) || !std::isfinite(input.position_m.z)) {
    record.status = PhotonStatus::invalid_input;
    return record;
  }
  Ray ray{input.position_m, *unit_direction};
  const auto first_hit = intersect_primary(ray, model);
  if (!first_hit) {
    record.status = PhotonStatus::missed_primary;
    record.final_direction = ray.direction;
    return record;
  }
  record.points_m[1] = first_hit->point_m;
  record.point_count = 2;
  record.path_length_m = first_hit->distance_m;
  record.incidence_primary_deg =
      std::acos(std::clamp(std::abs(dot(ray.direction, first_hit->unit_normal)), 0.0, 1.0)) *
      180.0 / std::numbers::pi;
  const auto after_primary = reflect(ray, *first_hit);
  if (!after_primary) {
    record.status = PhotonStatus::invalid_input;
    return record;
  }
  ray = *after_primary;
  if (model.secondary) {
    const auto second_hit = intersect_axisymmetric_mirror(ray, *model.secondary);
    if (!second_hit) {
      record.status = PhotonStatus::missed_screen;
      record.final_direction = ray.direction;
      return record;
    }
    record.points_m[2] = second_hit->point_m;
    record.point_count = 3;
    record.path_length_m += second_hit->distance_m;
    record.incidence_secondary_deg =
        std::acos(std::clamp(std::abs(dot(ray.direction, second_hit->unit_normal)), 0.0, 1.0)) *
        180.0 / std::numbers::pi;
    const auto after_secondary = reflect(ray, *second_hit);
    if (!after_secondary) {
      record.status = PhotonStatus::invalid_input;
      return record;
    }
    ray = *after_secondary;
  }
  const auto focal_distance = intersect_plane_z(ray, model.focal_plane_z_m);
  if (!focal_distance) {
    record.status = PhotonStatus::missed_screen;
    record.final_direction = ray.direction;
    return record;
  }
  const Vec3 focal_point = ray.position_m + ray.direction * *focal_distance;
  record.points_m[model.secondary ? 3 : 2] = focal_point;
  record.point_count = model.secondary ? 4 : 3;
  record.path_length_m += *focal_distance;
  record.incidence_focal_deg =
      std::acos(std::clamp(std::abs(ray.direction.z), 0.0, 1.0)) * 180.0 / std::numbers::pi;
  record.final_direction = ray.direction;
  record.status = std::hypot(focal_point.x, focal_point.y) <= model.focal_plane_radius_m
                      ? PhotonStatus::detected
                      : PhotonStatus::missed_screen;
  return record;
}

}  // namespace obdeect
