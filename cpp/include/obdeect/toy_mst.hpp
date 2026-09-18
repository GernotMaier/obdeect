#pragma once

#include "obdeect/interactions.hpp"
#include "obdeect/intersections.hpp"
#include "obdeect/photon_buffer.hpp"

#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

namespace obdeect {

// An intentionally small MST-inspired optical model.  It is not a CTAO
// production model: it has one continuous spherical mirror, a camera disk and
// four finite cylindrical mast legs.  Its purpose is a clear executable first
// vertical slice, including real pre-M1 structural shadowing.
struct ToyMstConfig {
  double mirror_radius_m{9.75};       // MST-like 12 m dish curvature scale.
  double mirror_aperture_radius_m{6.0};
  double focal_length_m{4.875};       // Spherical paraxial focus = R / 2.
  double screen_radius_m{2.0};
  double camera_radius_m{0.55};
  double camera_half_depth_m{0.25};
  double mast_radius_m{0.075};
  bool include_structure{true};
};

inline bool is_valid(const ToyMstConfig& config) {
  return std::isfinite(config.mirror_radius_m) && std::isfinite(config.mirror_aperture_radius_m) &&
         std::isfinite(config.focal_length_m) && std::isfinite(config.screen_radius_m) &&
         std::isfinite(config.camera_radius_m) && std::isfinite(config.camera_half_depth_m) &&
         std::isfinite(config.mast_radius_m) && config.mirror_radius_m > kEpsilon &&
         config.mirror_aperture_radius_m > kEpsilon &&
         config.mirror_aperture_radius_m <= config.mirror_radius_m && config.focal_length_m > kEpsilon &&
         config.screen_radius_m > kEpsilon && config.camera_radius_m >= 0.0 &&
         config.camera_half_depth_m >= 0.0 && config.camera_half_depth_m < config.focal_length_m &&
         std::isfinite(config.focal_length_m + config.camera_half_depth_m) && config.mast_radius_m >= 0.0;
}


inline std::vector<std::pair<Vec3, Vec3>> mast_legs(const ToyMstConfig& config) {
  std::vector<std::pair<Vec3, Vec3>> legs;
  constexpr double base_radius_m = 4.2;
  constexpr double camera_support_radius_m = 0.7;
  for (int leg = 0; leg < 4; ++leg) {
    const double angle = static_cast<double>(leg) * std::numbers::pi / 2.0;
    legs.push_back({{base_radius_m * std::cos(angle), base_radius_m * std::sin(angle), 0.30},
                    {camera_support_radius_m * std::cos(angle), camera_support_radius_m * std::sin(angle),
                     config.focal_length_m - config.camera_half_depth_m}});
  }
  return legs;
}

inline PathRecord trace_toy_mst(const Ray& input, std::uint64_t photon_id, const ToyMstConfig& config) {
  PathRecord record{};
  record.photon_id = photon_id;
  record.points_m[0] = input.position_m;
  record.point_count = 1;

  const auto direction = normalised_checked(input.direction);
  if (!is_valid(config) || !direction || !std::isfinite(input.position_m.x) ||
      !std::isfinite(input.position_m.y) || !std::isfinite(input.position_m.z)) {
    record.status = PhotonStatus::invalid_input;
    return record;
  }
  Ray ray{input.position_m, *direction};
  if (config.include_structure) {
    std::optional<std::pair<double, PhotonStatus>> nearest_obstruction;
    const auto consider_obstruction = [&nearest_obstruction](std::optional<double> candidate,
                                                               PhotonStatus status) {
      if (candidate && (!nearest_obstruction || *candidate < nearest_obstruction->first)) {
        nearest_obstruction = std::make_pair(*candidate, status);
      }
    };
    const double camera_front_z = config.focal_length_m + config.camera_half_depth_m;
    consider_obstruction(intersect_disk_z(ray, camera_front_z, config.camera_radius_m),
                         PhotonStatus::blocked_camera);
    for (const auto& [a, b] : mast_legs(config)) {
      consider_obstruction(intersect_finite_cylinder(ray, a, b, config.mast_radius_m),
                           PhotonStatus::blocked_mast);
    }
    if (nearest_obstruction) {
      record.status = nearest_obstruction->second;
      record.points_m[1] = ray.position_m + ray.direction * nearest_obstruction->first;
      record.point_count = 2;
      record.path_length_m = nearest_obstruction->first;
      return record;
    }
  }

  const Vec3 mirror_center{0.0, 0.0, config.mirror_radius_m};
  const auto mirror_t = intersect_lower_spherical_cap(ray, mirror_center, config.mirror_radius_m);
  if (!mirror_t) {
    record.status = PhotonStatus::missed_primary;
    return record;
  }
  const Vec3 mirror_hit = ray.position_m + ray.direction * *mirror_t;
  if (mirror_hit.x * mirror_hit.x + mirror_hit.y * mirror_hit.y >
      config.mirror_aperture_radius_m * config.mirror_aperture_radius_m) {
    record.status = PhotonStatus::missed_primary;
    record.points_m[1] = mirror_hit;
    record.point_count = 2;
    record.path_length_m = *mirror_t;
    return record;
  }

  const auto normal = normalised_checked(mirror_hit - mirror_center);
  if (!normal) {
    record.status = PhotonStatus::invalid_input;
    return record;
  }
  const auto reflected = normalised_checked(ray.direction - *normal * (2.0 * dot(ray.direction, *normal)));
  if (!reflected) {
    record.status = PhotonStatus::invalid_input;
    return record;
  }
  ray = {mirror_hit, *reflected};
  record.points_m[1] = mirror_hit;
  record.point_count = 2;
  record.path_length_m = *mirror_t;

  const auto screen_t = intersect_plane_z(ray, config.focal_length_m);
  if (!screen_t) {
    record.status = PhotonStatus::missed_screen;
    return record;
  }
  const Vec3 screen_hit = ray.position_m + ray.direction * *screen_t;
  record.points_m[2] = screen_hit;
  record.point_count = 3;
  record.path_length_m += *screen_t;
  record.status = screen_hit.x * screen_hit.x + screen_hit.y * screen_hit.y <=
                          config.screen_radius_m * config.screen_radius_m
                      ? PhotonStatus::detected
                      : PhotonStatus::missed_screen;
  return record;
}

// Low-discrepancy, uniform-area samples across the optical entrance pupil.
inline std::vector<Ray> parallel_blue_cherenkov_rays(std::size_t count, const ToyMstConfig& config,
                                                      double source_z_m = 20.0) {
  std::vector<Ray> rays;
  if (!is_valid(config) || !std::isfinite(source_z_m) || count == 0) return rays;
  rays.reserve(count);
  constexpr double golden_ratio_conjugate = 0.6180339887498948482;
  for (std::size_t index = 0; index < count; ++index) {
    const double u = (static_cast<double>(index) + 0.5) / static_cast<double>(count);
    const double r = config.mirror_aperture_radius_m * std::sqrt(u);
    const double phi = std::fmod(static_cast<double>(index) * golden_ratio_conjugate, 1.0) *
                       2.0 * std::numbers::pi;
    rays.push_back({{r * std::cos(phi), r * std::sin(phi), source_z_m}, {0.0, 0.0, -1.0}});
  }
  return rays;
}

}  // namespace obdeect
