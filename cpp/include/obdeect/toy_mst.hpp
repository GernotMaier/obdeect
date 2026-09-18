#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <string_view>
#include <vector>

namespace obdeect {

constexpr double kEpsilon = 1e-9;

struct Vec3 {
  double x{};
  double y{};
  double z{};

  constexpr Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
  constexpr Vec3 operator-(const Vec3& other) const { return {x - other.x, y - other.y, z - other.z}; }
  constexpr Vec3 operator*(double value) const { return {x * value, y * value, z * value}; }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3& value) { return std::sqrt(dot(value, value)); }
inline Vec3 normalised(const Vec3& value) { return value * (1.0 / norm(value)); }

struct Ray {
  Vec3 position_m;
  Vec3 direction;  // Must be unit length.
};

enum class PhotonStatus : std::uint8_t {
  detected,
  blocked_camera,
  blocked_mast,
  missed_primary,
  missed_screen,
};

inline std::string_view to_string(PhotonStatus status) {
  switch (status) {
    case PhotonStatus::detected: return "detected";
    case PhotonStatus::blocked_camera: return "blocked_camera";
    case PhotonStatus::blocked_mast: return "blocked_mast";
    case PhotonStatus::missed_primary: return "missed_primary";
    case PhotonStatus::missed_screen: return "missed_screen";
  }
  return "unknown";
}

struct PathRecord {
  std::uint64_t photon_id{};
  double wavelength_nm{400.0};
  PhotonStatus status{PhotonStatus::missed_primary};
  std::array<Vec3, 3> points_m{};
  std::uint8_t point_count{};
  double path_length_m{};
};

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

inline std::optional<double> intersect_sphere(const Ray& ray, const Vec3& center, double radius) {
  const Vec3 oc = ray.position_m - center;
  const double b = dot(oc, ray.direction);
  const double c = dot(oc, oc) - radius * radius;
  const double discriminant = b * b - c;
  if (discriminant < 0.0) return std::nullopt;
  const double root = std::sqrt(std::max(0.0, discriminant));
  const double first = -b - root;
  const double second = -b + root;
  if (first > kEpsilon) return first;
  if (second > kEpsilon) return second;
  return std::nullopt;
}

// The primary is the lower spherical cap (vertex z=0), not a closed sphere.
// Selecting this cap explicitly prevents an incoming ray from spuriously
// hitting the mathematically valid but physically absent upper hemisphere.
inline std::optional<double> intersect_lower_spherical_cap(const Ray& ray, const Vec3& center, double radius) {
  const Vec3 oc = ray.position_m - center;
  const double b = dot(oc, ray.direction);
  const double c = dot(oc, oc) - radius * radius;
  const double discriminant = b * b - c;
  if (discriminant < 0.0) return std::nullopt;
  const double root = std::sqrt(std::max(0.0, discriminant));
  for (const double t : {-b - root, -b + root}) {
    if (t > kEpsilon && (ray.position_m + ray.direction * t).z <= center.z + kEpsilon) return t;
  }
  return std::nullopt;
}

inline std::optional<double> intersect_plane_z(const Ray& ray, double z) {
  if (std::abs(ray.direction.z) < kEpsilon) return std::nullopt;
  const double t = (z - ray.position_m.z) / ray.direction.z;
  return t > kEpsilon ? std::optional<double>{t} : std::nullopt;
}

inline std::optional<double> intersect_disk_z(const Ray& ray, double z, double radius) {
  const auto t = intersect_plane_z(ray, z);
  if (!t) return std::nullopt;
  const Vec3 point = ray.position_m + ray.direction * *t;
  return point.x * point.x + point.y * point.y <= radius * radius ? t : std::nullopt;
}

// Side intersection of a finite circular cylinder from a to b.  The support
// legs have physical finite extent; caps are intentionally omitted because a
// grazing ray should be tested by the next stage rather than a zero-area cap.
inline std::optional<double> intersect_finite_cylinder(const Ray& ray, const Vec3& a, const Vec3& b,
                                                        double radius) {
  const Vec3 axis = b - a;
  const Vec3 offset = ray.position_m - a;
  const double axis2 = dot(axis, axis);
  const double d_axis = dot(ray.direction, axis);
  const double o_axis = dot(offset, axis);
  const double A = dot(ray.direction, ray.direction) - d_axis * d_axis / axis2;
  const double B = dot(ray.direction, offset) - d_axis * o_axis / axis2;
  const double C = dot(offset, offset) - o_axis * o_axis / axis2 - radius * radius;
  const double discriminant = B * B - A * C;
  if (A <= kEpsilon || discriminant < 0.0) return std::nullopt;
  const double root = std::sqrt(std::max(0.0, discriminant));
  for (const double t : {(-B - root) / A, (-B + root) / A}) {
    const double along = o_axis + t * d_axis;
    if (t > kEpsilon && along >= 0.0 && along <= axis2) return t;
  }
  return std::nullopt;
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

  Ray ray{input.position_m, normalised(input.direction)};
  if (config.include_structure) {
    const double camera_front_z = config.focal_length_m + config.camera_half_depth_m;
    if (const auto t = intersect_disk_z(ray, camera_front_z, config.camera_radius_m)) {
      record.status = PhotonStatus::blocked_camera;
      record.points_m[1] = ray.position_m + ray.direction * *t;
      record.point_count = 2;
      record.path_length_m = *t;
      return record;
    }
    for (const auto& [a, b] : mast_legs(config)) {
      if (const auto t = intersect_finite_cylinder(ray, a, b, config.mast_radius_m)) {
        record.status = PhotonStatus::blocked_mast;
        record.points_m[1] = ray.position_m + ray.direction * *t;
        record.point_count = 2;
        record.path_length_m = *t;
        return record;
      }
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

  const Vec3 normal = normalised(mirror_hit - mirror_center);
  ray = {mirror_hit, normalised(ray.direction - normal * (2.0 * dot(ray.direction, normal)))};
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
