#pragma once

#include <cmath>
#include <optional>

namespace obdeect {

constexpr double kEpsilon = 1e-9;
// Exact speed of light expressed in the time unit used by the public photon
// interfaces. Keep flight-time conversion in one place so all trace kernels
// use the same physical constant.
constexpr double kSpeedOfLightMPerNs = 0.299792458;

struct Vec3 {
  double x{};
  double y{};
  double z{};

  constexpr Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
  constexpr Vec3 operator-(const Vec3& other) const { return {x - other.x, y - other.y, z - other.z}; }
  constexpr Vec3 operator*(double value) const { return {x * value, y * value, z * value}; }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& value) { return std::sqrt(dot(value, value)); }

inline std::optional<Vec3> normalised_checked(const Vec3& value) {
  const double length = norm(value);
  if (!std::isfinite(length) || length <= kEpsilon) return std::nullopt;
  return value * (1.0 / length);
}

struct Ray {
  Vec3 position_m;
  Vec3 direction;  // Required to be unit length at a public kernel boundary.
};

}  // namespace obdeect
