#pragma once

#include <cmath>

namespace obdeect {

// Speed of light in vacuum, metres per nanosecond.
constexpr double kSpeedOfLightMPerNs = 0.299792458;
// Smallest accepted propagation distance; rejects self-intersections at the
// point a ray left.
constexpr double kMinDistanceM = 1e-9;

struct Vec3 {
  double x{}, y{}, z{};

  constexpr Vec3 operator+(const Vec3 &o) const { return {x + o.x, y + o.y, z + o.z}; }
  constexpr Vec3 operator-(const Vec3 &o) const { return {x - o.x, y - o.y, z - o.z}; }
  constexpr Vec3 operator-() const { return {-x, -y, -z}; }
  constexpr Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};

constexpr double dot(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
constexpr Vec3 cross(const Vec3 &a, const Vec3 &b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3 &a) { return std::sqrt(dot(a, a)); }
inline Vec3 normalized(const Vec3 &a) { return a * (1.0 / norm(a)); }
inline bool finite(const Vec3 &a) {
  return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}

// Ray with a unit direction.
struct Ray {
  Vec3 origin;
  Vec3 direction;
  [[nodiscard]] Vec3 at(double t) const { return origin + direction * t; }
};

// Specular reflection of a unit direction about a unit normal.
inline Vec3 reflect(const Vec3 &d, const Vec3 &n) { return d - n * (2.0 * dot(d, n)); }

} // namespace obdeect
