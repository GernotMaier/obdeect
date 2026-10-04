#pragma once

#include "obdeect/axisymmetric_optics.hpp"

#include <cstdint>
#include <numbers>
#include <optional>
#include <span>

namespace obdeect {

enum class AxisymmetricSegmentShape { hexagon, annular_sector };
struct AxisymmetricSegment {
  std::uint32_t id{};
  AxisymmetricSegmentShape shape{};
  Vec3 centre_m{};
  double diameter_m{};
  double rotation_rad{};
  double inner_radius_m{};
  double outer_radius_m{};
  double start_rad{};
  double span_rad{};
  double gap_m{};
  bool gap_at_start{};
};

// Ring boundary and gap convention matches sim_telarray/common/sim_imaging.c
// segment_hit case 4. Hexagons are tested in the tangent frame at their centre,
// preserving physical flat-to-flat diameter rather than projecting a circle.
[[nodiscard]] inline bool contains_axisymmetric_segment(const AxisymmetricSegment &segment,
                                                        const AxisymmetricMirror &mirror,
                                                        const Vec3 &point) {
  if (segment.shape == AxisymmetricSegmentShape::annular_sector) {
    const double radius = std::hypot(point.x, point.y);
    double angle = std::atan2(point.y, point.x) - segment.start_rad;
    angle -= 2 * std::numbers::pi * std::floor(angle / (2 * std::numbers::pi));
    const double distance_from_gap =
        segment.gap_at_start ? angle * radius : (segment.span_rad - angle) * radius;
    return radius >= segment.inner_radius_m && radius < segment.outer_radius_m &&
           angle < segment.span_rad && radius >= segment.gap_m &&
           distance_from_gap >= segment.gap_m;
  }
  const double radius = std::hypot(segment.centre_m.x, segment.centre_m.y);
  const double cx = radius > kEpsilon ? segment.centre_m.x / radius : 1;
  const double sx = radius > kEpsilon ? segment.centre_m.y / radius : 0;
  const double slope = mirror.surface.radial_slope(radius);
  const double cosine = 1 / std::sqrt(1 + slope * slope);
  const double sine = slope * cosine;
  const double dx = point.x - segment.centre_m.x, dy = point.y - segment.centre_m.y;
  const double dz = point.z - mirror.vertex_z_m - mirror.surface.sag(radius);
  const double radial = (dx * cx + dy * sx) * cosine + dz * sine;
  const double azimuthal = -dx * sx + dy * cx;
  const double x = radial * cx - azimuthal * sx, y = radial * sx + azimuthal * cx;
  const double u = x * std::cos(segment.rotation_rad) + y * std::sin(segment.rotation_rad);
  const double v = -x * std::sin(segment.rotation_rad) + y * std::cos(segment.rotation_rad);
  const double half = segment.diameter_m / 2;
  return std::abs(u) <= half && std::abs(0.5 * u + std::sqrt(3.0) / 2 * v) <= half &&
         std::abs(0.5 * u - std::sqrt(3.0) / 2 * v) <= half;
}

[[nodiscard]] inline std::optional<std::uint32_t>
axisymmetric_segment_id(std::span<const AxisymmetricSegment> segments,
                        const AxisymmetricMirror &mirror, const Vec3 &point,
                        std::uint32_t continuous_surface_id) {
  if (segments.empty())
    return continuous_surface_id;
  for (const auto &segment : segments)
    if (contains_axisymmetric_segment(segment, mirror, point))
      return segment.id;
  return std::nullopt;
}

[[nodiscard]] inline std::optional<AxisymmetricHit>
intersect_segmented_asphere(const Ray &ray, const AxisymmetricMirror &mirror,
                            std::span<const AxisymmetricSegment> segments, std::uint32_t surface_id,
                            bool *numerical_failure = nullptr) {
  return intersect_axisymmetric_mirror(
      ray, mirror, kEpsilon,
      [&](const Vec3 &point) {
        return axisymmetric_segment_id(segments, mirror, point, surface_id).has_value();
      },
      numerical_failure);
}

} // namespace obdeect
