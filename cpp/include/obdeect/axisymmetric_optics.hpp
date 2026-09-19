#pragma once

#include "obdeect/interactions.hpp"
#include "obdeect/intersections.hpp"

#include <array>
#include <cstddef>
#include <limits>

namespace obdeect {

// z(r) = sum_i coefficient_m[i] r^(2 i), expressed in the telescope frame.
// This covers a plane, paraboloid, spherical polynomial approximation and the
// even-polynomial CTAO SC primary/secondary descriptions.  The surface is
// deliberately independent of segmentation and material response.
struct EvenPolynomialSurface {
  std::array<double, 13> coefficient_m{};
  // One metre retains the conventional physical-power form.  SC imported
  // prescriptions may instead state powers of r / radial_scale_m.
  double radial_scale_m{1.0};

  [[nodiscard]] double sag(double radius_m) const {
    if (!std::isfinite(radial_scale_m) || radial_scale_m <= kEpsilon) return std::numeric_limits<double>::quiet_NaN();
    const double normalised_radius = radius_m / radial_scale_m;
    const double radius_squared = normalised_radius * normalised_radius;
    double value = 0.0;
    for (auto iterator = coefficient_m.rbegin(); iterator != coefficient_m.rend(); ++iterator) {
      value = value * radius_squared + *iterator;
    }
    return value;
  }

  [[nodiscard]] double radial_slope(double radius_m) const {
    if (!std::isfinite(radial_scale_m) || radial_scale_m <= kEpsilon) return std::numeric_limits<double>::quiet_NaN();
    const double normalised_radius = radius_m / radial_scale_m;
    const double radius_squared = normalised_radius * normalised_radius;
    double derivative_over_radius = 0.0;
    for (std::size_t index = coefficient_m.size(); index-- > 1;) {
      derivative_over_radius = derivative_over_radius * radius_squared +
                               2.0 * static_cast<double>(index) * coefficient_m[index];
    }
    return radius_m * derivative_over_radius / (radial_scale_m * radial_scale_m);
  }
};

struct AxisymmetricMirror {
  double vertex_z_m{};
  double inner_radius_m{};
  double outer_radius_m{};
  EvenPolynomialSurface surface{};
};

struct AxisymmetricHit {
  double distance_m{};
  Vec3 point_m{};
  Vec3 unit_normal{};
};

[[nodiscard]] inline bool is_valid(const AxisymmetricMirror& mirror) {
  if (!std::isfinite(mirror.surface.radial_scale_m) || mirror.surface.radial_scale_m <= kEpsilon)
    return false;
  for (const double coefficient : mirror.surface.coefficient_m)
    if (!std::isfinite(coefficient)) return false;
  return std::isfinite(mirror.vertex_z_m) && std::isfinite(mirror.inner_radius_m) &&
         std::isfinite(mirror.outer_radius_m) && mirror.inner_radius_m >= 0.0 &&
         mirror.outer_radius_m > kEpsilon && mirror.inner_radius_m < mirror.outer_radius_m;
}

// Newton solve for z(t) - vertex_z - sag(r(t)) = 0.  The caller supplies a
// directed stage and the solve returns only forward aperture-valid roots.
[[nodiscard]] inline std::optional<AxisymmetricHit> intersect_axisymmetric_mirror(
    const Ray& ray, const AxisymmetricMirror& mirror, double minimum_t_m = kEpsilon) {
  if (!is_valid(mirror) || !std::isfinite(minimum_t_m) || minimum_t_m < 0.0) return std::nullopt;
  const auto direction = normalised_checked(ray.direction);
  if (!direction || !std::isfinite(ray.position_m.x) || !std::isfinite(ray.position_m.y) ||
      !std::isfinite(ray.position_m.z)) {
    return std::nullopt;
  }
  if (std::abs(direction->z) < kEpsilon) return std::nullopt;

  double distance_m = (mirror.vertex_z_m + mirror.surface.sag(0.0) - ray.position_m.z) / direction->z;
  if (!std::isfinite(distance_m) || distance_m <= minimum_t_m) return std::nullopt;
  constexpr int maximum_iterations = 24;
  constexpr double residual_tolerance_m = 1.0e-12;
  for (int iteration = 0; iteration < maximum_iterations; ++iteration) {
    const Vec3 point = ray.position_m + *direction * distance_m;
    const double radius_m = std::hypot(point.x, point.y);
    const double residual = point.z - mirror.vertex_z_m - mirror.surface.sag(radius_m);
    if (!std::isfinite(residual)) return std::nullopt;
    if (std::abs(residual) <= residual_tolerance_m) {
      if (radius_m < mirror.inner_radius_m || radius_m > mirror.outer_radius_m) return std::nullopt;
      const double slope = mirror.surface.radial_slope(radius_m);
      const double radial_inverse = radius_m > kEpsilon ? 1.0 / radius_m : 0.0;
      const auto normal = normalised_checked(
          {-slope * point.x * radial_inverse, -slope * point.y * radial_inverse, 1.0});
      if (!normal) return std::nullopt;
      return AxisymmetricHit{distance_m, point, *normal};
    }
    const double current_radius_m = std::hypot(point.x, point.y);
    const double slope = mirror.surface.radial_slope(current_radius_m);
    const double radial_derivative = current_radius_m > kEpsilon
                                         ? (point.x * direction->x + point.y * direction->y) / current_radius_m
                                         : 0.0;
    const double derivative = direction->z - slope * radial_derivative;
    if (!std::isfinite(derivative) || std::abs(derivative) < kEpsilon) return std::nullopt;
    distance_m -= residual / derivative;
    if (!std::isfinite(distance_m) || distance_m <= minimum_t_m) return std::nullopt;
  }
  return std::nullopt;
}

[[nodiscard]] inline std::optional<Ray> reflect(const Ray& incident, const AxisymmetricHit& hit) {
  const auto direction = normalised_checked(incident.direction);
  const auto normal = normalised_checked(hit.unit_normal);
  if (!direction || !normal) return std::nullopt;
  const auto reflected = normalised_checked(*direction - *normal * (2.0 * dot(*direction, *normal)));
  return reflected ? std::optional<Ray>{{hit.point_m, *reflected}} : std::nullopt;
}

[[nodiscard]] inline EvenPolynomialSurface paraboloid(double focal_length_m) {
  EvenPolynomialSurface surface{};
  if (std::isfinite(focal_length_m) && focal_length_m > kEpsilon) {
    surface.coefficient_m[1] = 1.0 / (4.0 * focal_length_m);
  }
  return surface;
}

// CTAO SC records give coefficients in centimetre-based powers: z_cm =
// sum_i a_i r_cm^(2i). Convert once into SI before tracing:
// a_i[SI] = a_i[cm] * 0.01^(1 - 2i).
[[nodiscard]] inline EvenPolynomialSurface centimetre_even_polynomial_to_metres(
    const std::array<double, 13>& coefficient_cm) {
  EvenPolynomialSurface surface{};
  for (std::size_t index = 0; index < coefficient_cm.size(); ++index) {
    surface.coefficient_m[index] = coefficient_cm[index] * std::pow(0.01, 1.0 - 2.0 * index);
  }
  return surface;
}

// A general polynomial can specify centimetre sag coefficients in
// powers of r/R, where R is an explicit reference radius. Keep this distinct
// from centimetre_even_polynomial_to_metres(): conflating the conventions can
// produce physically nonsensical large-radius surfaces.
[[nodiscard]] inline EvenPolynomialSurface centimetre_normalised_radius_polynomial_to_metres(
    const std::array<double, 13>& coefficient_cm, double radial_scale_m) {
  EvenPolynomialSurface surface{};
  if (!std::isfinite(radial_scale_m) || radial_scale_m <= kEpsilon) return surface;
  surface.radial_scale_m = radial_scale_m;
  for (std::size_t index = 0; index < coefficient_cm.size(); ++index) {
    surface.coefficient_m[index] = coefficient_cm[index] * 0.01;
  }
  return surface;
}

// sim_telarray's reference-radius convention is z = R sum a_i (r/R)^(2i).
// R is in metres here. It scales BOTH sag and radius, independently of the
// physical aperture. See sim_config.c's parameter / R^(2i-1) conversion.
[[nodiscard]] inline EvenPolynomialSurface reference_radius_polynomial(
    const std::array<double, 13>& coefficients, double reference_radius_m) {
  EvenPolynomialSurface surface{};
  surface.radial_scale_m = reference_radius_m;
  for (std::size_t i = 0; i < coefficients.size(); ++i)
    surface.coefficient_m[i] = coefficients[i] * reference_radius_m;
  return surface;
}

}  // namespace obdeect
