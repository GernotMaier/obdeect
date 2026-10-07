#pragma once

#include "obdeect/interactions.hpp"
#include "obdeect/intersections.hpp"
#include "obdeect/polynomial_roots.hpp"

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
    if (!std::isfinite(radial_scale_m) || radial_scale_m <= kEpsilon)
      return std::numeric_limits<double>::quiet_NaN();
    const double normalised_radius = radius_m / radial_scale_m;
    const double radius_squared = normalised_radius * normalised_radius;
    double value = 0.0;
    for (auto iterator = coefficient_m.rbegin(); iterator != coefficient_m.rend(); ++iterator) {
      value = value * radius_squared + *iterator;
    }
    return value;
  }

  [[nodiscard]] double radial_slope(double radius_m) const {
    if (!std::isfinite(radial_scale_m) || radial_scale_m <= kEpsilon)
      return std::numeric_limits<double>::quiet_NaN();
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

[[nodiscard]] inline bool is_valid(const AxisymmetricMirror &mirror) {
  if (!std::isfinite(mirror.surface.radial_scale_m) || mirror.surface.radial_scale_m <= kEpsilon)
    return false;
  for (const double coefficient : mirror.surface.coefficient_m)
    if (!std::isfinite(coefficient))
      return false;
  return std::isfinite(mirror.vertex_z_m) && std::isfinite(mirror.inner_radius_m) &&
         std::isfinite(mirror.outer_radius_m) && mirror.inner_radius_m >= 0.0 &&
         mirror.outer_radius_m > kEpsilon && mirror.inner_radius_m < mirror.outer_radius_m;
}

// Restrict the directed ray to the finite aperture cylinder, then isolate all
// roots of its polynomial sag equation on that interval. This selects the
// nearest admissible root, including horizontal rays and grazing tangencies.
template <class Accept>
[[nodiscard]] inline std::optional<AxisymmetricHit>
intersect_axisymmetric_mirror(const Ray &ray, const AxisymmetricMirror &mirror, double minimum_t_m,
                              Accept accept, bool *numerical_failure = nullptr) {
  if (numerical_failure)
    *numerical_failure = false;
  if (!is_valid(mirror) || !std::isfinite(minimum_t_m) || minimum_t_m < 0.0)
    return std::nullopt;
  const auto direction = normalised_checked(ray.direction);
  if (!direction || !std::isfinite(ray.position_m.x) || !std::isfinite(ray.position_m.y) ||
      !std::isfinite(ray.position_m.z))
    return std::nullopt;
  const auto make_hit = [&](double distance) -> std::optional<AxisymmetricHit> {
    if (!std::isfinite(distance)) {
      if (numerical_failure)
        *numerical_failure = true;
      return std::nullopt;
    }
    if (distance <= minimum_t_m)
      return std::nullopt;
    const Vec3 point = ray.position_m + *direction * distance;
    const double radius = std::hypot(point.x, point.y);
    if (radius < mirror.inner_radius_m || radius > mirror.outer_radius_m)
      return std::nullopt;
    if (!accept(point))
      return std::nullopt;
    const double sag = mirror.vertex_z_m + mirror.surface.sag(radius);
    const double residual_tolerance = 1.e-10 * std::max(1.0, std::abs(sag));
    if (!std::isfinite(sag) || std::abs(point.z - sag) > residual_tolerance) {
      if (numerical_failure)
        *numerical_failure = true;
      return std::nullopt;
    }
    const double slope = mirror.surface.radial_slope(radius);
    const double inverse = radius > kEpsilon ? 1.0 / radius : 0.0;
    const auto normal =
        normalised_checked({-slope * point.x * inverse, -slope * point.y * inverse, 1.0});
    if (!normal && numerical_failure)
      *numerical_failure = true;
    return normal ? std::optional<AxisymmetricHit>{{distance, point, *normal}} : std::nullopt;
  };
  const long double transverse = static_cast<long double>(direction->x) * direction->x +
                                 static_cast<long double>(direction->y) * direction->y;
  if (transverse == 0) {
    if (direction->z == 0)
      return std::nullopt;
    const double radius = std::hypot(ray.position_m.x, ray.position_m.y);
    if (radius < mirror.inner_radius_m || radius > mirror.outer_radius_m)
      return std::nullopt;
    return make_hit((mirror.vertex_z_m +
                     mirror.surface.sag(std::hypot(ray.position_m.x, ray.position_m.y)) -
                     ray.position_m.z) /
                    direction->z);
  }
  const long double projection = static_cast<long double>(ray.position_m.x) * direction->x +
                                 static_cast<long double>(ray.position_m.y) * direction->y;
  const long double radial = static_cast<long double>(ray.position_m.x) * ray.position_m.x +
                             static_cast<long double>(ray.position_m.y) * ray.position_m.y;
  const long double radius_squared =
      static_cast<long double>(mirror.outer_radius_m) * mirror.outer_radius_m;
  const long double discriminant = projection * projection - transverse * (radial - radius_squared);
  if (discriminant < 0)
    return std::nullopt;
  const long double low = std::max(static_cast<long double>(minimum_t_m),
                                   (-projection - std::sqrt(discriminant)) / transverse);
  const long double high = (-projection + std::sqrt(discriminant)) / transverse;
  if (high < low)
    return std::nullopt;
  if (high == low)
    return make_hit(static_cast<double>(low));
  const long double span = high - low;
  const long double px = ray.position_m.x + low * direction->x;
  const long double py = ray.position_m.y + low * direction->y;
  const long double vx = span * direction->x, vy = span * direction->y;
  const long double scale_squared =
      static_cast<long double>(mirror.surface.radial_scale_m) * mirror.surface.radial_scale_m;
  const std::array<long double, 3> q{(px * px + py * py) / scale_squared,
                                     2 * (px * vx + py * vy) / scale_squared,
                                     (vx * vx + vy * vy) / scale_squared};
  detail::Polynomial polynomial{}, power{};
  power[0] = 1;
  for (std::size_t i = 0; i < mirror.surface.coefficient_m.size(); ++i) {
    for (std::size_t j = 0; j <= 2 * i; ++j)
      polynomial[j] -= mirror.surface.coefficient_m[i] * power[j];
    if (i + 1 == mirror.surface.coefficient_m.size())
      break;
    detail::Polynomial next{};
    for (std::size_t j = 0; j <= 2 * i; ++j)
      for (std::size_t k = 0; k < q.size(); ++k)
        next[j + k] += power[j] * q[k];
    power = next;
  }
  polynomial[0] += ray.position_m.z + low * direction->z - mirror.vertex_z_m;
  polynomial[1] += span * direction->z;
  for (const auto coefficient : polynomial)
    if (!std::isfinite(coefficient)) {
      if (numerical_failure)
        *numerical_failure = true;
      return std::nullopt;
    }
  const auto roots = detail::polynomial_roots_unit_interval(polynomial, 24);
  for (std::size_t i = 0; i < roots.count; ++i)
    if (auto hit = make_hit(static_cast<double>(low + span * roots.values[i])))
      return hit;
  return std::nullopt;
}

[[nodiscard]] inline std::optional<AxisymmetricHit>
intersect_axisymmetric_mirror(const Ray &ray, const AxisymmetricMirror &mirror,
                              double minimum_t_m = kEpsilon) {
  return intersect_axisymmetric_mirror(ray, mirror, minimum_t_m, [](const Vec3 &) { return true; });
}

[[nodiscard]] inline std::optional<Ray> reflect(const Ray &incident, const AxisymmetricHit &hit) {
  const auto direction = normalised_checked(incident.direction);
  const auto normal = normalised_checked(hit.unit_normal);
  if (!direction || !normal)
    return std::nullopt;
  const auto reflected =
      normalised_checked(*direction - *normal * (2.0 * dot(*direction, *normal)));
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
[[nodiscard]] inline EvenPolynomialSurface
centimetre_even_polynomial_to_metres(const std::array<double, 13> &coefficient_cm) {
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
[[nodiscard]] inline EvenPolynomialSurface
centimetre_normalised_radius_polynomial_to_metres(const std::array<double, 13> &coefficient_cm,
                                                  double radial_scale_m) {
  EvenPolynomialSurface surface{};
  if (!std::isfinite(radial_scale_m) || radial_scale_m <= kEpsilon)
    return surface;
  surface.radial_scale_m = radial_scale_m;
  for (std::size_t index = 0; index < coefficient_cm.size(); ++index) {
    surface.coefficient_m[index] = coefficient_cm[index] * 0.01;
  }
  return surface;
}

// sim_telarray's reference-radius convention is z = R sum a_i (r/R)^(2i).
// R is in metres here. It scales BOTH sag and radius, independently of the
// physical aperture. See sim_config.c's parameter / R^(2i-1) conversion.
[[nodiscard]] inline EvenPolynomialSurface
reference_radius_polynomial(const std::array<double, 13> &coefficients, double reference_radius_m) {
  EvenPolynomialSurface surface{};
  surface.radial_scale_m = reference_radius_m;
  for (std::size_t i = 0; i < coefficients.size(); ++i)
    surface.coefficient_m[i] = coefficients[i] * reference_radius_m;
  return surface;
}

} // namespace obdeect
