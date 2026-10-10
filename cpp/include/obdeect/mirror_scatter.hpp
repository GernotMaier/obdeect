#pragma once

#include "obdeect/interactions.hpp"
#include "obdeect/random.hpp"

#include <numbers>

namespace obdeect {

enum class MirrorScatterMethod { outgoing_angles, surface_slopes };
struct MirrorScatter {
  double sigma1_rad{};
  double fraction2{};
  double sigma2_rad{};
  MirrorScatterMethod method{};

  [[nodiscard]] bool is_valid() const {
    return std::isfinite(sigma1_rad) && sigma1_rad >= 0 && sigma1_rad < std::numbers::pi / 2 &&
           std::isfinite(sigma2_rad) && sigma2_rad >= 0 && sigma2_rad < std::numbers::pi / 2 &&
           std::isfinite(fraction2) && fraction2 >= 0 && fraction2 <= 1 &&
           (method == MirrorScatterMethod::outgoing_angles ||
            method == MirrorScatterMethod::surface_slopes);
  }
};

// The same identity, seed, surface and encounter always resolve the same two
// normal deviates. Box--Muller uses u in (0,1], so log never sees zero.
[[nodiscard]] inline std::optional<Vec3>
reflect_with_scatter(Vec3 incoming, Vec3 normal, Vec3 tangent, const MirrorScatter &scatter,
                     std::uint64_t photon_id, std::uint32_t surface_id,
                     std::uint64_t ray_tracing_seed, std::uint32_t encounter = 0) {
  const auto n = normalised_checked(normal);
  if (!n)
    return std::nullopt;
  const std::uint64_t dimension = 16 + 4 * static_cast<std::uint64_t>(encounter);
  const auto seed =
      ray_tracing_seed ^ (static_cast<std::uint64_t>(surface_id) * 0x9e3779b97f4a7c15ULL);
  const double sigma = source_uniform(photon_id, seed, dimension) < scatter.fraction2
                           ? scatter.sigma2_rad
                           : scatter.sigma1_rad;
  if (sigma == 0)
    return reflect_specular(incoming, *n);
  const double radius =
      sigma * std::sqrt(-2 * std::log(1 - source_uniform(photon_id, seed, dimension + 1)));
  const double phi = 2 * std::numbers::pi * source_uniform(photon_id, seed, dimension + 2);
  const double first_angle = radius * std::cos(phi), second_angle = radius * std::sin(phi);
  if (scatter.method == MirrorScatterMethod::surface_slopes) {
    // Polynomial prescription uses the local graph normal (-dz/dx,-dz/dy,1).
    // A caller supplies this frame's +z; vertical graph normals are unsupported.
    if (std::abs(n->z) <= kEpsilon)
      return std::nullopt;
    const Vec3 slope_normal{n->x / n->z + first_angle / 2, n->y / n->z + second_angle / 2, 1};
    return reflect_specular(incoming, slope_normal);
  }
  const auto reflected = reflect_specular(incoming, *n);
  auto first = normalised_checked(tangent - *n * dot(tangent, *n));
  if (!first)
    first = normalised_checked(cross(std::abs(n->z) < 0.9 ? Vec3{0, 0, 1} : Vec3{1, 0, 0}, *n));
  if (!reflected || !first)
    return std::nullopt;
  const Vec3 second = cross(*n, *first);
  const auto rotate = [](Vec3 d, Vec3 axis, double angle) {
    return d * std::cos(angle) + cross(axis, d) * std::sin(angle) +
           axis * (dot(axis, d) * (1 - std::cos(angle)));
  };
  return normalised_checked(rotate(rotate(*reflected, second, first_angle), *first, second_angle));
}

} // namespace obdeect
