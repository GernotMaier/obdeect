#pragma once

#include "obdeect/random.hpp"
#include "obdeect/sources.hpp"

#include <optional>
#include <span>

namespace obdeect {

[[nodiscard]] inline Vec3 sampled_pupil(std::uint64_t id, std::uint64_t seed, double radius_m) {
  const double radius = radius_m * std::sqrt(source_uniform(id, seed, 0));
  const double phi = 2 * std::numbers::pi * source_uniform(id, seed, 1);
  return {radius * std::cos(phi), radius * std::sin(phi), 0};
}

struct FiniteStarSource {
  // The source is a physical point above a horizontal entrance plane. Source
  // propagation up to that plane is included in the initial photon time.
  Vec3 position_m{0, 0, 10000};
  double entrance_z_m{50};
  double wavelength_nm{400};
  double emission_time_ns{};
};

[[nodiscard]] inline std::optional<OpticalPhoton>
sample_star(std::uint64_t id, std::uint64_t seed, double radius_m, const StarSource &source) {
  if (!std::isfinite(radius_m) || radius_m <= 0 || !is_valid_wavelength(source.wavelength_nm) ||
      !std::isfinite(source.source_plane_z_m))
    return std::nullopt;
  const auto direction =
      normalised_checked({std::tan(source.field_x_rad), std::tan(source.field_y_rad), -1});
  if (!direction || std::abs(direction->z) <= kEpsilon)
    return std::nullopt;
  const Vec3 centre{source.source_plane_z_m * direction->x / direction->z,
                    source.source_plane_z_m * direction->y / direction->z, source.source_plane_z_m};
  const Vec3 pupil = sampled_pupil(id, seed, radius_m);
  return OpticalPhoton{{centre + pupil, *direction},
                       id,
                       source.wavelength_nm,
                       dot(*direction, pupil) / kSpeedOfLightMPerNs,
                       1};
}

[[nodiscard]] inline std::optional<OpticalPhoton>
sample_star(std::uint64_t id, std::uint64_t seed, double radius_m, const FiniteStarSource &source) {
  if (!std::isfinite(radius_m) || radius_m <= 0 || !is_valid_wavelength(source.wavelength_nm) ||
      !std::isfinite(source.position_m.x) || !std::isfinite(source.position_m.y) ||
      !std::isfinite(source.position_m.z) || !std::isfinite(source.entrance_z_m) ||
      source.entrance_z_m < 0 || source.position_m.z <= source.entrance_z_m ||
      !std::isfinite(source.emission_time_ns))
    return std::nullopt;
  const Vec3 target = sampled_pupil(id, seed, radius_m);
  const Vec3 separation = target - source.position_m;
  const auto direction = normalised_checked(separation);
  if (!direction)
    return std::nullopt;
  const double travelled = (source.entrance_z_m - source.position_m.z) / direction->z;
  const Vec3 entrance = target + *direction * (source.entrance_z_m / direction->z);
  return OpticalPhoton{{entrance, *direction},
                       id,
                       source.wavelength_nm,
                       source.emission_time_ns + travelled / kSpeedOfLightMPerNs,
                       1};
}

// Normalization is explicit so independent blocks can represent one frozen
// emitted level. Changing the number of calls never silently changes weight.
[[nodiscard]] inline std::optional<OpticalPhoton>
sample_illuminator(std::uint64_t id, std::uint64_t seed, double radius_m,
                   const PointIlluminator &source, std::uint64_t normalization_samples) {
  if (!normalization_samples || !std::isfinite(radius_m) || radius_m <= 0 ||
      !is_valid_wavelength(source.wavelength_nm) || !std::isfinite(source.emitted_weight) ||
      source.emitted_weight < 0 || !std::isfinite(source.position_m.x) ||
      !std::isfinite(source.position_m.y) || !std::isfinite(source.position_m.z))
    return std::nullopt;
  const Vec3 separation = sampled_pupil(id, seed, radius_m) - source.position_m;
  const auto direction = normalised_checked(separation);
  if (!direction)
    return std::nullopt;
  const double distance_squared = dot(separation, separation);
  const double weight = source.emitted_weight * radius_m * radius_m * std::abs(direction->z) /
                        (4.0 * static_cast<double>(normalization_samples) * distance_squared);
  return OpticalPhoton{{source.position_m, *direction}, id, source.wavelength_nm, 0, weight};
}

[[nodiscard]] inline std::optional<OpticalPhoton> sample_laser(std::uint64_t id, std::uint64_t seed,
                                                               double beam_radius_m,
                                                               const LaserSource &source) {
  const auto central = normalised_checked(source.direction);
  if (!central || !std::isfinite(beam_radius_m) || beam_radius_m <= 0 ||
      !is_valid_wavelength(source.wavelength_nm) ||
      !std::isfinite(source.divergence_half_angle_rad) || source.divergence_half_angle_rad < 0 ||
      source.divergence_half_angle_rad >= std::numbers::pi / 2 ||
      !std::isfinite(source.origin_m.x) || !std::isfinite(source.origin_m.y) ||
      !std::isfinite(source.origin_m.z))
    return std::nullopt;
  const Vec3 reference = std::abs(central->z) < 0.9 ? Vec3{0, 0, 1} : Vec3{1, 0, 0};
  const auto first = normalised_checked(cross(*central, reference));
  if (!first)
    return std::nullopt;
  const Vec3 second = cross(*central, *first);
  const Vec3 pupil = sampled_pupil(id, seed, beam_radius_m);
  const Vec3 entrance = source.origin_m + *first * pupil.x + second * pupil.y;
  // Uniform solid angle within the declared cone, rather than the small-angle
  // theta*sqrt(u) approximation. Spatial and angular dimensions are independent.
  const double cosine =
      1 - source_uniform(id, seed, 2) * (1 - std::cos(source.divergence_half_angle_rad));
  const double sine = std::sqrt(std::max(0.0, 1 - cosine * cosine));
  const double phi = 2 * std::numbers::pi * source_uniform(id, seed, 3);
  const auto direction = normalised_checked(*central * cosine + *first * (sine * std::cos(phi)) +
                                            second * (sine * std::sin(phi)));
  if (!direction)
    return std::nullopt;
  return OpticalPhoton{{entrance, *direction}, id, source.wavelength_nm, 0, 1};
}

// Source adapter writes directly to reusable caller-owned storage. This is a
// batch boundary; transport sees only resolved OpticalPhoton data.
template <class Sampler>
[[nodiscard]] inline bool fill_source(std::span<OpticalPhoton> destination, std::uint64_t first_id,
                                      Sampler sampler) {
  if (!destination.empty() &&
      destination.size() - 1 > std::numeric_limits<std::uint64_t>::max() - first_id)
    return false;
  for (std::size_t index = 0; index < destination.size(); ++index) {
    const auto photon = sampler(first_id + index);
    if (!photon)
      return false;
    destination[index] = *photon;
  }
  return true;
}

} // namespace obdeect
