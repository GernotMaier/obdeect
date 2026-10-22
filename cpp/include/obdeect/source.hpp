#pragma once

#include "obdeect/random.hpp"
#include "obdeect/vec.hpp"

#include <cmath>
#include <cstdint>
#include <limits>

namespace obdeect {

// Light from a star (plane wave) or a point source at a finite distance,
// sampled as in sim_telarray's star-light set-up: target points uniform on a
// disk perpendicular to the optical axis, centred on the axis at height
// `disk_z_m`. As in sim_telarray, `distance_m` is measured from the disk
// centre. Each ray then starts where it crosses the plane `start_z_m`, above
// all telescope structure. The disk area is the sampled area for
// effective-area estimates.
struct PointSource {
  Vec3 towards_source{0.0, 0.0, 1.0}; // Unit vector from the telescope to the source.
  double distance_m{std::numeric_limits<double>::infinity()};
  double disk_radius_m{};
  double disk_z_m{};
  double start_z_m{};
  double wavelength_nm{400.0};
};

struct SourcePhoton {
  Ray ray;
  double time_ns{};
};

// Emission times are relative to the wavefront that passes the disk centre at t = 0.
[[nodiscard]] inline SourcePhoton sample_source(const PointSource &source, std::uint64_t photon_id,
                                                std::uint64_t seed) {
  const auto u = uniform4(seed, photon_id, Stream::source);
  const double r = source.disk_radius_m * std::sqrt(u[0]);
  const double phi = 6.283185307179586476925 * u[1];
  const Vec3 centre{0.0, 0.0, source.disk_z_m};
  const Vec3 target{r * std::cos(phi), r * std::sin(phi), source.disk_z_m};
  if (std::isinf(source.distance_m)) {
    const Vec3 direction = -source.towards_source;
    const Vec3 start = target + direction * ((source.start_z_m - target.z) / direction.z);
    return {{start, direction}, dot(start - centre, direction) / kSpeedOfLightMPerNs};
  }
  // A source below the start plane emits from its own position.
  const Vec3 position = centre + source.towards_source * source.distance_m;
  const Vec3 direction = normalized(target - position);
  const Vec3 start = position.z <= source.start_z_m
                         ? position
                         : target + direction * ((source.start_z_m - target.z) / direction.z);
  const double time = (norm(start - position) - source.distance_m) / kSpeedOfLightMPerNs;
  return {{start, direction}, time};
}

} // namespace obdeect
