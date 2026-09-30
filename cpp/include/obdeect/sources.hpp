#pragma once

#include "obdeect/math.hpp"

#include <cstdint>
#include <limits>
#include <numbers>
#include <string_view>
#include <vector>

namespace obdeect {

struct OpticalPhoton {
  Ray ray;
  std::uint64_t photon_id{};
  double wavelength_nm{400.0};
  double time_ns{};
  double weight{1.0};
  // Input adapters retain raw-bunch provenance here.  A wavelength of zero
  // means that the input did not specify a spectrum; it must be resolved by a
  // spectrum adapter before wavelength-dependent transport.
  std::uint64_t bunch_id{};
  double emission_height_m{std::numeric_limits<double>::quiet_NaN()};
  double emission_distance_m{std::numeric_limits<double>::quiet_NaN()};
};

enum class ArtificialSourceKind : std::uint8_t { star, illuminator, laser };

inline std::string_view to_string(ArtificialSourceKind kind) {
  switch (kind) {
    case ArtificialSourceKind::star: return "star";
    case ArtificialSourceKind::illuminator: return "illuminator";
    case ArtificialSourceKind::laser: return "laser";
  }
  return "unknown";
}

struct StarSource {
  // Telescope-frame apparent field direction. Positive x/y move the source
  // on the sky; photons consequently propagate in negative z.
  double field_x_rad{};
  double field_y_rad{};
  double source_plane_z_m{50.0};
  double wavelength_nm{400.0};
};

struct PointIlluminator {
  // A finite-distance flasher/illuminator. Weight follows 1/r^2 so rays are
  // physically comparable across the pupil before any telescope throughput.
  Vec3 position_m{0.0, 0.0, 50.0};
  double wavelength_nm{400.0};
  double emitted_weight{1.0};
};

struct LaserSource {
  // Beam direction is the central propagation direction, with divergence as
  // a hard half-angle sampling limit. A zero-divergence laser is collimated.
  Vec3 direction{0.0, 0.0, -1.0};
  // Centre of the launch plane. The sampled plane is normal to ``direction``.
  Vec3 origin_m{0.0, 0.0, 50.0};
  double divergence_half_angle_rad{};
  double wavelength_nm{400.0};
};

inline bool is_valid_wavelength(double wavelength_nm) {
  return std::isfinite(wavelength_nm) && wavelength_nm > 0.0;
}

inline std::vector<Vec3> fibonacci_pupil_points(std::size_t count, double radius_m, double z_m) {
  std::vector<Vec3> points;
  if (count == 0 || !std::isfinite(radius_m) || radius_m <= kEpsilon || !std::isfinite(z_m)) return points;
  points.reserve(count);
  constexpr double golden_ratio_conjugate = 0.6180339887498948482;
  for (std::size_t index = 0; index < count; ++index) {
    const double u = (static_cast<double>(index) + 0.5) / static_cast<double>(count);
    const double radius = radius_m * std::sqrt(u);
    const double phi = std::fmod(static_cast<double>(index) * golden_ratio_conjugate, 1.0) *
                       2.0 * std::numbers::pi;
    points.push_back({radius * std::cos(phi), radius * std::sin(phi), z_m});
  }
  return points;
}

inline std::vector<OpticalPhoton> star_photons(std::size_t count, double pupil_radius_m,
                                                const StarSource& source) {
  std::vector<OpticalPhoton> photons;
  if (!is_valid_wavelength(source.wavelength_nm)) return photons;
  const auto direction = normalised_checked(
      {std::tan(source.field_x_rad), std::tan(source.field_y_rad), -1.0});
  if (!direction) return photons;
  if (!std::isfinite(source.source_plane_z_m) || std::abs(direction->z) <= kEpsilon) return photons;
  // Aim the centre of the parallel bundle at the telescope optical axis.  A
  // horizontal launch plane at z alone would displace an off-axis star by
  // z*tan(field angle) at the primary; at astronomical source distances that
  // misses the telescope completely.
  const Vec3 launch_centre{source.source_plane_z_m * direction->x / direction->z,
                           source.source_plane_z_m * direction->y / direction->z,
                           source.source_plane_z_m};
  auto positions = fibonacci_pupil_points(count, pupil_radius_m, source.source_plane_z_m);
  for (auto& position : positions) {
    position.x += launch_centre.x;
    position.y += launch_centre.y;
  }
  photons.reserve(positions.size());
  // Samples on a horizontal launch plane have different phases for an
  // off-axis plane wave. Points farther along the propagation direction
  // cross a common wavefront later: t(r) = d dot (r - reference) / c.
  const Vec3 reference = launch_centre;
  for (std::size_t index = 0; index < positions.size(); ++index) {
    photons.push_back({{positions[index], *direction}, static_cast<std::uint64_t>(index),
                       source.wavelength_nm, dot(*direction, positions[index] - reference) /
                                                 kSpeedOfLightMPerNs,
                       1.0});
  }
  return photons;
}

inline std::vector<OpticalPhoton> illuminator_photons(std::size_t count, double pupil_radius_m,
                                                       const PointIlluminator& source) {
  std::vector<OpticalPhoton> photons;
  if (!is_valid_wavelength(source.wavelength_nm) || !std::isfinite(source.emitted_weight) ||
      source.emitted_weight < 0.0) {
    return photons;
  }
  const auto targets = fibonacci_pupil_points(count, pupil_radius_m, 0.0);
  const double sampled_area_m2 = std::numbers::pi * pupil_radius_m * pupil_radius_m;
  photons.reserve(targets.size());
  for (std::size_t index = 0; index < targets.size(); ++index) {
    const Vec3 separation = targets[index] - source.position_m;
    const double distance_m = norm(separation);
    const auto direction = normalised_checked(separation);
    if (!direction || distance_m <= kEpsilon) return {};
    // ``emitted_weight`` is the isotropically emitted photon count/weight.
    // Uniform-area pupil sampling represents the solid angle subtended by
    // each sample, including the projected pupil area.
    const double projected_area = sampled_area_m2 * std::abs(direction->z);
    const double sample_weight = source.emitted_weight * projected_area /
                                 (4.0 * std::numbers::pi * static_cast<double>(targets.size()) *
                                  distance_m * distance_m);
    photons.push_back({{source.position_m, *direction}, static_cast<std::uint64_t>(index), source.wavelength_nm,
                       0.0, sample_weight});
  }
  return photons;
}

inline std::vector<OpticalPhoton> laser_photons(std::size_t count, double pupil_radius_m,
                                                 const LaserSource& source) {
  std::vector<OpticalPhoton> photons;
  const auto central_direction = normalised_checked(source.direction);
  if (!central_direction || !is_valid_wavelength(source.wavelength_nm) ||
      !std::isfinite(source.divergence_half_angle_rad) || source.divergence_half_angle_rad < 0.0 ||
      source.divergence_half_angle_rad >= std::numbers::pi / 2.0) {
    return photons;
  }
  if (!std::isfinite(pupil_radius_m) || pupil_radius_m <= kEpsilon ||
      !std::isfinite(source.origin_m.x) || !std::isfinite(source.origin_m.y) ||
      !std::isfinite(source.origin_m.z)) {
    return photons;
  }
  photons.reserve(count);
  // Construct a stable transverse basis for a cone around the laser axis.
  const Vec3 seed = std::abs(central_direction->z) < 0.9 ? Vec3{0.0, 0.0, 1.0} : Vec3{1.0, 0.0, 0.0};
  const auto basis_x = normalised_checked(
      {central_direction->y * seed.z - central_direction->z * seed.y,
       central_direction->z * seed.x - central_direction->x * seed.z,
       central_direction->x * seed.y - central_direction->y * seed.x});
  if (!basis_x) return photons;
  const Vec3 basis_y{central_direction->y * basis_x->z - central_direction->z * basis_x->y,
                     central_direction->z * basis_x->x - central_direction->x * basis_x->z,
                     central_direction->x * basis_x->y - central_direction->y * basis_x->x};
  constexpr double golden_ratio_conjugate = 0.6180339887498948482;
  for (std::size_t index = 0; index < count; ++index) {
    // Position and direction use deliberately different low-discrepancy
    // dimensions, avoiding a radial/angle correlation in finite beams.
    const double position_u = (static_cast<double>(index) + 0.5) / static_cast<double>(count);
    const double position_radius = pupil_radius_m * std::sqrt(position_u);
    const double position_phi = std::fmod(static_cast<double>(index) * golden_ratio_conjugate, 1.0) *
                                2.0 * std::numbers::pi;
    const Vec3 position = source.origin_m + *basis_x * (position_radius * std::cos(position_phi)) +
                          basis_y * (position_radius * std::sin(position_phi));
    const std::size_t direction_index = (index * (count - 1U) + 17U) % count;
    const double u = (static_cast<double>(direction_index) + 0.5) / static_cast<double>(count);
    const double angle = source.divergence_half_angle_rad * std::sqrt(u);
    const double phi = std::fmod(static_cast<double>(direction_index) * golden_ratio_conjugate, 1.0) *
                       2.0 * std::numbers::pi;
    const auto direction = normalised_checked(*central_direction * std::cos(angle) +
                                              *basis_x * (std::sin(angle) * std::cos(phi)) +
                                              basis_y * (std::sin(angle) * std::sin(phi)));
    if (!direction) return {};
    photons.push_back({{position, *direction}, static_cast<std::uint64_t>(index), source.wavelength_nm,
                       0.0, 1.0});
  }
  return photons;
}

// Apply a deterministic top-hat emission profile without changing the source
// geometry or the plane-wave phase already carried by star photons.  The
// profile is useful for optical arrival-time studies; source-specific spectra
// and pulse shapes remain an adapter responsibility.
inline bool apply_top_hat_emission_times(std::vector<OpticalPhoton>& photons, double start_ns,
                                         double width_ns) {
  if (!std::isfinite(start_ns) || !std::isfinite(width_ns) || width_ns < 0.0) return false;
  if (photons.empty()) return true;
  const double count = static_cast<double>(photons.size());
  for (std::size_t index = 0; index < photons.size(); ++index) {
    photons[index].time_ns += start_ns + width_ns * (static_cast<double>(index) + 0.5) / count;
  }
  return true;
}

}  // namespace obdeect
