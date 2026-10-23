#pragma once

#include "obdeect/random.hpp"
#include "obdeect/vec.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

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

// Unnormalised probability density, linear between nodes and zero outside,
// sampled by inverting its cumulative distribution. An empty density is a
// point mass at zero.
class LinearDensity {
public:
  LinearDensity() = default;
  LinearDensity(std::vector<double> x, std::vector<double> y) : x_(std::move(x)), y_(std::move(y)) {
    if (x_.size() < 2 || x_.size() != y_.size())
      throw std::invalid_argument("a density needs at least two nodes and one value per node");
    cdf_.assign(x_.size(), 0.0);
    for (std::size_t i = 0; i < x_.size(); ++i) {
      if (!std::isfinite(x_[i]) || !std::isfinite(y_[i]) || y_[i] < 0.0 ||
          (i > 0 && !(x_[i] > x_[i - 1])))
        throw std::invalid_argument("density nodes must increase and values be finite and >= 0");
      if (i > 0)
        cdf_[i] = cdf_[i - 1] + 0.5 * (y_[i] + y_[i - 1]) * (x_[i] - x_[i - 1]);
    }
    if (!(cdf_.back() > 0.0))
      throw std::invalid_argument("a density must have a positive integral");
  }

  [[nodiscard]] double integral() const { return cdf_.back(); }
  [[nodiscard]] double maximum() const { return *std::max_element(y_.begin(), y_.end()); }

  [[nodiscard]] double at(double x) const {
    if (x < x_.front() || x > x_.back())
      return 0.0;
    const auto i = segment(x_, x);
    return y_[i] + (y_[i + 1] - y_[i]) * (x - x_[i]) / (x_[i + 1] - x_[i]);
  }

  // The value below which a fraction u in (0, 1) of the density lies.
  [[nodiscard]] double sample(double u) const {
    if (x_.empty())
      return 0.0;
    const double area = u * cdf_.back();
    const auto i = segment(cdf_, area);
    // Solve y_i s + k s^2 / 2 = a for s in the segment, stable for any slope k.
    const double a = area - cdf_[i], k = (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
    const double root = std::sqrt(std::max(0.0, y_[i] * y_[i] + 2.0 * k * a));
    const double s = y_[i] + root > 0.0 ? 2.0 * a / (y_[i] + root) : 0.0;
    return std::min(x_[i] + s, x_[i + 1]);
  }

private:
  // Segment [i, i + 1] of the nodes containing `value` of an increasing
  // column (x_ or cdf_), clamped to the table.
  [[nodiscard]] static std::size_t segment(const std::vector<double> &column, double value) {
    const auto upper = std::upper_bound(column.begin(), column.end(), value) - column.begin();
    return static_cast<std::size_t>(
        std::clamp<std::ptrdiff_t>(upper - 1, 0, static_cast<std::ptrdiff_t>(column.size()) - 2));
  }

  std::vector<double> x_, y_, cdf_;
};

// A point emitter such as a calibration light source. Its intensity per solid
// angle depends on u = 1 - cos(theta), theta measured from `axis`; emission
// times follow `pulse` [ns]. Only photons crossing a target disk (centre, unit
// normal, radius) are made: target points are uniform on the disk, and the
// weight of a photon is its fluence there relative to the largest possible,
//   w = I(u) / I_max * (h / r)^3,
// with h the distance of the emitter from the disk plane and r from the
// target point. Weights are in [0, 1]; keeping each photon with probability
// w gives the photons crossing the disk, so with a Poisson number of
// candidates (mean: `candidates_per_photon` times the photons emitted) the
// kept photons are Poisson distributed as well.
struct Emitter {
  Vec3 position;
  Vec3 axis{0.0, 0.0, 1.0};
  LinearDensity intensity; // Over u = 1 - cos(theta) in [0, 2].
  LinearDensity pulse;     // Over emission time [ns]; empty: all at 0.
  Vec3 disk_centre;
  Vec3 disk_normal{0.0, 0.0, 1.0};
  double disk_radius_m{};
  double wavelength_nm{400.0};

  [[nodiscard]] double disk_distance_m() const {
    return std::abs(dot(disk_centre - position, disk_normal));
  }

  // Mean candidates per emitted photon: the fluence bound times the disk area.
  [[nodiscard]] double candidates_per_photon() const {
    const double per_steradian =
        intensity.maximum() / (6.283185307179586476925 * intensity.integral());
    const double h = disk_distance_m();
    return per_steradian * 3.141592653589793238463 * disk_radius_m * disk_radius_m / (h * h);
  }
};

struct EmittedPhoton {
  Ray ray;
  double time_ns{}, weight{};
};

[[nodiscard]] inline EmittedPhoton sample_emitter(const Emitter &emitter, std::uint64_t photon_id,
                                                  std::uint64_t seed) {
  const auto u = uniform4(seed, photon_id, Stream::source);
  // An orthonormal frame on the disk.
  const Vec3 &n = emitter.disk_normal;
  const Vec3 e1 =
      normalized(std::abs(n.x) < 0.9 ? cross(n, {1.0, 0.0, 0.0}) : cross(n, {0.0, 1.0, 0.0}));
  const Vec3 e2 = cross(n, e1);
  const double r = emitter.disk_radius_m * std::sqrt(u[0]);
  const double phi = 6.283185307179586476925 * u[1];
  const Vec3 target = emitter.disk_centre + e1 * (r * std::cos(phi)) + e2 * (r * std::sin(phi));
  const Vec3 offset = target - emitter.position;
  const double distance = norm(offset);
  const Vec3 direction = offset * (1.0 / distance);
  const double cosine = std::clamp(dot(direction, emitter.axis), -1.0, 1.0);
  const double ratio = emitter.disk_distance_m() / distance;
  const double weight =
      emitter.intensity.at(1.0 - cosine) / emitter.intensity.maximum() * ratio * ratio * ratio;
  return {{emitter.position, direction}, emitter.pulse.sample(u[2]), weight};
}

} // namespace obdeect
