#pragma once

#include "obdeect/obscurers.hpp"
#include "obdeect/surfaces.hpp"
#include "obdeect/table.hpp"

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

namespace obdeect {

// Thrown for any invalid or unsupported optical model; for a model file the
// message names the JSON path.
struct ModelError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Random reflection error: with probability fraction2 the Gaussian width is
// sigma2, otherwise sigma1 (radians, per axis). `surface_slopes` perturbs the
// surface slope by half the angle; otherwise the reflected direction is rotated.
struct Scatter {
  double sigma1_rad{}, fraction2{}, sigma2_rad{};
  bool surface_slopes{};
};

struct Mirror {
  std::vector<Facet> facets; // Segmented mirror, or ...
  FacetGrid grid;
  std::optional<Asphere> asphere; // ... a continuous asphere with optional segments.
  std::vector<Segment> segments;
  std::vector<Segment::TangentFrame> segment_frames; // One per segment, from the asphere.
  SegmentGrid segment_grid;                          // Segments by x-y cell, sorted by id.
  std::optional<Response> reflectivity;
  std::optional<Scatter> scatter;
  // Relative reflectivity map, in facet map coordinates (Facet::map_u,
  // map_v) for facets and in telescope x-y for an asphere.
  std::optional<Map2D> degradation;
};

// Direction-dependent transmission T = on_axis / (1 + amplitude (sin(theta) /
// scale)^power)^outer_power, theta measured from the optical axis.
struct TelescopeTransmission {
  double on_axis{1.0}, amplitude{}, angular_scale_rad{1.0}, power{2.0}, outer_power{1.0};

  [[nodiscard]] double at(const Vec3 &direction) const {
    // Constant transmission: exact, as pow(1, p) == 1 and on_axis / 1.0 == on_axis.
    if (amplitude == 0.0)
      return on_axis;
    const double s = std::hypot(direction.x, direction.y);
    return on_axis /
           std::pow(1.0 + amplitude * std::pow(s / angular_scale_rad, power), outer_power);
  }
};

// Photons end on the focal surface: a plane z = plane_z or an asphere. The
// plane is unbounded, so the image is never truncated.
struct FocalSurface {
  std::uint32_t id{};
  double plane_z{};
  std::optional<Asphere> asphere;
};

struct OpticalModel {
  Mirror primary;
  std::optional<Mirror> secondary;
  // Without a dedicated shadow, light from the sky that meets the back of the
  // secondary mirror is absorbed there.
  bool secondary_blocks_incoming{};
  FocalSurface focal_surface;
  Obscurers obscurers;
  std::optional<TelescopeTransmission> telescope_transmission;
  double propagation_group_index{1.0};
};

} // namespace obdeect
