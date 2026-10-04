#pragma once

#include "obdeect/segmented_optical_model.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace obdeect {

// Finite detector entrance planes with an immutable, contiguous bounding-volume
// index. Compilation may allocate; traversal is stackless and never allocates.
class CompiledDetectorPlanes {
public:
  [[nodiscard]] const std::vector<ImportedDetectorSurface> &surfaces() const { return surfaces_; }

  [[nodiscard]] std::optional<DetectorSurfaceHit> intersect(const Ray &input,
                                                            double minimum_t_m = kEpsilon) const {
    const auto direction = normalised_checked(input.direction);
    if (!direction || !std::isfinite(input.position_m.x) || !std::isfinite(input.position_m.y) ||
        !std::isfinite(input.position_m.z) || !std::isfinite(minimum_t_m) || minimum_t_m < 0)
      return std::nullopt;
    const Ray ray{input.position_m, *direction};
    std::optional<DetectorSurfaceHit> nearest;
    std::size_t index = 0;
    while (index < nodes_.size()) {
      const auto &node = nodes_[index];
      const double maximum_t =
          nearest ? nearest->distance_m : std::numeric_limits<double>::infinity();
      if (!intersects_bounds(ray, node, minimum_t_m, maximum_t)) {
        index = node.escape;
        continue;
      }
      if (node.end > node.begin) {
        for (std::size_t surface = node.begin; surface < node.end; ++surface) {
          const auto hit =
              intersect_detector_surface_unchecked(ray, surfaces_[surface], minimum_t_m);
          if (hit &&
              (!nearest || hit->distance_m < nearest->distance_m ||
               (hit->distance_m == nearest->distance_m && hit->surface_id < nearest->surface_id)))
            nearest = hit;
        }
      }
      ++index;
    }
    return nearest;
  }

  [[nodiscard]] static std::optional<CompiledDetectorPlanes>
  compile(std::vector<ImportedDetectorSurface> surfaces) {
    if (surfaces.empty())
      return std::nullopt;
    std::unordered_set<std::uint32_t> ids;
    for (auto &surface : surfaces) {
      if (!is_valid(surface) || surface.shape > FacetShape::hexagon_flat_x ||
          surface.id == std::numeric_limits<std::uint32_t>::max() || !ids.insert(surface.id).second)
        return std::nullopt;
      surface.unit_normal = *normalised_checked(surface.unit_normal);
      if (surface.shape != FacetShape::circle)
        surface.unit_tangent_u = *normalised_checked(surface.unit_tangent_u);
      const double radius = bounding_radius(surface);
      for (std::size_t axis = 0; axis < 3; ++axis)
        if (!std::isfinite(coordinate(surface.centre_m, axis) - radius) ||
            !std::isfinite(coordinate(surface.centre_m, axis) + radius))
          return std::nullopt;
    }
    return CompiledDetectorPlanes(std::move(surfaces));
  }

private:
  struct Node {
    std::array<double, 3> lower{};
    std::array<double, 3> upper{};
    std::size_t begin{};
    std::size_t end{};
    std::size_t escape{};
  };

  explicit CompiledDetectorPlanes(std::vector<ImportedDetectorSurface> surfaces)
      : surfaces_(std::move(surfaces)) {
    nodes_.reserve(surfaces_.size() * 2);
    build(0, surfaces_.size());
  }

  [[nodiscard]] static double coordinate(const Vec3 &point, std::size_t axis) {
    return axis == 0 ? point.x : axis == 1 ? point.y : point.z;
  }

  [[nodiscard]] static double bounding_radius(const ImportedDetectorSurface &surface) {
    // Circumscribed circles cover every aperture orientation, including tilted
    // squares/hexagons. Padding includes the primitive's boundary tolerance.
    const double factor = surface.shape == FacetShape::square   ? std::sqrt(0.5)
                          : surface.shape == FacetShape::circle ? 0.5
                                                                : 1 / std::sqrt(3.0);
    return surface.diameter_m * factor + 2 * kEpsilon;
  }

  [[nodiscard]] static bool intersects_bounds(const Ray &ray, const Node &node, double lower_t,
                                              double upper_t) {
    for (std::size_t axis = 0; axis < 3; ++axis) {
      const double origin = coordinate(ray.position_m, axis);
      const double direction = coordinate(ray.direction, axis);
      if (direction == 0) {
        if (origin < node.lower[axis] || origin > node.upper[axis])
          return false;
      } else {
        const double first = (node.lower[axis] - origin) / direction;
        const double second = (node.upper[axis] - origin) / direction;
        lower_t = std::max(lower_t, std::min(first, second));
        upper_t = std::min(upper_t, std::max(first, second));
        if (lower_t > upper_t)
          return false;
      }
    }
    return true;
  }

  void build(std::size_t begin, std::size_t end) {
    Node node;
    node.lower.fill(std::numeric_limits<double>::infinity());
    node.upper.fill(-std::numeric_limits<double>::infinity());
    for (std::size_t index = begin; index < end; ++index) {
      const auto &surface = surfaces_[index];
      const double radius = bounding_radius(surface);
      for (std::size_t axis = 0; axis < 3; ++axis) {
        const double centre = coordinate(surface.centre_m, axis);
        node.lower[axis] = std::min(node.lower[axis], centre - radius);
        node.upper[axis] = std::max(node.upper[axis], centre + radius);
      }
    }
    const std::size_t index = nodes_.size();
    nodes_.push_back(node);
    if (end - begin <= 4) {
      nodes_[index].begin = begin;
      nodes_[index].end = end;
    } else {
      std::size_t axis = 0;
      for (std::size_t candidate = 1; candidate < 3; ++candidate)
        if (node.upper[candidate] - node.lower[candidate] > node.upper[axis] - node.lower[axis])
          axis = candidate;
      const std::size_t middle = begin + (end - begin) / 2;
      std::nth_element(surfaces_.begin() + static_cast<std::ptrdiff_t>(begin),
                       surfaces_.begin() + static_cast<std::ptrdiff_t>(middle),
                       surfaces_.begin() + static_cast<std::ptrdiff_t>(end),
                       [axis](const auto &left, const auto &right) {
                         const double first = coordinate(left.centre_m, axis);
                         const double second = coordinate(right.centre_m, axis);
                         return first < second || (first == second && left.id < right.id);
                       });
      build(begin, middle);
      build(middle, end);
    }
    nodes_[index].escape = nodes_.size();
  }

  std::vector<ImportedDetectorSurface> surfaces_;
  std::vector<Node> nodes_;
};

[[nodiscard]] inline std::optional<CompiledDetectorPlanes>
compile_detector_planes(std::vector<ImportedDetectorSurface> surfaces) {
  return CompiledDetectorPlanes::compile(std::move(surfaces));
}

} // namespace obdeect
