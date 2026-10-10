#pragma once

#include "obdeect/detector_planes.hpp"

#include <array>
#include <span>

namespace obdeect {

struct DetectorAssignmentGrid {
  std::uint32_t nx{}, ny{};
  double x_low_m{}, x_high_m{}, y_low_m{}, y_high_m{};
  Vec3 x_basis{}, y_basis{};
  std::optional<double> reference_plane_z_m{};
};

// Explicit projected candidate selection. Cell lists retain input surface order;
// the optical model chooses this policy instead of nearest physical intersection.
class CompiledDetectorAssignmentGrid {
public:
  [[nodiscard]] static std::optional<CompiledDetectorAssignmentGrid>
  compile(DetectorAssignmentGrid description, std::span<const ImportedDetectorSurface> surfaces) {
    const auto &g = description;
    const std::size_t cells = std::size_t(g.nx) * g.ny;
    if (!g.nx || !g.ny || cells / g.nx != g.ny || cells > 1000000 || surfaces.empty() ||
        surfaces.size() > std::numeric_limits<std::uint32_t>::max() || !std::isfinite(g.x_low_m) ||
        !std::isfinite(g.x_high_m) || !std::isfinite(g.y_low_m) || !std::isfinite(g.y_high_m) ||
        g.x_low_m >= g.x_high_m || g.y_low_m >= g.y_high_m || !normalised_checked(g.x_basis) ||
        !normalised_checked(g.y_basis) || std::abs(norm(g.x_basis) - 1) > kEpsilon ||
        std::abs(norm(g.y_basis) - 1) > kEpsilon ||
        std::abs(dot(g.x_basis, g.y_basis)) > kEpsilon ||
        (g.reference_plane_z_m && !std::isfinite(*g.reference_plane_z_m)))
      return {};
    CompiledDetectorAssignmentGrid result;
    result.surfaces_.assign(surfaces.begin(), surfaces.end());
    result.description_ = description;
    result.offsets_.resize(cells + 1);
    const auto visit_cells = [&](const auto &surface, const auto &visit) {
      const double x = dot(surface.centre_m, g.x_basis), y = dot(surface.centre_m, g.y_basis);
      const double r = surface.assignment_radius_m;
      const auto coordinate = [](double p, double low, double high, std::uint32_t count) {
        return int(std::clamp((p - low) / (high - low), 0.0, 1.0) * count);
      };
      const int x0 = std::max(0, coordinate(x - r, g.x_low_m, g.x_high_m, g.nx));
      const int x1 = std::min(int(g.nx) - 1, coordinate(x + r, g.x_low_m, g.x_high_m, g.nx));
      const int y0 = std::max(0, coordinate(y - r, g.y_low_m, g.y_high_m, g.ny));
      const int y1 = std::min(int(g.ny) - 1, coordinate(y + r, g.y_low_m, g.y_high_m, g.ny));
      for (int iy = y0; iy <= y1; ++iy)
        for (int ix = x0; ix <= x1; ++ix)
          visit(std::size_t(iy) * g.nx + std::size_t(ix));
    };
    for (const auto &surface : surfaces) {
      if (!is_valid(surface) || !std::isfinite(surface.assignment_radius_m) ||
          surface.assignment_radius_m <= 0)
        return {};
      const double x = dot(surface.centre_m, g.x_basis), y = dot(surface.centre_m, g.y_basis);
      if (x - surface.assignment_radius_m < g.x_low_m - kEpsilon ||
          x + surface.assignment_radius_m > g.x_high_m + kEpsilon ||
          y - surface.assignment_radius_m < g.y_low_m - kEpsilon ||
          y + surface.assignment_radius_m > g.y_high_m + kEpsilon)
        return {};
      visit_cells(surface, [&](std::size_t cell) { ++result.offsets_[cell + 1]; });
    }
    for (std::size_t cell = 0; cell < cells; ++cell)
      result.offsets_[cell + 1] += result.offsets_[cell];
    result.indices_.resize(result.offsets_.back());
    auto cursor = result.offsets_;
    for (std::size_t surface = 0; surface < surfaces.size(); ++surface)
      visit_cells(surfaces[surface], [&](std::size_t cell) {
        result.indices_[cursor[cell]++] = static_cast<std::uint32_t>(surface);
      });
    return result;
  }

  [[nodiscard]] const DetectorAssignmentGrid &description() const { return description_; }

  [[nodiscard]] std::optional<DetectorSurfaceHit> intersect(const Ray &ray,
                                                            const Vec3 &reference_point) const {
    const auto &g = description_;
    const double x = dot(reference_point, g.x_basis), y = dot(reference_point, g.y_basis);
    if (!std::isfinite(x) || !std::isfinite(y) || x < g.x_low_m || x > g.x_high_m ||
        y < g.y_low_m || y > g.y_high_m)
      return {};
    const auto ix = std::min(std::size_t(g.nx - 1),
                             std::size_t((x - g.x_low_m) / (g.x_high_m - g.x_low_m) * g.nx));
    const auto iy = std::min(std::size_t(g.ny - 1),
                             std::size_t((y - g.y_low_m) / (g.y_high_m - g.y_low_m) * g.ny));
    const std::size_t cell = iy * g.nx + ix;
    for (std::size_t entry = offsets_[cell]; entry < offsets_[cell + 1]; ++entry)
      if (const auto hit = intersect_detector_surface_unchecked(ray, surfaces_[indices_[entry]]))
        return hit;
    return {};
  }

private:
  DetectorAssignmentGrid description_{};
  std::vector<std::size_t> offsets_{};
  std::vector<std::uint32_t> indices_{};
  std::vector<ImportedDetectorSurface> surfaces_{};
};

} // namespace obdeect
