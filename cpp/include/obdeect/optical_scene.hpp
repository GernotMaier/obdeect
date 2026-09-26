#pragma once

#include "obdeect/frames.hpp"
#include "obdeect/segmented_scene.hpp"

#include <cstdint>
#include <limits>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace obdeect {

enum class SurfaceRole : std::uint8_t { mirror, detector, obscurer };
constexpr std::uint32_t kMaximumSceneInteractions = 64;

// The frame fixes the plane centre, normal, and in-plane orientation. The
// diameter defines a finite circle, square, or hexagon in that plane.
struct OpticalSurfaceRecord {
  std::uint32_t id{};
  RigidFrame frame{};
  FacetShape shape{FacetShape::circle};
  double diameter_m{};
  SurfaceRole role{SurfaceRole::obscurer};
  std::uint32_t material_id{std::numeric_limits<std::uint32_t>::max()};
};

struct ImportedOpticalScene {
  ModelProvenance provenance;
  std::vector<OpticalSurfaceRecord> surfaces;
  std::uint32_t max_interactions{8};
};

class CompiledOpticalScene {
 public:
  CompiledOpticalScene(CompiledOpticalScene&&) = default;
  CompiledOpticalScene(const CompiledOpticalScene&) = default;

  [[nodiscard]] const ModelProvenance& provenance() const { return provenance_; }
  [[nodiscard]] const std::vector<OpticalSurfaceRecord>& surfaces() const { return surfaces_; }
  [[nodiscard]] std::uint32_t max_interactions() const { return max_interactions_; }

 private:
  friend std::optional<CompiledOpticalScene> compile_optical_scene(const ImportedOpticalScene&);
  CompiledOpticalScene(ModelProvenance provenance, std::vector<OpticalSurfaceRecord> surfaces,
                       std::uint32_t max_interactions)
      : provenance_(std::move(provenance)), surfaces_(std::move(surfaces)),
        max_interactions_(max_interactions) {}

  ModelProvenance provenance_;
  std::vector<OpticalSurfaceRecord> surfaces_;
  std::uint32_t max_interactions_{};
};

[[nodiscard]] inline std::optional<CompiledOpticalScene> compile_optical_scene(
    const ImportedOpticalScene& input) {
  if (!has_valid_provenance(input.provenance) || input.surfaces.empty() || input.max_interactions == 0 ||
      input.max_interactions > kMaximumSceneInteractions)
    return std::nullopt;
  std::unordered_set<std::uint32_t> ids;
  ids.reserve(input.surfaces.size());
  for (const auto& surface : input.surfaces) {
    if (!surface.frame.is_valid() || !std::isfinite(surface.diameter_m) || surface.diameter_m <= kEpsilon ||
        surface.shape > FacetShape::hexagon_flat_x || surface.role > SurfaceRole::obscurer ||
        !ids.insert(surface.id).second)
      return std::nullopt;
  }
  return CompiledOpticalScene{input.provenance, input.surfaces, input.max_interactions};
}

struct OpticalSurfaceHit {
  std::uint32_t surface_id{};
  SurfaceRole role{SurfaceRole::obscurer};
  double distance_m{};
  Vec3 point_m{};
  Vec3 normal{};
};

[[nodiscard]] inline std::optional<OpticalSurfaceHit> intersect_nearest_surface(
    const Ray& ray, const CompiledOpticalScene& scene) {
  std::optional<OpticalSurfaceHit> nearest;
  for (const auto& surface : scene.surfaces()) {
    const ImportedDetectorSurface plane{surface.id, surface.frame.origin_m, surface.frame.z_axis,
                                        surface.diameter_m, surface.shape, surface.frame.x_axis};
    const auto hit = intersect_detector_surface_unchecked(ray, plane);
    if (hit && (!nearest || hit->distance_m < nearest->distance_m ||
                (hit->distance_m == nearest->distance_m && surface.id < nearest->surface_id)))
      nearest = OpticalSurfaceHit{surface.id, surface.role, hit->distance_m, hit->point_m,
                                  hit->unit_normal};
  }
  return nearest;
}

}  // namespace obdeect
