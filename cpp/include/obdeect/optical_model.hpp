#pragma once

#include "obdeect/axisymmetric_optics.hpp"
#include "obdeect/frames.hpp"
#include "obdeect/segmented_optical_model.hpp"
#include "obdeect/tables.hpp"

#include <cstdint>
#include <limits>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace obdeect {

enum class SurfaceRole : std::uint8_t { mirror, detector, obscurer, refractive_interface };
constexpr std::uint32_t kVacuumMediumId = std::numeric_limits<std::uint32_t>::max();
enum class InterfaceTransmissionSemantics { additional_coating, complete_interface };

struct MaterialCurve {
  std::vector<double> wavelength_nm;
  std::vector<double> value;
  [[nodiscard]] Table1DView view() const { return {wavelength_nm, value}; }
};

struct OpticalMaterialRecord {
  std::uint32_t id{};
  MaterialCurve phase_index;
  MaterialCurve group_index;
  MaterialCurve absorption_per_m;

  [[nodiscard]] bool is_valid() const {
    if (id == kVacuumMediumId || !phase_index.view().is_valid() || !group_index.view().is_valid() ||
        !absorption_per_m.view().is_valid())
      return false;
    for (const auto *curve : {&phase_index, &group_index, &absorption_per_m})
      for (const auto wavelength : curve->wavelength_nm)
        if (wavelength <= 0)
          return false;
    for (const auto *curve : {&phase_index, &group_index})
      for (const auto value : curve->value)
        if (value <= 0)
          return false;
    return true;
  }
};

struct ResolvedOpticalMaterial {
  double phase_index{};
  double group_index{};
  double absorption_per_m{};
};
constexpr std::uint32_t kMaximumOpticalModelInteractions = 64;

// The frame fixes the plane centre, normal, and in-plane orientation. The
// diameter defines a finite circle, square, or hexagon in that plane.
struct OpticalSurfaceRecord {
  std::uint32_t id{};
  RigidFrame frame{};
  FacetShape shape{FacetShape::circle};
  double diameter_m{};
  SurfaceRole role{SurfaceRole::obscurer};
  std::uint32_t material_id{std::numeric_limits<std::uint32_t>::max()};
  // Optional local z(r) surface with a circular/annular aperture. The frame
  // maps its vertex and local normal into the parent optical_model.
  std::optional<EvenPolynomialSurface> sag{};
  double inner_radius_m{};
  // Frame +z points into front_medium_id; -z points into back_medium_id.
  // The vacuum sentinel is explicit; every other ID must reference a material.
  std::uint32_t front_medium_id{kVacuumMediumId};
  std::uint32_t back_medium_id{kVacuumMediumId};
  std::optional<SpectralResponse> transmission{};
  std::optional<InterfaceTransmissionSemantics> transmission_semantics{};
  std::optional<SpectralResponse> reflectivity{};
};

struct ImportedOpticalModel {
  ModelProvenance provenance;
  std::vector<OpticalSurfaceRecord> surfaces;
  std::uint32_t max_interactions{8};
  std::vector<OpticalMaterialRecord> materials{};
  std::uint32_t entrance_medium_id{kVacuumMediumId};
};

class CompiledOpticalModel {
public:
  CompiledOpticalModel(CompiledOpticalModel &&) = default;
  CompiledOpticalModel(const CompiledOpticalModel &) = default;

  [[nodiscard]] const ModelProvenance &provenance() const { return provenance_; }
  [[nodiscard]] const std::vector<OpticalSurfaceRecord> &surfaces() const { return surfaces_; }
  [[nodiscard]] std::uint32_t max_interactions() const { return max_interactions_; }
  [[nodiscard]] std::uint32_t entrance_medium_id() const { return entrance_medium_id_; }
  [[nodiscard]] const std::vector<OpticalMaterialRecord> &materials() const { return materials_; }
  [[nodiscard]] std::optional<ResolvedOpticalMaterial> material_at(std::uint32_t id,
                                                                   double wavelength) const {
    if (id == kVacuumMediumId)
      return ResolvedOpticalMaterial{1, 1, 0};
    for (const auto &material : materials_) {
      if (material.id != id)
        continue;
      const auto phase = material.phase_index.view().interpolate_unchecked(wavelength);
      const auto group = material.group_index.view().interpolate_unchecked(wavelength);
      const auto absorption = material.absorption_per_m.view().interpolate_unchecked(wavelength);
      if (phase && group && absorption)
        return ResolvedOpticalMaterial{*phase, *group, *absorption};
      break;
    }
    return std::nullopt;
  }

private:
  friend std::optional<CompiledOpticalModel> compile_optical_model(const ImportedOpticalModel &);
  CompiledOpticalModel(ModelProvenance provenance, std::vector<OpticalSurfaceRecord> surfaces,
                       std::uint32_t max_interactions, std::vector<OpticalMaterialRecord> materials,
                       std::uint32_t entrance_medium_id)
      : provenance_(std::move(provenance)), surfaces_(std::move(surfaces)),
        max_interactions_(max_interactions), materials_(std::move(materials)),
        entrance_medium_id_(entrance_medium_id) {}

  ModelProvenance provenance_;
  std::vector<OpticalSurfaceRecord> surfaces_;
  std::uint32_t max_interactions_{};
  std::vector<OpticalMaterialRecord> materials_;
  std::uint32_t entrance_medium_id_{};
};

[[nodiscard]] inline std::optional<CompiledOpticalModel>
compile_optical_model(const ImportedOpticalModel &input) {
  if (!has_valid_provenance(input.provenance) || input.surfaces.empty() ||
      input.max_interactions == 0 || input.max_interactions > kMaximumOpticalModelInteractions)
    return std::nullopt;
  std::unordered_set<std::uint32_t> material_ids{kVacuumMediumId};
  for (const auto &material : input.materials)
    if (!material.is_valid() || !material_ids.insert(material.id).second)
      return std::nullopt;
  if (!material_ids.contains(input.entrance_medium_id))
    return std::nullopt;
  std::unordered_set<std::uint32_t> ids{kVacuumMediumId};
  ids.reserve(input.surfaces.size());
  for (const auto &surface : input.surfaces) {
    if (!surface.frame.is_valid() || !std::isfinite(surface.diameter_m) ||
        surface.diameter_m <= kEpsilon || surface.shape > FacetShape::hexagon_flat_x ||
        surface.role > SurfaceRole::refractive_interface || !ids.insert(surface.id).second ||
        surface.material_id != kVacuumMediumId)
      return std::nullopt;
    if (surface.transmission.has_value() != surface.transmission_semantics.has_value() ||
        (surface.transmission && !surface.transmission->is_valid()) ||
        (surface.reflectivity &&
         (surface.role != SurfaceRole::mirror || !surface.reflectivity->is_valid())))
      return std::nullopt;
    if (surface.role == SurfaceRole::refractive_interface) {
      if (surface.front_medium_id == surface.back_medium_id ||
          !material_ids.contains(surface.front_medium_id) ||
          !material_ids.contains(surface.back_medium_id))
        return std::nullopt;
      if (surface.transmission_semantics &&
          *surface.transmission_semantics != InterfaceTransmissionSemantics::additional_coating &&
          *surface.transmission_semantics != InterfaceTransmissionSemantics::complete_interface)
        return std::nullopt;
    } else if (surface.front_medium_id != kVacuumMediumId ||
               surface.back_medium_id != kVacuumMediumId || surface.transmission)
      return std::nullopt;
    if (!std::isfinite(surface.inner_radius_m) || surface.inner_radius_m < 0.0 ||
        (surface.sag && (surface.shape != FacetShape::circle ||
                         !is_valid(AxisymmetricMirror{0.0, surface.inner_radius_m,
                                                      surface.diameter_m / 2.0, *surface.sag}))) ||
        (!surface.sag && surface.inner_radius_m != 0.0))
      return std::nullopt;
  }
  return CompiledOpticalModel{input.provenance, input.surfaces, input.max_interactions,
                              input.materials, input.entrance_medium_id};
}

struct OpticalSurfaceHit {
  std::uint32_t surface_id{};
  SurfaceRole role{SurfaceRole::obscurer};
  double distance_m{};
  Vec3 point_m{};
  Vec3 normal{};
  std::size_t surface_index{};
};

[[nodiscard]] inline std::optional<OpticalSurfaceHit>
intersect_nearest_surface(const Ray &ray, const CompiledOpticalModel &optical_model) {
  std::optional<OpticalSurfaceHit> nearest;
  for (std::size_t surface_index = 0; surface_index < optical_model.surfaces().size();
       ++surface_index) {
    const auto &surface = optical_model.surfaces()[surface_index];
    std::optional<OpticalSurfaceHit> candidate;
    if (surface.sag) {
      const auto local_hit = intersect_axisymmetric_mirror(
          surface.frame.ray_from_parent(ray),
          AxisymmetricMirror{0.0, surface.inner_radius_m, surface.diameter_m / 2.0, *surface.sag});
      if (local_hit)
        candidate = OpticalSurfaceHit{surface.id, surface.role, local_hit->distance_m,
                                      surface.frame.point_to_parent(local_hit->point_m),
                                      surface.frame.direction_to_parent(local_hit->unit_normal)};
    } else {
      const ImportedDetectorSurface plane{surface.id,           surface.frame.origin_m,
                                          surface.frame.z_axis, surface.diameter_m,
                                          surface.shape,        surface.frame.x_axis};
      const auto hit = intersect_detector_surface_unchecked(ray, plane);
      if (hit)
        candidate = OpticalSurfaceHit{surface.id, surface.role, hit->distance_m, hit->point_m,
                                      hit->unit_normal};
    }
    if (candidate)
      candidate->surface_index = surface_index;
    if (candidate &&
        (!nearest || candidate->distance_m < nearest->distance_m ||
         (candidate->distance_m == nearest->distance_m && surface.id < nearest->surface_id)))
      nearest = candidate;
  }
  return nearest;
}

} // namespace obdeect
