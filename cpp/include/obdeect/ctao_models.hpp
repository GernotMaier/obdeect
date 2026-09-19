#pragma once

#include "obdeect/ctao_optical_specs.hpp"

#include <optional>
#include <string_view>

namespace obdeect {

// These are compact, redistributable *reference optical prescriptions*. They
// deliberately do not replace the licensed simulation-models repository:
// segment positions, alignment and camera pixels must be imported separately.
enum class CtaoTelescopeType : std::uint8_t { lst, mst, sst, sct };

[[nodiscard]] inline std::string_view to_string(CtaoTelescopeType type) {
  switch (type) {
    case CtaoTelescopeType::lst: return "LST";
    case CtaoTelescopeType::mst: return "MST";
    case CtaoTelescopeType::sst: return "SST";
    case CtaoTelescopeType::sct: return "SCT";
  }
  return "unknown";
}

struct CtaoReferenceModel {
  CtaoTelescopeType type{};
  TelescopeOpticalFamily family{};
  std::string_view identifier;
  std::string_view simulation_models_version{"6.3.0"};
  double primary_outer_radius_m{};
  double primary_inner_radius_m{};
  double primary_vertex_z_m{};
  // A spherical primary is used for the rotationally symmetric MST reference
  // baseline. Production modified Davies--Cotton tracing must replace it by
  // the imported facet list and facet normals.
  std::optional<double> primary_sphere_radius_m;
  EvenPolynomialSurface primary_surface{};
  std::optional<AxisymmetricMirror> secondary;
  double focal_plane_z_m{};
  double focal_plane_radius_m{};
};

[[nodiscard]] inline bool is_valid(const CtaoReferenceModel& model) {
  const bool common = !model.identifier.empty() && !model.simulation_models_version.empty() &&
                      std::isfinite(model.primary_outer_radius_m) &&
                      std::isfinite(model.primary_inner_radius_m) && std::isfinite(model.primary_vertex_z_m) &&
                      std::isfinite(model.focal_plane_z_m) && std::isfinite(model.focal_plane_radius_m) &&
                      model.primary_inner_radius_m >= 0.0 &&
                      model.primary_outer_radius_m > model.primary_inner_radius_m &&
                      model.focal_plane_radius_m > 0.0;
  if (!common) return false;
  if (model.primary_sphere_radius_m && (!std::isfinite(*model.primary_sphere_radius_m) ||
                                        *model.primary_sphere_radius_m <= kEpsilon)) return false;
  return !model.secondary || is_valid(*model.secondary);
}

[[nodiscard]] inline CtaoReferenceModel lst_reference_model() {
  return {CtaoTelescopeType::lst, TelescopeOpticalFamily::lst_parabolic, "LSTN-design", "6.3.0",
          11.5, 0.0, 0.0, std::nullopt, paraboloid(kLstNorthDesign.focal_length_m), std::nullopt,
          kLstNorthDesign.focal_length_m, 1.5};
}

[[nodiscard]] inline CtaoReferenceModel mst_reference_model() {
  // Central-facet curvature R=2 f is the validated first-order focus. This is
  // intentionally not a claim that a continuous sphere reproduces a segmented
  // modified Davies--Cotton telescope away from the optical axis.
  return {CtaoTelescopeType::mst, TelescopeOpticalFamily::mst_modified_davies_cotton,
          "MSTx-NectarCam", "6.3.0", 6.0, 0.0, 0.0, 2.0 * kMstNectarCam.focal_length_m,
          {}, std::nullopt, kMstNectarCam.focal_length_m, 1.0};
}

[[nodiscard]] inline CtaoReferenceModel sst_reference_model() {
  const auto reference = ssts_design_reference();
  auto secondary_surface = reference.secondary;
  secondary_surface.coefficient_m[0] = 0.0;
  return {CtaoTelescopeType::sst, TelescopeOpticalFamily::sst_sc, "SSTS-design", "6.3.0",
          reference.primary_diameter_m / 2.0, 0.0, 0.0, std::nullopt, reference.primary,
          AxisymmetricMirror{reference.secondary.coefficient_m[0], 0.0,
                              reference.secondary_diameter_m / 2.0, secondary_surface},
          reference.focal_length_m, 0.30};
}

[[nodiscard]] inline CtaoReferenceModel sct_reference_model() {
  const auto reference = scts_design_reference();
  auto secondary_surface = reference.secondary;
  secondary_surface.coefficient_m[0] = 0.0;
  return {CtaoTelescopeType::sct, TelescopeOpticalFamily::sct_sc, "SCTS-design", "6.3.0",
          reference.primary_diameter_m / 2.0, 0.0, 0.0, std::nullopt, reference.primary,
          AxisymmetricMirror{reference.secondary.coefficient_m[0], 0.0,
                              reference.secondary_diameter_m / 2.0, secondary_surface},
          reference.focal_length_m, 0.40};
}

[[nodiscard]] inline std::optional<CtaoReferenceModel> ctao_reference_model(std::string_view name) {
  if (name == "LST" || name == "lst" || name == "LSTN-design") return lst_reference_model();
  if (name == "MST" || name == "mst" || name == "MSTx-NectarCam") return mst_reference_model();
  if (name == "SST" || name == "sst" || name == "SSTS-design") return sst_reference_model();
  if (name == "SCT" || name == "sct" || name == "SCTS-design") return sct_reference_model();
  return std::nullopt;
}

}  // namespace obdeect
